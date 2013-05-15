#include <mysql/mysql.h> 
#include <curl/curl.h>
#include "exceptions.h"
#include "config.h"
#include <getopt.h>
#include <sqlite3.h>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <boost/filesystem/operations.hpp>
#include "list.h"
#include "utils.h"
#include "toolz.h"
#include <limits>
#include "gff_validation.h"
#include "boost/regex.hpp"

class int2msg {
	int status;
public:
	int2msg(int width) : status(width) {}
	friend std::ostream& operator<<(std::ostream& os, const int2msg & fw) {
		std::string s;
		switch(fw.status) {
            case FILE_NEW:              s = "new";                  break;
            case FILE_ASCI:             s = "non_asci";             break;
            case FILE_EMPTY:            s = "file_empty";           break;
            case FILE_MD5:              s = "md5_error";            break;
            case FILE_MD5IGNORE:        s = "md5_ignore";           break;
            case FILE_LOCALSTORE:       s = "local_copy";           break; // downloaded correctly, awaiting processing - need an sbm specific gene removal thing?!?
            case FILE_TYPE:             s = "file_type";            break;
            case FILE_GFF3VAL1:         s = "gff_val1";             break;
            case FILE_GFF3VAL2:         s = "gff_val2";             break;
            case FILE_GFF3NONEWGENES:   s = "no_new_genes";         break;
            case FILE_NONCDS:           s = "non_coding_cds";       break;
            case FILE_DONE:             s = "done";                 break; // saved with new gene models!?!
            case FILE_IGNORE:           s = "ignore";               break;
		}
		return os << s;
	}
};

inline int argcp2int(const char * cp) {
	std::string str(cp);
	if (str=="new")                 return FILE_NEW;
	else if (str=="asci")           return FILE_ASCI;
	else if (str=="empty")          return FILE_EMPTY;
	else if (str=="md5")            return FILE_MD5;
	else if (str=="md5_ignore")      return FILE_MD5IGNORE;
	else if (str=="local")          return FILE_LOCALSTORE;
	else if (str=="type")           return FILE_TYPE;
	else if (str=="gff_val1")       return FILE_GFF3VAL1;
	else if (str=="gff_val2")       return FILE_GFF3VAL2;
	else if (str=="no_new_genes")   return FILE_GFF3NONEWGENES;  // remove GFF3 from this part?!?
	else if (str=="stop_codons")    return FILE_NONCDS;
	else if (str=="done")           return FILE_DONE;
	else if (str=="ignore")         return FILE_IGNORE;
	else throw std::runtime_error("value '" + std::string(cp) + "' not recognised"); 
}

#define SCHEMA_20SEP2012 "CREATE TABLE `cap_files` (\
  `sbm_id` INTEGER PRIMARY KEY AUTOINCREMENT,\
  `submitter_name` varchar(40) NOT NULL,\
  `user` varchar(15) NOT NULL,\
  `user_email` varchar(40) NOT NULL,\
  `species` varchar(40) NOT NULL,\
  `file_name` varchar(255) NOT NULL,\
  `file_type` varchar(40) NOT NULL,\
  `file_md5` varchar(40) NOT NULL,\
  `file_desc` varchar(255) NOT NULL,\
  `file_size` int(11) NOT NULL,\
  `rts` timestamp DATETIME DEFAULT null,\
  `ts` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP,\
  `sbm_status` varchar(40) DEFAULT NULL\
);\
CREATE TABLE `cap_meta` (       `key` varchar(40) PRIMARY KEY NOT NULL,       `value` varchar(255) NOT NULL );\
CREATE TABLE `cap_species` (\
  `species` varchar(40) PRIMARY KEY NOT NULL,\
  `dsn`     varchar(255) NOT NULL,\
  `prefix`  varchar(10) NULL,\
  `max_id`  int(11) NULL\
  `email`   varchar(80) DEFAULT \"none\" \
);\
CREATE TABLE `cap_status` (\
  `status_id` INTEGER PRIMARY KEY AUTOINCREMENT,\
  `ts` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP,\
  `status` int NOT NULL\
);\
insert into cap_meta values ('cap.schema','20SEP2012');"

// 
#define E_CAP_PATCH_08MAY2013_1 "DROP TABLE IF EXISTS `canon_crc64`;"

#define E_CAP_PATCH_08MAY2013_2 "CREATE TABLE `canon_crc64` (\
    `crc64` varchar(50) NOT NULL,\
    PRIMARY KEY (`crc64`)\
) ENGINE=InnoDB DEFAULT CHARSET=latin1;"

#define E_CAP_PATCH_08MAY2013_3 "DROP TABLE IF EXISTS `source_rank`;"

#define E_CAP_PATCH_08MAY2013_4 "CREATE TABLE `source_rank` (\
    `cap_source` enum('cap','flybase','batch') NOT NULL,\
    `cap_rank` int(10) NOT NULL,\
    PRIMARY KEY (`cap_source`)\
) ENGINE=MyISAM DEFAULT CHARSET=latin1;"


#define E_CAP_PATCH_08MAY2013_5 "alter table gene add column `cap_sbm_id`  int(11) DEFAULT NULL;"
#define E_CAP_PATCH_08MAY2013_6 "alter table gene add column `cap_orig_stable_id` varchar(256);"
#define E_CAP_PATCH_08MAY2013_7 "alter table gene add column `cap_status` enum('INTEGRATED','OBSOLETE','CURRENT') DEFAULT 'CURRENT';"
#define E_CAP_PATCH_08MAY2013_8 "alter table gene add column `cap_source` enum('cap','flybase','batch') NOT NULL;"
#define E_CAP_PATCH_08MAY2013_9 "alter table gene add unique index (cap_status,cap_sbm_id);"

#define HELP_LIST "\nLists submitted files as of last check of file server.\n\nOptions to return partial list:\n\n\
        --help      -h      This little message\n\
        --sbmid     -s      Submission id\n\
        --status    -p      Submission status\n\
        --species   -o      Filter by species\n\
        --user      -u      File server username\n\
        --submitter -n      Submitter name (wildcarded)\n\
        --email     -e      email address associated with username (wildcarded)\n"

#define HELP_RESET "\nReset the status of an individual submission. Requires -sbmid and -status arguments.\n"

#define SMALL_INDENT 20

union ARGS {
    std::string _s;
    int _i;
};

using std::string;
using std::cout;
using std::cin;
using std::endl;
using std::runtime_error;

typedef boost::regex regex;

static inline void ignore_exit(std::string s) {
    cout << s;
    exit(0);
}

inline sqlite3* sqliteopen (char* capdb) { // what's this still doing here?

    std::cout << "Checking cap db" << std::endl;
    sqlite3 * db = 0;
    if(sqlite3_open(capdb, &db)) {
        sqlite3_close(db);
        throw Sqlite3Error("couldn't open database");
    }

    return db;
}

static inline void grab_value(string& s, regex& reg, const char* v, string d="") {
    string prnt = !d.empty() ? " [" + d + "]" : "";
    cout << v << prnt << ": ";
    char* buffer = (char*)malloc(STRING_LENGTH+1); // hmm - prolly a bit silly?!?
    assert(buffer!=0);
    for(;;) { // while(1) {
        cin.getline(buffer,STRING_LENGTH); // cin >> dbname;
        if(strlen(buffer)==0 && !d.empty()) s=d; 
        else if(regex_match(buffer,reg)) s=buffer;
        else { 
            cout << "please enter a sensible "<<v<<" : ";
            continue;
        }
        break;
    }
    free(buffer);
    cout << "using "<<v<<" : "<< s << "\n";
}

static inline void grab_value(unsigned int& s, const char* v, unsigned int d=0) { // bit nasty
    string prnt = d!=0 ? " [" + std::to_string(d) + "]" : "";
    cout << v << prnt << ": ";
    char* buffer = (char*)malloc(STRING_LENGTH+1);
    static regex reg_int("\\d{1,5}");
    assert(buffer!=0);
    for(;;) { // while(1) {
        cin.getline(buffer,STRING_LENGTH); 
        if(strlen(buffer)==0 && d!=0) s=d; 
        else if(regex_match(buffer,reg_int)) s=atoi(buffer);
        else { 
            cout << "please enter a sensible "<<v<<" : ";
            continue;
        }
        break;
    }
    free(buffer);
    cout << "using "<<v<<" : "<< s << "\n";
}

namespace opt {
    static char* user;   
    static char* submitter; 
    static char* email; 
    static int sbmid = 0; 
    static char* status; 
    static char* capdb; 
    static char* species; 
    static bool add = false;
    static bool edit = false;
    static bool wipe = false;
    static bool xcurl = false;
    static bool pstdout = false;
}

enum { OPT_CAPDB=1 };

static struct option long_options[] = {
    {"help",            no_argument,        0,      'h'},
    {"user",            required_argument,  0,      'u'},
    {"sbmid",           required_argument,  0,      's'},
    {"species",         required_argument,  0,      'o'},
    {"capdb",           required_argument,  0,      OPT_CAPDB},
    {"email",           required_argument,  0,      'e'},
    {"submitter",       required_argument,  0,      'n'},
    {"externalcurl",    no_argument,        0,      'c'},
    {"stdout",          no_argument,        0,      'd'},
    {"status",          required_argument,  0,      'p'},
    {0,                 0,                  0,      0}
};

void parse_list_opts(int argc, char **argv, bool) {

    for (int c; (c = getopt_long (argc, argv, "dcho:u:s:n:e:p:",long_options, NULL));) {

        if (c == -1) break; // unsafe way - can overflow

        switch (c) {
            case 'u':           opt::user       = optarg;           break;
            case 'o':           opt::species    = optarg;           break;
            case OPT_CAPDB:     opt::capdb      = optarg;           break;
            case 'e':           opt::email      = optarg;           break;
            case 'n':           opt::submitter  = optarg;           break;
            case 'c':           opt::xcurl      = true;             break;
            case 'd':           opt::pstdout    = true;             break;
            case 's':           opt::sbmid      = atoi(optarg);     break;
            case 'p':           opt::status     = optarg;           break;
            case 'h':  
                    cout << HELP_LIST << endl;
                    exit(0);
            break;
            default:    abort ();
            }
        }

    if(opt::capdb == NULL)  opt::capdb  =    const_cast<char*>(SQLITE_DB_NAME); 
       
    if (optind < argc) {
        printf ("non-option ARGV-elements: ");
        while (optind < argc)
          printf ("%s ", argv[optind++]);
        putchar ('\n');
    }

}

void parse_config_opts(int argc, char **argv, bool) {

    for (int c; (c = getopt_long (argc, argv, "aedo:",long_options, NULL));) {
        if (c == -1) break; 
            switch (c) {
                case 'e':   opt::edit       = true;         break;
                case 'a':   opt::add        = true;         break;
                case 'd':   opt::wipe       = true;         break;
                case 'o':   opt::species    = optarg;       break;
                default:    abort ();
            }
        }

    if(opt::capdb==0)  opt::capdb  =    const_cast<char*>(SQLITE_DB_NAME); 
       
    if (optind < argc) {
        printf ("non-option ARGV-elements: ");
        while (optind < argc)
          printf ("%s ", argv[optind++]);
        putchar ('\n');
    }
}

#define FILE_ATTRIBES "sbm_id, \
species, \
file_type, \
ts, \
rts, \
sbm_status, \
submitter_name, user, user_email, \
file_md5, \
file_size, \
file_desc, \
file_name"

#define SEP "--------------------------"

const static char * LINE = SEP SEP SEP;
const static char * SQL_PRINT_FILES = "select " FILE_ATTRIBES " from cap_files";

void List (int ac, char ** av) {

    parse_list_opts(ac,av,true);

    if (!toolz::file_exists(opt::capdb)) {
        cout << "Database "<<opt::capdb<< " doesn't exist. Bye." <<endl; 
        exit(0);
    } 

    toolz::SQLITE_ADAPTOR::connect(opt::capdb);
    string sbmdir =  toolz::SQLITE_ADAPTOR::get_instance()->meta_value("local.sbmdir");
    toolz::SQLITE_ADAPTOR::disconnect();
    string sql_string(SQL_PRINT_FILES);

    if (opt::sbmid!=0)              sql_string += " where sbm_id = " + std::to_string(opt::sbmid); 
    else if (opt::user!=0)          sql_string += " where user = '" + std::string(opt::user) + "'";
    else if (opt::submitter!=0)     sql_string += " where submitter_name like '%%" + string(opt::submitter) + "%%'"; 
    else if (opt::email!=0)         sql_string += " where user_email like '%%" + string(opt::email) + "%%'"; 
    else if (opt::status!=0)        sql_string += " where sbm_status = " + std::to_string(argcp2int(opt::status));
   
    if (opt::species!=0&&opt::user==0&&opt::submitter==0&&opt::email==0&&opt::status==0)
      sql_string += " where species like '%%" + string(opt::species) + "%%';"; 
    else if (opt::species!=0)
      sql_string += " and species like '%%" + string(opt::species) + "%%';";
    else sql_string += ";"; 

    sqlite3* db = sqliteopen(opt::capdb); // why is this still here?!

    cout << "Querying db " << opt::capdb << endl;

    char* buf = (char*)malloc((sql_string.size()+1)*sizeof(char)); // char buf[STRING_LENGTH];
    assert(buf!=0);
    snprintf(buf, sql_string.size()+1, sql_string.c_str(), FILE_NEW);

    sqlite3_stmt * statement;

    string species;
    int status=0;
    const char* a="'";
    int row =0;

    if(sqlite3_prepare_v2(db, buf, -1, &statement, 0) == SQLITE_OK) {

        int cols = sqlite3_column_count(statement);
        int result = 0;

        while(true) {

            result = sqlite3_step(statement);
            if(result == SQLITE_ROW) {

                if (row==0) cout << LINE << endl;
                row++;

                for(int col = 0; col < cols; col++) {
                
                    string s;
                    cout << std::setw(22);

                    // clean this mess!?
                    switch (col) {

                        case 0:     cout << std::right << "Submission ID : " << std::left << std::setw(4) << std::setfill('0') << 
                                    std::right << sqlite3_column_text(statement, col) << std::setfill(' ');                             break;
                        case 1:     
                            species = reinterpret_cast<const char*>(sqlite3_column_text(statement, col));
                            cout << std::right << "Species : " << std::left << species;                                                 break;
                        case 2:     cout << std::right << "File Type : " << std::left << sqlite3_column_text(statement, col);           break;
                        case 3:     cout << std::right << "Downloaded at : " << std::left << sqlite3_column_text(statement, col);       break;
                        case 4:     cout << std::right << "Uploaded at : " << std::left << sqlite3_column_text(statement, col);         break;
                        case 5:     
                            status = sqlite3_column_int(statement, col);
                            cout << std::right << "Status : " << std::left << int2msg(status);                                          break;
                        case 6:     cout << std::right << "Submitter Name : " << std::left << sqlite3_column_text(statement, col);      break;
                        case 7:     cout << std::right << "User : " << std::left << sqlite3_column_text(statement, col);                break;
                        case 8:     cout << std::right << "Email : " << std::left << sqlite3_column_text(statement, col);               break;
                        case 9:     cout << std::right << "Md5 Sum : " << std::left << sqlite3_column_text(statement, col);             break;
                        case 10:    cout << std::right << "File Size (Bytes) : " << std::left <<sqlite3_column_text(statement, col);    break;
                        case 11:    cout << std::right << "Description : " << std::left << sqlite3_column_text(statement, col);         break;
                        case 12: {  
                            string filename(reinterpret_cast<const char*>(sqlite3_column_text(statement, col)));
                            filepath(filename);
                            cout << std::right << "File Name : " << std::left << sbmdir << "/" << filename << ".slim";                  break;
                        }
                    }

                    cout << "\n";
                }

                if (opt::sbmid==0) cout << LINE << endl;

            } else break;  
        }

    } else throw runtime_error("that was odd: '" + sql_string + a );

    free(buf);
    sqlite3_finalize(statement);
    sqlite3_close(db);

    // can't be bothered to re-write above so just wait and now grab species info...
    if (opt::sbmid!=0&&(status==FILE_DONE||status==FILE_NONCDS)){ // !=nullptr
        toolz::SQLITE_ADAPTOR::connect(opt::capdb);
        std::vector<SPECIES> vs;
        toolz::SQLITE_ADAPTOR::get_instance()->query_species(vs,string(species));

        if(vs.size()==0) ignore_exit(string("cannot find anything like ") + species);
        else if (vs.size()>1) ignore_exit(string("too many species entries ") + species);

        char buf[STRING_LENGTH];
        snprintf(buf,STRING_LENGTH,"select biotype,count(1) from gene where cap_source = 'cap' and cap_sbm_id = %d group by biotype;",opt::sbmid);
        toolz::MYSQL_ADAPTOR::connect(vs[0].dbp);
        auto table = toolz::MYSQL_ADAPTOR::get_instance()->generic_query<TABLE>(buf);
        toolz::MYSQL_ADAPTOR::disconnect();

        for_each(table.begin(), table.end(), [](ROW r){ 
            cout<<std::setw(19) << std::right << r.at(0) << " : " <<std::left<<r.at(1)<<"\n";
        });
    }

    if (opt::sbmid!=0) cout << LINE << endl;
    cout << "Entries=" << row << endl;

    return;
}

void Reset (int ac, char ** av) {

    parse_list_opts(ac,av,false);
    string sql_string(SQL_PRINT_FILES);

    if (opt::sbmid==0 || opt::status==0){ 
        cout <<  "\nMust provide both --sbmid and --status arguments.\n" << endl;
        exit(0);
    }

    string STATUS(opt::status);
    if (STATUS!="new"&&STATUS!="local"&&STATUS!="ignore"&&STATUS!="done"&&STATUS!="md5_ignore"){ 
        cout << "status must be one of: download, local, ignore, md5_ignore, done. exiting." << endl;
        exit(0);
    } 

    toolz::SQLITE_ADAPTOR::connect(opt::capdb);

    char buf[STRING_LENGTH];
    sprintf(buf,"select species, sbm_status from cap_files where sbm_id = %d;",opt::sbmid);
    auto row = toolz::SQLITE_ADAPTOR::get_instance()->generic_query<ROW>(buf);

    std::stringstream silly;
    silly << int2msg(stoi(row.at(1)));
    string cs(silly.str());
    if(opt::status==cs) {
        cout << "Entry already has status " << cs << "\n";
        exit(0);
    }

    cout << "Will set status of submission " << opt::sbmid << " for species '" << row.at(0) 
      << "' from status '" << cs << "' to '" << opt::status  << "'\n";

    cout << "Enter [yY] to continue :";
    char c;
    cin >> c;

    if(c!='y'&&c!='Y') { 
        cout << "Exiting\n";
        exit(0);
    }

    sprintf(buf,SQL_UPDATE_SBM,argcp2int(opt::status),opt::sbmid);
    toolz::SQLITE_ADAPTOR::get_instance()->update(buf);

    toolz::SQLITE_ADAPTOR::disconnect();

    cout << "Updated entry\n";

    return;

}


void Pull (int ac, char ** av) {

    parse_list_opts(ac,av,false);
    string sql_string(SQL_PRINT_FILES);

    if (opt::sbmid==0){ 
        cout <<  "\nMust provide --sbmid argument.\n" << endl;
        exit(0);
    }

    toolz::SQLITE_ADAPTOR::connect(opt::capdb);

    char buf[STRING_LENGTH];
    sprintf(buf,"select file_name from cap_files where sbm_id = %d;",opt::sbmid);
    auto filename = toolz::SQLITE_ADAPTOR::get_instance()->generic_query<string>(buf);
    string url = toolz::SQLITE_ADAPTOR::get_instance()->meta_value("remote.jsonurl");

    std::stringstream strstrm;
    toolz::CURL_ADAPTOR::get_instance()->pull(url + "type=annot:return=file:value=" + filename, strstrm); 

    toolz::SQLITE_ADAPTOR::disconnect();

    if(opt::pstdout) {
        cout << "printing "<<filename<<" to stdout\n";
        cout << strstrm.str();
    } else {
        cout << "printing "<<filename<<" to stderr\n";
        std::cerr << strstrm.str();
    }

    return;

}

void Species (int ac, char ** av) {

    parse_config_opts(ac,av,true);

    static regex reg_dbname("[a-zA-Z_][\\w_]+");
    static regex reg_blah("[\\w_\\.-]+");
    static regex reg_sv("[A-Z][a-z]+\\s[a-z]+\\s\\S+");

    if(opt::edit+opt::add+opt::wipe>1) { 
        cout << "Do you want to edit, add or delete an entry?\n";
        exit(0);
    } else if (opt::edit+opt::add+opt::wipe==0) {

        toolz::SQLITE_ADAPTOR::connect(opt::capdb);

        std::vector<SPECIES> vs;
        if (opt::species) toolz::SQLITE_ADAPTOR::get_instance()->query_species(vs,string(opt::species));
        else toolz::SQLITE_ADAPTOR::get_instance()->query_species(vs);

        int size = vs.size();
        if (size>0) cout<<LINE<<"\n";

        int indent = size==1 ? 32 : SMALL_INDENT;

        if (size==1) cout << std::string(18,' ') << "cap e! db summary.\n\n";

        for_each(vs.begin(), vs.end(), [size,indent](SPECIES s){ 
            cout << std::setw(indent) << std::right   <<"species_version : "  <<std::left<<s.name<<"\n";

            if (!s.prefix.empty()) 
              cout <<std::setw(indent) << std::right   <<"prefix : "   <<std::left<<s.prefix<<"\n"
                <<  std::setw(indent) << std::right   <<"max_id : "   <<std::left<<s.max_id<<"\n";

            cout << std::setw(indent) << std::right   <<"dbname : "   <<std::left<<s.dbp.dbname()<<"\n"
              <<    std::setw(indent) << std::right   <<"host : "     <<std::left<<s.dbp.host()<<"\n"
              <<    std::setw(indent) << std::right   <<"user : "     <<std::left<<s.dbp.user()<<"\n"
              <<    std::setw(indent) << std::right   <<"pass : "     <<std::left<<s.dbp.pass()<<"\n"
              <<    std::setw(indent) << std::right   <<"port : "     <<std::left<<s.dbp.port()<<"\n";

            if (size==1) {

                std::stringstream ssm;
                ssm << "select sbm_status,count(1) from cap_files where species = '" << s.name << "' group by sbm_status;";

                auto table = toolz::SQLITE_ADAPTOR::get_instance()->generic_query<TABLE>(ssm.str());
                sort(table.begin(),table.end(),[](ROW a,ROW b) { return a.at(0)>b.at(0); });
                if(table.size()>0) cout << "\n" << std::string(18,' ') << "submission status summary.\n\n";
                for_each(table.begin(), table.end(), [](ROW r){ 
                    cout<<std::setw(29) << std::right <<int2msg(stoi(r.at(0))) << " : " <<std::left<<r.at(1)<<"\n";
                });

                table.clear();
                toolz::MYSQL_ADAPTOR::connect(s.dbp);
                table = toolz::MYSQL_ADAPTOR::get_instance()->generic_query<TABLE>("select biotype,count(1) from gene where cap_source = 'cap' and cap_status = 'CURRENT' group by biotype;");

                if(table.size()>0) cout << "\n" << std::string(18,' ') << "cap models with cap_status 'current' biotypes.\n\n";
                for_each(table.begin(), table.end(), [](ROW r){ 
                    cout<<std::setw(29) << std::right << r.at(0) << " : " <<std::left<<r.at(1)<<"\n";
                });

                table.clear();
                table = toolz::MYSQL_ADAPTOR::get_instance()->generic_query<TABLE>("select cap_status,count(1) from gene where cap_source ='cap' group by cap_status;");

                if(table.size()>0) cout << "\n" << std::string(18,' ') << "genes by cap_status summary.\n\n";
                for_each(table.begin(), table.end(), [](ROW r){ 
                    cout<<std::setw(29) << std::right << r.at(0) << " : " <<std::left<<r.at(1)<<"\n";
                });

                toolz::MYSQL_ADAPTOR::disconnect();

            }

            cout<<LINE<<"\n";
        });

        toolz::SQLITE_ADAPTOR::disconnect();

        cout << "Entries=" << vs.size() << endl;
        exit(0);

    } else if (opt::edit) {

        toolz::SQLITE_ADAPTOR::connect(opt::capdb);

        // can't be bothered - just use a vector and allow matching?!?
        std::vector<SPECIES> vs;

        if(opt::species) toolz::SQLITE_ADAPTOR::get_instance()->query_species(vs,string(opt::species));
        else ignore_exit("you need to give a species");

        string dbname, host, user, pass, prefix;
        unsigned int port = 0, max_id = 0;

        if(vs.size()==0) ignore_exit(string("cannot find anything like ") + opt::species);
        else if (vs.size()>1) {
            cout << "you need to be more specific:\n";
            for_each(vs.begin(),vs.end(),[](SPECIES s) { cout << "* "<<s.name << "\n"; });
        } else {

            cout << "editing entry for " << vs[0].name << "\n";

            grab_value(dbname,reg_dbname,"dbname",vs[0].dbp.dbname());
            grab_value(host,reg_blah,"hostname",vs[0].dbp.host());
            grab_value(user,reg_blah,"user",vs[0].dbp.user());
            grab_value(port,"port",vs[0].dbp.port());
            grab_value(pass,reg_blah,"pass",vs[0].dbp.pass());
            grab_value(prefix,reg_blah,"prefix",vs[0].prefix);
            grab_value(max_id,"max_id",vs[0].max_id);

        }

        cout << "\ntrying to access db " << dbname << " for '"<< vs[0].name << "' using host="
          <<host<<" user="<<user<<" port="<<port<<"\n";

        string dsn = "DBI:mysql:database="+dbname+";host="+host+";port="+std::to_string(port)+","+user+","+pass;

        DB_PARAMS dbp(dsn);

        try {
            toolz::MYSQL_ADAPTOR::connect(dbp);
        } catch(std::runtime_error e) {
            cout << "unable to query " << dbname << " : " <<e.what() <<"\n";
            throw;
        }

        auto baseuri = toolz::SQLITE_ADAPTOR::get_instance()->meta_value("remote.jsonurl"); 

        char* buffer = (char*)malloc(STRING_LENGTH); 
        
        // eeew, dangerous, will write till end and leave out null byte if it overflows?!?
        snprintf(buffer,STRING_LENGTH,"update cap_species set dsn = '%s', prefix = '%s', max_id = %d where species = '%s'",dsn.c_str(),prefix.c_str(),max_id,vs[0].name.c_str());

        cout << "updating cap db\n";
        toolz::SQLITE_ADAPTOR::get_instance()->update(buffer);
        toolz::SQLITE_ADAPTOR::disconnect();

        int has_sv = toolz::MYSQL_ADAPTOR::get_instance()->generic_query<int>("select count(1) from meta where meta_key = 'cap.species_version'");

        switch (has_sv) {
            case 0: snprintf(buffer,STRING_LENGTH,"insert into meta (species_id,meta_key,meta_value) values (1,'cap.species_version','%s')", vs[0].name.c_str()); break;
            case 1: snprintf(buffer,STRING_LENGTH,"update meta set meta_value = '%s' where meta_key = 'cap.species_version'", vs[0].name.c_str()); break;
            default: ignore_exit("it appears the meta table is a bit of a mess?!?"); break;
        }

        cout << "updating cap e! cap_species_version\n";

        toolz::MYSQL_ADAPTOR::get_instance()->no_return_query(buffer);

        int has_bu = toolz::MYSQL_ADAPTOR::get_instance()->generic_query<int>("select count(1) from meta where meta_key = 'cap.base_uri'");

        switch (has_bu) {
            case 0: snprintf(buffer,STRING_LENGTH,"insert into meta (species_id,meta_key,meta_value) values (1,'cap.base_uri','%s')", baseuri.c_str()); break;
            case 1: snprintf(buffer,STRING_LENGTH,"update meta set meta_value = '%s' where meta_key = 'cap.base_uri'", baseuri.c_str()); break;
            default: ignore_exit("it appears the meta table is a bit of a mess?!?"); break;
        }

        cout << "updating cap e! binding to cap frontend at " << baseuri << "\n";

        toolz::MYSQL_ADAPTOR::get_instance()->no_return_query(buffer);

        free(buffer);
        toolz::MYSQL_ADAPTOR::disconnect();

        cout << "done.\n";
        exit(0);
                
    } else if (opt::add) {

        string species, dbname, host, user, pass, prefix;
        unsigned int port = 0, max_id = 0;

        grab_value(species,reg_sv,"species_version (must be of form Genus species VERSION) ");
        grab_value(dbname,reg_dbname,"dbname");
        grab_value(host,reg_blah,"hostname");
        grab_value(user,reg_blah,"user");
        grab_value(port,"port");
        grab_value(pass,reg_blah,"pass");
        grab_value(prefix,reg_blah,"prefix");
        grab_value(max_id,"max_id");
        toolz::SQLITE_ADAPTOR::connect(opt::capdb);
        auto baseuri = toolz::SQLITE_ADAPTOR::get_instance()->meta_value("remote.jsonurl"); 

        string dsn = "DBI:mysql:database="+dbname+";host="+host+";port="+std::to_string(port)+","+user+","+pass;
        DB_PARAMS dbp(dsn); 

        try {
            toolz::MYSQL_ADAPTOR::connect(dbp);
        } catch(std::runtime_error e) {
            cout << "unable to query " << dbname << " : " <<e.what() <<"\n";
            throw;
        }

        char* buffer = (char*)malloc(STRING_LENGTH); 

        snprintf(buffer,STRING_LENGTH,"insert into cap_species values('%s','%s','%s',%d)",species.c_str(),dsn.c_str(),prefix.c_str(),max_id);

        cout << "updating cap db\n";
        toolz::SQLITE_ADAPTOR::get_instance()->update(buffer);
        toolz::SQLITE_ADAPTOR::disconnect();

        int has_sv = toolz::MYSQL_ADAPTOR::get_instance()->generic_query<int>("select count(1) from meta where meta_key = 'cap.species_version'");
        switch (has_sv) {
            case 0: snprintf(buffer,STRING_LENGTH,"insert into meta (species_id,meta_key,meta_value) values (1,'cap.species_version','%s')", species.c_str()); break;
            case 1: snprintf(buffer,STRING_LENGTH,"update meta set meta_value = '%s' where meta_key = 'cap.species_version'", species.c_str()); break;
            default: ignore_exit("it appears the meta table is a bit of a mess?!?"); break;
        }
        cout << "updating cap e! cap_species_version\n";
        toolz::MYSQL_ADAPTOR::get_instance()->no_return_query(buffer);

        int has_bu = toolz::MYSQL_ADAPTOR::get_instance()->generic_query<int>("select count(1) from meta where meta_key = 'cap.base_uri'");
        switch (has_bu) {
            case 0: snprintf(buffer,STRING_LENGTH,"insert into meta (species_id,meta_key,meta_value) values (1,'cap.base_uri','%s')", baseuri.c_str()); break;
            case 1: snprintf(buffer,STRING_LENGTH,"update meta set meta_value = '%s' where meta_key = 'cap.base_uri'", baseuri.c_str()); break;
            default: ignore_exit("it appears the meta table is a bit of a mess?!?"); break;
        }
        cout << "updating cap e! binding to cap frontend at " << baseuri << "\n";
        
        toolz::MYSQL_ADAPTOR::get_instance()->no_return_query(buffer);

        free(buffer);

        cout << "attempting to patch e! db\n";
        toolz::MYSQL_ADAPTOR::get_instance()->no_return_query(E_CAP_PATCH_08MAY2013_1);
        toolz::MYSQL_ADAPTOR::get_instance()->no_return_query(E_CAP_PATCH_08MAY2013_2);
        toolz::MYSQL_ADAPTOR::get_instance()->no_return_query(E_CAP_PATCH_08MAY2013_3);
        toolz::MYSQL_ADAPTOR::get_instance()->no_return_query(E_CAP_PATCH_08MAY2013_4);
        toolz::MYSQL_ADAPTOR::get_instance()->no_return_query(E_CAP_PATCH_08MAY2013_5);
        toolz::MYSQL_ADAPTOR::get_instance()->no_return_query(E_CAP_PATCH_08MAY2013_6);
        toolz::MYSQL_ADAPTOR::get_instance()->no_return_query(E_CAP_PATCH_08MAY2013_7);
        toolz::MYSQL_ADAPTOR::get_instance()->no_return_query(E_CAP_PATCH_08MAY2013_8);
        toolz::MYSQL_ADAPTOR::get_instance()->no_return_query(E_CAP_PATCH_08MAY2013_9);

        toolz::MYSQL_ADAPTOR::disconnect();
        exit(0);

    } else {

        toolz::SQLITE_ADAPTOR::connect(opt::capdb);
        std::vector<SPECIES> vs;
        if(opt::species) toolz::SQLITE_ADAPTOR::get_instance()->query_species(vs,string(opt::species));
        else ignore_exit("you need to give a species");
        if(vs.size()==0) ignore_exit(string("cannot find anything like ") + opt::species);
        else if (vs.size()>1) {
            cout << "you need to be more specific:\n";
            for_each(vs.begin(),vs.end(),[](SPECIES s) { cout << "* "<<s.name << "\n"; });
        } else {
            cout << "Will delete entry for " <<vs[0].name << ". Are you sure you want to continue [yes] : ";
            string hmm;
            cin >> hmm;
            if(hmm=="yes") {
                char buffer[STRING_LENGTH];
                snprintf(buffer,STRING_LENGTH,"delete from cap_species where species = '%s'",vs[0].name.c_str());
                cout << "Make sure "<<vs[0].name <<" is not enabled as an gene model upload organism within frontend!\n";
            } else cout << "Will not remove entry without 'yes'\n";
        }

        exit(0);

    }

    cout << "edit="<<opt::edit<<"\n";
    cout << "add="<<opt::add<<"\n";
    exit(0);
}

void Config (int ac, char ** av) {

    parse_config_opts(ac,av,true);

    if(opt::edit&&opt::add) { 
        cout << "Do you want to edit OR add an entry?\n";
        exit(0);
    } else if (!opt::edit&&!opt::add) {

        toolz::SQLITE_ADAPTOR::connect(opt::capdb);

        auto table = toolz::SQLITE_ADAPTOR::get_instance()->generic_query<TABLE>("select key,value from cap_meta");

        if (table.size()>0) { 

            sort(table.begin(), table.end(), [](ROW a, ROW b)->bool{ return a<b; }); 
            cout<<LINE<<"\n";

            for_each(table.begin(), table.end(), [](ROW r){ 
                cout<<std::setw(16)<<std::right<<r.at(0)<<" : "<<std::left<<r.at(1)<<"\n";
            });
            cout<<LINE<<"\n";
        }

        cout << "Entries=" << table.size() << endl;

        toolz::SQLITE_ADAPTOR::disconnect();

        exit(0);
    }

    cout << "not implemented. bye." << endl; 
    exit(0);

}

void Setup (int ac, char ** av) {

    cout << "NEED TO PUT IN CHECKS FOR ALL META ENTRIES e.g. modules, curl... PRIOR TO INSERTING THEM?!?\n";

    parse_list_opts(ac,av,true);

    if (toolz::file_exists(opt::capdb)) {
        cout << "File "<<opt::capdb<< " alredy exist. Will not over-write. Exiting." <<endl; 
        exit(0);
    }

    sqlite3 * db = 0; 
    if(sqlite3_open(opt::capdb, &db)) {
        sqlite3_close(db); 
        throw Sqlite3Error("couldn't open database");
    } 

    cout << "installing schema for " << opt::capdb << endl;
    char * zErrMsg = NULL;
    if(sqlite3_exec(db, SCHEMA_20SEP2012, NULL, NULL, &zErrMsg)!=SQLITE_OK) 
      throw std::runtime_error(zErrMsg);
    sqlite3_free(zErrMsg);

    set_cap_status(db,CAP_STATUS_READY);
    sqlite3_close(db);

    { using std::cin;
    namespace fs = boost::filesystem; 

    toolz::SQLITE_ADAPTOR::connect(opt::capdb);

    std::string response;
    while(1) {
        cout << "Where would you like to store 'slim'd gff files? ";
        cin >> response;
        if(fs::exists(response)) {
            cout<<"'"<<response<<"' already exists - will not use.\n";
        } else {
            cout << "Creating directory '" << response << "' for local gff file storage\n";
            fs::create_directories(response);
            toolz::SQLITE_ADAPTOR::get_instance()->meta_insert("local.sbmdir",response.c_str());
            break;
        }
    }

    cout << "Please enter base url for json retrieval e.g. 'http://vectorbase-cap.ensemblgenomes.org/?q=uploaded'? ";
    cin >> response;
    toolz::SQLITE_ADAPTOR::get_instance()->meta_insert("remote.jsonurl",response.c_str());

    cout << "Please enter location cookie file for json retrieval? ";
    cin >> response;
    toolz::SQLITE_ADAPTOR::get_instance()->meta_insert("remote.cookie",response.c_str());

    cout << "Please enter admin email address? ";
    cin >> response;
    toolz::SQLITE_ADAPTOR::get_instance()->meta_insert("email.admin",response.c_str());

    cout << "Please enter to email address e.g. maillist for notification of file uploads? ";
    cin >> response;
    toolz::SQLITE_ADAPTOR::get_instance()->meta_insert("email.to",response.c_str());

    cout << "Please enter from email address e.g. 'noreply@blah.org'? ";
    cin >> response;
    toolz::SQLITE_ADAPTOR::get_instance()->meta_insert("email.from",response.c_str());
    
    cout << "Please enter remote smtp server? ";
    cin >> response;
    toolz::SQLITE_ADAPTOR::get_instance()->meta_insert("remoate.smtp",response.c_str());
    
    cout << "Please enter location of perl bin? ";
    cin >> response;
    toolz::SQLITE_ADAPTOR::get_instance()->meta_insert("local.perlbin",response.c_str());

    cout << "Please enter location root location of Bio::Perl, EnsEMBL & GffDoc as comma-delimited list? ";
    cin >> response;
    toolz::SQLITE_ADAPTOR::get_instance()->meta_insert("local.perl5lib",response.c_str());

    toolz::SQLITE_ADAPTOR::get_instance()->meta_insert("capmon.loglevel","info");
    toolz::SQLITE_ADAPTOR::get_instance()->meta_insert("capmon.restart","900");
    toolz::SQLITE_ADAPTOR::get_instance()->meta_insert("capmon.sleep","60");
    toolz::SQLITE_ADAPTOR::get_instance()->meta_insert("capmon.xcurl","0");

    }

    toolz::SQLITE_ADAPTOR::disconnect();

    return;
}

void Ready (int ac, char ** av) {

    parse_list_opts(ac,av,false);

    sqlite3* db = sqliteopen(opt::capdb);

    cout << "Resetting db status" << endl;
    set_cap_status(db,CAP_STATUS_READY);
    sqlite3_close(db);
    return;
}


void Expt (int, char**) {

    cout << "experimental stuff atm?!?\n";

    toolz::SQLITE_ADAPTOR::connect(opt::capdb);
    toolz::SQLITE_ADAPTOR::get_instance()->meta_populate();
    toolz::SQLITE_ADAPTOR::disconnect();
    static std::map<std::string,std::string> meta;

    class META_DATA {

        std::map<std::string,std::string> meta;

    };

}

void Validate (int ac, char ** av) {

    parse_list_opts(ac,av,true);

    if (!toolz::file_exists(opt::capdb)) {
        cout << "Database "<<opt::capdb<< " doesn't exist. Bye." <<endl; 
        exit(0);
    } 

    toolz::SQLITE_ADAPTOR::connect(opt::capdb);
    string sbmdir =  toolz::SQLITE_ADAPTOR::get_instance()->meta_value("local.sbmdir");
    
    string sql_string;
    if (opt::sbmid!=0){ 
        std::stringstream ssm;
        ssm << "select file_name,species,sbm_status from cap_files where sbm_id = " << opt::sbmid;
        sql_string += ssm.str();
    } else { 
        cout << "You must give a submission id.";
        exit(0);
    }

    auto row = toolz::SQLITE_ADAPTOR::get_instance()->generic_query<ROW>(sql_string);

    string filename = row.at(0);
    string species = row.at(1);

    sql_string = "select dsn from cap_species where species = '" + species + "'";
    auto dsn = toolz::SQLITE_ADAPTOR::get_instance()->generic_query<std::string>(sql_string);
    toolz::SQLITE_ADAPTOR::disconnect();

    filepath(filename);
    filename = sbmdir + "/" + filename + ".slim"; 
    cout << "Running validation on submission : " << opt::sbmid << "\nLocal file : " << filename 
      << "\n\nCurrent file status : " << int2msg(stoi(row.at(2))) << "\n";

    if (!toolz::file_exists(filename.c_str())) {
        cout << "there is no file "<< filename << "\n";
        exit(0);
    } 

    DB_PARAMS dbp(dsn);

    cout << "Gff Validation check : " << std::boolalpha << standalone_gff_validation (filename.c_str(),&dbp);

}

void Rebind (int ac, char ** av) {

    parse_config_opts(ac,av,true);

    toolz::SQLITE_ADAPTOR::connect(opt::capdb);
    std::string atual = toolz::SQLITE_ADAPTOR::get_instance()->meta_value("remote.jsonurl");

    cout << "will rebind cap instance currently bound to '" << atual << "'\n"

      "please provide an http(s) address providing the minimal interface : ";

    std::cin >> atual;

    if(atual.substr(0,4)!="http") cout << "'" << atual << " does not look like an url. exiting\n", exit(0);

    std::vector<SPECIES> vs;
    toolz::SQLITE_ADAPTOR::get_instance()->query_species(vs);

    // cout << "will rebind cap instance " << opt::capdb << " to " << atual << "\n";
    cout << "will rebind cap instance and " << vs.size() << " patch dbs :\n";

    for_each(vs.begin(), vs.end(),[](SPECIES s){cout << s.name << " : " << s.dbp.dbname()<<"\n";});

    cout << "to '" << atual << "'\nDo you wish to continue [Ny] ";

    std::string sql;
    cin >> sql;
    if(sql!="y"&&sql!="Y") exit(0);

    cout << "updating cap db\n";
    toolz::SQLITE_ADAPTOR::get_instance()->update(("update cap_meta set value ='"+atual+"' where key = 'remote.jsonurl'").data());
    toolz::SQLITE_ADAPTOR::disconnect();

    sql="update meta set meta_value = '"+atual+"' where meta_key = 'cap.base_uri'";
    for_each(vs.begin(), vs.end(),[&sql](SPECIES s){ 
        toolz::MYSQL_ADAPTOR::connect(s.dbp);
        cout << "updating " << s.name << " patch db " << s.dbp.dbname() << "\n";
        toolz::MYSQL_ADAPTOR::get_instance()->no_return_query(sql.data());
        toolz::MYSQL_ADAPTOR::disconnect();
    });
    cout << "all done\n";

    return;

}

