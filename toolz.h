#ifndef CAPMON_TOOLZ_H
#define CAPMON_TOOLZ_H
#include <sstream>
#include <cassert>
#include <curl/curl.h>
#include <mysql/mysql.h>
#include <sqlite3.h> 
#include <string> 
#include <cstring>
#include <algorithm>
#include <vector>
// #include <type_traits>
#include <sys/stat.h> 
#include <stdexcept> // brought in by exceptions.h
#include "exceptions.h" 
#include "config.h"

typedef std::vector<std::string> ROW;
typedef std::vector<std::vector<std::string>> TABLE;

#ifndef NSUBMISSION
#include "json/json.h"
#endif 

#define STRING_LENGTH 250
#define SHORT_STRING 100
#define MAGIC_NUMBER 2500 

#define SEP1 "','"

// really must run this stuff through valgrind?!?

// this is really, really nasty!?
#ifndef NSUBMISSION
struct SUBMISSION { // really ought to make it at least a little private 

    int sbm_id;
    std::string submitter_name;
    std::string user;
    std::string user_email;
    std::string species;
    std::string file_name;
    std::string file_md5;
    std::string file_type;
    std::string file_desc;
    std::string timestamp;
    unsigned long file_size;
    int sbm_status; // ought to be an enum?!

    explicit SUBMISSION (const Json::Value& sbm) : 
      sbm_id(atoi(sbm["sbm_id"].asCString())),

      // postgres frontend changes...
      //user(sbm["user"].asString()), 
      //user_email(sbm["user_email"].asString()), 
      user(sbm["name"].asString()), 
      user_email(sbm["name_email"].asString()), 
      species(sbm["species_version"].asString()),
      file_name(sbm["file_name"].asString()), 
      file_md5(sbm["file_md5"].asString()),
      file_type(sbm["file_type"].asString()), 
      timestamp(sbm["timestamp"].asString()), 
      file_size(atoi(sbm["file_size"].asCString())) {

        std::string sbmname = sbm["submitter_name"].asString();
        sbmname.erase(std::remove(sbmname.begin(),sbmname.end(),'\''),sbmname.end());
        sbmname.erase(std::remove(sbmname.begin(),sbmname.end(),'"'),sbmname.end());
        this->submitter_name=sbmname;

        std::string desc = sbm["file_desc"].asString();
        desc.erase(std::remove(desc.begin(),desc.end(),'\''),desc.end());
        desc.erase(std::remove(desc.begin(),desc.end(),'"'),desc.end());
        file_desc = desc;

    }

    explicit SUBMISSION (ROW& r) : 
      sbm_id(std::stoi(r.at(2))), 
      file_name(r.at(0)), 
      file_md5(r.at(1)), 
      file_type(r.at(3)), 
      file_size(std::stol(r.at(4))),
      sbm_status(std::stol(r.at(5))) { /* construct */ }

    SUBMISSION (
        int _sbm_id,
        std::string _a,
        std::string _b,
        std::string _c,
        std::string _d,
        std::string _e
    ) : sbm_id(_sbm_id),
    submitter_name(_a),
    user_email(_b),
    species(_c),
    file_name(_d),
    file_type(_e) {};

    std::string to_sqlstring() { // horrible!?

        std:: stringstream ss;
        ss << "insert into cap_files (" 
            <<" sbm_id,              sbm_status,         submitter_name, "
            <<" user,                user_email,         species, "
            <<" file_name,           file_type,          file_md5, "
            <<" file_desc,           file_size,          rts "
            <<") values ("
            <<this->sbm_id<<","     <<FILE_NEW<<",'"            <<submitter_name<<SEP1
            <<user<<SEP1            <<user_email<<SEP1          <<species<<SEP1
            <<file_name<<SEP1       <<file_type<<SEP1           <<file_md5<<SEP1
            <<file_desc<<"',"       <<file_size<<",'"           <<timestamp<<"')";

        std::string blah(ss.str()); 
        return std::move(blah);

    }

};

typedef std::vector<SUBMISSION> SUBMISSIONLIST;

#endif

inline size_t my_curl_write(void * ptr, size_t, size_t, std::stringstream & stream) {
    std::string blah((char*)ptr);
    stream << blah;
    return blah.size(); // must return length or truncates
} 

inline void set_cap_status (sqlite3 * db, int status) {
    char buf[STRING_LENGTH];
    sprintf(buf,SQL_CAP_STATUS,status);
    char * zErrMsg = NULL;
    if(sqlite3_exec(db, buf, NULL, NULL, &zErrMsg)!=SQLITE_OK)
      throw std::runtime_error(zErrMsg);
    sqlite3_free(zErrMsg);
}

struct SPECIES {
    std::string name;
    DB_PARAMS dbp;
    std::string prefix;
    unsigned int max_id;
    SPECIES(std::string _n, std::string _d, std::string _p = "", unsigned int _m = 0) : name(_n), dbp(_d), prefix(_p), max_id(_m) {} 
    SPECIES()=delete;
};

namespace toolz {
    
inline bool file_exists(const char* dbname) {
    struct stat sts;
    bool x = (stat(dbname, &sts) == -1 && errno == ENOENT);
    return !x;
}

inline static int internal_meta_callback(void *blah, int, char **argv, char **){
    std::map<std::string,std::string>* _meta = static_cast<std::map<std::string,std::string>*>(blah);
    _meta->insert(std::pair<std::string,std::string>(argv[0],argv[1])); // _meta[argv[0]]=argv[1];
    return 0;
}

/* class ADAPTOR {
    virtual ~ADAPTOR = 0;
}; */

class MYSQL_ADAPTOR {

  // private:

    static MYSQL_ADAPTOR* me_ptr;
    MYSQL* conn;
    DB_PARAMS& dp;

    MYSQL_ADAPTOR(DB_PARAMS& _dp) : dp(_dp) { 

        conn = mysql_init(NULL);

        // very naughty...?!?
        const char * msg_str1 = "could not connect to species-specific e! cap database";
        size_t len1 = strlen(msg_str1);
        size_t len2 = strlen(dp.dbname());
        char * msg_full = (char*) malloc(len1 + len2 + 2);
        memcpy(msg_full, msg_str1, len1);
        msg_full[len1] = ' ';
        memcpy(msg_full+len1+1, dp.dbname(), len2+1);

        if(mysql_real_connect(conn,dp.host(),dp.user(),dp.pass(),dp.dbname(),dp.port(),NULL,0) == NULL) {
            std::string s(msg_full);
            free(msg_full);
            throw MySqlError(s);
        }

        free(msg_full);

    }

    ~MYSQL_ADAPTOR() {
        mysql_close(conn);
    }

    MYSQL_ADAPTOR(const MYSQL_ADAPTOR&) = delete;

  public:

    static MYSQL_ADAPTOR* get_instance() {
        if(!me_ptr) throw std::runtime_error("pointer to sqlite object is not initialised");
        return me_ptr;
    }
    
    static void connect(DB_PARAMS& dp) {
        if(!me_ptr) me_ptr = new MYSQL_ADAPTOR(dp);
    }

    static void disconnect() {
        if(me_ptr) delete me_ptr;
        me_ptr = NULL; 
        // return me_ptr;
    }

    int no_return_query(const char* q) {
        mysql_query(conn,q);
        return mysql_affected_rows(conn);
    }

    void func_helper2 (std::string q, int& x) { 

        MYSQL_RES *result;
        MYSQL_ROW row;

        if (mysql_query(conn,q.c_str()))
          throw std::runtime_error("unable to query mysql instance");

        if(!(result = mysql_store_result(conn)))
          throw std::runtime_error("result set pointer is null");
        else {
            row = mysql_fetch_row(result); 
            if (row == 0) throw std::runtime_error("no results for integer query " + q);
            x=atoi(row[0]);
        }
        mysql_free_result(result);
    }

    void func_helper2 (std::string q, std::string& x) { 

        MYSQL_RES *result;
        MYSQL_ROW row;

        if (mysql_query(conn,q.c_str()))
          throw std::runtime_error("unable to query mysql instance");

        if(!(result = mysql_store_result(conn)))
          throw std::runtime_error("result set pointer is null");
        else {
            row = mysql_fetch_row(result); 
            if (row == 0) throw std::runtime_error("no results for string query " + q);
            x=row[0];
        }
        mysql_free_result(result);
    }

    void func_helper2 (std::string q, TABLE& x) { 

        MYSQL_RES *result;
        MYSQL_ROW row;

        if (mysql_query(conn,q.c_str()))
          throw std::runtime_error("unable to query mysql instance");

        if(!(result = mysql_store_result(conn)))
          throw std::runtime_error("result set pointer is null");

        int num_fields = mysql_num_fields(result);

        while ((row = mysql_fetch_row(result))) {
            std::vector<std::string> s;
            for(int i = 0; i < num_fields; i++) s.push_back(row[i]);
            x.push_back(std::move(s));
        }

        mysql_free_result(result);
    }

    template <typename T> T generic_query(std::string y) {
        T x;
        func_helper2(y, x); 
        return std::move(x);
    }

};

class SQLITE_ADAPTOR { 

  // private:

    std::map<std::string,std::string> _meta;
    static SQLITE_ADAPTOR* me_ptr;
    sqlite3* database;
    char* db_name;
    sqlite3_stmt* statement;

    SQLITE_ADAPTOR(char* dbname) : db_name(dbname) { 
    
        if(!file_exists(dbname)) throw std::runtime_error("the db doesn't exist " + std::string(db_name));
        sqlite3_open(dbname, &database); 

    }

    ~SQLITE_ADAPTOR() { sqlite3_close(database); /* std::cout << "\ndtor! - just closed connection\n"; */ }
    SQLITE_ADAPTOR(const SQLITE_ADAPTOR&) = delete;

  public:

    static SQLITE_ADAPTOR* get_instance() {
        if(!me_ptr) throw std::runtime_error("pointer to sqlite object is not initialised");
        return me_ptr;
    }
    
    static void connect(char* dbname) {
        if(!me_ptr) me_ptr = new SQLITE_ADAPTOR(dbname);
        // return me_ptr;
    }

    static void disconnect() {
        if(me_ptr) delete me_ptr;
        me_ptr = NULL; 
        // return me_ptr;
    }

    void pull_files(const char*, int&) {
        std::cout << "this is the single int form\n";
    }

    void pull_files(const char*, std::string&) {
        std::cout << "this is the single string form\n";
    }

void func_helper2 (std::string q, double& x); 

void func_helper2 (std::string q, std::string& x) { 
    if(sqlite3_prepare_v2(database, q.c_str(), -1, &statement, 0) == SQLITE_OK) {
        if(sqlite3_step(statement)== SQLITE_ROW) x=reinterpret_cast<const char*>(sqlite3_column_text(statement, 0));
        else throw std::runtime_error("no rows!?!");
    } else {
        sqlite3_finalize(statement);
        throw std::runtime_error("error in response");
    }
    sqlite3_finalize(statement);
}

void func_helper2 (std::string q, int& x) { 

    if(sqlite3_prepare_v2(database, q.c_str(), -1, &statement, 0) == SQLITE_OK) {
        if(sqlite3_step(statement)== SQLITE_ROW) {
            x=sqlite3_column_int(statement, 0);
            // x=atoi(reinterpret_cast<const char*>(sqlite3_column_text(statement, 0)));
        } else throw std::runtime_error("no rows!?!");
    } else {
        sqlite3_finalize(statement);
        throw std::runtime_error("error in response");
    }
    sqlite3_finalize(statement);

}

void func_helper2 (std::string q, ROW& x) { 

    if(sqlite3_prepare_v2(database, q.c_str(), -1, &statement, 0) == SQLITE_OK) {
        int cols = sqlite3_column_count(statement);
        int result = sqlite3_step(statement);
        if(result == SQLITE_ROW) {
            std::vector<std::string> s;
            for(int col = 0; col < cols; col++) 
                x.push_back(reinterpret_cast<const char*>(sqlite3_column_text(statement, col)));
        } else throw std::runtime_error("error in response");
    } else {
        throw std::runtime_error("error in response");
    }

    sqlite3_finalize(statement);

}

void func_helper2 (std::string q, TABLE& x) {

    if(sqlite3_prepare_v2(database, q.c_str(), -1, &statement, 0) == SQLITE_OK) {
        int cols = sqlite3_column_count(statement);
        int result = 0;
        int row =0;
        while(true) {
            result = sqlite3_step(statement);
            if(result == SQLITE_ROW) {
                row++;
                std::vector<std::string> s;

                for(int col = 0; col < cols; col++)
                // basic_string::_S_construct NULL not valid : trying to create a string with a NULL char * pointer.
                  s.push_back(reinterpret_cast<const char*>(sqlite3_column_text(statement, col)));

                x.push_back(std::move(s));

            } else break;  
        }
    }

    sqlite3_finalize(statement);

}

template <typename T> T generic_query(std::string y) {
    T x;
    func_helper2(y, x); 
    return std::move(x);
}

template <typename T> T func_helper (T t, std::true_type) { 
    std::cout << t << " : true" << std::endl;
}
 
template <typename T> T func_helper (T t, std::false_type) { 
    std::cout << t << " : false" << std::endl;
}
 
template <typename T> T main_func(std::string y) { 
    T x;
    std::cout << "here we are = "<< y<<"\n";
    x = func_helper(x, typename std::is_integral<T>::type()); 
    return x;
}

void pull_files(const char*, std::vector<std::string>&) {

    if(sqlite3_prepare_v2(database, "select * from cap_files", -1, &statement, 0) == SQLITE_OK) {

        int cols = sqlite3_column_count(statement);
        int result = 0;
        int row =0;
        while(true) {
            result = sqlite3_step(statement);
            
            if(result == SQLITE_ROW) {
                row++;

                for(int col = 0; col < cols; col++) {
                    std::cout << "  [" << col<<"] " << sqlite3_column_text(statement, col);
                    std::string s = (char*)sqlite3_column_text(statement, col);
                }


                std::cout << "\n";

            } else break;  

            break;

        }
            
        sqlite3_finalize(statement);
    }
}

unsigned int check_capdb_status () {

    if (!me_ptr) throw std::runtime_error("you have to connect first!");
    int cap_status = me_ptr->generic_query<int>(SQL_CHECK_DB);

    if (cap_status==0) {
        std::cout << "Cap_status table appears to be empty\nexiting" << std::endl;
        exit(0);
    } else if (cap_status>0) {

        cap_status = me_ptr->generic_query<int>(SQL_QUERY_STATUS);

    } else throw Sqlite3Error("hmm");

    switch (cap_status) { 
    case CAP_STATUS_READY:
        std::cout << "Capdb has status 'ready'\n";
        std::cout << "Setting status to 'looping'\n";
        // get rid of this!?
        set_cap_status(database,CAP_STATUS_LOOPING);
        break;
    case CAP_STATUS_LOOPING:
        std::cout << "The database wasn't not shut down shut down properly - run the ready command" << std::endl;
        exit(0);
        break;
    case CAP_STATUS_LOCKED:
        std::cout << "Please reset database state with ready command" << std::endl;
        exit(0);
        break;
    default:
        std::cout << "Unknown status" << std::endl;
        exit(1);
    }

}

void update(const char* q) {
    char* zErrMsg = NULL;
    if(sqlite3_exec(database, q, NULL, NULL, &zErrMsg)!=SQLITE_OK) // throw std::runtime_error(zErrMsg);
      throw std::runtime_error(THROW_DEFAULT("unable to update sqlite3 db"));
    sqlite3_free(zErrMsg);

}

#ifndef NSUBMISSION

void get_local_file_list (SUBMISSIONLIST& sb) {

    char* buf = (char*)malloc((strlen(SQL_GET_LOCALFILES)+1));
    sprintf(buf,SQL_GET_LOCALFILES,FILE_LOCALSTORE); // not much point in snprintf?!?

    if(sqlite3_prepare_v2(database, buf, -1, &statement, 0) == SQLITE_OK) {

        int cols = sqlite3_column_count(statement);
        int result = 0;
        int row =0;
        while(true) {

            result = sqlite3_step(statement);
            if(result == SQLITE_ROW) {

                row++;
                int sbm_id = 0;
                std::string sbmname;
                std::string email;
                std::string species;
                std::string fname;
                std::string type;

                for(int col = 0; col < cols; col++) {

                    // if (col ==0) sbm_id = atoi((char*)sqlite3_column_text(statement, col));
                    // else if (col == 1) sbmname = (char*)sqlite3_column_text(statement, col);
                    switch (col) {
                        case 0:     sbm_id = atoi((char*)sqlite3_column_text(statement, col));  break;
                        case 1:     sbmname = (char*)sqlite3_column_text(statement, col);       break;
                        case 2:     email = (char*)sqlite3_column_text(statement, col);         break;
                        case 3:     species = (char*)sqlite3_column_text(statement, col);       break;
                        case 4:     fname = (char*)sqlite3_column_text(statement, col);         break;
                        case 5:     type = (char*)sqlite3_column_text(statement, col);          break;
                    }
                }

                SUBMISSION sbm(sbm_id, sbmname, email, species, fname, type);
                sb.push_back(sbm);

            } else break; 

        }

    }

    free(buf);
    sqlite3_finalize(statement);
}

#endif 

void update_file_status(int sbm_id, int status) {
    char* buf = (char*) malloc(STRING_LENGTH*sizeof(char));
    assert(buf!=0);
    snprintf(buf, STRING_LENGTH, SQL_SET_STATUS ,status,sbm_id);
    char * zErrMsg = NULL; 
    if(sqlite3_exec(database, buf, NULL, NULL, &zErrMsg)!=SQLITE_OK) throw std::runtime_error(zErrMsg);
    free(buf);
    sqlite3_free(zErrMsg);
}

void update_cap_status(int status) {
    char buf[STRING_LENGTH];
    sprintf(buf,SQL_CAP_STATUS,status);
    char * zErrMsg = NULL;
    if(sqlite3_exec(database, buf, NULL, NULL, &zErrMsg)!=SQLITE_OK)
    throw std::runtime_error(zErrMsg);
    sqlite3_free(zErrMsg);
}

void meta_insert(const char* k,const char* v) {
    char query[SHORT_STRING];
    sprintf(query,"insert into cap_meta values ('%s','%s');",k,v);
    char* zErrMsg = NULL;
    if(sqlite3_exec(database, query, NULL, NULL, &zErrMsg)!=SQLITE_OK) 
      throw Sqlite3Error(THROW_DEFAULT("unable to insert meta key/value pair into capdb"));
    sqlite3_free(zErrMsg);
}

void generic_insert(const std::string& q) { // it's really not appropriate to handle problems with normal flow - i.e. bool in this context - just throw?!?
    char* zErrMsg = NULL;
    if(sqlite3_exec(database, q.c_str(), NULL, NULL, &zErrMsg)!=SQLITE_OK) // throw std::runtime_error(zErrMsg);
        throw Sqlite3Error(THROW_DEFAULT("unable to insert into capdb"));
    sqlite3_free(zErrMsg);
}

void query_species(std::vector<SPECIES>& vs, std::string s="") { // just overload it - nah default args?!

    char buf[STRING_LENGTH];
    s = s.empty() ? "" : " where species like '%%" + s + "%%';";
    s = "select species,dsn,prefix,max_id from cap_species" + s;

    if(sqlite3_prepare_v2(database, s.c_str(), -1, &statement, 0) == SQLITE_OK) {

        while(sqlite3_step(statement)== SQLITE_ROW) {
            std::string species = (char*)sqlite3_column_text(statement,0);
            std::string dsn = (char*)sqlite3_column_text(statement,1);
        //try {
            if(sqlite3_column_type(statement,2)!=SQLITE_NULL) {
                std::string prefix = (char*)sqlite3_column_text(statement, 2);
                unsigned int max_id = (int)sqlite3_column_int(statement, 3);
                vs.push_back(SPECIES(species,dsn,prefix,max_id));
            } else  { 
                vs.push_back(SPECIES(species,dsn));
            }    
        // } catch (...) {
        // }

        } 
        
    } else throw std::runtime_error(buf);
    sqlite3_finalize(statement);
}

std::string meta_value(const char* key) {

    if(_meta.find(key)!=_meta.end()) // safe to use [] without generating entry
      return _meta[key];

    char q[SHORT_STRING];
    sprintf(q,"select value from cap_meta where key = '%s'",key);

    std::string local;
    if(sqlite3_prepare_v2(database, q, -1, &statement, 0) == SQLITE_OK) {

        if(sqlite3_step(statement)== SQLITE_ROW) {
            local=reinterpret_cast<const char*>(sqlite3_column_text(statement, 0));
        } else {
            // sqlite3_finalize(statement); // ?!?
            throw std::runtime_error("Cannot find entry for " + std::string(key)); // just use sprintf?!?
        }
            
    } else {
        sqlite3_finalize(statement);
        throw std::runtime_error("error in response");
    }

    // if(local=="local.perl5lib") std::for_each(perl5lib.begin(),perl5lib.end(),[](char& n){if(n==',')n=':';});

    sqlite3_finalize(statement);
    return std::move(local);

}

void meta_populate() {

    char *zErrMsg = 0; 
    int rc = sqlite3_exec(database, "select * from cap_meta", internal_meta_callback, &_meta, &zErrMsg);
    // int rc = sqlite3_exec(database, "select * from cap_meta", internal_meta_callback, 0, &zErrMsg);
    
    if(rc!=SQLITE_OK)
    fprintf(stderr, "SQL error: %s\n", zErrMsg);

}

bool species_info(std::string species, std::string& dsn, std::string& prefix, int& max_id) {

    char* buf = (char*)malloc(STRING_LENGTH); 
    snprintf(buf, STRING_LENGTH, "select dsn, prefix, max_id from cap_species where species = '%s'", species.c_str());
    
    if(sqlite3_prepare_v2(database, buf, -1, &statement, 0) == SQLITE_OK) {
        if(sqlite3_step(statement)== SQLITE_ROW) {
            dsn = (char*)sqlite3_column_text(statement, 0);
            if(sqlite3_column_type(statement,1)==SQLITE_NULL) prefix.clear();
            else  { 
                prefix = (char*)sqlite3_column_text(statement, 1);
                max_id = (int)sqlite3_column_int(statement, 2);
            }    
            // max_id = atoi((char*)sqlite3_column_text(statement, 2));
        } else return false; //throw std::runtime_error("no rows!?!");
    
    } else throw std::runtime_error(buf);

    free(buf);
    sqlite3_finalize(statement);

    return true;
}


};

class CURL_ADAPTOR {

    CURL* curl;
    bool xcurl;
    std::string cookie_file; 
    static CURL_ADAPTOR* me_ptr;

    CURL_ADAPTOR(bool xc,std::string ck) : xcurl(xc), cookie_file(ck) { 
        if(!xcurl) {
            curl_global_init(CURL_GLOBAL_ALL);
            curl = curl_easy_init();
        }
    }

    ~CURL_ADAPTOR() { 
        curl_easy_cleanup(curl); 
        if(me_ptr) delete me_ptr;
        me_ptr = NULL; 
    }

    CURL_ADAPTOR(const CURL_ADAPTOR&) = delete;

  public:

    void set_xcurl (bool b) { xcurl = b; }
    int get_xcurl () const { return xcurl; }

    static CURL_ADAPTOR* get_instance(bool xterm = false, std::string ck = COOKIE) {
        if(!me_ptr) me_ptr = new CURL_ADAPTOR(xterm, ck);
        return me_ptr;
    }

    void pull(std::string uri, std::stringstream& strstrm) {

        if (!xcurl) { // ND being a pain and have to supress SSL problems...?!
            curl_easy_setopt(this->curl,    CURLOPT_URL,                uri.c_str()); 
            curl_easy_setopt(curl,          CURLOPT_COOKIEFILE,         cookie_file.c_str()); // COOKIE);
            curl_easy_setopt(curl,          CURLOPT_COOKIEJAR,          cookie_file.c_str()); 
            curl_easy_setopt(curl,          CURLOPT_NOPROGRESS,         1);         
            curl_easy_setopt(curl,          CURLOPT_WRITEFUNCTION,      my_curl_write); 
            // curl_easy_setopt(curl,       CURLOPT_SSL_VERIFYPEER,      0); 
            curl_easy_setopt(curl,          CURLOPT_WRITEDATA,          &strstrm);      
            curl_easy_perform(curl);
        } else {
            char tmp[MAGIC_NUMBER];
            std::string cmd("curl -k -s -c " + cookie_file + " -b " + cookie_file + " ");
            cmd += uri;
            //std::cout << "using "<<cmd<<"\n";
            FILE * f2 = popen(cmd.c_str(), "r");
            while(( fgets( tmp, MAGIC_NUMBER, f2 )) != NULL ) strstrm << tmp;
            pclose(f2);
        }
    }

};


// odd multiple inclusion issue and not time to locate problem - make inline?!
inline int query_integer (char * query) {

    int db_return(0);
    sqlite3 * database;
    db_return = sqlite3_open("new_sqlite.db", &database);
    sqlite3_stmt * statement;

    if(sqlite3_prepare_v2(database, query, -1, &statement, 0) == SQLITE_OK) {
        if(sqlite3_step(statement)== SQLITE_ROW) {
            db_return = atoi((char*)sqlite3_column_text(statement, 0));
        } else throw std::runtime_error("no rows!?!");
    } else throw std::runtime_error("error in response");
    sqlite3_finalize(statement);

    return db_return;
}

}

#endif
