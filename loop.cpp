#include <fstream>
#include <getopt.h>
#include <openssl/md5.h> 
#include <boost/filesystem/operations.hpp>
#include "boost/regex.hpp"
#include <log4cpp/FileAppender.hh>
#include <log4cpp/OstreamAppender.hh>
#include <log4cpp/PatternLayout.hh>
#include "log4cpp/Category.hh"
#include <signal.h>   
#include <unistd.h>
#include <stdlib.h>
#include <sys/wait.h>
#include "email.h"
#include "utils.h"
#include "exceptions.h"
#include "toolz.h"
#include "levenshtein.h"
#include "json/json.h"
#include "time.h"
#include "gffdocwrap.h"  

/* MUST SET capmon TO INSERT HOSTNAME AND PROCESS ID OF PROCESS USING CAPDB - I.E. NOT JUST STATUS SO AS TO AVOID MULTIPLE INSTANCES WRITTING TO TABLES!?!? */

#define ANSI_COFF "\033[0m"
#define ANSI_REDB "\033[1;31m"
#define ANSI_RED "\033[0;31m"
#define ANSI_GREENB "\033[1;32m"
#define ANSI_GREEN "\033[0;32m"
#define ANSI_YELLOWB "\033[1;33m"
#define ANSI_YELLOW "\033[0;33m"
#define ANSI_BLUEB "\033[1;34m"
#define ANSI_BLUE "\033[0;34m"
#define ANSI_FUNKYB "\033[1;35m"
#define ANSI_FUNKY "\033[0;35m"
#define ANSI_CYANB "\033[1;36m" 
#define ANSI_CYAN "\033[0;36m" 
#define ANSI_WHITEB "\033[1;37m" 
#define ANSI_WHITE "\033[0;37m" 
#define CMD_GFFDOC "GffDoc.pl "
#define CMD_EFINGER "efingerprint "

#define CMD_GFFDOC_OPTS " -type contig=ignore -type match=ignore -type match_part=ignore \
-type pseudogenic_tRNA=mRNA:pseudogenic_tRNA -type ncRNA=mRNA:ncRNA \
-type tRNA=mRNA:tRNA -type miRNA=mRNA:miRNA -type rRNA=mRNA:rRNA \
-mRNA_biotype_dominant \
-non_protein_coding_types pseudogenic_tRNA \
-coordsystem toplevel -non_coding_cds -non_protein_coding -leafnonunique "

#define JSON_INCREMENT 100
#define LOGFILE "wrapper.log"
#define DEVLINE cout << "Dev : " << __LINE__ << endl;
#define SEP1 "','"

#define SQL_GET_FILES "select file_name, file_md5, sbm_id, file_type, file_size, sbm_status from cap_files where sbm_status == %d or sbm_status == %d"
#define SQL_GET_COUNT_BY_STATUS "select count(1) from cap_files where sbm_status == %d or sbm_status == %d"

using boost::filesystem::exists;
using std::runtime_error;
using std::stringstream;
using std::cout;
using std::endl;
using std::string;
using std::vector; 
// using std::pair; // using std::tuple; // using std::get; // tuple accessor?!?

typedef boost::regex regex;
typedef boost::smatch smatch; 

static const char *CORRECT_USAGE_MESSAGE =
"\nUsage: " EXEC " " COMMAND " [OPTIONS]\n"
"\n"
"      --help                       print this lovely message\n"
"\n"
"      --iterations     -i          set number of iterations to perform [default: none]\n"
"      --sleep          -s          turn on checking of das server\n"
"\n"
"\n"
"      -f, --from_addr              email address users receive from [default=\"cap_process_noreply@vectorbase.org\"]"
"\n"
"\n";

static unsigned char md5_result[MD5_DIGEST_LENGTH]; 
static std::map<std::string,std::string> meta;
static bool shchk = false;

enum { 
    OPT_UPLOADER = 1, OPT_SBMDIR, OPT_TEMPDIR, OPT_EXECDIR, 
    OPT_HELP, OPT_TEST, OPT_BASEURL, OPT_PERL, OPT_BIN, 
    OPT_LOGLEVEL, OPT_DAS, OPT_EMAILLIST, OPT_STARTEMAIL, OPT_CAPDB
};

static bool query_string(MYSQL*, const char*, std::string &);
static void generate_seqedits(DB_PARAMS&, std::string, int, std::string&);
static bool check_cap_db(log4cpp::Category&, DB_PARAMS&, SUBMISSIONLIST::iterator&, std::string&);
static bool check_cap_db_no_rows(log4cpp::Category&, DB_PARAMS&, char*, int);
static int get_html_table_for_genes(DB_PARAMS&, mailer*, std::string, int, std::string&, int);

void sig_action_function(int, siginfo_t *info, void*) { // inline...?!?
    std::cout << "\n---\nMessage received from child: " << (char*)info->si_value.sival_ptr<< "\n---\n" << std::endl; 
    // horrible?!?
    if (std::string((char*)info->si_value.sival_ptr)=="All done") {
        exit(0);
    } else if (std::string((char*)info->si_value.sival_ptr)=="Not safe to proceed") { 
        exit(1);
    }
}             

void filedir(std::string& filename, std::string dir) {

    static unsigned char md5_result[MD5_DIGEST_LENGTH];
    MD5((unsigned char*)filename.c_str(), filename.size(), md5_result); 
    char ca[3];
    std::string buf = std::move(dir);
    buf += "/";
    dir.clear();
    for(int i=0; i <2; i++) { 
        sprintf(ca,"%02X/",md5_result[i]);
        dir += ca;
    }                
    buf += dir;
    if (!exists(buf)) boost::filesystem::create_directories(buf);
    buf += filename;
    filename = std::move(buf);
    return;

}    

struct EMAILFNCTOR {
    mailer& mailobj;
    log4cpp::Category& log4;
    std::string subject;
    std::string body;
    EMAILFNCTOR(mailer& m, log4cpp::Category& l, std::string s, std::string b) : mailobj(m), log4(l),  subject(s), body(b) {};
    void operator() (std::string i) {
        log4.info("Copy email to : " + i);
        if (!mailobj.send_email(subject.c_str(), body.c_str(), i.c_str()))
          throw std::runtime_error(THROW_DEFAULT("couldn't send email"));
    }
};

inline static int meta_callback(void*, int, char **argv, char **){
    meta[argv[0]]=argv[1];
    return 0;
}

class LOOP {

    log4cpp::Category &log4;
    mailer &mailobj;
    int iterations;
    int nap;

  public:

    LOOP (log4cpp::Category &_log4, 
      mailer &_mailobj, 
      int _iterations, 
      int _nap
    ) : log4(_log4), 
      mailobj(_mailobj), 
      iterations(_iterations), 
      nap(_nap) {};

    void start ();
    void pull_remote_file_list ();
    void download_remote_files ();
    void process_local_file_list ();

};

bool update_file_status (sqlite3 * db, int sbm_id, int status);
std::string md5_sum2string(unsigned char* md);
void print_md5_sum(unsigned char* md);
bool fetch_dsn (sqlite3 * database, std::string species, std::string & dsn);
int query_integer (sqlite3 * database, const char * query);

void pull_new_files(CURL * curl_handle, sqlite3 * db, log4cpp::Category & log4);
void pull_file_list (CURL * curl_handle, sqlite3 * db, log4cpp::Category & log4);
void get_local_file_list (sqlite3 * db, SUBMISSIONLIST & sb);
void process_local_files(CURL * curl_handle, sqlite3 * db, log4cpp::Category & log4);

void checksqlite_forloop_populatemetamap (char*);

namespace opt {
    static int restart = 300;
    static int startemail = 0;
    static int realemail = 0;
    static int emaillist = 0;
    static int iterations = 0;
    static int max_retries = 10;
    static bool xcurl = false;
    static int forceretry = 0;
    static char *capdb;
    static char *loglevel;
    static char *xemails;
    std::vector<string> emailvector;
    static int sleep = 300;
}

namespace stash { 
    static string sbmdir;
    static string url;
    static string from_email;
    static string smtp;
    static string admin_email;
    static string bin;
    static string loglevel;
    //static int retstart = 900;
    //static int sleep = 60;
    //static bool xcurl = false;
}

namespace retry {
    static int retries = 0;
}

static struct option long_options[] = {

    {"realemail",           no_argument,        &opt::realemail,        1},
    {"forceretry",          no_argument,        &opt::forceretry,       1},
    {"externalcurl",        no_argument,        0,                      'c'},
    {"emaillist",           no_argument,        0,                      OPT_EMAILLIST},
    {"startemail",          no_argument,        0,                      OPT_STARTEMAIL}, 
    {"help",                no_argument,        0,                      OPT_HELP},
    {"extraemails",         required_argument,  0,                      'x'},
    {"iterations",          required_argument,  0,                      'i'},
    {"sleep",               required_argument,  0,                      's'},
    {"capdb",               required_argument,  0,                      OPT_CAPDB},
    {"bin",                 required_argument,  0,                      OPT_BIN},
    {"loglevel",            required_argument,  0,                      OPT_LOGLEVEL},
    {0,                     0,                  0,                      0}
};

void parse_opts(int argc, char **argv) {

    int index;
    for (int c; (c = getopt_long (argc, argv, "cx:f:edi:s:X:U:2:h:u:p:P:d:r:",long_options, &index));) {

        if (c == -1) break; // unsafe way - can overflow

        switch (c) {
            case 'i':               opt::iterations = atoi(optarg);     break;
            case 's':               opt::sleep      = atoi(optarg);     break;
            case 'c':               opt::xcurl      = true;             break;
            case 'x':               opt::xemails    = optarg;           break;
            // duh, always forget to add a mapping of value to option var?!?
            case OPT_LOGLEVEL:      opt::loglevel   = optarg;           break;
            case OPT_STARTEMAIL:    opt::startemail = 1;                break;
            case OPT_EMAILLIST:     opt::emaillist  = 1;                break;
            case OPT_HELP:
                std::cout << CORRECT_USAGE_MESSAGE;
                exit(0);
            default:    abort ();
        }
    }
   
    if(opt::loglevel == NULL)   opt::loglevel   =   const_cast<char*>("info");
    if(opt::xemails == NULL)    opt::xemails    =   const_cast<char*>(EMAIL_LISTSERVE);
    if(opt::capdb == NULL)      opt::capdb      =   const_cast<char*>(SQLITE_DB_NAME);

    if (optind < argc) {
        printf ("non-option ARGV-elements: ");
        while (optind < argc)
          printf ("'%s' ", argv[optind++]);
        putchar ('\n');
        std::cout << "exiting." << std::endl;
        exit(1);
    }
}

void check_perl_conf (char*) {

    toolz::SQLITE_ADAPTOR::connect(opt::capdb);
    auto dnslist = toolz::SQLITE_ADAPTOR::get_instance()->generic_query<TABLE>("select dsn from cap_species");

    std::vector<DB_PARAMS> dbs2;
    std::for_each(dnslist.begin(),dnslist.end(),[&dbs2](ROW r){ dbs2.push_back(r.at(0)); });

    std::string perl = "perl";
    std::cout << "Testing perl DBD::mysql and individual db connections to cap e! dbs"<<std::endl; 

    for (auto i = dbs2.begin() ; i!=dbs2.end() ; i++) { 

        toolz::MYSQL_ADAPTOR::connect(*i);
        std::string organism_name = toolz::MYSQL_ADAPTOR::get_instance()->generic_query<std::string>(SQL_GET_SPECIES);
        toolz::MYSQL_ADAPTOR::disconnect();

        std::string cmd(perl + " -MDBI -e 'print DBI->connect(q{DBI:mysql:" + i->dbname() + ";" + i->host() + ":"
          +  std::to_string(i->port()) + "},q{" + i->user() + "},q{" + i->pass() 
          + "})->selectrow_array(q{select meta_value from meta where meta_key =\"cap.species_version\"},undef)'"
        );

        if(system(cmd.c_str())!=0) 
          std::runtime_error("WE CANNOT CONTINUE!?! as the perl instance given is not able to query the db");

        FILE * f2 = popen(cmd.c_str(), "r"); 
        char tmp[STRING_SIZE];
        if (organism_name==std::string(fgets(tmp,STRING_SIZE,f2))) std::cout << " " << i->dbname() << "\n";
        else throw std::runtime_error("there is a problem with " + string(i->dbname()));
        pclose(f2);  

    }

    const char* modules[] = { "Bio::EnsEMBL::Analysis", "Bio::PrimarySeqI", "GffDoc::Exceptions" }; 

    string perl5lib = toolz::SQLITE_ADAPTOR::get_instance()->meta_value("local.perl5lib");

    for (const char** it=modules ; it!=modules+(sizeof(modules)/sizeof(char*)) ; it++) { 

        std::for_each(perl5lib.begin(),perl5lib.end(),[](char& n){if(n==',')n=':';}); 

        string cmd("PERL5LIB=" + perl5lib + " " + perl + " -M"+ *it +" -e 'print q{blah}' &> /dev/null");
        cout << "Checking " << *it << "\n";
        if(system(cmd.c_str())!=0) 
          throw std::runtime_error("There is a problem using module" + string(*it));
    }

    toolz::SQLITE_ADAPTOR::disconnect();
    return;
}

// this is horrible, horrible, horrible - need to use adaptor for all cases of curl/sqlite
inline void runthething(mailer& mailobj,log4cpp::Category& log4) {

    LOOP loop(log4,mailobj,opt::iterations,opt::sleep); 

    try { 

        log4.debug("starting loop");
        loop.start();

    } catch(runtime_error& e) { //} catch(MySqlConnError mse) {

        std::cout << "Houston we have a problem...\n" << "Threw a " << typeid(e).name() << " saying " << e.what() << "\nReseting db." << std::endl;
        toolz::SQLITE_ADAPTOR::connect(opt::capdb);
        toolz::SQLITE_ADAPTOR::get_instance()->update_cap_status(CAP_STATUS_READY);
        toolz::SQLITE_ADAPTOR::disconnect();
        std::cout << "Exiting." << std::endl;
        // too many emails?
        // mailobj.send_email(std::string("Houston we have a problem (" + std::string(typeid(e).name()) + "). Exiting.\r\n").c_str(),"blah\r\n", stash::admin_email.c_str());
        exit(0);

    } catch(std::logic_error& e) { 

        if(e.what()==std::string("basic_string::_S_create")) {
            std::cout << "oh dear : we have std::logic_error with basic_string::_S_create - seems we tried to generate a string longer than std::string::max_size()?!?\n";
            std::cout << "Exiting." << std::endl;
            mailobj.send_email("oh dear : we have std::logic_error with basic_string::_S_create - seems we tried to generate a string longer than std::string::max_size()?!?\r\n","blah\r\n", stash::admin_email.c_str());
        } else {

            std::cout << "WTF : " << "Threw a " << typeid(e).name() << " saying " << e.what() << "\nReseting db." << std::endl;
            toolz::SQLITE_ADAPTOR::connect(opt::capdb);
            toolz::SQLITE_ADAPTOR::get_instance()->update_cap_status(CAP_STATUS_READY);
            toolz::SQLITE_ADAPTOR::disconnect();
            std::cout << "Exiting." << std::endl;
            mailobj.send_email(std::string("WTF : (" + std::string(typeid(e).name()) + "). Exiting.\r\n").c_str(),"blah\r\n", stash::admin_email.c_str());

        }

        exit(0);

    } catch(std::exception& e) { 

        std::cout << "WTF : " << "Threw a " << typeid(e).name() << " saying " << e.what() << "\nReseting db." << std::endl;
        toolz::SQLITE_ADAPTOR::connect(opt::capdb);
        toolz::SQLITE_ADAPTOR::get_instance()->update_cap_status(CAP_STATUS_READY);
        toolz::SQLITE_ADAPTOR::disconnect();
        std::cout << "Exiting." << std::endl;
        mailobj.send_email(std::string("WTF : (" + std::string(typeid(e).name()) + "). Exiting.\r\n").c_str(),"blah\r\n", stash::admin_email.c_str());
        exit(0);

    } catch(...) { 

        std::cout << "Loop sent un-handled exception." << std::endl;
        exit(0);

    }

}

void Loop (int ac, char ** av) {

    parse_opts(ac,av);

    // why is this done like this?
    toolz::SQLITE_ADAPTOR::connect(opt::capdb);
    toolz::SQLITE_ADAPTOR::get_instance()->meta_populate();
    stash::admin_email = toolz::SQLITE_ADAPTOR::get_instance()->meta_value("email.admin");
    stash::sbmdir = toolz::SQLITE_ADAPTOR::get_instance()->meta_value("local.sbmdir");
    stash::url = toolz::SQLITE_ADAPTOR::get_instance()->meta_value("remote.jsonurl");
    stash::from_email = toolz::SQLITE_ADAPTOR::get_instance()->meta_value("email.from");
    stash::smtp = toolz::SQLITE_ADAPTOR::get_instance()->meta_value("remote.smtp");
    stash::bin = toolz::SQLITE_ADAPTOR::get_instance()->meta_value("local.bindir");
//    stash::loglevel = toolz::SQLITE_ADAPTOR::get_instance()->meta_value("capmon.loglevel");
//    stash::xcurl = stoi(toolz::SQLITE_ADAPTOR::get_instance()->meta_value("capmon.xcurl"));
//    stash::restart = stoi(toolz::SQLITE_ADAPTOR::get_instance()->meta_value("capmon.restart"));
//    stash::sleep = stoi(toolz::SQLITE_ADAPTOR::get_instance()->meta_value("capmon.sleep"));
    toolz::SQLITE_ADAPTOR::disconnect();

    mailer mailobj(stash::smtp.c_str(),stash::from_email.c_str());

    pid_t child_pid, wpid;
    int status = 0;

    struct sigaction act;
    memset (&act, '\0', sizeof(act));  
    act.sa_sigaction = sig_action_function;
    act.sa_flags = SA_SIGINFO;    
    sigaction(SIGUSR1, &act, 0); 
    sigaction(SIGUSR2, &act, 0); 

    check_perl_conf(opt::capdb);

    // eew
    char * shell;            
    shell = getenv ("SHELL");  
    if (shell!=NULL) {
        size_t len = strlen(shell);
        if(std::string(shell).substr(len-3,3)=="zsh" || std::string(shell).substr(len-4,4)=="bash") {
            std::cout << "Will use ANSI colours.\n";
            shchk = true;
        }
    }

    if(opt::realemail){
        if(shchk) std::cout << ANSI_WHITEB;
        std::cout << "\n\nYou have invoked with --realemail long-opt\n\nAre you sure [N|y]? ";
        char proceed = 'N';
        std::cin >> proceed;
        if (proceed=='y'||proceed=='Y') std::cout << "Starting in full mode." << std::endl;
        else {
            std::cout<< "exiting" << std::endl;
            exit(0);
        }
        std::cout << "Starting in 5s" << std::endl; 
        sleep(5);
        if(shchk) std::cout << ANSI_COFF;
    }    

    log4cpp::Appender* appender; 
    appender = new log4cpp::OstreamAppender("default", &std::cout);
    log4cpp::PatternLayout* patternLayout = new log4cpp::PatternLayout();

    log4cpp::Appender* appender_file; 
    appender_file = new log4cpp::FileAppender("FileAppender", LOG4_FILE);
    patternLayout->setConversionPattern("CAPMON : [%p] : %d{%Y/%m/%d %H:%M:%S} : %m\n"); //"%R %p %c %x: %m\n");
    appender->setLayout(patternLayout);

    log4cpp::Category& log4 = log4cpp::Category::getRoot();
    log4.addAppender(appender);
    log4.addAppender(appender_file);
    log4.setAppender(appender);

    if (std::string(opt::loglevel)=="debug") {
        std::cout << "Setting loglevel to debug" << std::endl;
        log4.setPriority(log4cpp::Priority::DEBUG);
    } else log4.setPriority(log4cpp::Priority::INFO);

    if(shchk) std::cout << ANSI_BLUEB;
                  
    stringstream ss(opt::xemails);
    std::string s;     
                  
    while (getline(ss, s, ',')) 
      opt::emailvector.push_back(s);
                  
    toolz::CURL_ADAPTOR::get_instance()->set_xcurl(opt::xcurl); 
    std::cout << "Sleep interval " << opt::sleep << "." << std::endl;
    std::cout << "External_Curl " << toolz::CURL_ADAPTOR::get_instance()->get_xcurl() << "." << std::endl;
    std::cout << "Additional email list " << opt::xemails << "." << std::endl;

    if(shchk) std::cout << ANSI_FUNKYB; 

    std::cout << "Storing files in '" << stash::sbmdir << "' directory." << std::endl;

    if(shchk) std::cout << ANSI_COFF;

    if(opt::iterations==1) {
        log4.info("Single iteration - will not fork");
        runthething(mailobj,log4);
        return;
    }

    while(1) { // for(;;) 

        child_pid = fork();
        if(shchk) std::cout << ANSI_WHITEB;
        std::cout << "Parent process pid=" << getpid() << " starting child fork" << std::endl; 

        switch (child_pid) { // if (child_pid==-1) {

            case -1:    std::cout << "Unable to fork!?! Exiting."<<endl;    exit(1);        break;

            case 0: { 

                if(shchk) std::cout << ANSI_COFF;
                std::cout << "Hello I'm child process pid=" << getpid() << std::endl; 

                runthething(mailobj,log4);

                char *messageText = const_cast<char*>("All done");
                union sigval signal_value;
                signal_value.sival_ptr = messageText;

                log4.info("Exiting Cap : Sending parent shut down signal");
                sigqueue(getppid(), SIGUSR2, signal_value); 

                return;

            } 

            default:        std::cout << "Parent fork pid="<<getpid()<<" generated child fork pid=" << child_pid << endl;  break;
        }

        while ((wpid = wait(&status)) > 0) {

            if (status>0) 
              log4.error("OUCH SEEMS WE EXCITED IN A MOST UNDIGNIFIED MANNER - SIGNAL : " + std::to_string(status)); 
            else 
              std::cout << "\nSeems child fork deliberately gave up.\n"<<endl;

            log4.warn("Waiting " + std::to_string((int)(opt::restart/60)) + "min before generating new child fork"); 
            sleep(opt::restart);
            mailobj.send_email("Creating new child fork.\r\n","If you receive many of these then you need to investigate!?!.\r\n", stash::admin_email.c_str()); 
        }
    }

    return;
}

bool query_string(MYSQL *conn, const char* query, std::string & stringy) {

    MYSQL_RES *result;
    MYSQL_ROW row;
    //MYSQL_FIELD *field;

    if (mysql_query(conn,query)) {
        stringy = "unable to query database : ";
        stringy += query;
        return false;
    }

    if(!(result = mysql_store_result(conn))) { 
        stringy = "resultset pointer is null : ";
        stringy += query;
        return false;
    } else {
        row = mysql_fetch_row(result);
        if (row == 0) {
            stringy = "row pointer is null : ";
            stringy += query;
            return false;
        }
        stringy = row[0];
    }

    mysql_free_result(result);

    return true;

}

size_t my_write_func2(void* ptr, size_t, size_t, stringstream & stream) {
    std::string blah((char*)ptr);
    stream << blah;
    return blah.size();
} 

void LOOP::start () {

    log4.info("Entering loop");

    int * date = day(); //y returns address of static int 
    int cnt = 0;
    int days = 0;

    for (;;) {

        if (cnt==0 && opt::startemail) {
            log4.debug("Sending startup email");
            if(!mailobj.send_email("Cap Monitor : started","You should receive one of the emails a day "
              "- else check the process is still running!?!\r\n", EMAIL_ADMIN)) 
              throw runtime_error("couldn't send email");
        } else if (*date != *(day())) {
            if(!mailobj.send_email(std::string("Daily I'm okay emails [" + std::to_string(++days) + "]").c_str(),
              "Just your daily reminder", EMAIL_ADMIN))
              throw runtime_error("couldn't send email");
        }

        // change this to the appropriate ostream manipulator?!?
        char msg[STRING_SIZE];
        sprintf(msg,"New file check - %05d",++cnt);
        log4.debug(msg);

        try {

            // wtf?!
            sqlite3 * db;
            int rc;
            rc = sqlite3_open(opt::capdb, &db);
            if(rc){ // why?!?
                fprintf(stderr, "Can't open database: %s\n", sqlite3_errmsg(db));
                sqlite3_close(db);
                throw Sqlite3Error("couldn't open database");
            }

            // why are you using this?!?
            pull_remote_file_list(); // this->pull_remote_file_list();
            (*this).download_remote_files();

            process_local_file_list();

        } catch(Sqlite3Error e) {

            log4.error(std::string("Sqlite3Error exception thrown : ") + e.what());
            ++retry::retries; 
            log4.error(std::string("Retry count=") + std::to_string(retry::retries));

            // let's not be throwing here!?!?
            if (retry::retries >= opt::max_retries) mailobj.send_email("Cap Monitor Giving Up", "sqlite problem.", EMAIL_ADMIN);
            if (retry::retries==1) mailobj.send_email("Cap Monitor Entering Retry Loop", "sqlite problem.", EMAIL_ADMIN);

            if (retry::retries >= opt::max_retries) throw; 
            log4.warn("Retry count below max retry value " + std::string(std::to_string(opt::max_retries)) + ". Will continue.");
        } catch (CapServerError e) {
            log4.error(std::string("CapServer exception thrown : ") + e.what());
            ++retry::retries; 
            log4.error(std::string("Retry count=") + std::to_string(retry::retries));

            if (retry::retries >= opt::max_retries) mailobj.send_email("Cap Monitor Giving Up", "capserver problem.", EMAIL_ADMIN);
            if (retry::retries==1) mailobj.send_email("Cap Monitor Entering Retry Loop", "capserver problem.", EMAIL_ADMIN);

            if (retry::retries >= opt::max_retries) throw; // rethrown
            log4.warn("Retry count below max retry value " + std::to_string(opt::max_retries) + ". Will continue.");
        } catch (JsonError e) {

            log4.error(std::string("JsonError exception thrown : ") + e.what());
            ++retry::retries;
            log4.error(std::string("Retry count=") + std::to_string(retry::retries));
            if (retry::retries >= opt::max_retries) mailobj.send_email("Cap Monitor Giving Up", "json problem.", EMAIL_ADMIN);
            if (retry::retries==1) mailobj.send_email("Cap Monitor Entering Retry Loop", "json problem.", EMAIL_ADMIN);
            if (retry::retries >= opt::max_retries) throw; // rethrown
            log4.warn("Retry count below max retry value " + std::to_string(opt::max_retries) + ". Will continue.");

        } catch (std::runtime_error e) {

            log4.error(std::string("Exception thrown : ") + e.what());
            throw; 

        }

        date = day();
        if ((iterations) && (iterations == cnt)) {
            log4.info("Final iteration");
            break;
        }

        if(retry::retries) {
            log4.warn(LAZY_SLEEP(RETRY_NAP));
            sleep(RETRY_NAP);
        } else sleep (nap);

    }

    log4.warn("Leaving loop");

}

void LOOP::pull_remote_file_list() {

    toolz::SQLITE_ADAPTOR::connect(opt::capdb);
    int local_max = toolz::SQLITE_ADAPTOR::get_instance()->generic_query<int>(SQL_GET_COUNT);
    // int sbmid_start = query_integer(sqlitedb, const_cast<char*>(SQL_GET_COUNT));

    if (local_max==0) log4.warn("We have an empty db");
    else if (local_max>0) {
        local_max = toolz::SQLITE_ADAPTOR::get_instance()->generic_query<int>(SQL_GET_MAX_SBM);
        std::string str("We have a db with " + std::to_string(local_max) + " entries and "); //?!?
        str += "max submission_id=" + std::to_string(local_max); 
        log4.info(str);
    } else throw std::runtime_error(THROW_DEFAULT("silly"));

    toolz::SQLITE_ADAPTOR::disconnect();

    std::stringstream thing;
    toolz::CURL_ADAPTOR::get_instance()->pull(stash::url + "type=annot:return=list:value=max", thing); 

    regex reg_max_sbm("\\{\"max\":\"(\\d+)\"\\}");
    smatch mo;
    
    if (!regex_search(thing.str(),mo,reg_max_sbm)) {
        mailobj.send_email("initial request response doesn't look like json.\r\n",thing.str().c_str(), EMAIL_ADMIN);
        throw runtime_error(THROW_DEFAULT("response doesn't look like json. are cookies aren't setup properly."));
    }

    int remote_max = atoi(std::string(mo[1]).c_str());

    SUBMISSIONLIST submissionlist;
    
    log4.debug("CURRENT LOCAL MAX = " + std::to_string(local_max));
    log4.debug("CURRENT REMOTE MAX = " + std::to_string(remote_max));

    regex reg_min(REGEX_MIN);
    regex reg_cookies(REGEX_ACCESS);

    // why?!
    for (int i=local_max?local_max:1,j=0 ; i-JSON_INCREMENT<remote_max&&i!=j; j=i,i=i+JSON_INCREMENT<remote_max?i+JSON_INCREMENT:remote_max) {
    
        if(j==0) continue;

        log4.debug("PULLING ENTRIES " + std::to_string(j==1?1:j+1) + "-" + std::to_string(i));

        stringstream ss;
        toolz::CURL_ADAPTOR::get_instance()->pull(stash::url + "type=annot:return=list:value=" + std::to_string(j==1?1:j+1) + "," + std::to_string(i), ss);

        std::string upload_json(ss.str()); 
        if (regex_match(upload_json,reg_cookies)) 
          throw runtime_error(THROW_DEFAULT("sub request response doesn't look like json."));

        if (regex_search(upload_json,reg_min)) log4.debug("Server looks healthy");
        else throw CapServerError(THROW_DEFAULT("Received unexpected response from server!?"));

        Json::Value root;
        Json::Reader reader;
        bool parsedSuccess = reader.parse(upload_json, root, false);
        
        if (!parsedSuccess)
          throw JsonError(THROW_DEFAULT("submission json wasn't as expected")); 
        
        int array_end = root.size();
        for ( int index = 0 ; index < array_end ; ++index ) { 

            //if (atoi(root[index]["sbm_id"].asCString()) <= local_max) continue; 
            const SUBMISSION rar(root[index]); 

            std::string str1("Cap server has new file '" + rar.file_name + "' (sbm_id=" + std::to_string(rar.sbm_id) 
              + ") from '" + rar.submitter_name + "' (" + rar.user_email + ") for species '" + rar.species + "'");
            log4.info(str1);
            std::string str2 = "Cap Monitor : New file uploaded to cap server for " + rar.species;

            if(!mailobj.send_email(str2.c_str(), str1.c_str(), EMAIL_ADMIN)) runtime_error("couldn't send email");
            submissionlist.push_back(std::move(rar)); 
        }
    }

    retry::retries = 0; // got proper response so reset retry count

    if(submissionlist.size()>0) { 
        
        toolz::SQLITE_ADAPTOR::connect(opt::capdb);
        std::for_each(submissionlist.begin(), submissionlist.end(), [&log4](SUBMISSION s) {
            log4.debug("Storing new entry sbmid=" + std::to_string(s.sbm_id) + " from " + s.user_email);
            toolz::SQLITE_ADAPTOR::get_instance()->generic_insert(s.to_sqlstring()); 
        });
        toolz::SQLITE_ADAPTOR::disconnect();
    }

    return;
}

void LOOP::download_remote_files() {

    char buf[STRING_LENGTH];
    sprintf(buf, SQL_GET_COUNT_BY_STATUS, FILE_NEW, FILE_MD5IGNORE);
    toolz::SQLITE_ADAPTOR::connect(opt::capdb);
    int sbmid = toolz::SQLITE_ADAPTOR::get_instance()->generic_query<int>(buf);

    if (sbmid==0) {
        log4.debug("There are no un-downloaded files.");
        return;
    }
        
    log4.info("There are " + std::to_string(sbmid) + " un-downloaded files.");

    sprintf(buf, SQL_GET_FILES, FILE_NEW, FILE_MD5IGNORE);

    auto files = toolz::SQLITE_ADAPTOR::get_instance()->generic_query<TABLE>(buf); 

    // why not generate directly or just iterate through the TABLE structure?!
    SUBMISSIONLIST submissionlist;
    for_each(files.begin(), files.end(), [&submissionlist](ROW r){ submissionlist.push_back(SUBMISSION(r)); });
    // submissionlist.push_back(sbm); // should be using std::move?!??!

    // seems i asked this before, why iterate through list again?!? - i.e. why not jsut process in first go?!
    for (SUBMISSIONLIST::iterator kit = submissionlist.begin(); kit != submissionlist.end(); kit++) { 

        std::string uri(stash::url + "/" + kit->file_name); 
        log4.debug("Downloading file " + uri);

        stringstream strstrm;
        toolz::CURL_ADAPTOR::get_instance()->pull(stash::url + "type=annot:return=file:value=" + kit->file_name, strstrm); 

        std::string file(strstrm.str()); 
        
        if(kit->file_size==0)  {
            std::string subject = "Cap : File Empty Error : " + kit->file_name;
            std::string str = std::string("There is a problem with file ") + kit->file_name + " (it's empty?!?) [sbm_id=" 
              + std::to_string(kit->sbm_id) + "/type=" + kit->file_type + "]";
            log4.error(str);
            log4.info("Will not store locally. Re-set status to FILE_NEW to re-attempt download");
            if (!mailobj.send_email(subject.c_str(), str.c_str(), EMAIL_ADMIN))
              throw runtime_error(THROW_DEFAULT("couldn't send email"));
            toolz::SQLITE_ADAPTOR::get_instance()->update_file_status(kit->sbm_id,FILE_EMPTY); 
            // update_file_status (sqlitedb,kit->sbm_id,FILE_EMPTY); 
            continue;
        }

        regex reg_body(".*?<body><pre>(.*)</pre></body></html>.*"); // regex reg_body(REGEX_HTML);
        smatch match_obj;

        // can't ignore non-asci as have nothing to work with - in all honestly there is little point in the md5 sums after this, but...
        if (!regex_match(file,match_obj,reg_body)){
            std::string subject = "Cap : File Format Error : " + kit->file_name;
            std::string str = std::string("There is a problem with file ") + kit->file_name + " (doesn't look like ASCI/download truncated?) [sbm_id=" 
              + std::to_string(kit->sbm_id) + "/type=" + kit->file_type + "]";
            if(shchk) std::cout << ANSI_REDB;
            log4.error(str);
            if(shchk) std::cout << ANSI_COFF;
            log4.info("Will not store locally. Re-set status to FILE_NEW to re-attempt download");

            if (!mailobj.send_email(subject.c_str(), str.c_str(), EMAIL_ADMIN))
              throw runtime_error(THROW_DEFAULT("couldn't send email"));

            toolz::SQLITE_ADAPTOR::get_instance()->update_file_status(kit->sbm_id,FILE_ASCI); 
            continue;
        }

        std::string file_contents(match_obj[1]);

        bool problem_flag = false;

        // switch (kit->sbm_status) {
        if(kit->sbm_status==FILE_MD5IGNORE) {
            if(shchk) std::cout << ANSI_REDB;
            log4.warn("Ignoring md5 checks.");
            if(shchk) std::cout << ANSI_COFF;
        } else {

            MD5((unsigned char*)file_contents.c_str(), file_contents.size(), md5_result); 

            if(kit->file_md5!=md5_sum2string(md5_result)) {

                problem_flag = true;

                if(shchk) std::cout << ANSI_REDB;
                std::string subject = "Cap : Md5 Error : md5 sums do not match : " + kit->file_name;
                log4.error(subject);
                if(shchk) std::cout << ANSI_COFF;

                std::string localmd5(md5_sum2string(md5_result));
                
                std::string str = std::string("There is a problem with file ") + kit->file_name + " (md5 sum doesn't match) [sbm_id="
                + std::to_string(kit->sbm_id) + "/type=" + kit->file_type + "]\n\n";
                    
                std::string s;
                if (kit->file_md5.empty() && kit->file_size==file_contents.size()) {
                    s = "However, the server-side md5 sub is empty and the file sizes are correct - try resetting status to FILE_LOCALCOPY";
                    log4.info(s);
                    str += "\n";
                } else if (kit->file_md5.empty()) {
                    s = "However, the server-side md5 sub is empty - try re-downloading by  resetting status to FILE_NEW";
                    log4.info(s);
                    str += "\n";
                } else str += "cap server md5 sum : '" + kit->file_md5 + "'\nlocal md5 sum : '" + localmd5 + "'\n"; 

                log4.info("Downloaded anyway. Reset to new to try again, or md5_ignore to just skip md5 checks.");

                toolz::SQLITE_ADAPTOR::get_instance()->update_file_status(kit->sbm_id,FILE_MD5); 
                if (!mailobj.send_email(subject.c_str(), str.c_str(), EMAIL_ADMIN))
                    throw runtime_error("couldn't send email");
            }

        }

        if(kit->file_type=="gff3"){ 

            stringstream gff_in(file_contents);
            stringstream slim_gff(stringstream::out);
            std::string sbuf;
            regex reg_ignore(TYPE_IGNORE);

            while(getline(gff_in,sbuf)) {
                if (regex_search(sbuf,reg_ignore)) continue;
                slim_gff << sbuf << "\n";
            }

            std::string slimfilename(kit->file_name);
            filedir(slimfilename,stash::sbmdir); // eek 
            slimfilename += LOCAL_GFF_SUFFIX;

            log4.info("storing as "+slimfilename);

            std::ofstream out(slimfilename.c_str());  
            if(out==0) throw runtime_error("problem opening gff file for writting"); 
            out << slim_gff.str();
            // out << slim_gff.rdbuf();
            out.close();

            {

            char* cmdtmp = (char*)malloc((MED_STRING_SIZE+slimfilename.size())*sizeof(char)); 
            assert(cmdtmp!=0);
            // scarry realisation that perl_run returns normally when invoked with -e if it can't open a file?!?
            snprintf(cmdtmp,MED_STRING_SIZE+slimfilename.size(),"sed -i 's/\\r//' %s",slimfilename.c_str());
            // snprintf(cmdtmp,MED_STRING_SIZE+slimfilename.size(),"perl -i -pe 's/\\r//' %s",slimfilename.c_str());
            assert(system(cmdtmp)==0);

            free(cmdtmp);
            
            }

            if (!problem_flag) toolz::SQLITE_ADAPTOR::get_instance()->update_file_status(kit->sbm_id,FILE_LOCALSTORE); 

        } else {
            log4.error("Illegal filetype - Cap system not implementing file type: " + std::string(kit->file_type));
            log4.info("To process as different filetype just update file_type and re-set status to FILE_NEW");
            toolz::SQLITE_ADAPTOR::get_instance()->update_file_status(kit->sbm_id,FILE_TYPE); 
            continue;
        }    
    }

    toolz::SQLITE_ADAPTOR::disconnect();

    return;
}

void LOOP::process_local_file_list() {
                     
    SUBMISSIONLIST submissionlist;

    toolz::SQLITE_ADAPTOR::connect(opt::capdb);
    toolz::SQLITE_ADAPTOR::get_instance()->get_local_file_list(submissionlist);
    toolz::SQLITE_ADAPTOR::disconnect();

    for (SUBMISSIONLIST::iterator kit = submissionlist.begin() ; kit != submissionlist.end() ; kit++ ) { 

        if(kit->file_type=="gff3") {

            std::string filename(kit->file_name);
            filepath(filename);
            filename = stash::sbmdir + "/" + filename + LOCAL_GFF_SUFFIX;

            if(!exists(filename)) throw runtime_error("fiel doesn't exist!?!");
            else log4.info("processing file " + filename);

            std::string dsn;
            std::string prefix;
            int max_id=0;
            toolz::SQLITE_ADAPTOR::connect(opt::capdb);
            bool has_prefix = toolz::SQLITE_ADAPTOR::get_instance()->species_info(kit->species, dsn, prefix, max_id);
            toolz::SQLITE_ADAPTOR::disconnect();

            if(has_prefix) {

                log4.debug(dsn); DB_PARAMS dbp(dsn); 

                if (check_cap_db(log4, dbp, kit, filename)) log4.debug("cap e! db looks fine");

                log4.info("running sbm_id=" + std::to_string(kit->sbm_id));

                if(!check_cap_db_no_rows(log4, dbp, const_cast<char*>("select count(1) from gene where cap_sbm_id = %d and cap_source = 'cap'"), kit->sbm_id)) {

                    log4.error("cap db " + std::string(dbp.dbname()) + " already has gene models with this sbmid!?!");
                    char *messageText = const_cast<char*>("Not safe to proceed");
                    union sigval signal_value;
                    signal_value.sival_ptr = messageText;
                    log4.info("Exiting Cap : Sending parent shut down signal");
                    sigqueue(getppid(), SIGUSR2, signal_value);
                    return;
                }

                // horrible, horrible, horrible...
                GFFPROC gffproc(
                  std::string(stash::bin + "/" + CMD_GFFDOC CMD_GFFDOC_OPTS).c_str(),  
                  std::string(stash::bin + "/" + CMD_EFINGER).c_str(), 
                  log4,
                  opt::capdb  /*, &mailobj*/
                ); 

                std::string output;
                std::string tmpmsg;
                int gff3val = gffproc.validate((char*)filename.c_str(), kit->sbm_id, &dbp, output, tmpmsg);

                std::string seqedit_report;
                if (gff3val==FILE_NONCDS) generate_seqedits(dbp, filename, kit->sbm_id, seqedit_report);

                std::string htmltable;
                max_id = get_html_table_for_genes(dbp, &mailobj, prefix, max_id, htmltable, kit->sbm_id); 

                {

                char buf[STRING_LENGTH];
                sprintf(buf, "update cap_species set max_id = %d where species = '%s'",max_id,kit->species.c_str());
                toolz::SQLITE_ADAPTOR::connect(opt::capdb);
                toolz::SQLITE_ADAPTOR::get_instance()->update(buf);
                toolz::SQLITE_ADAPTOR::get_instance()->update_file_status(kit->sbm_id,gff3val);
                toolz::SQLITE_ADAPTOR::disconnect();

                }
                
                std::string subject;
                stringstream ssm;
                
                // THEY GET ALL MESSAGES NOW
                bool email_user = true; // bool email_user = false; 

                switch (gff3val) {
                    case FILE_DONE:
                        subject =  kit->species + " : File Saved";
                        if(shchk) std::cout << ANSI_WHITEB;
                        log4.info(subject);
                        if(shchk) std::cout << ANSI_COFF;
                        output += htmltable;
                    break;
                    case FILE_NONCDS:
                        subject =  kit->species + " : File Saved has peptides with stop codons";
                        if(shchk) std::cout << ANSI_YELLOW;
                        log4.info(subject);
                        if(shchk) std::cout << ANSI_COFF;
                        output += htmltable + seqedit_report;
                    break;
                    case FILE_GFF3VAL1: 
                        subject =  kit->species + " : GFF3 Validation-1 Error";
                        if(shchk) std::cout << ANSI_REDB;
                        log4.error(subject);
                        if(shchk) std::cout << ANSI_COFF;
                        output = "Gff3 validation error at first validation step for file " + kit->file_name + "\n\n" + output;
                    break;
                    case FILE_GFF3VAL2:

                        subject =  "GFF3 Validation-2 Error (Experimental detailed report) : " + kit->species 
                          + " (" + kit->file_name + ").";

                        if(!mailobj.send_email(subject.c_str(), tmpmsg.c_str(), EMAIL_ADMIN))
                          runtime_error("couldn't send email");

                        subject =  kit->species + " : GFF3 Validation-2 Error";
                        if(shchk) std::cout << ANSI_FUNKYB;

                        log4.error(subject);

                        if(shchk) std::cout << ANSI_COFF;

                        output = "Gff3 validation error at second validation step (logged to file " + kit->file_name 
                          + SUFFIX_LOADING + ")." + output;

                    break;
                    case FILE_GFF3NONEWGENES:
                        subject =  kit->species + " : Gene Model Error (no new gene models detected)";
                        if(shchk) std::cout << ANSI_GREENB;
                        log4.error(subject);
                        if(shchk) std::cout << ANSI_COFF;
                        output = "Gff3 fingerprinting error - no new gene models detected (logged to file " + kit->file_name + SUFFIX_CLEANING + ").";
                        log4.error(output);
                    // default:
                    break;
                }

                subject += " : " + kit->file_name;

                output = "<p>Thankyou " + kit->user_email 
                  + ",</p><p>This is an automated email regarding the submission of file "
                  + kit->file_name + " to the Community Annotation Portal.\n\n</p>" 
                  + "<p>For further information see the <a href=\"http://vectorbase-cap.ensemblgenomes.org/?q=FAQ\">FAQ</a> or email " 
                  + EMAIL_LISTSERVE + "\n\n</p><h3>REPORT SUMMARY:</h3>\n\n" + output;

                if(!mailobj.send_email(subject.c_str(), output.c_str(), EMAIL_ADMIN))
                    runtime_error("couldn't send email");

                if(opt::realemail) {

                    std::for_each(opt::emailvector.begin(),opt::emailvector.end(),EMAILFNCTOR(mailobj,log4,subject,output));
                    // for_each(emailvector.begin(),emailvector.end(),pfunc); 

                    // get rid of this!?
                    if(email_user) {
                        std::string emailing = "Invoked with dev --realemail option and issue is relevant to user, emailing : ";
                        emailing += kit->user_email;
                        log4.info(emailing);
                        if(!mailobj.send_email(subject.c_str(), output.c_str(), kit->user_email.c_str()))
                          runtime_error("couldn't send email");
                    } else {
                        std::string emailing = "Invoked with dev --realemail option but issue not relevant to user - will not email.";
                        emailing += kit->user_email;
                        log4.warn(emailing);
                    }
                // #elif #undef...
                } else {
                    log4.warn("Invoked without --realemail so not emailing user");
                }
                // #endif

                #ifdef ONEFILE
                std::cout << "\nSilly=" << __FILE__ << ":" << __LINE__ << "\n" << std::endl;
                exit(0);
                #endif

            } else throw runtime_error (std::string("species unknown : ")+kit->species);

        } else throw runtime_error (THROW_DEFAULT("ilegal file format"));
    }

    #ifdef ONELOOP
    std::cout << "\nSilly=" << __FILE__ << ":" << __LINE__ << "\n" << std::endl;
    exit(0);
    #endif

    return;
}

std::string md5_sum2string(unsigned char* md) {
    int i;     
    char ca[2];
    std::string md5("");
    for(i=0; i <MD5_DIGEST_LENGTH; i++) {
        sprintf(ca,"%02x",md[i]);
        md5 += ca;
    }          
    return md5;
}    

void print_md5_sum(unsigned char* md) {
    int i;     
    for(i=0; i <MD5_DIGEST_LENGTH; i++) {
        printf("%02x",md[i]);
    }          
    printf("\n");
}    

size_t my_write_func(void *ptr, size_t size, size_t nmemb, FILE *stream) {
    std::cout << "here!!" <<std::endl;
    return fwrite(ptr, size, nmemb, stream);
} 

size_t my_write_func3(void *ptr, size_t size, size_t nmemb, FILE *stream) {

    // if forget to return it's called once!?!?!?
    int ret =  fwrite(ptr, size, nmemb, stream);
    std::cout<< "sizeof="<<sizeof((char*)ptr);
    std::string h((char*)ptr);
    std::cout << "STRING="<<h.size()<<std::endl;
    std::cout << "this is ret=" << ret<< " - this is=" << nmemb<< std::endl;

    return ret;

} 

bool check_cap_db_no_rows(log4cpp::Category&, DB_PARAMS & dbp, char* query, int value) {

    MYSQL *conn;
    conn = mysql_init(NULL);

    if(mysql_real_connect(conn,dbp.host(),dbp.user(),dbp.pass(),dbp.dbname(),dbp.port(),NULL,0) == NULL)
      throw MySqlError(THROW_DEFAULT("could not connect to species-specific e! cap database")); 

    std::string rs;
    char buf[STRING_LENGTH];
    sprintf(buf, query, value);

    if(!query_string(conn,buf,rs)) 
      throw runtime_error(THROW_CSTRING(query)); 

    mysql_close(conn);

    if(stoi(rs)==0) return true; //if(atoi(rs.c_str())==0) return true;

    return false;
}

void generate_seqedits(DB_PARAMS& dbp, std::string filename, int sbmid, std::string& report) {

    MYSQL *conn;
    conn = mysql_init(NULL); 

    if(mysql_real_connect(conn,dbp.host(),dbp.user(),dbp.pass(),dbp.dbname(),dbp.port(),NULL,0) == NULL)
      throw MySqlError(THROW_DEFAULT("could not connect to species-specific e! cap database")); 

    MYSQL_RES *result;
    MYSQL_ROW row;

    vector<string> update_strings;

    char q[STRING_LENGTH];
    sprintf(q,"SELECT gene_id FROM gene where biotype != 'protein_coding' and cap_sbm_id = %d order by gene_id DESC",sbmid);
    mysql_query(conn, q);
    result = mysql_store_result(conn);
    if (result==0) throw runtime_error("couldn't access gene table for id update");

    vector<int> non_coding_cds;
    while ((row = mysql_fetch_row(result))) non_coding_cds.push_back(atoi(row[0]));
    mysql_free_result(result);

    report = "<h3>SEQEDIT REPORT</h3>";
    vector<std::tuple<int,int,char>> transcripts_to_pull;

    for (auto it = non_coding_cds.begin() ; it != non_coding_cds.end() ; it++) {
        char x[250];
        sprintf(x,"select transcript_id,stable_id from transcript where gene_id = %d",*it);
        mysql_query(conn, x);
        
        result = mysql_store_result(conn);

        if (result==0) throw runtime_error("couldn't access gene table for id update");
        while ((row = mysql_fetch_row(result))) {
            std::cout << "NEED TO RUN SEQEDIT ROUTINE ON TRANSCRIPT "<< row[0] << " / " << row[1] << std::endl; 
            levenshtein_distance(dbp,filename,atoi(row[0]), row[1],report,conn);  // it's a temp. let's copy it?!?
        }
        mysql_free_result(result);
    }

    mysql_close(conn);

    return;

}

int get_html_table_for_genes(DB_PARAMS& dbp, mailer* mailptr, std::string prefix, int sqlite_gene_count, std::string& table, int sbmid) {

    MYSQL *conn;
    conn = mysql_init(NULL);

    if(mysql_real_connect(conn,dbp.host(),dbp.user(),dbp.pass(),dbp.dbname(),dbp.port(),NULL,0) == NULL)
      throw MySqlError(THROW_DEFAULT("could not connect to species-specific e! cap database")); // duh "'"

    MYSQL_RES *result;
    MYSQL_ROW row;

    if (prefix.empty()) {
        std::cout << "no organism prfix"<<std::endl; 
    } else {
        std::cout << "ORGANISM HAS A PREFIX VALUE ='"<<prefix<< "' THUS WILL GENERATE CAP IDS - CURRENT MAX="<<sqlite_gene_count<<std::endl; 

        vector<string> update_strings;

        // #define LIM " limit 20" // hhhrahrharehwajurhfsadfas, wtf?!
        char q[STRING_LENGTH];
        sprintf(q,"SELECT gene_id,stable_id  FROM gene where cap_sbm_id = %d",sbmid);
        mysql_query(conn, q);

        result = mysql_store_result(conn);
        if (result==0) throw runtime_error("couldn't access gene table for id update");

        vector<std::pair<int,int>> genes_to_pull;

        char sql[500];
        while ((row = mysql_fetch_row(result))) { 
            genes_to_pull.push_back(std::pair<int,int>(atoi(row[0]),sqlite_gene_count)); 
            sprintf(sql,"update gene set stable_id = '%s%06d', cap_orig_stable_id = '%s' where gene_id = %s",prefix.c_str(),sqlite_gene_count,row[1],row[0]);
            update_strings.push_back(sql);
            sqlite_gene_count++;
        }
        
        mysql_free_result(result);

        vector<std::tuple<int,int,char>> transcripts_to_pull;

        for (auto it = genes_to_pull.begin() ; it != genes_to_pull.end() ; it++) {
            char x[250];
            sprintf(x,"select transcript_id,stable_id from transcript where gene_id = %d",it->first);
            mysql_query(conn, x);
            result = mysql_store_result(conn);
            if (result==0) throw runtime_error("couldn't access gene table for id update");
            char suf = 'A';
            while ((row = mysql_fetch_row(result))) {
                int local_gene_count = std::get<1>(*it);
                transcripts_to_pull.push_back(std::tuple<int,int,char>(atoi(row[0]),local_gene_count,suf));// ->second,suf)); 
                sprintf(sql,"update transcript set stable_id = '%s%06d-R%c', cap_orig_stable_id = '%s' where transcript_id = %s",prefix.c_str(),local_gene_count,suf++,row[1],row[0]);
                update_strings.push_back(sql); 
            }
            mysql_free_result(result);
        }

        for (auto it = transcripts_to_pull.begin() ; it != transcripts_to_pull.end() ; it++) {
            char x[250];
            char trans_suf = std::get<2>(*it);
            sprintf(x,"select translation_id,stable_id from translation where transcript_id = %d",std::get<0>(*it)); // ->first);
            mysql_query(conn, x);
            result = mysql_store_result(conn);
            if (result==0) throw runtime_error("couldn't access gene table for id update");
            int u = 1;
            while ((row = mysql_fetch_row(result))) {
                string s(u==1?"":std::to_string(u));
                sprintf(sql,"update translation set stable_id = '%s%06d-P%c%s', cap_orig_stable_id = '%s' where translation_id = %s",
                prefix.c_str(),std::get<1>(*it),trans_suf,s.c_str(),row[1],row[0]);
                update_strings.push_back(sql);
            }
            mysql_free_result(result);
        }

        for (auto it = update_strings.begin() ; it != update_strings.end() ; it++ ) {
            mysql_query(conn,it->c_str());
            if(mysql_affected_rows(conn)==1) {
            } else if (mysql_affected_rows(conn)==0) {
                throw runtime_error("No rows were modified!?!"); 
            } else { 
                std::cout << "WTF="<<mysql_affected_rows(conn)<< std::endl; 
                throw runtime_error(THROW_CSTRING(it->c_str())); 
            }
        }
    }
    
    std::string rs;
    char buf[LONG_STRING_SIZE];
    snprintf(buf, LONG_STRING_SIZE, "select g.biotype, g.stable_id as cap_stable_id, g.cap_orig_stable_id as submitted_stable_id, "
      "g.gene_id as internal_gene_id, sr.name as scaffold_name, g.seq_region_start as gene_start, g.seq_region_end as gene_end " 
      "from gene g, seq_region sr where " 
      "sr.seq_region_id=g.seq_region_id and g.cap_sbm_id = %d order by g.gene_id", sbmid);

    if (mysql_query(conn,buf))
      throw runtime_error(THROW_DEFAULT("unable to check db"));

    MYSQL_FIELD *field;
    result = mysql_store_result(conn);

    int num_fields = mysql_num_fields(result);
    
    table = "<h3>SAVED GENE SUMMARY:</h3>"
    "<table style=\"width:585px;\"><tr>";

    // cos of that really strange incident when the server was hammered putting in a quite check?!
    bool all_unique=true;
    std::set<std::string> gene_names;

    while ((row = mysql_fetch_row(result))) {

        std::string bgc(std::string(*row)=="protein_coding"?"#e5f1ff":std::string(*row)=="terminal_stop_included"?"#FFB64E":"#FF927D"); 
        for(int i=0; i<num_fields; i++) {

            if (i==0) {
                while((field = mysql_fetch_field(result))) table += "<td BGCOLOR=\"#9FBC00\" style=\"padding:10px;\" >" + std::string(field->name) + "</td>";
                table += "</tr><tr>\n"; 
            }

            table += "<td BGCOLOR=\"" + bgc + "\" style=\"padding:10px;\" >";
            table += row[i] ? row[i] : "NULL";
            table += "</td>";

        }

        int name_index = prefix.empty() ? 1 : 2;
        if (gene_names.count(row[name_index])>0) all_unique=false;
        gene_names.insert(row[name_index]);

    }

    table += "</tr></table>";

    mysql_free_result(result);
    mysql_close(conn);

    if (!all_unique) 
      mailptr->send_email("VERY STRANGE : we have a non-unique id stored within a single submission\r\n","blah\r\n", EMAIL_ADMIN);

    return sqlite_gene_count;

}

bool check_cap_db(log4cpp::Category & log4, DB_PARAMS & dbp, SUBMISSIONLIST::iterator & kit, std::string & filename) {

    MYSQL *conn;
    conn = mysql_init(NULL);

    std::string str1("Connecting to e! cap database " + std::string(dbp.dbname()) + " for " + filename); // the rule of one-must-be...
    log4.debug(str1);

    if(mysql_real_connect(conn,dbp.host(),dbp.user(),dbp.pass(),dbp.dbname(),dbp.port(),NULL,0) == NULL)
      throw MySqlError(THROW_DEFAULT("could not connect to species-specific e! cap database")); // duh "'"

    std::string rs;
    if(!query_string(conn,SQL_GET_SPECIES,rs)) //throw runtime_error(rs);
      throw runtime_error(THROW_DEFAULT("unable to check e! cap meta table for meta_key 'cap.species_version'"));

    if (kit->species != rs) // if (kit->species != std::string(row[0])) 
      throw runtime_error(THROW_DEFAULT("cap db species doesn't match submission species"));

    if(!query_string(conn,"select meta_value from meta where meta_key = 'cap.base_uri'" ,rs)) //throw runtime_error(rs);
      throw runtime_error(THROW_DEFAULT("unable to check e! cap meta table for meta_key 'cap.base_uri'"));

    if (stash::url != rs) // if (kit->species != std::string(row[0])) 
      throw runtime_error(THROW_DEFAULT("cap db is locked to a different frontend (cap.base_uri)"));

    str1.clear(); // not sure why since just re-assining value?!?
    str1 = "e! cap db meta table species matches : " + rs;
    log4.debug(str1);

    if(!query_string(conn,SQL_CHECK_INIT,rs)) // throw runtime_error(rs);
      throw runtime_error(THROW_DEFAULT("unable to check crc64 table of e! cap db - has this been patched?"));

    if(stoi(rs)==0) { // if(atoi(rs.c_str())==0) {
        str1 += "The crc64 table is empty! - this database doesn't look initialised!";
        log4.warn(str1);
    } else {
        str1 += "Cap db has " + rs + " crc64 entries for canonical gene models";
        log4.debug(str1);
    }

    if(!query_string(conn,SQL_CHECK_READY,rs)) //throw runtime_error(rs);
      throw runtime_error(THROW_DEFAULT("unable to check e! cap gene table")); // make this a more specific exception type?!?

    if(stoi(rs)!=0) throw runtime_error(THROW_DEFAULT("There are unprocessed gene entries in e! cap db it isn't safe to continue!"));
    else log4.debug("Cap db is ready");

    mysql_close(conn);
    return true;

}

