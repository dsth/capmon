#include "boost/regex.hpp" 
#include <sstream>
#include "utils.h"
#include <openssl/md5.h>
#include <string>

#define SPC "<< \" \" <<"
#define DSN_PATTERN "The dsn string must be of the form: "

using std::string;

typedef boost::smatch::value_type subm;

class CAP_META {

    CAP_META ();
    std::string local_bindir_;
    std::string local_sbmdir_;
    std::string email_to_; 
    std::string email_from_; 
    std::string email_xtras_; 
    std::string capmon_loglevel_; 
    std::string remote_jsonurl_; 
    int capmon_restart_;
    int capmon_sleep_;
    bool capmon_xcurl_;

  public:

    CAP_META (char*);

    const char* local_bindir() { return local_bindir_.c_str(); } 
    const char* local_sbmdir() { return this->local_sbmdir_.c_str(); }
    const char* email_to() { return email_to_.c_str(); }
    const char* email_from() { return email_from_.c_str(); }
    const char* email_xtras() { return email_xtras_.c_str(); }
    const char* capmon_loglevel() { return capmon_loglevel_.c_str(); }
    const char* remote_jsonurl() { return remote_jsonurl_.c_str(); }
    unsigned int capmon_restart() { return capmon_restart_; }
    unsigned int capmon_sleep() { return capmon_sleep_; }
    unsigned int capmon_xcurl() { return capmon_xcurl_; }

};

std::ostream& operator<<(std::ostream & os, const subm & c) {
    // return os << "pink shit biyach!?!?!?!\n\n";
    std::string auto_tmp(c);
    return os << auto_tmp << " ";
}

DB_PARAMS::DB_PARAMS(std::string dsn) {

    boost::regex reg_dsn(DSN_REGEX);
    boost::smatch match_obj;
    if (boost::regex_match(dsn,match_obj,reg_dsn)){

        //// for some reason part assigning with temporaries just scambles stuff?!? - iterating?!?
        //// - try iterating through?!? - either pre-assign in batch to strings directly or use stringstream?!?

        std::stringstream pinky(std::stringstream::in|std::stringstream::out);
        pinky << match_obj[1] << match_obj[2] << match_obj[3] << match_obj[4] << match_obj[5];
        std::string db, h, u, p;
        int P;
        pinky >> db >> h >> P >> u >> p;

        /*
        // really shouldn't cast away const?!? - but it's just temporary so whatever
        dbname = (char*)db.c_str();
        host = (char*)h.c_str();
        */

        dbname_ = db;
        host_ = h;
        this->user_ = u;
        pass_ = p;
        port_ = P;

    } else throw std::runtime_error(DSN_PATTERN DSN_REGEX);
}

CAP_META::CAP_META(char*) { 
    std::cout << "let's grab those values"<<std::endl; 
    std::cout << "PINKERTON"<<std::endl; 
}

void filepath(string& filename) {
    static unsigned char md5_result[MD5_DIGEST_LENGTH];
    MD5((unsigned char*)filename.c_str(), filename.size(), md5_result); 
    char ca[3];
    string buf = std::move(filename);
    filename.clear();
    for(int i=0; i <2; i++) { 
        sprintf(ca,"%02X/",md5_result[i]);
        filename += ca;
    }                
    filename += buf;
    return;
}    
                                                             
