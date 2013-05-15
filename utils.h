#ifndef UTILS_H
#define UTILS_H

#define RETRY_NAP 180
#define LAZY_SLEEP(x) "In retry loop, sleeping for " STRINGIFY(x) " seconds."
#define COOKIE "cookie.txt"
#define STRING_SIZE 100
#define LONG_STRING_SIZE 500
#define REALLY_LONG_STRING_SIZE 2500
#define MED_STRING_SIZE 750 // should put check on file name length?!?
#define COMMAND "loop"

#include <string>

#define DSN_REGEX "DBI:mysql:database=(\\w+?);host=([\\w\\.-]+?);port=(\\d+?),(\\w+?),(\\w+?)"

class DB_PARAMS { // struct DB_PARAMS {

    std::string host_;
    std::string user_; 
    std::string pass_; 
    std::string dbname_;
    unsigned int port_;
    DB_PARAMS (); // prevent default construction - why?

public:

    DB_PARAMS (std::string _h, std::string _u, std::string _pw, int _p, std::string _db)
      : host_(_h), user_(_u), pass_(_pw), dbname_(_db), port_(_p) {}

    DB_PARAMS (std::string); // just pass by value - again why?

    const char * dbname() { return dbname_.c_str(); } 
    const char * host() { return this->host_.c_str(); }
    const char * user() { return user_.c_str(); }
    const char * pass() { return pass_.c_str(); }
    unsigned int port() { return port_; }

};

void filepath(std::string&);

#endif
