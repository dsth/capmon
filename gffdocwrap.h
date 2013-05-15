#ifndef GFFPROC_H
#define GFFPROC_H

#include "utils.h"
#include "log4cpp/Category.hh"
#include <string>

class GFFPROC {
    
    static const char summary_template[];
    std::string efingercmd;
    std::string gffdoccmd;
    log4cpp::Category& log4;
    DB_PARAMS* dbp;
    char* capdb; 

public:

    GFFPROC (std::string _ss, std::string _s, log4cpp::Category& _l, char* _cdb /*, mailer* _m*/) : 
      efingercmd(_s), gffdoccmd(_ss), log4(_l), capdb(_cdb)/*, mailptr(_m) */ {}
    int validate (char*, int, DB_PARAMS*, std::string&, std::string&);
    bool gffcheck (const char*, std::string&);
    int gffload (const char*, std::string&);
    std::pair<int,int> eclean (const char*, int);

};

#endif
