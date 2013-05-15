#include <iostream>
#include "config.h"
#include <getopt.h>
#include "batch.h"
#include "boost/lexical_cast.hpp"
#include "boost/filesystem.hpp"
#include "boost/regex.hpp"

#define HELP_LIST "\nLists submitted files as of last check of file server.\n\nOptions to return partial list:\n\n\
        --help      -h      This little message\n\
        --sbmid     -s      Submission id\n\
        --status    -p      Submission status\n\
        --species   -o      Filter by species\n\
        --user      -u      File server username\n\
        --submitter -n      Submitter name (wildcarded)\n\
        --email     -e      email address associated with username (wildcarded)\n"

#define HELP_RESET "\nReset the status of an individual submission. Requires -sbmid and -status arguments.\n"

#define STRING_LENGTH 250

using std::runtime_error;

namespace opt {
    static char * user;    
    static char * submitter; 
    static char * email; 
    static int sbmid = 0; 
    static char * status; 
    static char * capdb; 
    static char * species; 
}

enum { OPT_CAPDB=1 };

static struct option long_options[] = {
    {"help",       no_argument,  0,  'h'},
    {"user",       required_argument,  0,  'u'},
    {"sbmid",      required_argument,  0,  's'},
    {"species",      required_argument,0,  'o'},
    {"capdb",      required_argument,0,  OPT_CAPDB},
    {"email",      required_argument,  0,  'e'},
    {"submitter",  required_argument,  0,  'n'},
    {"status",      required_argument, 0,  'p'},
    {0, 0, 0, 0}
};

void parse_batch_opts(int argc, char **argv, bool) {

    for (int c; (c = getopt_long (argc, argv, "ho:u:s:n:e:p:",long_options, NULL));) {

        if (c == -1) break; // unsafe way - can overflow

        switch (c) {
            case 'u':           opt::user = optarg;             break;
            case 'o':           opt::species = optarg;          break;
            case OPT_CAPDB:     opt::capdb = optarg;            break;
            case 'e':           opt::email = optarg;            break;
            case 'n':           opt::submitter = optarg;        break;
            //y want shortopt form too so need to put in initialiser here rather than have it specified directly in long_options with required_arg, &opt::sbmid, 0}...
            case 's':           opt::sbmid = atoi(optarg);      break;
            case 'p':           opt::status = optarg;           break;
            case 'h':  
                std::cout << HELP_LIST << std::endl;
                exit(0);
            break; // why?
            default:    abort ();
            }
        }

    if(opt::capdb == NULL)  opt::capdb  =     const_cast<char*>(SQLITE_DB_NAME);
       
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
sbm_status, \
submitter_name, user, user_email, \
file_md5, \
file_size, \
file_desc, \
file_name"

#define SEP "--------------------------"

void Batch(int ac, char ** av) {

    parse_batch_opts(ac,av,true);

    std::cout << "this is in planning stages!?!?" << std::endl; 
    exit(0);

}





