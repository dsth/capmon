#include <signal.h>
#include "list.h"
#include "loop.h"
#include "batch.h"
#include "gffdocwrap.h"
#include "toolz.h"

#define DUMMY_SBM_ID 0

#define SQL_QUERY_STATUS "select status from cap_status where status_id \
  = (select max(status_id) from cap_status)"
#define SQL_CHECK_DB "select count(*) from cap_status"
#define SQL_INSERT_STATUS_READY "insert into cap_status (status) values ('ready')"
#define SQL_INSERT_STATUS_LOOPING "insert into cap_status (status) values ('looping')"

#define CAP_STATUS_READY 1
#define CAP_STATUS_LOOPING 2 
#define CAP_STATUS_LOCKED 3

#define SQL_CAP_STATUS "insert into cap_status (status) values (%d)"

#define EXEC "capmon"
#define VERSION "0.0.0"

// DROP TABLE IF EXISTS `source_rank`; 
// alter table gene drop column cap_sbm_id;
// alter table gene drop column cap_orig_stable_id;

void terminate (int sig, siginfo_t *info, void *ptr)  __attribute__((noreturn));

toolz::SQLITE_ADAPTOR* toolz::SQLITE_ADAPTOR::me_ptr = 0; // seems std::nullptr not yet implemented?!?  
toolz::MYSQL_ADAPTOR* toolz::MYSQL_ADAPTOR::me_ptr = 0;
toolz::CURL_ADAPTOR* toolz::CURL_ADAPTOR::me_ptr = 0;

using std::runtime_error;

const static char * CORRECT_USAGE = 
"Exec: " EXEC
"\nVersion: " VERSION " "
"\nBuild: " __DATE__ " @ " __TIME__
"\nWith gcc: " __VERSION__
"\nWith glibc: "  TOSTRING(__GLIBC__) "."  TOSTRING(__GLIBC_MINOR__)
"\n"
"\nDaniel S. T. Hughes. dsth@ebi.ac.uk, dsth@cantab.net\n"
"\nUsage: ./" EXEC " [COMMAND]\n"
"\n\tconfig\t\tbasic configuration admin"
"\n\tlist\t\tprint formatted file status list to stdout"
"\n\tloop\t\tstart monitor"
"\n\tpatch\t\tpatch a core e! db to a patch db"
"\n\tsetup\t\tsetup local cap file db instance"
"\n\tspecies\t\tcore e! cap database admin"
"\n\tready\t\treset cap db if loop exits improperly"
"\n\treset\t\treset the status of an individual file submission";

// this is all awful?!
static int my_pid;
void terminate (int, siginfo_t*, void*) { 
    if (my_pid==getpid()) 
      std::cout << "\n---\nSomeone sent a SIGINT\n---\n"; 
    exit(0);
}

int main(int ac, char **av) {

    my_pid = getpid();

    /*
    void (*prev_fn)(int);
    prev_fn = signal (SIGTERM,terminate); 
    prev_fn = signal (SIGINT,terminate); 
    */

    struct sigaction exit_act;
    memset (&exit_act, '\0', sizeof(exit_act));  
    exit_act.sa_sigaction = terminate; 
    exit_act.sa_flags = SA_SIGINFO;    
    sigaction(SIGTERM,&exit_act,0);
    sigaction(SIGINT,&exit_act,0);

    if (ac <= 1) {
        std::cout << CORRECT_USAGE << std::endl;
        return 0;
    }

    std::string mode(*(av+1));
    std::cout <<"Mode="<< mode << std::endl;

    if (mode=="single") {
    } else if (mode=="gff3") {
    } else if (mode=="list") {
        List(ac-1,av+1);
        return 0;
    } else if (mode=="validate") {
        Validate(ac-1,av+1);
        return 0;
    } else if (mode=="expt") {
        Expt(ac-1,av+1);
        return 0;
    } else if (mode=="pull") {
        Pull(ac-1,av+1);
        return 0;
    } else if (mode=="species") {
        Species(ac-1,av+1);
        return 0;
    } else if (mode=="config") {
        Config(ac-1,av+1);
        return 0;
    } else if (mode=="setup") {
        Setup(ac-1,av+1);
        return 0;
    } else if (mode=="reset") {
        Reset(ac-1,av+1);
        return 0;
    } else if (mode=="patch") {
        // std::cout << PATCH_30OCT2012<<"\n";
        std::cout << "not today...\n";
        return 0;
    } else if (mode=="unpatch") {
        // std::cout << UNPATCH_30OCT2012<<"\n";
        std::cout << "not any more\n";
        return 0;
    } else if (mode=="ready") {
        Ready(ac-1,av+1);
        return 0;
    } else if (mode=="loop") {
        Loop(ac-1,av+1);
    } else if (mode=="bload") {
        Batch(ac-1,av+1);
    } else if (mode=="rebind") {
        Rebind(ac-1,av+1);
    } else {
        std::cout << CORRECT_USAGE << std::endl;
        return 0;
    }


}

