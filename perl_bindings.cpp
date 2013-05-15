#include <EXTERN.h>
#include <perl.h>
#include "perl_bindings.h"

using namespace std;

EXTERN_C void xs_init (pTHX);
EXTERN_C void boot_DynaLoader (pTHX_ CV* cv);
EXTERN_C void xs_init(pTHX) {
	const char *file = __FILE__;
	dXSUB_SYS;
	newXS("DynaLoader::boot_DynaLoader", boot_DynaLoader, file);
}

static PerlInterpreter *my_perl;

std::string perl_grab_sequence_via_api(
    const char *host, 
    const char *user, 
    const char *pass, 
    const char *dbname, 
    int port,
    int transcid,
    bool translate
  ) {

    char *my_argv[] = { (char*)"", (char*)"perl_bindings.pl" }; // const_cast!?!
    my_perl = perl_alloc();
    perl_construct(my_perl);

    perl_parse(my_perl, xs_init, 2, my_argv, NULL);

    int retval;

    SV* value;

    HV * arghash = newHV();

    hv_store(arghash, "-dbname", 7, sv_2mortal(newSVpv(dbname,0)),0);
    hv_store(arghash, "-user", 5, sv_2mortal(newSVpv(user,0)),0);
    hv_store(arghash, "-pass", 5, sv_2mortal(newSVpv(pass,0)),0);
    hv_store(arghash, "-host", 5, sv_2mortal(newSVpv(host,0)),0);
    hv_store(arghash, "-port", 5, sv_2mortal(newSViv(port)),0);
    char trans_id[] = "transcript_id";
    hv_store(arghash, trans_id, strlen(trans_id), sv_2mortal(newSViv(transcid)),0); // hash keys are not null terminated?!?

    // assign to hashref?!?
    SV* hashref = newRV_inc((SV*)arghash);

    dSP;
    ENTER;
    SAVETMPS;
    PUSHMARK(sp);

    // push it on to perl stack?!?
    XPUSHs(sv_2mortal(hashref));
    PUTBACK;

    // feeling really lazy so ...
    std::string what(translate?"translatedCDS":"splicedExons");
    retval = perl_call_pv(what.c_str(), G_ARRAY);
    // retval = perl_call_pv("reorder", G_ARRAY);
    SPAGAIN;
    std::string pink;

    if (retval == 1) {
        value = POPs; 
        pink = (char*)SvPV_nolen(value);
    } else { 
        exit(1);
    }

    PUTBACK;
    FREETMPS;
    LEAVE;

    perl_destruct(my_perl);
    perl_free(my_perl);
    return pink;

}

