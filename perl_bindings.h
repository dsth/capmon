#ifndef PERL_GRAB_SEQS
#define PERL_GRAB_SEQS

#include <string>

extern std::string perl_grab_sequence_via_api(
    const char*,
    const char*,
    const char*,
    const char*,
    int,
    int,
    bool translate
);

#endif
