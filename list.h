#ifndef LIST_H
#define LIST_H

void List (int, char **);
void Expt (int, char**);
void Pull (int, char**);
void Reset (int, char**);
void Ready(int, char**);
void Species (int, char **) __attribute__((noreturn));
void Config(int, char**) __attribute__((noreturn));
void Setup (int, char**);
void Validate (int ac, char** av);
void Rebind(int, char**);

#endif
