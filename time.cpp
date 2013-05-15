#include<stdio.h>   
#include<time.h>   

int * day () {
    static int sec;
    time_t t;
    time(&t);
    tm * other = gmtime(&t);
    sec = (*other).tm_mday;
    return &sec;
}

