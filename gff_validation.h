#ifndef GFF_VALIDATION_H
#define GFF_VALIDATION_H
#include <string>
#include <map>
#include "feature.h"
#include "utils.h"

#define PERMITTED_TRANSCRIPT_TYPES_REGEX    "mRNA|tRNA|pseudogenic_tRNA|rRNA|miRNA|ncRNA"
enum PERMITTED_BIOTYPES                    { mRNA,tRNA,pseudogenic_tRNA,rRNA,miRNA,ncRNA };

struct feature_ext : public feature_min {

    // explicit feature_ext(const feature_min& fm) : feature_min(fm) {}
    feature_ext() = delete;
    feature_ext(const feature_min& fm) = delete;
    ~feature_ext () {}

    feature_ext(std::string chr, std::string p, uint a, uint b, unsigned char s, PERMITTED_BIOTYPES x) :
      feature_min(a,b), 
      _scfname(chr),
      _parent(p), 
      _strand(s), 
      _biotype(x) {} 
    
    std::string parent() const { return _parent; }   // don't ever return a non-const ref on a private member - i.e. you've just abolished privateness?!?
    unsigned char strand() const { return _strand; }
    PERMITTED_BIOTYPES biotype() const { return _biotype; }
    std::string scfname() const { return _scfname; }

private :
    std::string _scfname;
    std::string _parent;
    unsigned char _strand;
    PERMITTED_BIOTYPES _biotype;
};

struct gff_holder {
    //y overlap, fragmentation and names checks at..
    std::map<std::string,feature_ext> mrna_coords;
    std::multimap<std::string,feature_min> exonbymrna_coords;
    //y check cds have exons and/or check all mRNAs have CDS and exons
    std::multimap<std::string,feature_min> cdsbymrna_coords;
};

bool capmon_gff_validation (const char*, std::string&, DB_PARAMS* dbp=0);
bool validation_tests (std::string& report);
bool standalone_gff_validation (const char*, DB_PARAMS*);
bool parse_gff_for_loading (const char* filename, gff_holder*);

#endif
