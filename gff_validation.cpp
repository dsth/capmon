#include <fstream>
#include <iostream>
#include <sstream>
#include <cassert>
#include "boost/regex.hpp"
#include <set>
#include <algorithm>
#include "exceptions.h"
#include <bitset>
#include <iomanip> 
#include "utils.h"
#include "gff_validation.h"
#include <vector>

#ifdef CAPMON_EXT
#define NSUBMISSION
#include "toolz.h"
#endif
#define STAMPY(x) std::make_pair(#x,x)

#define IGNORE_TYPES_REGEX                  "contig|supercontig|match|match_part"
#define CORE_PERMITTED                      "gene|pseudogene|exon|CDS|"
#define ALL_PERMITTED                       CORE_PERMITTED PERMITTED_TRANSCRIPT_TYPES_REGEX
#define MAX_BAD_LINES 5
#define SPRINTF_STRING      "select count(1) from seq_region where name = '%s'"
#define SQL_STRING          "select * from seq_region where name = ''"
#define APOLLO_SCF_NAMES                            (1ull<<0) // (0x1<<0) 0x1
#define NAMES_HAVE_SPACES                           (1ull<<1)
#define LINES_WO_9COLS                              (1ull<<2)
#define LINE_ENDINGS                                (1ull<<3)
#define EXCESS_GENE_CONSISTENCY_PROB                (1ull<<5)
#define EXCESS_TRANS_REL_GENE_CONSISTENCY_PROB      (1ull<<6)
#define EXCESS_TRANS_REL_CDS_EXON_CONSISTENCY_PROB  (1ull<<7)
#define EXCESS_CDS_EXON_CONSISTENCY_PROB            (1ull<<8)
#define NON_PERMITTED_BIOTYPES                      (1ull<<9)
#define ID_WITHOUT_PARENT_NOT_GENE_PSEUDOGENE       (1ull<<10)
#define PARENT_WITHOUT_ID_NOT_CDS_EXON              (1ull<<11)

#define UNKNOWN_SCAFFOLD                            (1ull<<12)

#define NEGATIVE_COORDINATES                        (1ull<<13)
#define NON_PRINTING_X0D                            (1ull<<14)
#define BLANK_LINES                                 (1ull<<15)

#define NON_UNIQUE_ID                               (1ull<<16)

#define PSEUDOGENES_PRESENT                         (1ull<<17)
#define NO_CDS                                      (1ull<<18)

#define EMBL_FORMAT                                 (1ull<<19)
#define GFF_FASTA_HEADER                            (1ull<<20)
#define FASTA_HEADER                                (1ull<<21)

#define PARTIAL_MODEL                               (1ull<<22)

#define NO_FEATLINES                                (1ull<<23)
#define NO_GENES                                    (1ull<<24)
#define NO_EXON_CDS                                 (1ull<<25)
#define NO_TRANSCRIPTS                              (1ull<<26)

#define TRANSCRIPT_LACKS_EXONS                      (1ull<<27)
#define PROTEIN_CODING_LACKS_CDS                    (1ull<<28)
#define OVERLAPPING_EXONS                           (1ull<<29)
#define END_LESS_THAN_START                         (1ull<<30)

#define CDS_OUT_OF_TRANSCRIPT_BOUNDS                (1ull<<38)
#define CDS_MISSING_EXON                            (1ull<<39)

#define STAMP(x) std::make_pair(x,#x)

bool check_capmon_bool_test (const char*, std::string&, DB_PARAMS* dbp=0);

static void split (std::vector<std::string>& vs, std::string& blah, const char* token) {
    char* copy = strdup(blah.c_str());
    char *p = strtok((char*)copy,token);
    while (p) {
        vs.push_back(p);
        p = strtok(NULL,token);
    }
    free(copy);
}

using std::cout;
using std::endl;
using std::runtime_error;

typedef unsigned int uint;
typedef unsigned char uchar; 
typedef boost::regex regex;
typedef boost::smatch smatch;
typedef feature_min feature;
typedef unsigned long long unsignedll;
typedef std::map<unsigned long long,std::string> PEP;
typedef std::pair<std::string,feature_min> featpair;

unsignedll check_raw_bitflag_consistency_tests (const char*, std::string&, DB_PARAMS* dbp = 0);

struct name_holder {
    std::set<std::string> scaffolds;
    std::set<std::string> gene_ids;
    std::set<std::string> parents;
    std::set<std::string> gene_transcript_ids; // there's really not much reason for gene_transcript_ids to exist!?!?
};

class print_row {

    unsigned long long bf;
    unsigned long long bm;
    bool blah;
    std::string issue;
    std::string item;

  public:

	print_row(std::string i, std::string a, unsigned long long _bf, unsigned long long _bm) : 
      bf(_bf), 
      bm(_bm),
      issue(a), 
      item(i) {}

	print_row(std::string i, std::string a, bool b) :
      blah(b),
      issue(a), 
      item(i) {}

	friend std::ostream& operator<< (std::ostream& os, const print_row& bobj) {

        bool problem = (bobj.bf&bobj.bm); 

        std::string td, s; 
        if (bobj.bm==0) {
            td = bobj.blah ? "#FF927D" : "#e5f1ff";
            s =  bobj.blah ? "TRUE":"";
        } else { 
            td = problem ? "#FF927D" : "#e5f1ff"; 
            s =  problem ? "TRUE" : "FALSE"; 
        }

        td="<td BGCOLOR=\"" + td + "\" style=\"padding:10px;\" >";

		return os << " <tr> " << td 
          << bobj.item << " </td> " << td << bobj.issue 
          << " </td>" << td 
          << s << "</td></tr>\n"; 

	}

};

static std::map<std::string,PERMITTED_BIOTYPES> biotype_resolver = { 

/*
    std::make_pair("mRNA",mRNA),
    std::make_pair("pseudogenic_tRNA",pseudogenic_tRNA),
    std::make_pair("rRNA",rRNA),
    std::make_pair("miRNA",miRNA),
    std::make_pair("ncRNA",ncRNA)
    // std::make_pair("pseudogene",pseudogene) 
*/

    STAMPY(mRNA),
    STAMPY(pseudogenic_tRNA),
    STAMPY(rRNA),
    STAMPY(miRNA),
    STAMPY(ncRNA),

};

class print_em_pretty {
  
  public: 
    static PEP _pep;
    std::string grab_name(unsigned long long bf) {
        if(_pep.count(bf)==0) throw runtime_error("seems someone didn't populate the map?!?");
        _pep.find(bf);
        return _pep.find(bf)->second;
    }

};

PEP print_em_pretty::_pep = { 
    STAMP(CDS_MISSING_EXON), 
    STAMP(CDS_OUT_OF_TRANSCRIPT_BOUNDS),
    STAMP(END_LESS_THAN_START),
    STAMP(OVERLAPPING_EXONS),
    STAMP(PROTEIN_CODING_LACKS_CDS),
    STAMP(TRANSCRIPT_LACKS_EXONS),
    STAMP(NO_TRANSCRIPTS),
    STAMP(NO_EXON_CDS),
    STAMP(NO_GENES),
    STAMP(NO_FEATLINES),
    STAMP(PARTIAL_MODEL),
    STAMP(FASTA_HEADER),
    STAMP(GFF_FASTA_HEADER),
    STAMP(EMBL_FORMAT),
    STAMP(UNKNOWN_SCAFFOLD),
    STAMP(NO_CDS),
    STAMP(PSEUDOGENES_PRESENT),
    STAMP(NON_UNIQUE_ID),
    STAMP(APOLLO_SCF_NAMES),
    STAMP(NAMES_HAVE_SPACES),
    STAMP(LINES_WO_9COLS),
    STAMP(LINE_ENDINGS),
    STAMP(EXCESS_GENE_CONSISTENCY_PROB),
    STAMP(EXCESS_TRANS_REL_GENE_CONSISTENCY_PROB),
    STAMP(EXCESS_TRANS_REL_CDS_EXON_CONSISTENCY_PROB),
    STAMP(EXCESS_CDS_EXON_CONSISTENCY_PROB),
    STAMP(NON_PERMITTED_BIOTYPES),
    STAMP(ID_WITHOUT_PARENT_NOT_GENE_PSEUDOGENE),
    STAMP(PARENT_WITHOUT_ID_NOT_CDS_EXON),
    STAMP(UNKNOWN_SCAFFOLD),
    STAMP(NEGATIVE_COORDINATES),
    STAMP(NON_PRINTING_X0D),
    STAMP(BLANK_LINES),
    STAMP(NON_UNIQUE_ID)
};

static unsigned long long PROBLEM = ( EXCESS_GENE_CONSISTENCY_PROB | EXCESS_TRANS_REL_GENE_CONSISTENCY_PROB 
    | EXCESS_TRANS_REL_CDS_EXON_CONSISTENCY_PROB | EXCESS_CDS_EXON_CONSISTENCY_PROB | ID_WITHOUT_PARENT_NOT_GENE_PSEUDOGENE 
    | PARENT_WITHOUT_ID_NOT_CDS_EXON | UNKNOWN_SCAFFOLD|NON_PERMITTED_BIOTYPES | NEGATIVE_COORDINATES 
    | NON_UNIQUE_ID | EMBL_FORMAT | GFF_FASTA_HEADER | FASTA_HEADER | PARTIAL_MODEL
    | NO_FEATLINES | NO_GENES | NO_EXON_CDS | NO_TRANSCRIPTS | TRANSCRIPT_LACKS_EXONS 
    | OVERLAPPING_EXONS | LINES_WO_9COLS | PROTEIN_CODING_LACKS_CDS | END_LESS_THAN_START 
    | CDS_OUT_OF_TRANSCRIPT_BOUNDS | CDS_MISSING_EXON );

// PROBLEM = ~PROBLEM;
// static unsigned long long NO_CDS_AND_NO_PSEUDOGENES = ( CDS_PRESENT | PSEUDOGENES_PRESENT );

namespace gff {

    static regex reg_embl       ("FT\\s.*");
    static regex reg_gfffasta   ("#{1,2}FASTA.*");
    static regex reg_sequence   ("[actgnACTGN]+");
    static regex reg_fasta      (">\\w+.*");
    static regex reg_comment        ("\\s*#.*");
    static regex reg_line_ends  (".*[^;]\\s*");
    static regex reg_blank      ("\\s*");
    static regex reg_carriageret         ("\\r");

    static regex reg_feat       ("([^\\t]+)\\t[^\\t]+\\t([^\\t]+)\\t\\d+\\t\\d+\\t[^\\t]+\\t[+-\\.]+\\t[^\\t]+\\t([^\\t]+)");
    static regex reg_negstart   ("([^\\t]+)\\t[^\\t]+\\t([^\\t]+)\\t-\\d+\\t\\d+\\t[^\\t]+\\t[^\\t]+\\t[^\\t]+\\t([^\\t]+)");
    static regex reg_negend     ("([^\\t]+)\\t[^\\t]+\\t([^\\t]+)\\t-?\\d+\\t-\\d+\\t[^\\t]+\\t[^\\t]+\\t[^\\t]+\\t([^\\t]+)"); 

    static regex reg_scfname    ("(\\w+?\\d+?):\\d+?[:-]\\d+?");
    static regex reg_id         ("ID=([^;]+)");
    static regex reg_parent     ("Parent=([^;]+)");

    static regex reg_ignore_types(IGNORE_TYPES_REGEX);
    static regex reg_allow_types(ALL_PERMITTED);
    static regex reg_transcript_biotypes(PERMITTED_TRANSCRIPT_TYPES_REGEX);
}

namespace details {

static unsigned long long gff_basic_validation_1a_gff_parse (const char* filename, std::stringstream& strstrm, name_holder& nh, gff_holder& gh) {

    unsigned long long bitflag = ( NO_FEATLINES|NO_GENES|NO_CDS|NO_EXON_CDS|NO_TRANSCRIPTS );

    smatch match_obj;
    smatch match_obj2;

    std::ifstream in(filename);

    if(in == 0) throw runtime_error("problem opening gff file : " + std::string(filename));

    std::string s;
    unsigned char error_counter = 0;
    unsigned long long ln = 0;

    while(getline(in, s)) { 

        ++ln;

        if(regex_match(s,gff::reg_blank)) {
          bitflag |= BLANK_LINES;
          continue;
        } 

        if(regex_match(s,gff::reg_line_ends)) bitflag |= LINE_ENDINGS; 

        if(regex_match(s,gff::reg_gfffasta)) {
            bitflag |= GFF_FASTA_HEADER;
            strstrm << "<p>All lines must be valid gff feature lines or comments. Consequently, we do not accept submissions containing fasta sequence!</p>\n";
            continue;
        } else if(regex_match(s,gff::reg_fasta)) {
          bitflag |= FASTA_HEADER;
          continue;
        } else if(regex_match(s,gff::reg_embl)) {
          bitflag|=EMBL_FORMAT;
          continue;
        }
       
        // odd apollo thing of putting in negative coords 
        if (regex_match(s,gff::reg_negstart) || regex_match(s,gff::reg_negend)) {
            bitflag |= NEGATIVE_COORDINATES; 
            continue;
        }

        // do last or mis-diagnose the annoyance above
        if(regex_match(s,gff::reg_comment)) continue;

        if (!regex_match(s,match_obj,gff::reg_feat)) { 
            bitflag |= LINES_WO_9COLS;

            // be more helpful
            if(error_counter<MAX_BAD_LINES) { 
              strstrm << "<p>line " << ln << " does not have gff 9-col, tab-delimited format ";

                int tabs = count(s.begin(),s.end(),'\t');
                int spaces = count(s.begin(),s.end(),' ');

                if(tabs==8 && spaces < 7) strstrm << "- start & end must be numeric and strand must conform to '+', '-' or '.'";
                else if (tabs > 5 && tabs < 11) strstrm << "- please check all 9 columns are present";
                else if (spaces > 7 && spaces < 9 && tabs < 3) strstrm << "- is this line delimited with spaces and not tabs?";
                else if (regex_match(s,gff::reg_sequence)) strstrm << "- this line looks like sequence.";
                else strstrm << "- i'm really not sure what this is?";

                strstrm << " :</p>\n<pre>    " << s << "</pre>\n";

            } else if (error_counter==MAX_BAD_LINES) 
              strstrm<<"<p>Non 9-column, tab-delimited format errors truncated.</p>\n";

            error_counter++;
            continue;
        }

        /// 2 : we have a feature. at this stage we prolly ought to check that strand is overtly +/-?!?

        bitflag &= ~NO_FEATLINES;
        std::stringstream featurestrm(s); // lazy but we know it adheres to correct format 
        std::string scfname, score, strand, source, type, ignore, annot;
        int start = 0, end = 0; 

        featurestrm >> scfname >> source >> type >> start >> end >> score >> strand >> ignore;

        if(end<start) {
            strstrm << "<p>Seriously - end less than start?!?</p>\n";
            strstrm << "<pre>    " << s << "</pre>\n";
            bitflag |= END_LESS_THAN_START; // allow for 1 base cds
        }

        annot = match_obj[3]; 
        type = match_obj[2];

        if (regex_match(type,gff::reg_ignore_types)) continue; // ignore contigs etc.

        if (!regex_match(type,gff::reg_allow_types)) { // block other types?!?
            strstrm << "<p>Ilegal feature type='" << type << "'</p>\n";
            bitflag |= NON_PERMITTED_BIOTYPES;
        }
        
        // horrible apollo ids with coords inserted - temporarily clean them up?!?
        if (!scfname.empty() && regex_match(scfname,match_obj,gff::reg_scfname)) {
            nh.scaffolds.insert(match_obj[1]);
            bitflag|=APOLLO_SCF_NAMES; 
        } else if (!scfname.empty()) nh.scaffolds.insert(scfname); // completely unecessary conditional

        /// 2a : we have a 2o or 3o feature (have id and parent) : transcript or exon/CDS features

        if (regex_search(annot,match_obj,gff::reg_id)&&regex_search(annot,match_obj2,gff::reg_parent)){ 

            std::string id = match_obj[1];
            std::string parent = match_obj2[1];

            if(id.find(' ')!=std::string::npos || parent.find(' ')!=std::string::npos) bitflag|=NAMES_HAVE_SPACES;
                
            // really shouldn't do this without a reason - i.e. needs to be exon/cds and have ',' in parent id?!?
            std::vector<std::string> parents;

            if(type == "exon" || type == "CDS") split(parents,parent,",");
            else parents.push_back(parent);

            for (auto i = parents.begin(); i!=parents.end() ; i++) {

                if(type == "CDS") { 
                    bitflag &= ~NO_CDS; 
                    gh.cdsbymrna_coords.insert(std::pair<std::string,feature_min>(*i,feature_min(start,end))); 
                } else if(type == "exon") gh.exonbymrna_coords.insert(std::pair<std::string,feature_min>(*i,feature_min(start,end))); 

                if (type == "exon" || type == "CDS") { 
                    bitflag &= ~NO_EXON_CDS;
                    nh.parents.insert(*i);
                } else if (regex_match(type,gff::reg_transcript_biotypes)) {

                    bitflag &= ~NO_TRANSCRIPTS;
                    if(nh.gene_transcript_ids.count(id)!=0) {
                        bitflag |= NON_UNIQUE_ID;
                        strstrm << "<p>ID tags must be unique : i've seen the ID '"<<id<<"' before (at transcript/gene level no less!)</p>\n";
                    }

                    assert(biotype_resolver.find(type)!=biotype_resolver.end());
                    gh.mrna_coords.insert(std::pair<std::string,feature_ext>(id,feature_ext(scfname,*i,start,end,strand=="+"?1:0,biotype_resolver.find(type)->second))); 
                    nh.gene_transcript_ids.insert(id);

                } else {
                    bitflag |= NON_PERMITTED_BIOTYPES;
                    strstrm << "<p>Not accepting biotype '" << type << "' at transcript level</p>\n";
                    if(type.find(' ')!=std::string::npos) strstrm << "<pre>    We do not tolerate spaces in biotype names</pre>\n";
                }

            }

        /// 2b : we have 1o features - only have id

        } else if (regex_search(annot,match_obj,gff::reg_id)) {   

            std::string id = match_obj[1];

            if(id.find(' ')!=std::string::npos) bitflag|=NAMES_HAVE_SPACES;

            if(type == "pseudogene") bitflag|=PSEUDOGENES_PRESENT;

            if (type == "gene" || type == "pseudogene") {

                bitflag &= ~NO_GENES;
                if(nh.gene_transcript_ids.count(id)!=0) {
                    bitflag |= NON_UNIQUE_ID;
                    strstrm << "<p>ID tags must be unique : i've seen the ID '"<<id<<"' before (at transcript/gene level no less!)</p>\n";
                }

                nh.gene_ids.insert(id);
                nh.gene_transcript_ids.insert(id);

            } else {
                bitflag |= ID_WITHOUT_PARENT_NOT_GENE_PSEUDOGENE;
                strstrm << "<p>The following is not a gene/pseudogene and thus must have a Parent tag: </p>\n<pre>    "<<s<<"</pre>\n";
            }

        /// 2c : we have 3o feautres - only parent

        } else if (regex_search(annot,match_obj2,gff::reg_parent)) {

            std::string parent = match_obj2[1];
            if(parent.find(' ')!=std::string::npos) bitflag|=NAMES_HAVE_SPACES;

            std::vector<std::string> parents;

            if(type == "exon" || type == "CDS") split(parents,parent,",");
            else parents.push_back(parent);

            for (auto i = parents.begin(); i!=parents.end() ; i++) {

                if(type == "CDS") { // storing these now...
                    bitflag &= ~NO_CDS;
                    gh.cdsbymrna_coords.insert(std::pair<std::string,feature_min>(*i,feature_min(start,end))); 
                } else if(type == "exon") gh.exonbymrna_coords.insert(std::pair<std::string,feature_min>(*i,feature_min(start,end))); 

                ///y is this a very naughty feature?!?
                if (type == "exon" || type == "CDS") {
                    bitflag &= ~NO_EXON_CDS;
                    nh.parents.insert(*i);
                    //nh.parents.insert(parent);
                } else {
                    bitflag |= PARENT_WITHOUT_ID_NOT_CDS_EXON;
                    strstrm << "<p>The following is not a CDS/exon and thus must have an ID tag: </p>\n<pre>    "<<s<<"</pre>\n";
                }
            }

        } else { 
            // throw runtime_error("you seem to have broken file parsing quite severly?!? - is this something to do with merging of submission.h/tool.z and json.h - separate?!?"); 
            std::cout << "\nARGH!?!? check the regex!?!\n";
            sleep(5);
        }

    } 

    return bitflag;

}

static unsigned long long gff_basic_validation_1b_gff_name_checks (unsigned long long bitflag, std::stringstream& strstrm, name_holder& nh, gff_holder& gh,DB_PARAMS* dbp=0) {


    std::set<std::string> transcript_parents;
    for_each(gh.mrna_coords.begin(),gh.mrna_coords.end(),[&transcript_parents](std::pair<std::string,feature_ext> e){ transcript_parents.insert(e.second.parent()); });

    std::vector<std::string> top_level_problems;
    set_difference(nh.gene_ids.begin(),nh.gene_ids.end(),transcript_parents.begin(),transcript_parents.end(), std::back_inserter(top_level_problems));
    if (top_level_problems.size() > 0) {
        strstrm <<"<p>Top level set difference from genes_ids with no transcript_parents (missing transcript/excess gene features) = " <<top_level_problems.size()<<"</p>\n";
        bitflag |= EXCESS_GENE_CONSISTENCY_PROB;
    }
    for_each(top_level_problems.begin(),top_level_problems.end(),[&strstrm](std::string s) { strstrm << "<p> * Excess 'gene' ID=" << s << "</p>\n"; });
    top_level_problems.clear();

    set_difference(transcript_parents.begin(),transcript_parents.end(),nh.gene_ids.begin(),nh.gene_ids.end(),std::back_inserter(top_level_problems));
    if (top_level_problems.size() > 0) {
        strstrm << "<p>Top level set difference from transcript_parents with no gene_ids (excess transcript/missing gene features) = " <<top_level_problems.size() << "</p>\n";
        // "with no gene_ids (excess transcript/missing gene features) = " <<top_level_problems.size() << "</p>\n";
        bitflag |= EXCESS_TRANS_REL_GENE_CONSISTENCY_PROB;
    }
    for_each(top_level_problems.begin(),top_level_problems.end(),[&strstrm](std::string s) { strstrm << "<p> * No 'gene' corresponding to transcript Parent=" << s << "</p>\n"; });
    top_level_problems.clear();

    std::set<std::string> transcript_ids;
    for_each(gh.mrna_coords.begin(),gh.mrna_coords.end(),[&transcript_ids](std::pair<std::string,feature_ext> e){ transcript_ids.insert(e.first); });

    set_difference(transcript_ids.begin(),transcript_ids.end(),nh.parents.begin(),nh.parents.end(),std::back_inserter(top_level_problems));
    // std::inserter(top_level_problems, top_level_problems.end()) - was using set
    if (top_level_problems.size() > 0) {
        strstrm << "<p>Top level set difference from transcript_ids with no cds/exon_parents  (missing cds-exon/excess transcript features) = "
          << top_level_problems.size() << "</p>\n";
        bitflag |= EXCESS_TRANS_REL_CDS_EXON_CONSISTENCY_PROB;
    }
    for_each(top_level_problems.begin(),top_level_problems.end(),[&strstrm](std::string s) { strstrm << "<p> * Excess 'transcript' ID=" << s << "</p>\n"; } );
    top_level_problems.clear();

    set_difference(nh.parents.begin(),nh.parents.end(),transcript_ids.begin(),transcript_ids.end(),std::back_inserter(top_level_problems));
    // std::inserter(top_level_problems, top_level_problems.end()) - was using set
    if (top_level_problems.size() > 0) {
        strstrm << "<p>Top level set difference from cds/exon parents with no transcript_id (excess cds-exon/missing transcript features) = "<<top_level_problems.size() << "</p>\n";
        bitflag |= EXCESS_CDS_EXON_CONSISTENCY_PROB; 
    }
    for_each(top_level_problems.begin(),top_level_problems.end(),[&strstrm](std::string s) { strstrm << "<p> * No 'transcript' corresponding to cds/exon Parent=" << s << "</p>\n"; });

    if (dbp) { 
#ifdef CAPMON_EXT
        toolz::MYSQL_ADAPTOR::connect(*dbp);

        for (std::set<std::string>::iterator sit = nh.scaffolds.begin() ; sit != nh.scaffolds.end() ; sit++) { // prolly ought to be a functor...
            char qbuf[STRING_SIZE]; 
            sprintf(qbuf, SPRINTF_STRING, sit->c_str());
            if (toolz::MYSQL_ADAPTOR::get_instance()->generic_query<int>(qbuf) == 0) {
            // if (std::stoi(count) == 0) {
                strstrm << "<p>There is no scaffold named " << sit->c_str() << " in cap db</p>\n";
                bitflag |= UNKNOWN_SCAFFOLD;
            }
        }

        toolz::MYSQL_ADAPTOR::disconnect();
#else
        cout << "\n[1c] MUST compile with CAPMON_EXT for scaffold name checks aagainst mysql/e! db instance\n";
#endif
    }
    return bitflag;

}

static unsigned long long gff_basic_validation_1d_individual_model_checks (unsigned long long bitflag, std::stringstream& strstrm, gff_holder& gh) {

    for (auto it = gh.mrna_coords.begin(); it!=gh.mrna_coords.end(); ++it) {

        int count = gh.exonbymrna_coords.count(it->first);

        if (it->second.biotype()==mRNA && gh.cdsbymrna_coords.count(it->first)==0) {
            bitflag |= PROTEIN_CODING_LACKS_CDS;
            strstrm << "<p>Protein-coding transcript "<< it->first << " lacks CDS features (is it a pseudogene?).</p>\n";
        }

        switch (count) {

            case 0: 
                bitflag |= TRANSCRIPT_LACKS_EXONS;
                strstrm << "<p>Transcript "<< it->first << " has no exons!</p>\n";
            break;

            case 1: { // single exon case...
                // std::_Rb_tree_iterator<std::pair<const std::string, feature_min> > exon_it = exonbymrna_coords.find(it->first);
                std::multimap<std::string,feature_min>::iterator exon_it = gh.exonbymrna_coords.find(it->first);
                if ((it->second).gstart()!=(exon_it->second).gstart()||(it->second).gend()!=(exon_it->second).gend()) {
                    std::cout << "[single exon] FRAGMENTED MODEL "<<it->first << std::endl; 
                    strstrm << "<p>Transcript "<< it->first << " is fragmented (exon features do not fully account for transcript extension) :</p>\n"
                      <<"<pre>   transcript: "<<(it->second).gstart()<<"-"<<(it->second).gend()<<"\n    exon: "<<(exon_it->second).gstart()<<"-"<<(exon_it->second).gend()<<"</pre>\n";
                    bitflag |= PARTIAL_MODEL;
                }
            break; }

            default: /* declaring vars in switch!?! */ { 

                auto mmit = gh.exonbymrna_coords.equal_range(it->first);
                // auto buffer = mmit;
                unsigned int min_start = mmit.first->second.gstart();
                unsigned int max_end = mmit.first->second.gend();
                bool break_out = false;

                for (auto exon_it = ++mmit.first ; exon_it != mmit.second; ++exon_it) { 
                // for (auto exon_it = ++mmit.first ; exon_it != mmit.second && !break_out; ++exon_it) { 

                    // simple fragmentation due to positional extraction
                    if((exon_it->second).gstart()<min_start) min_start=(exon_it->second).gstart();
                    if((exon_it->second).gend()>max_end) max_end=(exon_it->second).gend();

                    // exon overlap
                    for (auto exon_it2 = mmit.first ; exon_it2 != mmit.second ; ++exon_it2) { 

                        if(exon_it==exon_it2) continue;
                        else if((exon_it2->second).gstart()<=(exon_it->second).gend()&&(exon_it2->second).gend()>=(exon_it->second).gstart()) {
                            bitflag |= OVERLAPPING_EXONS;
                            strstrm << "<p>exon " << exon_it->first << " (" << (exon_it->second).gstart() << "-" << (exon_it->second).gend()
                              << ") and " << exon_it2->first << " (" << (exon_it2->second).gstart() << "-" << (exon_it2->second).gend() 
                              << ") of transcript " << it->first << " overlap!</p>\n<p>Bypassing further overlap checks</p>\n"; 
                            break_out = true; 
                            break;
                        }
                        if (break_out) break;
                    }
                }

                if ((it->second).gstart()!=min_start||(it->second).gend()!=max_end) {
                    std::cout << "[multi exon] FRAGMENTED MODEL " << it->first << std::endl; 
                    strstrm << "<p>Transcript "<< it->first << " is fragmented (exon features do not fully account for transcript extension) :</p>\n"
                      <<"<pre>    transcript: "<<(it->second).gstart()<<"-"<<(it->second).gend()<<"\n    exons: "<<min_start<<"-"<<max_end<<"</pre>\n";
                    bitflag |= PARTIAL_MODEL;
                }

                // really wanted to avoid this - but it's just too much hassle not to - put a short-circuit to only run these tests with higher validation?!?

                auto cdsit = gh.cdsbymrna_coords.equal_range(it->first);
                std::vector<featpair> sillycds(cdsit.first, cdsit.second);

                // we could order by strand but not really much point?!?
                sort(sillycds.begin(),sillycds.end(),[](featpair a, featpair b) { return (a.second.gstart()<b.second.gstart()); });
                unsigned int c = 0;
                unsigned int lim = sillycds.size()-1;

                std::stringstream no_exon_strm;

                for(auto cit = sillycds.begin(); cit!=sillycds.end(); cit++,c++) {

                    if(cit->second.gstart()<it->second.gstart()||cit->second.gend()>it->second.gend()) {
                        strstrm << "<p>Transcript " << it->first << " has CDS that does not fall within bounds of transcript :</p>\n" 
                          << "<pre>    transcript: "<<(it->second).gstart()<<"-"<<(it->second).gend()<<"\n    CDS: " << cit->second.gstart() << "-" << cit->second.gend() << "</pre>\n";
                        bitflag |= CDS_OUT_OF_TRANSCRIPT_BOUNDS;

                    }

                    // seriously - we need to actually find the paired exons for them, it makes me want to cry?!?
                    std::stringstream silly_exon;
                    std::string argh;

                    bool no_exon = true;
                    auto blah = mmit.first;
                    for (auto exon_it = --blah ; exon_it != mmit.second; ++exon_it) { 

                        if(sillycds.size()==1&&cit->second.gstart()>=it->second.gstart()&&cit->second.gend()<=it->second.gend()) no_exon = false;
                     
                        if(c==0) {
                            argh = "lefthand flanking";
                            if(cit->second.gstart()>=exon_it->second.gstart()&&cit->second.gend()==exon_it->second.gend()) no_exon = false;
                        } else if(c>0&&c<lim) {
                            argh = "internal";
                            if(cit->second.gstart()==exon_it->second.gstart()&&cit->second.gend()==exon_it->second.gend()) no_exon = false;
                        } else if(c==lim) {
                            argh = "righthand flanking";
                            if(cit->second.gstart()==exon_it->second.gstart()&&cit->second.gend()<=exon_it->second.gend()) no_exon = false;
                        }

                        if(cit->second.gstart()<=exon_it->second.gend()&&cit->second.gend()>=exon_it->second.gstart()) 
                          silly_exon << "<pre>        > Candidate exon spans " << exon_it->second.gstart() << "-" << exon_it->second.gend() << "</pre>\n";

                    }

                    // std::string argh = c==0||c==lim ? "flanking" : "internal";
                    if (no_exon) {
                        /// should probably have a fragmneted gene model check by cds - but in lieu of that this should do the trick?!?
                        no_exon_strm << "<pre>    Cannot find exon corresponding to " << argh << " CDS with range " << cit->second.gstart() << "-" << cit->second.gend() 
                          << " for transcript " << it->first << "</pre>\n";
                        no_exon_strm << silly_exon.str();
                        bitflag |= CDS_MISSING_EXON;
                    }
                }

                if ((bitflag&CDS_MISSING_EXON)==CDS_MISSING_EXON) {

                    if (!no_exon_strm.str().empty()) 
                       strstrm << "<p>Problem locating matching exon for CDS. <b>Note :</b> All CDS must be paired with an exon as to establish the limits of any UTRs. "
                      "The flanking (first & last) CDS of a transcript must have an internal boundary that corresponds to its exon pair (the external boundary must "
                      "within the bounds of the exon), "
                      "while, internal CDS must match both the start and end of their corresponding exons! <i>This problem is often due to "
                      "modification of a CDS without concomitant modification of its paired exon</i> : </p>\n";
                    strstrm << no_exon_strm.str();
                }

            break; }
        };

    }

    return bitflag;

}

static void gff_basic_validation_1x_file_cleanup (unsigned long long bitflag,const char* filename,std::stringstream& strstrm) {

    // why?!?
    if ((bitflag&LINE_ENDINGS)==LINE_ENDINGS) {
        char cmdtmp[MED_STRING_SIZE];
        sprintf(cmdtmp,"perl -i -pe 's/(.*[^;]\\s*)\\n/${1};\\n/' %s",filename);
        // strstrm << "> fixing line endings in-place" << endl;
        if (system(cmdtmp) != 0) 
          throw runtime_error ("problem fixing line endings!");
    }

    if ((bitflag&APOLLO_SCF_NAMES)==APOLLO_SCF_NAMES) { 
        char cmdtmp[MED_STRING_SIZE];
        sprintf(cmdtmp,"perl -i -pe 's/^(\\w+\\d+):1-\\d+\\t/${1}\\t/' %s",filename);
        strstrm << "<p>Fixing apollo scaffold name endings in-place</p>\n";
        if (system(cmdtmp) != 0) 
          throw runtime_error("problem fixing apollo scaffold name endings!");
    }

    //r make sure we fix this last - i.e. after ';' etc..
    if ((bitflag&NAMES_HAVE_SPACES)==NAMES_HAVE_SPACES) { 

        char cmdtmp[MED_STRING_SIZE];
        sprintf(cmdtmp,"perl -i -pe 'if(/^(.*ID=)([ \\S]+?)(;.*)$/) { $c=$1;$d=$2;$e=$3;$d=~s/ /_/g; $_= qq{$c$d$e\\n} }' %s",filename);
        strstrm << "<p>Removing spaces from IDs in-place</p>\n";
        if (system(cmdtmp) != 0) 
          throw runtime_error("problem removing spaces from IDs!");
        sprintf(cmdtmp,"perl -i -pe 'if(/^(.*Parent=)([ \\S]+?)(;.*)$/) { $c=$1;$d=$2;$e=$3;$d=~s/ /_/g; $_= qq{$c$d$e\\n} }' %s",filename);
        strstrm << "<p>Removing spaces from Parents in-place</p>\n";
        if (system(cmdtmp) != 0) 
          throw runtime_error("problem removing spaces from Parents!");
    }

}

std::string capmon_html_table(unsigned long long bitflag,std::stringstream& strstrm) {

    std::stringstream strstrm2(std::stringstream::out);

    strstrm2 << "<p></p><table style=\"width:500px;\">\n"
      << print_row("A01","Contains EMBL qualifiers",bitflag,EMBL_FORMAT)
      << print_row("A02","Contains fasta entries",bitflag,(GFF_FASTA_HEADER|FASTA_HEADER))
      << print_row("A03","Does not contain gene features",bitflag,NO_GENES) 
      << print_row("A04","Does not contain transcript features",bitflag,NO_TRANSCRIPTS)
      << print_row("A05","Does not contain exon/CDS features",bitflag,NO_EXON_CDS)
      << print_row("A06","There are non-comment, non-9-col format lines present",bitflag,LINES_WO_9COLS)
      << print_row("A07","Gene id excess vs. transcript parent",bitflag,EXCESS_GENE_CONSISTENCY_PROB)
      << print_row("A08","Transcript parent excess vs. gene id",bitflag,EXCESS_TRANS_REL_GENE_CONSISTENCY_PROB)
      << print_row("A09","Transcript id excess vs. exon/CDS parent",bitflag,EXCESS_TRANS_REL_CDS_EXON_CONSISTENCY_PROB)
      << print_row("A10","Exon/CDS parent excess vs. transcript id",bitflag,EXCESS_CDS_EXON_CONSISTENCY_PROB)
      << print_row("A11","Unknown scaffold names",bitflag,UNKNOWN_SCAFFOLD)
      << print_row("A12","Non-permitted biotypes",bitflag,NON_PERMITTED_BIOTYPES)
      << print_row("A13","Negative coordinates - known Apollo bug",bitflag,NEGATIVE_COORDINATES)
      << print_row("A15","Non-gene features with ID but no Parent ",bitflag,ID_WITHOUT_PARENT_NOT_GENE_PSEUDOGENE)
      << print_row("A16","Non-CDS/exon features with Parent but no ID",bitflag,PARENT_WITHOUT_ID_NOT_CDS_EXON)
      << print_row("A17","Broken (partial) transcript(s)",bitflag,PARTIAL_MODEL)
      << print_row("A18","Non-unique ID at gene/transcript level",bitflag,NON_UNIQUE_ID)
      << print_row("A19","No gff feature lines",bitflag,NO_FEATLINES)
      << print_row("A20","Transcript lacks exons",bitflag,TRANSCRIPT_LACKS_EXONS)
      << print_row("A21","Gene has overlapping exons",bitflag,OVERLAPPING_EXONS)
      << print_row("A22","protein coding genes lacks CDS features",bitflag,PROTEIN_CODING_LACKS_CDS)
      << print_row("A23","feature end cannot be lower than feature start",bitflag,END_LESS_THAN_START)
      << print_row("A24","CDS does not have a corresponding exon",bitflag,CDS_MISSING_EXON)
      << print_row("A25","CDS falls outside of transcript bounds",bitflag,CDS_OUT_OF_TRANSCRIPT_BOUNDS)
      << "</table><p></p>";

    if (!strstrm.str().empty()) {
        std::string blah = "\n\n<h3>DETAILED REPORT:</h3>\n" + strstrm.str();
        strstrm2 << blah;
    }
    std::string report = strstrm2.str();

    return std::move(report);

}

} // end of namespace?!?

bool validation_tests (std::string& report) {

#ifdef CAPMON_EXT
{
    DB_PARAMS dbpnew("mysql-eg-devel-3.ebi.ac.uk","ensrw","scr1b3d3",4208,"dsth_CapDb_Mar07_anopheles_gambiae_core_13_66_3");
    assert( check_raw_bitflag_consistency_tests("UNKNOWN_SCAFFOLD",report, &dbpnew) == UNKNOWN_SCAFFOLD);
    assert( check_raw_bitflag_consistency_tests("UNKNOWN_SCAFFOLD_2KNOWNNAME",report, &dbpnew) == 0);
    assert ( check_capmon_bool_test ("UNKNOWN_SCAFFOLD",report, &dbpnew) == false );
    assert ( check_capmon_bool_test ("UNKNOWN_SCAFFOLD_2KNOWNNAME",report, &dbpnew) == true );

}

//    { 
//    unsigned long long bitflag = check_raw_bitflag_consistency_tests("CDS_MISSING_EXON",report);
//    cout << "return value =\n\t" << std::bitset<sizeof(long long)*8>(bitflag) << "\n";
//    for (unsigned long long i = 0 ; i < sizeof(long long)*8 ; i++) if (bitflag&(1ull<<i)) cout << " Active bit from return " << std::dec << i << "\n"; //  << " and " << x << "\n";
//    cout << report <<"\n";
//    }

    assert ( check_capmon_bool_test ("FINE",report) == true );
    assert ( check_capmon_bool_test ("NON_PERMITTED_BIOTYPES",report) == false );
    assert ( check_capmon_bool_test ("NAMES_HAVE_SPACES",report) == true );
    assert ( check_capmon_bool_test ("LINE_ENDINGS",report) == true );
    //     assert ( check_capmon_bool_test ("NON_PRINTING_X0D",report) == true );
    assert ( check_capmon_bool_test ("APOLLO_SCF_NAMES",report) == true );
    assert ( check_capmon_bool_test ("LINES_WO_9COLS",report) == false );
    assert ( check_capmon_bool_test ("EXCESS_GENE_CONSISTENCY_PROB",report) == false );
    assert ( check_capmon_bool_test ("EXCESS_TRANS_REL_GENE_CONSISTENCY_PROB",report) == false );
    assert ( check_capmon_bool_test ("EXCESS_TRANS_REL_CDS_EXON_CONSISTENCY_PROB",report) == false );
    //    assert ( check_capmon_bool_test ("EXCESS_CDS_EXON_CONSISTENCY_PROB",report) == false );
    assert ( check_capmon_bool_test ("NON_PERMITTED_BIOTYPES",report) == false );
    assert ( check_capmon_bool_test ("OVERLAPPING_EXONS",report) == false );
    assert ( check_capmon_bool_test ("ID_WITHOUT_PARENT_NOT_GENE_PSEUDOGENE",report) == false );
    assert ( check_capmon_bool_test ("NON_UNIQUE_ID",report) == false );

    assert ( check_capmon_bool_test ("PARENT_WITHOUT_ID_NOT_CDS_EXON",report) == false );
    assert ( check_capmon_bool_test ("NEGATIVE_COORDINATES",report) == false );
    assert ( check_capmon_bool_test ("BLANK_LINES",report) == true );
    assert ( check_capmon_bool_test ("PSEUDOGENES_PRESENT",report) == true );
    // assert ( check_capmon_bool_test ("NO_CDS",report) == true ); // not really relevant anymore?!?
    assert ( check_capmon_bool_test ("EMBL_FORMAT",report) == false );
    assert ( check_capmon_bool_test ("GFF_FASTA_HEADER",report) == false );
    assert ( check_capmon_bool_test ("FASTA_HEADER",report) == false );
    assert ( check_capmon_bool_test ("PARTIAL_MODEL",report) == false );
    cout << "here = "<<report<<"\n";
    assert ( check_capmon_bool_test ("NO_FEATLINES",report) == false );
    assert ( check_capmon_bool_test ("NO_GENES",report) == false );
    assert ( check_capmon_bool_test ("NO_EXON_CDS",report) == false );
    assert ( check_capmon_bool_test ("NO_TRANSCRIPTS",report) == false );
    assert ( check_capmon_bool_test ("TRANSCRIPT_LACKS_EXONS",report) == false );
    assert ( check_capmon_bool_test ("PROTEIN_CODING_LACKS_CDS",report) == false );
    assert ( check_capmon_bool_test ("OVERLAPPING_EXONS",report) == false );

    // these will modify the files
    if(system("cp testfiles/NAMES_HAVE_SPACES.gff_ORIG testfiles/NAMES_HAVE_SPACES.gff")!=0) throw runtime_error("couldn't copy file!?!");
    if(system("cp testfiles/APOLLO_SCF_NAMES.gff_ORIG testfiles/APOLLO_SCF_NAMES.gff")!=0) throw runtime_error("couldn't copy file!?!");
    if(system("cp testfiles/LINE_ENDINGS.gff_ORIG testfiles/LINE_ENDINGS.gff")!=0) throw runtime_error("couldn't copy file!?!");
    if(system("cp testfiles/NON_PRINTING_X0D.gff_ORIG testfiles/NON_PRINTING_X0D.gff")!=0) throw runtime_error("couldn't copy file!?!");
    if(system("cp testfiles/GFF_FASTA_HEADER.gff_ORIG testfiles/GFF_FASTA_HEADER.gff")!=0) throw runtime_error("couldn't copy file!?!");
    if(system("cp testfiles/FASTA_HEADER.gff_ORIG testfiles/FASTA_HEADER.gff")!=0) throw runtime_error("couldn't copy file!?!");
    if(system("cp testfiles/EMBL_FORMAT.gff_ORIG testfiles/EMBL_FORMAT.gff")!=0) throw runtime_error("couldn't copy file!?!");

#endif

    //  bitflag = check_raw_bitflag_consistency_tests("GFF_FASTA_HEADER",report); cout << "return value =\n\t" << std::bitset<sizeof(long long)*8>(bitflag) << "\n"; 
    //    for (int i = 0 ; i < sizeof(long long)*8 ; i++) if (int x = bitflag&(1ull<<i)) cout << " Active bit from return " << std::dec << i << "\n"; //  << " and " << x << "\n";
    //  unsigned long long bitflag = check_raw_bitflag_consistency_tests("LINES_WO_9COLS",report);
    //  cout << "return value =\n\t" << std::bitset<sizeof(long long)*8>(bitflag) << "\n";
    //    for (unsigned int i = 0 ; i < sizeof(long long)*8 ; i++) if (bitflag&(1ull<<i)) cout << " Active bit from return " << std::dec << i << "\n"; //  << " and " << x << "\n";
    //  cout << "report=\n\t"<<report<<"\n";

    assert( check_raw_bitflag_consistency_tests("FINE",report) == 0);
    assert( check_raw_bitflag_consistency_tests("CDS_MISSING_EXON",report) == CDS_MISSING_EXON );
    assert( check_raw_bitflag_consistency_tests("CDS_OUT_OF_TRANSCRIPT_BOUNDS",report) == (CDS_OUT_OF_TRANSCRIPT_BOUNDS|CDS_MISSING_EXON) );
    assert( check_raw_bitflag_consistency_tests("NON_PERMITTED_BIOTYPES",report) == (NON_PERMITTED_BIOTYPES|PARENT_WITHOUT_ID_NOT_CDS_EXON) );
    assert( check_raw_bitflag_consistency_tests("NON_PERMITTED_BIOTYPES",report) != (PARENT_WITHOUT_ID_NOT_CDS_EXON) );
    assert( check_raw_bitflag_consistency_tests("LINES_WO_9COLS",report) == LINES_WO_9COLS ); // 8 col
    assert( check_raw_bitflag_consistency_tests("LINES_WO_9COLS_2STRAND",report) == LINES_WO_9COLS); // 8 col
    assert( check_raw_bitflag_consistency_tests("LINES_WO_9COLS_3SEQ",report) == (LINES_WO_9COLS|LINE_ENDINGS) ); // 8 col
    assert( check_raw_bitflag_consistency_tests("OVERLAPPING_EXONS",report) == OVERLAPPING_EXONS );
    assert( check_raw_bitflag_consistency_tests("ID_WITHOUT_PARENT_NOT_GENE_PSEUDOGENE",report) == (NON_PERMITTED_BIOTYPES|ID_WITHOUT_PARENT_NOT_GENE_PSEUDOGENE) );
    assert( check_raw_bitflag_consistency_tests("ID_WITHOUT_PARENT_NOT_GENE_PSEUDOGENE",report) != (NON_PERMITTED_BIOTYPES) );
    assert( check_raw_bitflag_consistency_tests("ID_WITHOUT_PARENT_NOT_GENE_PSEUDOGENE",report) != (0x8) );

    assert( check_raw_bitflag_consistency_tests("NON_UNIQUE_ID",report) == NON_UNIQUE_ID ); // at mRNA level
    assert( check_raw_bitflag_consistency_tests("NON_UNIQUE_ID_2EXON",report) == OVERLAPPING_EXONS ); // duplicated an exon
    assert( check_raw_bitflag_consistency_tests("NON_UNIQUE_ID_3GENE",report) == NON_UNIQUE_ID ); // at gene level

    assert( check_raw_bitflag_consistency_tests("NEGATIVE_COORDINATES",report) == NEGATIVE_COORDINATES ); // at mRNA level
    assert( check_raw_bitflag_consistency_tests("NEGATIVE_COORDINATES_2END",report) == NEGATIVE_COORDINATES ); // at mRNA level
    assert( check_raw_bitflag_consistency_tests("NEGATIVE_COORDINATES_3BOTH",report) == NEGATIVE_COORDINATES ); // at mRNA level
    assert( check_raw_bitflag_consistency_tests("NEGATIVE_COORDINATES_3BOTH",report) != 0x80 );

    assert( check_raw_bitflag_consistency_tests("BLANK_LINES",report) == BLANK_LINES );

    if(system("cp testfiles/APOLLO_SCF_NAMES.gff_ORIG testfiles/APOLLO_SCF_NAMES.gff")!=0) throw runtime_error("couldn't copy file!?!");
    assert( check_raw_bitflag_consistency_tests("APOLLO_SCF_NAMES",report) == APOLLO_SCF_NAMES ); // a bit restrictive - uses scaffold names of form \w+\d and then the coords?!?

    if(system("cp testfiles/LINE_ENDINGS.gff_ORIG testfiles/LINE_ENDINGS.gff")!=0) throw runtime_error("couldn't copy file!?!");
    assert( check_raw_bitflag_consistency_tests("LINE_ENDINGS",report) == LINE_ENDINGS );

    assert( check_raw_bitflag_consistency_tests("NO_GENES",report) == (EXCESS_TRANS_REL_GENE_CONSISTENCY_PROB|NO_GENES) );
    assert( check_raw_bitflag_consistency_tests("NO_GENES",report) != (EXCESS_TRANS_REL_GENE_CONSISTENCY_PROB) );
    assert( check_raw_bitflag_consistency_tests("NO_GENES",report) != 0 );

    assert( check_raw_bitflag_consistency_tests("NO_CDS",report) == (NO_CDS|PROTEIN_CODING_LACKS_CDS) ); // this was only really used with pseudogene thing but not relevant now with PROTEIN_CODING_LACKS_CDS
    assert( check_raw_bitflag_consistency_tests("NO_CDS",report) != (PROTEIN_CODING_LACKS_CDS) ); 

    assert( check_raw_bitflag_consistency_tests("PARTIAL_MODEL",report) == PARTIAL_MODEL );
    assert( check_raw_bitflag_consistency_tests("PARTIAL_MODEL_2ENDPOS",report) == PARTIAL_MODEL );
    assert( check_raw_bitflag_consistency_tests("PARTIAL_MODEL_2ENDPOS",report) != 0x800 );

    assert( check_raw_bitflag_consistency_tests("GFF_FASTA_HEADER",report) == (LINE_ENDINGS|GFF_FASTA_HEADER) );

    assert( check_raw_bitflag_consistency_tests("FASTA_HEADER",report) == (FASTA_HEADER|LINE_ENDINGS) );

    assert( check_raw_bitflag_consistency_tests("EMBL_FORMAT",report) == (LINE_ENDINGS|EMBL_FORMAT) );

    assert( check_raw_bitflag_consistency_tests("NO_FEATLINES",report) == (NO_FEATLINES|LINES_WO_9COLS|NO_GENES|NO_CDS|NO_EXON_CDS|NO_TRANSCRIPTS) );

    assert( check_raw_bitflag_consistency_tests("PSEUDOGENES_PRESENT",report) == PSEUDOGENES_PRESENT );
    assert( check_raw_bitflag_consistency_tests("PSEUDOGENES_PRESENT_2TRANSCRIPT",report) == (NON_PERMITTED_BIOTYPES|EXCESS_CDS_EXON_CONSISTENCY_PROB) );
    assert( check_raw_bitflag_consistency_tests("EXCESS_GENE_CONSISTENCY_PROB",report) == EXCESS_GENE_CONSISTENCY_PROB );
    assert( check_raw_bitflag_consistency_tests("EXCESS_TRANS_REL_CDS_EXON_CONSISTENCY_PROB",report) == (EXCESS_TRANS_REL_CDS_EXON_CONSISTENCY_PROB|TRANSCRIPT_LACKS_EXONS|PROTEIN_CODING_LACKS_CDS) );
    assert( check_raw_bitflag_consistency_tests("NO_EXON_CDS",report) == (NO_EXON_CDS|NO_CDS|TRANSCRIPT_LACKS_EXONS|PROTEIN_CODING_LACKS_CDS|EXCESS_TRANS_REL_CDS_EXON_CONSISTENCY_PROB) );
    assert( check_raw_bitflag_consistency_tests("EXCESS_GENE_CONSISTENCY_PROB_2NOPROTEINCODING",report) == (EXCESS_GENE_CONSISTENCY_PROB|EXCESS_CDS_EXON_CONSISTENCY_PROB) );
    assert( check_raw_bitflag_consistency_tests("NO_TRANSCRIPTS",report) == (NO_TRANSCRIPTS|EXCESS_GENE_CONSISTENCY_PROB|EXCESS_CDS_EXON_CONSISTENCY_PROB) );
    assert( check_raw_bitflag_consistency_tests("PARENT_WITHOUT_ID_NOT_CDS_EXON",report) == (PARENT_WITHOUT_ID_NOT_CDS_EXON|EXCESS_CDS_EXON_CONSISTENCY_PROB) );
    assert( check_raw_bitflag_consistency_tests("PARENT_WITHOUT_ID_NOT_CDS_EXON_2MISNAMEDCDS",report) == (PARENT_WITHOUT_ID_NOT_CDS_EXON|NON_PERMITTED_BIOTYPES) );
    assert( check_raw_bitflag_consistency_tests("EXCESS_TRANS_REL_GENE_CONSISTENCY_PROB",report) == EXCESS_TRANS_REL_GENE_CONSISTENCY_PROB );
    assert( check_raw_bitflag_consistency_tests("TRANSCRIPT_LACKS_EXONS",report) == TRANSCRIPT_LACKS_EXONS );
    assert( check_raw_bitflag_consistency_tests("PROTEIN_CODING_LACKS_CDS",report) == PROTEIN_CODING_LACKS_CDS );

    assert( check_raw_bitflag_consistency_tests("NAMES_HAVE_SPACES",report) == NAMES_HAVE_SPACES );
    assert( check_raw_bitflag_consistency_tests("NAMES_HAVE_SPACES_2NOSPACES",report) == 0 );

    cout << "\nTESTS ARE FINE!?!?!\n";

    return 0; 

/*  
    const char* error_types[] = { 
        "APOLLO_SCF_NAMES",
        "NAMES_HAVE_SPACES",
        "LINES_WO_9COLS",
        "LINE_ENDINGS",
        "ILEGAL_FEAT_TYPES",
        "EXCESS_GENE_CONSISTENCY_PROB",
        "EXCESS_TRANS_REL_GENE_CONSISTENCY_PROB",
        "EXCESS_TRANS_REL_CDS_EXON_CONSISTENCY_PROB",
        "EXCESS_CDS_EXON_CONSISTENCY_PROB",
        "NON_PERMITTED_BIOTYPES",
        "ID_WITHOUT_PARENT_NOT_GENE_PSEUDOGENE",
        "PARENT_WITHOUT_ID_NOT_CDS_EXON",
        "UNKNOWN_SCAFFOLD",
        "NEGATIVE_COORDINATES",
        "NON_PRINTING_X0D",
        "BLANK_LINES",
        "NON_UNIQUE_ID",
        "PSEUDOGENES_PRESENT",
        "CDS_PRESENT",
        "EMBL_FORMAT",
        "GFF_FASTA_HEADER",
        "FASTA_HEADER",
        "PARTIAL_MODEL",
        "NO_FEATLINES",
        "NO_GENES",
        "NO_EXON_CDS",
        "NO_TRANSCRIPTS",
        "TRANSCRIPT_LACKS_EXONS",
        "PROTEIN_CODING_LACKS_CDS",
        "OVERLAPPING_EXONS"
    };

    for (unsigned int i = 0 ; i < sizeof(error_types) / sizeof(char*) ; i++) cout << "type = "<<*(error_types+i)<< " = " << error_types[i] << "\n";// ...
    for (const char** it = error_types ; it < error_types+(sizeof(error_types)/sizeof(char*)) ; it++) cout << "running test for " << *it << " with file\n";

unsignedll named_mask = NAMES_HAVE_SPACES|CDS_PRESENT;
unsignedll named_return = check_raw_bitflag_consistency_tests ("NAMES_HAVE_SPACES",report);
if(NAMES_HAVE_SPACES|CDS_PRESENT==check_raw_bitflag_consistency_tests ("NAMES_HAVE_SPACES",report)) {
if (named_mask==named_return) { 
if(NAMES_HAVE_SPACES|CDS_PRESENT==named_return) {
if(0x40002==check_raw_bitflag_consistency_tests ("NAMES_HAVE_SPACES",report)) {
if (check_raw_bitflag_consistency_tests ("NAMES_HAVE_SPACES",report)==named_mask) { 
unsignedll named_mask = (1ull<<1)|(1ull<<18);
unsignedll named_return = check_raw_bitflag_consistency_tests ("NAMES_HAVE_SPACES",report);
if((1ull<<1)|(1ull<<18)==check_raw_bitflag_consistency_tests ("NAMES_HAVE_SPACES",report)) {
if (named_mask==named_return) {
if((1ull<<1)|(1ull<<18)==named_return) {
if(0x40002==check_raw_bitflag_consistency_tests ("NAMES_HAVE_SPACES",report)) {
if (check_raw_bitflag_consistency_tests ("NAMES_HAVE_SPACES",report)==named_mask) {
*/

}

bool check_capmon_bool_test (const char* filename, std::string& report, DB_PARAMS* dbp) {

    std::string file(filename);
    file = "./testfiles/" + file + ".gff";
    cout << "Boolean check ";
    if (dbp!=0) cout << ": RUNNING SCF NAME CHECK ";
    else cout << ": NOT running scf name check ";
    bool x = capmon_gff_validation (file.c_str(), report, dbp);
    cout << "- using file " << file << " : " << std::boolalpha << x << "\n";
    return x;

}

unsignedll check_raw_bitflag_consistency_tests (const char* filename, std::string& report, DB_PARAMS* dbp) {

    std::string file(filename);
    file = "./testfiles/" + file + ".gff";
    cout << "Bitflag check - using file " << file << "\n";
    gff_holder dummy; // just let it go out of scope?!?
    std::stringstream strstrm(std::stringstream::out);
    unsignedll bitflag = 0;
    name_holder nh;
    bitflag = details::gff_basic_validation_1a_gff_parse (file.c_str(), strstrm, nh, dummy);
    bitflag = details::gff_basic_validation_1b_gff_name_checks (bitflag, strstrm, nh, dummy, dbp);
    
    bitflag = details::gff_basic_validation_1d_individual_model_checks(bitflag, strstrm, dummy);
    report = strstrm.str();
    return bitflag;
}

bool capmon_gff_validation (const char* filename, std::string& report, DB_PARAMS* dbp) {

    gff_holder dummy; 
    std::stringstream strstrm(std::stringstream::out);
    unsigned long long bitflag = 0;

{ // just avoiding persistance of name_holder for no reason?!?
    name_holder nh;

    try {
        bitflag = details::gff_basic_validation_1a_gff_parse (filename, strstrm, nh, dummy);
    } catch (std::length_error) {
        cout << "Seems that the error messages are just too long!?!\n";
        throw;
    }
    bitflag = details::gff_basic_validation_1b_gff_name_checks (bitflag, strstrm, nh, dummy, dbp);
    //y could put in fasta checking?!?

}

    bitflag = details::gff_basic_validation_1d_individual_model_checks(bitflag, strstrm, dummy);
    details::gff_basic_validation_1x_file_cleanup (bitflag,filename,strstrm);

    report = details::capmon_html_table(bitflag,strstrm);
    bool bitwise_stop = bitflag & PROBLEM; 

    if (bitwise_stop) return false;
    else return true;

}

bool standalone_gff_validation (const char* filename, DB_PARAMS* dbp) {

    gff_holder dummy; 
    std::stringstream strstrm(std::stringstream::out);
    unsigned long long bitflag = 0;
    name_holder nh;
    try {
        bitflag = details::gff_basic_validation_1a_gff_parse (filename, strstrm, nh, dummy);
    } catch (std::length_error) {
        cout << "Seems that the error messages are just too long!?!\n";
        throw;
    }
    bitflag = details::gff_basic_validation_1b_gff_name_checks (bitflag, strstrm, nh, dummy, dbp);
    bitflag = details::gff_basic_validation_1d_individual_model_checks(bitflag, strstrm, dummy);
    bool bitwise_stop = bitflag & PROBLEM; 
    print_em_pretty pep;
    if (bitwise_stop) cout << "\nFlags :\n";
    for (unsigned long long i = 0 ; i < sizeof(long long)*8 ; i++) 
      if (bitflag&(1ull<<i)) cout << " " << std::dec << pep.grab_name((1ull<<i)) << " ("<<i<<")\n"; //  << " and " << x << "\n";
    cout << "\nDetails :\n" << strstrm.str() << "\n";
    if (bitwise_stop) return false;
    else return true;
}

bool parse_gff_for_loading (const char* filename, gff_holder *blarp) {

    gff_holder &blah = *blarp;
    std::stringstream strstrm(std::stringstream::out);
    unsigned long long bitflag = 0;
    name_holder nh;
    bitflag = details::gff_basic_validation_1a_gff_parse (filename, strstrm, nh, blah);
    bitflag = details::gff_basic_validation_1b_gff_name_checks (bitflag, strstrm, nh, blah);
    bitflag = details::gff_basic_validation_1d_individual_model_checks(bitflag, strstrm, blah);
    bool bitwise_stop = bitflag & PROBLEM; 
    print_em_pretty pep;
    if (bitwise_stop) cout << "\nFlags :\n";
    for (unsigned long long i = 0 ; i < sizeof(long long)*8 ; i++) 
      if (bitflag&(1ull<<i)) cout << " " << std::dec << pep.grab_name((1ull<<i)) << " ("<<i<<")\n"; //  << " and " << x << "\n";
    cout << "\nDetails :\n" << strstrm.str() << "\n";
    if (bitwise_stop) return false;
    else return true;

}
