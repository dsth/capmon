#include <fstream>
#include "boost/regex.hpp"
#include "gffdocwrap.h"
#include <iomanip> 
#include "toolz.h"
#include "gff_validation.h"

#define REGEX_PROCESSED "CAPPROC : \\[INFO\\] : \\d+\\/\\d+\\/\\d+ \\d+:\\d+:\\d+ : Processing (\\d+) genes\\."
#define REGEX_REMOVED "CAPPROC : \\[INFO\\] : \\d+\\/\\d+\\/\\d+ \\d+:\\d+:\\d+ : Removed (\\d+) unprotected old gene models\\."
#define REGEX_PROTECTED "CAPPROC : \\[INFO\\] : \\d+\\/\\d+\\/\\d+ \\d+:\\d+:\\d+ : Protected (\\d+) new gene models\\."
#define REGEX_COM "\\s*#.*"
#define REGEX_LINE_ENDS ".*[^;]\\s*"
#define REGEX_BLANK "\\s*"
#define REGEX_CR "\\r"
#define REGEX_NEGSTART "([^\\t]+)\\t[^\\t]+\\t([^\\t]+)\\t-\\d+\\t\\d+\\t[^\\t]+\\t[^\\t]+\\t[^\\t]+\\t([^\\t]+)"
#define REGEX_NEGEND "([^\\t]+)\\t[^\\t]+\\t([^\\t]+)\\t\\d+\\t-\\d+\\t[^\\t]+\\t[^\\t]+\\t[^\\t]+\\t([^\\t]+)"
#define REGEX_FEAT "([^\\t]+)\\t[^\\t]+\\t([^\\t]+)\\t\\d+\\t\\d+\\t[^\\t]+\\t[^\\t]+\\t[^\\t]+\\t([^\\t]+)"
#define REGEX_SCFNAME "(\\w+?\\d+?):\\d+?[:-]\\d+?"
#define REGEX_ID "ID=([^;]+)"
#define REGEX_PARENT "Parent=([^;]+)"

#define REGEX_EMBL "FT\\s.*"
#define REGEX_GFFFASTA "#{1,2}FASTA.*"
#define REGEX_FASTA ">\\w+.*"

#define SPRINTF_STRING "select count(1) from seq_region where name = '%s'"
#define SQL_STRING "select * from seq_region where name = ''"

using std::setw;
using std::setfill;
using std::left;
using std::right;
using std::cout;
using std::endl;
using std::runtime_error;

// using boost::regex_match; // prolly shouldn't do this!?
typedef boost::regex regex;
typedef boost::smatch smatch;

#ifndef PERMITTED_TRANSCRIPT_TYPES_REGEX
#define PERMITTED_TRANSCRIPT_TYPES_REGEX "mRNA|tRNA|pseudogenic_tRNA|rRNA|miRNA|ncRNA|pseudogene"
#endif

#ifndef IGNORE_TYPES_REGEX
#define IGNORE_TYPES_REGEX "contig|supercontig|match|match_part"
#endif

int int_from_mysql_query (DB_PARAMS *, const char * query); // sheer lazyness

std::ostream & nl(std::ostream& os) { return os << " \\\n"; } // nullary maninulator, should probably use an effector

/* 
std::ostream & tb(std::ostream& os) { return os << " | ";   }
std::ostream & ts(std::ostream& os) { return os << "| ";    }
std::ostream & te(std::ostream& os) { return os << " |\n";  }

class bool_string {
	bool status; // implicit convertion to boo?!?
  public:
	bool_string(bool width) : status(width) {}
	friend std::ostream& operator<<(std::ostream& os, const bool_string & b_obj) {
		std::string s;
        s = b_obj.status ? "true" : "false";
		return os << s;
	}
};

class print_row {
	bool boolean; 
    std::string issue;
    std::string item;
public:
	print_row(std::string i, std::string a, bool b) : boolean(b), issue(a), item(i) {}
	friend std::ostream& operator<<(std::ostream& os, const print_row & b_obj) {
        std::string td(b_obj.boolean?"#FF927D":"#e5f1ff");
		std::string s(b_obj.boolean?"TRUE":"");
        td="<td BGCOLOR=\"" + td + "\" style=\"padding:10px;\" >";
		return os << " <tr> " << td 
          << b_obj.item << " </td> " << td << b_obj.issue 
          << " </td>" << td 
          << s << "</td></tr>\n"; ///r break-up table with '\n' to avoid that odd html email formating issue?!?
	}
};
*/

typedef std::set<std::string> FEAT_LIST; 

void list_from_singlecol_mysql_query (DB_PARAMS* dbp, const char* query, std::string& out) { // prefer string& param to return?!?

    try {

        MYSQL *conn;
        MYSQL_RES *result;
        MYSQL_ROW row;
        // MYSQL_FIELD *field;
        conn = mysql_init(NULL);

        if(mysql_real_connect(conn,dbp->host(),dbp->user(),dbp->pass(),dbp->dbname(),dbp->port(),NULL,0) == NULL)
          throw MySqlError("couldn't connect to database. exiting.");

            if (mysql_query(conn,query))
              throw MySqlError("unable to access e! cap db");

            if(!(result = mysql_store_result(conn))) { //y get the result set
                throw runtime_error ("unable to query e! cap database");
            } else {
                while ((row = mysql_fetch_row(result))) { // batch statement?!?
                    out += std::string(row[0]) + ", ";
                }
            // mysql_free_result(result);
        }

        mysql_free_result(result);
        mysql_close(conn);

    } catch (std::runtime_error e) {
        std::string prob("Exception thrown : ");
        throw; // propagate up...
    }

}

int int_from_mysql_query (DB_PARAMS * dbp, const char * query) {

    int query_int;

    try {

        MYSQL *conn;
        MYSQL_RES *result;
        MYSQL_ROW row;
        conn = mysql_init(NULL);

        if(mysql_real_connect(conn,dbp->host(),dbp->user(),dbp->pass(),dbp->dbname(),dbp->port(),NULL,0) == NULL)
          throw MySqlError("couldn't connect to database. exiting.");

        if (mysql_query(conn,query))
          throw MySqlError("unable to prepare query for int");

        if(!(result = mysql_store_result(conn))) { //y get the result set
            throw runtime_error ("unable to retrieve int from database");
        } else {
            row = mysql_fetch_row(result);
            if (row == 0) throw runtime_error ("null pointer returned on db query!");
            else query_int = atoi(row[0]);
        }

        mysql_free_result(result);
        mysql_close(conn);

    } catch (std::runtime_error e) {
        std::string prob("Exception thrown : ");
        throw; 
    }

    return query_int;
}

int GFFPROC::validate (char* filename, int sbmid, DB_PARAMS* db_params, std::string& email_body, std::string& tmpmsg) { 

    dbp = db_params; 
    char msg[MED_STRING_SIZE]; // this is very silly?!?
    snprintf(msg, MED_STRING_SIZE, "Running initial gff validation step on file '%s'",filename);
    log4.debug(msg);

    if(!capmon_gff_validation(filename,email_body,db_params)) return FILE_GFF3VAL1; 

    log4.debug("File passed first round validation");
    log4.debug("Running secondary validation and gene loading");

    email_body.clear(); 
    int saved_genes = (*this).gffload(filename,email_body); 

    // wtf is this doing here, was a quick hack to try and speed up gene model error oddities
    if(saved_genes==0) {
        tmpmsg = "<pre>";
        char tmp[2500];
        std::string cmd("./Validation2Probs.pl ");
        cmd += filename;
        cmd += " cap";
        FILE * f2 = popen(cmd.c_str(), "r");
        while(( fgets( tmp, 2500, f2 )) != NULL ) tmpmsg += tmp;
        pclose(f2);
        tmpmsg += "</pre>";
        return FILE_GFF3VAL2; 
    }

    log4.debug("File passed second round validation");
    sprintf(msg, "<p>Saved %d gene(s) (total number of gene models in gff submission).</p>", saved_genes);
    log4.debug(msg);
    email_body = msg;

    std::pair<int,int> cleanup = this->eclean(filename,sbmid);
    log4.debug("Fingerprinting gene(s) models for modified model detection");
    if(cleanup.first==0) return FILE_GFF3NONEWGENES;

    sprintf(msg, "<p>Protected %d new gene(s) models (total number of 'novel' gene models).</p>", cleanup.first); 
    log4.debug(msg);
    email_body += "\n";
    email_body += msg;

    sprintf(msg, "<p>Removed %d gene(s) (gene models identical to models in canonical set - to protect a gene model without changing it you must use the gene info form).</p>", cleanup.second);
    log4.debug(msg);
    email_body += "\n";
    email_body += msg;

    /*
    // just cos of that really freakish occurance when the servers went nuts!
    if(saved_genes<cleanup.first) 
      mailobj->send_email("SUBMISSION DOESN'T MAKE SENSE!\r\n","There are more protected models that saved models!?!.\r\n", EMAIL_ADMIN);
    */

    sprintf(msg,SQL_NONCDSNUM,sbmid);
    int non_cds = int_from_mysql_query(dbp,msg);
    if (non_cds>0) {
        sprintf(msg, "<p>Submission has %d gene(s) with stop codons (these will not be included in future patches to canonical gene set).</p>", non_cds);
        log4.debug(msg);
        email_body += "\n";
        email_body += msg;
        sprintf(msg,SQL_NONCDSNAMES,sbmid);
        std::string temp;
        list_from_singlecol_mysql_query (dbp,msg, temp);
        return FILE_NONCDS;
    }

    return FILE_DONE;
}

// such a mess?!?!? why's it even a class?!?
int GFFPROC::gffload (const char* filename ,std::string& exmsg) {
     
/*
#ifdef DBG2
    std::cout << "\nSilly=" << __FILE__ << ":" << __LINE__ << " : address of filename=" << &filename << "\n" << std::endl;
#endif
*/

    std::string log(filename);
    log += SUFFIX_LOADING;

    ///r make this a modifier!?!
    
    toolz::SQLITE_ADAPTOR::connect(const_cast<char*>(this->capdb)); // only to remind that its a data member?!?
        // string cmd("PERL5LIB=" + perl5lib + " " + perl + " -M"+ *it +" -e 'print q{blah}' &> /dev/null");
    std::string perl5lib = toolz::SQLITE_ADAPTOR::get_instance()->meta_value("local.perl5lib");
    toolz::SQLITE_ADAPTOR::disconnect();
    std::for_each(perl5lib.begin(),perl5lib.end(),[](char& n){if(n==',')n=':';});

    std::stringstream cmd(std::stringstream::out);

    // eew
    cmd << "PERL5LIB=" << perl5lib << " " << gffdoccmd << " \\\n"
      << " --dbname " << dbp->dbname() << nl
      << " --host " << dbp->host() << nl
      << " --user " << dbp->user() << nl
      << " --port " << dbp->port() << nl
      << " --pass " << dbp->pass() << nl
      << " --file " << filename << nl
      << "&> " << log;
    
/*
#ifdef DBG
    std::cout << "\nSilly=" << __FILE__ << ":" << __LINE__ << "\n" << std::endl;
    std::cout << cmd.str() << "\n" << std::endl;
#endif
*/

    if (system(cmd.str().c_str()) != 0)
      throw runtime_error ("Problem runnning GffDoc.pl. See : " + log);

    std::ifstream in(log.c_str());

    if(in == 0) 
      throw runtime_error ("Cannot find log file for gff load - something went wrong.");

    smatch match_obj;
    regex reg_saved("main:: > Saved (\\d+) genes\\.");
    // too lazy to re-write this 
    regex reg_ex_start("Class:   Exception"); 
    regex reg_ex_stop("Msg: "); 

    int genes_saved = 0;
    std::string s;
    bool exmsgflag = false;
    while(getline(in, s)) {
        if (regex_search(s,match_obj,reg_saved))
          genes_saved = atoi(std::string(match_obj[1]).c_str());
        if (regex_search(s,reg_ex_start)) exmsgflag = true;
        if (exmsgflag) exmsg += s + "\n";
        if (regex_search(s,reg_ex_stop)) exmsgflag = false;
    }

    if (!exmsg.empty()) exmsg = "<pre>\nPARSING EXCEPTION:\n" + exmsg + "</pre>";

    return genes_saved;
}

std::pair<int,int> GFFPROC::eclean (const char * filename, int sbmid) {
     
    std::string log(filename);
    log += SUFFIX_CLEANING;

    toolz::SQLITE_ADAPTOR::connect(const_cast<char*>((*this).capdb));
    std::string perl5lib = toolz::SQLITE_ADAPTOR::get_instance()->meta_value("local.perl5lib");
    toolz::SQLITE_ADAPTOR::disconnect();
    std::for_each(perl5lib.begin(),perl5lib.end(),[](char& n){if(n==',')n=':';});

    std::stringstream cmd(std::stringstream::out);
    cmd << "PERL5LIB=" << perl5lib << " " << efingercmd 
      << " --dbname " << dbp->dbname() << nl
      << " --host " << dbp->host() << nl
      << " --user " << dbp->user() << nl
      << " --port " << dbp->port() << nl
      << " --pass " << dbp->pass() << nl
      << " --sbmid " << sbmid << " \\\n" 
      << "&> " << log;

    std::stringstream strstrm(std::stringstream::out);

/*
#ifdef DBG
    std::cout << "\nSilly=" << __FILE__ << ":" << __LINE__ << "\n" << std::endl;
    std::cout << cmd.str() << "\n" << std::endl;
    exit(0);
#endif
*/

    if (system(cmd.str().c_str())!=0)
      throw runtime_error ("Problem database cleanup of old gene models. See : " + log);

    std::ifstream in(log.c_str());

    if(in == 0)
      throw runtime_error ("Cannot find log file for cap db cleanup - something went wrong.");

    smatch match_obj;
    regex reg_processed(REGEX_PROCESSED);
    regex reg_removed(REGEX_REMOVED);
    regex reg_protected(REGEX_PROTECTED);

    signed int genes_processed = -1;
    signed int genes_protected = -1;
    signed int genes_removed = -1;
    std::string s;

    while(getline(in, s)) {
        if (regex_search(s,match_obj,reg_protected))
          genes_protected = atoi(std::string(match_obj[1]).c_str());
        if (regex_search(s,match_obj,reg_processed))
          genes_processed = atoi(std::string(match_obj[1]).c_str());
        if (regex_search(s,match_obj,reg_removed))
          genes_removed = atoi(std::string(match_obj[1]).c_str());
    }

    if(genes_processed==-1||genes_protected==-1||genes_removed==-1)
      throw runtime_error ("was not able to parse efingerprint output!");

    if(genes_processed!=(genes_protected+genes_removed))
      throw runtime_error ("efingerprint cleanup numbers don't add up!");

    std::pair<int,int> tmp(genes_protected,genes_removed);
    return tmp;

}

/* 
const char GFFPROC::summary_template[] = "\n"
"====================================================================\n"
"| Validation 1 Report: issues marked with '*' are tolerated (atm!) |\n"
"====================================================================\n"
"| Does not contain gene features                           | %5s |\n" 
"| Does not contain transcript features                     | %5s |\n" 
"| Does not contain exon/CDS features                       | %5s |\n" 
"| There are non-comment, non-9-col format lines present    | %5s |\n"
"| Gene id excess vs. transcript parent                     | %5s |\n" 
"| Transcript parent excess vs. gene id                     | %5s |\n" 
"| Transcript id excess vs. exon/CDS parent                 | %5s |\n" 
"| Exon/CDS parent excess vs. transcript id                 | %5s |\n" 
"| Unknown scaffold names                                   | %5s |\n"
"| Non-permitted biotypes                                   | %5s |\n"
"| Negative coordinates                                     | %5s |\n"
"| Non-unique gene IDs                                      | %5s |\n"
"| Non-unique transcript IDs                                | %5s |\n"
"| Missing CDS                                              | %5s |\n"
"--------------------------------------------------------------------\n"
"| Blank lines *                                            | %5s |\n"
"| Non-printing carriage return *                           | %5s |\n"
"| There are lines that don't end with ';' character *      | %5s |\n" 
"| There are spaces in id or parent names *                 | %5s |\n"
"| Scaffold name format used apollo habit *                 | %5s |\n"
"--------------------------------------------------------------------\n\n";
*/

