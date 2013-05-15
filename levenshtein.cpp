#include <vector>
#include <cstdlib>
#include "matrix.h"
#include <tuple>
#include <fstream>
#include "perl_bindings.h"
#include <algorithm>
#include "boost/regex.hpp"
#include <mysql/mysql.h> 
#include "utils.h"
#include "levenshtein.h" 
#include "exceptions.h"
#include "config.h"

#define STRING_LENGTH 250

// using boost::regex_match;
typedef boost::regex regex;
typedef boost::smatch smatch;

struct SUBSTITUTION {
    int pos;
    char from;
    char too;
    SUBSTITUTION (int p ,char f,char t) : pos(p), from(f), too(t) {}
};

typedef std::pair<int,char> CHANGE;
typedef std::tuple<int,char,char> SUB;
typedef std::vector<CHANGE> CHANGES;
typedef std::vector<SUBSTITUTION> SUBSTITUTIONS;
typedef std::tuple<SUBSTITUTIONS,CHANGES,CHANGES> LIST_OF_CHANGES;
typedef std::pair<int,std::string> UPDATE;
typedef std::vector<UPDATE> UPDATES;

long update_insert_or_delete_from_mysql (const char* host, const char* user, const char* pass, const char* db, int port, const char* query) {

    MYSQL *conn; // MYSQL_RES *result; // MYSQL_ROW row; // MYSQL_FIELD *field;

    try {

        conn = mysql_init(NULL);

        if(mysql_real_connect(conn, host, user, pass, db, port, NULL,0) == NULL)
          throw std::runtime_error("couldn't connect to database. exiting.");
        if (mysql_query(conn,query))
          throw std::runtime_error("unable to access e! cap db");

        long affected = (long) mysql_affected_rows(conn);
        
        mysql_close(conn);
        return affected;

    } catch (std::runtime_error e) {
        mysql_close(conn);
        std::string prob("Exception thrown : ");
        throw;
    }

}

//using std::cout;
//using std::cin;
//using std::endl;
using namespace std;

int levenshtein_distance_and_backtracking(const string& s1, const string& s2, LIST_OF_CHANGES& list) {

    cout << "levenshtein_distance_and_backtracking called"<<endl;
    const size_t len1 = s1.size(), len2 = s2.size();
    matrix<int> d(len1+1,len2+1);

    d[0][0] = 0;
    for(unsigned int i = 1; i <= len1; ++i) d[i][0] = i;
    for(unsigned int i = 1; i <= len2; ++i) d[0][i] = i;

    for(unsigned int i = 1; i <= len1; ++i)
      for(unsigned int j = 1; j <= len2; ++j)
        d[i][j] = std::min( std::min(d[i - 1][j] + 1,d[i][j - 1] + 1),
        d[i - 1][j - 1] + (s1[i - 1] == s2[j - 1] ? 0 : 1) );

    cout << "generated matrix with final edit distance of "<<d[len1][len2]<<endl;

    SUBSTITUTIONS& subs = std::get<0>(list);
    CHANGES& inserts = std::get<1>(list);
    CHANGES& deletions = std::get<2>(list);
    unsigned int i = (unsigned int)len1;
    unsigned int j = (unsigned int)len2;

    int f = 0;
    while (1) {

        int current = d[i][j];
        int del = d[i-1][j];
        int ins = d[i][j-1];
        int diag = d[i-1][j-1];

        if (i <= 1 && j <= 1) {
            if (s1[i-1]!=s2[j-1]) { 
                std::cout << "> substitution at "<<i<<std::endl;
                subs.push_back(SUBSTITUTION(i-1,s1[i-1],s2[j-1])); // subs.push_back(std::make_tuple(i,s1[i-1],s2[j-1]));
            }
            break;
        }

        // if in doubt go top-left?!?
        if (diag <= ins && diag <= del) {
            if (i>1) i--;
            if (j>1) j--; 
            if (current==diag) { 
            } else if (current == diag+1) {
                std::cout << "substitution at "<<i<<std::endl;
                subs.push_back(SUBSTITUTION(i,s1[i],s2[j])); // subs.push_back(std::make_tuple(i,s1[i-1],s2[j-1]));
            } else { 
                cout << "\n(1) current="<<current<<" ins="<<ins<<" del="<<del<<" diag="<<diag<<" : ";
                throw "arse"; // throwing a string literal?!
            } 
        } else if (ins<=del) { 
            if (j>1) j--;
            std::cout << "insert at " << i << std::endl; 
            inserts.push_back(std::make_pair(i,s2[j]));
        } else if (ins>=del) { 
            if (i>1) i--;
            std::cout << "deletion at " << i <<std::endl; 
            deletions.push_back(std::make_pair(i,s1[i]));
        } else { 
            cout << "\n(2) current="<<current<<" ins="<<ins<<" del="<<del<<" diag="<<diag<<" : ";
            throw;
        }
        f++;
        if (f>=15000) break; // put limits on how far it bothers?!?
    }

    cout << "levenshtein_distance_and_backtracking returning"<<endl;
    return d[len1][len2];
}

void levenshtein_distance(DB_PARAMS& dbp, std::string& gfffile, int transcript_id, char* stable_id,std::string& summary,MYSQL* conn) {

    std::ifstream infs(gfffile.c_str());
    regex gff_regex("\\tmRNA\t.*ID=([^;]+)");
    regex seqedit_regex("Seq=([^;]+);");
    regex stopcodon_regex("\\*");
    smatch match;
    if (infs==0) throw std::runtime_error("cannot find gff file");

    std::string line;
    std::string gff_corrected_seq;
    bool found = false;
    while (getline(infs,line)) {
        if (boost::regex_search(line,match,gff_regex)) {
            if(match[1]==stable_id) {
                std::cout << line << std::endl; 
                std::string sid(match[1]);
                found=true;
                if (boost::regex_search(line,match,seqedit_regex)) {
                    summary += "<p>Transcript " + sid + " mRNA feature has Seq key.</p>\n";
                    gff_corrected_seq=match[1];
                } else {
                    summary += "<p>Transcript " + sid + " mRNA feature does not have Seq key with corrected sequence - will not generate seqedit(s).</p>\n";
                    return;
                }
            }
        }
    }

    if(!found) throw std::runtime_error("unable to find mRNA for transcript " + std::string(stable_id)); // throw?!? - i.e. shouldn't hapeen?!?

    if(gff_corrected_seq.empty()) { 
        summary += "<p> Will not process empty sequence</p>\n";
        return;
    }

    cout << "Wiped " << update_insert_or_delete_from_mysql (
        dbp.host(),dbp.user(),dbp.pass(),dbp.dbname(),dbp.port(),
        std::string("delete from transcript_attrib where transcript_id =  " + to_string(static_cast<long long int>(transcript_id)) +  " and attrib_type_id = 145").c_str() 
    ) 
      << " previously existing transcript_attrib features for transcript " << to_string(static_cast<long long int>(transcript_id)) << endl;

    summary += "<p> * Corrected sequence is " + to_string(gff_corrected_seq.size()) + " bases long : " + gff_corrected_seq + "</p>\n"; 

    string api_current_seq = perl_grab_sequence_via_api(dbp.host(),dbp.user(),dbp.pass(),dbp.dbname(),dbp.port(),transcript_id,false);
    summary += "<p> * Direct api translation of gff model sequence is " + to_string(api_current_seq.size()) + " bases long : " + api_current_seq + "</p>\n";
    
    LIST_OF_CHANGES things;

    if (api_current_seq==gff_corrected_seq) {
        summary += "<p> * Sequences are identical - no reason to proceed - should check there are no more stops?!?</p>\n";
        return;
    }

    int edit_dist = levenshtein_distance_and_backtracking(api_current_seq.c_str(),gff_corrected_seq.c_str(),things);
    
    if (edit_dist==0) {
        cout <<"edit distance is zero - shouldn't ever get here"<<endl; // silly to get here without checking?!?
        return;
    } else if (edit_dist>10) {
        summary += "<p> * Edit distance of " + to_string(edit_dist) + " is too high - will not generate seqedits.</p>\n";
        return;
    } else if (edit_dist>0) {
        summary += "<p> * Edit distance between gff_corrected_seq and api_current_seq = " + to_string(edit_dist) + "</p>\n";
        std::cout << "<p> * Edit distance between gff_corrected_seq and api_current_seq " << to_string(edit_dist) << "</p>\n";
    } else throw;

    UPDATES updates;
    int counter=0;

    for (auto it = std::get<0>(things).begin() ; it != std::get<0>(things).end() ; it++) { 
        cout << "# substitution @ "<<it->pos<< " "<<it->from<<"->"<<it->too<<endl;
        updates.push_back(make_pair(it->pos,"insert into transcript_attrib values ("+to_string(static_cast<long long int>(transcript_id))
          +",145,'"+to_string(static_cast<long long int>(it->pos+1))+" "+to_string(static_cast<long long int>(it->pos+1))+" "+it->too+"');"));
        summary += "<p>  * (" + to_string(++counter) + ") Substitution @ " + to_string(it->pos) + " = " + it->from + "->" + it->too + "</p>\n";
    }

    for (auto it = std::get<1>(things).begin() ; it != std::get<1>(things).end() ; it++) {
        cout << "# insertion @ "<<it->first<<" = "<<it->second<<endl;
        updates.push_back(make_pair(it->first,"insert into transcript_attrib values ("+to_string(static_cast<long long int>(transcript_id))
          +",145,'"+to_string(static_cast<long long int>(it->first+1))+" "+to_string(static_cast<long long int>(it->first))+" "+it->second+"');"));
        summary += "<p>  * (" + to_string(++counter) + ") Insertion @ " + to_string(it->first) + " = " + it->second + "</p>\n";
    }

    for (auto it = std::get<2>(things).begin() ; it != std::get<2>(things).end() ; it++) { 
        cout << "# deletion @ "<<it->first<<" = "<<it->second<<endl;
        updates.push_back(make_pair(it->first,"insert into transcript_attrib values ("+to_string(static_cast<long long int>(transcript_id))
          +",145,'"+to_string(static_cast<long long int>(it->first+1))+" "+to_string(static_cast<long long int>(it->first+1))+" ');"));
        summary += "<p>  * (" + to_string(++counter) + ") Deletion @ " + to_string(it->first) + " = " + it->second + "</p>\n";
    }

    std::sort(updates.begin(), updates.end());

    for (auto it = updates.begin() ; it != updates.end() ; it++) { 

        cout << "Inserted " << update_insert_or_delete_from_mysql (
            dbp.host(),dbp.user(),dbp.pass(),dbp.dbname(),dbp.port(),
            it->second.c_str()  
        ) 
          << " added seqedit for base " << to_string(static_cast<long long int>(it->first)) << endl;

        cout <<it->second<<" # "<< it->first << endl;

    }

    string api_seqedit_seq = perl_grab_sequence_via_api(dbp.host(),dbp.user(),dbp.pass(),dbp.dbname(),dbp.port(),transcript_id, false);

    const char* sep = "\n-------------------------------\n";
    if (api_seqedit_seq!=gff_corrected_seq) {

        cout << "seqedits have not worked!?!\ngff_corrected_seq:\n"<<gff_corrected_seq<<"\napi_current_seq:\n"<<api_current_seq<<"\napi_seqedit_seq:\n"<<api_seqedit_seq<<endl;

        LIST_OF_CHANGES blah;

        for (auto it = std::get<0>(blah).begin() ; it != std::get<0>(blah).end() ; it++) { 
            cout << "# substitution @ "<<it->pos<< " "<<it->from<<"->"<<it->too<<endl;
        }

        for (auto it = std::get<1>(blah).begin() ; it != std::get<1>(blah).end() ; it++) {
            cout << "# insertion @ "<<it->first<<" = "<<it->second<<endl;
        }

        for (auto it = std::get<2>(blah).begin() ; it != std::get<2>(blah).end() ; it++) { 
            cout << "# deletion @ "<<it->first<<" = "<<it->second<<endl;
        }

        exit(1);

    } else cout << sep << "SEQEDITS SEEM TO HAVE WORKED!?!"<< sep << endl;

    string translatedCDS = perl_grab_sequence_via_api(dbp.host(),dbp.user(),dbp.pass(),dbp.dbname(),dbp.port(),transcript_id, true);

    std::cout << translatedCDS<<std::endl; 
    summary += "<p> * Translation of coding sequence after application of seqedits to spliced sequence : " + translatedCDS + "</p>";

    if (boost::regex_search(translatedCDS,stopcodon_regex)) {
        summary += "<p> * Translation still contains stop codons - will not re-categorise as protein_coding.</p>\n";
    } else {
        summary += "<p> * Translation free of stop codons - will check all other transcripts.</p>\n";

        MYSQL_RES *result;
        MYSQL_ROW row;

        char q[STRING_LENGTH];
        sprintf(
          q,
          "select t.transcript_id, g.gene_id from gene g, transcript t "
            "where g.gene_id=t.gene_id and g.gene_id = (select gene_id "
            "from transcript where transcript_id = %d);",
          transcript_id
        );

        mysql_query(conn, q);
        // to retreive all results at once use mysql_store_result and to reteive row-by-row use mysql_use_result
        result = mysql_store_result(conn); 
        if (result==0) throw runtime_error("couldn't access e! cap db to retrieve gene and transcript ids");

        int gene_id=0; 
        vector<int> transcript_ids;

        while ((row = mysql_fetch_row(result))) { 
            if (!gene_id) gene_id=atoi(row[1]);
            if (transcript_id==atoi(row[0])) continue;

            std::cout << "gene_id "<< row[1] << " and " << row[0] << std::endl; 

            if (boost::regex_search(perl_grab_sequence_via_api(dbp.host(),dbp.user(),dbp.pass(),dbp.dbname(),dbp.port(),atoi(row[0]), true),stopcodon_regex)) {
                std::cout << "translation of transcript " << row[0] << " has stop codons - will not re-categorise gene model." <<std::endl; 
                summary += "<p>  * translation of transcript " + std::string(row[0]) + " has stop codons - will not re-categorise gene model.</p>\n"; 
                mysql_free_result(result);
                return;
            } else transcript_ids.push_back(atoi(row[0]));
        }
        
        mysql_free_result(result);
        transcript_ids.push_back(transcript_id);
        std::cout << "Checked all transcripts for gene " << gene_id << " will re-classify gene and all other transcripts as protein_coding" << std::endl; 
        summary += "<p>  * Checked all transcripts for gene " + to_string(gene_id) + " will re-classify gene and all other transcripts as protein_coding</p>\n"; 

        // getting really lazy!?! should use sprintf or stringstream?!?
        cout << "updated " << update_insert_or_delete_from_mysql ( 
            dbp.host(),dbp.user(),dbp.pass(),dbp.dbname(),dbp.port(),std::string("update gene set biotype = 'protein_coding' where gene_id = " + to_string(gene_id) + ";").c_str()) << " gene ids"<<std::endl; 
        std::for_each(transcript_ids.begin(),transcript_ids.end(),[&dbp](int i){ 
            cout << "updated " << update_insert_or_delete_from_mysql (
              dbp.host(),dbp.user(),dbp.pass(),dbp.dbname(),dbp.port(), 
              std::string("update transcript set biotype = 'protein_coding' where transcript_id = " + to_string(i) + ";").c_str()) << " transcripts ids"<<std::endl; 
        });

    }

    return;
}

