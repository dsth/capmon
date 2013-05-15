#ifndef CONFIG_H
#define CONFIG_H

#define EMAIL_ADMIN "dsth@cantab.net"
#define EMAIL_LISTSERVE "cap_qc@vectorbase.org" 
#define EMAIL_CAP "cap_process_noreply@vectorbase.org"

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)
#define THROW_DEFAULT(x) #x " @ " __FILE__ ":" TOSTRING(__LINE__) 

#define THROW_CSTRING(x) std::string(x) + " @ " __FILE__ ":" TOSTRING(__LINE__)

#define EXEC "capmon"

#define CAP_SPECIES_KEY "cap.species_version"

#define SUFFIX_CLEANING "_cleaning.log"
#define SUFFIX_LOADING "_loading.log"
#define LOCAL_GFF_SUFFIX ".slim"

#define LOG4_FILE "cap_proc.log"
#define SQLITE_DB_NAME ".cap.db"
#define BASE_URL "http://vectorbase-cap.ensemblgenomes.org/?q=uploaded"
#define TYPE_IGNORE "\\t(contig|protein_match|match|match_part|expressed_sequence_match|five_prime_utr|five_prime_UTR|three_prime_UTR|three_prime_utr)\\t"
#define DUMMY_SBM_ID 0

#define CAP_STATUS_READY 1
#define CAP_STATUS_LOOPING 2 
#define CAP_STATUS_LOCKED 3

#define SQL_CAP_STATUS "insert into cap_status (status) values (%d)"
#define SQL_UPDATE_SBM "update cap_files set sbm_status = %d  where sbm_id = %d"
#define COOKIE "cookie.txt"

#define PULL_SBM_INFO "select feature_no, new_mod_genes \
                       from cap_files \
                       where sbm_id = ?"
#define PULL_SBM_SQL "select feature_no, new_mod_genes \
                       from cap_files \
                       where sbm_id = %s"

#define SQL_QUERY "select sbm_id, user, uid, md5, species, file, ts, type from cap_files \
                        where status = 'new'"

#define SQL_NONCDSNUM "select count(1) from gene where biotype != 'protein_coding' and cap_sbm_id = %d"
#define SQL_NONCDSNAMES "select stable_id from gene where biotype != 'protein_coding' and cap_sbm_id = %d"
#define SQL_CHECK_READY "select count(1) from gene where cap_sbm_id is null"

#define SQL_QUERY_STATUS "select status from cap_status where status_id \
  = (select max(status_id) from cap_status)"

#define SQL_CHECK_DB "select count(*) from cap_status"

#define SQL_INSERT_STATUS_READY "insert into cap_status (status) values ('ready')"
#define SQL_INSERT_STATUS_LOOPING "insert into cap_status (status) values ('looping')"

#define PULL_SBM_INFO "select feature_no, new_mod_genes \
                       from cap_files \
                       where sbm_id = ?"
#define PULL_SBM_SQL "select feature_no, new_mod_genes \
                       from cap_files \
                       where sbm_id = %s"

#define SQL_QUERY "select sbm_id, user, uid, md5, species, file, ts, type from cap_files \
                        where status = 'new'"

#define SQL_GET_COUNT "select count(1) from cap_files"
#define SQL_GET_MAX_SBM "select max(sbm_id) from cap_files"

#define SQL_GET_SPECIES "select meta_value from meta where meta_key = '" CAP_SPECIES_KEY "';"
#define SQL_GET_DSN "select dsn from cap_species where species = '%s'"
#define SQL_SET_STATUS "update cap_files set sbm_status == %d where sbm_id == %d"
#define SQL_GET_LOCALFILES "select sbm_id, submitter_name, user_email, species, file_name, file_type from cap_files where sbm_status == %d"
#define SQL_CHECK_INIT "select count(1) from canon_crc64"
#define SQL_CHECK_READY "select count(1) from gene where cap_sbm_id is null"

#define REGEX_SBM "[\\[]?\\{\"sbm_id\":\"(\\d+)\",\"submitter_name\":\"(.*)\",\"user\
\":\"(.*)\",\"user_email\":\"(.*)\",\"uid\":\"\\d+\",\"ip\":\".*\",\"species\":\"(.*)\",\"file_name\
\":\"(.*)\",\"file_type\":\"(fasta|gff3|xls)\",\"file_md5\":\"(.*)\",\"file_size\
\":\"(\\d+)\",\"file_desc\":\"(.*)\"\\}?"

#define REGEX_MIN "\\[(.*)\\]"
#define REGEX_ACCESS ".*Access\\sdenied.*"
#define REGEX_HTML ".*?<body><pre>(.*)(</pre>)?</body></html>.*"

#define FILE_NEW 0 // newly downloaded
#define FILE_ASCI 1
#define FILE_MD5 2 // if problem wait ~30s?!? and try further ~2?!? times if still arsing 
#define FILE_LOCALSTORE 3 // file pulled and SFs stripped out
#define FILE_TYPE 4
#define FILE_GFF3VAL1 5
#define FILE_GFF3VAL2 6
#define FILE_GFF3NONEWGENES 7
#define FILE_NONCDS 8
#define FILE_DONE 9
#define FILE_IGNORE 10
#define FILE_EMPTY 11
#define FILE_MD5IGNORE 12

#define BLAH1 1
#define BLAH2 2

#endif

