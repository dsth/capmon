CC=g++-4.6

LDFLAGS= -L/homes/dsth/dev/NewCap/log4cpp/lib -L/homes/dsth/dev/NewCap/backend/monitor/sqlite/lib \
  -L/net/isilon3/production/panda/ensemblgenomes/development/dsth/boost_1_46_1/stage/lib/ `mysql_config --libs` \
  -L/homes/dsth/dev/curl/lib \
  `/usr/bin/perl -MExtUtils::Embed -e ldopts`

CFLAGS= -I/homes/dsth/dev/NewCap/log4cpp/include -I/homes/dsth/dev/NewCap/backend/monitor/sqlite/include \
  -I/homes/dsth/dev/NewCap/backend/monitor/boost_1_48_0 `~/dev/mysql/mysql-5.1.63/bin/mysql_config --cflags` \
  -I/homes/dsth/dev/curl/include \
  -I/homes/dsth/dev/NewCap/jsoncpp/include \
  `/usr/bin/perl -MExtUtils::Embed -e ccopts`

cc = ${CC}
d= ${D}
extra = -Wall -Werror -Wl,-Bstatic -llog4cpp -lmysqlclient -lboost_regex -lboost_system -lboost_program_options -lboost_filesystem -lsqlite3 -Wl,-Bdynamic -lz -lcrypt -lnsl -lm -lcurl json_static_compile/libjson.a

static = -static -static-libstdc++

w= -Wall -Werror -Wextra -pedantic -Wuninitialized -Winit-self -Wsequence-point -std=c++0x -O -Wno-long-long -Wpointer-arith  -Wmissing-noreturn 

%.o: %.cpp
	 $(cc) $(w) $(d) -g -c -fno-stack-protector $(CFLAGS) $< -o $@

OBJECTS = main.o loop.o email.o gffdocwrap.o \
		  time.o utils.o list.o batch.o perl_bindings.o levenshtein.o gff_validation.o

default: build_capmon

clean : 
	rm capmon $(OBJECTS) 

build_capmon : capmon main.o loop.o email.o gffdocwrap.o time.o utils.o list.o batch.o levenshtein.o gff_validation.o

capmon : $(OBJECTS)
	$(cc) $(w) $(d) -fno-stack-protector -o capmon \
	$(OBJECTS) \
	boost_1_48_0/lib/libboost_system.a \
	boost_1_48_0/lib/libboost_regex.a \
	boost_1_48_0/lib/libboost_filesystem.a \
	static_lib/libsqlite3.a \
	static_lib/liblog4cpp.a \
	static_lib/libcurl.a \
	static_lib/libmysqlclient.a \
	static_lib/libz.a  \
	openssl-1.0.1b/libcrypto.a \
	json_static_compile/libjson.a \
	-lpthread \
	-lrt \
	-static-libstdc++ \
	`/usr/bin/perl -MExtUtils::Embed -e ldopts`

main.o : main.cpp exceptions.h utils.h config.h matrix.h stride_iter.h  time.h toolz.h
loop.o : loop.cpp loop.h
email.o : email.cpp email.h
gffdocwrap.o : gffdocwrap.cpp gffdocwrap.h
time.o : time.cpp time.h
utils.o : utils.cpp utils.h
list.o : list.cpp list.h
batch.o : batch.cpp batch.h
perl_bindings.o : perl_bindings.cpp perl_bindings.h
levenshtein.o : levenshtein.cpp levenshtein.h
gff_validation.o : gff_validation.cpp gff_validation.h
	~/dev/gcc_4.6/bin/g++-4.6 -std=c++0x $(CFLAGS) -Wall -Wextra -Werror gff_validation.cpp -DCAPMON_EXT -lssl -O2 -c -o gff_validation.o

gfftest : test_clean test_extended test

test : tests.cpp 
	~/dev/gcc_4.6/bin/g++-4.6 -std=c++0x  -I/homes/dsth/dev/NewCap/jsoncpp/include -Wall -Wextra -Werror tests.cpp gff_validation.o -lboost_regex `mysql_config --libs` utils.o boost_1_48_0/lib/libboost_system.a boost_1_48_0/lib/libboost_filesystem.a -DCAPMON_EXT -lssl -O2 -o test

test_extended : gff_validation.cpp gff_validation.h 
	~/dev/gcc_4.6/bin/g++-4.6 -std=c++0x  -Wall -Wextra -Werror gff_validation.cpp -DCAPMON_EXT -lssl -O2 -c -o gff_validation.o

test_clean : 
	rm test gff_validation.o 

coord : gff_validation_mini.o converter.o
	~/dev/gcc_4.6/bin/g++-4.6 -Wall -Wextra -ansi -pedantic -Werror -g -static-libstdc++ -O -std=c++0x gff_validation_mini.o converter.o boost_1_48_0/lib/libboost_regex.a -o coord

converter.o : converter.cpp
	~/dev/gcc_4.6/bin/g++-4.6 $(CFLAGS)  -Wall -Wextra -ansi -pedantic -Werror -g -std=c++0x converter.cpp -c -o converter.o

gff_validation_mini.o : gff_validation.cpp
	/homes/dsth/dev/gcc_4.6/bin/g++-4.6 $(CFLAGS) -pedantic -ansi -std=c++0x -g -Wall -Wextra -Werror gff_validation.cpp -O2 -c -o gff_validation_mini.o

