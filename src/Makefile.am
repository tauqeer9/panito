panitoincludedir=$(includedir)/panito
panitoinclude_HEADERS=kseq.h gwani.h

bin_PROGRAMS = panito
panito_SOURCES = main.c
panito_LDADD = libpanito.la -lpthread

# libpanito.so is our library
lib_LTLIBRARIES = libpanito.la
libpanito_la_SOURCES = gwani.c
libpanito_la_LDFLAGS= -version-info 1:0

# "make check" target
TESTS = $(check_PROGRAMS)
check_PROGRAMS = run-all-tests
run_all_tests_SOURCES = \
	../tests/check-panito.c \
	../tests/run-all-tests.c
run_all_tests_CFLAGS = $(CHECK_CFLAGS) $(CFLAGS)
run_all_tests_LDADD = libpanito.la -lpthread $(CHECK_LIBS)
