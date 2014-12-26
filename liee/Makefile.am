AM_CPPFLAGS = \
    -I@with_boinc@/api \
    -I@with_boinc@/lib \
    -D"LOG_ENABLED=1" \
    -D"ALGLIB=1"

module_files = \
	module_config.hpp		module_config.cpp \
	module_factory.hpp \
	potential.hpp			potential.cpp \
	wave_function.hpp		wave_function.cpp \
	noumerov.hpp			noumerov.cpp \
	solver.hpp			solver.cpp \
	observer.hpp			observer.cpp

boinc_objects = \
	@with_boinc@/api/.libs/boinc_api.o \
	@with_boinc@/lib/libboinc.a \
	@with_boinc@/lib/libboinc_crypt.a

bin_PROGRAMS = liee_worker liee_tests

liee_worker_SOURCES =  liee_worker.cpp my_util.hpp my_util.cpp optimizer.hpp optimizer.cpp $(module_files)
liee_worker_LDADD = $(boinc_objects)

liee_tests_SOURCES =  liee_tests.cpp my_util.hpp my_util.cpp optimizer.hpp optimizer.cpp $(module_files)
liee_tests_LDADD = $(boinc_objects) -lboost_unit_test_framework
