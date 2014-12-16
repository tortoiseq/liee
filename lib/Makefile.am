# recompiles source here with -fPIC (position independent code for shared library)
# and deactivated logging (-D"LOG_DISABLED")

AM_CPPFLAGS = \
    -I../liee \
    -I@with_boinc@/lib \
    -I@with_boinc@/api \
    -D"LOG_DISABLED=1" \
    -D"ALGLIB=1" \
    -fPIC

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

liee_objects = \
	libliee.o \
	my_util.o \
	module_config.o \
	noumerov.o \
	observer.o \
	potential.o \
	solver.o \
	wave_function.o \
	optimizer.o

bin_PROGRAMS = testlib
testlib_SOURCES = $(module_files) libliee.hpp libliee.cpp testlib.cpp my_util.hpp my_util.cpp optimizer.hpp optimizer.cpp
testlib_LDADD = $(boinc_objects)

#	-lboost_filesystem -lboost_thread -lboost_system -lboost_serialization -lboost_regex
#	-lalglib -lfftw3 -ltinyxml -larchive -lbz2 -lz -llog4cxx -lpthread

libliee.so: $(liee_objects)
	g++ -shared -Wl,-soname,libliee.so -o libliee.so \
	    $(liee_objects) $(boinc_objects) \
	    -lboost_filesystem -lboost_thread -lboost_system -lboost_serialization -lboost_regex \
	    -lalglib -lfftw3 -ltinyxml -larchive -lbz2 -lz -llog4cxx -lpthread

buildlib: libliee.so

#-lcurl