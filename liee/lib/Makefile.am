# recompiles source here with -fPIC (position independent code for shared library)
# and deactivated logging (-D"LOG_DISABLED")

AM_CPPFLAGS = \
    -I.. \
    -I/opt/boinc-src/lib \
    -I/opt/boinc-src/api \
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

bin_PROGRAMS = testlib
testlib_SOURCES = testlib.cpp my_util.hpp my_util.cpp optimizer.hpp optimizer.cpp $(module_files) libliee.hpp libliee.cpp
testlib_LDADD = -lalglib -lfftw3 -lboinc_api -lboinc \
	-lboost_filesystem -lboost_thread -lboost_system -lboost_serialization -lboost_regex \
	-lcurl -ltinyxml -larchive -lbz2 -lz -llog4cxx -lpthread

libliee.so: libliee.o module_config.o my_util.o noumerov.o observer.o optimizer.o potential.o solver.o wave_function.o
	g++ -shared -Wl,-soname,libliee.so -o libliee.so libliee.o my_util.o module_config.o noumerov.o observer.o potential.o solver.o wave_function.o optimizer.o -lalglib -lfftw3 -lboinc_api -lboinc -lboost_filesystem -lboost_thread -lboost_system -lboost_serialization -lboost_regex -lcurl -ltinyxml -larchive -lbz2 -lz -llog4cxx -lpthread

buildlib: libliee.so
