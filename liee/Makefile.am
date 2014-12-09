SUBDIRS = lib

AM_CPPFLAGS = \
    -I/opt/boinc-src/api \
    -I/opt/boinc-src/lib \
    -D"LOG_ENABLED=1" \
    -D"ALGLIB=1"

# remider of include paths for server tools (liee_scheduler, ect.)
# (not needed here)
# -I/opt/boinc-src/sched
# -I/opt/boinc-src/db
# -I/usr/include/mysql

#server_side_files = \
#	optimizer.hpp			optimizer.cpp \
#	work_generator_helper.hpp	work_generator_helper.cpp

module_files = \
	module_config.hpp		module_config.cpp \
	module_factory.hpp \
	potential.hpp			potential.cpp \
	wave_function.hpp		wave_function.cpp \
	noumerov.hpp			noumerov.cpp \
	solver.hpp			solver.cpp \
	observer.hpp			observer.cpp

#alglib_objects = \
#	../alglib/alglibinternal.o \
#	../alglib/alglibmisc.o \
#	../alglib/ap.o \
#	../alglib/dataanalysis.o \
#	../alglib/diffequations.o \
#	../alglib/fasttransforms.o \
#	../alglib/integration.o \
#	../alglib/interpolation.o \
#	../alglib/linalg.o \
#	../alglib/optimization.o \
#	../alglib/solvers.o \
#	../alglib/specialfunctions.o \
#	../alglib/statistics.o


bin_PROGRAMS = liee_worker

liee_worker_SOURCES =  app_opti_wrap.cpp my_util.hpp my_util.cpp optimizer.hpp optimizer.cpp $(module_files)
liee_worker_LDADD = -lalglib -lboinc_api -lfftw3 -lboinc \
	-lboost_filesystem -lboost_thread -lboost_system -lboost_serialization -lboost_regex \
	-ltinyxml -larchive -lbz2 -lz -llog4cxx -lpthread
