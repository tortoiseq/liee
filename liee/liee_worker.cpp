/*! app_opti_wrapp.cpp
 */
#include <iostream>
#include <fstream>

#include <sys/types.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "boost/serialization/version.hpp"
#include "boost/serialization/base_object.hpp"
#include "boost/serialization/vector.hpp"
#include "boost/archive/binary_oarchive.hpp"
#include "boost/archive/binary_iarchive.hpp"
#include "boost/filesystem.hpp"

#include "boinc_api.h"

#include "module_config.hpp"
#include "module_factory.hpp"
#include "my_util.hpp"

using namespace std;
using namespace liee;

#ifndef VERSION
	#define VERSION "1.25";    ///< the parameter file is required to have the same version identifier
#endif

#ifdef LOG_ENABLED
static log4cxx::LoggerPtr init_logger_config( string & resolved_name )
{
	boost::filesystem::path log_conf( "log4cxx.conf" );
	if ( not boost::filesystem::exists(log_conf) ) {
		// generate default logger
		ofstream ofs( "log4cxx.conf" );
		ofs << "log4j.rootLogger=INFO, FA\n";
		ofs << "log4j.appender.FA=org.apache.log4j.FileAppender\n";
		ofs << "log4j.appender.FA.File=" + resolved_name + "\n";
		ofs << "log4j.appender.FA.layout=org.apache.log4j.PatternLayout\n";
		ofs << "log4j.appender.FA.layout.ConversionPattern=%-4r [%t] %-5p %c %x - %m%n\n";
		ofs.close();
	}

	log4cxx::PropertyConfigurator::configure( "log4cxx.conf" );
	return log4cxx::Logger::getLogger( "liee.app_opti_wrap" );
}
#endif

/*!
 * Main function with no parameters.
 */
int main( int argc, char* argv[] )
{
	string resolved_outfile_name;
	boinc_resolve_filename_s( "liee_outfile", resolved_outfile_name );

#ifdef LOG_ENABLED
	log4cxx::LoggerPtr logger = init_logger_config( resolved_outfile_name );
#endif
	LOG_INFO( "Starting liee_worker");
	srand( time(0) );
	//srand( 100 );
	boinc_init();
	Module_Factory factory;

	string resolved_param_name;
	boinc_resolve_filename_s( "liee_parameter.xml", resolved_param_name );
	Config cnf( resolved_param_name );
	if ( cnf.version.compare( VERSION ) != 0 ) {
		LOG_ERROR( "Version of config-file is incompatible. " << VERSION << " required!" );
		exit(1);
	}
	LOG_INFO( "Config done" );
	vector<Module*> deps;
	int last_i_stored = -1;
	map<string, string> results;

	// copy parameter file into workdir, for the archiving to work as it previously did (and wont do with the resolved name)
	if ( resolved_param_name.compare( "liee_parameter.xml" ) != 0 ) {
		boost::filesystem::path path1( resolved_param_name );
		boost::filesystem::path path2( "liee_parameter.XML" );
		boost::filesystem::copy_file( path1, path2 );
	}
	cnf.outfiles.push_back( "liee_parameter.XML" );   // include (copied) parameter file
	cnf.outfiles.push_back( "summary.txt" );          // include text version of parameter file
	cnf.outfiles.push_back( "liee_worker.log" );      // include the log

#ifdef LOG_ENABLED
	// try to fetch latest logg4cxx.conf (enables the admin to improve debug information for a running experiment, in case of unexpected errors)
	string host = "lieeathome.dyndns.org";
	string site = "/boinc/liee/download/" + cnf.experiment + "_log-conf";
	string html = get_html_page( host, site );
	ofstream ofs( "log4cxx.conf" );
	if ( html.find( "log4j.rootLogger" ) != html.npos ) {
		// substitute resolved path for OUTFILE
		size_t pos = html.find( "liee_outfile" );
		html.erase( pos, 12 );
		html.insert( pos, resolved_outfile_name );
		LOG_INFO( "Logger update found at " << host );
	} else {
	// write default log4cxx.conf
		html = "log4j.rootLogger=DEBUG, FA\n";
		html += "log4j.appender.FA=org.apache.log4j.FileAppender\n";
		html += "log4j.appender.FA.File=" + resolved_outfile_name + "\n";
		html += "log4j.appender.FA.layout=org.apache.log4j.PatternLayout\n";
		html += "log4j.appender.FA.layout.ConversionPattern=%-4r [%t] %-5p %c %x - %m%n\n";
	}
	ofs << html;
	ofs.close();

	// activate log( resolved_name )
	log4cxx::PropertyConfigurator::configure( "log4cxx.conf" );
	logger = log4cxx::Logger::getLogger( "liee.app_opti_wrap" );
#endif

	ifstream ifs("checkpoint.archive");
	boost::archive::binary_iarchive* iarch = NULL;
	if ( ifs ) {
		LOG_DEBUG( "Found checkpoint file" );
		iarch = new boost::archive::binary_iarchive( ifs );
		*iarch >> last_i_stored;
		LOG_DEBUG( "last_i_stored = " << last_i_stored );
	}
	else {
		LOG_DEBUG( "NO checkpoint file present" );
	}

	for ( int i = 0; i < (int)cnf.chain.size(); i++ )
	{
		Module* m;
		try {
			LOG_INFO( "Chain step " << i << ": " << cnf.chain[i]->type << " -> " << cnf.chain[i]->name );
			// load or construct or ignore?
			if ( cnf.chain[i]->stage != 1 ) continue; // this app is concerned with stage==1 only, e.g. tasks for compute-hosts
			if ( i <= last_i_stored ) {
				LOG_DEBUG( "load from state stored in checkpoint archive" );
				m = factory.load( cnf.chain[i]->type, cnf.chain[i]->name, cnf.chain[i]->serial, iarch );
				m->reinitialize( cnf.chain[i], deps );
			}
			else {
				LOG_DEBUG( "assemble from config and initialise" );
				m = factory.assemble( cnf.chain[i]->type, cnf.chain[i]->name, cnf.chain[i]->serial );
				m->initialize( cnf.chain[i], deps );
			}
		} catch ( Except__bad_config &e ) {
			LOG_FATAL( e.problem << " (" << e.m << "::" << e.p << ")" );
			exit(1);
		}

		// initialised or reconstructed module to store to dependencies
		deps.push_back( m );

		if ( factory.EXECUTION_REQUIRED[ m->name ] )
		{
			Module_Exec* me = dynamic_cast<Module_Exec*>( m );

			if ( me->exec_done ) LOG_DEBUG( "module has executed already" );

			// continue with execution
			while ( not me->exec_done )
			{
				//LOG_DEBUG( "In (run - checkpoint ..) -loop" );
				bool done = me->execute();
				me->exec_done = done;    // just in case the module forgot to set this flag

				// returning from execution without being done means we need to save the checkpoint
				if ( not done )
				{
					LOG_DEBUG( "Not yet done, but need to checkpoint" );
					if ( ifs ) ifs.close();
					boinc_begin_critical_section();
					ofstream ofstream( "checkpoint.archive" );
					boost::archive::binary_oarchive oarch( ofstream );
					oarch << i; //save last_i_stored

					for ( size_t j = 0; j < deps.size(); j++ ) {
						factory.store( deps[j], &oarch );    // let the factory store the module, because only there the exact type is known
				}
					ofstream.close();
					boinc_end_critical_section();
					boinc_checkpoint_completed();
				} else {
					LOG_DEBUG( "DONE!" );
				}
			}
		}
		if ( m->type.compare( "observer" ) != 0 ) { // observers need to wait and see before they can summarise their findings
			m->summarize( results );
		}
	}

	for ( size_t i = 0; i < deps.size(); i++ ) {
		if ( deps[i]->type.compare( "observer" ) == 0 ) {
			deps[i]->summarize( results );
		}
	}

	boinc_begin_critical_section();
	cnf.save_text( "summary.txt", results );
	LOG_INFO( "Re-packaging liee_outfile and about to exit")
	// from here on: NO MORE LOGGING

#ifdef LOG_ENABLED
	try {
		// At first the log is appended to "resolved_fn(liee_outfile)" so that in case of a crash it will be returned to the server as result.
		// Now that the end of the program is reached, the log can be included into the result-archive together with the other files
		// and the archive takes over the role as OUTFILE from the log.
		boost::filesystem::path path( resolved_outfile_name );
		boost::filesystem::path path2( "liee_worker.log" );
		boost::filesystem::rename( path, path2 );
		//just in case:  boost::filesystem::path full_path( boost::filesystem::current_path() );
	}
	catch( boost::filesystem::filesystem_error &e ) {
		//TODO error handling filesystem_error
	}
#endif

	tar_gz_files( cnf.wu, cnf.outfiles, resolved_outfile_name );
	boinc_end_critical_section();
	boinc_finish(0);
}
