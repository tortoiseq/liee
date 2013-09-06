/*
 * optimization_scheduler.cpp
 *
 *  Created on: 15-Jan-2012
 *      Author: quark
 */

#include <iostream>
#include <fstream>
#include <cstdio>
#include <sstream>

#include <boost/lexical_cast.hpp>
#include <boost/interprocess/sync/file_lock.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/foreach.hpp>
#include <boost/tokenizer.hpp>
#include <boost/filesystem.hpp>

#include "sched_config.h"
#include "sched_util.h"

#include "module_factory.hpp"
#include "module_config.hpp"
#include "downhill_simplex.hpp"
#include "optimizer.hpp"
#include "optimization_scheduler.hpp"

using namespace std;
using namespace liee;


Opti_Scheduler::Opti_Scheduler( string exp, int app_id )
{
	experiment = exp;
	state_file_name = "/var/www/boinc/liee/workspace/" + experiment + "_state";
	appid = app_id;
	opti = NULL;

	//some defaults
	priority = 0;
	delay_bound = 86400;
	min_quorum = 1;
	target_nresults = 1;
	max_error_results = 4;
	max_total_results = 8;
	max_success_results = 4;
	rsc_bandwidth_bound = 0;

	// fetch user configured boinc parameters from config
	string conf_filename = "/var/www/boinc/liee/workspace/" + experiment + "_conf.xml";  //TODO get project root from boinc config
	Config cnf( conf_filename );
	for ( size_t i = 0; i < cnf.chain.size(); i++ ) {
		if ( cnf.chain[i]->type.compare( "scheduler" ) == 0  &&  cnf.chain[i]->name.compare( "boinc" ) == 0  )
		{
			Conf_Module* boinc_param = cnf.chain[i];
			priority = (int)boinc_param->getParam("priority")->value;
			delay_bound = (int)boinc_param->getParam("delay_bound")->value;
			min_quorum = (int)boinc_param->getParam("min_quorum")->value;
			target_nresults = (int)boinc_param->getParam("target_nresults")->value;
			max_error_results = (int)boinc_param->getParam("max_error_results")->value;
			max_total_results = (int)boinc_param->getParam("max_total_results")->value;
			max_success_results = (int)boinc_param->getParam("max_success_results")->value;
			rsc_bandwidth_bound = (int)boinc_param->getParam("rsc_bandwidth_bound")->value;
		}
	}
}

Opti_Scheduler::~Opti_Scheduler() {}

void Opti_Scheduler::load()
{
	try // to obtain a lock on the state-file
	{
		boost::interprocess::file_lock flock( state_file_name.c_str() );
		// ^- This throws if the file does not exist or it can't open it with read-write access!
		if ( flock.try_lock() == false ) {
			cout << "\t try_lock failed \n";
			return;  //TODO error handling
		}

		if ( opti != NULL ) delete opti;
		// read optimiser state using serialisation
		ifstream ifs( state_file_name.c_str() );  //TODO handle more IO errors
		boost::archive::binary_iarchive iarch( ifs );

		{ // issues with abstract types and serialisation
			string type_name;
			iarch >> type_name;

			if ( type_name.compare( "down_simplex" ) == 0 ) {
				liee::opti::Downhill_Simplex* o;
				cout << "\t load Downhill_Simplex \n";
				iarch >> o;
				opti = o;
			}
			else if ( type_name.compare( "particle_swarm" ) == 0 ) {
				liee::opti::Particle_Swarm_Optimizer* o;
				cout << "\t load Particle_Swarm_Optimizer \n";
				iarch >> o;
				opti = o;
			}
			else if ( type_name.compare( "shot_gun" ) == 0 ) {
				liee::opti::Shot_Gun_Optimizer* o;
				cout << "\t load Shot_Gun_Optimizer \n";
				iarch >> o;
				opti = o;
			}
			else if ( type_name.compare( "rasterizer" ) == 0 ) {
				liee::opti::Rasterizer* o;
				cout << "\t load Rasterizer \n";
				iarch >> o;
				opti = o;
			}
			else {
				log_messages.printf( MSG_CRITICAL, "Unknown type of optimizer, don't know how to load! \n" );
				opti = NULL;
			}
		}

		ifs.close();
		flock.unlock();
	}
	catch ( boost::interprocess::interprocess_exception &ex ) {
		cout << ex.what() << endl;
		//TODO error handling
	}
}

void Opti_Scheduler::store( bool do_lock )
{
	if ( opti == NULL ) return; //TODO error handling

	try // to obtain a lock on the state-file
	{
		boost::interprocess::file_lock * flock;
		if ( do_lock ) {
			flock = new boost::interprocess::file_lock( state_file_name.c_str() );
			// ^- This throws if the file does not exist or it can't open it with read-write access!
			if ( flock->try_lock() == false ) {
				return;	//TODO error handling
			}
		}

		// write back changed optimiser state
		ofstream ofs( state_file_name.c_str() );
		boost::archive::binary_oarchive oarch( ofs );

		{ // issues with abstract types and serialisation
			oarch << opti->type_name;

			if ( opti->type_name.compare( "down_simplex" ) == 0 ) {
				liee::opti::Downhill_Simplex* o = dynamic_cast<liee::opti::Downhill_Simplex*>( opti );
				cout << "\t write as Downhill_Simplex \n";
				oarch << o;
			}
			else if ( opti->type_name.compare( "particle_swarm" ) == 0 ) {
				liee::opti::Particle_Swarm_Optimizer* o = dynamic_cast<liee::opti::Particle_Swarm_Optimizer*>( opti );
				cout << "\t write as Particle_Swarm_Optimizer \n";
				oarch << o;
			}
			else if ( opti->type_name.compare( "shot_gun" ) == 0 ) {
				liee::opti::Shot_Gun_Optimizer* o = dynamic_cast<liee::opti::Shot_Gun_Optimizer*>( opti );
				cout << "\t write as Shot_Gun_Optimizer \n";
				oarch << o;
			}
			else if ( opti->type_name.compare( "rasterizer" ) == 0 ) {
				liee::opti::Rasterizer* o = dynamic_cast<liee::opti::Rasterizer*>( opti );
				cout << "\t write as Rasterizer \n";
				oarch << o;
			}
			else {
				log_messages.printf( MSG_CRITICAL, "Unknown type of optimizer, don't know how to store! \n" );
			}
		}
		ofs.close();

		// release lock
		if ( do_lock ) {
			flock->unlock();
			delete flock;
		}
	}
	catch ( boost::interprocess::interprocess_exception &ex ) {
		cout << ex.what() << endl;
		return;  //TODO error handling
	}
}

void Opti_Scheduler::write_template( string & filename, vector<string> & infiles )
{
}

void Opti_Scheduler::initialize( int & flag )
{
	ifstream test_exist( state_file_name.c_str() );
	if ( test_exist ) {
		cout << "Error: The state-file " << state_file_name << " already exists! Won't overwrite.\n";
		flag = 1;
		return;
	}

	string conf_filename = "/var/www/boinc/liee/workspace/" + experiment + "_conf.xml";  //TODO get project root from boinc config
	Config cnf( conf_filename );
	infiles = cnf.infiles;

	// place additional infiles (other than the parameter file) into the appropriate download folder.
	for (size_t i = 1; i < cnf.infiles.size(); i++) {
		char path[512];
		dir_hier_path( cnf.infiles[i].c_str(), "/var/www/boinc/liee/download", 1024, path, true ); //TODO get uldl_dir_fanout from boinc config, here default value 1024 assumed, also get project root from config
		string dest_str( path );
		string srce_str( "/var/www/boinc/liee/workspace/" );
		srce_str += cnf.infiles[i];
		boost::filesystem::path dest( dest_str );
		boost::filesystem::path srce( srce_str );
		if ( not boost::filesystem::exists( dest ) ) {
			cout << cnf.infiles[i] << " not found in download path! (" << path << "). Trying to copy it from workspace...\n";
			boost::filesystem::copy_file( srce, dest );  //TODO error handling
		} else {
			cout << cnf.infiles[i] << " already in download path! (" << path << ")\n";
		}
	}

	{ // generate IN-template definition (to tell boinc about additional infiles)
		ofstream ofs( ("/var/www/boinc/liee/workspace/" + experiment + "_in").c_str() );
		ofs << "<input_template>\n";
		for (size_t i = 0; i < cnf.infiles.size(); i++) {
			ofs << "\t<file_info>\n";
			ofs << "\t\t<number>" << boost::lexical_cast<string>( i ) << "</number>\n";
			if ( i > 0 ) {
				ofs << "\t\t<sticky/> <no_delete/>\n";  // set all dependency files to sticky-no_delete, because they are probably shared among workunits
			}
			ofs << "\t</file_info>\n";
		}
		ofs << "\t<workunit>\n";
		for (size_t i = 0; i < cnf.infiles.size(); i++)	{
			ofs << "\t\t<file_ref>\n";
			ofs << "\t\t\t<file_number>" << boost::lexical_cast<string>( i ) << "</file_number>\n";
			ofs << "\t\t\t<open_name>" << cnf.infiles[i] << "</open_name>\n";
			ofs << "\t\t</file_ref>\n";
		}
		ofs << "\t\t<min_quorum>1</min_quorum>\n";  //TODO make configurable
		ofs << "\t\t<target_nresults>1</target_nresults>\n";  //TODO make configurable
		ofs << "\t</workunit>\n";
		ofs << "</input_template>\n";
		ofs.close();
	}

	{ // generate static OUT-template definition (just for symmetry and future variability of those parameters )
		ofstream ofs( ("/var/www/boinc/liee/workspace/" + experiment + "_out").c_str() );
		ofs << "<output_template>\n";
		ofs << "\t<file_info>\n";
		ofs << "\t\t<name><OUTFILE_0/></name>\n";
		ofs << "\t\t<generated_locally/>\n";
		ofs << "\t\t<upload_when_present/>\n";
		ofs << "\t\t<max_nbytes>1000000000</max_nbytes>\n";
		ofs << "\t\t<url><UPLOAD_URL/></url>\n";
		ofs << "\t</file_info>\n";

		ofs << "\t<result>\n";
		ofs << "\t\t<file_ref>\n";
		ofs << "\t\t<file_name><OUTFILE_0/></file_name>\n";
		ofs << "\t\t\t<open_name>liee_outfile</open_name>\n";
		ofs << "\t\t</file_ref>\n";
		ofs << "\t\t<report_immediately/>\n";
		ofs << "\t</result>\n";

		ofs << "</output_template>\n";
		ofs.close();
	}

	//fetch optimiser-config
	Conf_Module* oc = NULL;
	for ( size_t i = 0; i < cnf.chain.size(); i++ ) {
		if ( cnf.chain[i]->type.compare( "optimizer" ) == 0 ) {
			oc = cnf.chain[i];
			cout << " using Optimizer: " << oc->name << "\n";
			break;
		}
	}
	if ( oc == NULL ) {
		cout << "Error: No Optimizer specified in config's exec_chain. Exit \n";
		exit(1);
	}

	// Optimizer-Factory
	if ( opti != NULL ) delete opti;
	if ( oc->name.compare("down_simplex") == 0 ) {
		liee::opti::Downhill_Simplex* o = new liee::opti::Downhill_Simplex();
		opti = o;
	}
	else if ( oc->name.compare("particle_swarm") == 0 ) {
		int popul = (int) oc->getParam("population")->value;
		liee::opti::Particle_Swarm_Optimizer* o = new liee::opti::Particle_Swarm_Optimizer( popul );
		o->hood_sz = (int) oc->getParam("neighbours")->value;
		opti = o;
	}
	if ( oc->name.compare("shot_gun") == 0 ) {
		int bullets = (int) oc->getParam("num_bullets")->value;
		liee::opti::Shot_Gun_Optimizer* o = new liee::opti::Shot_Gun_Optimizer( );
		o->num_bullets_per_shot = bullets;
		opti = o;
	}
	if ( oc->name.compare("rasterizer") == 0 ) {
		liee::opti::Rasterizer* o;
		if ( oc->getParam("vector_num_samples")->values.size() == 0 ) {
			int unif = (int) oc->getParam("uniform_num_samples")->value;
			o = new liee::opti::Rasterizer( unif );
		}
		else {
			o = new liee::opti::Rasterizer( oc->getParam("vector_num_samples")->values );
		}
		opti = o;
	}
	else {}  //TODO error handling: unsupported optimizer

	opti->type_name = oc->name;
	opti->max_eval = (int)oc->getParam("max_eval")->value;
	opti->tolerance = oc->getParam("conv_tol")->value;
	opti->initialise( cnf.chain_params_merged );
	cout << "Initialization done.\n";
	store( false );
	flag  = 0;
}

void Opti_Scheduler::make_job( const liee::opti::Request & r, DB_WORKUNIT & job, Config & templ, int & flag )
{
	// make a unique name (for the job and its input file)
	char path[512];
	stringstream ss_name;
	ss_name << experiment;
	if ( opti->type_name.compare( "rasterizer" ) == 0 ) {
		for ( size_t i = 0; i < r.x.size(); i++ ) {
			ss_name << "_" << ( (int) r.x[i] );
		}
	}
	else {
		ss_name << "_" << r.id;
	}
	string name = ss_name.str();
	templ.wu = name;  // write wu-name into the config, so that the host-process can build the results-archive with a tailored path
	// Create the input file.
	// Put it at the right place in the download dir hierarchy
	flag = config.download_path( name.c_str(), path );
	if ( flag ) return;	// TODO error handling

	// insert requested values to template
	//TODO move this functionality into class Config, since we fiddle with its data members
	size_t param_index = 0;
	BOOST_FOREACH( Conf_Param* p, templ.chain_params_merged ) {
		if ( not p->fixed )
		{
			if ( param_index == r.x.size() ) {
				log_messages.printf( MSG_CRITICAL, "Wrong search space dimension: not enough variable parameters in template. \n" );
				exit( 2 );
			}
			cout << "\t change variable parameter: " << p->name << " <-- " << r.x_new[param_index] << "\n";
			p->value = r.x_new[param_index++];
			p->text = ""; // clear p->text in order to force usage of p->value
		}
	}
	templ.reevaluate_expressions();

	// let deployed modules estimate computational effort
	double flops = 0;
	double ram = 0;
	double disk = 0;
	vector<Module*> deps;
	BOOST_FOREACH( Conf_Module* cm, templ.chain ) {
		if ( cm->stage != 1 ) continue;  // the app is concerned with stage-1 only, e.g. tasks for compute-hosts
		Module* m = Module_Factory::assemble( cm->type, cm->name, cm->serial );
		m->estimate_effort( cm, flops, ram, disk );
	}
	//cout << "flops/ram/disk:" << flops << "\t" << ram << "\t" << disk << "\n";

	// write work unit
	string path_ = string( path );
	templ.save_file( path_ );

	// Fill in the job parameters
	job.clear();
	job.appid = appid;
	strcpy( job.name, name.c_str() );
	//TODO the estimated flops/ram/disk are currently much too low. --> multiply by 100 for estimate and 1000 for bound
	job.rsc_fpops_est = 100.0 * flops;
	job.rsc_fpops_bound = 1000.0 * flops;
	job.rsc_memory_bound = 10.0 * ram + 1e6;
	job.rsc_disk_bound = 10.0 * disk + 1e6;
	job.delay_bound = delay_bound;
	job.min_quorum = min_quorum;
	job.target_nresults = target_nresults;
	job.max_error_results = max_error_results;
	job.max_total_results = max_total_results;
	job.max_success_results = max_success_results;
	job.rsc_bandwidth_bound = rsc_bandwidth_bound;
}

void Opti_Scheduler::make_jobs( vector<DB_WORKUNIT> & jobs, int & flag )
{
	load();

	if ( opti == NULL ) {
		log_messages.printf( MSG_CRITICAL, "Couldn't restore optimizer from statefile \n" );
		exit( 1 );
	}

	vector<liee::opti::Request> work;
	opti->generate_requests( work );
	cout << "\t generate_requests() did succeed with " << work.size() << " new jobs \n";

	store( true );

	if ( work.size() == 0 ) return;

	//TODO instead of parsing the config in every loop, keep some state information in memory
	string conf_filename = "/var/www/boinc/liee/workspace/" + experiment + "_conf.xml";
	Config conf_template( conf_filename );
	infiles = conf_template.infiles;

	// process jobs: write input-files and create work units for registration
	for( size_t i = 0; i < work.size(); i++ ) {
		DB_WORKUNIT job;
		make_job( work[i], job, conf_template, flag );
		if ( flag ) return;  //TODO error handling
		jobs.push_back( job );
	}
	flag = 0;
}

int Opti_Scheduler::assimilate( int & flag )
{
	using boost::lexical_cast;
	using boost::bad_lexical_cast;
	int count;

	try // to obtain a lock on the result table
	{
		string table_name = "/var/www/boinc/liee/workspace/" + experiment + ".results";
		boost::interprocess::file_lock flock( table_name.c_str() );
		// ^- This throws if the file does not exist or it can't open it with read-write access!
		flock.lock();

		ifstream ifs_results( table_name.c_str() );
		string line;
		if ( ifs_results.is_open() )
		{
			vector<liee::opti::Request> v_res;
			// read whole table
			while ( ! ifs_results.eof() )
			{
				liee::opti::Request r;
				getline( ifs_results, line );
				size_t tab = line.find_first_of( "\t", 0 );
				if ( line.size() > 2  &&  tab != line.npos ) {
					r.id = lexical_cast<int>( line.substr( 0, tab ) );
					r.y = lexical_cast<double>( line.substr( tab + 1, strlen( line.c_str() ) ) );
					r.flag = 0;
					v_res.push_back( r );
				}
			}
			ifs_results.close();
			// empty table
			FILE * f = fopen( table_name.c_str(), "w" );
			fclose( f );
			flock.unlock();
			count = v_res.size();

			if ( v_res.size() > 0 )
			{
				load();
				if ( opti == NULL ) {
					log_messages.printf( MSG_CRITICAL, "Couldn't restore optimizer from statefile \n" );
					flag = 3;
					return 0;
				}

				log_messages.printf( MSG_DEBUG, "Sending %d results to the optimizer. \n", (int) v_res.size() );
				opti->assimilate_results( v_res );

				store( true );
			}
			else {
				cout << "\t nothing to assimilate \n";
			}
		}
		else {
			log_messages.printf( MSG_CRITICAL, "Can't open %s \n", table_name.c_str() );
			flag = 1;
			return 0;
		}
	}
	catch ( boost::interprocess::interprocess_exception &ex ) {
		cout << ex.what() << endl;
		flag = 2;
		return 0;
	}
	flag = 0;
	return count;
}

void Opti_Scheduler::get_state_of_affairs( int & evals, vector<double> & x, double & y, double & conv, int & flag )
{
	load();
	evals = opti->evaluations;
	x = opti->global_min_pos;
	y = opti->global_min;
	conv = opti->histo_var;
	flag = 0;
}
