// This file is based on a BOINC example.

// sample_work_generator.cpp: an example BOINC work generator.
// This work generator has the following properties
// (you may need to change some or all of these):
//
// - Runs as a daemon, and creates an unbounded supply of work.
//   It attempts to maintain a "cushion" of 100 unsent job instances.
//   (your app may not work this way; e.g. you might create work in batches)
// - Creates work for the application "uppercase".
// - Creates a new input file for each job;
//   the file (and the workunit names) contain a timestamp
//   and sequence number, so that they're unique.

#include <unistd.h>
#include <cstdlib>
#include <string>
#include <cstring>
#include <vector>
#include <iostream>
#include <fstream>

#include <boost/lexical_cast.hpp>

#include "boinc_db.h"
#include "error_numbers.h"
#include "backend_lib.h"
#include "parse.h"
#include "util.h"
#include "svn_version.h"

#include "sched_config.h"
#include "sched_util.h"
#include "sched_msgs.h"
#include "str_util.h"

#include "optimization_scheduler.hpp"
#include "optimizer.hpp"

// globals
char* wu_template;			// contents, not path
const char* result_template_filename;		// relative to project root
const char* result_template_filepath;		// absolute or relative to current dir
//
//DB_APP app;
char * experiment;
int APP_ID = 1; //TODO make configurable!



// create one new job
//
void make_job( DB_WORKUNIT & job, vector<string> & dep_files, int & flag ) 
{
    //TODO try using the API call instead of system call
    stringstream command;
    command << "/var/www/boinc/liee/bin/create_work --appname liee_worker --wu_name ";
    command << job.name; //wu_name
    command << " --wu_template ./workspace/" << experiment << "_in";
    command << " --result_template workspace/" << experiment << "_out";		// does not accept absolute paths
    command << " --priority " << job.priority;
    command << " --rsc_fpops_est " << job.rsc_fpops_est;
    command << " --rsc_fpops_bound " << job.rsc_fpops_bound;
    command << " --rsc_memory_bound " << job.rsc_memory_bound;
    command << " --rsc_disk_bound " << job.rsc_disk_bound;
    command << " --delay_bound " << job.delay_bound;
    command << " --min_quorum " << job.min_quorum;
    command << " --target_nresults " << job.target_nresults;
    command << " --max_error_results " << job.max_error_results;
    command << " --max_total_results " << job.max_total_results;
    command << " --max_success_results " << job.max_success_results;

    const char* infiles[ dep_files.size() ];
    for ( size_t i = 0; i < dep_files.size(); i++ ) {
	// char name[256];
	if ( i == 0 ) {
	    //strcpy( name, job.name );
	    command << " " << job.name;
	} else {
	    //strcpy( name, dep_files[i].c_str() );
	    command << " " << dep_files[i];
	}
	//infiles[i] = name; 
    }
    // Register the job with BOINC
    //flag = create_work( job, wu_template, result_template_filename, result_template_filepath, infiles, dep_files.size(), config );
    
    system( command.str().c_str() );
    //cout << command.str().c_str() << "\n";
}

void main_loop() 
{
    cout << "...in main_loop() \n";
    Opti_Scheduler sched( experiment, APP_ID ); 
    vector<DB_WORKUNIT> jobs;
    int flag;
    bool no_progress = true;

    //TODO timeouts adjustable at commandline
    while ( true ) 
    {
        sleep( 200 );
        //check_stop_daemons();
        
        // see if optimiser has new requests vector<Job_Record> & jobs
        sched.make_jobs( jobs, flag );
        
        if ( flag ) {
            log_messages.printf( MSG_DEBUG, "Work generation blocked or failed: %d\n", flag );
        } 
        else if ( jobs.size() == 0 ) {
            //log_messages.printf( MSG_DEBUG, "No more work available... \n" );
        } 
        else {
            log_messages.printf( MSG_DEBUG, "Making %d jobs\n", (int)jobs.size() );
        
            for ( size_t i = 0; i < jobs.size(); i++ ) 
            {
                make_job( jobs[i], sched.infiles, flag );
                if ( flag ) {
                    log_messages.printf( MSG_CRITICAL, "Can't make job (exiting): %d\n", flag );
                    exit( flag );
                } 
                else no_progress = false;
            }
            jobs.clear();
        }

        // check result table and return them to the optimizer
        int count = sched.assimilate( flag );
        
        // sleep some if both of the last two actions failed, otherwise hurry on
        if ( no_progress && (count == 0) ) {
            sleep( 1500 );
        }
    }
}

void usage(char *name) {
    fprintf(stderr, "This is an example BOINC work generator.\n"
        "This work generator has the following properties:\n"
        "\n"
        "- Creates work for the application example_app1\n"
        "- Runs as a daemon, and creates workunits whenever the optimizer\n"
        "  requests the evaluation of new points.\n"
        "\n"
        "Prerequisits:\n"
        "- The directory /var/www/boinc/liee/workspace must exist and be writable.\n"
        "- The optimizer must be initialized first using the -i command line\n"
        "  option!\n"
	"- A configuration file compatible to liee::Config has to be present\n"
	"  at location /var/www/boinc/liee/workspace/$exp_conf.xml\n"
        "- The complementary assimilator has to be running, which feeds\n"
        "  the results back to the optimizer, so that the optimizer\n"
        "  is able to make progress in its algorithm and request new work.\n"
        "\n"
        "- Creates a new input file for each job;\n"
        "  the file (and the workunit names) contain the experiment name\n"
        "  and a serial number.\n"
        "  Make sure to choose unique experiment names to avoid filename \n"
        "  collisions!\n"
        "\n"
        "Usage: %s [OPTION]...\n\n"
        "Options:\n"
        "  [ -i | --initialize ]       Creates a new optimizer with default\n"
        "                                parameters and its state-file.\n"
        "  [ -e X ]                    Sets the unique experiment name to X.\n"
        "  [ -d X ]                    Sets debug level to X.\n"
        "  [ -x ]                      Prints out: number of evaluations; so far\n"
        "                                best parameter-vector, its function value\n"
        "                                and convergence info. \n"
        "  [ -l ]                      Writes every 5 minutes the '-x' output to\n"
        "                                the logfile /var/www/boinc/liee/workspace/$exp.log\n"
        "  [ -h | --help ]             Shows this help text.\n"
        "  [ -v | --version ]          Shows version information.\n"
        "\n",
        name
    );
}

int print_status_info( ostream & os )
{
    int retval;
    Opti_Scheduler sched( experiment, APP_ID );
    vector<double> x;
    double y, conv;
    int n;
    sched.get_state_of_affairs( n, x, y, conv, retval );
    if ( retval ) return retval;
    
    os << n << "\t" << y << "\t" << conv << "\t" << "{";
    if ( x.size() > 0 ) os << x[0]; 
    for ( size_t i = 1; i < x.size(); i++ ) {
        os << "," << x[i];
    }
    os << "}" << endl;
    return retval;
}

int main(int argc, char** argv) {
    int i, retval;
    bool do_ini = false;
    bool do_log = false;
    bool print_info = false;

    for (i=1; i<argc; i++) {
        if (is_arg(argv[i], "e")) {
            if (!argv[++i]) {
                log_messages.printf(MSG_CRITICAL, "%s requires an argument\n\n", argv[--i]);
                usage(argv[0]);
                exit(1);
            }
            experiment = argv[i];
            
        } 
        else if (is_arg(argv[i], "d")) {
            if (!argv[++i]) {
                log_messages.printf(MSG_CRITICAL, "%s requires an argument\n\n", argv[--i]);
                usage(argv[0]);
                exit(1);
            }
            int dl = atoi(argv[i]);
            log_messages.set_debug_level(dl);
            if (dl == 4) g_print_queries = true;
        } 
        else if (is_arg(argv[i], "h") || is_arg(argv[i], "help")) {
            usage(argv[0]);
            exit(0);
        } 
        else if (is_arg(argv[i], "v") || is_arg(argv[i], "version")) {
            printf("%s\n", SVN_VERSION);
            exit(0);
        // changes
        } 
        else if (is_arg(argv[i], "i") || is_arg(argv[i], "initialize")) {
            do_ini = true;
        } 
        else if (is_arg(argv[i], "l") ) {
            do_log = true;
        } 
        else if (is_arg(argv[i], "x") ) {
            print_info = true;
        } 
        else {
            log_messages.printf(MSG_CRITICAL, "unknown command line argument: %s\n\n", argv[i]);
            usage(argv[0]);
            exit(1);
        }
    }
    
    if ( do_ini ) 
    {
        cout << "initializing " << experiment << "\n";
        Opti_Scheduler sched( experiment, APP_ID );
        sched.initialize( retval );
        // create empty table- and log file
        string exp = experiment;
        string table_name = "/var/www/boinc/liee/workspace/" + exp + ".results";
        FILE * f = fopen( table_name.c_str(), "w" );
        fclose( f );
        string log_name = "/var/www/boinc/liee/workspace/" + exp + ".log";
        f = fopen( log_name.c_str(), "w" );
        fclose( f );
        
        exit( retval );
    }

    if ( print_info ) {
        retval = print_status_info( cout );
        exit( retval );
    }
    
    if ( do_log ) {
        string exp = experiment;
        string logfile = "/var/www/boinc/liee/workspace/" + exp + ".log";
        ofstream ofs; 
        ofs.open( logfile.c_str(), fstream::app );
        while ( true ) 
        {
            check_stop_daemons();
            print_status_info( ofs );
            sleep( 300 );
        }
    }

    retval = config.parse_file();
    if ( retval ) {
        log_messages.printf(MSG_CRITICAL,
            "Can't parse config.xml: %s\n", boincerror(retval)
        );
        exit( 1 );
    }

    
    result_template_filename = ("workspace/" + string( experiment ) + "_out").c_str();
    result_template_filepath = config.project_path( result_template_filename );
    const char* wu_template_filename = config.project_path ( ("/workspace/" + string( experiment ) + "_in").c_str() );
    if ( read_file_malloc( wu_template_filename, wu_template) ) {      
        log_messages.printf(MSG_CRITICAL, "Can't read WU template\n");
        exit( 1 );
    }    

    retval = boinc_db.open(
        config.db_name, config.db_host, config.db_user, config.db_passwd
    );
    if ( retval ) {
        log_messages.printf(MSG_CRITICAL, "Can't open db\n");
        exit( 1 );
    }

    log_messages.printf(MSG_NORMAL, "Starting\n");
    
    main_loop();
}
