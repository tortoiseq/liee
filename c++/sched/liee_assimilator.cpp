// This file is part of BOINC.
// http://boinc.berkeley.edu
// Copyright (C) 2008 University of California
//
// BOINC is free software; you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License
// as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// BOINC is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with BOINC.  If not, see <http://www.gnu.org/licenses/>.

// A sample assimilator that:
// 1) if success, copy the output file(s) to a directory
// 2) if failure, append a message to an error log

#include <limits>
#include <vector>
#include <string>
#include <sstream>
#include <cstdlib>

#include <iostream>
#include <fstream>
#include <iomanip>

#include <boost/lexical_cast.hpp>
#include <boost/interprocess/sync/file_lock.hpp>
#include <boost/filesystem.hpp>

#include "boinc_db.h"
#include "error_numbers.h"
#include "filesys.h"
#include "sched_msgs.h"
#include "validate_util.h"
#include "sched_config.h"

using namespace std;

int write_error(char* p) {
    static FILE* f = 0;
    if (!f) {
        f = fopen(config.project_path("results/errors"), "a");
        if (!f) return ERR_FOPEN;
    }
    fprintf(f, "%s", p);
    fflush(f);
    return 0;
}

int assimilate_handler( WORKUNIT& wu, vector<RESULT>& /*results*/, RESULT& canonical_result ) 
{
    int retval;
    char buf[1024];

    //retval = boinc_mkdir( config.project_path( "results" ) );
    if ( retval ) return retval;

    if ( wu.canonical_resultid ) 
    {
        vector<OUTPUT_FILE_INFO> output_files;
        get_output_file_infos(canonical_result, output_files);

        // we only process one output_file
        {
            OUTPUT_FILE_INFO& fi = output_files[0];
            const char* result_dir = config.project_path("results");
            double y = 666;
	    
	    // test whether output_file is only the log (indicates error)
	    ifstream ifs( fi.path.c_str(), ifstream::in );
	    string line = "";
	    getline( ifs, line );
	    if ( line.find( "Starting liee_worker" ) != line.npos ) {
		// output_file is just the log. -> did not finish -> invalid result
		// TODO this needs to be checked in the validator, not here
		log_messages.printf( MSG_DEBUG, "found the first log-line in the outfile -> error-log\n" );
		
		// save log in results directory
		int i = 0;
		while ( i < 100 ) {
		    stringstream dest_name;
		    dest_name << result_dir << "/" << wu.name << "_error-log." << i;
		    boost::filesystem::path dest( dest_name.str() );
		    if ( not boost::filesystem::exists( dest ) ) {
			ofstream ofs( dest_name.str().c_str() );
			ofs << line << "\n";
			while ( getline( ifs, line ) ) {
			    ofs << line << "\n";
			}
			ofs.close();
			break;
		    } 
		    i++;
		}
		ifs.close();
		return 0;
	    }
	    ifs.close();
	    
            // extract the output files to the their directory
            stringstream tar_cmd;
            tar_cmd << "tar -C " << result_dir << " -xzf " << fi.path;
            int ret = system( tar_cmd.str().c_str() );
            
            stringstream result_file;
            result_file << result_dir << "/" << wu.name << "/summary.txt";
            line = "";
            char * endp;
            ifstream ifs_result( result_file.str().c_str() );
            string seek = "objective\t";
            while ( ifs_result.is_open() && ifs_result.good() )
            {
                getline( ifs_result, line );
                size_t pos = line.find( seek );
                if ( pos != line.npos ) { //found the right line
                    string objective = line.substr( pos + seek.length() );
		    try {
			 y = boost::lexical_cast<double>( objective );
		    } 
		    catch ( boost::bad_lexical_cast & e ) {
			log_messages.printf( MSG_DEBUG, "bad lexical cast when reading objective\n" );
		    }
                    break;
                }
            }
            ifs_result.close();
                
            // The first part of wu.name defines the experiment and therefore the results file.
            string name, path;
            name = wu.name;
            int spos = name.find_first_of( "_", 0 );
            path = "/var/www/boinc/liee/workspace/" + name.substr( 0, spos ) + ".results";
                
            boost::interprocess::file_lock * flock;
            try // to obtain a lock on the result table
            {
                flock = new boost::interprocess::file_lock( path.c_str() );
                // ^- This throws if the file does not exist or it can't open it with read-write access!
                flock->lock();
            }
            catch ( boost::interprocess::interprocess_exception &ex ) {
                cout << ex.what() << endl;
                return 666; //TODO error handling 
            }
            log_messages.printf( MSG_DEBUG, "return value = %g\n", y );
                
            ofstream ofs_table; 
            ofs_table.open( path.c_str(), fstream::app );
            ofs_table << setprecision( numeric_limits<double>::digits10 + 2 );
            ofs_table << name.substr( spos + 1, strlen( name.c_str() ) );   // the serial number of the wu
            ofs_table << "\t" << y << "\n";
            ofs_table.close();
             
            flock->unlock();
            delete flock;
            }
            
    } else {
        sprintf(buf, "%s: 0x%x\n", wu.name, wu.error_mask);
        return write_error( buf );
    }
    return 0;
}
