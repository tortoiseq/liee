/*! optimization_scheduler.hpp
 */

#ifndef WORK_GENERATOR_HELPER_HPP_
#define WORK_GENERATOR_HELPER_HPP_

#include <string>
#include <vector>

#include "boinc_db.h"

#include "module_config.hpp"
#include "optimizer.hpp"

using namespace std;

/*!
 *
 */
class Opti_Scheduler {
public:
	liee::opti::Asynch_Optimizer* opti;
	string experiment;
	string state_file_name;
	vector<string> infiles;
	int appid;

	int priority;
	int delay_bound;
	int min_quorum;
	int target_nresults;
	int max_error_results;
	int max_total_results;
	int max_success_results;
	int rsc_bandwidth_bound;

	Opti_Scheduler( string exp, int app_id );
	~Opti_Scheduler();

	/*!
	 * Sets up a fresh optimizer and writes its state-file to "workspace/experiment".
	 * Prerequisite: A configuration file compatible to class liee::Config has to be
	 * present at location "workspace/experiment_conf.xml" (experiment replaced by the
	 * actual name).
	 */
	void initialize( int & flag );

	/*!
	 * @see sample_work_generator--make_job()
	 */
	void make_jobs( vector<DB_WORKUNIT> & jobs, int & flag );

	/*!
	 * Reads the text-file  /var/boinc/workspace/$(experiment).results which is
	 * required to contain ($id - $value) pairs on each row as the results of
	 * previously requested evaluations.
	 * Those results are returned to the underlying optimization routine.
	 * @param flag	values other than zero indicate an error.
	 * @return the number of assimilated results
	 */
	int assimilate( int & flag );

	void get_state_of_affairs( int & evals, vector<double> & x, double & y, double & conv, int & flag );

private:
	void load();
	void store( bool do_lock );
	void make_job( const liee::opti::Request & x, DB_WORKUNIT & job, liee::Config & templ, int & flag );
	void write_template( string & filename, vector<string> & infiles );
};


#endif /* WORK_GENERATOR_HELPER_HPP_ */
