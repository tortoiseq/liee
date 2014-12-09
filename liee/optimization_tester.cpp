/*
 * optimization_tester.cpp
 *
 *  Created on: 24-Jan-2012
 *      Author: quark
 */

#include <string>
#include <vector>

//#include <filesys.h>
//#include <fstream>
//#include "boost/lexical_cast.hpp"
//#include "boinc_api.h"

#include "my_util.hpp"
#include "potential.hpp"

#include "optimizer.hpp"
#include "downhill_simplex.hpp"
#include "objective_functions.hpp"

using namespace std;
using namespace liee;


/*!
 * Main function with (currently) no parameters.
 */
#ifdef SERVER_BUILD
int main( int argc, char* argv[] )
#else
int main__( int argc, char* argv[] )
#endif
{
	//cout << "Testing Particle_Swarm_Optimizer... \n";
	cout << "Testing Downhill_Simplex... \n";
	cout.precision( 6 );

	int dim = 10;
	int N = 25;
	int MAX_ITER = 100000;

	vector<double> lo;
	vector<double> hi;
	for ( int i = 0; i < dim; i++ ) {
		lo.push_back( -10 );
		hi.push_back(  10 );
	}
	opti::Particle_Swarm_Optimizer optim( N );
	//opti::Downhill_Simplex optim;
	//opti::Shot_Gun_Optimizer optim;
	optim.initialise( lo, hi );
	//optim.temperature = 0;

	vector<opti::Request> jobs;
	for ( int iter = 0; iter < MAX_ITER; iter++ )
	{
		// do next iteration
		optim.generate_requests( jobs );
		if ( jobs.size() > 0 ) {
			for ( size_t i = 0; i < jobs.size(); i++ )
			{
				jobs[i].y = rosenbrock( jobs[i].x_new );
				//jobs[i].y = paraboloid( jobs[i].x_new );
				//jobs[i].y = rastrigin( jobs[i].x_new );
				jobs[i].flag = 0;
			}
		}
		else {
			cout << "oops \n";
			break;
		}
		optim.assimilate_results( jobs );

		// gather some statistics
		double xm = 0;
		double vm = 0;
		//for ( int i = 0; i < N; i++ )
		for ( int i = 0; i < dim + 1; i++ )
		{
			double xi = 0;
			double vi = 0;
			for ( int d = 0; d < dim; d++ )	{
				//xi += optim.simplex[i].x[d] * optim.simplex[i].x[d];
				//xi += optim.swarm[i].x[d] * optim.swarm[i].x[d];
				//vi += optim.swarm[i].v[d] * optim.swarm[i].v[d];
			}
			xm += sqrt( xi );
			vm += sqrt( vi );
		}
		xm /= N;
		vm /= N;

		// log status
		if ( iter % 1 == 0 ) {
			//optim.temperature = 0.025 * abs( optim.yhi - optim.ylo );
			cout << iter << "\t" << optim.global_min << "\t" << optim.current_min << "\t" << optim.histo_var << "\n";
		}
		//cout << optim.evaluations << "\t" << optim.ylo << "\t" << optim.yhi << "\n";
	}
	return 0;
}
