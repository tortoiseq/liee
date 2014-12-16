/*! optimizer.cpp
 *
 *  Implementation of sample optimiser(s).
 */

#include <time.h>

#include "optimizer.hpp"

#ifdef SERVER
#include "sched_config.h"  // only for server-sided  logging
#endif

using namespace std;

namespace liee {

void opti::Asynch_Optimizer::convergence_test( const vector<opti::Request> & population )
{
	current_min = 6.66e300;

	for ( size_t i = 0; i < population.size(); i++ ) {
		if ( population[i].flag == 0 )
		{
			if ( population[i].y < global_min ) {
				global_min = population[i].y;
				global_min_pos = population[i].x;
			}
			if ( population[i].y < current_min ) {
				current_min = population[i].y;
				current_min_pos = population[i].x;
			}
		}
	}

	history[history_i % history_sz] = current_min;
	if ( history_i >= history_sz ) {
		histo_var = abs( current_min - history[0] );
		histo_var = sqrt( variance( history ) );
	}
	else if ( history_i > 2 ) {
		vector<double> sub_hist( history.begin(), history.begin() + history_i );
		histo_var = sqrt( variance( sub_hist ) );
	}
	else { histo_var = 6.66e300; }
	history_i++;
}

//----------------------------------------------------------------------------------------------------------

int opti::Shot_Gun_Optimizer::initialise( const vector<double> & lower, const vector<double> & upper )
{
	unsigned long long timestamp = time( NULL );
	random = new Ranq1( timestamp );
	if ( lower.size() != upper.size() ) return 1;
	dim = lower.size();
	lower_bounds = lower;
	upper_bounds = upper;
	for ( int i = 0; i < num_bullets_per_shot; i++ ) {
		in_flight[i].id = -1;
		in_flight[i].flag = -1;
		in_flight[i].x.resize( dim );
		in_flight[i].x_new.resize( dim );
	}
	return 0;
}

int opti::Shot_Gun_Optimizer::generate_requests( vector<Request> & work )
{
	if ( evaluations >= max_eval ) return 1;  // end condition

	work.clear();
	for ( int i = 0; i < num_bullets_per_shot; i++ ) {
		if ( in_flight[i].id == -1 )  // free slot?
		{
			in_flight[i].id = next_id++;
			// draw random point
			for ( int d = 0; d < dim; d++ ) {
				in_flight[i].x_new[d] = random->doub() * ( upper_bounds[d] - lower_bounds[d] ) + lower_bounds[d];
			}
			work.push_back( in_flight[i] );
			evaluations++;
		}
	}
	return 0;
}

int opti::Shot_Gun_Optimizer::assimilate_results( const vector<Request> & result )
{
	for ( size_t ri = 0; ri < result.size(); ri++ ) {
		for ( int bi = 0; bi < num_bullets_per_shot; bi++ )	{
			if ( in_flight[bi].id == result[ri].id )
			{
			    in_flight[bi].id = -1;  // freeing the slot
				if ( result[ri].flag == 0 )
				{
					in_flight[bi].x = in_flight[bi].x_new;
					in_flight[bi].y = result[ri].y;
					in_flight[bi].flag = 0;
				}
				else {
					//TODO error handling: something for corrupt results
				}
				break;
			}
		}
	}
	convergence_test( in_flight );  // this is silly, it wont converge unless the objective is flat
	return 0;
}

//-------------------------------------------------------------------------------------------------------------

int opti::Particle_Swarm_Optimizer::initialise( const vector<double> & lower, const vector<double> & upper )
{
	unsigned long long timestamp = time( NULL );
	random = new Ranq1( timestamp );
	if ( lower.size() != upper.size() ) return 1;
	dim = lower.size();
	lower_bounds = lower;
	upper_bounds = upper;
	swarm.resize( swarm_sz );
	global_min_pos.resize( dim );

	for ( int i = 0; i < swarm_sz; i++ ) {
		swarm[i].id = -1;
		swarm[i].flag = -1;
		swarm[i].x.resize( dim );
		swarm[i].v.resize( dim );
		swarm[i].x_new.resize( dim );
		swarm[i].x_min.resize( dim );
		swarm[i].min_value = 6.66e300;

		for ( int d = 0; d < dim; d++ ) {
			swarm[i].x[d] = lower_bounds[d] + random->doub() * ( upper_bounds[d] - lower_bounds[d] );
			swarm[i].v[d] = 0.01 * random->doub() * ( upper_bounds[d] - lower_bounds[d] );
			swarm[i].x_min[d] = swarm[i].x[d];
		}
	}

	for ( int d = 0; d < dim; d++ ) {
		global_min_pos[d] = 0.0;
	}
	return 0;
}

int opti::Particle_Swarm_Optimizer::generate_requests( vector<Request> & work )
{
	double MAX_V = 10.0 * exp( -evaluations / swarm_sz / 10000.0 );  //TODO make configurable

	work.clear();
	if ( evaluations >= max_eval ) return 1;  // give up end condition
	if ( histo_var <= tolerance ) return 0;  // success end condition

	for ( int i = 0; i < swarm_sz; i++ )
	{
		if ( swarm[i].id == -1 )  // free slot?
		{
			// find the hood's leader
			int j_best = (i + hood_sz) % swarm_sz;  // default is the last of the neighbourhood, which isn't touched in the for-loop
			for ( int j = i - hood_sz; j < i + hood_sz; j++ )
			{
				int j_ = ( j + swarm_sz ) % swarm_sz;  // add swarm size to make sure dividend is positive
				if ( swarm[j_].min_value < swarm[j_best].min_value ) {
					j_best = j_;
				}
			}

			// adapt particle's velocity
			double v = 0;
			for ( int d = 0; d < dim; d++ )
			{
				swarm[i].v[d] = inertia * swarm[i].v[d]
					+ cognition * random->doub() * ( swarm[  i   ].x_min[d] - swarm[i].x[d] )
					+ coherence * random->doub() * ( swarm[j_best].x_min[d] - swarm[i].x[d] );
				v += swarm[i].v[d] * swarm[i].v[d];
			}
			// limit velocity
			if ( abs( v ) >  MAX_V )
			{
				double scale = abs( MAX_V / v );
				for ( int d = 0; d < dim; d++ ) {
					swarm[i].v[d] *= scale;
				}
			}
			// change particle's (new) position
			for ( int d = 0; d < dim; d++ ) {
				swarm[i].x_new[d] = swarm[i].x[d] + swarm[i].v[d];
				//TODO test bounds
			}

			swarm[i].id = next_id++;
			work.push_back( swarm[i] );
			evaluations++;
		}
	}
	return 0;
}

int opti::Particle_Swarm_Optimizer::assimilate_results( const vector<Request> & result )
{
	//process results
	for ( size_t r = 0; r < result.size(); r++ ) {
		for ( int i = 0; i < swarm_sz; i++ ) {
			if ( swarm[i].id == result[r].id )
			{
				swarm[i].id = -1;  // in any case: make available for next round
				if ( result[i].flag == 0 )
				{
					swarm[i].x = swarm[i].x_new;
					swarm[i].y = result[r].y;
					swarm[i].flag = 0;

					if ( swarm[i].y < swarm[i].min_value )
					{
						swarm[i].min_value = swarm[i].y;
						swarm[i].x_min = swarm[i].x;
						if ( swarm[i].y < global_min )
						{
							global_min = swarm[i].y;
							global_min_pos = swarm[i].x;
						}
					}
				}
				else { } //TODO implement
				//TODO time-out for late results
			}
		}
	}
	//elect new leader
	leader = 0;
	for ( int i = 1; i < swarm_sz; i++ ) {
		if ( swarm[i].y < swarm[leader].y ) {
			leader = i;
		}
	}
	convergence_test( swarm );
	return 0;
}

void opti::Particle_Swarm_Optimizer::convergence_test( const vector<Particle> & population )
{
	//TODO remove dirty'n'slow hack and play by the book of polymorphism
	vector<Request> cast_back;
	for (size_t i = 0; i < population.size(); i++) {
		cast_back.push_back( population[i] );
	}
	Asynch_Optimizer::convergence_test( cast_back );
}

//------------------------------------------------------------------------------------------------------------------
int opti::Rasterizer::initialise( const vector<double> & lower, const vector<double> & upper )
{
	unsigned long long timestamp = time( NULL );
	random = new Ranq1( timestamp );
	dim = lower.size();
	lower_bounds = lower;
	upper_bounds = upper;
	if ( ((int)logscale.size()) != dim ) {
		logscale.resize( dim );
		for ( int d = 0; d < dim; d++ ) { logscale[d] = false; }
	}
	discharged = false;
	if ( ((int)num_samples.size()) != dim ) {  // apply uniform N to all dimensions
		int N = num_samples[0];
		num_samples.resize( dim );
		for ( size_t i = 0; i < num_samples.size(); i++ ) {
			num_samples[i] = N;
		}
	}
	return 0;
}

int opti::Rasterizer::initialise( const vector<Conf_Param*> & params )
{
	vector<double> lo, hi;
	BOOST_FOREACH( Conf_Param* p, params ) {
		if ( not p->fixed ) {
			if ( p->logscale )
			{
				if ( p->min == 0.0 || p->max == 0.0 ) {
					p->logscale = false;
					cout <<  "WARNING: Bound of " << p->name << " set to zero despite log-scale. Logarithmic scaling DISABLED! \n";
				} else {
					cout <<  "Using logarithmic scaling for " << p->name << "\n";
					p->min = log( p->min );
					p->max = log( p->max );
				}
			}
			logscale.push_back( p->logscale );
			cout << "adding param " << p->name << " to the variable dimensions \n";
			lo.push_back( p->min );
			hi.push_back( p->max );
		}
	}
	return initialise( lo, hi );
}

/*!
 * Recursively iterate over all dimension's ranges using backtracking.
 * @param d dimension as backtrack pointer
 * @param r is the working record
 * @param work is the vector to which all requests finally get added (in last stage)
 */
void opti::Rasterizer::iterate( int d, Request r, vector<Request> & work )
{
	if ( d < 0 ) {
		work.push_back( r );
	}
	else {
		for ( double i = 0; i < num_samples[d]; i++ )
		{
			r.x[d] = i;	// store index to give a hint at the work-unit's naming
			r.x_new[d] =  i / ( num_samples[d] - 1 ) * ( upper_bounds[d] - lower_bounds[d] ) + lower_bounds[d];
			if ( logscale[d] ) {
				r.x_new[d] = exp( r.x_new[d] );
			}
			iterate( d - 1, r, work );
		}
	}
}

int opti::Rasterizer::generate_requests( vector<Request> & work )
{
	cout << "...opti::Rasterizer::generate_requests() \n";
	if ( discharged ) {
		//TODO 1) use logger in #else case 2) really critical or maybe just warning level?
		#ifdef SERVER
			log_messages.printf( MSG_CRITICAL, "This Rasterizer created all its jobs already. Re-initialize if you want to play again! \n" );
		#else
			cout << "This Rasterizer created all its jobs already. Re-initialize if you want to play again! \n";
		#endif
		return 1;
	}
	work.clear();
	Request req;
	req.x.resize( dim );
	req.x_new.resize( dim );
	iterate( dim - 1, req, work );
	discharged = true;
	return 0;
}

double opti::Golden_Section_Search::minimize( boost::function<double (double)> func, const double a, const double b, const double c )
{
	double xa = a;
	double xb = b;
	double xc = c;
	double xd;
	double fa = func( xa );
	double fb = func( xb );
	double fc = func( xc );
	double fd;
	double gs = ( 3.0 - sqrt( 5.0 ) ) / 2.0;

	if (not ((a < b) && (b < c) && (fb < fa) && (fb < fc)) ) {
		throw Except__Preconditions_Fail( __LINE__ );
	}

	while ( abs( xc - xa ) > tol )
	{
		double ab = abs( xb - xa );
		double bc = abs( xb - xc );
		if ( ab > bc ) {  // bisect [ab] for it is the bigger interval
			xd = xb - gs * ab;
			fd = func( xd );
			if ( fd < fb ) {  // test point d is the lowest --> exclude c
				xc = xb; fc = fb;
				xb = xd; fb = fd;
			}
			else {  // b remains the lowest --> pull in a to d
				xa = xd; fa = fd;
			}
		}
		else {  // bisect interval [bc]
			xd = xb + gs * bc;
			fd = func( xd );
			if ( fd < fb ) {  // test point d is the lowest --> exclude a
				xa = xb; fa = fb;
				xb = xd; fb = fd;
			}
			else {  // b remains the lowest --> pull in c to d
				xc = xd; fc = fd;
			}
		}
	}
	fmin = fb;
	return xb;
}

} // namespace liee
