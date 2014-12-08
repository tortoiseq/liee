#include <iostream>

#include "filesys.h"
#include "boinc_api.h"

#include "boost/math/special_functions/erf.hpp"
#include "boost/thread.hpp"
#include "boost/bind.hpp"

#include "potential.hpp"
#include "wave_function.hpp"
#include "solver.hpp"

using namespace std;
namespace liee {

void Solver::renormalize()
{
	double norm_factor = 1.0 / sqrt( integrate_psi_sqr( 0.0, r_range ) );
	if ( isinf( norm_factor ) || isnan( norm_factor ) ) {
		LOG_FATAL("Normalization of WaveFunction failed. (check that initial WF is non-zero)");
		exit(1);
	}
	for ( size_t i = 0; i < Nr; i++ ) {
		psi[i] *= norm_factor;
	}
}

double Solver::integrate_psi_sqr()
{
	return integrate_psi_sqr( 0, r_range );
}

/*!
 * Integrates the squared Wavefunction amplitude over the interval [a,b]
 * Caution! The domain for the bounds is [0, r_range] and not [r_start, r_start+r_range]
 *
 * using extended trapeziodal rule with all available data points of the WF
 */
double Solver::integrate_psi_sqr( double a, double b )
{
	size_t ia = a / dr;
	size_t ib = b / dr;
	if ( ib >= psi.size() ) { ib = psi.size()-1; }

	Kahan_Summator<double> sum;
	sum.add( -0.5 * real( psi[ia] * conj( psi[ia] ) ) - 0.5 * real( psi[ib] * conj( psi[ib] ) ) );  //end-points have only half the weight
	for ( size_t i = ia; i <= ib; i++ ) {
		sum.add( real( psi[i] * conj( psi[i] ) ) );
	}
	return dr * sum.get_sum();
}

bool Solver::execute()
{
	const double WORK = Nr * t_end / dt;  // samples * frames
	const time_t start_time = time(0);
	time_t last_report_time = start_time;
	time_t last_checkpoint_time = start_time;
	double last_progress = t / t_end;
	int count_tic = 100;

	while ( t < t_end )
	{
		boinc_begin_critical_section();
		for( size_t i = 0; i < obs.size(); i++ ) {
			obs[i]->observe( this );
		}
		boinc_end_critical_section();

		evolve_1step();  //<<<<<<<<<<<<<<< this is where the action is
		t += dt;

		if ( ++count % count_tic == 0 ) {
			time_t now = time(0);
			boinc_fraction_done( t / t_end );
			if ( boinc_time_to_checkpoint() || ( now - last_checkpoint_time > 300 ) ) {  // in any case checkpoint every 5 min
				last_checkpoint_time = now;
				return false;
			}

			// progress info
			double wait = now - last_report_time;
			if ( wait > 20 ) {
				last_report_time = now;
				double leap = t / t_end - last_progress;
				last_progress = t / t_end;

				double current_speed = leap / wait;
				double samples_per_s = current_speed * WORK;
				double minutes_remain = (1.0 - last_progress) / current_speed / 60.0;
				count_tic = 1 + (int)abs( 10.0 * samples_per_s / Nr );  // check progress about once in 10s
				INFO_SHOW4( count, last_progress, samples_per_s, minutes_remain );
			}
		}
	}
	// observe the final state
	t = t_end;
	exec_done = true;

	boinc_begin_critical_section();
	for( size_t i = 0; i < obs.size(); i++ ) {
		obs[i]->observe( this );
	}
	// optionally save final state
	if ( outfile.length() > 0 ) {
		FILE* f = boinc_fopen( outfile.c_str(), "w" );
		double r_start = potential->get_r_start();
		for ( size_t i = 0; i < psi.size(); i++ ) {
			fprintf( f, "%1.16g\t%1.16g\t%1.16g\n", r_start + i * dr , real( psi[i] ), imag( psi[i] ) );
		}
		fclose( f );
	}
	boinc_end_critical_section();

	return exec_done;
}

//--------------------------------------------------------------------------------------------

void Crank_Nicholson::register_dependencies( vector<Module*> dependencies )
{
	//TODO when we have more than just the Crank-Nicholson, this code should move to Solver::register_dependencies
	for ( size_t i = 0; i < dependencies.size(); i++ ) {
		if ( dependencies[i]->type.compare("potential") == 0 ) {
			LOG_DEBUG( "found potential to use: " << dependencies[i]->name );
			potential = dynamic_cast<Potential*>( dependencies[i] );
		}
		else if ( count == 0  &&  dependencies[i]->type.compare("initial_wf") == 0 ) {
			LOG_DEBUG( "found initial WV to use: " << dependencies[i]->name );
			Wave_Function* wf = dynamic_cast<Wave_Function*>( dependencies[i] );

			psi.resize( Nr );
			for( size_t i = 0; i < Nr; i++ ) {
				if ( i < wf->psi.size() ) {
					psi[i] = wf->psi[i];
				}
				else {
					psi[i] = dcmplx( 0, 0 );
				}
			}
		}
		else if ( dependencies[i]->type.compare("observer") == 0 ) {
			Observer* observer = dynamic_cast<Observer*>( dependencies[i] );
			if ( observer->target == my_id ) {
				obs.push_back( observer );
				LOG_DEBUG( "registered observer: " << dependencies[i]->name );
			}
		}
	}
}


void Crank_Nicholson::initialize( Conf_Module* config, vector<Module*> dependencies )
{
	GET_LOGGER( "liee.Crank_Nicholson" );
	my_id = config->serial;
	dr = config->get_double("dr") / CONV_au_nm;
	t = 0;
	t_end = config->get_double("t_range") / CONV_au_fs;
	r_range = config->get_double("r_range") / CONV_au_nm;
	dt = config->get_double("dt") / CONV_au_fs;
	dt_ = dt;
	Nr = 1 + r_range / dr;
	LOG_INFO( "Problem size (Nr x Nt): " << Nr << " x " << (int)(t_end / dt) << " = " << (Nr * t_end / dt) );
	c1 = dcmplx( 0.0, -dt / ( 8.0 * SQR(dr) ) );  //(44)
	c3 = dcmplx( 0.0,  dt / 4.0 );
	gamma.resize(Nr);
	g.resize(Nr);
	count = 0;
	exec_done = false;
	register_dependencies( dependencies );
	renormalize();
	potential->set_grid_CN( dr, dt, Nr );
	outfile = config->get_string("OUTFILE");
	multithreaded = false;
	if (multithreaded) {
		d.resize(Nr);
		mutex_d_buffer.lock();
			producer_pos = 0;
			consumer_pos = 0;
			t_threadsafe = 0;
			producer_sleeps = false;
			consumer_sleeps = false;
		mutex_d_buffer.unlock();
		boost::thread producer_fred( &Crank_Nicholson::fill_d_thread, this );
	}
}

void Crank_Nicholson::reinitialize( Conf_Module* config, vector<Module*> dependencies )
{
	GET_LOGGER( "liee.Crank_Nicholson" );
	gamma.resize(Nr);
	g.resize(Nr);
	register_dependencies( dependencies );
	potential->set_grid_CN( dr, dt, Nr );
	if (multithreaded) {
		d.resize(Nr);
		mutex_d_buffer.lock();
			producer_pos = 0;
			consumer_pos = 0;
			t_threadsafe = t;
			producer_sleeps = false;
			consumer_sleeps = false;
		mutex_d_buffer.unlock();
		boost::thread producer_fred( &Crank_Nicholson::fill_d_thread, this );
	}
}

void Crank_Nicholson::estimate_effort( Conf_Module* config, double & flops, double & ram, double & disk )
{
	double N = config->get_double("r_range") / config->get_double("dr");
	double Nt = config->get_double("t_range") / config->get_double("dt");

	flops += abs( Nt * N * ( 6*6 + 5*2 + 3*4 + 3) );  // 6 complex(*) + 5 complex(+) + 3 double(/^2) + 3 double(+)  ... for every cell and every step
	ram += abs( 3 * sizeof( dcmplx ) * N + 256 );
	disk += 0;
}


void Crank_Nicholson::fill_d_thread()
{
	mutex_d_buffer.lock();
		size_t my_pos = producer_pos;
	mutex_d_buffer.unlock();

	for (;;) {
		// update state
		mutex_d_buffer.lock();
			producer_pos = my_pos;
			size_t yo_pos = consumer_pos;
			double my_t = t_threadsafe;
		mutex_d_buffer.unlock();

		// if almost full --> sleep
		if ( (my_pos - yo_pos) > 0.9 * Nr ) {
			boost::mutex::scoped_lock lock( mutex_d_buffer );
			producer_sleeps = true;
			poke_producer.wait( lock );
			producer_sleeps = false;
			continue;  // re-read new state
		}

		// fill as much as possible
		for ( ; my_pos < (yo_pos + Nr); my_pos++ ) {
			size_t i = my_pos % Nr;
			// update every 20%
			if ( my_pos % (Nr/5) == 0 ) {
				mutex_d_buffer.lock();
					producer_pos = my_pos;
				mutex_d_buffer.unlock();
			}

			// wake up consumer when half full
			if ( consumer_sleeps && (my_pos - yo_pos) > Nr/2 ) {
				mutex_d_buffer.lock();
					producer_pos = my_pos;
					poke_consumer.notify_one();
				mutex_d_buffer.unlock();
			}

			double tt = (my_pos % Nr < yo_pos % Nr) ? my_t+dt : my_t;
			d[i] = c3 * potential->V_indexed( i, tt );  //#(45) // all this infrastructure for this single line of computation
		}
	}
}

void Crank_Nicholson::evolve_1step_MT()
{
	mutex_d_buffer.lock();
		t_threadsafe = t;
		size_t my_pos = consumer_pos;
		size_t yo_pos = producer_pos;
	mutex_d_buffer.unlock();

	size_t start = my_pos;
	for ( ; my_pos < start + Nr; my_pos++ ) {
		size_t i = my_pos % Nr;
		// check if data available, otherwise --> sleep
		while ( yo_pos == my_pos ) {
			mutex_d_buffer.lock();
				yo_pos = producer_pos;
			mutex_d_buffer.unlock();
			if ( yo_pos != my_pos ) { break; }

			boost::mutex::scoped_lock lock( mutex_d_buffer );
			consumer_sleeps = true;
			poke_consumer.wait( lock );
			consumer_sleeps = false;
			// woken up --> check again
		}

		// update every 20%
		if ( my_pos % (Nr/5) == 0 ) {
			mutex_d_buffer.lock();
				consumer_pos = my_pos;
				yo_pos = producer_pos;
			mutex_d_buffer.unlock();
		}

		// wake up producer if half empty
		if ( producer_sleeps && (yo_pos - my_pos) < Nr/2 ) {
			mutex_d_buffer.lock();
				consumer_pos = my_pos;
				poke_producer.notify_one();
				mutex_d_buffer.unlock();
		}

		dcmplx alfa = d[i] - c1 * gamma[i-1];  //#(50)
		gamma[i] = c1 / alfa;  //#(49)
		g[i] = ( psi[i] - c1 * g[i-1] ) / alfa;  //#(53)
	}
	// done here, set consumer_pos for the next time-step
	mutex_d_buffer.lock();
		consumer_pos = my_pos;
	mutex_d_buffer.unlock();

	dcmplx alfa = d[0];
	gamma[0] = c1 / alfa;
	g[0] = psi[0] / alfa;  //#(52)
	for (size_t j=1; j < Nr; j++) {
		alfa = d[j] - c1 * gamma[j-1];  //#(50)
		gamma[j] = c1 / alfa;  //#(49)
		g[j] = ( psi[j] - c1 * g[j-1] ) / alfa;  //#(53)
	}

	dcmplx phi = g[Nr-1];  //#(55)
	psi[Nr-1] = phi - psi[Nr-1];  //#(39)
	for (int j=Nr-2; j >= 0; j--) {
		phi = g[j] - gamma[j] * phi;  //#(56)
		psi[j] = phi - psi[j];  //#(39)
	}
}

void Crank_Nicholson::evolve_1step()
{
	//if ( multithreaded ) { evolve_1step_MT(); return; }

	dcmplx alfa = c3 * potential->V_indexed( 0, t );
	gamma[0] = c1 / alfa;
	g[0] = psi[0] / alfa;  //#(52)
	for (size_t j=1; j < Nr; j++) {
		alfa = c3 * potential->V_indexed( j, t ) - c1 * gamma[j-1];  //#(45),#(50)
		gamma[j] = c1 / alfa;  //#(49)
		g[j] = ( psi[j] - c1 * g[j-1] ) / alfa;  //#(53)
	}

	dcmplx phi = g[Nr-1];  //#(55)
	psi[Nr-1] = phi - psi[Nr-1];  //#(39)
	for (int j=Nr-2; j >= 0; j--) {
		phi = g[j] - gamma[j] * phi;  //#(56)
		psi[j] = phi - psi[j];  //#(39)
	}
}

} //namespace
