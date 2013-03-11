#include <iostream>

#include "boinc_api.h"

#include "potential.hpp"
#include "wave_function.hpp"
#include "solver.hpp"

using namespace std;
namespace liee {

void Solver::renormalize()
{
	double norm_factor = 1.0 / sqrt( integrate_psi_sqr( 0.0, r_range ) );
	for ( size_t i = 0; i < Nr; i++ ) {
		psi[i] *= norm_factor;
	}
}

double Solver::integrate_psi_sqr()
{
	return integrate_psi_sqr( 0, r_range );
}

/*!
 *	Integrates the squared Wavefunction amplitude over the interval [a,b]
 *	Caution! The domain for the bounds is [0, r_range] and not [r_start, r_start+r_range]
 *
 *	using extended trapeziodal rule with all available data points of the WF
 */
double Solver::integrate_psi_sqr( double a, double b )
{
	size_t ia = a / dr;
	size_t ib = b / dr;
	if ( ib >= psi.size() ) { ib = psi.size()-1; }

	double sum = -0.5 * real( psi[ia] * conj( psi[ia] ) ) - 0.5 * real( psi[ib] * conj( psi[ib] ) );  //end-points have only half the weight
	for ( size_t i = ia; i <= ib; i++ ) {
		sum += real( psi[i] * conj( psi[i] ) );
	}
	return dr * sum;
}

bool Solver::execute()
{
	const double WORK = Nr * t_end / dt; // samples * frames
	const time_t start_time = time(0);
	time_t last_report_time = start_time;
	time_t last_checkpoint_time = start_time;
	double last_progress = t / t_end;

	while ( t < t_end )
	{
		// first observe (to catch the starting condition as well)
        boinc_begin_critical_section();
		for( size_t i = 0; i < obs.size(); i++ ) {
			obs[i]->observe( this );
		}
        boinc_end_critical_section();

        evolve_1step();

		if ( ++count % 1 == 0 ) { //TODO hmm??
			time_t now = time(0);
			boinc_fraction_done( t / t_end );
			//if ( t / t_end > 0.5 ) return false; //TODO: test if a resumed computation yields the same result
			if ( boinc_time_to_checkpoint() || ( now - last_checkpoint_time > 300 ) ) { 	// in any case checkpoint every 5 min
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
				double hours_remain = (1.0 - last_progress) / current_speed / 3600;
				INFO_SHOW4( count, last_progress, samples_per_s, hours_remain );
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
    boinc_end_critical_section();

    return exec_done;
}

//--------------------------------------------------------------------------------------------

void Crank_Nicholson::register_dependencies( vector<Module*> dependencies )
{
	//TODO when we have more than just the Crank-Nicholson, this code needs to move to Solver::register_dependencies
	for ( size_t i = 0; i < dependencies.size(); i++ ) {
		if ( dependencies[i]->type.compare("potential") == 0 ) {
			LOG_INFO( "found potential to use");;
			potential = dynamic_cast<Potential*>( dependencies[i] );
		}
		else if ( count == 0  &&  dependencies[i]->type.compare("initial_wf") == 0 ) {
			LOG_INFO( "found initial WV to use");;
			Wave_Function* wf = dynamic_cast<Wave_Function*>( dependencies[i] );

			this->psi.resize( Nr );
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
			//TODO check, if observer's->target matches this->serial
			LOG_INFO( "found an observer to use");;
			obs.push_back( dynamic_cast<Observer*>( dependencies[i] ) );
		}
	}
}


void Crank_Nicholson::initialize( Conf_Module* config, vector<Module*> dependencies )
{
	GET_LOGGER( "liee.Crank_Nicholson" );
	dr = config->getParam("dr")->value * 1e-9 / CONV_au_m;
	t_end = config->getParam("t_range")->value * 1e-15 / CONV_au_s;
	r_range = config->getParam("r_range")->value * 1e-9 / CONV_au_m;
	dt = config->getParam("dt")->value * 1e-15 / CONV_au_s;
	Nr = 1 + r_range / dr;
	DEBUG_SHOW2( Nr, t_end/dt );
	jb = Nr - 1;
	c = dcmplx( 0.0, -dt / ( 8.0 * dr*dr ) );   //(44)
	alfa.resize(Nr);
	gamma.resize(Nr);
	g.resize(Nr);
	phi.resize(Nr);
	d.resize(Nr);
	count = 0;
	exec_done = false;
	register_dependencies( dependencies );
	renormalize();
}

void Crank_Nicholson::reinitialize( Conf_Module* config, vector<Module*> dependencies )
{
	GET_LOGGER( "liee.Crank_Nicholson" );
	alfa.resize(Nr);
	gamma.resize(Nr);
	g.resize(Nr);
	phi.resize(Nr);
	d.resize(Nr);
	register_dependencies( dependencies );
}

void Crank_Nicholson::estimate_effort( Conf_Module* config, double & flops, double & ram, double & disk )
{
	double N = config->getParam("r_range")->value / config->getParam("dr")->value;
	double Nt = config->getParam("t_range")->value / config->getParam("dt")->value;

	flops += abs( Nt * N * ( 6*6 + 5*2 + 3*4 + 3) ); // 6 complex(*) + 5 complex(+) + 3 double(/^2) + 3 double(+)  ... for every cell and every step
	ram += abs( 7 * sizeof( dcmplx ) * N + 256 );
	disk += 0;
}


void Crank_Nicholson::evolve_1step()
{
	// update potential and recalculate 'd'
    for (size_t j=0; j < Nr; j++) {
   		d[j] = dcmplx( 0.5, dt / (4.0 * pow(dr, 2)) ) + dcmplx(0.0, dt / 4.0) * potential->V( potential->r_start + dr * j, t ); //#(45)
    }

    alfa[0] = d[0];
    gamma[0] = c / alfa[0];
    for (int j=1; j < jb+1; j++) {
    	alfa[j] = d[j] - c * gamma[j-1]; //#(50)
		gamma[j] = c / alfa[j]; //#(49)
    }

    g[0] = psi[0] / alfa[0]; //#(52)
    for (int j=1; j < jb+1; j++) {
    	g[j] = ( psi[j] - c * g[j-1] ) / alfa[j]; //#(53)
    }

    phi[jb] = g[jb]; //#(55)
    for (int j=jb-1; j > -1; j--) {
    	phi[j] = g[j] - gamma[j] * phi[j+1]; //#(56)
    }

    for (int j=0; j <= jb; j++) {
    	psi[j] = phi[j] - psi[j]; //#(39)
    }

    t += dt;
    count++;
}

} //namespace
