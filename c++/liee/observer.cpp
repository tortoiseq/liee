#include <fstream>

#include "filesys.h"
#include "boinc_api.h"

#include <boost/function.hpp>
#include <boost/bind.hpp>

#include "my_util.hpp"
#include "solver.hpp"
#include "observer.hpp"
#include "wave_function.hpp"

using namespace std;

namespace liee {

void Obs_Snapshot_WF::initialize( Conf_Module* config, vector<Module*> dependencies )
{
	GET_LOGGER( "liee.Obs_Snapshot_WF" );
	// prepare spatial downsampling
	double r_range = config->getParam("r_range")->value / CONV_au_nm;
	double dr = config->getParam("dr")->value / CONV_au_nm;

	double r0 = 0.0;
	double r1 = r_range;
	if ( not config->getParam("obs_r0")->textual ) { r0 = config->getParam("obs_r0")->value / CONV_au_nm; }
	if ( not config->getParam("obs_r1")->textual ) { r1 = config->getParam("obs_r1")->value / CONV_au_nm; }
	ir0 = (int) (r0 / dr);
	ir1 = (int) (r1 / dr);
	DEBUG_SHOW2(ir0, ir1);

	int Nr = 1.0 + ( r1 - r0 ) / dr;
	step_r = floor( Nr / config->getParam("r_samples")->value );
	if ( step_r < 1 ) { step_r = 1; }
	if ( step_r > Nr ) { step_r = Nr; }
	config->getParam("r_samples")->value = Nr;	// save the actual number of samples

	// prepare temporal downsampling
	t_range = config->getParam("t_range")->value / CONV_au_fs;
	dt = config->getParam("dt")->value / CONV_au_fs;
	t0 = 0.0;
	t1 = t_range;
	if ( not config->getParam("obs_t0")->textual ) { t0 = config->getParam("obs_t0")->value / CONV_au_fs; }
	if ( not config->getParam("obs_t1")->textual ) { t1 = config->getParam("obs_t1")->value / CONV_au_fs; }
	counter = 0;
	int Nt = 1.0 + ( t1 - t0 ) / dt;
	step_t = floor( Nt / config->getParam("t_samples")->value );
	if ( step_t < 1 ) { step_t = 1; }
	if ( step_t > Nt ) { step_t = Nt; }
	config->getParam("t_samples")->value = Nt;	// save the actual number of samples

	do_square = config->getParam("square")->text.compare("true") == 0;
	stringstream ss;
	int digits = int( config->getParam("precision")->value );
	if (do_square) {
		ss << "%1." << digits << "g\t";
	} else {
		ss << "%1." << digits << "g,%1." << digits << "g\t";
	}
	format = ss.str();

	boinc_resolve_filename_s( config->getParam("OUTFILE")->text.c_str(), filename );
	do_normalize = config->getParam("normalize")->text.compare("true") == 0;
	// delete previous snapshot file
	FILE* f = boinc_fopen( filename.c_str(), "w" );
	fclose( f );
}

void Obs_Snapshot_WF::reinitialize( Conf_Module* config, vector<Module*> dependencies )
{
	GET_LOGGER( "liee.Obs_Snapshot_WF" );
}

void Obs_Snapshot_WF::estimate_effort( Conf_Module* config, double & flops, double & ram, double & disk )
{
	double N = config->getParam("r_range")->value / config->getParam("dr")->value;
	double Nt = config->getParam("t_range")->value / config->getParam("dt")->value;
	double Nr = config->getParam("r_samples")->value;
	int samples = (int) config->getParam("t_samples")->value;
	bool sqr = config->getParam("square")->text.compare("true") == 0;

	flops += abs( 8 * Nt + samples * N * 20 );	// called Nt times without saving + sample times with saving
	ram += 1024;
	double dsk = (samples + 2) * Nr * 2 * 23 * 1.5; // factor: complex:2, double2str:23, outfile duplicated in gz archive: 1.5
	if ( sqr ) dsk /= 2;
	disk += abs( dsk );
}

void Obs_Snapshot_WF::observe( Module* state )
{
	Solver* s = dynamic_cast<Solver*>( state );
	if ( s->t < t0  ||  s->t > t1 ) return;		// not in observation window
	if ( counter++ % step_t != 0 ) return;		// downsample

   	double fac = 1.0;
   	if ( do_normalize ) {
   		fac = 1.0 / s->integrate_psi_sqr();
   	}

    FILE *file;
	file = boinc_fopen( filename.c_str(), "a" );

   	for ( size_t i = ir0; i <= ir1 && i < s->Nr; i += step_r ) {
		double re, im;
    	if ( do_square ) {
    		re = fac * real( s->psi[i] * conj( s->psi[i] ) );
    		fprintf( file, format.c_str(), re );
    	}
    	else {
    		re = fac * real( s->psi[i] );
    		im = fac * imag( s->psi[i] );
    		fprintf( file, format.c_str(), re, im );
    	}
	}
	fprintf( file, "\n" );
   	fclose( file );
}

//-------------------------------------------------------------------------------------------------------------

void Obs_Tunnel_Ratio::initialize( Conf_Module* config, vector<Module*> dependencies )
{
	GET_LOGGER( "liee.Obs_Tunnel_Ratio" );
	Potential* pot;
	for ( size_t i = 0; i < dependencies.size(); i++ ) {
		if ( dependencies[i]->type.compare( "potential" ) == 0 ) {
			pot = dynamic_cast<Potential*>( dependencies[i] );
		}
	}

	filename = config->getParam("OUTFILE")->text;
	is_objective = config->getParam("is_objective")->text.compare( "true" ) == 0;
	ra = config->getParam("r_a")->value / CONV_au_nm;
	rb = config->getParam("r_b")->value / CONV_au_nm;
	if ( ra < 0 ) ra = 0;
	if ( rb > pot->get_r_range() ) { rb = pot->get_r_range(); }

	// temporal downsampling
	t_range = config->getParam("t_range")->value / CONV_au_fs;
	dt = config->getParam("dt")->value / CONV_au_fs;
	counter = 0;
	int Nt = 1.0 + t_range / dt;
	step_t = floor( Nt / config->getParam("t_samples")->value );
	dt *= step_t;
	if ( step_t < 1 ) { step_t = 1; }
	if ( step_t > Nt ) { step_t = Nt; }
	t_samples = 1 + Nt / step_t;
	config->getParam("t_samples")->value = t_samples;	// save the actual number of samples
}

void Obs_Tunnel_Ratio::reinitialize( Conf_Module* config, vector<Module*> dependencies )
{
	GET_LOGGER( "liee.Obs_Tunnel_Ratio" );
}

void Obs_Tunnel_Ratio::estimate_effort( Conf_Module* config, double & flops, double & ram, double & disk )
{
	//TODO account for (ra..rb) integration
	double Nr = config->getParam("r_range")->value / config->getParam("dr")->value;
	double Nt = config->getParam("t_range")->value / config->getParam("dt")->value;
	double N =  config->getParam("t_samples")->value;

	flops += abs( 8 * Nt + N * Nr * 8 );	// called Nt times without saving + t_sample times to integrate (psi psi*)
	ram += 1024;
	disk += 30 * N;
}

void Obs_Tunnel_Ratio::summarize( map<string, string> & results )
{
	results["tunnel_ratio"] = doub2str( tunnel_ratio ); //TODO boost conversion
	if ( is_objective ) {
		results["objective"] = doub2str( tunnel_ratio );
	}
}

void Obs_Tunnel_Ratio::observe( Module* state )
{
	if ( counter++ % step_t != 0 ) return;
	Solver* s = dynamic_cast<Solver*>( state );
	psi_sqr.push_back( s->integrate_psi_sqr( ra, rb ) );

	// save after getting the last sample
	if ( (int)psi_sqr.size() >= t_samples-1 ) {	// for the sake of robustness we are willing to write the outfile at last-1 step and overwrite at last step, to make sure it gets written at least once.
											//TODO: check log output, if t_samples calculation is always spot on -> two writes, and remove the additional write out at last-1
		LOG_DEBUG( "Writing Tunnel data at time-step " << psi_sqr.size() << " of " << t_samples );
    	tunnel_ratio = 1.0  -  psi_sqr.back() / psi_sqr.front();
    	FILE* f = boinc_fopen( filename.c_str(), "w" );
    	double last = psi_sqr.front();
    	for ( size_t i = 0; i < psi_sqr.size(); i++ )
    	{
    		double j = ( last - psi_sqr[i] ) / dt;
    		last = psi_sqr[i];
    		fprintf( f, "%1.6g\t%1.15g\t%1.10g\n", i * dt , psi_sqr[i], j );
    	}
    	fprintf( f, "\n" );
    	fclose( f );
    }
}

//-------------------------------------------------------------------------------------------------------------

void Obs_JWKB_Tunnel::initialize( Conf_Module* config, vector<Module*> dependencies )
{
	GET_LOGGER( "liee.Obs_JWKB_Tunnel" );
	//register dependencies
	Wave_Function* wf;
	for ( size_t i = 0; i < dependencies.size(); i++ ) {
		if ( dependencies[i]->type.compare( "potential" ) == 0 ) {
			V = dynamic_cast<Potential*>( dependencies[i] );
		}
		else if ( dependencies[i]->type.compare( "initial_wf" ) == 0 ) {
			wf = dynamic_cast<Wave_Function*>( dependencies[i] );
		}
	}
	filename = config->getParam("OUTFILE")->text;
	is_objective = config->getParam("is_objective")->text.compare( "true" ) == 0;
	dr = config->getParam("dr")->value / CONV_au_nm;
	N = (int)config->getParam("int_samples")->value / CONV_au_nm;

	//calculate E from curvature of WF: E = V - hbar^2/(2m Psi) d^2/dr^2 Psi
	double rmin = V->getPot_const()->get_Vmin_pos();
	double Vmin = V->getPot_const()->V(rmin);
	Vp0 = V->getPot_const()->V(dr);
	double r0 = V->get_r_start();
	int i_middle = (int)( (rmin - r0) / dr );
	E = Vmin - ( wf->psi[i_middle-1].real() - 2 * wf->psi[i_middle].real() + wf->psi[i_middle+1].real() ) / ( pow( dr, 2.0 ) * 2.0 * wf->psi[i_middle].real() );
	DEBUG_SHOW(E);
	g = 4.0 * CONST_PI * sqrt( 2.0 );

	// temporal downsampling
	// (code duplication with tunnel ratio observer)
	t_range = config->getParam("t_range")->value / CONV_au_fs;
	dt = config->getParam("dt")->value / CONV_au_fs;
	counter = 0;
	int Nt = 1.0 + t_range / dt;
	step_t = floor( Nt / config->getParam("t_samples")->value );
	dt *= step_t;
	if ( step_t < 1 ) { step_t = 1; }
	if ( step_t > Nt ) { step_t = Nt; }
	t_samples = 1 + Nt / step_t;
	config->getParam("t_samples")->value = t_samples;	// save the actual number of samples
}

void Obs_JWKB_Tunnel::reinitialize( Conf_Module* config, vector<Module*> dependencies )
{
	GET_LOGGER( "liee.Obs_JWKB_Tunnel" );
	for ( size_t i = 0; i < dependencies.size(); i++ ) {
		if ( dependencies[i]->type.compare( "potential" ) == 0 ) {
			V = dynamic_cast<Potential*>( dependencies[i] );
		}
	}
}

void Obs_JWKB_Tunnel::estimate_effort( Conf_Module* config, double & flops, double & ram, double & disk )
{
	double Nr = config->getParam("r_range")->value / config->getParam("dr")->value;
	double Nt = config->getParam("t_range")->value / config->getParam("dt")->value;
	double N =  config->getParam("t_samples")->value;

	flops += abs( 8 * Nt + N * Nr * 10 );	// called Nt times without saving + t_sample times to integrate over the squareroot of barrier hight
	ram += 1024;
	disk += 20 * N;
}

void Obs_JWKB_Tunnel::summarize( map<string, string> & results )
{
	results["jwkb_J"] = doub2str( sum_j ); //TODO boost conversion
	if ( is_objective ) {
		results["objective"] = doub2str( sum_j );
	}
}

void Obs_JWKB_Tunnel::observe( Module* state )
{
	if ( counter++ % step_t != 0 ) return;
	Solver* s = dynamic_cast<Solver*>( state );
	boost::function<double(double)> deltaE;
	deltaE = boost::bind( &Potential::deltaV, V, _1, s->t, E );

	// estimate actual field strength from the change in potentials height compared with t0
	double F = ( Vp0 - V->V( dr, s->t ) ).real() / dr;
	if (F <= 0 ) {
		// barrier goes up +++ assume slow varying F --> infinite barrier --> j=0
		j.push_back( 0.0 );
	}
	else {
		// find left turning point starting at r==0
		double r1 = 0;
		double d = dr;
		double a = deltaE( r1 );
		if ( a < 0 ) { d = -dr; } // r1 already inside barrier -->, direction backwards
		while ( deltaE( r1 + d ) * a > 0 ) {
			d *= 2;
			if ( isinf( d ) ) {	throw Except__Too_Far_Out( __LINE__ ); }
		}
		r1 = find_root( deltaE, r1, r1 + d, 1e-12 );

		//TODO resolve code duplication with Pot_const::get_outer_turningpoints() by writing a robust and general root-finding routine
		// find right turning point
		// assume constant F +++ assume V approaching 0 for r > 0
		double r2 = -E / F;

		if ( abs( deltaE( r2 ) ) > abs( deltaE( last_r2 ) ) ) {
			// last turning-point was better than new estimate
			r2 = last_r2;
		}
		d = dr;
		a = deltaE( r2 );
		if ( a > 0 ) { d = -dr; } // r1 already outside barrier --> direction backwards
		while ( deltaE( r2 + d ) * a > 0 ) {
			d *= 2;
			if ( isinf( d ) ) {	throw Except__Too_Far_Out( __LINE__ ); }
		}
		r2 = find_root( deltaE, r2, r2 + d, 1e-12 );

		DEBUG_SHOW2(r1, r2);

		// integrate: dr sqrt( V(r) - E )
		d = ( r2 - r1 ) / ( N - 1.0 );
		vector<Point> samples;
		for ( double r = r1; r <= r2; r += d ) {
			samples.push_back( Point( r, sqrt( -deltaE(r) ) ) );
		}

		double A = 0;
	#ifdef ALGLIB
		alglib::spline1dinterpolant spline = to_cubic_spline( samples );
		double A0_xa = alglib::spline1dintegrate( spline, r1 );
		double A0_xb = alglib::spline1dintegrate( spline, r2 );
		A = A0_xb - A0_xa;
	#else
		A = simple_integrate( samples, r1, r2 );
	#endif
		DEBUG_SHOW( A );
		j.push_back( exp( -g * A ) );
		DEBUG_SHOW( j.back() );
	}

	// save after getting the last sample (again some code duplication with tunnel observer)
	if ( (int)j.size() >= t_samples-1 ) {	// for the sake of robustness we are willing to write the outfile at last-1 step and overwrite at last step, to make sure it gets written at least once.
											//TODO: check log output, if t_samples calculation is always spot on -> two writes, and remove the additional write out at last-1
		LOG_DEBUG( "Writing Tunnel data at time-step " << j.size() << " of " << t_samples );
		sum_j = sum(j);
    	FILE* f = boinc_fopen( filename.c_str(), "w" );
    	for ( size_t i = 0; i < j.size(); i++ )	{
    		fprintf( f, "%1.6g\t%1.10g\n", i * dt , j[i] );
    	}
    	fprintf( f, "\n" );
    	fclose( f );
    }
}

} // namespace liee
