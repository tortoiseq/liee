#include <fstream>

#include "filesys.h"
#include "boinc_api.h"

#include <fftw3.h>

#include <boost/foreach.hpp>
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
	counter = 0;
	writtenLns = 0;
	foundLns = 0;
	do_square = config->getParam("square")->text.compare("true") == 0;
	do_fourier = config->getParam("fourier")->text.compare("true") == 0;
	do_normalize = config->getParam("normalize")->text.compare("true") == 0;
	rel_change = config->getParam("rel_change")->text.compare("true") == 0;
	N = (int)( 1.0 + config->getParam("r_range")->value / config->getParam("dr")->value );

	stringstream ss;
	int digits = int( config->getParam("precision")->value );
	if ( do_square ) {
		ss << "%1." << digits << "g\t";
	} else {
		ss << "%1." << digits << "g\t%1." << digits << "g\t";
	}
	format = ss.str();

	double dr = config->getParam("dr")->value / CONV_au_nm;
	if ( not do_fourier ) {
		// prepare spatial downsampling
		double r_range = config->getParam("r_range")->value / CONV_au_nm;

		double r0 = ( config->getParam("obs_r0")->textual )  ?  0       :  config->getParam("obs_r0")->value / CONV_au_nm;
		double r1 = ( config->getParam("obs_r1")->textual )  ?  r_range :  config->getParam("obs_r1")->value / CONV_au_nm;
		ir0 = (int) (r0 / dr);
		ir1 = (int) (r1 / dr);

		size_t Nr = 1.0 + ( r1 - r0 ) / dr;
		step_r = floor( Nr / config->getParam("r_samples")->value );
		if ( step_r < 1 ) { step_r = 1; }
		if ( step_r > Nr ) { step_r = Nr; }
		num_r = 1 + Nr / step_r;
		config->getParam("r_samples")->value = num_r;
	} else {
		// prepare spectral downsampling
		double k_range = 0.5 / ( config->getParam("dr")->value / CONV_au_nm );  // k_max = 1/(2 dr)
		double dk = 2.0 * k_range / ( config->getParam("r_range")->value / config->getParam("dr")->value );

		double E0 = ( config->getParam("obs_E0")->textual )  ?  - pow( k_range / CONST_PI, 2.0 ) / 8.0  :  config->getParam("obs_E0")->value / CONV_au_eV;
		double E1 = ( config->getParam("obs_E1")->textual )  ?  + pow( k_range / CONST_PI, 2.0 ) / 8.0  :  config->getParam("obs_E1")->value / CONV_au_eV;
		double k0 = sign(E0) * 2.0 * CONST_PI * sqrt( 2.0 * abs(E0) );
		double k1 = sign(E1) * 2.0 * CONST_PI * sqrt( 2.0 * abs(E1) );
		if ( k0 < -k_range ) { k0 = -k_range; }
		if ( k1 < -k_range ) { k1 = -k_range; }
		if ( k0 > +k_range ) { k0 = +k_range; }
		if ( k1 > +k_range ) { k1 = +k_range; }
		ir0 = (int) ( (k0 + k_range) / dk );
		ir1 = (int) ( (k1 + k_range) / dk );

		size_t Nk = 1.0 + ( k1 - k0 ) / dk;
		step_r = floor( Nk / config->getParam("E_samples")->value );
		if ( step_r < 1 ) { step_r = 1; }
		if ( step_r > Nk ) { step_r = Nk; }
		num_r = 1 + Nk / step_r;
		config->getParam("E_samples")->value = num_r;

		psi_fft = (fftw_complex*) fftw_malloc( sizeof(fftw_complex) * N );
		plan = fftw_plan_dft_1d( N, psi_fft, psi_fft, FFTW_FORWARD, FFTW_ESTIMATE );
	}

	// prepare temporal downsampling
	t_range = config->getParam("t_range")->value / CONV_au_fs;
	dt = config->getParam("dt")->value / CONV_au_fs;
	t0 = ( config->getParam("obs_t0")->textual )  ?  0        :  config->getParam("obs_t0")->value / CONV_au_fs;
	t1 = ( config->getParam("obs_t1")->textual )  ?  t_range  :  config->getParam("obs_t1")->value / CONV_au_fs;
	size_t Nt = 1.0 + ( t1 - t0 ) / dt;
	step_t = floor( Nt / config->getParam("t_samples")->value );
	if ( step_t < 1 ) { step_t = 1; }
	if ( step_t > Nt ) { step_t = Nt; }
	num_t = 1 + Nt / step_t;
	config->getParam("t_samples")->value = num_t;

	valrec.resize( num_r );
	if ( rel_change ) {
		valrec_prev.resize( num_r );
	}
	rel_change_ready = false;

	boinc_resolve_filename_s( config->getParam("OUTFILE")->text.c_str(), filename );
	FILE* f = boinc_fopen( filename.c_str(), "w" );  // delete previous snapshot file
	fclose( f );
}

void Obs_Snapshot_WF::reinitialize( Conf_Module* config, vector<Module*> dependencies )
{
	GET_LOGGER( "liee.Obs_Snapshot_WF" );
	// save the actual number of samples again
	config->getParam("t_samples")->value = num_t;

	if ( do_fourier ) {
		config->getParam("k_samples")->value = num_r;
		psi_fft = (fftw_complex*) fftw_malloc( sizeof(fftw_complex) * N );
		plan = fftw_plan_dft_1d( N, psi_fft, psi_fft, FFTW_FORWARD, FFTW_ESTIMATE );
	} else {
		config->getParam("r_samples")->value = num_r;
	}
	valrec.resize( num_r );
	if ( rel_change ) {
		valrec_prev.resize( num_r );
	}
	rel_change_ready = false;

	// count the lines actually written to outfile (there might have been some appended between checkpoint-writing and program-exit)
	FILE* f = boinc_fopen( filename.c_str(), "r" );
	foundLns = 0;
	int ch;
	while ( EOF != ( ch = fgetc(f) ) ) {
		if ( ch == '\n' ) { ++foundLns; }
	}
	fclose( f );
}

void Obs_Snapshot_WF::estimate_effort( Conf_Module* config, double & flops, double & ram, double & disk )
{
	double N = config->getParam("r_range")->value / config->getParam("dr")->value;
	double Nt = config->getParam("t_range")->value / config->getParam("dt")->value;
	double Nr = config->getParam("r_samples")->value;
	int samples = (int) config->getParam("t_samples")->value;
	bool sqr = config->getParam("square")->text.compare("true") == 0;

	flops += abs( 8 * Nt + samples * N * 20 );  // called Nt times without saving + sample times with saving
	ram += 1024;
	double dsk = (samples + 2) * Nr * 2 * 23 * 1.5;  // factor: complex:2, double2str:23, outfile duplicated in gz archive: 1.5
	if ( sqr ) dsk /= 2;
	disk += abs( dsk );
}

void Obs_Snapshot_WF::observe( Module* state )
{
	Solver* s = dynamic_cast<Solver*>( state );
	if ( s->t < t0  ||  s->t > t1 ) return;  // not in observation window
	if ( counter++ % step_t != 0  &&  !rel_change_ready ) return;  // down-sample
	if ( foundLns >= ++writtenLns ) return;  // already written previously

	if ( do_fourier ) {
		for ( size_t i = 0; i < N; i++ ) {
			psi_fft[i][0] = real( s->psi[i] );
			psi_fft[i][1] = imag( s->psi[i] );
		}
		fftw_execute( plan );
	}

	double fac = 1.0;
	if ( do_normalize  &&  !rel_change  &&  !do_fourier ) {
		fac = 1.0 / s->integrate_psi_sqr();
	}

	for ( size_t i = ir0, j = 0; i <= ir1 && j < valrec.size(); i += step_r, j++ ) {
		if ( do_fourier ) {
			size_t i_ = ( i + N/2 ) % N;  // ring-looping necessary because fftw places positive frequencies in the first half before the negative frequencies in output-array.
			valrec[j] = dcmplx( psi_fft[i_][0], psi_fft[i_][1] );
		} else {
			if ( i >= s->Nr ) {
				LOG_WARN( "unexpectedly ran out of bounds in snapshot observer (N != s->Nr) ?" );
				break;
			}
			valrec[j] = fac * s->psi[i];
		}
	}

	if ( rel_change &&  !rel_change_ready ) {
		// remember last state for calculation of deviation
		valrec_prev = valrec;
		rel_change_ready = true;
		return;
	}

	FILE *file;
	file = boinc_fopen( filename.c_str(), "a" );
	for ( size_t i = 0; i < valrec.size(); i++ ) {
		dcmplx value = rel_change  ?  ( valrec[i] - valrec_prev[i] )  :  valrec[i];
		if ( do_square ) {
			value = value * conj( value );
			fprintf( file, format.c_str(), real( value ) );
		}
		else {
			fprintf( file, format.c_str(), real( value ), imag( value ) );
		}
	}
	fprintf( file, "\n" );
	fclose( file );

	if ( rel_change ) { rel_change_ready = false; }
}

//-------------------------------------------------------------------------------------------------------------

void Obs_Wigner_Distribution::initialize( Conf_Module* config, vector<Module*> dependencies )
{
	GET_LOGGER( "liee.Obs_Wigner_Distribution" );
	counter = 0;
	writtenFrames = 0;
	foundFrames = 0;

	// prepare temporal downsampling
	t_range = config->getParam("t_range")->value / CONV_au_fs;
	dt = config->getParam("dt")->value / CONV_au_fs;
	t0 = ( config->getParam("obs_t0")->textual )  ?  0        :  config->getParam("obs_t0")->value / CONV_au_fs;
	t1 = ( config->getParam("obs_t1")->textual )  ?  t_range  :  config->getParam("obs_t1")->value / CONV_au_fs;
	size_t Nt = 1.0 + ( t1 - t0 ) / dt;
	step_t = floor( Nt / config->getParam("t_samples")->value );
	if ( step_t < 1 ) { step_t = 1; }
	if ( step_t > Nt ) { step_t = Nt; }
	// save the actual number of samples
	num_t = 1 + Nt / step_t;
	config->getParam("t_samples")->value = num_t;

	double r_start = 0;
	for ( size_t i = 0; i < dependencies.size(); i++ ) {
		if ( dependencies[i]->type.compare("potential") == 0 ) {
			Potential* pot = dynamic_cast<Potential*>( dependencies[i] );
			r_start = pot->get_r_start();
		}
	}
	double r_end = r_start + config->getParam("r_range")->value / CONV_au_nm;
	double dr = config->getParam("dr")->value / CONV_au_nm;
	double kmax = CONST_PI / dr;
	r0 = ( config->getParam("obs_r0")->textual )  ?  r_start  :  config->getParam("obs_r0")->value / CONV_au_nm;  //TODO regard simulation bounds
	r1 = ( config->getParam("obs_r1")->textual )  ?  r_end    :  config->getParam("obs_r1")->value / CONV_au_nm;
	k0 = ( config->getParam("obs_k0")->textual )  ?  -kmax    :  config->getParam("obs_k0")->value * CONV_au_nm;
	k1 = ( config->getParam("obs_k1")->textual )  ?  +kmax    :  config->getParam("obs_k1")->value * CONV_au_nm;

	num_r = (int)config->getParam("r_samples")->value;
	num_k = (int)config->getParam("k_samples")->value;

	stringstream ss;
	int digits = int( config->getParam("precision")->value );
	ss << "%1." << digits << "g\t";
	format = ss.str();

	boinc_resolve_filename_s( config->getParam("OUTFILE")->text.c_str(), filename );
	FILE* f = boinc_fopen( filename.c_str(), "w" );  // delete previous snapshot file
	fclose( f );
}

void Obs_Wigner_Distribution::reinitialize( Conf_Module* config, vector<Module*> dependencies )
{
	GET_LOGGER( "liee.Obs_Wigner_Distribution" );
	// save the actual number of samples again
	config->getParam("t_samples")->value = num_t;

	// count the frames actually written to outfile (there might have been some appended between checkpoint-writing and program-exit)
	foundFrames = 0;
	FILE* f = boinc_fopen( filename.c_str(), "r" );
	int ch;
	while ( EOF != ( ch = fgetc(f) ) ) {
		if ( ch == '#' ) { ++foundFrames; }
	}
	fclose( f );
}

void Obs_Wigner_Distribution::estimate_effort( Conf_Module* config, double & flops, double & ram, double & disk )
{
	//int Nrr = (int) round( config->getParam("r_range")->value / config->getParam("dr")->value );
	//int Ntt = (int) round( config->getParam("t_range")->value / config->getParam("dt")->value );
	//int Nr = (int)config->getParam("r_samples")->value;
	//int Nk = (int)config->getParam("k_samples")->value;
	//int Nt = (int)config->getParam("t_samples")->value;
	//TODO estimate
}

void Obs_Wigner_Distribution::observe( Module* state )
{
	Solver* s = dynamic_cast<Solver*>( state );
	if ( s->t < t0  ||  s->t > t1 ) return;  // not in observation window
	if ( counter++ % step_t != 0 ) return;  // down-sample
	if ( foundFrames >= ++writtenFrames ) return;  // already written previously

	double r_start = s->potential->get_r_start();
	double dr = ( r1 - r0 ) / ( num_r - 1 );
	double dk = ( k1 - k0 ) / ( num_k - 1 );
	dcmplx two_i = dcmplx( 0, 2 );

	FILE *file;
	file = boinc_fopen( filename.c_str(), "a" );
	fprintf( file, "#i=%d\n#t=%g\n#rmin=%g\n#kmin=%g\n#dr=%g\n#dk=%g\n#Nr=%d\n#Nk=%d\n",
			(int)( counter / step_t ), s->t*CONV_au_fs, r0*CONV_au_nm, k0/CONV_au_nm,
			dr*CONV_au_nm, dk/CONV_au_nm, (int)num_r, (int)num_k );

	for ( size_t ri = 0; ri < num_r; ri++ )
	{
		double r = r0 + ri * dr;
		size_t rii = (size_t) round( (r - r_start) / s->dr );
		size_t delta_max = min( rii, s->Nr - 1 - rii );
		if ( rii + delta_max >= s->Nr  ||  rii - delta_max < 0 ) {
			DEBUG_SHOW4("Wigner out of range!", rii, delta_max, s->Nr );
		}

		for ( size_t ki = 0; ki < num_k; ki++ )
		{
			dcmplx kx2i = two_i * ( k0 + ki * dk );
			//TODO implement bucket summator to reduce numeric errors from adding small and large values
			dcmplx sum = 0;
			for ( size_t delta_i = 0; delta_i <= delta_max; delta_i++ ) {
				double delta_r = s->dr * delta_i;
				sum += conj( s->psi[rii + delta_i] ) * s->psi[rii - delta_i] * exp( kx2i * delta_r );
			}
			fprintf( file, format.c_str(), real( sum / CONST_PI ) );
		}
		fprintf( file, "\n" );
	}

	fprintf( file, "\n" );
	fclose( file );
	// to split the file before gnu-plotting use: awk '/##/{n++}{print >"frame" n ".dat" }' phasespace.dat
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

	flops += abs( 8 * Nt + N * Nr * 8 );  // called Nt times without saving + t_sample times to integrate (psi psi*)
	ram += 1024;
	disk += 30 * N;
}

void Obs_Tunnel_Ratio::summarize( map<string, string> & results )
{
	// save after getting the last sample
	tunnel_ratio = 1.0  -  psi_sqr.back() / psi_sqr.front();
	FILE* f = boinc_fopen( filename.c_str(), "w" );
	double last = psi_sqr.front();
	for ( size_t i = 0; i < psi_sqr.size(); i++ )
	{
		double j = ( last - psi_sqr[i] ) / dt;
		last = psi_sqr[i];
		fprintf( f, "%1.6g\t%1.16g\t%1.16g\n", i * dt , psi_sqr[i], j );
	}
	fprintf( f, "\n" );
	fclose( f );

	results["tunnel_ratio"] = doub2str( tunnel_ratio );  //TODO boost conversion
	if ( is_objective ) {
		results["objective"] = doub2str( tunnel_ratio );
	}
}

void Obs_Tunnel_Ratio::observe( Module* state )
{
	if ( counter++ % step_t != 0 ) return;
	Solver* s = dynamic_cast<Solver*>( state );
	psi_sqr.push_back( s->integrate_psi_sqr( ra, rb ) );
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
	r_end = V->get_r_phys_end();

	//calculate E from curvature of WF: E = V - hbar^2/(2m Psi) d^2/dr^2 Psi
	double rmin = V->getPot_const()->get_Vmin_pos();
	double Vmin = V->getPot_const()->V(rmin);
	Vp0 = V->getPot_const()->V(dr);
	double r0 = V->get_r_start();
	int i_middle = (int)( (rmin - r0) / dr );
	E = Vmin - ( wf->psi[i_middle-1].real() - 2 * wf->psi[i_middle].real() + wf->psi[i_middle+1].real() ) / ( pow( dr, 2.0 ) * 2.0 * wf->psi[i_middle].real() );
	DEBUG_SHOW(E);
	g = 2.0 * sqrt( 2.0 );
	last_r2 = numeric_limits<double>::quiet_NaN();
	burst = false;

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
	config->getParam("t_samples")->value = t_samples;  // save the actual number of samples
}

void Obs_JWKB_Tunnel::reinitialize( Conf_Module* config, vector<Module*> dependencies )
{
	GET_LOGGER( "liee.Obs_JWKB_Tunnel" );
	config->getParam("t_samples")->value = t_samples;  // save the actual number of samples again
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

	flops += abs( 8 * Nt + N * Nr * 10 );  // called Nt times without saving + t_sample times to integrate over the squareroot of barrier hight
	ram += 1024;
	disk += 20 * N;
}

void Obs_JWKB_Tunnel::summarize( map<string, string> & results )
{
	sum_j = 0;
	BOOST_FOREACH( double x, j ) {
		if ( not isinf(x) ) {  // its no use to invalidate the aggregated J_tot just because some infinities have been picked up along the way
			sum_j += x;
		}
	}

	// save to file after getting the last sample (again some code duplication with tunnel observer)
	FILE* f = boinc_fopen( filename.c_str(), "w" );
	for ( size_t i = 0; i < j.size(); i++ )	{
		fprintf( f, "%1.6g\t%1.10g\t%1.6g\t%1.6g\t%1.8g\n", i * dt , j[i], rt[i].x, rt[i].y, A[i] );
	}
	fprintf( f, "\n" );
	fclose( f );

	results["jwkb_J"] = doub2str( sum_j ); //TODO boost conversion
	results["jwkb_burst"] = burst==true ? "true" : "false";
	if ( is_objective ) {
		results["objective"] = doub2str( sum_j );
	}
}

void Obs_JWKB_Tunnel::observe( Module* state )
{
	if ( counter++ % step_t != 0 ) return;
	Solver* s = dynamic_cast<Solver*>( state );

	boost::function<double(double)> deltaE;
	deltaE = boost::bind( &Potential::deltaV, V, _1, s->t, E );  // deltaE(r) is mapped to deltaV(r, s->t, E)

	double F = ( Vp0 - V->V( dr, s->t ) ).real() / dr;  // estimate actual field strength from the change in potentials height compared with t0
	double r2 = -E / F;  // assume constant F +++ assume V approaching 0 for r > 0
	if (F <= 0  ||  r2 > 2 * r_end ) {
		j.push_back( 0.0 );  // barrier up or barely down +++ assume slow varying F --> infinite barrier --> j=0
		rt.push_back( Point(0,0) );
		A.push_back( 6e66 );
	}
	else {
		// find left turning point starting at r==0
		double r1 = 0;
		double d = dr;
		double a = deltaE( r1 );
		if ( a < 0 ) { d = -dr; }  // r1 already inside barrier --> direction backwards
		while ( deltaE( r1 + d ) * a > 0 ) {
			if ( r2 + d > r_end ) {  // no barrier at all --> j=infinity
				j.push_back( numeric_limits<double>::infinity() );
				burst = true;
				rt.push_back( Point(0,0) );
				A.push_back( 6e66 );
				return;
			}
			d *= 2;
			if ( isinf( d ) ) { throw Except__Too_Far_Out( __LINE__ ); }  // this shouldn't happen
		}
		r1 = find_root( deltaE, r1, r1 + d, 1e-12 );

		//TODO resolve code duplication with Pot_const::get_outer_turningpoints() by writing a robust and general root-finding routine
		// find right turning point
		if ( abs( deltaE( r2 ) ) > abs( deltaE( last_r2 ) ) && not isnan(last_r2) ) {
			r2 = last_r2;  // last turning-point was better than new estimate
		}
		d = dr;
		a = deltaE( r2 );
		if ( a > 0 ) { d = -dr; }  // r2 already outside barrier --> direction backwards
		while ( deltaE( r2 + d ) * a > 0 ) {
			if ( r2 + d > r_end ) {
				j.push_back( 0.0 );  // barrier extends farther than simulation range --> j=0
				rt.push_back( Point(0,0) );
				A.push_back( 6e66 );
				return;
			}
			d *= 2;
			if ( isinf( d ) ) { throw Except__Too_Far_Out( __LINE__ ); }  // this shouldn't happen
		}
		r2 = find_root( deltaE, r2, r2 + d, 1e-12 );
		last_r2 = r2;

		// integrate: dr sqrt( V(r) - E ) |r1..r2
		d = ( r2 - r1 ) / ( N - 1 );
		vector<Point> samples;
		for ( double r = r1; r <= r2; r += d ) {
			samples.push_back( Point( r, sqrt( abs( deltaE( r ) ) ) ) );
		}

		double A = 0;
	#ifdef ALGLIB
		alglib::spline1dinterpolant spline = to_cubic_spline( samples );
		A = alglib::spline1dintegrate( spline, r2 );
	#else
		A = simple_integrate( samples, r1, r2 );
	#endif
		j.push_back( exp( -g * A ) );
		rt.push_back( Point( r1, r2 ) );
		this->A.push_back( A );
	}
}

} // namespace liee
