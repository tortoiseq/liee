#include <fstream>
#include <iomanip>
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
	do_square = config->get_bool("square");
	do_fourier = config->get_bool("fourier");
	do_normalize = config->get_bool("normalize");
	rel_change = config->get_bool("rel_change");
	do_average = config->get_bool("average");
	if ( do_average ) { rel_change = false; }

	string digits = config->get_string("precision");  // untested if "precision" is integer
	format = "%1." + digits + "g\t";
	if ( not do_square ) {
		format = format + format;  // twice the format-string for real- and imaginary part
	}
	stringstream range_info;
	range_info << std::scientific << std::setprecision(8);

	double dr = config->get_double("dr") / CONV_au_nm;
	double r0 = config->get_double("obs_r0") / CONV_au_nm;
	double r1 = config->get_double("obs_r1") / CONV_au_nm;
	ir0 = (int) (r0 / dr);
	ir1 = (int) (r1 / dr);

	if ( not do_fourier ) {
		// prepare spatial downsampling
		size_t Nr = 1.0 + ( r1 - r0 ) / dr;
		step_r = floor( Nr / config->get_double("r_samples") );
		if ( step_r < 1 ) { step_r = 1; }
		if ( step_r > Nr ) { step_r = Nr; }
		num_r = (Nr % step_r)  ?  Nr / step_r + 1  :  Nr / step_r;  // ceil( Nr/step_r )
		config->set_int("r_samples", num_r);
		valrec.resize( num_r );
		range_info << "##\t" << "r0=" << r0 * CONV_au_nm << ";\n";
		range_info << "##\t" << "r1=" << r1 * CONV_au_nm << ";\n";
		range_info << "##\t" << "Nr=" << num_r << ";\n";
	} else {
		// prepare spectral downsampling
		double k_range = 0.5 / dr;  // previously: k_max = 1/(2 dr)
		k_range *= 2.0 * CONST_PI;  // values were off by 2pi
		double dk = 2.0 * k_range / ( (r1 - r0) / dr );
		double k0 = config->param_is_nan("obs_k0")  ?  -k_range  :  config->get_double("obs_k0") * CONV_au_nm;
		double k1 = config->param_is_nan("obs_k1")  ?  +k_range  :  config->get_double("obs_k1") * CONV_au_nm;
		if ( k0 < -k_range ) { k0 = -k_range; }
		if ( k1 < -k_range ) { k1 = -k_range; }
		if ( k0 > +k_range ) { k0 = +k_range; }
		if ( k1 > +k_range ) { k1 = +k_range; }
		//TODO check indices ik0, ik1! Experiments at low k-resolution indicate that expected k(ik) is actually k(ik-1), i.e. the index might be offset by 1
		ik0 = (int) ( (k0 + k_range) / dk );
		ik1 = (int) ( (k1 + k_range) / dk );

		size_t Nk = 1.0 + ik1 - ik0;
		step_k = floor( Nk / config->get_double("k_samples") );
		if ( step_k < 1 ) { step_k = 1; }
		if ( step_k > Nk ) { step_k = Nk; }
		num_k = (Nk % step_k)  ?  Nk / step_k + 1  :  Nk / step_k;
		config->set_int("k_samples", num_k);
		valrec.resize( num_k );
		range_info << "##\t" << "k0=" << k0 / CONV_au_nm << ";\n";
		range_info << "##\t" << "k1=" << k1 / CONV_au_nm << ";\n";
		range_info << "##\t" << "Nk=" << num_k << ";\n";

		Nfou = 1 + ir1 - ir0;
		psi_fft = (fftw_complex*) fftw_malloc( sizeof(fftw_complex) * Nfou );  //TODO free this memory when done
		plan = fftw_plan_dft_1d( Nfou, psi_fft, psi_fft, FFTW_FORWARD, FFTW_ESTIMATE );
	}

	// prepare temporal downsampling
	t_range = config->get_double("t_range") / CONV_au_fs;
	dt = config->get_double("dt") / CONV_au_fs;
	t0 = config->get_double("obs_t0") / CONV_au_fs;
	t1 = config->get_double("obs_t1") / CONV_au_fs;
	size_t Nt = 1.0 + ( t1 - t0 ) / dt;
	step_t = floor( Nt / config->get_double("t_samples") );
	if ( step_t < 1 ) { step_t = 1; }
	if ( step_t > Nt ) { step_t = Nt; }
	num_t = (Nt % step_t)  ?  Nt / step_t + 1  :  Nt / step_t;
	config->set_int("t_samples",  num_t);
	double v_g = num_t / (t1 - t0) * (r1 - r0);  // only relevant for spectrometer calibration by a runtime-correction-factor v/v_g {v_g being the velocity at which (on average) just one sample is taken while the wave-packet is passing the detector's length. Slower packets are oversampled by a factor of v_g/v. As long v remains below relativistic, v_g might be surpass the speed of light without causing problems.}

	range_info << "##\t" << "t0=" << t0 * CONV_au_fs << ";\n";
	range_info << "##\t" << "t1=" << t1 * CONV_au_fs << ";\n";
	range_info << "##\t" << "Nt=" << num_t << ";\n";
	range_info << "##\t" << "v_g=" << v_g * CONV_au_nm/CONV_au_fs << ";\n";

	if ( rel_change || do_average ) {
		valrec_prev.resize( valrec.size(), dcmplx(0.0, 0.0) );
	}
	rel_change_ready = false;

	boinc_resolve_filename_s( config->get_string("OUTFILE").c_str(), filename );
	FILE* f = boinc_fopen( filename.c_str(), "w" );  // delete previous snapshot file
	fprintf( f, "%s", range_info.str().c_str() );  // write range-information as comments
	fclose( f );
}

void Obs_Snapshot_WF::reinitialize( Conf_Module* config, vector<Module*> dependencies )
{
	GET_LOGGER( "liee.Obs_Snapshot_WF" );
	// save the actual number of samples again
	config->set_int("t_samples", num_t);

	if ( do_fourier ) {
		config->set_int("k_samples", num_k);
		psi_fft = (fftw_complex*) fftw_malloc( sizeof(fftw_complex) * Nfou );
		plan = fftw_plan_dft_1d( Nfou, psi_fft, psi_fft, FFTW_FORWARD, FFTW_ESTIMATE );
	} else {
		config->set_int("r_samples", num_r);
	}
	valrec.resize( num_r );
	if ( rel_change ) {
		valrec_prev.resize( num_r );
	}
	rel_change_ready = false;

	if (do_average) return;
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
	double N = config->get_double("r_range") / config->get_double("dr");
	double Nt = config->get_double("t_range") / config->get_double("dt");
	int Nr = config->get_int("r_samples");
	int samples = config->get_int("t_samples");
	bool sqr = config->get_bool("square");

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
		for ( size_t i = ir0; i <= ir1; i++ ) {
			psi_fft[i-ir0][0] = real( s->psi[i] );
			psi_fft[i-ir0][1] = imag( s->psi[i] );
		}
		fftw_execute( plan );
	}

	double fac = 1.0;
	if ( do_normalize  &&  !rel_change  &&  !do_fourier ) {
		fac = 1.0 / sqrt( s->integrate_psi_sqr() );
	}

	//TODO optional smoothing of data before sub-sampling
	if ( do_fourier ) {
		//fac = sqrt( s->dr / Nfou) );
		fac = s->dr / sqrt( 2.0 * CONST_PI );
		for ( size_t i = ik0, j = 0; i <= ik1 && j < valrec.size(); i += step_k, j++ ) {
			// ring-looping necessary because fftw places positive frequencies in the first half before the negative frequencies in output-array.
			size_t i_ = ( i + Nfou/2 ) % Nfou;
			valrec[j] = fac * dcmplx( psi_fft[i_][0], psi_fft[i_][1] );
		}
	}
	else for ( size_t i = ir0, j = 0; i <= ir1 && j < valrec.size(); i += step_r, j++ ) {
		valrec[j] = fac * s->psi[i];
	}

	if ( rel_change &&  !rel_change_ready ) {
		// remember last state for calculation of deviation
		valrec_prev = valrec;
		rel_change_ready = true;
		return;
	}
	if ( rel_change ) { rel_change_ready = false; }

	if ( do_average ) {
		// use valrec_prev to store the sum, (rel_change and do_average are exclusive)
		for ( size_t i = 0; i < valrec.size(); i++ ) {
			if ( do_square ) {
				valrec_prev[i] += valrec[i] * conj( valrec[i] );
			} else {
				valrec_prev[i] += valrec[i];
			}
		}
	}
	else {
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
	}
}

void Obs_Snapshot_WF::summarize( map<string, string> & results )
{
	if ( ! do_average ) return;

	// save after getting the last sample
	FILE* f = boinc_fopen( filename.c_str(), "a" );
	for ( size_t i = 0; i < valrec_prev.size(); i++ ) {
		if ( do_square ) {
			fprintf( f, "%1.16g\n", real( valrec_prev[i] )/num_t );
		}
		else {
			fprintf( f, "%1.16g\t%1.16g\n", real( valrec_prev[i] )/num_t, imag( valrec_prev[i] )/num_t );
		}
	}
	fclose( f );
}
//-------------------------------------------------------------------------------------------------------------

void Obs_Wigner_Distribution::initialize( Conf_Module* config, vector<Module*> dependencies )
{
	GET_LOGGER( "liee.Obs_Wigner_Distribution" );
	counter = 0;
	writtenFrames = 0;
	foundFrames = 0;

	// prepare temporal downsampling
	t_range = config->get_double("t_range") / CONV_au_fs;
	dt = config->get_double("dt") / CONV_au_fs;
	t0 = config->get_double("obs_t0") / CONV_au_fs;
	t1 = config->get_double("obs_t1") / CONV_au_fs;
	size_t Nt = 1.0 + ( t1 - t0 ) / dt;
	step_t = floor( Nt / config->get_double("t_samples") );
	if ( step_t < 1 ) { step_t = 1; }
	if ( step_t > Nt ) { step_t = Nt; }
	// save the actual number of samples
	num_t = (Nt % step_t)  ?  Nt / step_t + 1  :  Nt / step_t;
	config->set_int("t_samples", num_t);

	double r_start = 0;
	for ( size_t i = 0; i < dependencies.size(); i++ ) {
		if ( dependencies[i]->type.compare("potential") == 0 ) {
			Potential* pot = dynamic_cast<Potential*>( dependencies[i] );
			r_start = pot->get_r_start();
		}
	}
	double r_end = r_start + config->get_double("r_range") / CONV_au_nm;
	double dr = config->get_double("dr") / CONV_au_nm;
	double kmax = CONST_PI / dr;
	r0 = ( config->param_is_nan("obs_r0") )  ?  r_start  :  config->get_double("obs_r0") / CONV_au_nm;  //TODO regard simulation bounds
	r1 = ( config->param_is_nan("obs_r1") )  ?  r_end    :  config->get_double("obs_r1") / CONV_au_nm;
	k0 = ( config->param_is_nan("obs_k0") )  ?  -kmax    :  config->get_double("obs_k0") * CONV_au_nm;
	k1 = ( config->param_is_nan("obs_k1") )  ?  +kmax    :  config->get_double("obs_k1") * CONV_au_nm;

	num_r = config->get_int("r_samples");
	num_k = config->get_int("k_samples");

	stringstream range_info;
	range_info << std::scientific << std::setprecision(8);
	range_info << "##\t" << "t0=" << t0 * CONV_au_fs << ";\n";
	range_info << "##\t" << "t1=" << t1 * CONV_au_fs << ";\n";
	range_info << "##\t" << "Nt=" << num_t << ";\n";
	range_info << "##\t" << "r0=" << r0 * CONV_au_nm << ";\n";
	range_info << "##\t" << "r1=" << r1 * CONV_au_nm << ";\n";
	range_info << "##\t" << "Nr=" << num_r << ";\n";
	range_info << "##\t" << "k0=" << k0 / CONV_au_nm << ";\n";
	range_info << "##\t" << "k1=" << k1 / CONV_au_nm << ";\n";
	range_info << "##\t" << "Nk=" << num_k << ";\n";

	format = "%1." + config->get_string("precision") + "g\t";

	boinc_resolve_filename_s( config->get_string("OUTFILE").c_str(), filename );
	FILE* f = boinc_fopen( filename.c_str(), "w" );  // delete previous snapshot file
	fprintf( f, "%s", range_info.str().c_str() );  // write range-information as comments
	fclose( f );
}

void Obs_Wigner_Distribution::reinitialize( Conf_Module* config, vector<Module*> dependencies )
{
	GET_LOGGER( "liee.Obs_Wigner_Distribution" );
	// save the actual number of samples again
	config->set_int("t_samples", num_t);

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
			// TODO use pairwise summation or Kahan summation to reduce numeric errors from adding small and large values
			double sum = 0;
			for ( size_t delta_i = 0; delta_i <= delta_max; delta_i++ ) {
				double delta_r = s->dr * delta_i;
				sum += real( conj( s->psi[rii + delta_i] ) * s->psi[rii - delta_i] * exp( kx2i * delta_r ) );
			}
			fprintf( file, format.c_str(), sum / CONST_PI );
		}
		fprintf( file, "\n" );
	}
	fprintf( file, "\n" );
	fclose( file );
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
	filename = config->get_string("OUTFILE");
	is_objective = config->get_bool("is_objective");
	dr = config->get_double("dr") / CONV_au_nm;
	N = config->get_int("int_samples");
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
	//TODO lots of code duplication with other observers
	t_range = config->get_double("t_range") / CONV_au_fs;
	dt = config->get_double("dt") / CONV_au_fs;
	counter = 0;
	int Nt = 1.0 + t_range / dt;
	step_t = floor( Nt / config->get_double("t_samples") );
	dt *= step_t;
	if ( step_t < 1 ) { step_t = 1; }
	if ( step_t > Nt ) { step_t = Nt; }
	num_t = (Nt % step_t)  ?  Nt / step_t + 1  :  Nt / step_t;
	config->set_int("t_samples", num_t);  // save the actual number of samples
}

void Obs_JWKB_Tunnel::reinitialize( Conf_Module* config, vector<Module*> dependencies )
{
	GET_LOGGER( "liee.Obs_JWKB_Tunnel" );
	config->set_int("t_samples", num_t);  // save the actual number of samples again
	for ( size_t i = 0; i < dependencies.size(); i++ ) {
		if ( dependencies[i]->type.compare( "potential" ) == 0 ) {
			V = dynamic_cast<Potential*>( dependencies[i] );
		}
	}
}

void Obs_JWKB_Tunnel::estimate_effort( Conf_Module* config, double & flops, double & ram, double & disk )
{
	double Nr = config->get_double("r_range") / config->get_double("dr");
	double Nt = config->get_double("t_range") / config->get_double("dt");
	int N =  config->get_int("t_samples");

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

//-------------------------------------------------------------------------------------------------------------

void Obs_Probability_Current::initialize( Conf_Module* config, vector<Module*> dependencies )
{
	GET_LOGGER( "liee.Probability_Current" );
	double r_range;
	for ( size_t i = 0; i < dependencies.size(); i++ ) {
		if ( dependencies[i]->type.compare( "potential" ) == 0 ) {
			r_range = dynamic_cast<Potential*>( dependencies[i] )->get_r_range();
		}
	}

	filename = config->get_string("OUTFILE");
	is_objective = config->get_bool("is_objective");
	r_detect = config->get_double("r_detect") / CONV_au_nm;
	if ( r_detect < 0 ) r_detect = 0;
	if ( r_detect > r_range ) { r_detect = r_range; }
	double dr = config->get_double("dr") / CONV_au_nm;
	ri = (int)( 0.5 + r_detect / dr );
	prefac = dcmplx( 0.0, 0.25 / dr );

	// temporal downsampling
	t_range = config->get_double("t_range") / CONV_au_fs;
	dt = config->get_double("dt") / CONV_au_fs;
	counter = 0;
	int Nt = 1.0 + t_range / dt;
	step_t = floor( Nt / config->get_double("t_samples") );
	dt *= step_t;
	if ( step_t < 1 ) { step_t = 1; }
	if ( step_t > Nt ) { step_t = Nt; }
	num_t = (Nt % step_t)  ?  Nt / step_t + 1  :  Nt / step_t;
	config->set_int("t_samples", num_t);  // save the actual number of samples
}

void Obs_Probability_Current::reinitialize( Conf_Module* config, vector<Module*> dependencies )
{
	GET_LOGGER( "liee.Probability_Current" );
	config->set_int("t_samples", num_t); // save the actual number of samples again
}

void Obs_Probability_Current::estimate_effort( Conf_Module* config, double & flops, double & ram, double & disk )
{
	double Nt = config->get_double("t_range") / config->get_double("dt");
	int N =  config->get_int("t_samples");

	flops += abs( 8 * Nt + N * 4 );  // called Nt times without saving + t_sample times to save a single complex number
	ram += 8 * N;
	disk += 30 * N;
}

void Obs_Probability_Current::summarize( map<string, string> & results )
{
	// save after getting the last sample
	FILE* f = boinc_fopen( filename.c_str(), "w" );
	double sumJ = 0.0;
	for ( size_t i = 0; i < j.size(); i++ )
	{
		sumJ += j[i];  // TODO use pairwise summation or Kahan summation
		fprintf( f, "%1.6g\t%1.16g\n", i * dt * CONV_au_fs , j[i] );
	}
	fprintf( f, "\n" );
	fclose( f );

	double J = sumJ * dt;
	results["sum_over_j(t)"] = doub2str( J );  //TODO boost conversion
	if ( is_objective ) {
		results["objective"] = doub2str( J );
	}
}

void Obs_Probability_Current::observe( Module* state )
{
	if ( counter++ % step_t != 0 ) return;
	Solver* s = dynamic_cast<Solver*>( state );
	dcmplx l = s->psi[ri-1];
	dcmplx m = s->psi[ri];
	dcmplx r = s->psi[ri+1];
	j.push_back( real( prefac * ( m * conj(r-l) - conj(m) * (r - l) ) ) );
}


} // namespace liee
