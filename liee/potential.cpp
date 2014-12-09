#include <cmath>

#include "filesys.h"
#include "boinc_api.h"

#include "specialfunctions.h"  // from alglib
#include "fftw3.h"

#include "boost/foreach.hpp"
#include "boost/lexical_cast.hpp"
#include "boost/math/special_functions/erf.hpp"
#include "boost/function.hpp"
#include "boost/bind.hpp"

#include "potential.hpp"
#include "optimizer.hpp"

using namespace std;
namespace liee {


double Pot_const::get_Vmin_pos()
{
	if (been_there) return r_min;

	double b = 0;
	double Vb = V(b);
	// bracket minimum at b==0 by pushing a and b further along the negative/positive axis
	double a = -1e-40;	while ( V(a) <= Vb && !isinf(a) ) { a *= 2; }
	double c = 1e-40;	while ( V(c) <= Vb && !isinf(c) ) { c *= 2; }

	boost::function<double(double)> funct;
	funct = boost::bind( &Pot_const::V, this, _1 );

	opti::Golden_Section_Search minimizer( 1e-12 );
	r_min = minimizer.minimize( funct, a, b, c );

	been_there = true;
	return r_min;
}

void Pot_const::get_outer_turningpoints( double E, double & leftmost, double & rightmost )
{
	boost::function<double(double)> deltaE;
	deltaE = boost::bind( &Pot_const::deltaV, this, _1, E );

	double E_V0 = deltaE( r_min );
	if ( E_V0 < 0 ) {
		// E must not be lower than minimum
		leftmost = numeric_limits<double>::quiet_NaN();
		rightmost = numeric_limits<double>::quiet_NaN();
		return;
	}

	vector<double> result;
	result.resize(2);
	for ( int i = 0; i <= 1 ; i++ )  // left and right turning-point
	{
		double d = pow( -1, double( i + 1 ) );  // -1 and 1

		while ( deltaE( r_min + d ) * E_V0 > 0 )   {
			d *= 2;
			if ( isinf( d ) ) {
				result[i] = numeric_limits<double>::quiet_NaN();
				goto skip_find_root;
			}
		}
		result[i] = find_root( deltaE, r_min, r_min + d, 1e-12 );
		skip_find_root:;
	}
	leftmost = result[0];
	rightmost = result[1];
}

//----------------------------------------------------------------------------------------------------------

void Pot_Round_Well_wImage::initialize( Conf_Module* config, vector<Module*> dependencies )
{
	GET_LOGGER( "liee.Module.Pot_Round_Well_wImage" );
	been_there = false;
	width = config->get_double("width") / CONV_au_nm;
	depth = config->get_double("depth") / CONV_au_eV;
	expo = config->get_double("boxness") / ( width / 2.0 );

	a = depth / cosh( expo * width / 2.0 );
	shift_cosh = log( (8.0 * depth + expo - sqrt( expo*expo + 16.0 * depth * expo ) ) / 4.0 / a ) / expo;
	shift_mirror = 1.0 / ( 4.0 * depth - 2.0 * a * exp( expo * shift_cosh ) );
}
void Pot_Experimental::initialize( Conf_Module* config, vector<Module*> dependencies )
{
	GET_LOGGER( "liee.Module.Experimental" );
	been_there = false;
	width = config->get_double("width") / CONV_au_nm;
	depth = config->get_double("depth") / CONV_au_eV;
	expo = config->get_double("boxness") / ( width / 2.0 );
	a = depth / cosh( expo * width / 2.0 );
	shift_cosh = log( (8.0 * depth + expo - sqrt( expo*expo + 16.0 * depth * expo ) ) / 4.0 / a ) / expo;
	shift_mirror = 1.0 / ( 4.0 * depth - 2.0 * a * exp( expo * shift_cosh ) );

	v_scale = config->get_double("v_scale");
	h_scale = config->get_double("h_scale");
	power = config->get_double("power");
	b.resize(6);

	b[0] = 2.3 * width;	// left-sigma / left-width
	b[1] = config->get_double("boxness") / ( b[0] / 2.0 );  // left-expo
	b[2] = depth / cosh( b[1] * b[0] / 2.0 );  //left-a
	/*
	b[0] = -2164.61;
	b[1] = -422.862;
	b[2] = -32.9888;
	b[3] = -1.27608;
	b[4] = -0.0244718;
	b[5] = -0.000186276;
	*/
}

void Pot_Round_Well_wImage::get_outer_turningpoints( const double E, double & leftmost, double & rightmost )
{
	// -delta_switch is the center of the box-potential
	double exponential = 1.0 / expo * log( 2.0 * (E + depth) / a );
	leftmost = -shift_cosh - exponential + shift_mirror;
	double mirror = -shift_mirror - 0.25 / E;
	if ( mirror > 0 ) {
		rightmost = mirror + shift_mirror;
	} else {
		rightmost = -shift_cosh + exponential + shift_mirror;
	}
}

double Pot_Round_Well_wImage::get_Vmin_pos()
{
	return -shift_cosh + shift_mirror;
}

/*!
 * The stationary part of the potential is a potential well with given
 * width and constant depth. For r > width the potential is governed by
 * the image charge effect.
*/
double Pot_Round_Well_wImage::V( double r )
{
	double r_ = r - shift_mirror;
	if ( r_ < 0 )
		return -depth + a * cosh( expo * (r_ + shift_cosh) );  // potential well (right)
	else
		return -0.25 / r ;  // mirror charge
}

double Pot_Experimental::V( double r )
{
	double r_ = r - shift_mirror;
	if ( r_ < 0 ) {
		double x = r_ + shift_cosh;
		double wall = -depth + a * cosh( expo * x );

		double a1 = 0.0001; // vertical scale of ROOT
		double a2 = 1.0;    // horizontal squeeze of ROOT
		double a3 = 0.0;    // shift ROOT rightwards (in ratio of width)
		double a4 = 1.2;    // exponent of r-power

		if ( (-x + a3 * width) < 0 ) {
			return wall;
		} else {
			double root = -depth + a1 * pow( a2 * ( -x + a3 * width ), a4 );

			double a5 = 0.25; // shift midpoint of weighting Gauss leftwards (in ratio of width) relative to the middle of the well
			double a6 = 0.1;  // right-sigma of weighting Gauss (in ratio of width)
			double a7 = 2.3;  // left-sigma of weighting Gauss (in ratio of width)


			double gx = x + a5 * width;
			if ( gx > 0 ) {
				gx = gx / ( a6 * width );
			} else {
				gx = gx / ( a7 * width );
				wall = -depth + b[2] * cosh( b[1] * x );  // use wider left exponential wall to avoid too sudden onset
			}
			double weight = exp( -pow( gx, 2.0 ) );
			return ( 1 - weight ) * wall  +  weight * root;
		}

		/*
		double poli = b[0];
		for ( size_t i = 1; i < b.size(); i++ ) {
			poli += b[i] * pow( (r * CONV_au_nm)-7.0, (double)i );  // evaluate polynomial well (left)
		}
		poli =  poli / CONV_au_eV;

		double weight = 0.5 * ( 1.0 + boost::math::erf<double>( 0.01 * (r_ + shift_cosh + 5*19) ) ); // manually adjusted
		double wall = -depth + a * cosh( expo * (r_ + shift_cosh) ); // for morphing former potential shape into the polynomial
		return weight * wall + ( 1 - weight ) * poli;
		*/
	}
	else
		return -0.25 / r ;  // mirror charge (right)
}

//----------------------------------------------------------------------------------------------------------

void Gaussian_Pulse::initialize( Conf_Module* config, vector<Module*> dependencies )
{
	t_ofs = config->get_double("t_center") / CONV_au_fs;
	fwhm = config->get_double("FWHM") / CONV_au_fs;
	ga = 2.0 * log( 2.0 ) / pow( this->fwhm, 2.0 );
	F0 = config->get_double("amplitude") / CONV_au_V * CONV_au_m ;
	omega0 = 2 * PI * 137.0 / ( config->get_double("wavelength") / CONV_au_nm );  // 2 pi c / lambda
	phi0 = config->get_double("ce_phase");
	theta0 = config->get_double("chirp") * CONV_au_fs * CONV_au_fs;
}

double Gaussian_Pulse::electric_field( double t )
{
	double t_ = t - t_ofs;
	double t2 = pow(t_, 2.0);
	double gauss = exp( -ga * t2 ) * cos( theta0 * t2 + omega0 * t_ + phi0 );
	return F0 * gauss;
}

//----------------------------------------------------------------------------------------------------------

void Potential::register_dependencies( vector<Module*> dependencies )
{
	bool not_well = true;
	for ( size_t i = 0; i < dependencies.size(); i++ ) {
		if ( dependencies[i]->type.compare("pot_const") == 0  &&  dependencies[i]->serial == reg_serials[0] ) {
			LOG_INFO( "found a static potential to use: " + dependencies[i]->name)
			well = dynamic_cast<Pot_const*>( dependencies[i] );
			not_well = false;
		}
		else if ( dependencies[i]->type.compare("pulse") == 0 ) {
			for( size_t id = 1; id < reg_serials.size(); id++ ) {
				if ( reg_serials[id] == dependencies[i]->serial ) {
					LOG_INFO( "found a laser pulse to use: " + dependencies[i]->name )
					pulses.push_back( dynamic_cast<Laser_Field*>( dependencies[i] ) );
				}
			}
		}
	}
	if (not_well) {
		LOG_FATAL("Required dependency not found. At least need one pot_const module.");
		exit(1);
	}
}

void Potential::initialize( Conf_Module* config, vector<Module*> dependencies )
{
	GET_LOGGER( "liee.Module.Potential" );
	// save the relevant module id's of constant potential and pulses and register them
	reg_serials.push_back( config->get_int("my_pot_const") );
	BOOST_FOREACH( double serial, config->get_array("my_pulses") ) {
		reg_serials.push_back( (int) serial );
	}
	register_dependencies( dependencies );

	r_range = config->get_double("r_range") / CONV_au_nm;
	double inner_cutoff = config->get_double("inner_cutoff") / CONV_au_eV;
	double unused;
	well->get_outer_turningpoints( inner_cutoff, r_start, unused );  // set r_start where the inner_cutoff-Energy is reached
	if ( isnan(r_start) ) {
		r_start = 0;
		LOG_INFO( "Potential remains below cut-off value --> r_start = 0" );
	} else {
		LOG_INFO( "r_start: " << r_start << "\t( @ " << inner_cutoff << ")" );
	}

	F_dc = -1 * config->get_double("F_dc") / (CONV_au_V / CONV_au_nm);
	gamma = config->get_double("near_amplf");
	double s = config->get_double("near_width") / CONV_au_nm;
	s2 = pow( s , 2.0);
	wcap  = config->get_double("wcap") / CONV_au_nm;
	r_cap = r_start + r_range - wcap;
	int_samples  = config->get_int("int_samples");
	t_on = config->get_double("t_DC_on") / CONV_au_fs;
	t_full = config->get_double("t_DC_full") / CONV_au_fs;
	ini_erf = config->get_bool("ini_erf_DC");
	deltaDC = config->get_double("delta_DC") / CONV_au_nm;
	deltaAC = config->get_double("delta_AC") / CONV_au_nm;

	// set spatial positions for the evaluation of the laser-pulse-field.
	Pulse_samples.resize( int_samples );
	dx_sample = r_range / (int_samples - 1);
	for( size_t i = 0; i < int_samples; i++ ) {
		Pulse_samples[i] = Point( r_start + i * dx_sample, 0.0 );
	}
	t_now = -1;
}

void Potential::reinitialize( Conf_Module* config, vector<Module*> dependencies )
{
	GET_LOGGER( "liee.Module.Pot_Gauss" );
	register_dependencies( dependencies );
	// set spatial positions for the evaluation of the laser-pulse-field.
	Pulse_samples.resize( int_samples );
	dx_sample = r_range / (int_samples - 1);
	for( size_t i = 0; i < int_samples; i++ ) {
		Pulse_samples[i] = Point( r_start + i * dx_sample, 0.0 );
	}
	t_now = -1;
}

void Potential::estimate_effort( Conf_Module* config, double & flops, double & ram, double & disk )
{
	double N = config->get_double("r_range") / config->get_double("dr");
	flops += abs( 23 * N );
	ram += sizeof( Potential );
	disk += 0;
}

void Potential::summarize( map<string, string> & results )
{
	results["r_start"] = doub2str( r_start * CONV_au_nm );
}

dcmplx Potential::V( double r, double t )
{
	double z = (r - r_cap ) / wcap;  // CAP coordinate
	if (z > 0.0 && z <= 1.0)  // r is in CAP territory -> all other potential contributions kept constatnt
	{
		return dcmplx( well->V( r_cap ) +  V_Fdc( r_cap, t ) + V_pulse( r_cap, t ), V_cap( z ) );
	}
	else {
		return dcmplx( well->V( r ) +  V_Fdc( r, t ) + V_pulse( r, t ), 0.0 );
	}
}

dcmplx Potential::V_indexed( size_t ri, double t )
{
	// return without Vdc before it is switched on
	if ( t < t_on ) { return cache_Vwell[ri] + V_pulse( grid_r[ri], t ); }

	// while ramping up V_dc, it has to be re-calculated for each time-step
	if ( t <= t_full ) {
		return cache_Vwell[ri] + V_Fdc( grid_r[ri], t ) + V_pulse( grid_r[ri], t );
	}

	// otherwise use fully charged constant potential from cache
	return cache_Vconst[ri] + V_pulse( grid_r[ri], t );
}

void Potential::set_grid( double dr, size_t N )
{
	grid_dr = dr;
	grid_r.clear();
	cache_Vwell.clear();
	cache_Vconst.clear();
	for ( size_t i = 0; i < N; i++ ) {
		double r = r_start + i * dr;
		if ( r > r_cap ) { r = 1.000000001 * r_cap; }  // in the CAP-range, all normal contributions to the potential are kept constant at their value at r_cap
		grid_r.push_back( r );

		double z = ( (r_start + i * dr) - r_cap ) / wcap;
		if (z > 0.0  &&  z <= 1.0)  // r is in CAP territory
		{
			cache_Vconst.push_back( dcmplx( well->V( r_cap ) + V_Fdc( r_cap, t_full ), V_cap(z) ) );
			cache_Vwell.push_back( dcmplx( well->V( r_cap ), V_cap(z) ) );
		}
		else {
			cache_Vconst.push_back( dcmplx( well->V(r) + V_Fdc( r, t_full ), 0.0 ) );
			cache_Vwell.push_back( dcmplx( well->V(r), 0.0 ) );
		}
	}
}

void Potential::set_grid_CN( double dr, double dt, size_t N )
{
	set_grid( dr, N );
	dcmplx c = dcmplx( 0.5, dt / (4.0 * SQR(dr)) )  /  dcmplx( 0.0, dt / 4.0 );
	for ( size_t i = 0; i < N; i++ ) {
		cache_Vwell[i]  += c;
		cache_Vconst[i] += c;
	}
}

/*!
 * Returns only the real part, e.g. without the complex absorbing potential.
 */
double Potential::Vr( double r, double t )
{
	return well->V( r ) +  V_Fdc( r, t ) + V_pulse( r, t );
}


/*!
 * Absorbing Boundary Condition with quartic imaginary potential.
 * If the width defined by parameter 'wcap' is sufficiently large,
 * then portions of the wafefunction entering the absorbing region,
 * will be damped to zero to avoid reflection at the simulation edge.
 * @param z the coordinate spanning the CAP region, range: 0..1
*/
double Potential::V_cap( double z )
{
	return -6.13899 * pow(z, 4);
}

/*!
 * old version: V = int( dr' U / (k  r') = U / k * ln r', r' = r + r0
 * new version: constant F, linear V
 * new new version: exponential decay into negative range, adiabatic activation redux
 */
double Potential::V_Fdc( double r, double t )
{
	if ( t < t_on ) return 0;

	double V = ( r < 0 )  ?  F_dc * deltaDC * exp( r / deltaDC )  :  F_dc * ( deltaDC + r );  // exponential raise till r=0 then linear with r

	if ( t < t_full ) {
		// adiabatic activation
		double amp = (t - t_on) / (t_full - t_on);
		if ( ini_erf ) {
			amp = ( 1.0 + boost::math::erf<double>( 2 * exp(1.0) * (amp - 1.0) + exp(1.0) ) ) / 2.0;  // activation function: (erf(2e*t+e)+1)/2 |(-1 .. 0) --> (6e-5 .. 1 - 6e5)
		}
		return amp * V;
	}
	return V;
}

double Potential::F_pulse( double r, double t )
{
	double F = 0;
	double t_ = t + r / 137.0;
	for ( size_t i = 0; i < pulses.size(); i++ ) {
		F += pulses[i]->electric_field( t_ );
	}
	double amp = ( gamma > 0 )  ?  1.0 + gamma * s2 / (s2 + pow(r, 2.0))  :  1.0;
	if ( r >= 0 ) return amp * F;
	return amp * F * exp( r / deltaAC );
}

double Potential::V_pulse( double r, double t )
{
	if ( pulses.size() == 0 ) { return 0.0; }
	if ( t != t_now ) {
		// new time-step -> need to recalculate the integral samples using Simpson's rule
		Pulse_samples[0].y = 0;
		double fa = F_pulse( Pulse_samples[0].x, t );
		double fb, fm;
		double fac = (Pulse_samples[1].x - Pulse_samples[0].x) / 6.0;  // assume homogeneous grid

		for ( size_t i = 1; i < int_samples; i++ ) {
			fb = F_pulse( Pulse_samples[i].x, t );
			fm = F_pulse( 0.5*(Pulse_samples[i-1].x + Pulse_samples[i].x), t );
			Pulse_samples[i].y = Pulse_samples[i-1].y + fac * (fa + 4*fm + fb);
			fa = fb;
		}
		t_now = t;
	}

	double r_ = r - r_start;
	if ( r_ < 0 ) { return 0.0; }
	size_t i = floor( r_ / dx_sample );
	if (i >= int_samples) {
		return Pulse_samples[int_samples-1].y;
	}
	// linear interpolation
	float weight = r_ / dx_sample - i;
	return (1 - weight) * Pulse_samples[i].y  +  weight * Pulse_samples[i+1].y;
}

//------------------------------------------------------------------------------------------------------------

double Spatial_Light_Modificator::electric_field( double t )
{
	return 0;
}

void Spatial_Light_Modificator::register_dependencies( vector<Module*> dependencies )
{
	for ( size_t i = 0; i < dependencies.size(); i++ ) {
		if ( dependencies[i]->type.compare("pulse") == 0 ) {
			for( size_t id = 1; id < pulse_serials.size(); id++ ) {
				if ( pulse_serials[id] == dependencies[i]->serial ) {
					LOG_INFO( "found a laser pulse to use: " + dependencies[i]->name )
					pulses.push_back( dynamic_cast<Laser_Field*>( dependencies[i] ) );
				}
			}
		}
	}
}

void Spatial_Light_Modificator::initialize( Conf_Module* config, vector<Module*> dependencies )
{
	GET_LOGGER( "liee.Module.Pot_GaussSLM" );
	BOOST_FOREACH( double serial, config->get_array("orig_pulses") ) {
		pulse_serials.push_back( (int) serial );
	}
	register_dependencies( dependencies );

	//TODO improve the following sandbox code!
	/*
	double th = config->get_double("FWHM")->value / CONV_au_fs;
	double sw = 2.0 * sqrt( 2.0 * log(2.0) ) / th;
	st = 1.0 / sw;
	L = 10 * sw;
	DEBUG_SHOW4(omega0, st, sw, L);

	int N_slm = 1024;
	vector<double> slm_Amp( N_slm );
	vector<double> slm_Phi( N_slm );
	double T_sqr = 2.5 * th;
	for ( int k = 0; k < N_slm; k++ ) {
		double wk = omega0 - L/2  +  (k + 0.5) * L / N_slm;
		slm_Amp[k] = sin( (wk-omega0) * T_sqr ) / ( (wk-omega0) * T_sqr );
		slm_Phi[k] = 0;
	}

	// first we test an empty SLM to see if the FFT transforms the frequency-version of our gaussian pulse to the expected time-version
	double omegaC = 100.0/2 * omega0;
	double delta = 1.0 / 2.0 / (omegaC/2/PI);
	dcmplx z = dcmplx( 1.0/st/st, 2*theta0 );
	int N = 50 * (int) (ceil( omegaC * th) );	// T = 50 FWHM
	DEBUG_SHOW3( omegaC, delta, N );

	fftw_plan plan8, plan9;
	fftw_complex *pulse_t = (fftw_complex*) fftw_malloc( sizeof(fftw_complex) * N );
	fftw_complex *pulse_f = (fftw_complex*) fftw_malloc( sizeof(fftw_complex) * N );
	plan8 = fftw_plan_dft_1d( N, pulse_t, pulse_f, FFTW_FORWARD, FFTW_ESTIMATE );
	plan9 = fftw_plan_dft_1d( N, pulse_f, pulse_t, FFTW_BACKWARD, FFTW_ESTIMATE );

	pulse_f[0][0] = 0;
	pulse_f[0][1] = 0;
	FILE* f = boinc_fopen( "F_of_w.dat", "w" );
	fprintf( f, "%1.6g\t%1.6g\t%1.6g\n", 0.0, pulse_f[0][0], pulse_f[0][1] );
	for ( int i = 1; i <= N/2; i++ ) {
		double w = -i * omegaC / (N/2);
		int k = (int)floor( (w - (omega0 - L/2))  / (L / N_slm) -0.5);
		dcmplx F;
		if ( k < 0  ||  k >= N_slm ) {
			//transparent
			F = F0 / sqrt( z ) * exp( -pow( w - omega0, 2.0 ) / (2.0 * z) - dcmplx(0, phi0) );
		} else {
			F = slm_Amp[k] * F0 / sqrt( z ) * exp( -pow( w - omega0, 2.0 ) / (2.0 * z) - dcmplx(0, phi0) );
		}
		//F = F0 / sqrt( z ) * exp( -pow( w - omega0, 2.0 ) / (2.0 * z) - dcmplx(0, phi0) );
		pulse_f[N-i][0] = real(F);
		pulse_f[N-i][1] = imag(F);
		//w = w / CONV_au_fs;
		if ( k > 0  &&  k < N_slm ) {
			fprintf( f, "%1.6g\t%1.6g\t%1.6g\t%1.6g\t%d\n", w, pulse_f[N-i][0], pulse_f[N-i][1], slm_Amp[k], k );
		} else {
			fprintf( f, "%1.6g\t%1.6g\t%1.6g\t%s\t%d\n", w, pulse_f[N-i][0], pulse_f[N-i][1], "NaN", k );
		}

		w = i * omegaC / (N/2);
		k = (int)floor( (w - (omega0 - L/2)) / (L / N_slm) -0.5);
		if ( k < 0  ||  k >= N_slm ) {
			F = F0 / sqrt( z ) * exp( -pow( w - omega0, 2.0 ) / (2.0 * z) - dcmplx(0, phi0) );
		} else {
			F = slm_Amp[k] * F0 / sqrt( z ) * exp( -pow( w - omega0, 2.0 ) / (2.0 * z) - dcmplx(0, phi0) );
		}
		//F = F0 / sqrt( z ) * exp( -pow( w - omega0, 2.0 ) / (2.0 * z) - dcmplx(0, phi0) );
		pulse_f[i][0] = real(F);
		pulse_f[i][1] = imag(F);
		//w = w / CONV_au_fs;
		if ( k > 0  &&  k < N_slm ) {
			fprintf( f, "%1.6g\t%1.6g\t%1.6g\t%1.6g\t%d\n", w, pulse_f[i][0], pulse_f[i][1], slm_Amp[k], k );
		} else {
			fprintf( f, "%1.6g\t%1.6g\t%1.6g\t%s\t%d\n", w, pulse_f[i][0], pulse_f[i][1], "NaN", k );
		}
	}
	fclose( f );

	fftw_execute(plan9);

	f = boinc_fopen( "F_of_t.dat", "w" );
	for ( int i = 0; i < N; i++ ) {
		//double t = ( i - N / 2.0 ) / N * th * 10;
		//dcmplx F = F0 * exp( dcmplx( -0.5 * pow( t/st, 2.0 ), -(theta0*t*t + 0*omega0*t + phi0) ) );
		//pulse_t[i][0] = real(F);
		//pulse_t[i][1] = imag(F);
		double t = (i >= N/2) ? (i-N)*delta : i*delta;
		fprintf( f, "%1.6g\t%1.6g\t%1.6g\n", t, pulse_t[i][0], pulse_t[i][1] );
	}
	fclose( f );
	exit(0);
	*/
}

void Spatial_Light_Modificator::estimate_effort( Conf_Module* config, double & flops, double & ram, double & disk )
{
	//TODO implement
}

void Spatial_Light_Modificator::summarize( map<string, string> & results )
{
	//TODO implement
}

//------------------------------------------------------------------------------------------------------------

/*!
 * parabolic potential
 */
double Pot_Harm_Oscillator::V( double r )
{
	return 0.5 * k * pow( r - shift, 2.0 );
}

void Pot_Harm_Oscillator::get_outer_turningpoints( const double E, double & leftmost, double & rightmost )
{
	leftmost  = shift - sqrt( 2.0 * E / k );
	rightmost = shift + sqrt( 2.0 * E / k );
}

double Pot_Harm_Oscillator::get_Vmin_pos()
{
	return shift;
}

void Pot_Harm_Oscillator::initialize( Conf_Module* config, vector<Module*> dependencies )
{
	GET_LOGGER( "liee.Module.Pot_Harm_Oscillator" );
	been_there = false;
	k = config->get_double("k") / 1e-30 * CONV_au_s * CONV_au_s;
	w = sqrt( k );
	shift = config->get_double("shift") / CONV_au_nm;
}

/*!
 * not normalised
 */
double Pot_Harm_Oscillator::analytic_eigenfunction( int n, double x )
{
	double H = alglib::hermitecalculate( n, sqrt( w ) * (x - shift) );
	double expo = exp( -0.5 * w * pow( x - shift, 2.0 ) );
	//double F0 = pow( w / PI, 0.25 ) / sqrt( pow( 2, n ) / n! );
	return H * expo;
}

//------------------------------------------------------------------------------------------------------------

double Pot_Piecewise::V( double r )
{
	if ( r < segs.front().start_r ) {
		r = segs.front().start_r;
	}
	if ( r > segs.back().end_r ) {
		r = segs.back().end_r;
	}
	for ( size_t i = 0; i < segs.size(); i++ ) {
		if ( r >= segs[i].start_r  &&  r <= segs[i].end_r ) {
			if ( segs[i].sign == 0 ) {
				return segs[i].m * r + segs[i].n;  // line segment
			}
			if ( SQR( segs[i].radius ) - SQR( r - segs[i].cr ) < 0.0 ) {
				return segs[i].cV;
			}
			return segs[i].cV + segs[i].sign * sqrt( SQR( segs[i].radius ) - SQR( r - segs[i].cr ) ) / scale_y;  // ellipse segment
		}
	}
	LOG_ERROR( "Pot_Piecewise is not continuous at r = " << r )
	exit(1);
}

void Pot_Piecewise::get_outer_turningpoints( const double E, double & leftmost, double & rightmost )
{
	//TODO adapt to rounding procedure
	leftmost  = numeric_limits<double>::quiet_NaN();
	rightmost = numeric_limits<double>::quiet_NaN();

	for ( size_t i = 0; i < X.size()-1; i++ ) {
		if ( (E >= Y[i] && E < Y[i+1]) || (E <= Y[i] && E > Y[i+1]) )
		{
			double turningpoint = X[i] + ( E - Y[i] ) / ( Y[i+1] - Y[i] ) * ( X[i+1] - X[i] );

			if ( isnan( leftmost ) ) {
				leftmost = turningpoint;
			}
			rightmost = turningpoint;
		}
	}
}

double Pot_Piecewise::get_Vmin_pos()
{
	//TODO adapt to rounding procedure
	size_t i_min = 0;
	for ( size_t i = 1; i < X.size(); i++ ) {
		if ( Y[i] < Y[i_min] ) {
			i_min = i;
		}
	}
	return X[i_min];
}


void Pot_Piecewise::Segment::line( double x1, double x2, double y1, double y2 )
{
	sign = 0;
	start_r = x1;
	end_r = x2;
	m = (y2 - y1) / (x2 - x1);
	n = y1 - m * x1;
	ortho_X = -m  / sqrt( 1.0 + SQR(m) );
	ortho_Y = 1.0 / sqrt( 1.0 + SQR(m) );
}

void Pot_Piecewise::initialize( Conf_Module* config, vector<Module*> dependencies )
{
	GET_LOGGER( "liee.Module.Pot_Piecewise" );
	double epsilon = config->get_double("epsilon") / CONV_au_nm;
	rounding = config->get_array("ellipse")[0] / CONV_au_nm;
	scale_y = rounding / ( config->get_array("ellipse")[1] / CONV_au_eV );
	X = config->get_array("r_list");
	Y = config->get_array("V_list");
	if ( X.size() != Y.size() || X.size() == 0 ) {
		LOG_ERROR( "r_list and V_list are required to have the same size with at least 1 element." )
		exit(1);
	}
	if ( X.size() == 1 ) {  // constant potential
		X.push_back( X[0] + config->get_double("r_range") );
		Y.push_back( Y[0] );
	}

	for ( size_t i = 0; i < X.size(); i++ ) {
		X[i] /= CONV_au_nm;
		Y[i] /= CONV_au_eV;
		if ( i > 0  &&  X[i-1] > X[i] ) {
			LOG_ERROR( "r_list is required to be sorted in ascending order. exiting." )
			exit(1);
		}
		if ( i > 0  &&  X[i-1] == X[i] ) {
			X[i] += epsilon;
		}
	}

	Segment first;
	first.line( X[0], X[1], Y[0], Y[1] );
	segs.push_back( first );

	for ( size_t i = 1; i < X.size()-1; i++ )
	{
		Segment left, right;
		left = segs.back();
		right.line( X[i], X[i+1], Y[i], Y[i+1] );

		if ( rounding > 0  &&  left.m != right.m ) {
			// note: tangential-circle calculation is done in a y-stretched system in order to
			// work with a circle instead of an ellipse. The x-positions are valid in both systems,
			// while the circle centre-V is scaled back to the actual ellipse-centre
			// (V(r) draws the ellipse according to scale_y)
			for ( double rr = rounding; rr > 1e-4 * rounding; rr /= 1.5 )
			{
				Segment sleft, sright, circle;
				double S = scale_y;
				sleft.line( X[i-1], X[i], S*Y[i-1], S*Y[i] );
				sright.line( X[i], X[i+1], S*Y[i], S*Y[i+1] );

				circle.sign = ( left.m > right.m )  ?  +1  : -1;  // sign=+1: circle under the roof-like kink, positive "square-root"-circle is touching tangents and opposite for sign=-1
				double dist = -circle.sign * rr;  // the possible centre-positions for the tangential circle are on a parallel line, distanced 1 rounding-radius in orthogonal direction, either above or below
				double shifted_n_l = (S*Y[i] + dist * sleft.ortho_Y)  - sleft.m  * (X[i] + dist * sleft.ortho_X);  // ordinate of the left-side parallel
				double shifted_n_r = (S*Y[i] + dist * sright.ortho_Y) - sright.m * (X[i] + dist * sright.ortho_X);
				circle.cr = (shifted_n_l - shifted_n_r) / (sright.m - sleft.m);  // circle centre at the intersection of the two parallels
				circle.cV = ( sleft.m * circle.cr + shifted_n_l ) / S;
				circle.start_r = circle.cr - dist * sleft.ortho_X;  // from circle centre back on the normal-vector to the respective lines
				circle.end_r   = circle.cr - dist * sright.ortho_X;
				circle.radius = rr;

				// only accept rounding circle if it does not exceed the middle of neighbouring segments
				if (    circle.start_r > 0.5 * (segs.back().start_r + segs.back().end_r)
				     && circle.end_r   < 0.5 * (right.start_r + right.end_r) )
				{
					segs.back().end_r = circle.start_r;  // reduce left segment extend
					right.start_r = circle.end_r;
					segs.push_back( circle );
					break;
				} // else continue loop with a smaller radius
			}
		}
		segs.push_back( right );
	}
}

void Pot_Piecewise::reinitialize( Conf_Module* config, vector<Module*> dependencies )
{
	GET_LOGGER( "liee.Module.Pot_Piecewise" );
}

//------------------------------------------------------------------------------------------------------------

/*!
 * the potential
 */
double Chulkov_Image_Potential::V( double r )
{
	r = abs( r );
	double V = 0;

	if ( r <= D ) {
		V = A10 + A1 * cos( 2.0 * PI / as * r );
	}
	else if ( r > D && r <= z1 ) {
		V = -A20 + A2 * cos( beta * (r - D) );
	}
	else if ( r > z1 && r <= zim ) {
		V = A3 * exp( -alpha * (r - z1) );
	}
	else {  // if ( r > zim )
		V = ( exp( -lambda * (r - zim) ) - 1.0 ) / ( 4.0 * (r - zim) );
	}

	return V;
}

void Chulkov_Image_Potential::get_outer_turningpoints( const double E, double & leftmost, double & rightmost )
{
	double E4 = E * 4;
	rightmost = +zim + lambert_w( exp( lambda / E4 ) * lambda / E4, true ) / lambda  -  1 / E4;
	if ( rightmost-zim < as/1e9 ) {
		// not in the r>zim zone, try z1<r<=zim
		rightmost = z1 - log( E / A3 ) / alpha;
		if ( rightmost-z1 < as/1e9 ) {
			// not in the z1<r<=zim zone, try D<r<=z1
			double offset = acos( (E + A20) / A2 );
			if ( not isnan(offset) ) {
				double r_last = z1;
				for ( int k = 0; k < 3 * nlay; k++ ) {  //TODO(perf.) instead of stupid iteration from the beginning, start from the end
					for ( int sign = -1; sign <= 1; sign += 2 ) {
						double r = D + ( sign * offset + k * 2 * PI ) / beta;
						if ( r > z1 ) {
							// now exceeded target zone --> r_last was the wanted solution on the final up-slope
							rightmost = r_last;
							leftmost = -rightmost;
							return;
						}
						r_last = r;
					}
				}
			}
			else {
				// in the r<D zone
				//TODO(code dupl.) with previous case
				double offset = acos( (E - A10) / A1 );
				if ( not isnan(offset) ) {
					double r_last = 0;
					for ( int k = 0; k < 3 * nlay; k++ ) {  //TODO(perf.) instead of stupid iteration from the beginning, start from the end
						for ( int sign = -1; sign <= 1; sign += 2 ) {
							double r = ( sign * offset + k * 2 * PI ) * as /2.0 /PI;
							if ( r > D ) {
								// now exceeded target zone --> r_last was the wanted solution on the final up-slope
								rightmost = r_last;
								leftmost = -rightmost;
								return;
							}
							r_last = r;
						}
					}
				} else {}  //TODO error-handling
			}
		}
	}
	leftmost = -rightmost;
}

double Chulkov_Image_Potential::get_Vmin_pos()
{
	return as / 2;
}

void Chulkov_Image_Potential::initialize( Conf_Module* config, vector<Module*> dependencies )
{
	GET_LOGGER( "liee.Module.Chulkov_Image_Potential" );
	been_there = false;
	//this->r_end = config->get_double("r_end") * 1e-9 / CONV_au_m;

	// Das hier und die Funktion im Anhang erzeugen das Potential für Beryllium(0001)
	// mit 41*2 Cos-Perioden im inneren und einer Wandstärke von 40*2 Perioden.
	as = 3.387;
	nlay = 41;
	//double dslab = 2 * as * nlay  +  40 * as;
	shift = 0;
	A10 = -18.750 / CONV_au_eV;
	A1 = 6.2 / CONV_au_eV;
	A2 = 4.1354 / CONV_au_eV;
	beta = 4.2630;
	// moved from function V to initialisation
	A20 = A2 - A10 - A1;
	D = nlay * as;
	z1 = 5 * PI / (4 * beta) + D;
	A3 = -A20 + A2 * cos( beta * (z1 - D) );
	alpha = A2 / A3 * beta * sin( beta * (z1 - D) );
	lambda = 2 * alpha;
	zim = z1 - 1 / alpha * log( -lambda / (4 * A3) );

	/*
	FILE* f = fopen( "initial_potential.dat", "w" );
	//double dr = r_end / 2000;
	//for ( double x = 0; x < r_end; x += dr ) {
	for ( double x = -dslab/2; x < dslab/2; x += dslab/10000 ) {
		fprintf( f, "%1.8g\t%1.8g\n", x, Vr( x, 0 ) );
	}
	fclose( f );
	*/
}

} //namespace
