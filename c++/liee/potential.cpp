#include <cmath>

#include "filesys.h"
#include "boinc_api.h"

#include "../alglib/specialfunctions.h"

#include <fftw3.h>

#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>

#include "potential.hpp"

using namespace std;
namespace liee {


void Pot_Round_Well_wImage::initialize( Conf_Module* config, vector<Module*> dependencies )
{
	GET_LOGGER( "liee.Module.Pot_Round_Well_wImage" );
	width = config->getParam("width")->value / CONV_au_nm ;
	depth = config->getParam("depth")->value / CONV_au_eV;
	expo = config->getParam("boxness")->value / ( width / 2.0 );

	a = depth / cosh( expo * width / 2.0 );
	shift_cosh = log( (8.0 * depth + expo - sqrt( expo*expo + 16.0 * depth * expo ) ) / 4.0 / a ) / expo;
	shift_mirror = 1.0 / ( 4.0 * depth - 2.0 * a * exp( expo * shift_cosh ) );
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
inline double Pot_Round_Well_wImage::V( double r )
{
	double r_ = r - shift_mirror;
	if ( r_ < 0 ) {
		double v = -depth + a * cosh( expo * (r_ + shift_cosh) ); // potential well
		return v;
	}
	else {
		return -1.0 / ( 4.0 * ( r_ + shift_mirror ) );	// mirror charge
	}
}

//----------------------------------------------------------------------------------------------------------

void Gaussian_Pulse::initialize( Conf_Module* config, vector<Module*> dependencies )
{
	t_ofs = config->getParam("t_envl_center")->value / CONV_au_fs;
	fwhm = config->getParam("FWHM")->value / CONV_au_fs;
	ga = 2.0 * log( 2.0 ) / pow( this->fwhm, 2.0 );
	F0 = config->getParam("amplitude")->value / CONV_au_V * CONV_au_m ;
	omega0 = 2 * PI * 137.0 / ( config->getParam("wavelength")->value / CONV_au_nm ); // 2 pi c / lambda
	phi0 = config->getParam("ce_phase")->value;
	theta0 = config->getParam("chirp")->value * CONV_au_fs * CONV_au_fs;
}

inline double Gaussian_Pulse::electric_field( double t )
{
	//if (r <= 0) return 0;
	double t_ = t - t_ofs;	// + r / 137.0;
	double t2 = pow(t_, 2.0);
	double lorentz = 1.0; // + gamma * s2 / ( s2 + pow(r, 2.0) );
	double gauss = exp( -ga * t2 ) * cos( theta0 * t2 + omega0 * t_ + phi0 );
	return F0 * lorentz * gauss;
}

//----------------------------------------------------------------------------------------------------------

void Potential::register_dependencies( vector<Module*> dependencies )
{
	for ( size_t i = 0; i < dependencies.size(); i++ ) {
		if ( dependencies[i]->type.compare("pot_const") == 0  &&  dependencies[i]->serial == reg_serials[0] ) {
			LOG_INFO( "found a static potential to use: " + dependencies[i]->name);;
			well = dynamic_cast<Pot_const*>( dependencies[i] );
		}
		else if ( dependencies[i]->type.compare("pulse") == 0 ) {
			for( size_t id = 1; id < reg_serials.size(); id++ ) {
				if ( reg_serials[id] == dependencies[i]->serial ) {
					LOG_INFO( "found a laser pulse to use: " + dependencies[i]->name );;
					pulses.push_back( dynamic_cast<Laser_Field*>( dependencies[i] ) );
				}
			}
		}
	}
}

void Potential::initialize( Conf_Module* config, vector<Module*> dependencies )
{
	GET_LOGGER( "liee.Module.Potential" );
	r_range = config->getParam("r_range")->value / CONV_au_nm;
	double inner_cutoff = config->getParam("inner_cutoff")->value / CONV_au_eV;
	double dummy;
	well->get_outer_turningpoints( inner_cutoff, r_start, dummy );	// set r_start where the inner_cutoff-Energy is reached

	r0 = config->getParam("r0")->value / CONV_au_nm;
	k_geom = config->getParam("k_geom")->value;
	F_dc = -1 * config->getParam("U_dc")->value / CONV_au_V / k_geom / r0;
	t_charge = config->getParam("t_charge")->value / CONV_au_fs;
	gamma = config->getParam("near_amplf")->value;
	double s = config->getParam("near_width")->value / CONV_au_nm;
	s2 = pow( s , 2.0);
	wcap  = config->getParam("wcap")->value / CONV_au_nm;
	int_samples  = (int) config->getParam("int_samples")->value;

	// set spatial positions for the evaluation of the laser-pulse-field.
	Pulse_samples.resize( int_samples );
	double pos_range = r_start + r_range;
	dx_sample = pos_range / (int_samples - 1);
	for( size_t i = 0; i < int_samples; i++ ) {
		Pulse_samples[i] = Point( i * dx_sample, 0.0 );
	}

	// save the relevant module id's of constant potential and pulses and register them
	reg_serials.push_back( (int)config->getParam("my_pot_const")->value );
	BOOST_FOREACH( double serial, config->getParam("my_pulses")->values ) {
		reg_serials.push_back( (int) serial );
	}
	register_dependencies( dependencies );
}

void Potential::estimate_effort( Conf_Module* config, double & flops, double & ram, double & disk )
{
	double N = config->getParam("r_range")->value / config->getParam("dr")->value;
	flops += abs( 23 * N );
	ram += sizeof( Potential );
	disk += 0;
}

void Potential::summarize( map<string, string> & results )
{
	results["r_start"] = doub2str( r_start * CONV_au_nm ); //TODO boost convert
}

dcmplx Potential::V( double r, double t )
{
    double z = (r - (r_start + r_range - wcap) ) / wcap;	// CAP coordinate
	if (z > 0.0 && z <= 1.0) 								// r is in CAP territory
	{
		double r_ = r_start + r_range - wcap;
		return dcmplx( well->V( r_ ) +  V_Fdc( r_, t ) + V_pulse( r_, t ), V_cap( z ) );
	}
	else {
		return dcmplx( well->V( r ) +  V_Fdc( r, t ) + V_pulse( r, t ), 0.0 );
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
inline double Potential::V_cap( double z )
{
	return -6.13899 * pow(z, 4);
}

/*!
 * old version: V = int( dr' U / (k  r') = U / k * ln r', r' = r + r0
 * new version: constant F, linear V
 */
inline double Potential::V_Fdc( double r, double t )
{
	if ( r < 0 ) return 0; // inside the metal, the external electric field is shielded
	double V = F_dc * r;
	if ( t < t_charge ) { return t / t_charge * V; }
		else			{ return V; }
}

inline double Potential::F_pulse( double r, double t )
{
	double F = 0;
	double t_ = t + r / 137.0;
	BOOST_FOREACH( Laser_Field* p, pulses ) {
		F += p->electric_field( t_ );
	}
	return F;
}

inline double Potential::V_pulse( double r, double t )
{
	if ( t != t_current ) {
		// new time-step -> need to recalculate the integral samples using Simpson's rule
		Pulse_samples[0].y = 0;
		double fa = F_pulse( Pulse_samples[0].x, t );
		double fb, fm;
		double fac = (Pulse_samples[1].x - Pulse_samples[0].x) / 6.0;	// assume homogeneous grid

		for ( size_t i = 1; i < int_samples; i++ ) {
			fb = F_pulse( Pulse_samples[i].x, t );
			fm = F_pulse( 0.5*(Pulse_samples[i-1].x + Pulse_samples[i].x), t );
			Pulse_samples[i].y = Pulse_samples[i-1].y + fac * (fa + 4*fm + fb);
			fa = fb;
		}
		t_current = t;
	}

	if (r <= 0) {
		return 0;
	}
	size_t i = floor( r / dx_sample );
	if (i >= int_samples) {
		return Pulse_samples[int_samples-1].y;
	}
	// linear interpolate
	float weight = r / dx_sample - i;
	return (1 - weight) * Pulse_samples[i].y  +  weight * Pulse_samples[i+1].y;
}

//------------------------------------------------------------------------------------------------------------

inline double Spatial_Light_Modificator::electric_field( double t )
{
	return 0;
}

void Spatial_Light_Modificator::register_dependencies( vector<Module*> dependencies )
{
	for ( size_t i = 0; i < dependencies.size(); i++ ) {
		if ( dependencies[i]->type.compare("pulse") == 0 ) {
			for( size_t id = 1; id < pulse_serials.size(); id++ ) {
				if ( pulse_serials[id] == dependencies[i]->serial ) {
					LOG_INFO( "found a laser pulse to use: " + dependencies[i]->name );;
					pulses.push_back( dynamic_cast<Laser_Field*>( dependencies[i] ) );
				}
			}
		}
	}
}

void Spatial_Light_Modificator::initialize( Conf_Module* config, vector<Module*> dependencies )
{
	GET_LOGGER( "liee.Module.Pot_GaussSLM" );
	BOOST_FOREACH( double serial, config->getParam("orig_pulses")->values ) {
		pulse_serials.push_back( (int) serial );
	}
	register_dependencies( dependencies );

	//TODO improve the following sandbox code!
	/*
	double th = config->getParam("FWHM")->value / CONV_au_fs;
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
	//TODO
}

void Spatial_Light_Modificator::summarize( map<string, string> & results )
{
	//TODO
}

//------------------------------------------------------------------------------------------------------------

/*!
 * parabolic potential
 */
inline double Pot_Harm_Oscillator::V( double r )
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
	this->k = config->getParam("k")->value / 1e-30 * CONV_au_s * CONV_au_s;
	this->w = sqrt( k );
	this->shift = config->getParam("shift")->value  * 1e-9 / CONV_au_m;
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

/*!
 * the potential
 */
inline double Chulkov_Image_Potential::V( double r )
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
	else {		// if ( r > zim )
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
				for ( int k = 0; k < 3 * nlay; k++ ) {	//TODO(perf.) instead of stupid iteration from the beginning, start from the end
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
					for ( int k = 0; k < 3 * nlay; k++ ) {	//TODO(perf.) instead of stupid iteration from the beginning, start from the end
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
				} else { /*TODO(errorhandling)*/ }
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
	//this->r_end = config->getParam("r_end"]->value * 1e-9 / CONV_au_m;

	// Das hier und die Funktion im Anhang erzeugen das Potential für Beryllium(0001)
	// mit 41*2 Cos-Perioden im inneren und einer Wandstärke von 40*2 Perioden.
	as = 3.387;
	nlay = 41;
	double dslab = 2 * as * nlay  +  40 * as;
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
