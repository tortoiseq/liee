#ifndef POTENTIAL_H_
#define POTENTIAL_H_

#include <boost/serialization/version.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>

#include <limits>

#include "my_util.hpp"
#include "module_config.hpp"

using namespace std;
namespace liee
{

double const PI = 4.0 * atan( 1.0 );
double const sqr_pi_quater = PI * PI / 4.0;
double const sqrt_2_half = sqrt( 2.0 ) / 2.0;

/*!
 * Definition of potential interface.
 */
class Potential : public Module
{
public:
	//! potential energy at position r and time t
	virtual dcmplx V( double r, double t ) = 0;
	virtual double Vr( double r, double t ) = 0;

	//! exposes the classical turning-points for a particle with energy E
	virtual void get_outer_turningpoints( const double E, double & leftmost, double & rightmost ) = 0;

	//! returns the position r for which the potential has its global minimum
	virtual double get_Vmin_pos() = 0;

	/*! total range including CAP */
	double r_range;
	double r_start;
};

/*!
 * For testing purpose. V(r):=0 for all r.
 */
class Zero_Pot : public Potential
{
public:
	virtual dcmplx V( double r, double t ) { return dcmplx(0,0); }
	virtual double Vr( double r, double t ) { return 0; }
	virtual void get_outer_turningpoints( const double E, double & leftmost, double & rightmost ) {
		leftmost = rightmost = std::numeric_limits<double>::quiet_NaN();
	}
	virtual double get_Vmin_pos() { return 0; }
	// empty interface declarations
	virtual void initialize( Conf_Module* config, vector<Module*> dependencies ) {
		GET_LOGGER( "liee.Module.Zero_Pot" );
	}
	virtual void reinitialize( Conf_Module* config, vector<Module*> dependencies ) {
		GET_LOGGER( "liee.Module.Zero_Pot" );
	}
	virtual void estimate_effort( Conf_Module* config, double & flops, double & ram, double & disk ) {} // whatever
	virtual void summarize( map<string, string> & results ) {}

	friend class boost::serialization::access;
    template<class Archive>
    void serialize( Archive & ar, const unsigned int version ) {}
};

/*!
 * Potential well with memory charge, adiabatic activated dc-field, absorbing potential on the right side and chirped
 * Gaussian Laser pulse, experiencing resonant amplification in the near-field.
 */
class Pot_Gauss : public Potential
{
public:
	Pot_Gauss() {}
	virtual dcmplx V( double r, double t );
	virtual double Vr( double r, double t );
	virtual void initialize( Conf_Module* config, vector<Module*> dependencies );
	virtual void reinitialize( Conf_Module* config, vector<Module*> dependencies ) {
		GET_LOGGER( "liee.Module.Pot_Gauss" );
	}
	virtual void estimate_effort( Conf_Module* config, double & flops, double & ram, double & disk );
	virtual void summarize( map<string, string> & results );
	virtual void get_outer_turningpoints( const double E, double & leftmost, double & rightmost );
	virtual double get_Vmin_pos();

	/*!
	 * write a matrixes to "fileprefix"-F.dat and "fileprefix"-V.dat containing sampled electric field strength
	 * and potential for the given time-spatial range and resolution.
	 */
	void debug_output( string filename, double ra, double rb, int Nr, double ta, double tb, int Nt );

//private:
	double width;
	double depth;
	double expo;
	double inner_cutoff;
	double wcap;
	double F_dc;
	double r0;
	double k_geom;
	double t_charge;
	double t_ofs;
	double fwhm;
	double ga;
	double F0;
	double omega0;
	double phi0;
	double theta0;
	double gamma;
	double s2;
	double a;
	double shift_cosh;
	double shift_mirror;
	double t_current;
	double dx_sample;
	size_t int_samples;
	vector<Point> Pulse_samples;

	friend class boost::serialization::access;
    template<class Archive>
    void serialize( Archive & ar, const unsigned int version )
    {
        ar & r_range;
    	ar & r_start;
        ar & width;
        ar & depth;
        ar & expo;
        ar & inner_cutoff;
        ar & wcap;
        ar & F_dc;
        ar & r0;
        ar & k_geom;
        ar & t_charge;
        ar & t_ofs;
        ar & fwhm;
        ar & ga;
        ar & F0;
        ar & omega0;
        ar & phi0;
        ar & theta0;
        ar & gamma;
        ar & s2;
    	ar & a;
    	ar & shift_cosh;
    	ar & shift_mirror;
    	ar & int_samples;
    }

	virtual void ini_common( Conf_Module* config, vector<Module*> dependencies );
    inline double V_cap( double r );
	inline double V_stat( double r );
	inline double V_Fdc( double r, double t );
	virtual double F_pulse( double r, double t );
	inline double V_pulse( double r, double t );
	inline double integral_gauss_x_oscill( double t_a, double t_b );
};

/*!
 * Gaussian Potential (see Pot_Gauss) with additional pulse shaping by spatial light modificators (SLM)
 */
class Pot_GaussSLM : public Pot_Gauss
{
public:
	Pot_GaussSLM() {}
	virtual void initialize( Conf_Module* config, vector<Module*> dependencies );
	virtual void reinitialize( Conf_Module* config, vector<Module*> dependencies ) {
		GET_LOGGER( "liee.Module.Pot_GaussSLM" );
	}
	virtual void estimate_effort( Conf_Module* config, double & flops, double & ram, double & disk );
	virtual void summarize( map<string, string> & results );

	double L;
	double st;

	friend class boost::serialization::access;
    template<class Archive>
    void serialize( Archive & ar, const unsigned int version )
    {
    	ar & boost::serialization::base_object<Pot_Gauss>( *this );
        ar & L;
        ar & st;
    }

private:
	virtual double F_pulse( double r, double t );
};

/*!
 * Parabolic potential of the harmonic oscillator to compare results with the known analytical solutions.
 * (mass=1)
 */
class Pot_Harm_Oscillator : public Potential
{
public:
	Pot_Harm_Oscillator() {}
	virtual dcmplx V( double r, double t );
	virtual double Vr( double r, double t );
	virtual void initialize( Conf_Module* config, vector<Module*> dependencies );
	virtual void get_outer_turningpoints( const double E, double & leftmost, double & rightmost );
	virtual double get_Vmin_pos();

	double analytic_eigenfunction( int n, double x );

	// empty interface declarations
	virtual void reinitialize( Conf_Module* config, vector<Module*> dependencies ) {
		GET_LOGGER( "liee.Module.Pot_Harm_Oscillator" );
	}
	virtual void estimate_effort( Conf_Module* config, double & flops, double & ram, double & disk ) {} // whatever
	virtual void summarize( map<string, string> & results ) {} // no relevant information to summarise

private:
	/*! force constant. equals omega^2, since m=1 */	double k;
	/*! eigen frequency in (atomic time unit)^-2 */		double w;
	/*! shift parabola by this amount to the right (to be able to stay in the positive domain for moderate quantum numbers)*/
														double shift;

	friend class boost::serialization::access;
    template<class Archive>
    void serialize( Archive & ar, const unsigned int version )
    {
        ar & r_range;
        ar & r_start;
    	ar & k;
    	ar & w;
    	ar & shift;
    }
};

class Chulkov_Image_Potential : public Potential
{
public:
	Chulkov_Image_Potential() {}
	virtual dcmplx V( double r, double t );
	virtual double Vr( double r, double t );
	virtual void initialize( Conf_Module* config, vector<Module*> dependencies );
	virtual void get_outer_turningpoints( const double E, double & leftmost, double & rightmost );
	virtual double get_Vmin_pos();

	// empty interface declarations
	virtual void reinitialize( Conf_Module* config, vector<Module*> dependencies ) {
		GET_LOGGER( "liee.Module.Chulkov_Image_Potential" );
	}
	virtual void estimate_effort( Conf_Module* config, double & flops, double & ram, double & disk ) {} // whatever
	virtual void summarize( map<string, string> & results ) {} // no relevant information to summarise

private:
	int nlay;
	double as;
	double A10;
	double A1;
	double A2;
	double beta;

	double A20;
	double D;
	double z1;
	double A3;
	double alpha;
	double lambda;
	double zim;

	/*! shift parabola by this amount to the right (to be able to stay in the positive domain for moderate quantum numbers)*/
	double shift;

	friend class boost::serialization::access;
    template<class Archive>
    void serialize( Archive & ar, const unsigned int version )
    {
    }
};

} //namespace liee

#endif /* POTENTIAL_H_ */
