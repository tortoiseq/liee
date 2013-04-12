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
 * Interface for the static portion of the potential.
 */
class Pot_const : public Module
{
public:
	//! potential energy at position r
	virtual inline double V( double r ) = 0;

	/*!
	 * exposes the classical turning-points for a particle with energy E
	 *
	 * limited default implementation provided:
	 * 		- uses root-finding by bracketing and bisection from pivot point get_Vmin_pos() outwards
	 * 		- is not guaranteed to reach the _outer_ turning points, yet more likely the inner ones (for a periodic potential)
	 * 		- much slower than evaluation of a moderate expression
	 * --> overwrite with analytic solution E=V(r), r minimal (maximal) for anything but testing
	 */
	virtual void get_outer_turningpoints( const double E, double & leftmost, double & rightmost );

	/*!
	 * returns the position r for which the potential has its global minimum
	 *
	 * limited default implementation provided:
	 * 		- uses minimum search by bracketing starting at r = 0
	 * 		- is not guaranteed to reach the lowest V(r) for weird periodic potentials,
	 * 		  but should do the trick in most practical cases
	 * 		- result is cached in case the method is called many times
	 */
	virtual double get_Vmin_pos();

protected:
	inline double deltaV( double r, double E ) { return E - V(r); }
	bool been_there;
	double r_min;
};

class Pot_Round_Well_wImage : public Pot_const
{
public:
	virtual inline double V( double r );
	virtual void initialize( Conf_Module* config, vector<Module*> dependencies );
	virtual void reinitialize( Conf_Module* config, vector<Module*> dependencies ) {
		GET_LOGGER( "liee.Module.Pot_Round_Well_wImage" );
	}
	virtual void estimate_effort( Conf_Module* config, double & flops, double & ram, double & disk ) {}	// effort depends on callers
	virtual void summarize( map<string, string> & results ) {}	// nothing to write home about

	virtual void get_outer_turningpoints( const double E, double & leftmost, double & rightmost );
	virtual double get_Vmin_pos();

protected:
	double width;
	double depth;
	double expo;
	double a;
	double shift_cosh;
	double shift_mirror;

	friend class boost::serialization::access;
    template<class Archive>
    void serialize( Archive & ar, const unsigned int version )
    {
        ar & width;
        ar & depth;
        ar & expo;
    	ar & a;
    	ar & shift_cosh;
    	ar & shift_mirror;
    }
};

/*!
 * For experimenting with the left potential wall's influence on the eigenvalue spectrum,
 * the strong exponential repulsion from Pot_Round_Well_wImage is substituded by a more
 * curvy polynomial.
 */
class Pot_Experimental : public Pot_const
{
public:
	virtual inline double V( double r );
	virtual void initialize( Conf_Module* config, vector<Module*> dependencies );
	virtual void reinitialize( Conf_Module* config, vector<Module*> dependencies ) {
		GET_LOGGER( "Pot_Experimental" );
	}
	virtual void estimate_effort( Conf_Module* config, double & flops, double & ram, double & disk ) {}	// effort depends on callers
	virtual void summarize( map<string, string> & results ) {}	// nothing to write home about

private:
	double width;
	double depth;
	double expo;
	double a;
	double shift_cosh;
	double shift_mirror;

	double v_scale;
	double h_scale;
	double power;
	vector<double> b;

	friend class boost::serialization::access;
    template<class Archive>
    void serialize( Archive & ar, const unsigned int version )
    {
        ar & width;
        ar & depth;
        ar & expo;
    	ar & a;
    	ar & shift_cosh;
    	ar & shift_mirror;
        ar & v_scale;
        ar & h_scale;
        ar & power;
        ar & b;
    }
};

//----------------------------------------------------------------------------------------------------------

/*!
 * Interface for the Laser Pulse.
 */
class Laser_Field : public Module
{
public:
	//! potential energy at position r and time t
	virtual inline double electric_field( double t ) = 0;
};


class Gaussian_Pulse : public Laser_Field
{
public:
    virtual inline double electric_field( double t );
	virtual void initialize( Conf_Module* config, vector<Module*> dependencies );
	virtual void reinitialize( Conf_Module* config, vector<Module*> dependencies ) {
		GET_LOGGER( "liee.Module.Gaussian_Pulse" );
	}
	virtual void estimate_effort( Conf_Module* config, double & flops, double & ram, double & disk ) {}	// effort depends on callers
	virtual void summarize( map<string, string> & results ) {}	// nothing to write home about

private:
	double t_ofs;
	double fwhm;
	double ga;
	double F0;
	double omega0;
	double phi0;
	double theta0;

protected:
	friend class boost::serialization::access;
    template<class Archive>
    void serialize( Archive & ar, const unsigned int version )
    {
        ar & t_ofs;
        ar & fwhm;
        ar & ga;
        ar & F0;
        ar & omega0;
        ar & phi0;
        ar & theta0;
    }

};


/*!
 * Additional pulse shaping by spatial light modificators (SLM) done to a Laser_Field dependency
 */
class Spatial_Light_Modificator : public Laser_Field
{
public:
	virtual inline double electric_field( double t );
	void register_dependencies( vector<Module*> dependencies );
	virtual void initialize( Conf_Module* config, vector<Module*> dependencies );
	virtual void reinitialize( Conf_Module* config, vector<Module*> dependencies ) {
		GET_LOGGER( "liee.Module.Pot_GaussSLM" );
		register_dependencies( dependencies );
	}
	virtual void estimate_effort( Conf_Module* config, double & flops, double & ram, double & disk );
	virtual void summarize( map<string, string> & results );

private:
	vector<Laser_Field*> pulses;
	vector<int> pulse_serials;		///< remember which dependency modules to register
	double L;
	double st;

protected:
	friend class boost::serialization::access;
    template<class Archive>
    void serialize( Archive & ar, const unsigned int version )
    {
    	ar & pulse_serials;
        ar & L;
        ar & st;
    }
};

//----------------------------------------------------------------------------------------------------------

/*!
 * General potential interface.
 */
class Potential : public Module
{
public:
	//! potential energy at position r and time t
	dcmplx V( double r, double t );
	//! real part of V(r,t)
	double Vr( double r, double t );

	//! grant access to the static potential
	Pot_const* getPot_const() { return well; }

	//! total range including CAP
	double get_r_range() { return r_range; }
	//! at which r to start the simulation
	double get_r_start() { return r_start; }
	//! starting-position of CAP is the physical end of the experiment
	double get_r_phys_end() { return r_range + r_start - wcap; }

	inline double deltaV( double r, double t, double E ) { return E - V(r, t).real(); }

	void register_dependencies( vector<Module*> dependencies );
	virtual void initialize( Conf_Module* config, vector<Module*> dependencies );
	virtual void reinitialize( Conf_Module* config, vector<Module*> dependencies ) {
		GET_LOGGER( "liee.Module.Pot_Gauss" );
		register_dependencies( dependencies );
	}
	virtual void estimate_effort( Conf_Module* config, double & flops, double & ram, double & disk );
	virtual void summarize( map<string, string> & results );

private:
	Pot_const* well;
	vector<Laser_Field*> pulses;
	vector<int> reg_serials;		///< remember which dependency modules to register


	double r_range;
	double r_start;
	double wcap;
	double F_dc;
	double r0;
	double k_geom;
	double t_charge;
	double gamma;
	double s2;
	double t_current;
	double dx_sample;
	size_t int_samples;
	vector<Point> Pulse_samples;

    inline double V_cap( double r );
	inline double V_Fdc( double r, double t );
	inline double F_pulse( double r, double t );
	inline double V_pulse( double r, double t );
	inline double integral_gauss_x_oscill( double t_a, double t_b );

protected:
	friend class boost::serialization::access;
    template<class Archive>
    void serialize( Archive & ar, const unsigned int version )
    {
    	ar & reg_serials;
    	ar & r_range;
    	ar & r_start;
        ar & wcap;
        ar & F_dc;
        ar & r0;
        ar & k_geom;
        ar & t_charge;
        ar & gamma;
        ar & s2;
    	ar & int_samples;
    }
};


/*!
 * Parabolic potential of the harmonic oscillator to compare results with the known analytical solutions.
 * (mass=1)
 */
class Pot_Harm_Oscillator : public Pot_const
{
public:
	virtual inline double V( double r );
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
	double k;		///< force constant. equals omega^2, since m=1
	double w;		///< eigen frequency in (atomic time unit)^-2
	double shift;	///< shift to the right (to be able to stay in the positive domain for moderate quantum numbers)

	friend class boost::serialization::access;
    template<class Archive>
    void serialize( Archive & ar, const unsigned int version )
    {
    	ar & k;
    	ar & w;
    	ar & shift;
    }
};

class Chulkov_Image_Potential : public Pot_const
{
public:
	virtual inline double V( double r );
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

	double shift;	///< shift to the right (to be able to stay in the positive domain for moderate quantum numbers)

	friend class boost::serialization::access;
    template<class Archive>
    void serialize( Archive & ar, const unsigned int version )
    {
    }
};

/*!
 * For testing purpose. V(r):=0 for all r.
 */
class Zero_Pot : public Pot_const
{
public:
	virtual inline double V( double r ) { return 0; }
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
} //namespace liee

#endif /* POTENTIAL_H_ */
