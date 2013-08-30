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

	/*! for repeated evaluations on a fixed grid, performance is improved by caching the constant
	 *  portion	of the potential and access the combined potential via indexed.
	 *  set_grid() has to be called before using V_indexed()!
	 */
	dcmplx V_indexed( size_t ri, double t );

	/*! sets the grid transformation and stores the constant potential for repeated usage
	 *  important! has to be rerun in reinitialize() because the cache is not stored in the
	 *  checkpoint archive!
	 */
	void set_grid( double dr, size_t N );

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
	virtual void reinitialize( Conf_Module* config, vector<Module*> dependencies );
	virtual void estimate_effort( Conf_Module* config, double & flops, double & ram, double & disk );
	virtual void summarize( map<string, string> & results );

private:
	Pot_const* well;
	vector<Laser_Field*> pulses;
	vector<int> reg_serials;		///< remember which dependency modules to register
	vector<double> grid_r;			///< all r-positions of the grid (optional) (!not saved to archive!)
	vector<dcmplx> cache_Vconst;	///< indexed constant potential (V_well+V_dc+V_cap) (!not saved to archive!) (!not saved to archive!)
	vector<dcmplx> cache_Vwell;		///< indexed constant potential without V_dc, needed for the non-constant ramping-up period (!not saved to archive!)
	double grid_dr;					///< spacing of the requested grid

	double t_on;					///< activation start-time of the constant field, F_dc(t < t_on) = 0
	double t_full;					///< time when activation of the constant field is completed
	bool ini_erf;					///< if "true", the activation follows an erf()-like smooth swing, linear otherwise
	double deltaDC;					///< distance at which the static field has decreased to 1/e inside the tip (r<0)
	double deltaAC;					///< distance at which the Laser field has decreased to 1/e inside the tip (r<0)
	double r_range;					///< spatial range of simulation
	double r_start;					///< start position defined by where the negative-side repulsion potential reaches a certain threshold
	double wcap;					///< width of complex absorbing potential at the right side of the simulation range
	double r_cap;					///< start of complex absorbing potential, short for r_start + r_range - wcap
	double F_dc;					///< electric field at tip apex from applying a constant voltage
	double gamma;					///< resonant amplification factor
	double s2;						///< resonant amplification range
	//following members are not saved to checkpoint archive!
	double t_now;					///< as long as t==t_current the cache is valid
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
    	ar & t_on;
    	ar & t_full;
    	ar & ini_erf;
    	ar & deltaDC;
    	ar & deltaAC;
    	ar & r_range;
    	ar & r_start;
        ar & wcap;
        ar & r_cap;
        ar & F_dc;
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
    void serialize( Archive & ar, const unsigned int version ) //TODO
    {
    }
};

/*!
 * A potential composed of concatenated linear sections
 */
class Pot_Piecewise : public Pot_const
{
public:
	virtual inline double V( double r );
	virtual void get_outer_turningpoints( const double E, double & leftmost, double & rightmost );
	virtual double get_Vmin_pos();

	virtual void initialize( Conf_Module* config, vector<Module*> dependencies );
	virtual void reinitialize( Conf_Module* config, vector<Module*> dependencies );
	virtual void estimate_effort( Conf_Module* config, double & flops, double & ram, double & disk ) {} // TODO
	virtual void summarize( map<string, string> & results ) {}// TODO

	friend class boost::serialization::access;
    template<class Archive>
    void serialize( Archive & ar, const unsigned int version )
    {
    	ar & r_;
    	ar & V_;
    }

private:
    vector<double>	r_;		///< list of r-coordinates
    vector<double>	V_;		///< V_list[ r_list[i] {+epsilon} ]: list of potential energies at the positions in r_list
};
} //namespace liee

#endif /* POTENTIAL_H_ */
