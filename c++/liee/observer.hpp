/*
 * observer.hpp
 *
 *  Created on: 30-Jan-2012
 *      Author: quark
 */

#ifndef OBSERVER_HPP_
#define OBSERVER_HPP_

#include <boost/serialization/version.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>

#include "my_util.hpp"
#include "module_config.hpp"

using namespace std;
namespace liee {

/*!
 * Base class for Observers of iterative Experiments (Solver).
 * The interface function observe() will be called every iteration and the
 * specialised implementation accesses the state to evaluate the properties
 * interested in. The data members can optionally be used to keep track of
 * sampling rate.
 */
class Observer : public Module
{
public:
	int 	N;			///< total number of observations
	double 	t_range;	///< simulation end-time
	string 	filename;	///< write the result to file
	double 	t_last;		///< time of last measurement
	int		counter;
	double 	dt;			///< time resolution of experiment

	friend class boost::serialization::access;
    template<class Archive>
    void serialize( Archive & ar, const unsigned int version )
    {
        ar & N;
        ar & t_range;
        ar & counter;
        ar & filename;
        ar & t_last;
        ar & dt;
    }

	//! state should point to a Solver
	virtual void observe( Module* state ) = 0;
};

/*!
 * This Observer can be used together with a Solver of time-dependent Schroedinger equation.
 * It records the whole wave function (WV) N times during the course of the simulation and writes it
 * to a text file. The file contains the WV-samples as rows, with spatial samples in tab-separated columns.
 * If complex values are stored, real and imaginary part are separated by a comma. //TODO better interlaced: row-a real, row-b imag, ...
 * To reduce the size of the output the user can decrease the spatial resolution by the factor
 * downsample (//TODO downsample is obsolete)or ignore the complex nature of the WV by setting do_square to true.
 */
class Obs_Snapshot_WF : public Observer
{
public:
	size_t	num_r;		///< number of spatial samples to produce
	size_t	num_t;		///< number of temporal samples to produce
	size_t 	step_r;		///< number of spatial samples to skip over between those being saved
	size_t 	step_t;		///< number of temporal samples to skip over between the saved ones
	string	format;		///< cache format string for fprintf
	bool 	do_square;	///< only save the square of the complex wf
	bool	do_normalize;///< if true, normalise integral to Psi Psi* == 1
	size_t	ir0;		///< index of start-pos of recording window
	size_t	ir1;		///< index of end-pos of recording window
	double	t0;			///< start-time of recording window
	double	t1;			///< end-time of recording window
	size_t  writtenLns;	///< number of lines written to the outfile
	size_t	foundLns;	///< number of lines found in outfile (can be bigger than writtenLns after load from checkpoint)

    virtual void observe( Module* state );
	virtual void initialize( Conf_Module* config, vector<Module*> dependencies );
	virtual void reinitialize( Conf_Module* config, vector<Module*> dependencies );
	virtual void estimate_effort( Conf_Module* config, double & flops, double & ram, double & disk );
	virtual void summarize( map<string, string> & results ){}

	friend class boost::serialization::access;
    /*! When the class Archive corresponds to an output archive, the
     *  & operator is defined similar to <<.  Likewise, when the class Archive
     *  is a type of input archive the & operator is defined similar to >>. */
    template<class Archive>
    void serialize( Archive & ar, const unsigned int version )
    {
    	ar & boost::serialization::base_object<Observer>( *this );
    	ar & num_r;
    	ar & num_t;
        ar & step_r;
        ar & step_t;
        ar & format;
        ar & do_square;
        ar & do_normalize;
        ar & ir0;
        ar & ir1;
        ar & t0;
        ar & t1;
        ar & writtenLns;
    }
};


/*!
 * This Observer can be used together with a Solver of time-dependent Schroedinger equation.
 * It integrates the Psi^2 for the first and last time step. The result is a single valued:		//TODO describe extended functionality
 * the relative portion of tunnelled (electron) density.
 */
class Obs_Tunnel_Ratio : public Observer
{
public:
	vector<double>	psi_sqr;	///< total probability over integration region for each time sample
	int 	step_t;
	int		t_samples;
	double 	tunnel_ratio;		///< portion of probability lost during simulation
	double	ra;
	double	rb;
	bool 	is_objective;

    virtual void observe( Module* state );
	virtual void initialize( Conf_Module* config, vector<Module*> dependencies );
	virtual void reinitialize( Conf_Module* config, vector<Module*> dependencies );
	virtual void estimate_effort( Conf_Module* config, double & flops, double & ram, double & disk );
	virtual void summarize( map<string, string> & results );

	friend class boost::serialization::access;
    template<class Archive>
    void serialize( Archive & ar, const unsigned int version )
    {
    	ar & boost::serialization::base_object<Observer>( *this );
        ar & psi_sqr;
        ar & step_t;
        ar & t_samples;
        ar & tunnel_ratio;
        ar & ra;
        ar & rb;
        ar & is_objective;
    }
};

/*!
 * The Obs_JWKB_Tunnel observer can be used in conjunction with a solver of the time-dependent Schroedinger
 * equation, but the WF-evolution is ignored. Instead the time-dependent potential is used for the JWKB
 * approximation of the quasi-stationary tunnel rate by integrating over the barrier. The resulting tunnel
 * rate in comparison with the dynamic rate (Obs_Tunnel_Ratio) gives insights on the ratio of
 * absorption+over-the-barrier versus through-the-barrier emission.
 */
class Obs_JWKB_Tunnel : public Observer
{
public:
	vector<double>	j;	///< total probability over integration region for each time sample
	vector<Point>	rt;	///< debug: record turning points
	vector<double>	A;	///< debug: record integral of square-root(V-E) over barrier
	double 	sum_j;		///< portion of probability lost during simulation
	Potential* V;		///< potential for integration
	double E;			///< energy of initial state, reverse engineered from curvature of WF at the potentials minimum
	double dr;			///< global dx, used as small step into the barrier to reverse engineer the effective Field strength at the surface
	double Vp0;			///< cache V(0+dx, 0) potential without field one step into the barrier (assuming field is deactivated at t=0)
	double last_r2;		///< cache latest second turning point
	int N;				///< number of samples for numerical integration over the barrier
	double g;			///< JWKB constant
	double r_end;		///< for r > r_end, the potential stops to be useful (because CAP takes over)
	bool burst;			///< true if: at least at one sampling time the barrier was pushed bellow the states initial energy
	int 	step_t;
	int		t_samples;
	bool 	is_objective;

    virtual void observe( Module* state );
	virtual void initialize( Conf_Module* config, vector<Module*> dependencies );
	virtual void reinitialize( Conf_Module* config, vector<Module*> dependencies );
	virtual void estimate_effort( Conf_Module* config, double & flops, double & ram, double & disk );
	virtual void summarize( map<string, string> & results );

	friend class boost::serialization::access;
    template<class Archive>
    void serialize( Archive & ar, const unsigned int version )
    {
    	ar & boost::serialization::base_object<Observer>( *this );
    	ar & j;
    	ar & rt;
    	ar & A;
    	ar & sum_j;
    	ar & E;
    	ar & dr;
    	ar & Vp0;
    	ar & last_r2;
    	ar & N;
    	ar & g;
    	ar & r_end;
    	ar & burst;
        ar & step_t;
        ar & t_samples;
        ar & is_objective;
    }
};

} // namespace liee

#endif /* OBSERVER_HPP_ */
