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
	int 	step_r;
	int 	step_t;
	string	format;
	bool 	do_square;	///< only save the square of the complex wf
	bool	do_normalize;

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
        ar & step_r;
        ar & step_t;
        ar & do_square;
        ar & do_normalize;
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
	//bool	written;
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
        //ar & written;
        ar & tunnel_ratio;
        ar & ra;
        ar & rb;
        ar & is_objective;
    }
};

} // namespace liee

#endif /* OBSERVER_HPP_ */
