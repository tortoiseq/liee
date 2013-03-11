/*!
 * optimizer.hpp
 *
 * Definition of generic structures and interfaces used by various optimisation algorithms.
 *
 *  Created on: 10-Jan-2012
 *      Author: quark
 */

#ifndef OPTIMIZER_H_
#define OPTIMIZER_H_

#include <limits>
#include <vector>

#include <boost/serialization/version.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/function.hpp>
#include <boost/foreach.hpp>

#include "my_util.hpp"
#include "module_config.hpp"


using namespace std;
namespace liee {
namespace opti {

/*!
 * Encapsulates the data structure for one evaluation request of the objective function.
 */
class Request {
public:
	/*! coordinate in the optimisation domain */	vector<double> 		x;
	/*! the new position for which an evaluation was
	    requested but no result has yet arrived */	vector<double> 		x_new;
	/*! value of objective function at x */			double 				y;
	/*! unique identifier of the request made */	int 				id;
	/*! indicates success: 0 or fail: otherwise */	int 				flag;

	friend class boost::serialization::access;
    /*! When the class Archive corresponds to an output archive, the
     *  & operator is defined similar to <<.  Likewise, when the class Archive
     *  is a type of input archive the & operator is defined similar to >>. */
    template<class Archive>
    void serialize( Archive & ar, const unsigned int version )
    {
        ar & id;
        ar & flag;
        ar & x;
        ar & x_new;
        ar & y;
    }
};

/*!
 * Common interface for asynchronous optimisation algorithms.
 */
class Asynch_Optimizer
{
public:
	/*! number of dimensions */					int 				dim;
	/*! lower limits indexed by dimension */	vector<double> 		lower_bounds;
	/*! upper limits indexed by dimension */	vector<double> 		upper_bounds;
	/*! number of _requested_ evaluations */	int					evaluations;
	/*! end condition if not converging */		int					max_eval;

	/*! next id (serial to match request and results) to send out  */
												int					next_id;

	/*! the so-far lowest value found */		double				global_min;
	/*! lowest value in current set */			double				current_min;
	/*! domain vector for global_min */			vector<double> 		global_min_pos;
	/*! domain vector for current_min */		vector<double> 		current_min_pos;
	/*!	*/										vector<double>		history;
	/*! */										int					history_sz;
	/*! */										int					history_i;

	/*! convergence tolerance */				double				tolerance;
	/*! convergence */							double				histo_var;

	/*! to identify the specific type */		string				type_name;

    friend class boost::serialization::access;
    /*! When the class Archive corresponds to an output archive, the
      * & operator is defined similar to <<.  Likewise, when the class Archive
      * is a type of input archive the & operator is defined similar to >>. */
    template<class Archive>
    void serialize( Archive & ar, const unsigned int version )
    {
        ar & dim;
        ar & lower_bounds;
        ar & upper_bounds;
        ar & tolerance;
        ar & histo_var;
        ar & evaluations;
        ar & max_eval;
        ar & next_id;
        ar & global_min;
        ar & current_min;
        ar & global_min_pos;
        ar & current_min_pos;
        ar & history;
        ar & history_sz;
        ar & history_i;
        ar & type_name;
    }

    Asynch_Optimizer()
    {
    	dim = 1;
    	tolerance = 1e-10;
    	histo_var = 6.66e300; //numeric_limits<double>::max();
    	evaluations = 0;
    	max_eval = (int) 2e9;
    	next_id = 1;
    	global_min = 6.66e300; //numeric_limits<double>::max();
    	current_min = 6.66e300; //numeric_limits<double>::max();
    	history_sz = 200;
    	history.resize( history_sz );
    	history_i = 0;
    }

    /*! Always have a virtual destructor in base classes. */
    virtual ~Asynch_Optimizer() {}

	/*! Initialises the optimisation problem for the first time.
	 *  Caller has to specify the bounds. The size of the bounds vectors determines the dimensionality.
	 *
	 *  Algorithm specific parameters are assumed to have been set by the constructor. See the algorithm of
	 *   choice for details. */
	virtual	int	initialise( const vector<double> & lower, const vector<double> & upper ) = 0;

	/*! Initialises the optimisation problem for the first time.
	 * 	Using richer data structure of Conf_Param which might include additional information
	 * 	about the desired probability distribution of variables, besides upper and lower bound.
	 * 	Not all algorithms will support those optional settings, so the default implementation
	 * 	here just extracts the bounds and calls the above initialize() method.
	 * 	Optimizers making use of the specifics of Conf_Param need to overwrite this method.
	 */
	virtual	int	initialise( const vector<Conf_Param*> & params )
	{
		vector<double> lo, hi;
		BOOST_FOREACH( Conf_Param* p, params ) {
			if ( not p->fixed ) {
				lo.push_back( p->min );
				hi.push_back( p->max );
			}
		}
		return initialise( lo, hi );
	}

	/*! In order for the optimisation to advance, new evaluations of the objective function are requested by this method.
	 *  The previous content of work is erased. If (work.size()==0 && result==0), then there is no work,
	 *  and the algorithm is waiting to assimilate more requested results.
	 *  @result	non-zero value indicates a non-recoverable error, therefore no more waiting is necessary. */
	virtual int generate_requests( vector<Request> & work ) = 0;

	/*! Feeds the results of previous requests from generate_work. If a result is faulty, this
	 *  should be marked by Request.flag!=0, so that the algorithm can recover through re-requesting. */
	virtual int assimilate_results( const vector<Request> & result ) = 0;
protected:
	virtual void convergence_test( const vector<Request> & population );
	virtual void convergence_test_( const vector<Request*> & population ) { return; }
};
//BOOST_CLASS_VERSION(Asynch_Optimizer, 0) TODO


/*!
 * Example implementation for Asynch_Optimiser.
 * One of the worst optimisation algorithms which is guaranteed to never converge.
 * It simply picks random domain vectors, equally distributed between lower and upper bounds.
 * Delivers the (so far) best lucky hit (lowest value of the objective function),
 * which slowly approaches the global minimum in the very long run.
 */
class Shot_Gun_Optimizer : public Asynch_Optimizer
{
public:
	/*! number of simultaneously requested evaluations */	int 				num_bullets_per_shot;
	/*! evaluations in progress */							vector<Request>		in_flight;
	/*! random number generator */							Ranq1 * 			random;

	/*! Default constructor sets number of parallel requests to 10. */
	Shot_Gun_Optimizer() : num_bullets_per_shot( 20 )
	{
		random = new Ranq1( 42 );
		in_flight.resize( num_bullets_per_shot );
	}

	friend class boost::serialization::access;
    /*! When the class Archive corresponds to an output archive, the
     *  & operator is defined similar to <<.  Likewise, when the class Archive
     *  is a type of input archive the & operator is defined similar to >>. */
    template<class Archive>
    void serialize( Archive & ar, const unsigned int version )
    {
    	// serialize base class information
    	ar & boost::serialization::base_object<Asynch_Optimizer>( *this );
        ar & num_bullets_per_shot;
        ar & random->v;
        ar & in_flight;
    }

	virtual int initialise( const vector<double> & lower, const vector<double> & upper );
	virtual int generate_requests( vector<Request> & work );
	virtual int assimilate_results( const vector<Request> & result );
};
//BOOST_CLASS_VERSION(Shot_Gun_Optimizer, 0) TODO


class Particle : public Request {
public:
	/*! velocity of particle */				vector<double> 		v;
	/*! best position of this particle */	vector<double> 		x_min;
	/*! lowest y on this trajectory */		double 				min_value;

	friend class boost::serialization::access;
    /*! When the class Archive corresponds to an output archive, the
     *  & operator is defined similar to <<.  Likewise, when the class Archive
     *  is a type of input archive the & operator is defined similar to >>. */
    template<class Archive>
    void serialize( Archive & ar, const unsigned int version )
    {
    	// serialize base class information
    	ar & boost::serialization::base_object<Request>( *this );
        ar & v;
        ar & x_new;
        ar & x_min;
        ar & min_value;
    }
};

/*!
 * Implementation of "A Modified Particle Swarm Optimizer" by Yuhui Shi et.al.
 *
 */
class Particle_Swarm_Optimizer : public Asynch_Optimizer
{
public:
	/*! population size */						int 		swarm_sz;
	/*! inertia weight (> 1: exploration; < 1: local search) */
												double 		inertia;
	/*! cognition weight (particle memory) */	double 		cognition;
	/*! social weight (follow the leader) */	double 		coherence;
	/*! index of momentarily best particle */	int			leader;
	/*! number of hops for which adjacent particles are considered
	 *  part of the neighbourhood of the centre
	 *  particle (in a ring topology) */		int			hood_sz;
	/*! random number generator */				Ranq1 * 	random;
//private:
	/*! particle population */			vector<Particle>	swarm;


public:
	/*! Default constructor sets population to 25 and other parameters to their generally recommended values */
	Particle_Swarm_Optimizer() :
		swarm_sz ( 25   ),
		inertia  ( 1.0 ),
		cognition( 2.0  ),
		coherence( 2.0  ),
		leader 	 ( 0 	),
		hood_sz  ( 2    )
	{
		random = new Ranq1( 4684 );
	}

	/*! Constructor with choice of population size. Change additional parameters (weights) manually! */
	Particle_Swarm_Optimizer( int population_size ) :
		swarm_sz ( population_size ),
		inertia  ( 1.0 ),
		cognition( 0.2  ),
		coherence( 0.3  ),
		leader 	 ( 0 	),
		hood_sz  ( 1    )
	{
		random = new Ranq1( 4684 );
	}

	friend class boost::serialization::access;
    /*! When the class Archive corresponds to an output archive, the
     *  & operator is defined similar to <<.  Likewise, when the class Archive
     *  is a type of input archive the & operator is defined similar to >>. */
    template<class Archive>
    void serialize( Archive & ar, const unsigned int version )
    {
    	// Serialise base class information
    	ar & boost::serialization::base_object<Asynch_Optimizer>( *this );
        ar & swarm_sz;
        ar & inertia;
        ar & cognition;
        ar & coherence;
        ar & leader;
        ar & hood_sz;
        ar & random->v;
        ar & swarm;
    }

	virtual int initialise( const vector<double> & lower, const vector<double> & upper );
	virtual int generate_requests( vector<Request> & work );
	virtual int assimilate_results( const vector<Request> & result );
	virtual void convergence_test( const vector<Particle> & population );
};
//BOOST_CLASS_VERSION(Shot_Gun_Optimizer, 0) TODO


/*!
 * The Rasterizer scans the search space in a regular fashion, its intended use is to
 * explore a region of interest and gather data for plotting. Like the Shot_Gun_Optimizer,
 * this Class can't be used effectively for optimisation, despite implementing the
 * Asynch_Optimizer interface.
 *
 * All work is requested in a single batch, make sure not to ask for more samples than you can handle!
 * Possible improvements reluctantly envisaged:
 * - smaller configurable batches
 * - random (or interlaced) ordering of requests to make early overview-plots over the whole range possible
 */
class Rasterizer : public Asynch_Optimizer
{
public:
	/*! resolution for each dimension */		vector<int> 	num_samples;
	/*! indicator for each dimensions scaling*/	vector<bool> 	logscale;
	/*! random number generator */				Ranq1 * 		random;
	/*! indicates that the one and only job creation was already done */
												bool			discharged;

	/*! boost serialization seems to need default constructor */
	Rasterizer()
	{
		this->random = new Ranq1( 42 );
		this->discharged = false;
	}

	/*! sets the number of samples uniformly for each dimension */
	Rasterizer( int num_samples_uni )
	{
		this->num_samples.resize( 1 );
		this->num_samples[0] = num_samples_uni;
		this->random = new Ranq1( 42 );
		this->discharged = false;
	}

	/*! sets the number of samples for every dimension */
	Rasterizer( vector<int> & num_samples )
	{
		this->num_samples = num_samples;
		this->random = new Ranq1( 42 );
		this->discharged = false;
	}

	/*! for convenience */
	Rasterizer( vector<double> & num_samples )
	{
		this->num_samples.resize( num_samples.size() );
		for ( size_t i = 0; i < num_samples.size(); i++ ) {
			this->num_samples[i] = (int) num_samples[i];
		}
		this->random = new Ranq1( 42 );
		this->discharged = false;
	}

	friend class boost::serialization::access;
    /*! When the class Archive corresponds to an output archive, the
     *  & operator is defined similar to <<.  Likewise, when the class Archive
     *  is a type of input archive the & operator is defined similar to >>. */
    template<class Archive>
    void serialize( Archive & ar, const unsigned int version )
    {
    	// serialize base class information
    	ar & boost::serialization::base_object<Asynch_Optimizer>( *this );
        ar & num_samples;
        ar & logscale;
        ar & random->v;
        ar & discharged;
    }

	virtual int initialise( const vector<double> & lower, const vector<double> & upper );
	virtual	int	initialise( const vector<Conf_Param*> & params );
	virtual int generate_requests( vector<Request> & work );
	virtual int assimilate_results( const vector<Request> & result ) { return 0; }
private:
	void iterate( int d, Request r, vector<Request> & work );

};

/*!
 * Basic minimisation method in one dimension.
 * Given three points a < b < c with f(b) < min( f(a), f(b) ) and a convergence tolerance,
 * the method will keep bisecting the interval with ratio of the golden section ( 3 - sqrt(5) ) / 2
 * by exclusion of the high point until the outer points are less then tolerance apart.
 */
class Golden_Section_Search
{
public:
	double tol;
	double fmin;
	Golden_Section_Search( double tolerance ) : tol( tolerance ) {}

	double minimize( boost::function<double (double)> func, const double a, const double b, const double c );
};

} // namespace opti
} // namespace liee


#endif /* OPTIMIZER_H_ */
