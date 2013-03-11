/*!
 * downhill_simplex.hpp
 *
 *  Created on: 24-Jan-2012
 *      Author: quark
 *
 * Implementation of the downhill simplex optimization method
 * by Nelder and Mead  algorithm along the lines of the implementation
 * from Numerical Recipes Webnote No. 15, Rev.1.
 *
 * TODO request all 3 possible relocations of the worst point at once to exploit parallel (partially wasted) execution for higher cadence. (then psum cache is of no use also)
 */

#ifndef DOWNHILLSIMPLEX_H_
#define DOWNHILLSIMPLEX_H_

#include <boost/serialization/version.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>

#include "optimizer.hpp"

using namespace std;
namespace liee {
namespace opti {

/*!
 * Multidimensional minimisation by the downhill simplex method of Nelder and Mead.
 * Also featuring Simulated Annealing as suggested by Numerical Recipes Webnote No. 15.
 */
class Downhill_Simplex : public Asynch_Optimizer {
public:
	/*! The number of points of the simplex. */			int 	mpts;

	/*! reflection coefficient */						double 	alfa;
	/*! expansion coefficient */						double 	gamma;
	/*! contraction coefficient */						double 	roh;
	/*! shrink coefficient */							double 	sigma;
	/*! Temperature to tune fluctuations for simulated annealing */
														double 	temperature;

	/*! Current simplex */						vector<Request> simplex;
	/*! to try a single reflection point */				Request trial;

	/*! random number generator */						Ranq1 * random;

	/*! Index of highest (worst) vertex. */				int 	ihi;
	/*! Index of next-highest vertex. */				int 	inhi;
	/*! Index of lowest (best) vertex. */				int 	ilo;


	/*! lowest function value incl. fluctuation) */		double 	ylo;
	/*! highest function value incl. fluctuation) */	double 	yhi;
	/*! 2nd-highest function value incl. fluctuation)*/	double 	ynhi;
	/*! to remember a temporary result between calls */	double 	ysave;

	/*! interim result: psum[d] = sum of x[d] over all vertices, needed for every reflection. */
												 vector<double>	psum; //TODO This saves roughly dim^2 additions per iteration, is it worth the trouble of storing?

	//! This flag is set after reflection-request to ensure possible swapping of high-value after the new point got evaluated.
														bool 	checkReflection;

	//! Indicates where to continue after return for the next request.
														int 	returnAddr;


	friend class boost::serialization::access;
    /*! When the class Archive corresponds to an output archive, the
     *  & operator is defined similar to <<.  Likewise, when the class Archive
     *  is a type of input archive the & operator is defined similar to >>. */
    template<class Archive>
    void serialize( Archive & ar, const unsigned int version )
    {
    	// serialize base class information
    	ar & boost::serialization::base_object<Asynch_Optimizer>( *this );
        ar & mpts;
        ar & alfa;
        ar & gamma;
        ar & roh;
        ar & sigma;
        ar & temperature;
        ar & simplex;
        ar & trial;
        ar & random->v;
        ar & ihi;
        ar & inhi;
        ar & ilo;
        ar & ylo;
        ar & yhi;
        ar & ynhi;
        ar & ysave;
        ar & psum;
        ar & checkReflection;
        ar & returnAddr;
    }

    //! Default constructor
	Downhill_Simplex();

	/*! Constructor
	*   @param maxEval	The maximum number of function evaluations before giving up
	*   @param convTol 	The fractional convergence tolerance to be achieved in the function value.
	*   @param temp		The temperature to start with; leads to limited upphill movement due to fluctuations.
	*/
	Downhill_Simplex( int maxEval, double convTol, double temp );

	//! Destructor
	~Downhill_Simplex();

	virtual int initialise( const vector<double> & lower, const vector<double> & upper );
	virtual int generate_requests( vector<Request> & work );
	virtual int assimilate_results( const vector<Request> & result );

	/*!
	 * Redefines the initial simplex (Do only use directly after initialise() !).
	 * Initialise() sets the initial simplex to the maximum extent of the search space.
	 * Use to probe a smaller region or provide randomised starting conditions.
	 */
	void set_simplex(const vector<double> & startingPoint, const double displacement );

	void set_simplex(const vector<double> & startingPoint, const vector<double> & displacements );

	bool is_waiting_for_results();

private:
	void request_reflected_vertex( double factor, vector<Request> & work );
	void recalculate_psum();
	void fluctuate_extrema();
	bool has_converged();
	void blubb();
};

} // namespace liee
} // namespace opti
#endif /* DOWNHILLSIMPLEX_H_ */
