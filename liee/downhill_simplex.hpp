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

#include "boost/serialization/version.hpp"
#include "boost/serialization/base_object.hpp"
#include "boost/serialization/vector.hpp"

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
	int     mpts;       ///< The number of points of the simplex.
	double  alfa;       ///< reflection coefficient
	double  gamma;      ///< expansion coefficient
	double  roh;        ///< contraction coefficient
	double  sigma;      ///< shrink coefficient
	double  temperature;///< Temperature to tune fluctuations for simulated annealing

	vector<Request> simplex; ///< Current simplex
	Request trial;      ///< to try a single reflection point
	Ranq1 * random;     ///< random number generator

	int     ihi;        ///< index of highest (worst) vertex
	int     inhi;       ///< index of next-highest vertex
	int     ilo;        ///< index of lowest (best) vertex

	double  ylo;        ///< lowest function value incl. fluctuation)
	double  yhi;        ///< highest function value incl. fluctuation)
	double  ynhi;       ///< 2nd-highest function value incl. fluctuation)
	double  ysave;      ///< to remember a temporary result between calls

	vector<double> psum; ///< interim result: psum[d] = sum of x[d] over all vertices, needed for every reflection.

	bool    checkReflection; ///< This flag is set after reflection-request to ensure possible swapping of high-value after the new point got evaluated.
	int     returnAddr; ///< Indicates where to continue after return for the next request.

	SERIALIZE( boost::serialization::base_object<Asynch_Optimizer>( *this )
	          & mpts & alfa & gamma & roh & sigma & temperature & simplex & trial & random->v
	          & ihi & inhi & ilo & ylo & yhi & ynhi & ysave & psum & checkReflection & returnAddr )

	//! Default constructor
	Downhill_Simplex();

	/*! Constructor
	*   @param maxEval  The maximum number of function evaluations before giving up
	*   @param convTol  The fractional convergence tolerance to be achieved in the function value.
	*   @param temp     The temperature to start with; leads to limited upphill movement due to fluctuations.
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
