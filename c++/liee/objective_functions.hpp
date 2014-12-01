/*!
 * objective_functions.hpp
 *
 *  Created on: 6-Jan-2012
 *      Author: quark
 *
 * A collection of standard test problems for optimisation problems.
 */

#ifndef TESTOBJECTIVE_H_
#define TESTOBJECTIVE_H_

using namespace std;
namespace liee {

/*!
* Just the square of x. (for readability)
*/
inline double sqr( double x ) { return x * x; }

/*!
* The generalised Rosenbrock's Function is a test for multivariant optimisation.
* The minimum is y_min(1, ..., 1) = 0.
* The coordinates should be bound by the interval [-30;30].
* @param x      point in search space as vector of double
* @return       the value of the objective function at x
*/
inline double rosenbrock( const vector<double> & x )
{
	double result = 0;
	for ( size_t i = 0; i < x.size()-1; i++ ) {
		result += 100 * sqr( x[i+1] - sqr( x[i] ) ) + sqr( x[i] - 1 );
	}
	return result;
}

/*!
* The generalised Rastrigin Function is a test for multivariant optimisation.
* This function is a fairly difficult problem due to its large search space and its large number of local minima.
* The minimum is y_min(1, ..., 1) = 0.
* The coordinates should be bound by the interval [-5.12;5.12].
* @param x      point in search space as vector of double
* @return       the value of the objective function at x
*/
inline double rastrigin( const vector<double> & x )
{
	double result = 10.0 * x.size();
	for ( size_t i = 0; i < x.size()-1; i++ ) {
		result += sqr( x[i] ) - 10.0 * cos( 2 * PI * x[i] );
	}
	return result;
}

/*!
* Sum of all coordinates-squared and shifted to (1, ..., 1)
* The minimum is y_min(1, ..., 1) = 0.
* @param x      point in search space as vector of double
* @return       the value of the objective function at x
*/
inline double paraboloid( const vector<double> & x )
{
	double result = 0;
	for ( size_t i = 0; i < x.size(); i++ ) {
		result += sqr( 1 - x[i] );
	}
	return result;
}

inline double equation_solving( const vector<double> & x )
{
	double err = 0.0;

	double A  = x[0];
	double B  = x[1];
	double C  = x[2];
	double D  = x[3];
	double W  = x[4];
	double x1 = x[5];
	double x2 = x[6];

	err += pow( (1/C/x1 - D) - (B*cosh(A*(x1-W/2))), 2.0);
	err += pow( (B*cosh(A*(x2-W/2))) - (-1/4/(x2-W)), 2.0);
	err += pow( (-1/(C*x1*x1)) - (A*B*sinh(A*(x1-W/2))), 2.0);
	err += pow( (A*B*sinh(A*(x2-W/2))) - (1/(4*pow(W-x2, 2.0))), 2.0);

	return err;
}

} //namespace liee

#endif /* TESTOBJECTIVE_H_ */
