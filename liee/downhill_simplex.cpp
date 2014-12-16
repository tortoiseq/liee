/*! downhill_simplex.cpp
 *
 * Implementation of the downhill simplex optimization method
 * by Nelder and Mead  algorithm along the lines of the implementation
 * from Numerical Recipes Webnote No. 15, Rev.1.
 */

#include "downhill_simplex.hpp"

using namespace std;
using namespace liee;

namespace liee {
namespace opti {

string vec2str( vector<double> v ) {
	stringstream ss;
	ss.precision( 3 );
	ss << "(";
	for ( size_t i = 0; i < v.size(); i++ ) {
		if ( i > 0 ) { ss << "|"; }
		ss << v[i];
	}
	ss << ")";
	return ss.str();
}

Downhill_Simplex::Downhill_Simplex()
{
	random = new Ranq1( 654646ULL );

	max_eval = 1e6;
	tolerance = 1e-6;
	temperature = 0;

	alfa = 1;
	gamma = 2;
	roh = 0.5;
	sigma = 0.5;

	checkReflection = false;
}

Downhill_Simplex::Downhill_Simplex(int maxEval, double convTol, double temp)
{
	random = new Ranq1( 654646ULL );

	max_eval = maxEval;
	tolerance = convTol;
	temperature = temp;

	alfa = 1;
	gamma = 2;
	roh = 0.5;
	sigma = 0.5;

	checkReflection = false;
}

Downhill_Simplex::~Downhill_Simplex() {}

int Downhill_Simplex::initialise( const vector<double> & lower, const vector<double> & upper )
{
	if ( lower.size() != upper.size() ) return 1;
	dim = lower.size();
	lower_bounds = lower;
	upper_bounds = upper;
	mpts = dim + 1;
	psum.resize( mpts );

	simplex.resize( mpts );
	for ( int v = 0; v < mpts; v++ ) {
		simplex[v].flag = -1;
		simplex[v].id = -1;
		simplex[v].x.resize( dim );
		simplex[v].x_new.resize( dim );
		for ( int d=0; d < dim; d++ ) {
			if ( d == v ) simplex[v].x[d] = upper_bounds[d];
			else          simplex[v].x[d] = lower_bounds[d];
		}
	}
	trial.flag = -1;
	trial.id = -1;
	trial.x.resize( dim );
	trial.x_new.resize( dim );
	return 0;
}

void Downhill_Simplex::set_simplex(const vector<double> & startingPoint, const double displacement)
{
	//TODO implement
}

void Downhill_Simplex::set_simplex(const vector<double> & startingPoint, const vector<double> & displacements)
{
	//TODO implement
}

int Downhill_Simplex::generate_requests( vector<Request> & work )
{
	//cout << "starting request-generation with returnAddr = " << returnAddr << "\n";
	work.clear();

	if ( this->is_waiting_for_results() ) return 0;

	if (evaluations == 0) {
		// During the first run, the function needs to be evaluated at all vertices of the simplex
		for ( int v = 0; v < mpts; v++ ) {
			simplex[v].id = next_id++;
			simplex[v].x_new = simplex[v].x;
			work.push_back( simplex[v] );
			evaluations++;
		}
		returnAddr = 1;
		return 0;
	}

	// We added a thermal fluctuation to all the current vertices, but we subtract it here,
	// so as to give the simplex a thermal Brownian motion: It "likes" to accept any suggested change.
	double yflu = trial.y + temperature * log( random->doub() );	// (mind the negative log)
	if ( checkReflection && yflu < yhi ) {
		// If reflected trial point is better than the current worst, than replace the highest.
		//cout << "trial accepted \n";
		simplex[ihi].y = trial.y;
		for ( int d = 0; d < dim; d++ ) {
			psum[d] += trial.x[d] - simplex[ihi].x[d];  //alter cached sum accordingly {ptry[d] - p[ihi*ndim+d];}
			simplex[ihi].x[d] = trial.x[d];
		}
		fluctuate_extrema();
	}
	trial.y = yflu;
	checkReflection = false;

	if ( returnAddr == 1 )
	{
		newIteration:
		recalculate_psum();
		fluctuate_extrema();
		if ( this->has_converged() ) {
			cout << "has converged\n";
			blubb();
			return 0;
		}
		if ( evaluations >= max_eval ) {
			cout << "has exceeded\n";
			return -2;
		}
		// First extrapolate by a factor -1 through the face of the
		// simplex across from the high point, i.e., reflect the simplex from the high point.
		//cout << "\tREFLECT\n";
		request_reflected_vertex( -alfa, work );
		returnAddr = 2;
		return 0;
	}
	else if ( returnAddr == 2 )
	{
		if ( trial.y <= ylo ) {
			// Gives a result better than the best point, so try an additional expansion by a factor 2.
			//cout << "\tEXPAND\n";
			request_reflected_vertex( gamma, work );
			returnAddr = 1;
			return 0;
		}
		else {
			if ( trial.y >= ynhi ) {
				// The reflected point is worse than the second-highest, so look for an
				// intermediate lower point, i.e., do a one-dimensional contraction.
				//cout << "\tPULL BACK\n";
				ysave = yhi;
				request_reflected_vertex( roh, work );
				returnAddr = 3;
				return 0;
			}
		}
		// Reflection accepted
		goto newIteration;
	}
	else if ( returnAddr == 3 )
	{
		if ( trial.y >= ysave ) {
			// Can't seem to get rid of that high point. Better contract around the lowest (best) point.
			//cout << "\tCONTRACT\n";
			for ( int v = 0; v < mpts; v++) {
				if (v != ilo)
				{
					for ( int d = 0; d < dim; d++ ) {
						simplex[v].id = next_id++;
						simplex[v].x_new[d] = sigma * ( simplex[v].x[d] + simplex[ilo].x[d] );
					}
					work.push_back( simplex[v] );
					evaluations++;
				}
			}
			returnAddr = 1;
			return 0;
		}
		// Pull back accepted
		returnAddr = 1;
		goto newIteration;
	}
	return -1; // should never reach this line
}

int Downhill_Simplex::assimilate_results( const vector<Request> & result )
{
	for ( size_t r = 0; r < result.size(); r++ )
	{
		if ( trial.id == result[r].id ) {
			trial.id = -1;  // free slot!
			if ( result[r].flag == 0 )
			{
				trial.x = trial.x_new; //TODO avoid code duplication by integrating 'trial' into simplex data structure
				trial.y = result[r].y;
				trial.flag = 0;
				//cout << "trial returned " << vec2str( trial.x ) << "->" << trial.y << "\n";
			}
			else {
				cout << "DESASTER \n";	//TODO error handling: something for corrupt results
			}
			break;
		}

		for ( int v = 0; v < mpts; v++ )	{
			if ( simplex[v].id == result[r].id )
			{
				simplex[v].id = -1; // free slot!
				if ( result[r].flag == 0 )
				{
					simplex[v].x = simplex[v].x_new;
					simplex[v].y = result[r].y;
					simplex[v].flag = 0;
				}
				else {
					cout << "DESASTER \n"; //TODO error handling:  something for corrupt results
				}
				break;
			}
		}
	}
	return 0;
}

bool Downhill_Simplex::has_converged()
{
	convergence_test( simplex );
	//histo_var = 2.0 * abs( simplex[ihi].y - simplex[ilo].y ) / ( simplex[ihi].y + simplex[ilo].y + 1e-10 );
	if ( histo_var < tolerance )
		return true;
	else
		return false;
}

void Downhill_Simplex::request_reflected_vertex( double fac, vector<Request> & work )
{
	double fac1 = (1.0 - fac) / dim;
	double fac2 = fac1 - fac;
	for ( int d = 0; d < dim; d++ ) {
		trial.x_new[d] = psum[d] * fac1 - simplex[ihi].x[d] * fac2;
	}
	trial.id = next_id++;
	work.push_back( trial );
	evaluations++;
	checkReflection = true;
	//cout << "trying to reflect )" << vec2str( trial.x_new ) << "\n";
}

void Downhill_Simplex::recalculate_psum()
{
	for ( int d = 0; d < dim; d++ )
	{
		double sum = 0.0;
		for ( int v = 0; v < mpts; v++ ) {
			sum += simplex[v].x[d];
		}
		psum[d] = sum;
	}
}

void Downhill_Simplex::fluctuate_extrema()
{
	// the thermal fluctuations only appear in the ylo/yhi/ynhi/ytry variables; simplex[v].y remains unchanged
	ilo = 0;
	ylo = simplex[ilo].y - temperature * log( random->doub() );  // adding fluctuation (mind the negative log)
	ihi = 1;
	yhi = simplex[ihi].y - temperature * log( random->doub() );
	if ( ylo > yhi ) {
		ihi = 0;
		ilo = 1;
		ynhi = yhi;
		yhi = ylo;
		ylo = ynhi;
	}

	for ( int v = 2; v < mpts; v++ )
	{
		double yflu = simplex[v].y - temperature * log( random->doub() );
		if (yflu <= ylo)	{
			ilo = v;
			ylo = yflu;
		}
		if (yflu > yhi) {
			ynhi = yhi;
			ihi = v;
			yhi = yflu;
		}
		else if (yflu > ynhi && v != ihi) {
			ynhi = yflu;
		}
	}
	//blubb();
}

bool Downhill_Simplex::is_waiting_for_results()
{
	//TODO replace by counter for outstanding results
	if ( trial.id != -1 ) return true;
	for ( int v = 0; v < mpts; v++ ) {
		if ( simplex[v].id != -1 ) return true;
	}
	return false;
}

void Downhill_Simplex::blubb()
{
	for ( size_t i = 0; i < simplex.size(); i++ ) {
		if ( i > 0 ) { cout << "\t"; }
		for ( size_t d = 0; d < dim; d++ ) {
			cout << simplex[i].x[d] << "|";
		}
		cout << simplex[i].y;
	}
	cout << "\n";
}

} // namespace liee
} // namespace opti
