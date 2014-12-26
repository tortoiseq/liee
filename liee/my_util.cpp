/*! my_util.cpp
 *
 * Contains a rather unstructured collection of utility routines, which are
 * both simple and generic to help reuse and avoid duplication.
 */

#include <math.h>
#include <sys/stat.h>
#include <archive.h>
#include <archive_entry.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <string.h>

#include "boost/foreach.hpp"
#include "boost/lexical_cast.hpp"
#include "boost/asio.hpp"
#include "boost/regex.hpp"

#include "my_util.hpp"

using namespace std;
using namespace liee;
namespace liee {

ExpressionParser::ExpressionParser()
{
	operators.push_back( Operator("exp", 3, true,  &expOP) );
	operators.push_back( Operator("ln",  2, true,   &lnOP) );
	operators.push_back( Operator("sin", 3, true,  &sinOP) );
	operators.push_back( Operator("cos", 3, true,  &cosOP) );
	operators.push_back( Operator("tan", 3, true,  &tanOP) );
	operators.push_back( Operator("^",   1, false, &powOP) );
	operators.push_back( Operator("*",   1, false, &mulOP) );
	operators.push_back( Operator("/",   1, false, &divOP) );
	operators.push_back( Operator("-",   1, false, &subOP) );
	operators.push_back( Operator("+",   1, false, &addOP) );
}

double ExpressionParser::evaluate( const string &expression, const std::map<string, double> &vars )
{
	string ex = expression;
	strip_white( ex );
	return eval_expression( ex, vars );
}

void ExpressionParser::strip_white( string &s )
{
	size_t pos = s.find_first_of(" \t\n");
	while ( pos != s.npos ) {
		s.erase( pos, 1 );
		pos = s.find_first_of(" \t\n");
	}
}

double ExpressionParser::eval_expression( string ex, const std::map<string, double> &vars )
{
	// first, resolve nested brackets
	int bracket_level = 0;
	size_t start_pos = ex.npos;
	for ( size_t i = 0; i < ex.length(); i++ )
	{
		if ( ex[i] == '(' ) {
			if ( bracket_level == 0 ) { start_pos = i; }
			bracket_level++;
			continue;
		}

		if ( ex[i] == ')' )
		{
			bracket_level--;
			if ( bracket_level < 0 ) {
				throw new Except__Preconditions_Fail( 101 );  //TODO error handling: have more specific exception
			}

			if ( bracket_level == 0 ) {
				// final closing bracket for this operand. cut the sub-expression and evaluate it recursively
				string subex = ex.substr( start_pos+1, i-start_pos-1 );
				ex.erase( start_pos, i-start_pos+1 );
				double x = eval_expression( subex, vars );
				i = start_pos; // reset loop-index to where the removed bracketed expression was
				operands.push_back( x );
				string ref = "$" + boost::lexical_cast<string>( operands.size() - 1 );
				ex.insert( start_pos, ref );  // replace by reference to the operand
			}
		}
	}
	if ( bracket_level != 0 ) { // brackets don't match
		throw new Except__Preconditions_Fail( 102 );  //TODO error handling: have more specific exception
	}
	// now, there should be no more brackets left

	// for the special case of a sign-operator at the start of the expression:  insert a "0"
	if ( (ex.find_first_of("+-") == 0) ) {
		ex = "0" + ex;
	}

	// find first scientific notations, so they don't introduce spurious +/- operators with positive/negative exponents
	boost::regex scinum_regex = boost::regex("(?:0|[1-9]\\d*)(?:\\.\\d*)?(?:[eE][+\\-]?\\d+)");  // does not match numbers without exponent
	vector<string> matches;
	boost::sregex_iterator it( ex.begin(), ex.end(), scinum_regex );
	boost::sregex_iterator end;
	for (; it != end; ++it) {
		matches.push_back( it->str() );
	}
	for ( size_t i = 0; i < matches.size(); i++ ) {
		size_t pos = ex.find( matches[i] );
		double value = eval_operand( matches[i], vars );
		operands.push_back( value );
		string ref = "$" + boost::lexical_cast<string>( operands.size() - 1 );
		ex.erase( pos, matches[i].length() );
		ex.insert( pos, ref );  // replace by reference to the operand
	}

	double result;
	// repeat to evaluate operation-substrings according to given ranking order of operators
	bool noop = true;
	for ( size_t o = 0; o < operators.size(); o++ ) {
		while ( ex.find( operators[o].symbol ) != ex.npos )
		{
			size_t pos = ex.find( operators[o].symbol );
			string operand_str;
			size_t r_pos, l_pos = pos;
			r_pos = copy_adjacent_operand( ex, operand_str, pos + operators[o].sy_len, true, o );
			double r_value = eval_operand( operand_str, vars );

			if ( operators[o].is_unary ) {
				result = (*operators[o].function)( 0.0, r_value );
			} else {
				l_pos = copy_adjacent_operand( ex, operand_str, pos, false, o );
				double l_value = eval_operand( operand_str, vars );
				result = (*operators[o].function)( l_value, r_value );
			}
			ex.erase( l_pos, r_pos - l_pos + 1 );
			operands.push_back( result );
			string ref = "$" + boost::lexical_cast<string>( operands.size() - 1 );
			ex.insert( l_pos, ref );  // replace by reference to the operand
		}
	}
	if (noop) return eval_operand( ex, vars );
	else      return result;
}

size_t ExpressionParser::copy_adjacent_operand( const string &ex, string &copyto, const size_t pos, const bool rightward, const size_t from_op )
{
	size_t op_pos;
	size_t farthest_operand_char;
	if ( rightward ) { farthest_operand_char = ex.length()-1; }
	else             { farthest_operand_char = 0; }

	// find the next operator-string
	for ( size_t o = from_op; o < operators.size(); o++ )
	{
		size_t sy_len = operators[o].sy_len;

		if ( rightward && ( ex.length() - pos > sy_len ) ) {
			op_pos = ex.find( operators[o].symbol, pos );
			if ( op_pos != ex.npos  &&   op_pos-1 < farthest_operand_char ) {
				farthest_operand_char = op_pos-1;
			}
		}
		else if ( not rightward && ( pos > sy_len ) ) {
			op_pos = ex.rfind( operators[o].symbol, pos - sy_len );
			if ( op_pos != ex.npos  &&  op_pos + sy_len > farthest_operand_char ) {
				farthest_operand_char = op_pos + sy_len;
			}
		}
	}

	if ( rightward ) { copyto = ex.substr( pos, farthest_operand_char - pos + 1 ); }
	else             { copyto = ex.substr( farthest_operand_char, pos - farthest_operand_char ); }

	return farthest_operand_char;
}

double ExpressionParser::eval_operand( string &ex, const std::map<string, double> &vars )
{
	// test for $-escaped reference
	if ( ex.find("$") == 0 ) {
		ex.erase(0, 1);
		size_t i = ex.npos;
		try { i = boost::lexical_cast<int>( ex ); } catch ( boost::bad_lexical_cast& e ) {
			// not a previously extracted operand, try out named variables
			if ( vars.find(ex) != vars.end() ) {
				return vars.find(ex)->second;
			}
			// ...negative -> undefined variable
			throw new Except__Preconditions_Fail( 106 );  //TODO error handling: have more specific exception
		}
		if ( i < ex.npos ) {
			if ( i >= operands.size() ) throw new Except__Preconditions_Fail( 107 );
			return operands[i];
		}
	}
	// must be numeral then
	try { return boost::lexical_cast<double>( ex ); } catch ( boost::bad_lexical_cast& e ) {
		throw new Except__Preconditions_Fail( 108 );
	}
}


alglib::spline1dinterpolant to_cubic_spline( vector<Point> & data )
{
	alglib_impl::ae_state my_state;
	alglib_impl::ae_state_init( &my_state );
	alglib_impl::ae_vector vec_x;
	alglib_impl::ae_vector vec_y;
	alglib_impl::ae_vector_init( &vec_x, 0, alglib_impl::DT_REAL, &my_state, true );
	alglib_impl::ae_vector_init( &vec_y, 0, alglib_impl::DT_REAL, &my_state, true );
	alglib_impl::ae_vector_set_length( &vec_x, data.size(), &my_state);
	alglib_impl::ae_vector_set_length( &vec_y, data.size(), &my_state );

	//stringstream ss;
	//ss << "debug-" << time(0) << "-" << clock();
	//FILE* f = fopen( ss.str().c_str(), "w" );
	for( size_t i = 0; i < data.size(); i++ ) {
		vec_x.ptr.p_double[i] = data[i].x;
		vec_y.ptr.p_double[i] = data[i].y;
		//fprintf( f, "%1.8g\t%1.8g\n", data[i].x, data[i].y );
		if ( isnan( data[i].x * data[i].y ) || isinf( data[i].y ) || isinf( data[i].x ) ) {
			/*debug emergency out
			stringstream ss;
			ss << "debug-" << time(0) << "-" << clock();
			FILE* f = fopen( ss.str().c_str(), "w" );
			for( size_t i = 0; i < data.size(); i++ ) {
				fprintf( f, "%1.8g\t%1.8g\n", data[i].x, data[i].y );
			}
			fclose( f );
			//\debug emergency out */
			throw Except__Preconditions_Fail( __LINE__ );
		}
	}
	//fclose( f );

	alglib::real_1d_array arr_x( &vec_x );
	alglib::real_1d_array arr_y( &vec_y );

	int _2ND_DERIVATIVE = 2;
	int N = data.size();

	alglib::spline1dinterpolant spline;
	alglib::spline1dbuildcubic( arr_x, arr_y, N, _2ND_DERIVATIVE, 0.0, _2ND_DERIVATIVE, 0.0, spline );
	alglib_impl::ae_state_clear( &my_state );
	return spline;
}

vector<double> cerf_cache;
dcmplx cerf( dcmplx z, size_t max_iter )
{
	if ( max_iter > cerf_cache.size() ) {
		cerf_cache.resize( max_iter );
		for ( size_t k = 1; k < max_iter; k++ ) {
			cerf_cache[k] = -(2*k - 1) / ( k * (2*k + 1) );
		}
	}
	dcmplx z_ = z;
	dcmplx z__ = 0;
	dcmplx z2 = z * z;

	for ( size_t n = 0; n < max_iter; n++ )
	{
		z_ = z;
		for ( double k = 1; k <= n; k++ ) {
			z_ *= cerf_cache[k] * z2;
		}
		dcmplx test = z__ + z_;
		if ( test ==  z__ ) break;
		z__ = test;
	}
	return M_2_SQRTPI * z__;
}

double me3eep(double y) {
	double x;
	size_t N = (5 + ( (size_t)( 10* drand48() ) ) ) * 1000*1000;
	for ( unsigned long i = 0; i < N; i++ ) {
		x = sin(x + y);
	}
	return x;
}

/*!
 * implementation after http://www.whim.org/nebula/math/lambertw.html
 */
double lambert_w( double x, bool principal )
{
	const double TINY = 1e-10;
	const int MAX_ITER = 1000;
	// set initial value
	double w;

	if ( principal ) {
		if ( x >= -1.0 / exp(1.0)  &&  x <= 10 ) {
			w = 0;
		}
		else if ( x > 10 ) {
			w = log(x) - log( log(x) );
		}
		else {
			throw Except__Preconditions_Fail( __LINE__ );
		}
	}
	else {
		if ( x >= -1.0 / exp(1.0)  &&  x <= -0.1 ) {
			w = -2;
		}
		else if ( x > 0.1  &&  x < 0 ) {
			w = log(-x) - log( -log(-x) );
		}
		else {
			throw Except__Preconditions_Fail( __LINE__ );
		}
	}

	// iterate until converged to difference less than TINY
	double w_new;
	int iter = 0;

	while ( true )
	{
		w_new = ( x * exp(-w) + w * w ) / (w + 1);
		iter++;

		if ( abs(w_new - w) < TINY ) {
			break;
		}
		if ( iter > MAX_ITER ) {
			throw Except__Convergence_Fail( __LINE__ );
		}
		w = w_new;
	}

	return w_new;
}

void linear_fit( const vector<Point>& data, double& m, double& n, double& dm, double& dn )
{
	double x = 0.0;
	double y = 0.0;
	double xx = 0.0;
	double xy = 0.0;
	int N = data.size();

	if (N == 2) {
		m = ( data[1].y - data[0].y ) / ( data[1].x - data[0].x );
		n = data[0].y - m * data[0].x;
		dm = dn = 0;  // this error margins does not include floating-point inaccuracies
		return;
	}
	if (N < 2) {
		m = n = dm = dn = 0;
		return;
	}

	// calculate sums
	BOOST_FOREACH( Point p, data ) {
		x += p.x;
		y += p.y;
		xx += pow( p.x, 2.0 );
		xy += p.x * p.y;
	}

	m = ( N * xy - x * y ) / ( N * xx - x*x );
	n = ( xx * y - x * xy ) / ( N * xx - x*x );

	// calculate error
	double sqrErr = 0;
	BOOST_FOREACH( Point p, data ) {
		sqrErr += pow( m * p.x + n - p.y, 2.0 );
	}

	dn = sqrt( 1.0 / (N - 2) * sqrErr );
	dm = sqrt( dn*dn * N / ( N * xx - x*x ) );
}

//-----------------------------------------------------------------------------------

void tar_gz_files( const string& dir_prefix, const vector<string>& files, const string& archive_name )
{
	struct archive *a;
	struct archive_entry *entry;
	struct stat st;
	char buff[8192];
	int len;
	int fd;

	a = archive_write_new();
	//archive_write_set_compression_gzip( a );
	//archive_write_set_compression_lzma( a );
	archive_write_add_filter_lzma( a );
	archive_write_set_format_pax_restricted( a );
	archive_write_open_filename( a, archive_name.c_str() );

	for ( size_t i = 0; i < files.size(); i++ )
	{
		string fit = dir_prefix + "/" + files[i];
		char filename[256];
		char file_in_tar[256]; // include the wu-name as encompassing directory
		strcpy( filename, files[i].c_str() );
		strcpy( file_in_tar, fit.c_str() );
		stat( filename, &st );
		entry = archive_entry_new();
		archive_entry_set_pathname( entry, file_in_tar );
		archive_entry_set_size( entry, st.st_size );
		archive_entry_set_filetype( entry, AE_IFREG );
		archive_entry_set_perm( entry, 0644 );
		archive_entry_set_atime( entry, st.st_atime, 0 );
		archive_entry_set_ctime( entry, st.st_ctime, 0 );
		archive_entry_set_mtime( entry, st.st_mtime, 0 );
		archive_write_header( a, entry );
		fd = open( filename, O_RDONLY );
		len = read( fd, buff, sizeof( buff ) );
		while ( len > 0 ) {
			archive_write_data( a, buff, len );
			len = read( fd, buff, sizeof( buff ) );
		}
		close( fd );
		archive_entry_free( entry );
	}
	archive_write_close( a );
	//archive_write_free( a ); not declared
}

string doub2str( double d, size_t precision ){
	ostringstream ss;
	ss << std::scientific << std::setprecision( precision ) << d;
	return ss.str();
}


void append_array_literal( const string & str, vector<double> & v_out )
{
	size_t last = 0;
	char * endp;
	const char * tok;
	size_t next = str.find_first_of( ",", last+1 );
	while (next != string::npos)
	{
		tok = str.substr( last+1, next-last-1 ).c_str();
		double value = strtod( tok, &endp );
		// test validity
		if (tok != endp && *endp == '\0') {
			v_out.push_back( value );
		}

		last = next;
		next = str.find_first_of( ",", last+1 );
	}
	tok = str.substr( last+1, str.length()-last-2 ).c_str();
	double value = strtod( tok, &endp );
	// test validity
	if (tok != endp && *endp == '\0') {
		v_out.push_back( value );
	}
}

void parse_datafile( const string filename, vector<double> & data )
{
	ifstream ins( filename.c_str() );
	data.clear();
	string line;

	if (ins.is_open())
	{
		while ( ins.good() )
		{
			getline( ins, line );

			// ignor comments
			int found = line.find_first_of( "#" );
			if ( found != -1 ) {
				line = line.substr( 0, found );
			}
			if ( line.length() == 0 ) {
				continue;
			}

			// tokenize
			char *cstr;
			cstr = new char[line.size() + 1];
			strcpy( cstr, line.c_str() );
			char *tok, *endp;
			double value;
			tok = strtok( cstr, " ,\t" );
			while ( tok != NULL )
			{
				value = strtod( tok, &endp );
				// test validity
				if ( tok != endp && *endp == '\0' ) {
					data.push_back(value);
				}
				// next token
				tok = strtok( NULL, " ,\t" );
			}
		}
		ins.close();
	}
}

//! sums up all doubles in vector
double sum( const vector<double> & x )
{
	double sum = 0;
	for ( size_t i = 0; i < x.size(); i++ ) {
		sum += x[i];
	}
	return sum;
}

double arithmetic_mean( const vector<double> & x )
{
	return sum( x ) / x.size();
}

size_t max_pos( const vector<double> & x )
{
	double max = numeric_limits<double>::min();
	size_t max_i = 0;
	for ( size_t i = 0; i < x.size(); i++ ) {
		if ( x[i] > max ) {
			max = x[i];
			max_i = i;
		}
	}
	return max_i;
}

size_t max_pos( const vector<dcmplx> & x )
{
	double max = numeric_limits<double>::min();
	size_t max_i = 0;
	for ( size_t i = 0; i < x.size(); i++ ) {
		double sqr = real( x[i] * conj( x[i] ) );
		if ( sqr > max ) {
			max = sqr;
			max_i = i;
		}
	}
	return max_i;
}

/*!
 * Calculates the statistic property variance of a given set of numbers.
 */
double variance( const vector<double> & x )
{
	double var = 0;
	double mean = arithmetic_mean( x );
	for ( size_t i = 0; i < x.size(); i++ ) {
		var += pow( x[i] - mean, 2 );
	}
	return var;
}

double find_root( boost::function<double (double)> func, double a, double b, double tol )
{
	double fa = func(a);
	double fb = func(b);
	double m;

	while ( abs( a - b ) / ( abs(a) + abs(b) ) > tol )
	{
		m = 0.5 * (a + b);
		double fm = func(m);
		if ( fa * fm < 0 ) {
			b = m;
			fb = fm;
		}
		else if ( fb * fm < 0 ) {
			a = m;
			fa = fm;
		}
		else if ( fm == 0 ) {
			break;
		}
	}
	return m;
}

double offgrid( double pos, double step )
{
	double c = abs(pos - ceil(pos / step) * step);
	double f = abs(pos - floor(pos / step) * step);
	if ( c < f ) {
		return c;
	}
	else {
		return f;
	}
}

double simple_integrate( const vector<Point> data, double a, double b )
{
	double A = 0;
	for(size_t i = 1; i < data.size(); i++ ) {
		if ( data.at(i).x > a && data.at(i).x < b ) {
			A += 0.5 * ( data.at(i).y + data.at(i-1).y ) * ( data.at(i).x - data.at(i-1).x );  // trapezoid rule
		}
		else if ( data.at(i).x >= b ) break;
	}
	return A;
}


} //namespace liee
