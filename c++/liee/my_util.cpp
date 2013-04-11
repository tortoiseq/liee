/*!
 * my_util.cpp
 *
 *  Created on: 6-Jan-2012
 *      Author: quark
 *
 * Contains a rather unstructured collection of utility routines, which are
 * both simple and generic to help reuse and avoid duplication.
 */

#include <math.h>

#include <sys/stat.h>
#include <archive.h>
#include <archive_entry.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string.h>

#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/asio.hpp>
#include <boost/asio/ip/tcp.hpp>

#include "my_util.hpp"

using namespace std;
using namespace liee;
namespace liee {

Ranq1::Ranq1(unsigned long long seed)
{
	v = 4101842887655102017LL;
	v ^= seed;
	v = int64();
}

BracketedInfixParser::BracketedInfixParser( std::map<string, double> * directory_of_variables )
{
	vars = directory_of_variables;
}

double BracketedInfixParser::evaluate( string expression )
{
	size_t a = expression.find_first_of( '(' );
	size_t b = expression.find_last_of( ')' );
	return eval( strip_white( expression.substr( a+1, b-a-1 ) ) );
}

double BracketedInfixParser::eval( string ex )
{
	double operand[2] = {};

	// first, resolve nested expressions
	int bracket_level = 0;
	size_t start_pos = ex.npos;
	for ( size_t i = 0; i < ex.length(); i++ )
	{
		if ( ex[i] == '(' ) {
			if ( bracket_level == 0 )	{ start_pos = i; }
			bracket_level++;
		}

		if ( ex[i] == ')' )
		{
			bracket_level--;
			if ( bracket_level < 0 ) {
				throw new Except__Preconditions_Fail( 101 );	//TODO give it exceptions of its own
			}

			if ( bracket_level == 0 ) {
				// final closing bracket for this operand. cut the sub-expression and evaluate it recursively
				double x = eval( ex.substr( start_pos+1, i-start_pos-1 ) );
				ex.erase( start_pos, i-start_pos+1 );
				i = start_pos; // reset loop-index to where the removed bracketed expression was

				if ( start_pos == 0 ) {
					operand[0] = x;
					ex.insert( start_pos, "$1" );	//replace by reference to operand-1
				}
				else {
					operand[1] = x;
					ex.insert( start_pos, "$2" );
				}
			}
		}
	}
	if ( bracket_level != 0 ) { // brackets don't match
		throw new Except__Preconditions_Fail( 102 );	//TODO give it exceptions of its own
	}

	// now, there should be exactly one operator left in the string.
	// is it an unary operator?
	if ( ex.find("-") == 0 ) {
		return - eval( ex.substr( 1, ex.npos ), operand );
	}
	if ( ex.find("exp") == 0 ) {
		return exp( eval( ex.substr( 3, ex.npos ), operand ) );
	}
	if ( ex.find("ln") == 0 ) {
		return log( eval( ex.substr( 2, ex.npos ), operand ) );
	}
	if ( ex.find("sin") == 0 ) {
		return sin( eval( ex.substr( 3, ex.npos ), operand ) );
	}
	if ( ex.find("cos") == 0 ) {
		return cos( eval( ex.substr( 3, ex.npos ), operand ) );
	}

	// must be binary then
	size_t op_pos = ex.find_first_of("+-*/^");
	if ( op_pos == ex.npos ) {
		// did not find any supported operator	//TODO support "(sin(10.3))" where the inside bracket has no operator but obviously can be evaluated
		throw new Except__Preconditions_Fail( 103 );	//TODO give it exceptions of its own
	}
	operand[0] = eval( ex.substr( 0, op_pos ), operand );
	operand[1] = eval( ex.substr( op_pos+1, ex.npos ), operand );

	if ( ex.find("-") != ex.npos ) {
		return operand[0] - operand[1];
	}
	if ( ex.find("+") != ex.npos ) {
		return operand[0] + operand[1];
	}
	if ( ex.find("*") != ex.npos ) {
		return operand[0] * operand[1];
	}
	if ( ex.find("/") != ex.npos ) {
		return operand[0] / operand[1];
	}
	if ( ex.find("^") != ex.npos ) {
		return pow( operand[0], operand[1] );
	}

	return 0; // never reach this line
}

/*!
 * expressions evaluated by this method are either numerals or references to variables/operands ($...)
 * might throw boost::bad_lexical_cast
 */
double BracketedInfixParser::eval( string ex, double* operand )
{
	if ( ex[0] == '$' ) {
		ex.erase( 0, 1 );

		if ( ex[0] == '1' )	return operand[0];
		if ( ex[0] == '2' )	return operand[1];

		if ( vars->find( ex ) == vars->end() ) {
			// undefined variable
			throw new Except__Preconditions_Fail( 104 );	//TODO give it exceptions of its own
		}
		else {
		  return vars->operator[]( ex );
		}

	}

	// must be numeral then
	return boost::lexical_cast<double>( ex );
}

string BracketedInfixParser::strip_white( string s )
{
	size_t pos = s.find_first_of(" \t\n");
	while ( pos != s.npos ) {
		s.erase( pos, 1 );
		pos = s.find_first_of(" \t\n");
	}
	return s;
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

double const PI = 4.0 * atan( 1.0 );
double const CONST_2_OVER_SQRT_PI = 2.0 / sqrt( PI );

dcmplx cerf( dcmplx z )
{
	//TODO (perf.) cache the products to sum over
	//TODO (perf.) break summation if last addend is TINY relative to the sum
	int MAX_ITER = 100;
	//double TINY = 1e-15;

	dcmplx z_ = z;
	dcmplx z__ = 0;
	dcmplx z2 = z * z;
	cout << "...in cerf() \t" << MAX_ITER << "\t" << z << "\n";

	for ( int n = 0; n < MAX_ITER; n++ )
	{
		z_ = z;
		for ( double k = 1; k <= n; k++ )
		{
			double a = -(2*k - 1) / ( k * (2*k + 1) );
			z_ *= a * z2;
		}
		z__ += z_;
	}
	return CONST_2_OVER_SQRT_PI * z__;
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

    if (N < 3) {
    	m = n = dm = dn = 0; //TODO for N=2 can give the line parameters with zero error
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
	archive_write_set_compression_gzip( a );
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

string get_html_page( const string& host, const string& item, long timeout )
{
	boost::asio::ip::tcp::iostream stream;
	stream.expires_from_now( boost::posix_time::seconds( timeout ) );
	stream.connect( host, "http" );
	stream << "GET " << item << " HTTP/1.0\r\n";
	stream << "Host: " << host << "\r\n";
	stream << "Accept: */*\r\n";
	stream << "Connection: close\r\n\r\n";
	stream.flush();
	stringstream stuff;
	string content = "";
	if ( stream.bad() ) return content;
	stuff << stream.rdbuf();
	content = stuff.str();
	size_t header_end_pos = content.find( "\r\n\r\n" );
	if ( content.find( "200 OK" ) < header_end_pos ) {
		content.erase( 0, header_end_pos + 4 );
	} else { // not "200 OK"
		content = "";
	}
	return content;
}

string doub2str( double d )
{
	char buff[25];
	sprintf( buff, "%1.15g", d );
	std::string buffAsStdStr = buff;
	return buffAsStdStr;
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

//!
double arithmetic_mean( const vector<double> & x )
{
	return sum( x ) / x.size();
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
			A += 0.5 * ( data.at(i).y + data.at(i-1).y ) * ( data.at(i).x - data.at(i-1).x );	// trapezoid rule
		}
		else if ( data.at(i).x >= b ) break;
	}
	return A;
}


} //namespace liee
