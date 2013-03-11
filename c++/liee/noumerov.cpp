#include <limits>

#include "filesys.h"
#include "boinc_api.h"

#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>

#include "my_util.hpp"
#include "noumerov.hpp"

using std::vector;
namespace liee
{

void Noumerov1d::register_dependencies( vector<Module*> dependencies )
{
	for ( size_t i = 0; i < dependencies.size(); i++) {
		if ( dependencies[i]->type.compare( "potential" ) == 0 ) {
			LOG_INFO( "found a potential to use");
			pot = dynamic_cast<Potential*>( dependencies[i] );
			double V_min = pot->Vr( pot->get_Vmin_pos(), 0.0 );
			if ( Q_lo < V_min ) {
				Q_lo = V_min;
			}
		}
	}
}

void Noumerov1d::initialize( Conf_Module* config, vector<Module*> dependencies )
{
	GET_LOGGER( "liee.Module.Noumerov1d" );
	epsilon = config->getParam("n_eps")->value;
	tail_tiny = config->getParam("tail_tiny_range")->values[0];
	ttoo_tiny = config->getParam("tail_tiny_range")->values[1];
	N_min = (int)config->getParam("num_sample_range")->values[0];
	N_max = (int)config->getParam("num_sample_range")->values[1];
	lvl_lo = (int)config->getParam("energy_levels")->values[0];
	lvl_up = (int)config->getParam("energy_levels")->values[1];
	Q_lo = config->getParam("search_range")->values[0] / CONV_au_eV;
	Q_up = config->getParam("search_range")->values[1] / CONV_au_eV;
	retries = (int)config->getParam("num_retries")->value;
	max_iterations = (int)config->getParam("max_iterations")->value;
	filename = config->getParam("OUTFILE")->text;
	is_objective = config->getParam("is_objective")->text.compare( "true" ) == 0;
	if ( is_objective ) {
		target_E = config->getParam("target_E")->text;
	}
	exec_done = false;
	iteration = 0;
	limbo = 0;
    register_dependencies( dependencies );
}

void Noumerov1d::reinitialize( Conf_Module* config, vector<Module*> dependencies )
{
	GET_LOGGER( "liee.Module.Noumerov1d" );
	register_dependencies( dependencies );
}

void Noumerov1d::estimate_effort( Conf_Module* config, double & flops, double & ram, double & disk )
{
	//double N = config->getParam("r_end"]->value / config->getParam("dr"]->value;
	flops += 0; //TODO(impl.)
	ram += 0; //TODO(impl.)
}

void Noumerov1d::summarize( map<string, string> & results )
{
	results["init_E0"] = doub2str( spectrum[0].E );
	if ( is_objective ) {
		if ( target_E.size() > 0 ) {
			double target = boost::lexical_cast<double>( target_E );
			results["objective"] = doub2str( abs( target - spectrum[0].E ) );
		} else {
			results["objective"] = doub2str( spectrum[0].E );
		}
	}
}

/*!
 * Is the given floating point number d closer than epsilon to the nearest integer?
 */
bool Noumerov1d::is_integer( double d, double epsilon )
{
	double df = abs( floor( d ) - d );
	double dc = abs( ceil( d ) - d );

	if ( df < epsilon  ||  dc < epsilon ) {
		return true;
	}
	else return false;
}

/*!
 * Integrates from inner point a to outer b for estimated energy eigenvalue E
 * (negative dx interchanges inner and outer point)
 */
void Noumerov1d::noumerovate( double Q, double a, double b, vector<Point> & solution )
{
//	DEBUG_BEACON("%");
//	DEBUG_SHOW3( Q, a, b );
	double dx = abs(a - b) / N_min;
	b = a + dx * ceil( (b - a) / dx ); // make sure a and b are separated by dx times an integer
	double cumulative_error;
	double u_next;
	double steps_needed = N_min / 2;
	int steps_took;

	do { // retry loop
		solution.clear();
		cumulative_error = 0;
		steps_needed *= 2.0;
		steps_took = 1;
		double h;
		double h_ = pow( epsilon / 7.2 / steps_needed / pow( (pot->Vr( a, 0 ) - Q) / 6, 3.0 ), 1.0/6 );		//(3.2)

		h = dx; // maximal step size
		while ( h > h_ ) {	// h has to satisfy: dx / 2^i   ;i being an integer
			h /= 2.0;
		}
		if ( b < a ) h *= -1;	// reverse direction with negative step size
		double h2o6 = h * h / 6.0;
		double T_thresh = pow( epsilon / ( 7.2 * steps_needed ), 1.0/3 );		//(3.2)

		double u_last, u_this, T_last, T_this, T_next;
		double x = a + h; //points at 'u_next'
		double homerun = abs( a - b );

		T_this = h2o6 * (pot->Vr( a, 0 ) - Q);        //(2.7);
		u_this = 1.0;
		T_next = h2o6 * (pot->Vr( x, 0 ) - Q);        //(2.7);
		u_next = exp( sqrt( 3.0 * T_this ) + sqrt( 3.0 * T_next ) );    //(7.3)
		cumulative_error += 7.2 * pow( abs( T_next ), 3 );
		solution.push_back( Point( a, u_this ) );
		solution.push_back( Point( x, u_next ) );

		do { // the integration loop
			u_last = u_this;
			T_last = T_this;
			u_this = u_next;
			T_this = T_next;
			x += h; //points at 'u_next'
			steps_took++;

			T_next = h2o6 * (pot->Vr( x, 0 ) - Q);        //(2.7);
			u_next = ((2.0 + 10.0 * T_this) * u_this - (1.0 - T_last) * u_last) / (1.0 - T_next);  //(2.8)
			cumulative_error += 7.2 * pow( abs( T_next ), 3 );
			solution.push_back( Point( x, u_next ) );

			homerun = abs( x - b );
			// save distance from finish     && |T| much smaller than thresh. && not exceeding max step size
			if ( homerun > abs( 4.0 * h ) && abs( T_next ) < T_thresh / 10.0  &&  abs(2.0 * h) <= dx ) {
				// double step size
				h *= 2.0;
				h2o6 *= 4.0;
				T_thresh = pow( epsilon / ( 7.2 * steps_needed ), 1.0/3 );		//(3.2)
				u_this = u_last;
				T_this = 4.0 * T_last;
				T_next *= 4.0;
			}
			// if |T| exceeds save value
			else if ( abs( T_next ) > T_thresh ||
			// or         close to finish      && current step size would miss exact finish point, (but don't get too tiny)
					( (homerun < abs( 2.0 * h) && (not is_integer( homerun / h, 5e-3 )) && (dx / h < 1024) ) ) )
			{
				// half step size
				double x0 = x - 0.5 * h;
				h *= 0.5;
				h2o6 *= 0.25;
				T_thresh = pow( epsilon / ( 7.2 * steps_needed ), 1.0/3 );		//(3.2)
				T_next /= 4.0;
				T_this = h2o6 * (pot->Vr( x0, 0 ) - Q);
				u_this = ( (1.0 - T_next) * u_next + (1.0 - h2o6 * (pot->Vr( x0 - h, 0 ) - Q)) * u_this ) / (2.0 + 10.0 * T_this); //(3.5)
			}

		} while ( homerun > abs( h / 2.0 ) ); //we want to end the loop with x->next ~=~ b and two previous samples
		// here was the slope calculation
	}
	// we _really_ want no error bigger than epsilon, so repeat the whole integration with smaller
	// thresholds until summed up error is under control or the maximum allowed number of samples is reached.
	while ( cumulative_error > epsilon  &&  2 * steps_took < N_max );
}
/* BACKUP of the slope - calculation, which has been dumped in favour of the mean square error minimisation
// here we end with u_next==u(b), u_this==u(b-h), u_last==u(b-2h)
double u_plus_h = ( (2.0 + 10.0 * h2 * T(b) ) * u_next - ( 1.0 - h2 * T(b-h) ) * u_this ) / ( 1.0 - h2*T(b+h) );  		//(2.8)
double u_plus_2h = ( (2.0 + 10.0 * h2 * T(b+h) ) * u_plus_h - ( 1.0 - h2 * T(b) ) * u_next ) / ( 1.0 - h2*T(b+2*h) );  	//(2.8)
double A1 = 0.5 * (u_plus_h - u_this);         								//(4.4)
double A2 = 0.5 * (u_plus_2h - u_last);        								//(4.7)
double B1 = T(x + h) * u_plus_h - T(x - h) * u_this;						//(4.5)
double B2 = T(x + 2.0*h) * u_plus_2h - T(x - 2.0*h) * u_last;				//(4.8)
slope = 16.0/21.0 * ( -A1 + 37.0/32.0*A2 - 37.0/5.0*B1 - 17.0/40.0*B2 ) / h;    //(4.9)
*/


double Noumerov1d::penetrate_border( double Q, double x_turn, double d, double not_less_than )
{
    double A; // area under sqrt(V) beyond turning points
    do {
    	d *= 2.0; //increase penetration to the left until A is sufficiently large (triangle estimation)
    	if ( ( pot->Vr( x_turn + d, 0 ) - Q ) < 0 )
    		A = 0;
    	else
    		A = 0.5 * abs( d ) * sqrt( 2.0 * ( pot->Vr( x_turn + d, 0 ) - Q ) );

    } while ( A < -log( tail_tiny ) );

    double x = x_turn;
    // start integrating (sum) all over again with a step-size of 1/1000 of the estimated penetration depth
    for ( A = 0.0; ; x += d / 1000.0 ) {
    	if ( real( pot->V( x, 0 ) ) < Q ) continue;
    	if ( A > -log( tail_tiny ) ) {
    		if ( d > 0  &&  x > not_less_than ) break;
    		if ( d < 0  &&  x < not_less_than ) break;
    	}
    	A += abs( d ) / 1000.0 * sqrt( 2.0 * ( pot->Vr( x, 0 ) - Q ) );
    }
    // to fulfil the not_less_than condition, outward integration would go further than the allowed underflow ttoo_tiny,
    if ( A > -log( ttoo_tiny ) ) {
    	throw Except__Too_Far_Out( __LINE__ );
    }
    return x;
}

/*!
 * Blends the two solutions so that the exponential tail is taken from only
 * the solution which starts on the same side (has minimal error).
 * The middle is mixed 50%-50%.
 */
void Noumerov1d::blend_wf( Integration_Rec& rec )
{
	alglib::spline1dinterpolant sp_r = to_cubic_spline( rec.rightwards );
	alglib::spline1dinterpolant sp_l = to_cubic_spline( rec.leftwards );

	vector<double> sample_positions;
	BOOST_FOREACH( Point p, rec.leftwards ) {
		sample_positions.push_back( p.x );
	}
	BOOST_FOREACH( Point p, rec.rightwards ) {
		sample_positions.push_back( p.x );
	}

	rec.blend.clear();
	double sigma = 0.1 * abs(rec.xb - rec.xa);

	BOOST_FOREACH( double x, sample_positions ) {
    	double y;
    	double yr = alglib::spline1dcalc( sp_r, x );
    	double yl = alglib::spline1dcalc( sp_l, x );
    	double weight = exp( -pow( (x - rec.middle) / sigma , 2.0) );

    	if ( x < rec.middle && x > rec.xa ) {
    		y = yr * (1.0 - weight) + yl * weight;
    	}
    	else if ( x > rec.middle && x < rec.xb ) {
    		y = yl * (1.0 - weight) + yr * weight;
    	}
    	else if ( x < rec.xa ) {
    		y = yr;
    	}
    	else if ( x > rec.xb ) {
    		y = yl;
    	}

    	rec.blend.push_back( Point( x, y ) );
    }
	std::sort( rec.blend.begin(), rec.blend.end(), less_than_point_x() );
}

/*!
 * Normalise the area between the classical turning points to equal 1.
 */
void Noumerov1d::normalize( vector<Point> & wf, double xa, double xb, int sign )
{
	vector<Point> wf_sqr;
	for ( size_t i = 0; i < wf.size(); i++ ) {
		wf_sqr.push_back( Point( wf[i].x, pow( wf[i].y, 2.0 ) ) );
	}
	alglib::spline1dinterpolant spline = to_cubic_spline( wf_sqr );
	double A0_xa = alglib::spline1dintegrate( spline, xa );
	double A0_xb = alglib::spline1dintegrate( spline, xb );
	double A = A0_xb - A0_xa;

	double norm_factor = sign / sqrt( A );
	for ( size_t i = 0; i < wf.size(); i++ ) {
		wf[i].y *= norm_factor;
	}
}

double Noumerov1d::mean_square_error( Integration_Rec& rec )
{
	alglib::spline1dinterpolant sp1 = to_cubic_spline( rec.leftwards );
	alglib::spline1dinterpolant sp2 = to_cubic_spline( rec.rightwards );

	vector<double> sample_positions;
	BOOST_FOREACH( Point p, rec.leftwards ) {
		if ( p.x > rec.xa  &&  p.x < rec.xb ) {
			sample_positions.push_back( p.x );
		}
	}
	BOOST_FOREACH( Point p, rec.rightwards ) {
		if ( p.x > rec.xa  &&  p.x < rec.xb ) {
			sample_positions.push_back( p.x );
		}
	}

	double result = 0;
	double sigma = 0.2 * abs(rec.xb - rec.xa);

	// weight in favour of centring points
	BOOST_FOREACH( double x, sample_positions ) {
		double weight = exp( -pow( (x - rec.middle) / sigma , 2.0) );
		result += weight * pow( alglib::spline1dcalc( sp1, x ) - alglib::spline1dcalc( sp2, x ), 2.0 );
	}
	return result;
}

double Noumerov1d::evaluate_energy( Integration_Rec& ir )
{
	if ( not ir.fixed_bounds )
	{
		pot->get_outer_turningpoints( ir.E, ir.xa, ir.xb );
		ir.a = penetrate_border( ir.E, ir.xa, 1e-10 * (ir.xa - ir.xb), ir.xa );
		ir.b = penetrate_border( ir.E, ir.xb, 1e-10 * (ir.xb - ir.xa), ir.xb );
		ir.middle = 0.5 * (ir.xa + ir.xb);
	}
	ir.rightwards.clear();
	ir.leftwards.clear();
	int sign = 1;
	noumerovate( ir.E, ir.a, ir.xb, ir.rightwards );
	normalize( ir.rightwards, ir.xa, ir.xb, sign );
	if ( ir.rightwards.back().y < 0 ) {
		sign = -1;
	}
	noumerovate( ir.E, ir.b, ir.xa, ir.leftwards );
	std::sort( ir.leftwards.begin(), ir.leftwards.end(), less_than_point_x() );
	normalize( ir.leftwards, ir.xa, ir.xb, sign );

	blend_wf( ir );
	ir.level = count_nodes( ir.blend );
	ir.werror = mean_square_error( ir );
	return ir.werror;
}

/*!
 * Tries to find integration bounds for both energies Q1 and Q2 such that
 * the wave function's exponential decay on both sides ends up in the interval
 * between ttoo_tiny and tail_tiny. If possible, those bounds are set and
 * kept fixed for the remaining convergence. Otherwise each evaluation will
 * set its own bound to fulfil the tail_tiny requirement.
 */
void Noumerov1d::try_fixate_bounds( Integration_Rec& ir1, Integration_Rec& ir2, Integration_Rec& ir3, Integration_Rec& ir4 )
{
	// if already fixed, keep it that way
	if ( not ir1.fixed_bounds )
	{
		if ( ir1.E > ir3.E || ir2.E > ir3.E ) {
			throw Except__Preconditions_Fail( __LINE__ );
		}
		double x1a, x3a, x1b, x3b;
		double a1, a3, b1, b3;
		// first find integration bounds for E3 which is higher and therefore reach out farther
		pot->get_outer_turningpoints( ir3.E, x3a, x3b );
		a3 = penetrate_border( ir3.E, x3a, 1e-10 * (x3a - x3b), x3a );
		b3 = penetrate_border( ir3.E, x3b, 1e-10 * (x3b - x3a), x3b );

		// then set integration bounds for Q1 under the condition to reach at least as far as for Q2
		pot->get_outer_turningpoints( ir1.E, x1a, x1b );
		try {
			a1 = penetrate_border( ir1.E, x1a, 1e-10 * (x1a - x1b), a3 );
			b1 = penetrate_border( ir1.E, x1b, 1e-10 * (x1b - x1a), b3 );
		}
		catch ( Except__Too_Far_Out & e ) {
			// failure: E1 and E3 are too far apart to have common bounds (and to fixate them)
			return;
		}
		// success
		ir1.set_bounds( a3, b3, x1a, x1b, true );
		ir2.set_bounds( a3, b3, x1a, x1b, true );
		ir3.set_bounds( a3, b3, x1a, x1b, true );
		ir4.set_bounds( a3, b3, x1a, x1b, true );
	}
}

/*!
 * Counts the number of alternating signs of neighbouring points in the given vector.
 * It is assumed that the points are already sorted in regard to x.
 */
size_t Noumerov1d::count_nodes( vector<Point> & wf )
{
	size_t num_nodes = 0;
	for ( size_t i = 0; i < wf.size() - 1; i++ ) {
		if ( wf[i].y * wf[i+1].y < 0 ) {
			num_nodes++;
		}
	}
	return num_nodes;
}

void Noumerov1d::save_graph( string & filename, vector<Point>& data )
{
	FILE* f = boinc_fopen( filename.c_str(), "w" );
	BOOST_FOREACH( Point p, data ) {
		fprintf( f, "%1.16g\t%1.16g\n", p.x, p.y );
	}
	fclose( f );
}

void Noumerov1d::save_results( string & filename )
{
	vector<string> files;
	for ( size_t i = lvl_lo; i <= lvl_up; i++ )
	{
		if ( spectrum[i].num_trial > 0 ) {
			evaluate_energy( spectrum[i - lvl_lo] );

			string si = boost::lexical_cast<string>( i );
			string fi = filename + "_" + string( 4 - si.size(), '0') + si + ".dat"; 	// or printf("%04d", i) ?
			//save_graph( ss0.str() + "_l.dat", spectrum[i - lvl_lo].leftwards );
			//save_graph( ss0.str() + "_r.dat", spectrum[i - lvl_lo].rightwards );
			files.push_back( fi );
			save_graph( fi , spectrum[i - lvl_lo].blend );
			spectrum[i - lvl_lo].clear();
		}
	}

	string fspec = filename + "_spectrum.dat";
	files.push_back( fspec );
	FILE* f = boinc_fopen( fspec.c_str(), "w" );

	for ( size_t i = lvl_lo; i <= lvl_up; i++ ) {
		if ( spectrum[i  -lvl_lo].num_trial > 0 ) {
			fprintf( f, "%d\t%1.16g\t%1.16g\n", (int)i, spectrum[i - lvl_lo].E, spectrum[i - lvl_lo].werror );
		}
	}
	fclose( f );
	tar_gz_files( ".", files, filename );
}

/*!
 * Runs a number of random-initialised minimisations of the mean square error between left- and rightward integration.
 */
bool Noumerov1d::execute()
{
	LOG_INFO( "Searching for Eigenstates in interval [" << Q_lo << " .. " << Q_up << "]" );;
	Noum_Golden_Section_Search golden( 1e-15 );
	size_t lvls = lvl_up - lvl_lo + 1;
	spectrum.resize( lvls );
	double Q_bottom = Q_lo;
	double Q_top = Q_up;

	while ( true )
	{
		Integration_Rec Qa, Qb, Qc;
		Qa.E = Q_bottom + ((double)rand() / (double)RAND_MAX) * (Q_top - Q_bottom);
		Qb.E = Q_bottom + ((double)rand() / (double)RAND_MAX) * (Q_top - Q_bottom);
		Qc.E = Q_bottom + ((double)rand() / (double)RAND_MAX) * (Q_top - Q_bottom);
		sort3( Qa, Qb, Qc);
		try {
			golden.minimize( this, Qa, Qb, Qc );
		}
		catch ( Except__Preconditions_Fail & e ) {
			// evaluation of the three energies could not bracket a minimum --> retry
			continue;
		}

		LOG_DEBUG( "Converged at level " << Qb.level << "with Energy " << Qb.E << " \t remaining difference: " << Qb.werror );;

		if ( Qb.level <= lvl_up  && Qb.level >= lvl_lo )
		{
			int num_trial = spectrum[Qb.level - lvl_lo].num_trial;
			if ( ( spectrum[Qb.level - lvl_lo].num_trial == 0 )	|| spectrum[Qb.level - lvl_lo].werror > Qb.werror )	{
				// found a new one or better one
				spectrum[Qb.level - lvl_lo] = Qb;
				DEBUG_SHOW4( num_trial, Qb.level, Qb.E, Qb.werror );
			}
			DEBUG_SHOW4( iteration, num_trial, spectrum[Qb.level-lvl_lo].werror, Qb.werror);
			spectrum[Qb.level - lvl_lo].num_trial = num_trial + 1;
		}
		else if ( Qb.level > lvl_up ) {
			Q_top = Q_up = Qb.E;
		}
		else if ( Qb.level < lvl_lo ) {
			Q_bottom = Q_lo = Qb.E;
		}

		if ( (not target_missed_levels( Q_bottom, Q_top )) || ++iteration > max_iterations ) {
			// nothing left to target or iteration limit reached
			break;
		}
		//DEBUG_SHOW( iteration );
		DEBUG_SHOW2( Q_bottom, Q_top );
	}

	// save wave functions and free memory
	if ( filename.size() > 0 ) {
		save_results( filename );
	}
	BOOST_FOREACH( Integration_Rec r, spectrum ) {
		r.leftwards.resize( 0 );
		r.rightwards.resize( 0 );
		r.blend.resize( 0 );
	}
	return true;
}

bool Noumerov1d::target_missed_levels( double & Q_bottom, double & Q_top )
{
	//DEBUG_BEACON("*");
	size_t void_size = 0;
	int void_start = -1;
	size_t empty_lvl_count = 0;
	size_t lvls = lvl_up - lvl_lo + 1;

	retarget:;
	for ( size_t i = 0; i < lvls; i++ )
	{
		if ( spectrum[i].num_trial <= limbo ) {
			empty_lvl_count++;
		}
		else if ( empty_lvl_count > 0 ) {
			if ( empty_lvl_count > void_size ) {
				void_size = empty_lvl_count;
				void_start = i - empty_lvl_count;
			}
			empty_lvl_count = 0;
		}
	}
	if ( empty_lvl_count > void_size ) {
		void_size = empty_lvl_count;
		void_start = lvls - empty_lvl_count;
	}

	if ( void_start == -1 ) {
		if ( spectrum[0].num_trial > limbo ) {
			// all levels got at least limbo+1 trials
			if ( ++limbo > retries ) 	return false;
			else 						goto retarget;
		}
		Q_bottom = Q_lo;
		Q_top = Q_up;
	}
	else {
		Q_bottom = ( void_start == 0 ) ? Q_lo : spectrum[void_start - 1].E;
		Q_top = ( void_start + void_size == lvls ) ? Q_up : spectrum[void_start + void_size].E;
		//DEBUG_SHOW4( void_start, void_size, Q_bottom, Q_top );
	}
	//DEBUG_BEACON("#");
	return true;
}

//-----------------------------------------------------------------------------------

void Noum_Golden_Section_Search::minimize( Noumerov1d* provider, Integration_Rec& a, Integration_Rec& b, Integration_Rec& c )
{
//	DEBUG_BEACON("*");
	Integration_Rec d( a );
//	DEBUG_BEACON("?");
	provider->evaluate_energy( a );
//	DEBUG_BEACON(".");
	a.clear();							// clear the data, since for the moment we are only interested in E
//	DEBUG_BEACON("?");
	provider->evaluate_energy( b );
//	DEBUG_BEACON(".");
	b.clear();							// ^-- ...and werror. copy-assigning filled vectors would be wasted effort.
//	DEBUG_BEACON("?");
	provider->evaluate_energy( c );
//	DEBUG_BEACON(".");
	c.clear();

//	DEBUG_BEACON("?");
	provider->try_fixate_bounds( a, b, c, d );
//	DEBUG_BEACON(".");
	double gs = ( 3.0 - sqrt( 5.0 ) ) / 2.0;

	if (not ((a.E < b.E) && (b.E < c.E) && (b.werror < a.werror) && (b.werror < c.werror)) ) {
//		DEBUG_BEACON("-");
		throw Except__Preconditions_Fail( __LINE__ );
	}

//	DEBUG_BEACON("^");
	while ( abs( c.E - a.E ) > tol )
	{
//		DEBUG_SHOW( abs( c.E - a.E ) );
		double ab = abs( b.E - a.E );
		double bc = abs( b.E - c.E );
		if ( ab > bc ) {	// bisect [ab] for it is the bigger interval
			d.E = b.E - gs * ab;
//			DEBUG_BEACON("?");
			provider->evaluate_energy( d );
//			DEBUG_BEACON(".");
			d.clear();
			if ( d.werror < b.werror ) { // test point d is the lowest --> exclude c
				c = b;
				b = d;
			}
			else { // b remains the lowest --> pull in a to d
				a = d;
			}
		}
		else { // bisect interval [bc]
			d.E = b.E + gs * bc;
//			DEBUG_BEACON("?");
			provider->evaluate_energy( d );
//			DEBUG_BEACON(".");
			d.clear();
			if ( d.werror < b.werror ) { // test point d is the lowest --> exclude a
				a = b;
				b = d;
			}
			else { // b remains the lowest --> pull in c to d
				c = d;
			}
		}
//		DEBUG_BEACON("?");
		provider->try_fixate_bounds( a, b, c, d );
//		DEBUG_BEACON(".");
	}
//	DEBUG_BEACON("#");
}

} //namespace liee
