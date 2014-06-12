#include <cstddef>  //needed in windows build
#include "filesys.h"
#include "boinc_api.h"

#include <boost/foreach.hpp>

#include "my_util.hpp"
#include "wave_function.hpp"


using namespace std;
namespace liee {

void WF_Reader::initialize( Conf_Module* config, vector<Module*> dependencies )
{
	GET_LOGGER( "liee.Module.WF_Reader" );
	Potential* pot;
	for ( size_t i = 0; i < dependencies.size(); i++ ) {
		if ( dependencies[i]->type.compare("potential") == 0 ) {
			pot = dynamic_cast<Potential*>( dependencies[i] );
		}
	}

	string resolved_name;
	boinc_resolve_filename_s( config->get_string("INFILE").c_str(), resolved_name );

	double dr = config->get_double("dr") / CONV_au_nm;
	double r_range = config->get_double("r_range") / CONV_au_nm;
	double r_start = pot->get_r_start();
	int Nr = 1 + r_range / dr;

	//bool complex = config->getParam("complex"]->text.compare("true") == 0;  //TODO(depricated)
	vector<double> parsed;
	vector<Point> wave_funct;
	parse_datafile( resolved_name, parsed );
	if ( parsed.size() == 0 ) {
		LOG_FATAL( "Can't read WF data (file not found): " + resolved_name );
		Except__Preconditions_Fail( __LINE__ );
	}
	double r_min = parsed[0];
	double r_max = parsed[ parsed.size() - 2 ];
	DEBUG_SHOW2(r_min, r_max);
	LOG_INFO( "Read " << parsed.size() / 2 << "WF samples from file" )

	for( size_t i = 0; i < parsed.size()-2; i++ ) {
		if ( i % 2 == 0 ) {
			wave_funct.push_back( Point( parsed[i], parsed[i+1] ) );
		}
	}

	// prepare left tail extrapolation
	double slope_end = abs( (wave_funct[0].y - wave_funct[1].y) / (wave_funct[0].x - wave_funct[1].x) );
	vector<Point> tail;
	Point l = wave_funct[1];
	BOOST_FOREACH( Point p, wave_funct )
	{
		tail.push_back( Point( p.x, log(p.y) ) );
		if ( p.x > l.x ) {
			double slope = abs( (l.y - p.y) / (l.x - p.x) );
			if ( slope > 2 * slope_end  &&  tail.size() > 10 ) break;
		}
		l = p;
	}
	double m_l, n_l, dm_l, dn_l;
	linear_fit( tail, m_l, n_l, dm_l, dn_l);
	DEBUG_SHOW3( tail.size(), m_l, n_l );

	// prepare right tail extrapolation
	int e = wave_funct.size() - 1;
	slope_end = abs( (wave_funct[e].y - wave_funct[e-1].y) / (wave_funct[e].x - wave_funct[e-1].x) );
	tail.clear();
	l = wave_funct[e-1];
	int sign_r_tail = 1;
	if ( wave_funct.back().y < 0 ) {
		sign_r_tail = -1;
	}
	BOOST_REVERSE_FOREACH( Point p, wave_funct )
	{
		tail.push_back( Point( p.x, log( abs(p.y) ) ) );
		if ( p.x < l.x ) {
			double slope = abs( (l.y - p.y) / (l.x - p.x) );
			if ( slope > 2 * slope_end  &&  tail.size() > 10 ) break;
		}
		l = p;
	}

	double m_r, n_r, dm_r, dn_r;
	linear_fit( tail, m_r, n_r, dm_r, dn_r);
	DEBUG_SHOW3( tail.size(), m_r, n_r );
	tail.clear();

	alglib::spline1dinterpolant spline = to_cubic_spline( wave_funct );

	psi.clear();
	for ( int i = 0; i < Nr; i++ ) {
		double r = r_start + i * dr;
		if ( r >= r_min  &&  r <= r_max  ) {
			psi.push_back( dcmplx( alglib::spline1dcalc( spline, r ), 0.0 ) );
		}
		else {  // extrapolate exponential decay of left and right tail
			if ( r < r_min ) {
				psi.push_back( dcmplx( exp( m_l * r + n_l ), 0.0 ) );
			}
			else if ( r > r_max ) {
				psi.push_back( dcmplx( sign_r_tail * exp( m_r * r + n_r ), 0.0 ) );
			}
		}
	}
}

void WF_Reader::reinitialize( Conf_Module* config, vector<Module*> dependencies )
{
	GET_LOGGER( "liee.Module.WF_Reader" );
	// nothing more to do after reload from checkpoint file
}

void WF_Reader::estimate_effort( Conf_Module* config, double & flops, double & ram, double & disk )
{
	double N = config->get_double("r_range") / config->get_double("dr");
	flops += abs( 20 * N );
	ram += abs( sizeof( dcmplx ) * N );
	disk += abs( sizeof( dcmplx ) * N );
}

//-------------------------------------------------------------------------------------------------------------

void WF_Gauss_Packet::initialize( Conf_Module* config, vector<Module*> dependencies )
{
	GET_LOGGER( "liee.Module.WF_Gauss_Packet" );
	size_t l = config->get_array("r0").size();
	for ( size_t i=0; i<l; i++ ) {
		r0.push_back( config->get_array("r0")[i] / CONV_au_nm );
		double E0 = config->get_array("E0")[i] / CONV_au_eV;
		k0.push_back( sign(E0) * sqrt( 2.0 * abs(E0) ) );
		sigma.push_back( config->get_array("sigma")[i] / CONV_au_nm );
	}
	dr = config->get_double("dr") / CONV_au_nm;
	N = 1 + (int)( config->get_double("r_range") / config->get_double("dr") );
	compute_wf();
}

void WF_Gauss_Packet::reinitialize( Conf_Module* config, vector<Module*> dependencies )
{
	GET_LOGGER( "liee.Module.WF_Gauss_Packet" );
	compute_wf();
}

void WF_Gauss_Packet::estimate_effort( Conf_Module* config, double & flops, double & ram, double & disk )
{
	double N = config->get_double("r_range") / config->get_double("dr");
	flops += 50 * N;
	ram += abs( sizeof( dcmplx ) * N );
	disk += abs( sizeof( dcmplx ) * N );
}

void WF_Gauss_Packet::compute_wf()
{
	psi.resize( N, dcmplx( 0.0, 0.0 ) );
	for (size_t j=0; j < k0.size(); j++ ) {
		double normfac = 1.0 / sqrt( sigma[j] * sqrt( CONST_PI ) );
		for ( size_t i = 0; i < N; i++ ) {
			psi[i] += normfac * exp( -0.5 * pow( (i*dr - r0[j]) / sigma[j], 2.0 ) ) * exp( dcmplx( 0.0, k0[j] * i*dr ) );
		}
	}
}


}//namespace liee
