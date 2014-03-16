/*
 * module_factory.hpp
 *
 *  Created on: 29-Jan-2012
 *      Author: quark
 */

#ifndef MODULE_FACTORY_HPP_
#define MODULE_FACTORY_HPP_

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include <map>

#include "my_util.hpp"
#include "module_config.hpp"
#include "potential.hpp"
#include "wave_function.hpp"
#include "noumerov.hpp"
#include "observer.hpp"
#include "solver.hpp"

using namespace std;
namespace liee {

#define STORE_MODULE( _class ) \
		_class* instance = dynamic_cast<_class*>( module ); \
		if ( instance == NULL ) LOG_ERROR( "DYNAMIC CAST FAILED"); \
		*oarch << instance;

#define LOAD_MODULE( _class ) \
		_class* instance = new _class(); \
		if ( iarch != NULL ) \
		*iarch >> instance; \
		module = instance;

#define SWITCH_MODULES( _load_or_store, _type, _name ) \
		if ( _type.compare( "potential" ) == 0 ) \
		{ \
			if      ( _name.compare( "main_pot" ) == 0 )        { _load_or_store( Potential ); } \
		} \
		else if ( _type.compare( "pulse" ) == 0 ) \
		{ \
			if      ( _name.compare( "gaussian" ) == 0 )        { _load_or_store( Gaussian_Pulse ); } \
			else if ( _name.compare( "slm" ) == 0 )             { _load_or_store( Spatial_Light_Modificator ); } \
		} \
		else if ( _type.compare( "pot_const" ) == 0 ) \
		{ \
			if      ( _name.compare( "round_well" ) == 0 )      { _load_or_store( Pot_Round_Well_wImage ); } \
			else if ( _name.compare( "pot_lol" ) == 0 )         { _load_or_store( Pot_Experimental ); } \
			else if ( _name.compare( "harmonic_osci" ) == 0 )   { _load_or_store( Pot_Harm_Oscillator ); } \
			else if ( _name.compare( "metal_surface" ) == 0 )   { _load_or_store( Chulkov_Image_Potential ); } \
			else if ( _name.compare( "piecewise_linear" ) == 0 ){ _load_or_store( Pot_Piecewise ); } \
		} \
		else if ( _type.compare( "initial_wf" ) == 0 ) \
		{ \
			if      ( _name.compare( "wf_reader" ) == 0 )       { _load_or_store( WF_Reader ); } \
			else if ( _name.compare( "wave_packet" ) == 0 )     { _load_or_store( WF_Gauss_Packet); } \
		} \
		else if ( _type.compare( "executable" ) == 0 ) \
		{ \
			if      ( _name.compare( "numerov" ) == 0 )         { _load_or_store( Noumerov1d ); } \
		} \
		else if ( _type.compare( "observer" ) == 0 ) \
		{ \
			if      ( _name.compare( "wf_snapshots" ) == 0 )    { _load_or_store( Obs_Snapshot_WF ); } \
			else if ( _name.compare( "jwkb_tunnel" ) == 0 )     { _load_or_store( Obs_JWKB_Tunnel ); } \
			else if ( _name.compare( "wigner_phasespace" ) == 0 ){ _load_or_store( Obs_Wigner_Distribution ); } \
			else if ( _name.compare( "probability_current" ) == 0 ){ _load_or_store( Obs_Probability_Current ); } \
		} \
		else if ( _type.compare( "solver" ) == 0 ) \
		{ \
			if      ( _name.compare( "crank_nicholson" ) == 0 ) { _load_or_store( Crank_Nicholson ); } \
		} \
		else { \
			throw Except__Wrong_Module_Type( __LINE__ ); \
		}


struct Except__Wrong_Module_Type {
	int specification_code;
	Except__Wrong_Module_Type( int c ) : specification_code( c ) {}
};

/*!
 *
 */
class Module_Factory
{
public:
	map<string, bool> EXECUTION_REQUIRED;

	Module_Factory() {
		EXECUTION_REQUIRED["main_pot"]          = false;
		EXECUTION_REQUIRED["gaussian"]          = false;
		EXECUTION_REQUIRED["slm"]               = false;
		EXECUTION_REQUIRED["round_well"]        = false;
		EXECUTION_REQUIRED["pot_lol"]           = false;
		EXECUTION_REQUIRED["harmonic_osci"]     = false;
		EXECUTION_REQUIRED["metal_surface"]     = false;
		EXECUTION_REQUIRED["piecewise_linear"]  = false;
		EXECUTION_REQUIRED["wf_reader"]         = false;
		EXECUTION_REQUIRED["wave_packet"]       = false;
		EXECUTION_REQUIRED["wf_snapshots"]      = false;
		EXECUTION_REQUIRED["jwkb_tunnel"]       = false;
		EXECUTION_REQUIRED["wigner_phasespace"] = false;
		EXECUTION_REQUIRED["probability_current"] = false;
		EXECUTION_REQUIRED["crank_nicholson"]   = true;
		EXECUTION_REQUIRED["numerov"]           = true;
	}

	static Module* assemble( string & type, string & name, int serial )
	{
		return load( type, name, serial, NULL );
	}

	static void store( Module*  module, boost::archive::binary_oarchive* oarch )
	{
		DECLARE_LOGGER;
		GET_LOGGER( "liee.Module_Factory" );
		SWITCH_MODULES( STORE_MODULE, module->type, module->name );
	}


	static Module* load( string & type, string & name, int serial, boost::archive::binary_iarchive* iarch )
	{
		DECLARE_LOGGER
		GET_LOGGER( "liee.Module_Factory" );
		Module* module;
		SWITCH_MODULES( LOAD_MODULE, type, name );
		module->type = type;
		module->name = name;
		module->serial = serial;
		return module;
	}

};

} //namespace liee

#endif /* MODULE_FACTORY_HPP_ */
