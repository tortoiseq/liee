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

/*!
 *
 */
class Module_Factory
{
public:
	map<string, bool> EXECUTION_REQUIRED;

	Module_Factory() {
		EXECUTION_REQUIRED["main_pot"] 			= false;
		EXECUTION_REQUIRED["gaussian"]			= false;
		EXECUTION_REQUIRED["SLM"] 				= false;
		EXECUTION_REQUIRED["round_well"] 		= false;
		EXECUTION_REQUIRED["pot_lol"] 			= false;
		EXECUTION_REQUIRED["harmonic_osci"] 	= false;
		EXECUTION_REQUIRED["metal_surface"] 	= false;
		EXECUTION_REQUIRED["zero_pot"] 			= false;

		EXECUTION_REQUIRED["wf_reader"] 		= false;
		EXECUTION_REQUIRED["wf_snapshots"] 		= false;
		EXECUTION_REQUIRED["tunnel_ratio"] 		= false;
		EXECUTION_REQUIRED["JWKB_tunnel"] 		= false;
		EXECUTION_REQUIRED["crank_nicholson"] 	= true;
		EXECUTION_REQUIRED["noumerov"] 			= true;
	}

	static Module* assemble( string & type, string & name, int serial )
	{
		return load( type, name, serial, NULL );
	}

	static void store( Module*  m, boost::archive::binary_oarchive* oarch )
	{
		DECLARE_LOGGER
		GET_LOGGER( "liee.Module_Factory" );
		if ( m->type.compare( "potential" ) == 0 && m->name.compare( "main_pot" ) == 0 ) {
			Potential* p = dynamic_cast<Potential*>( m );
			if ( p == NULL ) LOG_ERROR( "DYNAMIC CAST FAILED");
			LOG_DEBUG( "writing as Potential");
			*oarch << p;
		}
		else if ( m->type.compare( "pulse" ) == 0 && m->name.compare( "gaussian" ) == 0 ) {
			Gaussian_Pulse* p = dynamic_cast<Gaussian_Pulse*>( m );
			if ( p == NULL ) LOG_ERROR( "DYNAMIC CAST FAILED");
			LOG_DEBUG( "writing as Gaussian_Pulse");
			*oarch << p;
		}
		else if ( m->type.compare( "pulse" ) == 0 && m->name.compare( "SLM" ) == 0 ) {
			Spatial_Light_Modificator* p = dynamic_cast<Spatial_Light_Modificator*>( m );
			if ( p == NULL ) LOG_ERROR( "DYNAMIC CAST FAILED");
			LOG_DEBUG( "Spatial_Light_Modificator");
			*oarch << p;
		}
		else if ( m->type.compare( "pot_const" ) == 0 && m->name.compare( "round_well" ) == 0 ) {
			Pot_Round_Well_wImage* p = dynamic_cast<Pot_Round_Well_wImage*>( m );
			if ( p == NULL ) LOG_ERROR( "DYNAMIC CAST FAILED");
			LOG_DEBUG( "Pot_Round_Well_wImage");
			*oarch << p;
		}
		else if ( m->type.compare( "pot_const" ) == 0 && m->name.compare( "pot_lol" ) == 0 ) {
			Pot_Experimental* p = dynamic_cast<Pot_Experimental*>( m );
			if ( p == NULL ) LOG_ERROR( "DYNAMIC CAST FAILED");
			LOG_DEBUG( "Pot_Round_Well_wImage");
			*oarch << p;
		}
		else if ( m->type.compare( "pot_const" ) == 0 && m->name.compare( "harmonic_osci" ) == 0 ) {
			Pot_Harm_Oscillator* p = dynamic_cast<Pot_Harm_Oscillator*>( m );
			if ( p == NULL ) LOG_ERROR( "DYNAMIC CAST FAILED");
			LOG_DEBUG( "writing as Pot_Harm_Oscillator");
			*oarch << p;
		}
		else if ( m->type.compare( "pot_const" ) == 0 && m->name.compare( "metal_surface" ) == 0 ) {
			Chulkov_Image_Potential* p = dynamic_cast<Chulkov_Image_Potential*>( m );
			if ( p == NULL ) LOG_ERROR( "DYNAMIC CAST FAILED");
			LOG_DEBUG( "writing as Chulkov_Image_Potential");
			*oarch << p;
		}
		else if ( m->type.compare( "pot_const" ) == 0 && m->name.compare( "zero_pot" ) == 0 ) {
			Zero_Pot* p = dynamic_cast<Zero_Pot*>( m );
			if ( p == NULL ) LOG_ERROR( "DYNAMIC CAST FAILED");
			LOG_DEBUG( "writing as Zero_Potential");
			*oarch << p;
		}
		else if ( m->type.compare( "initial_wf" ) == 0 && m->name.compare( "wf_reader" ) == 0 ) {
			WF_Reader* p = dynamic_cast<WF_Reader*>( m );
			if ( p == NULL ) LOG_ERROR( "DYNAMIC CAST FAILED");
			LOG_DEBUG( "writing as WF_Reader");
			*oarch << p;
		}
		else if ( m->type.compare( "executable" ) == 0 && m->name.compare( "noumerov" ) == 0 ) {
			Noumerov1d* p = dynamic_cast<Noumerov1d*>( m );
			if ( p == NULL ) LOG_ERROR( "DYNAMIC CAST FAILED");
			LOG_DEBUG( "writing as Noumerov");
			*oarch << p;
		}
		else if ( m->type.compare( "observer" ) == 0 && m->name.compare( "wf_snapshots" ) == 0 ) {
			Obs_Snapshot_WF* p = dynamic_cast<Obs_Snapshot_WF*>( m );
			if ( p == NULL ) LOG_ERROR( "DYNAMIC CAST FAILED");
			LOG_DEBUG( "writing as Snapshot_WF");
			*oarch << p;
		}
		else if ( m->type.compare( "observer" ) == 0 && m->name.compare( "tunnel_ratio" ) == 0 ) {
			Obs_Tunnel_Ratio* p = dynamic_cast<Obs_Tunnel_Ratio*>( m );
			if ( p == NULL ) LOG_ERROR( "DYNAMIC CAST FAILED");
			LOG_DEBUG( "writing as Obs_Tunnel_Ratio");
			*oarch << p;
		}
		else if ( m->type.compare( "observer" ) == 0 && m->name.compare( "JWKB_tunnel" ) == 0 ) {
			Obs_JWKB_Tunnel* p = dynamic_cast<Obs_JWKB_Tunnel*>( m );
			if ( p == NULL ) LOG_ERROR( "DYNAMIC CAST FAILED");
			LOG_DEBUG( "writing as Obs_JWKB_Tunnel");
			*oarch << p;
		}
		else if ( m->type.compare( "solver" ) == 0 && m->name.compare( "crank_nicholson" ) == 0 ) {
			Crank_Nicholson* p = dynamic_cast<Crank_Nicholson*>( m );
			if ( p == NULL ) LOG_ERROR( "DYNAMIC CAST FAILED");
			LOG_DEBUG( "writing as Crank_Nicholson");
			*oarch << p;
		}
	}


	static Module* load( string & type, string & name, int serial, boost::archive::binary_iarchive* iarch )
	{
		DECLARE_LOGGER
		GET_LOGGER( "liee.Module_Factory" );
		Module* m;

		if ( type.compare( "potential" ) == 0 && name.compare( "main_pot" ) == 0 ) {
			Potential* p = new Potential();
			if ( iarch != NULL ) *iarch >> p;
			m = p;
		}
		else if ( type.compare( "pulse" ) == 0 && name.compare( "gaussian" ) == 0 ) {
			Gaussian_Pulse* p = new Gaussian_Pulse();
			if ( iarch != NULL ) *iarch >> p;
			m = p;
		}
		else if ( type.compare( "pulse" ) == 0 && name.compare( "SLM" ) == 0 ) {
			Spatial_Light_Modificator* p = new Spatial_Light_Modificator();
			if ( iarch != NULL ) *iarch >> p;
			m = p;
		}
		else if ( type.compare( "pot_const" ) == 0 && name.compare( "round_well" ) == 0 ) {
			Pot_Round_Well_wImage* p = new Pot_Round_Well_wImage;
			if ( iarch != NULL ) *iarch >> p;
			m = p;
		}
		else if ( type.compare( "pot_const" ) == 0 && name.compare( "pot_lol" ) == 0 ) {
			Pot_Experimental* p = new Pot_Experimental;
			if ( iarch != NULL ) *iarch >> p;
			m = p;
		}
		else if ( type.compare( "pot_const" ) == 0 && name.compare( "harmonic_osci" ) == 0 ) {
			Pot_Harm_Oscillator* p = new Pot_Harm_Oscillator();
			if ( iarch != NULL ) *iarch >> p;
			m = p;
		}
		else if ( type.compare( "pot_const" ) == 0 && name.compare( "metal_surface" ) == 0 ) {
			Chulkov_Image_Potential* p = new Chulkov_Image_Potential();
			if ( iarch != NULL ) *iarch >> p;
			m = p;
		}
		else if ( type.compare( "pot_const" ) == 0 && name.compare( "zero_pot" ) == 0 ) {
			Zero_Pot* p = new Zero_Pot();
			if ( iarch != NULL ) *iarch >> p;
			m = p;
		}
		else if ( type.compare( "initial_wf" ) == 0 && name.compare( "wf_reader" ) == 0 ) {
			WF_Reader* p = new WF_Reader();
			if ( iarch != NULL ) *iarch >> p;
			m = p;
		}
		else if ( type.compare( "executable" ) == 0 && name.compare( "noumerov" ) == 0 ) {
			Noumerov1d* p = new Noumerov1d();
			if ( iarch != NULL ) *iarch >> p;
			m = p;
		}
		else if ( type.compare( "observer" ) == 0 && name.compare( "wf_snapshots" ) == 0 ) {
			Obs_Snapshot_WF* p = new Obs_Snapshot_WF();
			if ( iarch != NULL ) *iarch >> p;
			m = p;
		}
		else if ( type.compare( "observer" ) == 0 && name.compare( "tunnel_ratio" ) == 0 ) {
			Obs_Tunnel_Ratio* p = new Obs_Tunnel_Ratio();
			if ( iarch != NULL ) *iarch >> p;
			m = p;
		}
		else if ( type.compare( "observer" ) == 0 && name.compare( "JWKB_tunnel" ) == 0 ) {
			Obs_JWKB_Tunnel* p = new Obs_JWKB_Tunnel();
			if ( iarch != NULL ) *iarch >> p;
			m = p;
		}
		else if ( type.compare( "solver" ) == 0 && name.compare( "crank_nicholson" ) == 0 ) {
			Crank_Nicholson* p = new Crank_Nicholson();
			if ( iarch != NULL ) *iarch >> p;
			m = p;
		}

		m->type = type;
		m->name = name;
		m->serial = serial;
		return m;
	}

};

} //namespace liee

#endif /* MODULE_FACTORY_HPP_ */
