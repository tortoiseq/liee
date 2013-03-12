#include <fstream>

#include "filesys.h"
#include "boinc_api.h"

#include "my_util.hpp"
#include "solver.hpp"
#include "observer.hpp"

using namespace std;

namespace liee {

void Obs_Snapshot_WF::initialize( Conf_Module* config, vector<Module*> dependencies )
{
	GET_LOGGER( "liee.Obs_Snapshot_WF" );
	// prepare spatial downsampling
	double r_range = config->getParam("r_range")->value / CONV_au_nm;
	double dr = config->getParam("dr")->value / CONV_au_nm;

	double r0 = 0.0;
	double r1 = r_range;
	if ( not config->getParam("obs_r0")->textual ) { r0 = config->getParam("obs_r0")->value / CONV_au_nm; }
	if ( not config->getParam("obs_r1")->textual ) { r1 = config->getParam("obs_r1")->value / CONV_au_nm; }
	ir0 = (int) (r0 / dr);
	ir1 = (int) (r1 / dr);
	DEBUG_SHOW2(ir0, ir1);

	int Nr = 1.0 + ( r1 - r0 ) / dr;
	step_r = floor( Nr / config->getParam("r_samples")->value );
	if ( step_r < 1 ) { step_r = 1; }
	if ( step_r > Nr ) { step_r = Nr; }
	config->getParam("r_samples")->value = Nr;	// save the actual number of samples

	// prepare temporal downsampling
	t_range = config->getParam("t_range")->value / CONV_au_fs;
	dt = config->getParam("dt")->value / CONV_au_fs;
	t0 = 0.0;
	t1 = t_range;
	if ( not config->getParam("obs_t0")->textual ) { t0 = config->getParam("obs_t0")->value / CONV_au_fs; }
	if ( not config->getParam("obs_t1")->textual ) { t1 = config->getParam("obs_t1")->value / CONV_au_fs; }
	counter = 0;
	int Nt = 1.0 + ( t1 - t0 ) / dt;
	step_t = floor( Nt / config->getParam("t_samples")->value );
	if ( step_t < 1 ) { step_t = 1; }
	if ( step_t > Nt ) { step_t = Nt; }
	config->getParam("t_samples")->value = Nt;	// save the actual number of samples

	do_square = config->getParam("square")->text.compare("true") == 0;
	stringstream ss;
	int digits = int( config->getParam("precision")->value );
	if (do_square) {
		ss << "%1." << digits << "g\t";
	} else {
		ss << "%1." << digits << "g,%1." << digits << "g\t";
	}
	format = ss.str();

	boinc_resolve_filename_s( config->getParam("OUTFILE")->text.c_str(), filename );
	do_normalize = config->getParam("normalize")->text.compare("true") == 0;
	// delete previous snapshot file
	FILE* f = boinc_fopen( filename.c_str(), "w" );
	fclose( f );
}

void Obs_Snapshot_WF::reinitialize( Conf_Module* config, vector<Module*> dependencies )
{
	GET_LOGGER( "liee.Obs_Snapshot_WF" );
}

void Obs_Snapshot_WF::estimate_effort( Conf_Module* config, double & flops, double & ram, double & disk )
{
	double N = config->getParam("r_range")->value / config->getParam("dr")->value;
	double Nt = config->getParam("t_range")->value / config->getParam("dt")->value;
	double Nr = config->getParam("r_samples")->value;
	int samples = (int) config->getParam("t_samples")->value;
	bool sqr = config->getParam("square")->text.compare("true") == 0;

	flops += abs( 8 * Nt + samples * N * 20 );	// called Nt times without saving + sample times with saving
	ram += 1024;
	double dsk = (samples + 2) * Nr * 2 * 23 * 1.5; // factor: complex:2, double2str:23, outfile duplicated in gz archive: 1.5
	if ( sqr ) dsk /= 2;
	disk += abs( dsk );
}

void Obs_Snapshot_WF::observe( Module* state )
{
	Solver* s = dynamic_cast<Solver*>( state );
	if ( s->t < t0  ||  s->t > t1 ) return;		// not in observation window
	if ( counter++ % step_t != 0 ) return;		// downsample

   	double fac = 1.0;
   	if ( do_normalize ) {
   		fac = 1.0 / s->integrate_psi_sqr();
   	}

    FILE *file;
	file = boinc_fopen( filename.c_str(), "a" );

   	for ( size_t i = ir0; i <= ir1 && i < s->Nr; i += step_r ) {
		double re, im;
    	if ( do_square ) {
    		re = fac * real( s->psi[i] * conj( s->psi[i] ) );
    		fprintf( file, format.c_str(), re );
    	}
    	else {
    		re = fac * real( s->psi[i] );
    		im = fac * imag( s->psi[i] );
    		fprintf( file, format.c_str(), re, im );
    	}
	}
	fprintf( file, "\n" );
   	fclose( file );
}

//-------------------------------------------------------------------------------------------------------------

void Obs_Tunnel_Ratio::initialize( Conf_Module* config, vector<Module*> dependencies )
{
	GET_LOGGER( "liee.Obs_Tunnel_Ratio" );
	Potential* pot;
	for ( size_t i = 0; i < dependencies.size(); i++ ) {
		if ( dependencies[i]->type.compare( "potential" ) == 0 ) {
			pot = dynamic_cast<Potential*>( dependencies[i] );
		}
	}

	filename = config->getParam("OUTFILE")->text;
	is_objective = config->getParam("is_objective")->text.compare( "true" ) == 0;
	ra = config->getParam("r_a")->value / CONV_au_nm;
	rb = config->getParam("r_b")->value / CONV_au_nm;
	if ( ra < 0 ) ra = 0;
	if ( rb > pot->r_range ) { rb = pot->r_range; }

	// temporal downsampling
	t_range = config->getParam("t_range")->value / CONV_au_fs;
	dt = config->getParam("dt")->value / CONV_au_fs;
	counter = 0;
	int Nt = 1.0 + t_range / dt;
	step_t = floor( Nt / config->getParam("t_samples")->value );
	dt *= step_t;
	if ( step_t < 1 ) { step_t = 1; }
	if ( step_t > Nt ) { step_t = Nt; }
	t_samples = 1 + Nt / step_t;
	config->getParam("t_samples")->value = t_samples;	// save the actual number of samples
//	written = false;
}

void Obs_Tunnel_Ratio::reinitialize( Conf_Module* config, vector<Module*> dependencies )
{
	GET_LOGGER( "liee.Obs_Tunnel_Ratio" );
}

void Obs_Tunnel_Ratio::estimate_effort( Conf_Module* config, double & flops, double & ram, double & disk )
{
	//TODO account for (ra..rb) integration
	double Nr = config->getParam("r_range")->value / config->getParam("dr")->value;
	double Nt = config->getParam("t_range")->value / config->getParam("dt")->value;
	double N =  config->getParam("t_samples")->value;

	flops += abs( 8 * Nt + N * Nr * 8 );	// called Nt times without saving + t_sample times to integrate (psi psi*)
	ram += 1024;
	disk += 20 * N;
}

void Obs_Tunnel_Ratio::summarize( map<string, string> & results )
{
	results["tunnel_ratio"] = doub2str( tunnel_ratio ); //TODO boost conversion
	if ( is_objective ) {
		results["objective"] = doub2str( tunnel_ratio );
	}
}

void Obs_Tunnel_Ratio::observe( Module* state )
{
	if ( counter++ % step_t != 0 ) return;
	Solver* s = dynamic_cast<Solver*>( state );
	psi_sqr.push_back( s->integrate_psi_sqr( ra, rb ) );

	// save after getting the last sample
	if ( (int)psi_sqr.size() >= t_samples-1 ) {	// for the sake of robustness we are willing to write the outfile at last-1 step and overwrite at last step, to make sure it gets written at least once.
											//TODO: check log output, if t_samples calculation is always spot on -> two writes, and remove the additional write out at last-1
		LOG_DEBUG( "Writing Tunnel data at time-step " << psi_sqr.size() << " of " << t_samples );
    	tunnel_ratio = 1.0  -  psi_sqr.back() / psi_sqr.front();
    	FILE* f = boinc_fopen( filename.c_str(), "w" );
    	double last = psi_sqr.front();
    	for ( size_t i = 0; i < psi_sqr.size(); i++ )
    	{
    		double j = ( last - psi_sqr[i] ) / dt;
    		last = psi_sqr[i];
    		fprintf( f, "%1.6g\t%1.15g\t%1.10g\n", i * dt , psi_sqr[i], j );
    	}
    	fprintf( f, "\n" );
    	fclose( f );
    	//written = true;
    }
}

} // namespace liee
