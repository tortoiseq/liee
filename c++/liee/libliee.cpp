#include <iostream>
#include <fstream>

#include "libliee.hpp"
#include "module_factory.hpp"

using namespace std;
using namespace liee;

namespace lib_liee {

void read_config()
{
	if (not config_initialised) {
		try {
			cnf = new Config( "liee_parameter.xml" );
		}
		catch ( Except__Preconditions_Fail & e )
		{
			if ( e.specification_code == 6544 ) {
				cnf = new Config( "liee_parameter.XML" );  // try if the auto-renamed parameter file is present instead of the default one
			}
			else { throw e; }
		}
		config_initialised = true;
	}
}

extern "C" void init_potential()
{
	read_config();
	// go over the chain like in app_opti_wrap.cpp to assemble the required Modules (but don't execute any tasks)
	Module_Factory factory;
	vector<Module*> deps;
	for ( int i = 0; i < (int)cnf->chain.size(); i++ )
	{
		if ( cnf->chain[i]->type.compare("potential") == 0       // so far, only concerned with using potential components
				|| cnf->chain[i]->type.compare("pot_const") == 0
				|| cnf->chain[i]->type.compare("pulse") == 0 )
		{
			Module* m;
			m = factory.assemble( cnf->chain[i]->type, cnf->chain[i]->name, cnf->chain[i]->serial );
			m->initialize( cnf->chain[i], deps );
			deps.push_back( m );

			if ( cnf->chain[i]->type.compare( "potential" ) == 0 ) {
				my_pot = dynamic_cast<Potential*>( m );
			}
		}
	}
}

extern "C" void calc_potential( double r, double t, double & V_real, double & V_imag )
{
	dcmplx z = my_pot->V( r, t );
	V_real = z.real();
	V_imag = z.imag();
}

extern "C" void export_params( char *outfile_name, int module_id )
{
	read_config();
	Module_Factory factory;
	ofstream myfile;
	myfile.open ( outfile_name );
	vector<Module*> deps;
	for ( int i = 0; i < (int)cnf->chain.size(); i++ )
	{
		try {
			Module* m;
			m = factory.assemble( cnf->chain[i]->type, cnf->chain[i]->name, cnf->chain[i]->serial );
			m->initialize( cnf->chain[i], deps );
			deps.push_back( m );
		} catch ( Except__Wrong_Module_Type &e ) {}  // if unknown module type then just don't build it

		if ( cnf->chain[i]->serial == module_id ) {
			map<string, Conf_Param*>::iterator iter;
			for ( iter = cnf->chain[i]->param.begin(); iter != cnf->chain[i]->param.end(); ++iter ) {
				myfile << iter->first << "=" << iter->second->value << "\n";
			}
		}
	}
	myfile.close();
}

extern "C" void get_param( int module_id, char *param_name, char *text, double &value )
{
	read_config();
	for ( int i = 0; i < (int)cnf->chain.size(); i++ ) {

	}

}


}
