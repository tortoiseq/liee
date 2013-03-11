#include "libliee.hpp"
#include "module_config.hpp"
#include "module_factory.hpp"

using namespace std;
using namespace liee;

namespace lib_liee {


/*!
 * Precondition: a valid parameter-file "liee_parameter.xml" must be present in the working directory!
 * The potential will be setup as defined in the parameter-file.
 * @returns a pointer to the complex number (two doubles) in which the result of a calc_potential call gets written.
 */
extern "C" void init_potential()
{
    Module_Factory factory;
    Config* config;
    try {
    	config = new Config( "liee_parameter.xml" );
    }
    catch ( Except__Preconditions_Fail & e )
    {
    	if ( e.specification_code == 6544 ) {
    		config = new Config( "liee_parameter.XML" );	// try if the auto-renamed parameter file is present instead of the default one
    	}
    	else { throw e; }
    }

    for ( int i = 0; i < (int)config->chain.size(); i++ )
    {
    	// find the potential, assemble and initialise
    	if ( config->chain[i]->type.compare( "potential" ) == 0 ) {
    		Module* m = factory.assemble( config->chain[i]->type, config->chain[i]->name );
    		vector<Module*> deps;	// empty since not used in potential
    		m->initialize( config->chain[i], deps );
    		my_pot = dynamic_cast<Potential*>( m );
    	}
    }
    delete config;
}

extern "C" void calc_potential( double r, double t, double & V_real, double & V_imag )
{
	dcmplx z = my_pot->V( r, t );
	V_real = z.real();
	V_imag = z.imag();
}

}
