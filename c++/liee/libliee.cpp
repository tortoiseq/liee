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
    Config* cnf;
    try {
    	cnf = new Config( "liee_parameter.xml" );
    }
    catch ( Except__Preconditions_Fail & e )
    {
    	if ( e.specification_code == 6544 ) {
    		cnf = new Config( "liee_parameter.XML" );	// try if the auto-renamed parameter file is present instead of the default one
    	}
    	else { throw e; }
    }

    // go over the chain like in app_opti_wrap.cpp to assemble the required Modules (but don't execute any tasks)
    vector<Module*> deps;
    for ( int i = 0; i < (int)cnf->chain.size(); i++ )
    {
    	if ( cnf->chain[i]->stage != 1 ) continue; // only concerned with stage-1 only, e.g. tasks for compute-hosts

    	Module* m;
   		m = factory.assemble( cnf->chain[i]->type, cnf->chain[i]->name, cnf->chain[i]->serial );
		m->initialize( cnf->chain[i], deps );
    	deps.push_back( m );

    	if ( cnf->chain[i]->type.compare( "potential" ) == 0 ) {
    		my_pot = dynamic_cast<Potential*>( m );
    	}
    }

    delete cnf;
}

extern "C" void calc_potential( double r, double t, double & V_real, double & V_imag )
{
	dcmplx z = my_pot->V( r, t );
	V_real = z.real();
	V_imag = z.imag();
}

}
