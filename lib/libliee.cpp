#include <iostream>
#include <fstream>
#include <iomanip>

#include "boost/filesystem/operations.hpp"
#include "boost/filesystem/path.hpp"
namespace fs = boost::filesystem;

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

extern "C" void export_params( int module_id )
{
	read_config();
	Module_Factory factory;
	ofstream myfile;
	myfile.open ("set_vars.gp");
	myfile << std::scientific << std::setprecision(8);
	vector<Module*> deps;
	for ( int i = 0; i < (int)cnf->chain.size(); i++ )
	{
		if ( cnf->chain[i]->serial == module_id ) {
			map<string, Conf_Param*>::iterator iter;
			myfile << "serial=" << module_id << "\n";
			for ( iter = cnf->chain[i]->param.begin(); iter != cnf->chain[i]->param.end(); ++iter ) {
				if ( iter->second->textual && iter->second->text[0] != '(' ) {
					myfile << iter->first << "=\"" << iter->second->text << "\"\n";
				} else {
					myfile << iter->first << "=" << iter->second->value << "\n";
				}
			}
		}
	}
	myfile.close();
}



extern "C" void plot( int module_id )
{
	read_config();

	for ( int i = 0; i < (int)cnf->chain.size(); i++ ) {
		if ( cnf->chain[i]->type.compare("plot") == 0  &&  ( module_id == -1 || module_id == cnf->chain[i]->serial ) ) {
			// see if either template_name or template_name.gp exists in the local directory, otherwise try the default directory. if no template is found, then just continue
			string template_name = cnf->chain[i]->param["template"]->text;
			fs::path f1( template_name.c_str() );
			fs::path f2( template_name.append(".gp").c_str() );
			fs::path d1( fs::current_path() );
			fs::path d2( fs::path("/usr/local/src/gnuplot") );
			string command = "gnuplot '";

			if      ( fs::exists( d1/f1 ) ) { command.append( (d1/f1).c_str() ); }
			else if ( fs::exists( d1/f2 ) ) { command.append( (d1/f2).c_str() ); }
			else if ( fs::exists( d2/f1 ) ) { command.append( (d2/f1).c_str() ); }
			else if ( fs::exists( d2/f2 ) ) { command.append( (d2/f2).c_str() ); }
			else continue;

			export_params( cnf->chain[i]->serial );
			system( command.append("'").c_str() );
		}
	}
}

}
