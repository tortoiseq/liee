#include <sstream>
#include <iomanip>
#include <iostream>

#include "filesys.h"

#include <boost/foreach.hpp>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>

#include "module_config.hpp"
#include "my_util.hpp"

using namespace std;
namespace liee {

Conf_Param::Conf_Param(TiXmlElement* pParamNode)
{
	GET_LOGGER( "liee.Conf_Param" );
	// fetch attributes
	name = pParamNode->Attribute( "name" );
	text = pParamNode->Attribute( "value" );
	pParamNode->QueryDoubleAttribute( "lower" , &min );
	pParamNode->QueryDoubleAttribute( "upper" , &max );

	fixed = true;
	logscale = false;
	if ( pParamNode->Attribute( "variable" ) != NULL && strcmp( pParamNode->Attribute( "variable" ), "YES" ) == 0 ) {
		fixed = false;
		if ( pParamNode->Attribute( "logscale" ) != NULL && strcmp( pParamNode->Attribute( "logscale" ), "true" ) == 0 ) {
			logscale = true;
		}
	}

	if ( text.size() == 0 ) {
		textual = true;
		return;
	}

	int check = pParamNode->QueryDoubleAttribute( "value", &value );
	if ( check == TIXML_WRONG_TYPE ) {
		// parameter might be an array (we are only interested in array of doubles)
		if ( text[0] == '{' ) {
			append_array_literal( text, values );
			textual = false;
		} else {
			textual = true;
			// note: also bracketed expressions like '(dr/6)' get only the textual flag at this stage
		}
	} else {
		textual = false;
	}
}

/*!
 * Generates an XML representation for the parameter, which omits ranges and comments.
 * Its supposed to go into the work-unit definition where only name->value pairs are of concern.
 */
TiXmlElement* Conf_Param::minimal_xml_element()
{
	TiXmlElement* param = new TiXmlElement( "param" );
	param->SetAttribute( "name", this->name.c_str() );

	if (textual) {
		param->SetAttribute( "value", this->text.c_str() );
	}
	else {
		if ( this->values.size() > 0 ) {
			//TODO put the conversion in a to_string() function
			ostringstream ss;
			ss << "{" << setprecision( 15 );	//TODO investigate more into portable and exact string representation of doubles
			for ( size_t i = 0; i < values.size(); i++ ) {
				if ( i == values.size() - 1 ) {
					ss << values[i] << "}";
				}
				else {
					ss << values[i] << ",";
				}
			}
			param->SetAttribute( "value", ss.str().c_str() );
		}
		else {
			param->SetDoubleAttribute( "value", this->value );
		}
	}
	return param;
}

Conf_Module::Conf_Module(TiXmlElement* pModuleNode)
{
	GET_LOGGER( "liee.Conf_Module" )
	// fetch attributes
	pModuleNode->QueryIntAttribute( "serial", &serial );
	type = pModuleNode->Attribute( "type" );
	name = pModuleNode->Attribute( "name" );
	pModuleNode->QueryIntAttribute( "stage", &stage );

	// gather parameters
	TiXmlHandle handle = TiXmlHandle( pModuleNode );
	TiXmlElement* pParamNode;
	for( pParamNode = handle.FirstChild( "param" ).Element() ; pParamNode; pParamNode = pParamNode->NextSiblingElement() )
	{
		if ( strcmp(pParamNode->Value(), "param") == 0 ) {
			Conf_Param* p = new Conf_Param( pParamNode );
			//LOG4CXX_DEBUG(logger, "found parameter: " << p->name << "\t" << p->value )
			param[pParamNode->Attribute( "name" )] = p;
		}
	}
}

Conf_Param* Conf_Module::getParam( const string & id )
{
	if ( param.count( id ) == 1 ) {
		return param[id];
	} else {
		LOG_FATAL( "Parameter " << id << " not found in Module " << this->name << "! Exiting..." );
		exit(3);
	}
}

Conf_Param* Conf_Module::getParam( const char* id )
{
	return getParam( string(id) );
}


TiXmlElement* Conf_Module::minimal_xml_element() //TODO the name minimal is misleading: refacture
{
	TiXmlElement* module = new TiXmlElement( "module" );
	module->SetAttribute( "serial", this->serial );
	module->SetAttribute( "type", this->type.c_str() );
	module->SetAttribute( "name", this->name.c_str() );
	module->SetAttribute( "stage", this->stage );

	map<string, Conf_Param*>::iterator iter;
	for ( iter = this->param.begin(); iter != this->param.end(); ++iter ) {
		module->LinkEndChild( iter->second->minimal_xml_element() );
	}
	return module;
}

void Conf_Module::evaluate_expressions()
{
	map<string, Conf_Param*>::iterator iter;
	map<string, double> vars;	// a directory of the actual parameter values

	size_t num_unresolved;
	size_t num_unres_last = param.size();	// to pretend there was an improvement from the last loop

	// repeat maximum as many times as there are parameters, to be able to resolve cascaded dependencies
	for ( size_t i = 0; i < param.size(); i++ )
	{
		vars.clear();

		// find ready-to-use parameters
		for ( iter = param.begin(); iter != param.end(); ++iter ) {
			Conf_Param* p = iter->second;
			if ( (not p->textual) && p->values.size() == 0 ) {	// arrays are not supported
				// add the parameter's value to the variables, since it is ready to use
				vars[iter->first] = p->value;
			}
		}

		// find expressions and try to evaluate
		BracketedInfixParser parser( &vars );
		num_unresolved = 0;
		for ( iter = param.begin(); iter != param.end(); ++iter )
		{
			Conf_Param* p = iter->second;
			if ( p->textual && p->text[0] == '(' ) {
				try {
					double v = parser.evaluate( p->text );
					// parsing succeeded
					p->textual = false;
					p->value = v;		// next loop this parameter is available in vars too
				}
				catch ( Except__Preconditions_Fail& e ) {
					num_unresolved++;
				}
			}
		}

		// check end conditions
		if ( num_unresolved == 0 ) {
			// done with numeric evaluations, now lastly check for strings with variable suffix, for e.g. numbered filenames
			for ( iter = param.begin(); iter != param.end(); ++iter )
			{
				Conf_Param* p = iter->second;
				// identify the hints for the syntax: state_$lvl
				size_t opos = p->text.find("_$");
				if ( p->textual  &&  opos != p->text.npos ) {
					string key = p->text.substr( opos + 2, p->text.npos );
					string prefix = p->text.substr( 0, opos );
					if ( vars.find( key ) == vars.end() ) continue;		// suffix is not a known variable --> ignore
					p->text = prefix + boost::lexical_cast<string>( (int)vars[key] );
					DEBUG_SHOW3( prefix, key, p->text );
				}
			}
			break;
		}
		if ( not num_unresolved < num_unres_last ) {
			LOG_FATAL( "Could not evaluate all given expressions" );
			exit(1);
		}
		num_unres_last = num_unresolved;
	}
}

Config::Config( string filename )
{
	GET_LOGGER( "liee.Config" );
	infiles.push_back( "liee_parameter.xml" );
	//infiles.push_back( "log4cxx.conf" )
	TiXmlDocument doc( filename.c_str() );

	if ( doc.LoadFile() == false ) {
		LOG_FATAL( "Couldn't load work-unit file: " << filename );
		throw Except__Preconditions_Fail( 6544 );	//TODO have more specific exception
	}
	TiXmlHandle hDoc( &doc );
	TiXmlElement* pModuleNode;
	TiXmlHandle hRoot( 0 );

	pModuleNode = hDoc.FirstChildElement().Element();
	if ( not pModuleNode ) {	// should always have a valid root
		LOG_FATAL( "work-unit is malformed XML: " << filename );
		throw ( Except__Preconditions_Fail( 6545 ) );	//TODO have more specific exception
	}

	this->version = pModuleNode->Attribute( "version" );
	this->project = pModuleNode->Attribute( "project" );
	this->experiment = pModuleNode->Attribute( "experiment" );
	this->wu = pModuleNode->Attribute( "wu" );
	this->exec_chain = pModuleNode->Attribute( "exec_chain" );

	vector<double> chain_members;
	append_array_literal( exec_chain, chain_members );
	this->chain.resize( chain_members.size() );

	hRoot = TiXmlHandle( pModuleNode );

	for( pModuleNode = hRoot.FirstChild( "module" ).Element() ; pModuleNode; pModuleNode = pModuleNode->NextSiblingElement() )
	{
		if ( strcmp( pModuleNode->Value(), "module" ) == 0 )
		{
			Conf_Module* m = new Conf_Module( pModuleNode );
			LOG_INFO( "found module " << m->serial << ": " << m->name );

			if ( m->serial == 0 && m->type.compare("global") == 0 ) {
				globals = m->param;

				// add the global parameters (once) to the merged-perspective
				map<string, Conf_Param*>::iterator iter;
				for ( iter = globals.begin(); iter != globals.end(); ++iter ) {
					chain_params_merged.push_back( iter->second );
				}
			}
			else {
				// add the module to the execution chain if requested so
				for ( size_t i = 0; i < chain_members.size(); i++ ) {
					if ( ( (int)chain_members[i] ) == m->serial ) {
						this->chain[i] = m;

						// add the parameters of module m to the merged-perspective (before adding the globals)
						map<string, Conf_Param*>::iterator iter;
						for ( iter = m->param.begin(); iter != m->param.end(); ++iter ) {
							chain_params_merged.push_back( iter->second );
							//also keep account of IN and OUT-FILES
							if ( iter->second->text.length() > 0 ) {
								if ( iter->second->name.compare("INFILE") == 0 ) {
									infiles.push_back( iter->second->text );
								}
								else if ( iter->second->name.compare("OUTFILE") == 0 ) {
									outfiles.push_back( iter->second->text );
								}
							}
						}
					}
				}

				// add the global parameters to each module
				map<string, Conf_Param*>::iterator iter;
				for ( iter = globals.begin(); iter != globals.end(); ++iter ) {
					m->param[iter->first] = iter->second;
				}

				m->evaluate_expressions();
			}
		}
	}

}

void Config::save_file( string & filename )
{
	TiXmlDocument doc; //TODO(perf.) cache the assembled TiXmlDocument for multiple wu-creation with only necessary param.-changes
	TiXmlDeclaration * decl = new TiXmlDeclaration( "1.0", "UTF-8", "" );
	doc.LinkEndChild( decl );

	TiXmlElement* plan9 = new TiXmlElement( "plan9" );
	plan9->SetAttribute( "version", version.c_str() );
	plan9->SetAttribute( "project", project.c_str() );
	plan9->SetAttribute( "experiment", experiment.c_str() );
	plan9->SetAttribute( "wu", wu.c_str() );
	plan9->SetAttribute( "exec_chain", exec_chain.c_str() );

	for ( size_t i = 0; i < chain.size(); i++) {
		plan9->LinkEndChild( chain[i]->minimal_xml_element() );
	}

	doc.LinkEndChild( plan9 );
	doc.SaveFile( filename.c_str() );
}

void Config::save_text( string filename, map<string, string> & results )
{
	FILE* f = boinc_fopen( filename.c_str(), "w" );

	// write results
	fprintf( f, "%s\n", "Results" );
	map<string, string>::iterator iter_res;
	for ( iter_res = results.begin(); iter_res != results.end(); ++iter_res ) {
		fprintf( f, "%20s\t%s\n", iter_res->first.c_str(), iter_res->second.c_str() );
	}

	// write global parameters
	fprintf( f, "\n%s\n", "Globals" );
	map<string, Conf_Param*>::iterator iter;
	for ( iter = globals.begin(); iter != globals.end(); ++iter ) {
		if ( iter->second->text[0] == '(' ) { //TODO code duplication -> if (do write)...
			// print evaluation result instead of expression
			fprintf( f, "%20s\t%1.10g\n", iter->first.c_str(), iter->second->value );
		} else {
			fprintf( f, "%20s\t%s\n", iter->first.c_str(), iter->second->text.c_str() );
		}
	}

	// write modules in chain
	for ( size_t i = 0; i < chain.size(); i++ )
	{
		fprintf( f, "\n%s\n", chain[i]->name.c_str() );

		map<string, Conf_Param*>::iterator iter_param;
		for ( iter_param = chain[i]->param.begin(); iter_param != chain[i]->param.end(); ++iter_param )
		{
			// exclude global parameters again
			bool do_write = true;
			map<string, Conf_Param*>::iterator iter_glob;
			for ( iter_glob = globals.begin(); iter_glob != globals.end(); ++iter_glob )
			{
				if ( iter_glob->first.compare( iter_param->first ) == 0 ) {
					do_write = false;
					break;
				}
			}

			if ( do_write ) {
				if ( iter_param->second->textual || iter_param->second->values.size() > 0 ) {
					fprintf( f, "%20s\t%s\n", iter_param->first.c_str(), iter_param->second->text.c_str() );
				} else {
					fprintf( f, "%20s\t%1.10g\n", iter_param->first.c_str(), iter_param->second->value );
				}
			}
		}
	}
	fclose( f );
}

void Config::reevaluate_expressions()
{
	// first reset textual flag to true for the expressions
	//TODO bad style! coding with side-effects: Conf_Module::evaluate_expression() is setting some Conf_Param::textual to false for internal housekeeping, which has to be undone here. better avoid this side effect in the first place!
	BOOST_FOREACH( Conf_Param* p, chain_params_merged ) {
		if ( p->text[0] == '(' ) {
			p->textual = true;
		}
	}

	BOOST_FOREACH( Conf_Module* m, chain ) {
		m->evaluate_expressions();
	}
}


}// namespace liee
