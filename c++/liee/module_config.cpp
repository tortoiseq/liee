#include <iomanip>
#include <time.h>

#include <filesys.h>

#include <boost/foreach.hpp>
#include "boost/tokenizer.hpp"
#include "boost/lexical_cast.hpp"

#include "module_config.hpp"
#include "my_util.hpp"

using namespace std;
namespace liee {

Conf_Param::Conf_Param(TiXmlElement* pParamNode, int parent)
{
	GET_LOGGER( "liee.Conf_Param" );
	// fetch attributes
	parent_id = parent;
	evaluated = false;
	name = pParamNode->Attribute( "name" );
	text = pParamNode->Attribute( "value" );

	if ( text.length() == 0 ) {
		text = pParamNode->Attribute( "default" );
	}
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

	if ( text.length() == 0 ) {
		textual = true;
		return;
	}

	try {
		value = boost::lexical_cast<double>( text );
		textual = false;
	} catch ( boost::bad_lexical_cast& e ) {
		// parameter might be an array (we are only interested in array of doubles)
		if ( text[0] == '{' ) {
			append_array_literal( text, values );
			textual = false;
		} else {
			textual = true;
			// note: also bracketed expressions like '($dr/6)' get flagged as textual here
		}
	}

	// alter expressions such that references "$var_name" to local variables are expanded to the
	// "$parentid::var_name" -syntax, because later on, the locality context is lost in the merged parameter map.
	if ( text[0] == '(' ) {
		size_t i = 0;
		while ( text.find("$", i) != text.npos ) {
			i = text.find("$", i);
			size_t j = text.find("::", i);
			try {
				boost::lexical_cast<int>( text.substr( i+1, j-i-1 ) );
			}
			catch ( const boost::bad_lexical_cast& e ) {
				text.insert( i+1, boost::lexical_cast<string>(parent_id) + "::" );
			}
			i++;
		}
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
			ss << "{" << setprecision( 15 );    //TODO investigate more into portable and exact string representation of doubles
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
	pModuleNode->QueryStringAttribute( "include", &includes );
	pModuleNode->QueryIntAttribute( "stage", &stage );

	// gather parameters
	TiXmlHandle handle = TiXmlHandle( pModuleNode );
	TiXmlElement* pParamNode;
	for( pParamNode = handle.FirstChild( "param" ).Element() ; pParamNode; pParamNode = pParamNode->NextSiblingElement() )
	{
		if ( strcmp(pParamNode->Value(), "param") == 0 ) {
			Conf_Param* p = new Conf_Param( pParamNode, serial );
			//LOG4CXX_DEBUG(logger, "found parameter: " << p->name << "\t" << p->value )
			param[pParamNode->Attribute( "name" )] = p;
		}
	}
}

void Conf_Module::check_param_exists( const char* id  )
{
	if ( param.count( id ) == 0 ) throw Except__bad_config( "Parameter not found in module", name, id );
	if ( param.count( id ) > 1 )  throw Except__bad_config( "Parameter appears more than once in module", name, id );
}

double Conf_Module::get_double( const char* id  )
{
	check_param_exists(id);
	if ( param[id]->textual && param[id]->evaluated == false ) {
		throw Except__bad_config( "Parameter's value is not a double", name, id );
	}
	return param[id]->value;
}

int Conf_Module::get_int( const char* id  )
{
	check_param_exists(id);
	if ( param[id]->textual && param[id]->evaluated == false ) {
		throw Except__bad_config( "Parameters value is not an integer", name, id );
	}
	//TODO sanity check if value is a reasonable integer
	return (int)param[id]->value;
}

void Conf_Module::set_int( const char* id, const int x)
{
	check_param_exists(id);
	param[id]->value = x;
	param[id]->text = boost::lexical_cast<string>(x);
	param[id]->textual = false;
}

string Conf_Module::get_string( const char* id  )
{
	check_param_exists(id);
	return param[id]->text;
}

bool Conf_Module::get_bool( const char* id  )
{
	check_param_exists(id);
	if ( param[id]->text.compare("true") == 0 ) return true;
	if ( param[id]->text.compare("false") == 0 ) return false;

	throw Except__bad_config( "Boolean parameter is neither \"true\" nor \"false\"", name, id );
	return false;
}

bool Conf_Module::param_is_nan( const char* id )
{
	check_param_exists(id);
	return param[id]->textual;
}
vector<double>& Conf_Module::get_array( const char* id )
{
	check_param_exists(id);
	return param[id]->values;
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

Config::Config( string filename )
{
	GET_LOGGER( "liee.Config" );
	infiles.push_back( "liee_parameter.xml" );
	//infiles.push_back( "log4cxx.conf" )
	TiXmlDocument doc( filename.c_str() );

	if ( doc.LoadFile() == false ) {
		LOG_FATAL( "Couldn't load work-unit file: " << filename );
		throw Except__Preconditions_Fail( 6544 );  //TODO have more specific exception
	}
	TiXmlHandle hDoc( &doc );
	TiXmlElement* pModuleNode;
	TiXmlHandle hRoot( 0 );

	pModuleNode = hDoc.FirstChildElement().Element();
	if ( not pModuleNode ) {  // should always have a valid root
		LOG_FATAL( "work-unit is malformed XML: " << filename );
		throw ( Except__Preconditions_Fail( 6545 ) );  //TODO have more specific exception
	}

	this->version = pModuleNode->Attribute( "version" );
	this->project = pModuleNode->Attribute( "project" );
	this->experiment = pModuleNode->Attribute( "experiment" );
	this->wu = pModuleNode->Attribute( "wu" );
	this->exec_chain = pModuleNode->Attribute( "exec_chain" );

	vector<double> chain_members;
	append_array_literal( exec_chain, chain_members );
	for ( size_t i = 0; i < chain_members.size(); i++) { chain.push_back( new Conf_Module() ); }

	hRoot = TiXmlHandle( pModuleNode );
	num_variables = 0;

	for( pModuleNode = hRoot.FirstChild( "module" ).Element() ; pModuleNode; pModuleNode = pModuleNode->NextSiblingElement() )
	{
		if ( strcmp( pModuleNode->Value(), "module" ) == 0 )
		{
			Conf_Module* m = new Conf_Module( pModuleNode );
			LOG_INFO( "Processing config-parameters of #" << m->serial );

			// include parameters from other modules
			if ( m->includes.length() > 0 ) {
				boost::char_separator<char> sep(",");
				boost::tokenizer<boost::char_separator<char> > tokens( m->includes, sep );
				BOOST_FOREACH( const string& tok, tokens ) {
					bool found = false;
					int include_id = boost::lexical_cast<int>(tok);
					for ( size_t i = 0; i < chain_members.size(); i++ ) {
						if ( ( (int)chain_members[i] ) == include_id ) {
							// add included parameters
							found = true;
							map<string, Conf_Param*>::iterator iter;
							for ( iter = chain[i]->param.begin(); iter != chain[i]->param.end(); ++iter ) {
								m->param[iter->first] = iter->second;
							}
						}
					}
					if ( not found ) throw ( Except__Preconditions_Fail( 6546 ) );  //TODO have more specific exception
				}
			}

			// add the module to the specific slot in execution chain
			for ( size_t i = 0; i < chain_members.size(); i++ ) {
				if ( ( (int)chain_members[i] ) == m->serial ) {
					this->chain[i] = m;

					// add the parameters of module m to the merged-perspective TODO prevent multiple additions of included parameters
					map<string, Conf_Param*>::iterator iter;
					for ( iter = m->param.begin(); iter != m->param.end(); ++iter ) {
						Conf_Param* p = iter->second;
						string var_name = boost::lexical_cast<string>( p->parent_id ) + "::" + p->name;
						merged[ var_name ] = p;
						// also keep account of IN and OUT-FILES
						if ( p->text.length() > 0 ) {
							if ( p->name.compare("INFILE") == 0 ) {
								infiles.push_back( p->text );
							}
							else if ( p->name.compare("OUTFILE") == 0 ) {
								outfiles.push_back( p->text );
							}
						}
						if ( not p->fixed ) { num_variables++; }
					}
				}
			}
		}
	}
	for ( size_t i = 0; i < chain.size(); i++ ) {
		if ( chain[i]->serial == -1 ) {
			LOG_FATAL("The "<< i <<"th member of the execution chain was not found in the config file.");
			exit(1);
		}
	}

	evaluate_expressions();
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
	fprintf( f, "%s\n", "Info" );	//TODO stop using c-style printf
	fprintf( f, "%20s\t%s\n", "version", this->version.c_str() );
	fprintf( f, "%20s\t%s\n", "project", this->project.c_str() );
	fprintf( f, "%20s\t%s\n", "experiment", this->experiment.c_str() );
	fprintf( f, "%20s\t%s\n", "workunit", this->wu.c_str() );
	fprintf( f, "%20s\t%g\n", "completed", (double)(time(0)/3600.0/24.0) );

	// write results
	fprintf( f, "%s\n", "Results" );
	map<string, string>::iterator iter_res;
	for ( iter_res = results.begin(); iter_res != results.end(); ++iter_res ) {
		fprintf( f, "%20s\t%s\n", iter_res->first.c_str(), iter_res->second.c_str() );
	}

	// write modules in chain
	for ( size_t i = 0; i < chain.size(); i++ )
	{
		fprintf( f, "\n%d::%s\n", chain[i]->serial, chain[i]->name.c_str() );

		map<string, Conf_Param*>::iterator iter_param;
		for ( iter_param = chain[i]->param.begin(); iter_param != chain[i]->param.end(); ++iter_param )
		{
			if ( iter_param->second->textual || iter_param->second->values.size() > 0 ) {
				fprintf( f, "%20s\t%s\n", iter_param->first.c_str(), iter_param->second->text.c_str() );
			} else {
				fprintf( f, "%20s\t%1.10g\n", iter_param->first.c_str(), iter_param->second->value );
			}
		}
	}
	fclose( f );
}

void Config::evaluate_expressions()
{
	map<string, Conf_Param*>::iterator iter;
	map<string, double> vars;  // a directory of the actual parameter values

	size_t num_unresolved;
	size_t num_unres_last = merged.size();  // circumvent stop condition (pretend there was an improvement from the last loop)

	// reset evaluated-flags
	for ( iter = merged.begin(); iter != merged.end(); ++iter ) {
		iter->second->evaluated = false;
	}

	// repeat maximum as many times as there are parameters, to be able to resolve cascaded dependencies
	for ( size_t i = 0; i < merged.size(); i++ )
	{
		vars.clear();

		// find ready-to-use parameters
		for ( iter = merged.begin(); iter != merged.end(); ++iter ) {
			Conf_Param* p = iter->second;
			if ( not p->textual ) {
				// arrays are only supported as inputs with the syntax "serial::name[int]"
				for ( size_t j; j < p->values.size(); j++ ) {
					string key = iter->first + "[" + boost::lexical_cast<string>(j) + "]";
					vars[key] = p->values[j];  // dump every constant from the array
				}
				if ( p->values.size() == 0 ) {
					vars[iter->first] = p->value; // add the parameter's scalar value to the variables, since it's ready to use
				}
			}
		}

		// find expressions and try to evaluate
		ExpressionParser parser;
		num_unresolved = 0;
		for ( iter = merged.begin(); iter != merged.end(); ++iter )
		{
			Conf_Param* p = iter->second;
			if ( p->textual && p->text[0] == '(' && p->evaluated == false ) {
				try {
					double v = parser.evaluate( p->text, vars );
					// parsing succeeded
					p->evaluated = true;
					p->value = v;  // next loop this parameter is available in vars too
				}
				catch ( Except__Preconditions_Fail& e ) {
					num_unresolved++;
				}
			}
		}

		if ( num_unresolved == 0 ) { // check end condition
			// done with numeric evaluations, now lastly check for strings pointing to other strings or have a variable integer-suffix
			for ( iter = merged.begin(); iter != merged.end(); ++iter )
			{
				Conf_Param* p = iter->second;
				if ( p->textual  && p->text[0] == '$' ) {  // without brackets, a pointer to another parameter can only refer to another string
					string ref = p->text.substr(1, p->text.length() );
					for ( map<string, Conf_Param*>::iterator it = merged.begin(); it != merged.end(); ++it ) {
						if ( ref.compare( it->first ) == 0 ) {
							p->text = it->second->text;
						}
					}
				}

				size_t opos = p->text.find("_$");  // find the syntax: string_$refInteger
				if ( p->textual  &&  opos != p->text.npos ) {
					string key = p->text.substr( opos + 2, p->text.npos );
					if ( p->text.find("::") == p->text.npos ) {
						key = boost::lexical_cast<string>( p->parent_id ) + "::" + key;  // adding the optional "serail::" as required with vars-directory
					}
					string prefix = p->text.substr( 0, opos );
					if ( vars.find( key ) == vars.end() ) continue;  // suffix is not a known variable --> ignore
					p->text = "" + prefix + boost::lexical_cast<string>( (int)vars[key] );
				}
			}
			break;  // the loop for recursive evaluations, since num_unresolved is zero
		}
		if ( not num_unresolved < num_unres_last ) {
			LOG_FATAL( "Could not evaluate all given expressions" );
			exit(1);
		}
		num_unres_last = num_unresolved;
	}
}

vector<Conf_Param*> Config::get_Variables()
{
	vector<Conf_Param*> result;
	map<string, Conf_Param*>::iterator iter;
	for ( iter = merged.begin(); iter != merged.end(); ++iter ) {
		if ( iter->second->fixed == false ) {
			result.push_back( iter->second );
		}
	}
	return result;
}

}// namespace liee
