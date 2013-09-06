#ifndef CONFIG_H_
#define CONFIG_H_


#include <string>
#include <vector>
#include <map>
#ifdef LOG_ENABLED
	#include "log4cxx/logger.h"
#endif

#include <tinyxml.h>
//#include "../tinyxml/tinyxml.h" OS_X BUILD
#include "my_util.hpp"


using namespace std;
namespace liee {

class Conf_Param
{
protected:
	DECLARE_LOGGER;
public:
	string  name;
	double  value;
	bool    fixed;
	double  min;
	double  max;
	bool    logscale;
	vector<double> values;
	string  text;
	bool    textual;

	/*! Constructor */
	Conf_Param( TiXmlElement * pParamNode );

	/*!
	 * Assembles a minimal xml-element from the stored data, containing only name-value-pairs
	 * and ignores attributes: fixed, lower, higher, because those are not required for a simple evaluation.
	 */
	TiXmlElement* minimal_xml_element();
};

class Conf_Module {
protected:
	DECLARE_LOGGER;
public:
	int serial;
	int stage;
	string type;
	string name;
	map<string, Conf_Param*> param;

	/*! Constructor */
	Conf_Module( TiXmlElement * pmoduleNode );

	/*!
	 * Generates an XML representation for the module, which omits ranges and comments in (children) parameter definitions.
	 * Its supposed to go into the work-unit where only name->value pairs are of concern.
	 */
	TiXmlElement* minimal_xml_element();
	Conf_Param* getParam( const string & id );
	Conf_Param* getParam( const char* id );

	/*!
	 * Parameters can be declared dependent of other parameters within Module+global scope using math expressions.
	 * This Method tries to evaluate all such expressions and set the double Conf_Param.value accordingly.
	 * Should be called after the global parameters have been added to the Module.
	 */
	void evaluate_expressions();
};

class Config {
protected:
	DECLARE_LOGGER;
public:
	map<string, Conf_Param*> globals;
	vector<Conf_Param*>      chain_params_merged;
	vector<Conf_Module*>     chain;

	string                   version;
	string                   project;
	string                   experiment;
	string                   wu;
	string                   exec_chain;
	vector<string>           infiles;
	vector<string>           outfiles;


	/*!
	 * Loads the configuration from the specified XML file.
	 * @see documentation to liee-config.xml syntax.
	 */
	Config(string filename);

	/*!
	 * Writes the configuration with (possibly) altered parameters to a file.
	 * This file is supposed to serve as work-unit for remote solvers.
	 * Comments and unused modules are skipped to save a little bandwidth.
	 */
	void save_file( string & filename );

	/*!
	 * Writes the current configuration to a simple text digest in order to
	 * spare the quick glance from XML syntax and other details.
	 * @param filename	name of the textfile to write to
	 * @param results	optionally add some observer results as name-value-pairs
	 */
	void save_text( string filename, map<string, string> & results );

	/*!
	 * Call this method after altering parameter-values to reevaluate all expressions.
	 * Some of them might depend on those changes and need to change too.
	 */
	void reevaluate_expressions();
};

/*!
 * Interface class for program modules.
 */
class Module
{
protected:
	DECLARE_LOGGER;
public:
	string name;
	string type;
	int serial;  //< id defined in parameter file

	virtual ~Module() {}
	virtual void initialize( Conf_Module* config, vector<Module*> dependencies ) = 0;
	virtual void reinitialize( Conf_Module* config, vector<Module*> dependencies ) = 0;
	virtual void estimate_effort( Conf_Module* config, double & flops, double & ram, double & disk ) = 0;

	/*! Appends name-value pairs of the most _relevant_ data to the result-summary.
	 *  Should be limited to a few values (or none), not extensive data.
	 */
	virtual void summarize( map<string, string> & results ) = 0;
};

/*!
 * Interface class for executable program modules.
 * The executable modules are differentiated from the simpler ones by the fact that they typically take more
 * than a few seconds to reach their final state and therefore might be caught by a compute-suspension in the
 * middle of their work. They need to be prepared to checkpoint and restore in cooperation with the function
 * boinc_time_to_checkpoint().
 */
class Module_Exec : public Module
{
public:
	bool exec_done;

	/*!
	 * @return true if execution was completed, false if boinc-checkpoint-signal was received
	 */
	virtual bool execute() = 0;
};

} // namespace liee

#endif /* CONFIG_H_ */
