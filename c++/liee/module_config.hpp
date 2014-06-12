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

struct Except__bad_config {
	string p;
	string m;
	string problem;
	Except__bad_config( string problem, string module="", string parameter="") : p(parameter), m(module), problem(problem) {}
};

class Conf_Param
{
protected:
	DECLARE_LOGGER;
public:
	string  name;
	string  text;
	double  value;
	bool    fixed;
	double  min;
	double  max;
	bool    logscale;
	vector<double> values;
	bool    textual;
	bool    evaluated;      ///< flags successful evaluation of an expression
	int     parent_id;      ///< point to the parent-module to allow cross module references like $[0]wavelength

	/*! Constructor */
	Conf_Param( TiXmlElement * pParamNode, int parent );

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

	/*! Default constructor to create an uninitialised Conf_Module that will mark missing modules*/
	Conf_Module() : serial(-1), stage(-1) {}

	/*! Constructor */
	Conf_Module( TiXmlElement * pmoduleNode );

	/*!
	 * Generates an XML representation for the module, which omits ranges and comments in (children) parameter definitions.
	 * Its supposed to go into the work-unit where only name->value pairs are of concern.
	 */
	TiXmlElement* minimal_xml_element();

	void check_param_exists( const char* id );
	double get_double( const char* id );
	int get_int( const char* id );
	string get_string( const char* id );
	bool get_bool( const char* id );
	vector<double>& get_array( const char* id );
	bool param_is_nan( const char* id );

	void set_int( const char* id, const int x );
};

class Config {
protected:
	DECLARE_LOGGER;
public:
	map<string, Conf_Param*> merged;         ///< parameters of modules in the chain merged and indexed by "module_serial::parameter_name"
	vector<Conf_Module*>     chain;

	string                   version;
	string                   project;
	string                   experiment;
	string                   wu;
	string                   exec_chain;
	vector<string>           infiles;
	vector<string>           outfiles;
	size_t                   num_variables;  ///< number of variable (not fixed) parameters


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
	 * Call this method after altering parameter-values to re-evaluate all expressions.
	 * Some of them might depend on those changes and need to change too.
	 *
	 * Parameters can be declared dependent of other parameters within using math expressions.
	 * This Method tries to evaluate all such expressions and sets the Conf_Param.value accordingly.
	 */
	void evaluate_expressions();

	/*!
	 * returns the not-fixed parameters of this->merged
	 */
	vector<Conf_Param*> get_Variables();
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
