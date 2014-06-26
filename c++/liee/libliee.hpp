/*
 * libliee.hpp
 *
 *  Created on: Oct 31, 2012
 *      Author: quark
 *
 *  Interface to provide some essential stuff to external programs like the Python plotters.
 */

#ifndef LIBLIEE_HPP_
#define LIBLIEE_HPP_

#include "my_util.hpp"
#include "potential.hpp"
#include "module_config.hpp"


namespace lib_liee {

liee::Potential* my_pot;
liee::Config* cnf;
bool config_initialised = false;

void read_config();

extern "C" {

	/*!
	 * Precondition: a valid parameter-file "liee_parameter.xml" must be present in the working directory!
	 * The potential will be setup as defined in the parameter-file.
	 * @returns a pointer to the complex number (two doubles) in which the result of a calc_potential call gets written.
	 */
	void init_potential();

	void calc_potential( double r, double t, double & V_real, double & V_imag );

	/*!
	 * This function provides access to the parameter file's values in a gnuplot-friendly form.
	 * For convenience it writes the data as name=value pairs to a file of given name. Normally this
	 * file is used as gnuplot-macro, which provides variable initialisation to the main routine.
	 *
	 * @param[outfile_name]         the filename to write to
	 * @param[module_id]            specify the module by its "serial" attribute in liee_parameter.xml
	 */
	void export_params( char *outfile_name, int module_id );

	void get_param( int module_id, char *param_name, char *text, double &value );
}

}

#endif /* LIBLIEE_HPP_ */
