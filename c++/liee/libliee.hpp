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

namespace lib_liee {

liee::Potential* my_pot;

extern "C" {
	void init_potential();
	void calc_potential( double r, double t, double & V_real, double & V_imag );
}

}

#endif /* LIBLIEE_HPP_ */
