#ifndef WAVEFUNCTION_H_
#define WAVEFUNCTION_H_

#include "my_util.hpp"
#include "potential.hpp"
#include "module_config.hpp"

using namespace std;
namespace liee {

/*!
 * Definition of wave function interface.
 * The only requirement is to grant access to psi.
 */
class Wave_Function
{
public:
	vector<dcmplx> psi;
};


/*!
 * Potential well with memory charge, dc-field and absorbing potential on the right side.
 */
class WF_Reader : public Wave_Function, public Module
{
public:
	WF_Reader() {}
	virtual void initialize( Conf_Module* config, vector<Module*> dependencies );
	virtual void reinitialize( Conf_Module* config, vector<Module*> dependencies );
	virtual void estimate_effort( Conf_Module* config, double & flops, double & ram, double & disk );
	virtual void summarize( map<string, string> & results ){}

	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & ar, const unsigned int version)
	{
		ar & psi;
	}
};

} //namespace liee

#endif /* WAVEFUNCTION_H_ */
