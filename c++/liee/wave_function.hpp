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
 * Reads the WF and interpolates to the grid of the experiment.
 */
class WF_Reader : public Wave_Function, public Module
{
public:
	virtual void initialize( Conf_Module* config, vector<Module*> dependencies );
	virtual void reinitialize( Conf_Module* config, vector<Module*> dependencies );
	virtual void estimate_effort( Conf_Module* config, double & flops, double & ram, double & disk );
	virtual void summarize( map<string, string> & results ){}

	SERIALIZE( psi )
};


/*!
 * For using an Gaussian shaped electron wave packet as initial state.
 * Initial uncertainty: delta_r0 = sigma/sqrt(2)		delta_k0 = 1/sigma/sqrt(2)
 */
class WF_Gauss_Packet : public Wave_Function, public Module
{
public:
	double r0;     ///< starting expectation value for position
	double k0;     ///< starting expectation value for wavenumber
	double sigma;  ///< width parameter of Gaussian
	size_t N;      ///< number of spatial samples
	double dr;     ///< spatial grid spacing
	virtual void initialize( Conf_Module* config, vector<Module*> dependencies );
	virtual void reinitialize( Conf_Module* config, vector<Module*> dependencies );
	virtual void estimate_effort( Conf_Module* config, double & flops, double & ram, double & disk );
	virtual void summarize( map<string, string> & results ){}
	void compute_wf();

	SERIALIZE( psi & r0 & k0 & sigma & N & dr )
};

} //namespace liee

#endif /* WAVEFUNCTION_H_ */
