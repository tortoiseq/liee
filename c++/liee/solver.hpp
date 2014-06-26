#ifndef SOLVER_H_
#define SOLVER_H_

#include <vector>

#include "boost/serialization/version.hpp"
#include "boost/serialization/base_object.hpp"
#include "boost/serialization/vector.hpp"

#include "my_util.hpp"
#include "module_config.hpp"
#include "potential.hpp"
#include "observer.hpp"

using namespace std;
namespace liee {

class Solver : public Module_Exec
{
public:
	virtual void evolve_1step() = 0;

	virtual bool execute();
	virtual void renormalize();
	virtual double integrate_psi_sqr();
	virtual double integrate_psi_sqr( double a, double b);

	Potential* potential;
	vector<dcmplx> psi;    ///< wave function
	vector<Observer*> obs; ///< registered observers

	size_t Nr;             ///< number of positions on the grid
	double r_range;        ///< spatial range
	double dr;             ///< step size for r
	double t;              ///< simulation-time
	double t_end;          ///< the end of simulation-time
	double dt;             ///< step size for t
	double dt_;            ///< backup of the actual dt while changing time-steps for adiabatic activation
	int count;
	string outfile;        ///< optionally write final state to this file

	SERIALIZE( psi & Nr & r_range & dr & t & t_end & dt & dt_ & count & exec_done & outfile )
};

class Crank_Nicholson : public Solver
{
public:
	virtual void evolve_1step();
	void register_dependencies( vector<Module*> dependencies );
	virtual void initialize( Conf_Module* config, vector<Module*> dependencies );
	virtual void reinitialize( Conf_Module* config, vector<Module*> dependencies );
	virtual void estimate_effort( Conf_Module* config, double & flops, double & ram, double & disk );
	virtual void summarize( map<string, string> & results ){}
private:
	// compute variables
	int jb;
	vector<dcmplx> alfa;
	vector<dcmplx> gamma;
	vector<dcmplx> g;
	vector<dcmplx> phi;
	vector<dcmplx> d;
	dcmplx c;

	SERIALIZE( boost::serialization::base_object<Solver>( *this ) & jb & c )
};

//------------------------------------------------------------------------------------------------

} //namespace liee
#endif /* SOLVER_H_ */
