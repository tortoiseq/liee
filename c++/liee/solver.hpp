#ifndef SOLVER_H_
#define SOLVER_H_

#include <vector>

#include <boost/serialization/version.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>

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

	Potential* 			potential;
	vector<dcmplx> 		psi;		///< wave function
	vector<Observer*> 	obs;		///< registered observers

	size_t Nr;						///< number of positions on the grid
	int Nt_adiab;					///< number of time steps spend on the adiabatic activation
	double t_adiab;					///< simulation time in a.u. the adiabatic activation evolves. the window [-1..0) is scaled by dt's accordingly
	double r_range;					///< spatial range
	double dr;						///< step size for r
	double t;						///< simulation-time
	double t_end;					///< the end of simulation-time
	double dt;						///< step size for t
	double dt_;						///< backup of the actual dt while changing time-steps for adiabatic activation
	size_t max_i;					///< for debug, index of initially highest WF^2
	int count;

	friend class boost::serialization::access;
    /*! When the class Archive corresponds to an output archive, the
     *  & operator is defined similar to <<.  Likewise, when the class Archive
     *  is a type of input archive the & operator is defined similar to >>. */
    template<class Archive>
    void serialize( Archive & ar, const unsigned int version )
    {
        ar & psi;
        ar & Nr;
        ar & Nt_adiab;
        ar & t_adiab;
        ar & r_range;
        ar & dr;
        ar & t;
        ar & t_end;
        ar & dt;
        ar & dt_;
        ar & max_i;
        ar & count;
        ar & exec_done;
    }
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

	friend class boost::serialization::access;
    /*! When the class Archive corresponds to an output archive, the
     *  & operator is defined similar to <<.  Likewise, when the class Archive
     *  is a type of input archive the & operator is defined similar to >>. */
    template<class Archive>
    void serialize( Archive & ar, const unsigned int version )
    {
    	// serialize base class information
    	ar & boost::serialization::base_object<Solver>( *this );
    	ar & jb;
    	ar & c;
    }
};


//------------------------------------------------------------------------------------------------

} //namespace liee
#endif /* SOLVER_H_ */
