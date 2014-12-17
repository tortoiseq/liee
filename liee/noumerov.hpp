/*! noumerov.hpp
 */
#ifndef NOUMEROV_HPP_
#define NOUMEROV_HPP_

#include "my_util.hpp"
#include "potential.hpp"
#include "module_config.hpp"
#include "wave_function.hpp"

using namespace std;
namespace liee {

struct Integration_Rec
{
	double E;
	double werror;
	size_t level;
	double xa;              ///< Left classical turning point
	double xb;              ///< Right classical turning point
	double a;               ///< Left integration bound
	double b;               ///< Right integration bound
	double dx;              ///< maximum step size
	double middle;
	bool fixed_bounds;      ///< if true, then keep bounds fixed to assure final convergence
	vector<Point> blend;
	vector<Point> rightwards, leftwards;
	int num_trial;

	SERIALIZE( E & werror & level & xa & xb & a & b & dx & middle & fixed_bounds & num_trial )

	// default constructor
	Integration_Rec( ) : E(0), werror(0), level(0), xa(0), xb(0), a(0), b(0), dx(0),
		middle(0), fixed_bounds(false), num_trial(0) {}

	void set_bounds( double a, double b, double xa, double xb, double dx, bool fixed )
	{
		this->a = a; this->b = b; this->xa = xa; this->xb = xb; this->dx = dx;
		this->fixed_bounds = fixed;
		middle = 0.5 * (xa + xb);
	}

	void clear()
	{
		blend.clear();
		rightwards.clear();
		leftwards.clear();
	}

	bool operator> ( const Integration_Rec& that ) {
		return ( this->E  > that.E );
	}
};

class Noumerov1d : public Module_Exec
{
public:
	Pot_const* pot;
	vector<Integration_Rec> spectrum;
	double Q_up;
	double Q_lo;
	double epsilon;
	double tail_tiny;
	double ttoo_tiny;
	int N_min;
	int N_max;
	int retries;
	int limbo;
	int iteration;
	int max_iterations;
	size_t lvl_lo;
	size_t lvl_up;
	string filename;
	bool is_objective;
	string target_E;

	SERIALIZE( Q_up & Q_lo & epsilon & N_min & N_max & retries & limbo & iteration & max_iterations
	          & tail_tiny & ttoo_tiny & lvl_lo & lvl_up & filename & is_objective & target_E )

	void register_dependencies( vector<Module*> dependencies );
	virtual void initialize( Conf_Module* config, vector<Module*> dependencies );
	virtual void reinitialize( Conf_Module* config, vector<Module*> dependencies );
	virtual void estimate_effort( Conf_Module* config, double & flops, double & ram, double & disk );
	virtual void summarize( map<string, string> & results );
	virtual bool execute();

	double evaluate_energy( Integration_Rec& ir );
	void try_fixate_bounds( Integration_Rec& ir1, Integration_Rec& ir2, Integration_Rec& ir3, Integration_Rec& ir4 );

private:
	void noumerovate( double Q, double a, double m, double b, double dx, vector<Point> & solution );
	double penetrate_border( double Q, double x_turn, double d, double not_less_than );

	void save_graph( string & filename, vector<Point>& data );
	void save_results( string & filename );
	void blend_wf( Integration_Rec& ir );
	void drop_overly_dense_samples( vector<Point> & wf, double portion_to_drop );
	void normalize_area( vector<Point> & wf, double xa, double xb, int sign );
	void scale_to_matching_midpoint( Integration_Rec& ir );
	double mean_square_error( Integration_Rec& ir );
	size_t count_nodes( vector<Point>& wf );
	bool target_missed_levels( double& Q_bottom, double& Q_top );
};

/*!
 * This is a slight adaption of liee::opti::Golden_Section_Search, which adds
 * nothing but the call to try_fixate_bounds() after each iteration. The generic
 * functor is swapped with Noumerov1d object which provides this functionality
 * as well as the objective function evaluate_energy().
 */
class Noum_Golden_Section_Search
{
public:
	double tol;
	Noum_Golden_Section_Search( double tolerance ) : tol( tolerance ) {}

	void minimize( Noumerov1d* provider, Integration_Rec& a, Integration_Rec& b, Integration_Rec& c );
};

}//namespace liee

#endif /* NOUMEROV_HPP_ */
