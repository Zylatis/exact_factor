/////////////////////////////////////////////////////////////////////////////////////
// GLOBAL VARIABLE AND COMPLEX DOUBLE DEFINITIONS
/////////////////////////////////////////////////////////////////////////////////////

#include <omp.h> //for parallelization
#include <iostream>                                                             
#include <fstream>    
#include <ctime>
#include <numeric>
#include <vector>
#include <string>
#include <math.h>
#include <sstream>
#include <set>
#include <sys/time.h>
#include <complex>
#include <cstdlib>
#include <functional>
#include <iomanip>
#include <sys/types.h>
#include <sys/stat.h>

/////////////////////////////////////////////////////////////////////////////////////
// GLOBAL VARS
/////////////////////////////////////////////////////////////////////////////////////

using namespace std;
/////////////////////////////////////////////////////////////////////////////////////
// define complex double
const   complex<double> ii(0.0,1.0);   	
typedef complex<double> dcomp; 

/////////////////////////////////////////////////////////////////////////////////////
// initialize global variables ( values set elsewhere)
double dt, tmax, dr, dR, r_min, r_max, R_min, R_max, init_nd_norm, outer_tol, inner_tol, global_t;

int xx, nr, nR, LS_n, l, deriv_err_order, step_gap;
int perc_done, prev_perc;
int n_states;
int total_steps;
// 1D vector of doubles
vector<double> r_points, R_points, nd_arr;

// 1D vector of ints
vector<int> nd_edges_outer, nd_edges_inner, nd_edges_inner_smooth, nd_edges_outer_smooth;;
vector<int> BO_basis;
vector<int>  grad_coefs_c, grad_coefs_l, grad_coefs_r; // higher stencil order deriv stuff - not used
vector<int>  grad2_coefs_c, grad2_coefs_l, grad2_coefs_r;

// 1D vector of complexes
vector<dcomp> empty_vec_nR, empty_vec_nr;
vector<dcomp> unity_vec_nR, unity_vec_nr;
vector<dcomp> prev_dt_phi_global;


vector<vector< double > > nd_norm_timeseries;
vector<vector< dcomp > > grad_eps, grad_sq_eps;

vector< vector< vector<dcomp> > >  NACV_1, NACV_2, gradNACV_1;

bool have_phi;
bool dump_state;

vector<vector< dcomp > > BOPES;
// string for file stuff
string outputFolder; 
string submodel;
string model_loc;

double m_N;// nuclear mass
/////////////////////////////////////////////////////////////////////////////////////
// initialize global variables ( values set here)
int FILE_P = 15; // precision of outputs
double LAGtol = pow(10.,-45.); // BICSTAB error tol
double t = 0.;
int n_threads;
int step = 0;
int file_n = 0;
vector<double> global_mask;

vector<dcomp> SG_coeffs;
int SGwidth, SGorder;

struct config{

	vector<int> basis;
	vector<double> edge_tols, weights;
	string model, sub_model, sim_name;
	int step_gap, n_threads, init_state;
	
	// false/0 things
	double tmax, dt, R0, k0, sigsq, mask_width, nuc_mass, bdry_mask_edge;
	bool have_phi, A_re, sub_imUen, mask_coeffs, grad_approx, mask_uen, mask_GCOC, mask_dt, ctmqc_uen;
	bool subGauge, free_chi,phi_BE, phi_CN, ft_deriv;
	int fixed_chi;
	
	// true/1 things
	double dir_uen_tdpes;
	bool phi_RK45, re_tdpes, mask_A;
	
	
	// default constructor to initialise stuff
	config() : 
		tmax(0.), dt(0.), R0(0.), k0(0.), sigsq(0.), mask_width(0.), nuc_mass(0.), bdry_mask_edge(0.),
		have_phi(0), A_re(0), sub_imUen(0), mask_coeffs(0), grad_approx(0), mask_uen(0), mask_GCOC(0), mask_dt(0), ctmqc_uen(0), phi_BE(0), phi_CN(0), subGauge(0), free_chi(0),phi_RK45(0),ft_deriv(0),
		
		dir_uen_tdpes(1.), mask_A(1), re_tdpes(0)
	{}
};

config in;


/////////////////////////////////////////////////////////////////////////////////////
//END OF FILE
/////////////////////////////////////////////////////////////////////////////////////
