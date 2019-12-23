/////////////////////////////////////////////////////////////////////////////////////
// EF Solver
// G. H. Gossel 2017/18
// graeme.gossel@gmail.com
// https://aip.scitation.org/doi/10.1063/1.5090802
/////////////////////////////////////////////////////////////////////////////////////

#include "vars.h"
#include "derivs.h"
//~ #include "fftw_deriv.h"
#include "ints.h"
#include "masks.h"

#include "functions.h"
#include "compNACV.h"
#include "linsol.h"
#include "state_objs.h"
#include "RK45.h"
#include "io.h"

/////////////////////////////////////////////////////////////////////////////////////
//MAIN
/////////////////////////////////////////////////////////////////////////////////////
int main ( int argc, char *argv[] ){

	// Read terminal inputs
	if(argc!=2){
		cout<<"No config file specified, exiting."<<endl;
		exit(0);
	}
	
	SGorder = 2;
	SGwidth = 11;
	// Reads in config file and applies config values
	get_config( in, argv[1]);
	
	// Allocate memory and import initial wfn info and other var declarations
	do_setup();

	n_threads = in.n_threads;
	omp_set_num_threads( n_threads );
	
	// Time integration settings
	dt  = in.dt; //time step (atomic units)
	tmax =  in.tmax;// max time (atomic units)
	step_gap = in.step_gap; // number of iterations between data dumps - IMPORTANT FOR TIME STAMP OF FILES
	LS_n = 1000; // Max number of iterations for linear solver
	m_N = in.nuc_mass;
	//////////////////////////////////////////////////////////////

	// Create state objects
	phi_state phiS;
	chi_state chiS;
	
	// Call initialization routines, passing config info
	phiS.init( in );
	chiS.init( in );
	
	
	if( dt == 0. || tmax == 0. || step_gap == 0. || model_loc == "" ){
		cout<<"Error, invalid/incomplete config file - exiting."<<endl;
		exit(0);
	}
	
	vector<dcomp> chiT, phiT;
	print("#################################################");
	print("		Start simulation:");
	print("#################################################");
	cout<<"Step #"<<"\t"<<"t(au)"<<"\t"<<"nuc. norm dev"<<endl;
	
	update_states( phiS, chiS );
	save_output( phiS, chiS );
	
	vector< vector<int> > mask_edge_table;
	chiS.density();
	mask_edge_table.push_back(nd_edges_outer_smooth);
	
	if(!(in.phi_BE || in.phi_CN)){
		in.phi_RK45 = 1;
	} else {
		in.phi_RK45 = 0;
	}
	
	double wt = get_wall_time();
	int n_flat_points = n_states*nR;
	while(t < tmax){	
		global_t = t;
		
		// solve nuclear system with CN
		if(in.fixed_chi == 0 || in.free_chi == 1){
			chiT = algebraic_solve(chiS, dt, t, chiS.size );
			chiS.vec = chiT;
		}
	
		chiS.density();	
		
		// Toggle between using RK45 or a backward euler/crank nicholson scheme for
		// electronic subsystem (choice between those two done in the phiS class in state_objs.h)
		if(in.phi_RK45){
			RK_solve( phiS, dt );		
		} else {
			phiT = algebraic_solve(phiS, dt, t, n_flat_points );
			phiS.vec = phiT;
			phiS.unroll();				
		}
		update_states(phiS, chiS);
			
		// advance steps
		t+=dt;
		step++;	
		
		// print bits to screen
		cout<<".";
		cout.flush();
		if(step % step_gap == 0 ){
			cout<<endl;
			// write stuff to file
			save_output( phiS, chiS );
			mask_edge_table.push_back(nd_edges_outer_smooth);			
		}
	
	// end time loop
	}
	// write time series of nuclear norm to file
	write_vector2D(nd_norm_timeseries, "nd_norm_series");
	write_vector2D(mask_edge_table, "mask_edge_table");
	// print duration of simulation to screen 
	cout<<get_wall_time()-wt<<endl;
}
/////////////////////////////////////////////////////////////////////////////////////
//END OF FILE
/////////////////////////////////////////////////////////////////////////////////////
