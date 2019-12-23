/////////////////////////////////////////////////////////////////////////////////////
// DEFINE ELECTRONIC AND NUCLEAR STATE CLASSES
/////////////////////////////////////////////////////////////////////////////////////
//~ #include "fftw_deriv.h"

/////////////////////////////////////////////////////////////////////////////////////
// Forward declarations stuff
// This function is used by the nuclear state to fit 'in the wings' to get smooth GCOC
vector<double> do_fit( double ** data);
	

/////////////////////////////////////////////////////////////////////////////////////
// Electronic wfn state class
class phi_state {
	
	public:
		// name of state for debugging if passed as un-named state
		string name;
		
		// 1D vectors for the vector potential, its derivative, grad chi/chi, the TDPES
		vector<dcomp> Avec, divAvec, GCOC, TDPES, G2COC, divAtemp, stability_track;
			
		// 2D vectors to store the C_i(R) vectors, their spatial derivatives, and their time derivatives
		vector<vector< dcomp> > coeff_set, grad_coeff_set, grad_sq_coeff_set, dt_C_set;
		
		// 2D vectors to store the BOPES and norms of each state
		vector< vector<double> > BO_eps, coeff_sq;
		
		// 1D vector to store the norm at each 'R'
		vector<double> elec_norm, local_mask, Acheck;
		
		vector<dcomp> vec, b;

		vector<dcomp> exp_Uen, gauge_cond, ft_grad, ft_grad_sq;
		// 2D vectors to store Uen*C_i and some intermediate vectors for debugging
		vector<vector< dcomp> >  UenC,T1_vec,T2_vec, T3_vec, UenC_trunc;
		
		// Length of all 1D vectors being considerd (nR basically)
		int size;
		
		double LAGtol;
		
		
		///////////////////////
		// Constructor with attribute initialzations - no init here because of RK45 objects
		phi_state(): Avec(nR,0) { }
		
		///////////////////////
		void init( config in ){
			cout<<Avec.size()<<endl;
			name = "phi";
			size = nR; // remember this is the length of the 1D vectors and so is number of nuclear grid points (no nr here, in BO basis!)
			
			
			// Appropriately size (allocate memory) for vectors of length nR
			vector<dcomp> row;
			row.resize(nR); // define row of length nR to shove everywhere
			fill(row.begin(),row.end(),0.); // fill with zeros just cuz
			vec.resize( n_states * nR );
			b.resize( n_states * nR );
			// make a list of addresses to the vectors we want to resize
			vector<vector<dcomp>* > vecs1D = {&Avec, &divAvec, &GCOC, &TDPES, &G2COC, &divAtemp, &exp_Uen, &stability_track, &ft_grad, &ft_grad_sq };
			vector< vector<vector<dcomp> >* > vecs2D = { &coeff_set, &grad_coeff_set, &grad_sq_coeff_set, &dt_C_set,  &UenC, &T1_vec, &T2_vec, &T3_vec, &UenC_trunc  };
			
			// go through and resize all the bits. Don't really need to use addresses here, but useful to fold into external function later (pass by reference)
			for( int i = 0; i<vecs1D.size(); i++){
				(*vecs1D[i]) = row;
			}
			for( int i = 0; i<vecs2D.size(); i++){
					(*(vecs2D[i])).resize( n_states );		
				for(int n = 0; n < n_states; n++){
					(*vecs2D[i])[n] = row;
				}
			}
		
			// put 1. in starting state coefficient, leave the rest zero
			fill(coeff_set[in.init_state].begin(), coeff_set[in.init_state].end(), 1.);
			
			// Calculate initial gradients (should be zero) and Avec (also should be zero);
			calc_Gradients();
			calc_A();
			
			for(int i = 0; i < n_states; i++){
				calc_Uen_Ci(i);
			}
			
			calc_TDPES();	
			
			gauge_cond = empty_vec_nR;
			local_mask.resize(nR);
			Acheck.resize(nR);
		}
			
		///////////////////////
		void calc_Gradients(){
			vector<dcomp> blank = empty_vec_nR;
			for(int i = 0; i<n_states; i++){
				
				// call gradient and laplacian functions
				// function has the following form: 
				// calc_grad_vec( vector<dcomp> &vec, int dim, int l_boundary, int r_boundary )
				// So a lot of this stuff is vestigial but gives freedom if we want to faff with boundaries/masks/2D stuff
				
				// 23 Dec 2019
				// For the purposes of this fresh git upload, the fourier derivs here have been removed
				// The linkage to fftw remains in the makefile but the header, fftw_deriv.h, is moved to 'deprecated' folder
				// as using this approach added boundary artefacts and didn't make it into the paper. Also not convinved it was 100% correct,
				// though running it on simple functions and comparing with analytic derivs seemed to work, but spidey senses still tingled...

				// BUT! Should definitely be revisited and tested everywhere
				
				// if(!in.ft_deriv){
					grad_coeff_set[i] = calc_grad_vec( coeff_set[i], 1,  0, nR-1);
					grad_sq_coeff_set[i] =  calc_grad_sq_vec(  coeff_set[i], 1,0, nR-1);
				// }  else {
				// 	fourier_deriv_dcomp( coeff_set[i], ft_grad, ft_grad_sq, dR, 10 );
				// 	grad_coeff_set[i] = ft_grad;
				// 	grad_sq_coeff_set[i] = ft_grad_sq;
				// }
			}	
		}	

		///////////////////////
		void calc_A(){
			
			// See Eqn 16 in notes
			// Re-zeroing the vectors we build - probably (hopefully!) not needed, can check later.			
			fill(Avec.begin(), Avec.end(), 0.);
			fill(divAvec.begin(), divAvec.end(), 0.);
			fill(divAtemp.begin(), divAtemp.end(), 0.);
			
			// naughty naughty, loops within loops...
			for( int R = 0; R < nR; R++ ){
				for( int i = 0; i < n_states; i++ ){
					dcomp ci = coeff_set[i][R];
					dcomp dci = grad_coeff_set[i][R];
					dcomp d2ci = grad_sq_coeff_set[i][R];
					
					Avec[R] += -ii*(conj(ci)*dci);
					
					divAtemp[R] += real(-ii*(conj(ci)*d2ci+dci*conj(dci)));
					
					for( int j = 0; j < n_states; j++ ){
						
						dcomp cj = coeff_set[j][R];
						dcomp dcj = grad_coeff_set[j][R];
						dcomp dij = NACV_1[i][j][R];
						
						Avec[R] += -ii*( conj(ci)*cj*dij);
						divAtemp[R] += (-ii*(conj(ci)*(cj*gradNACV_1[i][j][R]+ NACV_1[i][j][R]*dcj) + cj*NACV_1[i][j][R]*conj(dci)));
					}	
				}
			}
			
			// this toggle comes from the configuration file
			// if we want to output some kind of tracker of the imaginary part as an indicator of fucked-ness, we can do that here
			if(in.A_re){
				// force keeping real part only
				// 10 JUL - fixed the real part of divA here ( was done on other branch but not back-merged)
				for( int c = 0; c < nR; c++ ){
					Avec[c] = real(Avec[c]);
					divAtemp[c] = real(divAtemp[c]);
				}
			}
			
			// Toggle to force a mask on the Avec stuff, seems to help in benchmark case
			if(in.mask_A){
				mask_vec( Avec, in.bdry_mask_edge,nR - in.bdry_mask_edge );
				mask_vec( divAtemp,in.bdry_mask_edge,nR - in.bdry_mask_edge );
			}
			divAvec = divAtemp;	
		}
		
	
		///////////////////////
		void calc_Uen_Ci( int i ){ // i is projected state!
			// See Eqn 15 in AdiabaticEF notes
			// Define some local vars and vectors for temporary computation
			vector<dcomp> UenCi, UenCi_trunc, UenCi_ctmqc;
			UenCi.resize( nR );
			UenCi_trunc.resize( nR );
			UenCi_ctmqc.resize( nR );
			dcomp T1, T2, T3;
			dcomp T1_t, T2_t, T3_t;
			dcomp T1_ctmqc, T2_ctmqc, T3_ctmqc;
			
			// Definitely need to zero these, they are the summed values!
			T1 = 0.;
			T2 = 0.;
			//~ T3 = 0.;

			T1_t = 0.;
			T2_t = 0.;
			//~ T3_t = 0.;
			
			T1_ctmqc = 0.;
			T2_ctmqc = 0.;
		
			// The projected Uen in all it's gory.
			for(int R = 0; R < nR; R++){
				T1 = (0.5*ii*divAvec[R]*in.weights[0] - 0.5*pow(Avec[R], 2.)*in.weights[1] + GCOC[R]*ii*Avec[R]*in.weights[2])*coeff_set[i][R];
				T2 =  -0.5*grad_sq_coeff_set[i][R]*(in.weights[3])- GCOC[R]*grad_coeff_set[i][R]*in.weights[4];	
				T3 = 0.;
				
			
				T1_t = (  Avec[R]*Avec[R] + ii*divAvec[R] )*coeff_set[i][R];
				T2_t = 2.*ii*Avec[R]*grad_coeff_set[i][R] - grad_sq_coeff_set[i][R];
				T3_t = 0.;
				
				T1_ctmqc = (-ii*GCOC[R] + Avec[R]);
				T2_ctmqc = T1_ctmqc*(grad_coeff_set[i][R] - Avec[R]*coeff_set[i][R]);
				T3_ctmqc = 0.;
	
				// Sum over 'j' terms given by last term in Eqn 15
				for(int j = 0; j<n_states; j++){
					T3 += 0.5*coeff_set[j][R]*NACV_2[i][j][R]*in.weights[5]	+ NACV_1[i][j][R]*(grad_coeff_set[j][R]*in.weights[6] + coeff_set[j][R]*GCOC[R]*in.weights[7] );
					T3_t +=  2.*NACV_1[i][j][R]*( ii*Avec[R]*coeff_set[j][R]- grad_coeff_set[j][R] ) - coeff_set[j][R]*NACV_2[i][j][R] ;
					T3_ctmqc += -ii*coeff_set[j][R]*NACV_1[i][j][R]*T1_ctmqc;
				}
				
				T1_vec[i][R] = T1;
				T2_vec[i][R] = T2;
				T3_vec[i][R] = T3;

				UenCi[R] = (T1 + T2 - T3)/m_N;
				UenCi_trunc[R] = 0.5*(T1_t + T2_t + T3_t)/m_N;
				UenCi_ctmqc[R] = (T2_ctmqc + T3_ctmqc)/m_N;
			}
			
			// config toggle for masking Uen
			// Here the lack of specification of mask points means it is taken from the nuclear density, i.e. comoving
			// If we want to fix the position as some kind of test, do it here
			// Could also add this to the config file for testing
			if(in.mask_uen){
				//~ mask_vec(UenCi,in.bdry_mask_edge,nR-in.bdry_mask_edge);	
				local_mask = mask_vec(UenCi);	
			}
			
			if(in.ctmqc_uen){
				UenC[i] = UenCi_ctmqc;
			} else {
				UenC[i] = UenCi;
				UenC_trunc[i] = UenCi_trunc;
			}
		}
		
	///////////////////////
	void calc_exp_Uen(){
			
			// One of the many hopefully not necessary re-zeroings.
			// Need to comb through code and one by one remove these to make sure no cross pollination
			fill( exp_Uen.begin(), exp_Uen.end(), 0. );
			
			// Is a sum over states so need that to be the inner loop
			for( int R = 0; R < nR; R++ ){
				for(int i = 0; i<n_states; i++){
					exp_Uen[R] += conj(coeff_set[i][R])*UenC[i][R];
				}
			}
			
			// here we redefine Uen to have the imaginary part of <Uen> subtracted
			// NOTE - requires un-commenting of term in calc_dt_C as well ( so that this is actually called!)
			if(in.sub_imUen){
				for( int R = 0; R<nR; R++){
					for(int i = 0; i < n_states; i++){
						UenC[i][R] = UenC[i][R] - ii*imag(exp_Uen[R])*coeff_set[i][R];
					}				
				}
			}
		}
		
		///////////////////////
		void calc_TDPES(){ 
			
			// Assumes gradients and Uen already calculated
			// 4Jul18 commented out, is this needed?
			
			// 11 Jul now using this again, see below, to ensure no nasty business coming from bits in Uen that
			// shouldn't contribute to TDPES - question is how to deal with this in backward Euler case?
			// should be okay, i.e. no need to freeze it? Especially given it should in principle be idential with
			// what is going on here, i.e. frozen stuff. Check on that branch.
			calc_exp_Uen();
			
			// default for in.dir_uen_tdpes is 1! See vars.h.
			
			// Ensure we are starting the summation from zero
			fill(TDPES.begin(), TDPES.end(), 0.);
			
			// Defined in Eqn 18 in notes
			for(int R = 0; R< nR; R++){
				for(int i = 0; i < n_states; i++){
					dcomp Ci = coeff_set[i][R];	 // shorthand
					TDPES[R] += BOPES[i][R]*pow(abs(Ci),2)  + in.dir_uen_tdpes*UenC[i][R]*conj(Ci);
					
					// IMPORTANT NOTE 11 AUGUST 17: if instead of having UenCi here we have the expectation value, get slightly different answer
					// Not entirely surprising given we need to update the expUen too
				}
				
				TDPES[R] = TDPES[R] + exp_Uen[R]*(1.- in.dir_uen_tdpes);
				if(in.re_tdpes){
					TDPES[R] = real(TDPES[R]);
				}
			}			
		}
		
		///////////////////////
		vector< vector<dcomp> > calc_dt_C( ){
			if(in.mask_coeffs){
				
				mask_vec(coeff_set[0],in.bdry_mask_edge,nR-in.bdry_mask_edge, 0.);
				mask_vec(coeff_set[1],in.bdry_mask_edge,nR-in.bdry_mask_edge, 1.);
			}

			// This should be the root call of these functions
			calc_Gradients();
			calc_A();
			//~ // Use these gradients/Avec (fixed GCOC from prev step) to calc Uen C_i
			for(int i = 0; i < n_states; i++){
				calc_Uen_Ci(i);
			}
			
			// if subtracting Im<Uen>, need to make this correction here before computing TDPES
			if(in.sub_imUen){
				calc_exp_Uen();
			}
			
			calc_TDPES();
	
			for(int R = 0; R < nR; R++){
				for(int i = 0; i < n_states; i++){
					dt_C_set[i][R] = -ii*( coeff_set[i][R]*(BOPES[i][R]-TDPES[R] ) +  UenC[i][R]);
				}
			}	
			
			// possibly inconsistent with the enforcement of gauge condition? i.e. no mask on the coeffs themselves
			if(in.mask_dt){
				for(int i = 0; i < n_states; i++){
					mask_vec(dt_C_set[i]);
				}
			}
			
			if(in.subGauge){
				checkGauge();
				for (int R = 0; R < nR; R++ ){
					for(int i = 0; i < n_states; i++){
						dt_C_set[i][R] = dt_C_set[i][R]  - gauge_cond[R]*coeff_set[i][R];  
					}
				}	
			}
			
			
		
			return dt_C_set;
		}
	
		///////////////////////
		void calc_norms(){
			// Sums over states to get electronic normsd
			// Also outputs |C_i|^2, as well as d/dt C_i and some spatial gradients
			// Norm is output as the actual norm as well as the norm deviation (1-PNC(R))
			vector<double> elec_norm, elec_norm_dev;
			vector<vector< double > > coeff_sq;

			elec_norm.resize(nR);
			elec_norm_dev.resize(nR);
			coeff_sq.resize( n_states );

			for(int i = 0; i<n_states; i++){
				coeff_sq[i].resize(nR);
				for(int R = 0; R<nR; R++){
					coeff_sq[i][R] = pow(abs(coeff_set[i][R]),2.);
					elec_norm[R] += coeff_sq[i][R];
				}
				write_vector1D(coeff_sq[i], toString(file_n) + "C"+toString(i+1)+"_sq");
				write_vector1D(grad_coeff_set[i], toString(file_n) + "grad_C"+toString(i+1)+"");
				write_vector1D(grad_sq_coeff_set[i], toString(file_n) + "grad_sq_C"+toString(i+1)+"");
				//~ write_vector1D(dt_C_set[i], toString(file_n) + "dt_C_set"+toString(i+1)+"");
				
				
			}
			for(int R = 0; R<nR; R++){
				elec_norm_dev[R] = 1.-elec_norm[R];
			}
			//~ write_vector1D(elec_norm, toString(file_n) +"phi_norm");
			//~ write_vector1D(elec_norm_dev, toString(file_n) +"phi_norm_dev");
		}
		
		///////////////////////
		void calc_B( double &dt ){
			int c = 0;
			
			for( int i = 0; i < n_states; i++ ){
					for(int R = 0; R < nR; R++){
						vec[c] = coeff_set[i][R];
						c++;
					}
				}

			// choose between CN and BE
			if(in.phi_BE){
				b = vec;
			} else if( in.phi_CN ) {
				
				vector<dcomp> H_Cvec;
				H_Cvec.resize(n_states*nR);
				
				vector< vector< dcomp > > dt_set = calc_dt_C();
				
				c = 0;
				for( int i = 0; i < n_states; i++ ){
					for(int R = 0; R < nR; R++){
						H_Cvec[c] = ii*dt_set[i][R];
						b[c] =  vec[c] - ii*dt*H_Cvec[c]/2.;	
						c++;
					}
				}
			}
		}
		///////////////////////
		// Calculate A.x = (1 + i*dt*H/2)|phi> for Crank Nicholson solvercalcaasd
		// This function is called repeatedly in the linear algebra solver
		vector< dcomp > AmatrixVec(  double dt ){
			
			// choose between CN and BE
			double factor(0);
			if(in.phi_BE){
				factor = 1.;
			} else if(in.phi_CN){
				factor = 2.;
			}
			
			LAGtol = pow(10,-15);
			vector<dcomp> H_Cvec, stepped_Cvec;
			H_Cvec.resize( n_states*nR );
			stepped_Cvec.resize( n_states*nR );
			
			int c = 0;
			for( int i = 0; i < n_states; i++ ){
				for(int R = 0; R < nR; R++){
					coeff_set[i][R] = vec[c];
					c++;
				}
			}

			vector< vector< dcomp > > dt_set = calc_dt_C();
			
			c = 0;
			for( int i = 0; i < n_states; i++ ){
				for(int R = 0; R < nR; R++){
					H_Cvec[c] = ii*dt_set[i][R];
					c++;
				}
			}
			
			for( int c = 0; c < n_states*nR; c++ ){
				stepped_Cvec[c] = vec[c] + ii*dt*H_Cvec[c]/factor;	
			}
			return stepped_Cvec;
			
		}
		
		///////////////////////
		void checkGauge(){
			gauge_cond = empty_vec_nR;
			for(int i = 0; i < n_states; i++){
				for(int R = 0; R < nR; R++){
					gauge_cond[R] += conj(coeff_set[i][R])*dt_C_set[i][R];
				}
			}	
		}
		
		///////////////////////
		void unroll(){
			for(int i = 0; i < n_states; i++){
				for(int R = 0; R < nR; R++){
					coeff_set[i][R] = vec[i*nR + R];
				}
			}
		}
		
	// end public
};	 // End of phi class

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
// Nuclear wfn state class
class chi_state {
	public:
		// variables local to chi object
		string name;
		int size, R0pos;
		
		// Store chi (vec), gradients, b (used in BICSTAB algorithm), grad chi/chi (GCOC), and H|chi>
		vector<dcomp> vec, grad_R, grad_sq_R, b, GCOC, H_chi, G2COC, dGCOC;
		
		// Store externally calculated things, Avec, divAvec, surface (TDPES)
		vector<dcomp> Avec, divAvec, surface, nuc_phase, nuc_phase_raw, grad_nuc_phase;
		
		vector<dcomp>  ft_grad, ft_grad_sq;
		
		// Declare internally used tolerances, and the norm
		double nd_edge_tol_outer, nd_edge_tol_inner, norm, R0, k0, sigsq;
		
		// Nuclear density
		vector<double> nd;
		
		vector<dcomp> ft_output;
		
		double LAGtol;
		
		///////////////////////
		void init( config in ){
			name = "chi";
			cout<<"Making chi state:"<<endl;

			size = nR;
			// These tolerances define the values of the density at which we start and stop fitting for GCOC
			// and also sometimes the mask on parts of the electronic potential
			// i.e. this gets the points between which the density spans this tolerance range and uses that data
			// to fit a 2nd order polynomial and extrapolate from there to the boundary for GCOC
			
			nd_edge_tol_outer = in.edge_tols[0];
			nd_edge_tol_inner = in.edge_tols[1];
			
			nd_edges_inner_smooth.resize(2);
			nd_edges_outer_smooth.resize(2);
			
			R0 = in.R0;
			k0 = in.k0;
			sigsq = in.sigsq;
			
			// Setup initial wfn (assumes gaussian initial state)
			vec.resize( size );
			for(int c = 0; c< size;c++){		
				double Rval = R_points[c];
				vec[c] = exp(-pow( Rval - R0,2.)/(sigsq))*exp(ii*k0*(Rval-R0));
			}
			
			// Allocate memory for various vectors	
			nd.resize( size );
			H_chi.resize( size );
			grad_R.resize( size );
			grad_sq_R.resize( size );
			GCOC.resize( size );
			G2COC.resize( size );
			dGCOC = GCOC;
			
			ft_grad_sq.resize( size );
			ft_grad.resize( size );
			
			ft_output.resize( size );
			
			// Normalize
			nd_norm();
			
			
			for(int c = 0; c<size; c++){
				vec[c] = vec[c]/sqrt(norm);
			}	
			
			// Calculate initial gradients and GCOC
			calcGradients();		
			fit_GCOC();	
			density();
			
			Avec = empty_vec_nR;
			divAvec = empty_vec_nR;
			surface = empty_vec_nR;
			calcB_CN( dt ); // Added 18 Oct 17, seems should have been in here from start?
				
			LAGtol = pow(10,-45);
		}
		
		///////////////////////
		void calcGradients(){
			for(int c = 0; c < size ; c++){
				grad_R[c] = deriv(vec, c, 1, 1, 0, nR-1);
				grad_sq_R[c] = laplacian(vec, c, 1, 1, 0, nR-1);	
			}	
			
			
			// set boundary conditions
			vec[0] = 0.;
			vec[nR-1] = 0.;
			
			grad_R[0] = 0.;
			grad_R[nR-1] = 0.;
			
			grad_sq_R[0] = 0.;
			grad_sq_R[nR-1] = 0.;
		}
		///////////////////////
		void fit_GCOC(){ // will come back and tidy this up once I work out which cases I want to test
			density();
			int max_pos = 0.;
			int R = 0;
			
			if(in.fixed_chi == 0){
				for(int c =0 ; c<size; c++){
					GCOC[c] = grad_R[c]/vec[c];	
					if(isnan(abs(GCOC[c])) || isinf(abs(GCOC[c]))){ 
						GCOC[c] = 0.;
					}
				}
			}
			max_pos = get_max_pos( nd );
			if(in.fixed_chi != 0){
				
				for(int R = 0; R<nR; R++){
					
					if(in.fixed_chi == 1){
						
						GCOC[R] = (ii*k0 -2.*dR*(R-max_pos)/sigsq) ; // comoving density GCOC shit (not really, needs free chi needs to be on)
						
					} else if( in.fixed_chi == 2 ){
						GCOC[R] = ii*k0 -2.*dR*(R-R0-0)/sigsq ; // Fixed at centre start point
					
					} else if( in.fixed_chi == 3 ){
						GCOC[R] = ii*k0 -2.*dR*(R-R0-9)/sigsq ;//  right
				
					}  else if( in.fixed_chi == 4 ){
						
						GCOC[R] = ii*k0 -2.*dR*(R-R0-9)/sigsq +100. ;// left
					
					} else if( in.fixed_chi == 5 ){
				
						GCOC[R] = 10.;
						
					} else if( in.fixed_chi == 6 ){
						
						GCOC[R] =  50.*nd[R];//*exp(pow(R-R0,2))/sigsq); // gaussian at centre val 2
					}
				}
			}	
			
			if(in.mask_GCOC){	
				//~ mask_vec(GCOC, 400, 1150);
				mask_vec(GCOC);
			}
			//, +20 for left
		}	
		//~ void calc_nuc_phase(){
			
			//~ vector<dcomp> temp_vec = vec;
		
			//~ nuc_phase_raw = empty_vec_nR;
			//~ nuc_phase = total_phase(temp_vec,nuc_phase_raw);
		
			//~ grad_nuc_phase = calc_grad_vec(nuc_phase,1, 0, nR-1);
			
		//~ }
		///////////////////////
		// calculate H|chi>
		// assumes gradients already calculated (see calcB_CN() and AmatrixVec below )
		void calc_H_chi(){
			vector<dcomp> kinval, el_pot, s2;
			kinval.resize(size);
			el_pot.resize(size);
			
			if(in.free_chi){
				Avec=empty_vec_nR;
				divAvec=empty_vec_nR;
			}
			
			for( int c = 0; c <size; c++){ 
				if(in.free_chi){
					surface[c]=BOPES[1][c];
				}
				kinval[c] = (1/(2.*m_N))*(-1.*grad_sq_R[c]);
				el_pot[c] = (1/(2.*m_N))*(- 2.*ii*(Avec[c])*grad_R[c]- ii*vec[c]*(divAvec[c]) + vec[c]*pow((Avec[c]),2.) ) + real(surface[c])*vec[c]; 
				H_chi[c] = kinval[c] + el_pot[c];
			}		
		}
		
		///////////////////////
		// calculate b = (1-i*dt*H/2)|chi> in A.x = b for asymmetric linear system
		void calcB_CN( double &dt  ){
			calcGradients();
			b.resize( size );
			calc_H_chi(  );
		
			for(int c = 0; c<size;c++){
				b[c] = vec[c] - ii*dt*H_chi[ c ]/2.;
			}
		}
		
		
		///////////////////////
		// Calculate A.x = (1 + i*dt*H/2)|chi> for Crank Nicholson solvercalcaasd
		// This function is called repeatedly in the linear algebra solver
		vector< dcomp > AmatrixVec(  double dt ){
			vector< dcomp > stepped_chi;
			stepped_chi.resize(size);

			calcGradients();
			calc_H_chi( );
			
			for(int c = 0; c < size ;c++){
				stepped_chi[c] = vec[c] + ii*dt*H_chi[c]/2.;	
			}
			
			return stepped_chi;
		}
		
		///////////////////////
		// calculate nuclear density and find edge points according to chosen tolerances
		vector<double> density(){
		
			for(int c = 0; c<size; c++){
				nd[c] = pow(abs(vec[c]),2.);
			}
			int max_pos = get_max_pos( nd );
			
			// save edge values globally (see vars.h)
			nd_edges_outer = find_edges( nd, max_pos, nd_edge_tol_outer);
			nd_edges_inner = find_edges( nd, max_pos, nd_edge_tol_inner);
	
			// hack for interpolating to get smooth mask development
			// just doing for 'outer' atm
			double p1, p2;
			double nd1, nd2;
		
			p1 = (double) nd_edges_outer[0] - 1;
			p2 = (double) nd_edges_outer[0] + 1;
			
			nd1 = nd[p1];
			nd2 = nd[p2];
	
			nd_edges_outer_smooth[0] = p1 + 2.*(nd1-nd_edge_tol_outer)/(nd1-nd2);
			p1 = (double) nd_edges_outer[1] - 1;
			p2 = (double) nd_edges_outer[1] + 1;

			nd1 = nd[p1];
			nd2 = nd[p2];
			nd_edges_outer_smooth[1] = p1 + 2.*(nd1-nd_edge_tol_outer)/(nd1-nd2);
			
			return nd;
		
		}
		///////////////////////
		// calculate norm of chi
		double nd_norm(){
			for(int c = 0; c<size; c++){
				nd[c] = pow(abs(vec[c]),2.);
			}
			
			norm = 0.;
			for(int c = 0; c<size; c++){
				norm += dR*nd[c];
			}
			
			return norm;
		}
	
};	// End of chi class


/////////////////////////////////////////////////////////////////////////////////////
void update_states( phi_state &phiS, chi_state &chiS){
	
	// calc things based on current state
	phiS.calc_Gradients();
	phiS.calc_A();
	
	chiS.calcGradients();
	chiS.density();
	chiS.fit_GCOC();
	
	// distribute info to electronic state
	phiS.GCOC = chiS.GCOC;
	
	// calc TDPES and Uen contributions
	for(int i = 0; i<n_states; i++){
		phiS.calc_Uen_Ci(i);
		
	}	
	
	phiS.calc_TDPES();

	for(int R = 0; R<nR; R++){
		phiS.stability_track[R] = phiS.divAvec[R]/m_N + 2.*phiS.Avec[R]*real(phiS.GCOC[R]);
	}
	
	// send info to chi state
	chiS.Avec = phiS.Avec;
	chiS.divAvec = phiS.divAvec;
	chiS.surface = phiS.TDPES;

	chiS.calc_H_chi();
	chiS.calcB_CN( dt );
	if(!in.phi_RK45){
		phiS.calc_B( dt );
	}
}


//~ struct rparams
  //~ {
    //~ double a;
    //~ double b;
    //~ chi_state state;
  //~ };

//~ int rosenbrock_f (const gsl_vector_complex * x, void *params, 
              //~ gsl_vector_complex * f )
//~ {
  //~ double a = ((struct rparams *) params)->a;
  //~ double b = ((struct rparams *) params)->b;

  //~ gsl_complex  x0 = gsl_vector_complex_get (x, 0);
  //~ const double x1 = gsl_vector_complex_get (x, 1);

  //~ const double y0 = a * (1. - x0);
  //~ const double y1 = b * (x1 - x0 * x0);

  //~ chi_state chiS = ((struct rparams *) params)->state;
  //~ cout<<chiS.size<<endl;
  //~ exit(0);
  //~ gsl_vector_complex_set (f, 0, y0);
  //~ gsl_vector_complex_set (f, 1, y1);

  //~ return GSL_SUCCESS;
//~ }



//~ void test(chi_state chiS ){
	
	//~ const gsl_multiroot_fsolver_type *T;
  //~ gsl_multiroot_fsolver *s;

  //~ int status;
  //~ size_t i, iter = 0;

  //~ const size_t n = 2;
  //~ struct rparams p = {1.0, 10.0, chiS};
  //~ gsl_multiroot_function f = {&rosenbrock_f, n, &p};
  //~ dcomp x_init[2] = {-10.0, -5.0};
  //~ gsl_vector_complex *x = gsl_vector_complex_alloc (n);

  //~ gsl_vector_complex_set (x, 0, x_init[0]);
  //~ gsl_vector_complex_set (x, 1, x_init[1]);
	
  //~ T = gsl_multiroot_fsolver_hybrids;
  //~ s = gsl_multiroot_fsolver_alloc (T, 2);
  //~ gsl_multiroot_fsolver_set (s, &f, x);

//~ }
/////////////////////////////////////////////////////////////////////////////////////
//END OF FILE
/////////////////////////////////////////////////////