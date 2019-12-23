/////////////////////////////////////////////////////////////////////////////////////
// FUNCTIONS TO CALCULATE THE NON-ADIABATIC COUPLING VECTORS
// 
// Not used if 'have phi' set to 0 (will read in previously computed NACVs from folder)
// It is this code that requires the keeping of the 2d derivatives in derivs.h (phi_R(r,t) is 2d fn)
/////////////////////////////////////////////////////////////////////////////////////

///////////////////////
void NACVs( int i, int j, int size, vector<vector< dcomp > > &flat_BO_vecs ){
	// calculate first and second order NACV
	vector<dcomp> ket_1, ket_2, integrand_1, integrand_2, out_1, out_2;
	integrand_1.resize( size );
	integrand_2.resize( size );
	
	
	// Here gradient acts on j-state so  d_ij = <phi_i \grad phi_j >
	ket_1 =  calc_grad_vec( flat_BO_vecs[j], 2, 0, nR-1 ); // taking deriv of 2D state (r,R in R direction)
	ket_2 =  calc_grad_sq_vec( flat_BO_vecs[j], 2, 0, nR-1 ); // taking deriv of 2D state (r,R in R direction)
	
	for(int c = 0; c<size;c++){
		integrand_1[c] = conj(flat_BO_vecs[i][c])*ket_1[c]; // first order
		integrand_2[c] = conj(flat_BO_vecs[i][c])*ket_2[c]; // second order
	}
	
	// take expectation values over electronic coordinate using intergrate_r
	out_1 = integrate_r(integrand_1, dr );
	out_2 = integrate_r(integrand_2, dr );
	
	// See Eqn 13 in notes: d_ij = <phi_i \grad phi_j >
	// This is important when considering indicies in Eqn 15: the 2nd index of the NACV should be the one being summed over!
	NACV_1[i][j] = out_1;
	NACV_2[i][j] = out_2;		
}
	
///////////////////////
void calc_NACVs(){
	// read in adiabatic basis vectors 
	// (computed from Ali's code for ShinMetiu, for Tully NACVs are computed in mathematica file)
	vector< vector< vector<double> >  > BO_states;
	BO_states.resize( n_states );
	
	// NOTE MAKE CHECK OF SIZE OF BOPES VS INIT WFN OTHERWISE SHITTY FUCKING SEGFAULTS
	for( int c = 0; c<n_states;c++){
		BO_states[c] = FileRead(model_loc+"initwfn/"+toString(BO_basis[c])+"_reboewf"+submodel+".dat");
	}
	
	// get info about grid (# points, spacing, etc)
	int size = BO_states[0].size();
	r_points = get_col(BO_states[0],0);
	get_unique(r_points);
	nr = r_points.size();
	dr = abs(r_points[1]-r_points[0]);
	
	// Read y-values of BO wfns into 1d lists (following same format as BOwfns)
	vector<vector< dcomp> >  flat_BO_vecs;
	flat_BO_vecs.resize(n_states);
	for(int i = 0; i<n_states; i++){
		flat_BO_vecs[i].resize(nr*nR);
		for(int c = 0; c<nr*nR; c++){
			flat_BO_vecs[i][c] = BO_states[i][c][2];
		}
	}
	
	// Calculate first and 2nd order NACVs (not currently imposing any conditions here)
	for(int i = 0; i<n_states; i++){
		for(int j = 0; j<n_states; j++){
			NACVs(i,j, size, flat_BO_vecs);
		}
		fill(NACV_1[i][i].begin(), NACV_1[i][i].end(), 0.);
	}	
	
	
	// Write to file
	cout<<"Writing NACV to output:"<<endl;
	for(int i = 0; i<n_states; i++){
		for(int j = 0; j<n_states; j++){
			write_vector1D(NACV_1[i][j], "NACV1-"+toString( i + 1 )+toString( j + 1 )+submodel, model_loc+"NACV/" );
			write_vector1D(NACV_2[i][j], "NACV2-"+toString( i + 1 )+toString( j + 1 )+submodel, model_loc+"NACV/" );
		}
	}
}
	
/////////////////////////////////////////////////////////////////////////////////////
//END OF FILE
/////////////////////////////////////////////////////////////////////////////////////

