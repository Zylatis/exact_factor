/////////////////////////////////////////////////////////////////////////////////////
// DEFINE RK45 ROUTINE
/////////////////////////////////////////////////////////////////////////////////////

vector<vector< dcomp > > add_RK_vecs( vector<vector< dcomp > > init_coefs, vector<vector< dcomp > > k_i , double h){
	vector<vector< dcomp > > out;
	int ns = init_coefs.size();
	int l = init_coefs[0].size();
	out.resize(ns);
	for(int i = 0; i<ns; i++){
		out[i].resize(l);

		for(int R = 0; R<l; R++){
			out[i][R] = init_coefs[i][R] + h*k_i[i][R];
		}
	}
	return out;
	
}

/////////////////////////////////////////////////////////////////////////////////////
// RK45 solver routine for phi
// Recalculates Avec at each multi-step point!
void RK_solve( phi_state &state, double dt ){
	phi_state k1s, k2s, k3s, k4s;
	vector< vector<dcomp> > k1, k2, k3, k4, final;
	vector< vector<dcomp> > k2_in, k3_in, k4_in;
	int s = state.size;
	
	vector<vector<dcomp> > y0 = state.coeff_set;
	
	k1s = state;
	k2s = state;
	k3s = state;
	k4s = state;
	
	k1s.coeff_set = y0;
	k1 = k1s.calc_dt_C();
	k2_in = add_RK_vecs( y0, k1, dt/2. );
	
	k2s.coeff_set = k2_in;
	k2 = k2s.calc_dt_C();
	k3_in = add_RK_vecs( y0, k2, dt/2. );
	
	k3s.coeff_set = k3_in;
	k3 = k3s.calc_dt_C();
	k4_in = add_RK_vecs( y0, k3, dt );
	
	k4s.coeff_set = k4_in;
	k4 = k4s.calc_dt_C();
	
	final = y0;
	
	vector<vector< dcomp > > dt_vec;
	dt_vec.resize(n_states);
	
	for(int i = 0; i<n_states; i++){
		vector<dcomp> temp;
		temp.resize(nR);
		fill(temp.begin(), temp.end(), 0.);
		//~ #pragma omp parallel for
		for(int R = 0; R<nR; R++){
			temp[R] = (dt/6.)*(k1[i][R] + 2.*k2[i][R] + 2.*k3[i][R] + k4[i][R]);
		}
	
		if(in.mask_dt){
			mask_vec(temp);
		}

		for(int R = 0; R<nR; R++){
			final[i][R] += temp[R];
		}
		
	}
	
	state.coeff_set = final;
}

/////////////////////////////////////////////////////////////////////////////////////

// Archived attempt at (not really but sorta) predictor-corrector code
// Just does some mixing of forward and backward solution
//~ void pred_cor( phi_state &state, double dt ){
	//~ vector<vector< dcomp > > dt_vecP, dt_vecC, local_coeff_set, y0_state;
	//~ phi_state PCstate = state;
	//~ y0_state = state.coeff_set;
	
	//~ // initial step
	//~ dt_vecP = RK_solve( PCstate, dt );
	
	//~ for(int j = 0; j<10; j++){
		//~ for(int i = 0; i < n_states; i++){
			//~ for(int R = 0; R < nR; R++){
				//~ PCstate.coeff_set[i][R] = y0_state[i][R] + dt_vecP[i][R];
			//~ }
		//~ }
		
		//~ // get d/dt at new point
		//~ dt_vecC = RK_solve( PCstate, dt );
		
		//~ // do mixing, i.e. project d/dt backwards and mix
		//~ for(int i = 0; i < n_states; i++){
			//~ for(int R = 0; R < nR; R++){
				//~ dt_vecP[i][R] = (1./2.)*(dt_vecP[i][R] + dt_vecC[i][R]);
			//~ }
		//~ }
	//~ }
	
	//~ state.coeff_set = PCstate.coeff_set;
//~ }

/////////////////////////////////////////////////////////////////////////////////////
//END OF FILE
/////////////////////////////////////////////////////////////////////////////////////
