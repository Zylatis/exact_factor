/////////////////////////////////////////////////////////////////////////////////////
// CONTAINS LINEAR ALGEBRA SOLVER FOR CRANK NICHOLSON STEP AS WELL AS
// CRANK NICHOLSON STEPPER FUNCTION FOR EXACT FACTOR CASE

// See BiCGStabAlg.pdf for algorithm used.
/////////////////////////////////////////////////////////////////////////////////////

#ifndef linsol // takes care of circular dependencies
#define linsol

/////////////////////////////////////////////////////////////////////////////////////
// A BUNCH OF FUNCTIONS FOR LINEAR ALGEBRA STUFF 
/////////////////////////////////////////////////////////////////////////////////////

// Compute dot product of two complex vectors
// Here dot product represents inner product on vector space, hence the conjugate
// Concerning that the OpenMP part here doesn't seem to lead to more core usage, need better profiling!
dcomp  dot_product( vector< dcomp  > &v1, vector<dcomp> &v2, int &l){
	dcomp sum = 0.;
	vector<dcomp> par_sum;
	par_sum.resize(n_threads);
	#pragma omp parallel for
	for(int i = 0; i<l;i++){
		int id = omp_get_thread_num();
		par_sum[ id ] += conj( v1[i] )*v2[i];
	}
	
	//do serial summation of parallel bits
	for(int c = 0;c<n_threads; c++){
		sum += par_sum[c];
	}
	
	return sum;
}
	
/////////////////////////////////////////////////////////////////////////////////////
// Subtract two 1D vectors
vector< dcomp > subtract_vector( vector< dcomp > &v1, vector< dcomp > &v2, int &l){
	vector< dcomp > temp;
	temp.resize(l);
	#pragma omp parallel for
	for(int i = 0; i<l;i++){
		temp[i] = v1[i] - v2[i];
	}
	return temp;
}

/////////////////////////////////////////////////////////////////////////////////////
// Add two 1D vectors
vector< dcomp > add_vector( vector< dcomp > &v1, vector< dcomp > &v2,  double a, double b, int &l){
	vector< dcomp > temp;
	temp.resize(l);
	#pragma omp parallel for
	for(int i = 0; i<l;i++){
		temp[i] = a*v1[i] + b*v2[i];
	}
	return temp;
}

/////////////////////////////////////////////////////////////////////////////////////
// Add two 1D vectors
vector< dcomp > sum_vec( vector< dcomp > &v1, vector< dcomp > &v2 ){
	vector< dcomp > temp;
	int s = v1.size();
	temp.resize(s);
	#pragma omp parallel for
	for(int i = 0; i<s;i++){
		temp[i] = v1[i] + v2[i];
	}
	return temp;
}

/////////////////////////////////////////////////////////////////////////////////////
// conjugate vector
vector< dcomp > conj_vec( vector< dcomp > &v1 ){
	vector<dcomp> temp;
	int s = v1.size();
	temp.resize( s );
	for(int c = 0; c<s; c++){
		temp[c] = conj(v1[c]);
	} 
	return temp;
}

/////////////////////////////////////////////////////////////////////////////////////
// Multiply element-wise two vectors
vector< dcomp > product_of_vectors( vector< vector< dcomp> > &vec_set){
	int n_vecs = vec_set.size();
	int s = vec_set[0].size();
	
	vector<dcomp> out_vec;
	out_vec.resize( s );
	
	#pragma omp parallel for
	for(int c = 0; c<s; c++){
		out_vec[c] = 1.;
		for(int n = 0; n<n_vecs; n++){
			out_vec[c] = out_vec[c]*vec_set[n][c];
		}
	}
	return out_vec;
}
/////////////////////////////////////////////////////////////////////////////////////
// Calculate the norm of a 1D vector
double vec_norm( vector< dcomp  >  &vec ){
	double sum = 0;
	for(int c = 0; c<l; c++){
		sum += (double) pow(abs(vec[c]),2.);
	}
	return sum*dr*dR;
}

/////////////////////////////////////////////////////////////////////////////////////
// Calc p_n for linear solver
vector<dcomp> calc_pn( vector<dcomp> &rprev, vector<dcomp> &pprev, vector<dcomp> &vprev, dcomp &beta, dcomp &omegaprev, int &l){
	vector< dcomp > temp;
	temp.resize(l);
	#pragma omp parallel for
	for(int i=0; i<l; i++){
		temp[i] = rprev[i] + beta*(pprev[i] - omegaprev*vprev[i]);
	}
	return temp;
}

/////////////////////////////////////////////////////////////////////////////////////
// Calculate s_n for linear solver
vector<dcomp> calc_sn( vector<dcomp> &rprev, vector<dcomp> &vn, dcomp &alpha, int &l){
	vector< dcomp > temp;
	temp.resize(l);
	#pragma omp parallel for
	for(int i=0; i<l; i++){
		temp[i] = rprev[i] - alpha*vn[i];
	}
	return temp;
}

/////////////////////////////////////////////////////////////////////////////////////
// Calculate x_n for linear solver
vector<dcomp> calc_xn( vector<dcomp> &xprev, vector<dcomp> &pn, vector<dcomp> &sn, dcomp &alpha, dcomp &omegan,  const int &l){
	vector< dcomp > temp;
	temp.resize(l);
	#pragma omp parallel for
	for(int i=0; i<l; i++){
		temp[i] = xprev[i] + alpha*pn[i] + omegan*sn[i];
	}
	return temp;
}

/////////////////////////////////////////////////////////////////////////////////////
// Calculate r_n for linear solver
vector<dcomp> calc_rn( vector<dcomp> &sn, vector<dcomp> &tn,  dcomp &omegan,  int &l){
	vector< dcomp > temp;
	temp.resize(l);
	#pragma omp parallel for
	for(int i=0; i<l; i++){
		temp[i] = sn[i] - omegan*tn[i];
	}
	return temp;
}

/////////////////////////////////////////////////////////////////////////////////////
// Run algebraic solver to compute Crank Nicholson step
// Specifically, solve (1+i*dt*H/2)psi_n+1 = (1-i*dt*H/2)psi_n for psi_n+1, given H and psi_n
template <class T> 
vector<dcomp> algebraic_solve( T &state,  double &dt, double &t, int &len){
	
	dcomp rho0, rhon, rhoprev;
	dcomp alpha;
	dcomp omega0, omegan, omegaprev;
	dcomp beta;
	
	vector<dcomp> r0, rn, rprev, xF, rhat0;
	vector<dcomp> v0,vn, vprev;
	vector<dcomp> p0, pn, pprev;
	vector<dcomp> sn;
	vector<dcomp> tn;
	vector<dcomp> x0, xn, xprev;
	vector<dcomp> diff1;
	
	xF.resize( len );
	x0.resize( len );

	p0.resize( len );
	v0.resize( len );
	rhat0.resize( len );

	x0 = state.b; //  initial guess - alg pdf suggests 0, Lionel suggests 'b'
	fill(p0.begin(), p0.end(), 0.);
	fill(v0.begin(), v0.end(), 0.);
	
	rho0 = 1.;
	alpha = 1.;
	omega0 = 1.;
	omegaprev = omega0;
	rhoprev = rho0;

	xprev = x0;
	pprev = p0;
	vprev = v0;
	
	state.vec = x0;
	vector<dcomp> state_b = state.b;  // must be already calculated
	vector<dcomp> state_AmatrixVec = state.AmatrixVec( dt );
	r0 = subtract_vector(state_b, state_AmatrixVec, len);	

	rn = r0;
	rprev = r0;
	rhat0 = r0;
	vector<dcomp> tempv;
	vector<dcomp> state_save = state.vec;
	
	double err;
	for(int i = 0; i < 	LS_n + 1; i++ ){
		
		rhon = dot_product( r0 , rprev, len);
		
		beta = (rhon/rhoprev)*(alpha/omegaprev);
		pn = calc_pn(rprev, pprev, vprev, beta, omegaprev, len);	
		
		state.vec = pn;
		vn = state.AmatrixVec(dt);
		
	
		alpha = rhon/dot_product(r0,vn, len);
		sn = calc_sn(rprev, vn, alpha, len);

		state.vec = sn;
		tn = state.AmatrixVec(dt);
		
		omegan = dot_product(tn,sn, len)/dot_product(tn,tn, len);
		xn = calc_xn(xprev, pn, sn, alpha, omegan, len);
		
		rn = calc_rn(sn,tn,omegan, len);
		err = abs( dot_product(rn, rn, len)/dot_product(xn , xn ,len) );
	
	
		// if relative error less than given tolerance, keep this solution
		if(err < state.LAGtol){
			// ~ cout<<"LA solver: "<<state.name<<" converged to "<<err<<" in "<<i<<" steps."<<	endl;
			return xn;
		}
		// catch NaNs
		if(isnan(err)){
			cout<<"Isnan err" <<state.name<<" i = "<<i<<endl;
			cout<<"rhon"<<"\t"<<rhon<<endl;	
			cout<<"beta"<<"\t"<<beta<<endl;	
			cout<<"alpha \t"<<alpha<<endl;
			cout<<"omegan \t"<<omegan;
			write_vector1D(pn, "failed_pn");
			write_vector1D(vn, "failed_vn");
			write_vector1D(tn, "failed_tn");
			write_vector1D(xn, "failed_xn");
			write_vector1D(rn, "failed_rn");
			write_vector1D(state.b, "failed_b");
			write_vector1D(r0, "failed_r0");
		
			dump_state = true;
			exit(0);
			return xn;
		}

		xprev = xn;
		rprev = rn;
		omegaprev = omegan;
		rhoprev = rhon;
		pprev = pn;
		vprev = vn;
    }
    //If we get to this point, pretty fucked!
	cout<<"Wfn "+state.name + " failed to converge in " + toString(LS_n) + " iterations!"<<endl;
	write_vector1D(xn, "unconverged_xn");
	cout<<err<<"\t"<<state.LAGtol<<endl;
	exit(0);
	return xn;
}


#endif
/////////////////////////////////////////////////////////////////////////////////////
//END OF FILE
/////////////////////////////////////////////////////////////////////////////////////
