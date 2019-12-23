////////////////////////////////////////////////////////////////////////////////////
// DEFINITION OF MASK AND SMOOTHING FUNCTIONS
/////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////
// Single sided mask based on tanh function - similar to logistic function
// but has nicer derivative properties (though probably not important here, thought 
// it might have been previously)
double mask_LL( int x, double c, double w, double s ){
	// s = 1, step goes up, s = - 1, step goes down (you can't explain that?!)
	
	double x1 = c - w/2;
	double x2 = c + w/2;
	
	double Lval, Rval;
	Lval = 0.5*(1. - s);
	Rval = 0.5*(1. - s*(-1.));
	
	
	
	if( x <= x1 ){
		return Lval;
 	} else if (x >= x2){
		return Rval;
	} else {
		double k1 = ((double) w)/abs(x-x1);
		double k2 = ((double) w)/abs(x-x2);
		double mv = 0.5*(1. - s*tanh( k1-k2 ));
		return mv;
	}		
}


/////////////////////////////////////////////////////////////////////////////////////
// n-point moving average for (possibly) complex vector
// np must be odd (duh! - is centred by default)
vector<dcomp> SMA( vector<dcomp> &vec, int np ){
	if (!(np % 2)) {cout<<"np is even, exiting:"<<endl; exit(0); }
	vector<dcomp> smoothed;
	int s = vec.size();
	smoothed.resize( s );
	double npp = (double) (np);
	//LHS edge
	for(int i = 0; i < np-2; i++){
		dcomp val = 0.;
		for(int j = 0; j < np-1; j++){
			val += vec[i+j];
		}
		smoothed[i] = vec[i];//val/(npp-1);
	}
	
	//RHS edge
	for(int i = s - np  + 2; i < s; i++){
		dcomp val = 0.;
		for(int j = 0; j < np-1; j++){
			val += vec[i-j];
		}
		smoothed[i] = vec[i];//val/(npp-1);
	}

	//centre
	for(int i = np - 2; i < s - np + 2; i++){
		dcomp val = 0.;
		for(int j = -(np-1)/2; j <= (np-1)/2; j++){
			val += vec[i+j];
		}
		smoothed[i] = val/npp;
	}
	return smoothed;
}

/////////////////////////////////////////////////////////////////////////////////////
// n-point moving average for (possibly) complex vector
//~ // np must be odd (duh! - is centred by default)
//~ vector<dcomp> SMA2( vector<dcomp> &vec, int np ){
	//~ if (!(np % 2)) {cout<<"np is even, exiting:"<<endl; exit(0); }
	//~ vector<dcomp> smoothed;
	//~ int s = vec.size();
	//~ smoothed.resize( s );
	//~ double npp = (double) (np);
	//~ vector<int> default_points;
	//~ default_points.resize( np );
	
	//~ for(int i = 0; i<np; i++){
		//~ default_points[i] = i-(np-1)/2;
	//~ }
	
	//~ for( int i = 0; i < s; i++){
		//~ int ns = (np-1)/2;
		//~ vector< int > points;
		//~ points.resize( np );
		//~ int y;
		
		//~ if( i - ns < 0){ //LH edge
			//~ y = i - ns;
		//~ } else if( i + ns > s-1 ){ //RH edge
			//~ y = i - ns - s - 1;
		//~ } else { // rest
			//~ y = 0;
		//~ }
		
		//~ for(int j = 0; j < np; j++){
			//~ int cp = j- ns + y;
			//~ dcomp sum 
		//~ }
		
		
		
	//~ }
	//~ return smoothed;
//~ }
/////////////////////////////////////////////////////////////////////////////////////
// Applies given mask to a vector and saves mask values in 'global_mask' for plotting
vector<double> mask_vec( vector<dcomp> &vec, double l =  nd_edges_outer_smooth[0], double r =  nd_edges_outer_smooth[1], double target = 0. ){
	double lw, rw;
	lw = in.mask_width; // mask width
	rw = lw;
	
	for(int R = 0; R<vec.size(); R++){
		double mv1 = mask_LL( R, l, lw, 1 );
		double mv2 = mask_LL( R, r, rw, -1 );
		
		double mv1p = mask_LL( R, l, lw, -1 );
		double mv2p = mask_LL( R, r, rw, 1 );
		
		global_mask[R] = mv1*mv2;
		vec[R] = vec[R]* mv1*mv2 + target*(mv1p+mv2p);
	}
	return global_mask;
}




/////////////////////////////////////////////////////////////////////////////////////
vector<dcomp> total_phase( vector<dcomp> vec, vector<dcomp> &raw_vec){
	int s = vec.size();	
	vector<dcomp> out;
	out.resize( s );
	fill(out.begin(), out.end(), 0.);
	out[0] = arg(vec[0]);
	//~ mask_vec(vec);
	int n = 0;
	double shift = 0.;
	for(int R = 1; R < s; R++){
		raw_vec[R] = arg(vec[R]);
		double phase = arg( vec[R] );
		double phase_diff = arg(vec[R]) -arg(vec[R-1]);
		
		if(phase_diff>M_PI){
			shift -= 2.*M_PI;
		} else if(phase_diff<-M_PI){
			shift += 2.*M_PI;
		}
		out[R] = phase + shift;
		
	}
	out[s-1] = arg( vec[s-2] );

	return out;
}

/////////////////////////////////////////////////////////////////////////////////////
// polynomial flat step mask

void step_mask(vector<dcomp> &vec, int r2_i, int w, int dir ){
	//~ write_vector1D(vec, "TEST_raw");
	int r1_i, r3_i;
	if(dir == 1){
		r1_i = r2_i - w;
		r3_i = r2_i + w;
	} else if(dir == -1){
		r1_i = r2_i + w;
		r3_i = r2_i - w;
	}
	
	dcomp df1 = deriv( vec, r1_i, 1, 0, 0, nR-1)*dR;
	dcomp f1 = vec[r1_i];
	dcomp f2 = vec[r2_i];
	
	dcomp r1 = (double) r1_i;
	dcomp r2 = (double) r2_i;
	dcomp r3 = (double) r3_i;
	
	dcomp a = (f2*pow(r1,2.)*(r1 - 3.*r3) + pow(r3,2.)*(f1*(3.*r1 - r3) + df1*r1*(-r1 + r3)))/pow(r1 - r3,3.);
	dcomp b = (r3*(6.*f1*r1 - 6.*f2*r1 + df1*(-2.*pow(r1,2.) + r1*r3 + pow(r3,2.))))/pow(-r1 + r3,3.);
	dcomp c = (3.*f1*(r1 + r3) - 3.*f2*(r1 + r3) - df1*(r1 - r3)*(r1 + 2.*r3))/pow(r1 - r3,3.);
	dcomp d = (-2.*f1 + 2.*f2 + df1*(r1 - r3))/pow(r1 - r3,3.);
	
	if(dir == 1){
		for(int R = r1_i; R<=r3_i; R++){
			dcomp Rv = (double) R;

			vec[R] = a + b*Rv + c*pow(Rv,2.) + d*pow(Rv,3.);
		}
		for(int R = r3_i; R<=nR; R++){
			vec[R] = f2;
		}
	} else if (dir == -1){
		for(int R = r1_i; R>=r3_i; R--){
			dcomp Rv = (double) R;
			vec[R] = a + b*Rv + c*pow(Rv,2.) + d*pow(Rv,3.);
		}
		for(int R = r3_i; R>=0; R--){
			vec[R] = f2;
		}
	}
	
}

/////////////////////////////////////////////////////////////////////////////////////
vector<dcomp> SG_filter_vec( vector<dcomp> &vec){
	int s = vec.size();
	int side_l = (SGwidth - 1)/2;
	
	vector<dcomp> padded;
	int padded_s = s + 2*side_l;
	padded.reserve( padded_s );
	
	vector<dcomp> LH, RH;
	LH.resize(side_l);
	RH.resize(side_l);
	
	
	dcomp LH_edge_val = vec[0];
	dcomp RH_edge_val = vec[s-1];

	for(int i = 0; i< side_l; i++){
		dcomp LHdiff, RHdiff, LHval,RHval;
		LHval = vec[side_l - 1 -i+1 ];
		RHval = vec[ s-1-i-1];
		
		LHdiff = LHval - LH_edge_val;
		RHdiff = -(RHval - RH_edge_val);

		LH[i] = LH_edge_val - LHdiff;
		RH[i] = RH_edge_val + RHdiff;

	}
	
	padded.insert( padded.end() , LH.begin(), LH.end() );
	padded.insert( padded.end() , vec.begin(), vec.end() );
	padded.insert( padded.end() , RH.begin(), RH.end() );
	
	vector<dcomp> out;
	out.resize(s);
	fill(out.begin(), out.end(), 0.+0.*ii);
	for(int R = 0; R<s;R++){
		vector<dcomp> subVec(padded.begin() + R, padded.begin() + R + SGwidth);
		for(int i = 0; i<SGwidth;i++){
			subVec[i] = subVec[i]*SG_coeffs[i];
		}
		out[R] = accumulate(subVec.begin(), subVec.end(), 0.+ 0.*ii);
			
	}

	vec = out;
	return padded;
}

/////////////////////////////////////////////////////////////////////////////////////
//END OF FILE
/////////////////////////////////////////////////////////////////////////////////////
