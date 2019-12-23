/////////////////////////////////////////////////////////////////////////////////////
// CONTAINS FIRST AND SECOND DERIVATIVE FUNCTIONS

// Uses template function to allow to pass double or complex double types
// 'dir' corresponds to direction: 0 for r, 1 for R (i.e. d/dr or d/dR)
// dim is dimension of target vector i.e. phi(r,R) is dim = 2, chi(R) is dim = 1
// TODO: modularise this much more! Could be made far more compact.
/////////////////////////////////////////////////////////////////////////////////////

int calc_c_pos( int rpos, int Rpos) {
	return Rpos + nR*rpos ;
}

/////////////////////////////////////////////////////////////////////////////////////
template <typename T> 
T deriv( vector<T> &vec, int c ,int dir, int dim, int l_boundary, int r_boundary){
	
	int c_l, c_r,c_ll,c_rr, c_lll,c_rrr,c_llll,c_rrrr, xpos, ypos;

	xpos = floor(c/nR);	
	ypos = c-xpos*nR;
	int pos, pm;
	T h;
	if(dim == 1){
		pos = c;
		h = dR;
		pm = nR;
		c_l  = c - 1;
		c_ll = c - 2;
		c_lll = c - 3;
		c_llll = c - 4;
		
		c_r  = c + 1;
		c_rr = c + 2;
		c_rrr = c + 3;
		c_rrrr = c + 4;
		
	} else { // very cumbersome atm :/
		if(dir == 0){
			pos = xpos;
			h = dr;
			pm = nr;
			
			c_l  = ypos + nR*(xpos - 1);
			c_ll = ypos + nR*(xpos - 2);
			c_lll = ypos + nR*(xpos - 3);
			c_llll = ypos + nR*(xpos - 4);
			
			c_r  = ypos + nR*(xpos + 1);
			c_rr = ypos + nR*(xpos + 2);
			c_rrr = ypos + nR*(xpos + 3);
			c_rrrr = ypos + nR*( xpos + 4 );
		}
		else if(dir == 1){
			h = dR;
			pm = nR;
			pos = ypos;
			c_l  = ( ypos - 1 ) + nR*xpos;
			c_ll = ( ypos - 2 ) + nR*xpos;
			c_lll = ( ypos - 3 ) + nR*xpos;
			c_llll = ( ypos - 4 ) + nR*xpos;
			
			c_r  = ( ypos + 1 ) + nR*xpos;
			c_rr = ( ypos + 2 ) + nR*xpos;
			c_rrr = ( ypos + 3 ) + nR*xpos;
			c_rrrr = ( ypos + 4 ) + nR*xpos;
			
		}
	}
	
	if( pos  == l_boundary ){ // left boundary
		return  -(25.*vec[c] - 48.*vec[c_r] + 36.*vec[c_rr] -16.*vec[c_rrr] +3.*vec[c_rrrr] )/(12.*h);

	} else if(pos == r_boundary ) {	// right boundary
		return  (25.*vec[c] - 48.*vec[c_l] + 36.*vec[c_ll] -16.*vec[c_lll] + 3.*vec[c_llll] )/(12.*h);
	
	} else if(pos == l_boundary + 1 ) {// near left boundary
		return (-3.*vec[c_l] -10.*vec[c] + 18.*vec[c_r] -6.*vec[c_rr] + vec[c_rrr] )/(12.*h);

	} else if(pos == r_boundary - 1 ) { // near right boundary
		return -(vec[c_lll] -6.*vec[c_ll] + 18.*vec[c_l] - 10.*vec[c] - 3.*vec[c_r])/(12.*h);
		 
	} else if(pos >l_boundary && pos < r_boundary) {
		return  (-vec[ c_rr] + 8.*vec[c_r] - 8.*vec[c_l] + vec[c_ll] )/(12.*h);
		//~ return (vec[c_r] - vec[c_l])/(2.*h);
	} else if(pos <l_boundary || pos>r_boundary){
		return 0.;
	}
}


/////////////////////////////////////////////////////////////////////////////////////
template <typename T> 
T laplacian( vector<T> &vec, int c ,int dir, int dim , int l_boundary, int r_boundary){
	
	int c_l, c_r,c_ll,c_rr, c_lll,c_rrr,c_llll,c_rrrr, xpos, ypos;
	T f, fl, fll, flll, fllll, fr, frr, frrr, frrrr;
	xpos = floor(c/nR);	
	ypos = c-xpos*nR;

	int pos, pm;
	T h;
	if(dim == 1){
		pos = c;
		h = dR;
		pm = nR;
	
		c_l  = c-1;
		c_ll = c-2;
		c_lll = c-3;
		c_llll = c-4;
		
		c_r  = c+1;
		c_rr = c+2;
		c_rrr = c+3;
		c_rrrr = c+4;
		
	} else {
		
		
		if(dir == 0){
			pos = xpos;
			h = dr;
			pm = nr;
			
			c_l  = ypos + nR*(xpos - 1);
			c_ll = ypos + nR*(xpos - 2);
			c_lll = ypos + nR*(xpos - 3);
			c_llll = ypos + nR*(xpos - 4);
			
			c_r  = ypos + nR*(xpos + 1);
			c_rr = ypos + nR*(xpos + 2);
			c_rrr = ypos + nR*(xpos + 3);
			c_rrrr = ypos + nR*( xpos + 4 );
			
			
		}
		else if(dir == 1){
			h = dR;
			pm = nR;
			pos = ypos;
			c_l  = ( ypos - 1 ) + nR*xpos;
			c_ll = ( ypos - 2 ) + nR*xpos;
			c_lll = ( ypos - 3 ) + nR*xpos;
			c_llll = ( ypos - 4 ) + nR*xpos;
			
			c_r  = ( ypos + 1 ) + nR*xpos;
			c_rr = ( ypos + 2 ) + nR*xpos;
			c_rrr = ( ypos + 3 ) + nR*xpos;
			c_rrrr = ( ypos + 4 ) + nR*xpos;
			
		}
	}
	
	if( pos  == l_boundary ){ // left boundary
		return (35.*vec[c]-104.*vec[ c_r]+114.*vec[ c_rr]-56.*vec[ c_rrr]+11.*vec[ c_rrrr ])/(12.*pow(h,2.)); 
	} 
	 else if(pos == r_boundary) { // right boundary 
		return  (35.*vec[c]-104.*vec[ c_l]+114.*vec[ c_ll]-56.*vec[ c_lll]+11.*vec[ c_llll ])/(12.*pow(h,2.));	
	}
	 else if(pos ==  l_boundary + 1 ) { // near left boundary
		return -(vec[c_rrr] - 4.*vec[c_rr] -6.*vec[c_r] + 20.*vec[c] -11.*vec[c_l])/(12.*pow(h,2.));
	} 
	else if(pos == r_boundary - 1) { // near right boundary
		return -(vec[c_lll] - 4.*vec[c_ll] -6.*vec[c_l] + 20.*vec[c] -11.*vec[c_r])/(12.*pow(h,2.));
		
	} else if(pos >l_boundary && pos < r_boundary) {
		return ( -vec[c_rr] + 16.*vec[c_r] - 30.*vec[c] +16.*vec[c_l] -vec[c_ll])/(12.*pow(h,2.));
	
	} else if(pos <l_boundary || pos>r_boundary){
		return 0.;
	}

}


/////////////////////////////////////////////////////////////////////////////////////
//END OF FILE
/////////////////////////////////////////////////////////////////////////////////////
