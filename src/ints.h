/////////////////////////////////////////////////////////////////////////////////////
// CONTAINS INTEGRALS (lazy for over R, simpson 3/8 for over r - no good reason for this)
/////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////
// Integrate over a simple 1D vector (in this case probably chi(R) )
template <typename T> 
vector<T> integrate_R( vector<T> &integrand , double h) {  
	vector< T > out;
	int size = integrand.size();
	out.resize( size );
	T val;
	out[0] = 0.;
	
	for(int Rend = 2; Rend<size; Rend++){
		out[Rend] = out[Rend-1]+integrand[Rend]*h;;
	}
	return out;
}

/////////////////////////////////////////////////////////////////////////////////////
// Integrate out the r coord i.e. to calculate A = -i<phi| grad_R | phi>
// Previously tested higher order integral i.e. Simpson, didn't alter answer significantly
// probably should add back in anyway, though too higher order can cause spurrious oscillations

template <typename T> 
vector<T> integrate_r( vector<T> &integrand, double h ) {  
	vector< T > out;
	out.resize( nR );
	
	#pragma omp parallel for shared(out)
	for(int Rpos = 0; Rpos < nR; Rpos++){
		T evens, odds, first, last;
		int c;
		
		out[ Rpos ] = 0.;
		evens = 0.;
		odds = 0.;
		
		// evens
		for(int rpos = 0; rpos < nr; rpos = rpos + 2){
			c = Rpos + nR*rpos; // convert from 1d index (c) to 2d index (r,R)
			evens += integrand[c];
		}	
		// odds
		for(int rpos = 1; rpos < nr; rpos = rpos + 2){
			c = Rpos + nR*rpos;
			odds += integrand[c];
		}	
		first = integrand[Rpos];
		last = integrand[  Rpos + nR*(nr-1) ];
		out[ Rpos ] = (h/3.)*( first + last + 4.*evens + 2.*odds );
	}
	return out;
}


/////////////////////////////////////////////////////////////////////////////////////
//END OF FILE
/////////////////////////////////////////////////////////////////////////////////////
