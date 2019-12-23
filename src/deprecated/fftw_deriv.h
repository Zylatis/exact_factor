// WIP FILE, PROBABLY WON'T WORK!
// For reference only (came back to this after a looong time and it looks like it was left in quite a state...)

#ifndef fftw_deriv // takes care of circular dependencies
#define fftw_deriv

#include <memory.h>
#include <fftw3.h>
/////////////////////////////////////////////////////////////////////////////////////


// First derivative using FFTW3
// NOTE: assumes # points is even
void ft_deriv( vector<dcomp> &vec , vector<dcomp> & first, vector<dcomp> &second ){
	
	int n = vec.size();
	double T = n;
	double df = 1./T;
	
	vector<dcomp> kleft, kright, kvec;
	
	kvec.resize(n);
	fill(kvec.begin(),kvec.end(),0.);
	dcomp ci;
	dcomp pref = ii*2.*M_PI*df;
	
	
	for(int c = 0; c<n/2; c++){
		kvec[c] = pref*((dcomp) c);
	}
	
	for(int c = n/2+1; c<n; c++){
		kvec[c] = pref*((dcomp) (c-n));
	}
	
	kvec[ n/2 ] = 0.;
	
	vector<dcomp> out_vec;
	
	out_vec.resize( n );

	fftw_complex raw[n];
	fftw_complex raw_ft[n];
	
	fftw_complex ft_grad[n];
	fftw_complex grad[n];
	
	fftw_complex ft_grad_sq[n];
	fftw_complex grad_sq[n];
	
	for(int c = 0; c < n; c++){
		raw[c][0] = real( vec[c] );
		raw[c][1] = imag( vec[c] );
	}
	
	fftw_plan forward = fftw_plan_dft_1d( n, raw, raw_ft, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute( forward );
	

	dcomp temp1, temp2;
	for(int c = 0; c < n; c++){
		temp1 = kvec[c]*(raw_ft[c][0] + ii*raw_ft[c][1]);
		temp2 = (kvec[c]*kvec[c])*(raw_ft[c][0] + ii*raw_ft[c][1]);
		
		ft_grad[c][0] = real(temp1);
		ft_grad[c][1] = imag(temp1);
		
		ft_grad_sq[c][0] = real(temp2);
		ft_grad_sq[c][1] = imag(temp2);
	}	
	ft_grad[n/2][0] = 0.;
	ft_grad[n/2][1] = 0.;
	
	
	fftw_plan inverse_grad = fftw_plan_dft_1d( n, ft_grad, grad, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_plan inverse_grad_sq = fftw_plan_dft_1d( n, ft_grad_sq, grad_sq, FFTW_BACKWARD, FFTW_ESTIMATE);
	
	
	fftw_execute( inverse_grad );
	fftw_execute( inverse_grad_sq );
	
	for(int c = 0; c < n; c++){
		first[c] = (grad[c][0] + grad[c][1]*ii)/((double)n)/dR;
		second[c] = (grad_sq[c][0] + grad_sq[c][1]*ii)/((double)n)/(dR*dR);
	}
	
	fftw_destroy_plan( forward );
	fftw_destroy_plan( inverse_grad );
	fftw_destroy_plan( inverse_grad_sq );
	fftw_cleanup();

}

void fourier_deriv( const vector<double> &vec, vector<double> &first, vector<double> &second, double dx, int n_copies, vector<double> &pad_output ){
	int n = vec.size();
	vector<double> padded_vec;
	padded_vec.reserve(n*n_copies);
	pad_output.resize(n*n_copies);
	for(int i = 0; i<n_copies; i++){
		padded_vec.insert( padded_vec.begin() + i*n, vec.begin(), vec.end() );
	}
	n = n_copies*n;
	double *in = padded_vec.data();
	
	fftw_complex out[n], ft_first[n], ft_second[n];
	double kvec[n], first_out[n], second_out[n];
	
	fftw_plan forward = fftw_plan_dft_r2c_1d( n, in, out,FFTW_ESTIMATE);
	fftw_execute( forward );
	
	
	// come back and refactor dumb shit like this
	for(int c = 0; c<n / 2; c++){
		kvec[c] = c;
	}
		
	for(int c = n/2+1; c<n; c++){
		kvec[c] = c-n;
	}	
	
	dcomp ft_temp, pref(ii*2.*M_PI/(1.*n));
	kvec[n/2] = 0.;
	
	for(int c = 0; c<n;c++){
		ft_temp = ( out[c][0] + ii*out[c][1] );
		ft_first[c][0] = real(pref*ft_temp)*kvec[c];
		ft_first[c][1] = imag(pref*ft_temp)*kvec[c];

		ft_second[c][0] = real(ft_temp*pow(pref*kvec[c],2));
		ft_second[c][1] = imag(ft_temp*pow(pref*kvec[c],2));
	
	}
	
	fftw_plan inverse_first = fftw_plan_dft_c2r_1d( n, ft_first, first_out, FFTW_ESTIMATE);
	fftw_plan inverse_second = fftw_plan_dft_c2r_1d( n, ft_second, second_out,FFTW_ESTIMATE);
	fftw_execute( inverse_first );
	fftw_execute( inverse_second );
	
	cout<<padded_vec.size()<<endl;
	for(int c = 0; c<vec.size(); c++){
		first[c] = first_out[c]/(1.*n*dx);
		second[c] = second_out[c]/(1.*n*dx*dx);
	}
	pad_output = padded_vec;
	
	
}

void fourier_deriv_dcomp( const vector<dcomp> &vec, vector<dcomp> &first, vector<dcomp> &second, double dR, int n_copies ){
	vector<double> vec_re(nR, 0.), vec_im(nR, 0), re_first(nR, 0), re_second(nR, 0), im_first(nR, 0), im_second(nR, 0);
	for( int i = 0; i<nR; i++){
		vec_re[i] = real(vec[i]);
		vec_im[i] = imag(vec[i]);
	}
	
	fourier_deriv(vec_re, re_first, re_second, dR, n_copies); 
	fourier_deriv(vec_im, im_first, im_second, dR, n_copies); 
	
	for( int i = 0; i<nR; i++){
		first[i] = re_first[i] + ii*im_first[i];
		second[i] = re_second[i] + ii*im_second[i];
	}
}

#endif
