#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multimin.h>	

/////////////////////////////////////////////////////////////////////
// least squares cost function
double compute_SSR(const gsl_vector *v,  void * data){

	int fit_par_n = 2;
	double SSR = 0.;
	
	double * local_dat = (double * ) data;
	int size = local_dat[0];
	vector<double> pars;
	
	pars.resize(fit_par_n);
	pars[0] = gsl_vector_get(v, 0);
	pars[1] = gsl_vector_get(v, 1);
	//~ pars[2] = gsl_vector_get(v, 2);
	//~ pars[3] = gsl_vector_get(v, 3);
	//~ pars[3] = 0.;

	for(int c = 0; c < size; c++){
		SSR = SSR + pow(poly_fn( c , pars ) - local_dat[ c + 1 ], 2.) ;
	}	

	//~ cout<<SSR<<endl;
	return SSR;
	
}

/////////////////////////////////////////////////////////////////////
// fit function currently global - ideal to have in class but meh.
vector<double> do_fit( double * data ){
	// Using non-rando Nelder-Mead without analytic Jacobian
	const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
	gsl_multimin_fminimizer *s = NULL;
	gsl_vector *ss, *x;
	gsl_multimin_function minex_func;
	
	int nn = 2;
	size_t iter = 0;
	int status;
	double size;
	/* Starting point */
	x = gsl_vector_alloc(nn);
	gsl_vector_set (x, 0, 1.);
	gsl_vector_set (x, 1, 1.);
	//~ gsl_vector_set (x, 2, 1.);
	//~ gsl_vector_set (x, 3, pow(10.,-8));
	ss = gsl_vector_alloc(nn);

	gsl_vector_set(ss, 0, 1.);
	gsl_vector_set(ss, 1, 1.);
	//~ gsl_vector_set(ss, 2, 1.);
	//~ gsl_vector_set(ss, 3, 1.);
	//~ /* Initialize method and iterate */
	minex_func.n = nn;
	minex_func.f = compute_SSR;
	minex_func.params = data;
	
	s = gsl_multimin_fminimizer_alloc (T, nn);

	gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

	do
	{
	  iter++;
	  status = gsl_multimin_fminimizer_iterate(s);
	  
	  if (status) 
		break;

	  size = gsl_multimin_fminimizer_size (s);
	  status = gsl_multimin_test_size (size, 1e-7);

	 
	} while (status == GSL_CONTINUE && iter < pow(10.,5));
	
	
	vector<double> pars;
	pars.resize( nn );
	for(int c = 0; c<nn;c++){
		pars[c] = gsl_vector_get (s->x, c);
	}

	gsl_vector_free(x);
	gsl_vector_free(ss);
	gsl_multimin_fminimizer_free (s);

	return pars;
}
