#include <math.h>
#include <complex.h>
#include <fftw3.h>
/* gsl has a good interface for Gaussian random numbers */
#include <gsl/gsl_randist.h>
/* to seed the RNG */
#include <time.h>

#include "../util.h"
#include "field_gen.h"


double spec_flat(double k)
{
	/* just a uniform power spectrum for testing purposes */
	return 8.0;
}
double spec_gaussian(double k)
{
	return  exp(- k * k * pow(1.0e2, 2)) * 1e8;
}

double spec_delta(double k)
{
	if (k >= 1e-3)
		return 0.0;
	return 1.0e9;

}

double spec_linear(double k)
{
	return k;
}


static double bbks_f(double x);
static const double EPSILON = 1e-8;

double spec_bbks(double k)
{
	/* units of h/Mpc in k, (Mpc/h)^3 in A */
	const double A = 5.1e4;
	const double k_eq = 0.01;
	const double n = 0.967;

	return A * pow(k / k_eq, n) * pow(bbks_f(k / k_eq), 2);

}


static double bbks_f(double x)
{
	double Q = log(1 + 0.171 * x) / (0.171 * x + EPSILON);
	double R = 1 + 0.274 * x + pow(1.18 * x, 2) + pow(0.399 * x, 3) + pow(0.49 * x, 4);
	return Q * pow(R, -0.25);
}


void gen_field(complex double *field_buffer, size_t KX, double mode_spacing,
		double (*power_spectrum) (double))
{
	/* Volume element in k-space is just 1/(number of points in real space * vol of rsp)*/
	/* This is also (delta k)^3 / ( (2pi)^3 (n points)) */
	double dV = pow(mode_spacing / (2 * M_PI ), 3);

	/* Set up the RNG. Type is luxury random numbers, which have good
	 * decorrelation.*/
	gsl_rng *rng = gsl_rng_alloc(gsl_rng_ranlxd2);
	/* seed the RNG. Wants a 32-bit float*/
	gsl_rng_set(rng, time(NULL));
	
	/* Step through the modes, get the power spectrum, and generate the field */

	for (size_t l = 0; l < KX; ++l) {
		for (size_t m = 0; m < KX; ++m) {
			for (size_t n = 0; n < KX/2 + 1; ++n) {

				double k = index_to_k(l, m, n, KX, mode_spacing);
				/* evaluate power spectrum at k, then draw real and imaginary
				 * parts from corresponding Gaussian distributions */
				double ps = power_spectrum(k);
				double sigma = sqrt(ps / (2.0 * dV));
				//sigma = sqrt(ps /2.0);

				double real_part = gsl_ran_gaussian(rng, sigma);
				double imag_part = gsl_ran_gaussian(rng, sigma);



				/* we are stepping through in row-major order */
				field_buffer[field_index(l, m, n, KX)] = real_part + I * imag_part;
				
				// enforce real-ness of field
				field_buffer[field_index((KX - l) % KX, (KX - m) % KX, (KX - n) % KX, KX)]
					= real_part - I * imag_part;
				// special cases where -k and k modes are equal
				if ((l == 0 || (l == KX/2 && ! (KX % 2))) &&
						(m == 0 || (m == KX/2 && ! (KX % 2))) &&
						(n == 0 || (n == KX/2 && ! (KX % 2))) ){

					field_buffer[field_index(l, m, n, KX)] = real_part;
				}

				// extra-special case: zero mean density contrast
				if (!(l || m || n)) {
					field_buffer[0] = 0;
				}

				/* ********* *
				 * DEBUGGING *
				 * ********* **/
				#if PARAM_FG_DEBUG
				if (!(l ||m ||n)) {
					eprintf("Debugging field mode");
				}
				/*
				if (( l == KX/4 ) ) {
					real_part = KX * KX * KX * KX;
				}*/
				complex double field_val = cexp(-k * k * 1e4 - KX/2 * I * l) * KX * KX * KX;


				/* we are stepping through in row-major order */
				field_buffer[field_index(l, m, n, KX)] = field_val;
				
				// enforce real-ness of field
				//eprintf("%ld %ld %ld   " , l, m, n);
				field_buffer[field_index((KX - l) % KX, (KX - m) % KX, (KX - n) % KX, KX)]
					= conj(field_val);
				if ((l == 0 || (l == KX/2 && ! (KX % 2))) &&
						(m == 0 || (m == KX/2 && ! (KX % 2))) &&
						(n == 0 || (n == KX/2 && ! (KX % 2))) ){
					field_buffer[field_index(l, m, n, KX)] = creal(field_val);
				}
				#endif

			}
		}
	}

	/* Free the RNG to avoid a memory leak! */
	gsl_rng_free(rng);



}


