#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>


#include "thread_pool.h"
#include "util.h"


static double index_to_k_helper(size_t l, size_t KX, double mode_spacing);

int eprintf(const char *restrict format, ...)
{
	
	va_list ap;
	va_start(ap, format);
	int ret =  vfprintf(stderr, format, ap);
	va_end(ap);
	return ret;
}


double index_to_k(size_t l, size_t m, size_t n, size_t KX, double mode_spacing)
{
	double kl = index_to_k_helper(l, KX, mode_spacing);
	double km = index_to_k_helper(m, KX, mode_spacing);
	double kn = index_to_k_helper(n, KX, mode_spacing);


	return sqrt(kl * kl + km * km + kn * kn);
}

static double index_to_k_helper(size_t l, size_t KX, double mode_spacing)
{
	/* Accounts for the ordering of the modes in FFTW output */
	int mode_number = (int) l;
	if (l > KX / 2) 
		mode_number = ((int) l) - (int) KX;
	
	return mode_spacing * mode_number;
}

void index_to_vec_k(size_t l, size_t m, size_t n, size_t KX,
		double mode_spacing, double *out)
{
	out[0] = index_to_k_helper(l, KX, mode_spacing);
	out[1] = index_to_k_helper(m, KX, mode_spacing);
	out[2] = index_to_k_helper(n, KX, mode_spacing);
}




/* Converts from [l][m][n] notation to [l + ()m + ()()n] */
size_t field_index(size_t l, size_t m, size_t n, size_t KX)
{
	return n + KX * m + KX * KX * l;
}
size_t field_rsp_index(size_t l, size_t m, size_t n, size_t X)
{
	return n + X * m + X * X * l;
}

void index_to_coords(size_t index, size_t KX, size_t *out)
{
	out[2] = index % KX;
	index -= index % KX;
	index /= KX;

	out[1] = index % KX;
	index -= index % KX;
	index /= KX;


	out[0] = index % KX;
	index -= index % KX;
	index /= KX;

}



complex double discrete_ksp_gradient(size_t l, size_t m, size_t n, size_t component,
		size_t KX, double mode_spacing, double real_spacing)
{
	double k_i;
	switch (component) {
		case 0:
			k_i =  index_to_k_helper(l, KX, mode_spacing);
			break;
		case 1:
			k_i =  index_to_k_helper(m, KX, mode_spacing);
			break;
		case 2:
			k_i =  index_to_k_helper(n, KX, mode_spacing);
			break;
		default:
			// this should never happen...
			return 0;
	}
	//return 1.0 / real_spacing * (cexp(I * k_i * real_spacing) - 1);
	return I * k_i;
}

complex double discrete_ksp_laplacian(size_t l, size_t m, size_t n, size_t KX,
		double mode_spacing, double real_spacing)
{
	complex double laplacian = 0.0;
	for (size_t i = 0; i < 3; ++i) {
		complex double cpt = discrete_ksp_gradient(l, m, n, i, KX, mode_spacing, real_spacing);
		laplacian += cpt * cpt;
	}
	return laplacian;
}



void normalise_k_to_r(void *in, void *out, size_t index, void *arg)
{
	struct normalisation_arg *a = (struct normalisation_arg *) arg;
	*(complex double *)out = *(complex double *)in / (a->real_dV * a->N);

}
void normalise_r_to_k(void *in, void *out, size_t index, void *arg)
{
	struct normalisation_arg *a = (struct normalisation_arg *) arg;
	*(complex double *)out = *(complex double *)in * a->real_dV;

}


/* conveniently wrap the threading interface for normalisation calls */

void normalise_k_to_r_t(complex double *rsp, double real_dV, size_t N, size_t N_THREADS)
{
	struct normalisation_arg nn_arg = {.real_dV = real_dV, .N = N};
	thread_map(&normalise_k_to_r, (void *) &nn_arg, (void *) rsp, sizeof(complex double),
			(void *) rsp, sizeof(complex double), N, N_THREADS);
}

void normalise_r_to_k_t(complex double *ksp, double real_dV, size_t N, size_t N_THREADS)
{

	struct normalisation_arg nn_arg = {.real_dV = real_dV, .N = N};
	thread_map(&normalise_r_to_k, (void *) &nn_arg, (void *) ksp, sizeof(complex double),
			(void *) ksp, sizeof(complex double), N, N_THREADS);

}
