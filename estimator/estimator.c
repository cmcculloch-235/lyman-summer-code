#include <stdlib.h>
#include <complex.h>
#include <fftw3.h>
#include <math.h>

#include "config.h"
#include "util.h"
#include "estimator.h"


/* TODO: write a more general field correlator function */
/* It will be a fairly simple copy-paste of this one */

void power_spectrum(complex double *field, size_t KX, double mode_spacing,
		double *k_buffer, double *bin_buffer, size_t *n_buffer, double k_min,
		double k_max, size_t n_bins)
{
	/* This is so we can find the mean k and power in each bin */
	//int *n_buffer = calloc(n_bins, sizeof(int));
	//actually, n_buffer is now provided to us
	
	/* For buffer assignment */
	double l_k_min = log(k_min);
	double l_k_max = log(k_max);
	double l_range = l_k_max - l_k_min;


	/* mode volume in k-space/density of states: if we get PS by doing
	 * an integral, equivalent in DFT is multiplying by V where it's nonzero*/
	double dV = pow(mode_spacing / (2 * M_PI), 3);

	/* Step through the field components in row-major order */
	for (size_t l = 0; l < KX; ++l) {
		for (size_t m = 0; m < KX; ++m){
			for (size_t n = 0; n < KX; ++n) {
				/* Find out what k is */
				double k = index_to_k(l, m, n, KX, mode_spacing);

				/* Find out which bin this k goes in. They are uniform
				 * in log-space. */
				double l_k = log(k);
				int bin = (int) trunc(n_bins * (l_k - l_k_min) / l_range);
				/* if k is out of the desired range, don't record it */

				if (bin < 0 || bin >= (int) n_bins)
				{
					continue;
				}

				/* Record the power, wavevector, and where it ended up */
				n_buffer[bin] += 1;

				complex double f_val = field[field_index(l, m, n, KX)];

				//complex double f_val_2 = field[field_index((KX - l) % KX, (KX - m) % KX, (KX -n) % KX, KX)];
				//complex double power = f_val * f_val_2;
				//eprintf("%f+%fi  ", creal(power), cimag(power));
				double power = f_val * conj(f_val) * dV;
				
				bin_buffer[bin] += power;

				k_buffer[bin] += k;
				/*
				if (k < 5.0) {
					eprintf("%f %f\n", k, f_imag);
				}
				*/

			}
		}
	}

	/* Step through the bins and take means of k and power */
	for (size_t i = 0; i < n_bins; ++i) {
		bin_buffer[i] /= n_buffer[i];
		k_buffer[i] /= n_buffer[i];
	}

	//free(n_buffer);
}


static double window_deconv_helper(double x);

void correlator(complex double *field_1, complex double *field_2, size_t KX,
		double mode_spacing, double h, double *k_buffer, double *bin_buffer,
		size_t *n_buffer, double k_min, double k_max, size_t n_bins)
{
	/* This is so we can find the mean k and power in each bin */
	//int *n_buffer = calloc(n_bins, sizeof(int));
	//actually, n_buffer is now provided to us
	
	/* For buffer assignment */
	double l_k_min = log(k_min);
	double l_k_max = log(k_max);
	double l_range = l_k_max - l_k_min;


	/* mode volume in k-space/density of states: if we get PS by doing
	 * an integral, equivalent in DFT is multiplying by V where it's nonzero*/
	double dV = pow(mode_spacing / (2 * M_PI), 3);

	/* Step through the field components in row-major order */
	for (size_t l = 0; l < KX; ++l) {
		for (size_t m = 0; m < KX; ++m){
			for (size_t n = 0; n < KX/2 + 1; ++n) {
				/* Find out what k is */
				double k = index_to_k(l, m, n, KX, mode_spacing);

				/* Find out which bin this k goes in. They are uniform
				 * in log-space. */
				double l_k = log(k);
//				int bin = (int) round(n_bins * (l_k - l_k_min) / l_range);
				int bin = (int) round(n_bins * (k - k_min) / (k_max - k_min));
				/* if k is out of the desired range, don't record it */

				if (bin < 0 || bin >= (int) n_bins)
				{
					double kmax = index_to_k(KX/2, KX/2, KX/2, KX, mode_spacing);
					if (k >= kmax)
						eprintf("skip: k = %e (max %e)\n", k, kmax );
					continue;
				}


				complex double f1_val = field_1[field_index(l, m, n, KX)];
				complex double f2_val = field_2[field_index(l, m, n, KX)];

				double corr = f1_val * conj(f2_val) * dV;


				/* debug */
				/*
				double q, r;
				if ((l >= KX/2 + 1  && m >= KX/2 + 1)) {
					signed int L = ((signed int) l) - KX;
					signed int M = ((signed int) m) - KX;
					q = sqrt(L * L + M * M + n*n);
				}
				else
					continue;
				corr = q;
				*/

				/* apply estimator corrections */
				/* window deconvolution */
				double ks[3] = {0.0};
				index_to_vec_k(l, m, n, KX, mode_spacing, ks);
				/* Nyquist wavenumber */
				double k_Ny = mode_spacing * KX / 2.0;

				double window_deconv = 1.0;
				window_deconv *= window_deconv_helper(ks[0] * M_PI_2 / k_Ny);
				window_deconv *= window_deconv_helper(ks[1] * M_PI_2 / k_Ny);
				window_deconv *= window_deconv_helper(ks[2] * M_PI_2 / k_Ny);



				/* shot noise: 
				* just the inverse of the  total density of particles
				* particle number defined in config.h */
				/* In fact, dV = (mode_spacing/2pi)^3 = 1/L^3 = 1/vol */
				complex double s_noise = 1.0 / (dV * N_PARTICLES);


				
				/* Record the correlator, wavevector, and where it ended up */

				bin_buffer[bin] += corr / (window_deconv * window_deconv) - s_noise;
				/* debug */
				//bin_buffer[bin] += corr;
				n_buffer[bin] += 1;
				k_buffer[bin] += k;

			}
		}
	}

	/* Step through the bins and take means of k and correlator */
	for (size_t i = 0; i < n_bins; ++i) {
		bin_buffer[i] /= n_buffer[i];
		k_buffer[i] /= n_buffer[i];

		/* NOTE: lengths are in cMpc/h, so right now, correlator is in
		 * (cMpc/h)^3 (since we do an integral d^3 x and square to get pspec * delta
		 * and integrate again to get delta to go away).
		 * multiply by h^3 to get power in cMpc^6, which doesn't get extra
		 * time-dependent scaling and is in the coordinates in which perturbations
		 * grow in a simple manner. */
#ifndef USE_MPC_H
		/* convert to cMpc^3 */
		bin_buffer[i] *= pow(h, 3);
#endif

	//	bin_buffer[i] = 1.0;
		//k_buffer[i] = 0.0;
	}

}


static double window_deconv_helper(double x)
{
	if (x > 0.1)
		return sin(x)/x;
	return -pow(x,6.)/5040.0 + pow(x,4.)/120.0 - pow(x,2.)/6.0 + 1.0;
}
