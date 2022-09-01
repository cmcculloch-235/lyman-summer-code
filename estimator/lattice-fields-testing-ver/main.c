#include <stdlib.h>
#include <complex.h>
#include <fftw3.h>
#include <math.h>
#include <string.h>

#include "../util.h"
#include "../estimator.h"
#include "../thread_pool.h"
#include "../transformations.h"

#include "config.h"
#include "field_gen.h"

int main(int argc, char *argv[])
{




	/* *********** *
	 * Set up FFTW * 
	 * *********** */
	/* Based on `info fftw` and Jeong's thesis */

	const int N_THREADS = 8;


	/* Documentation says to enable multi-threading before calling ANY fftw
	 * functions. Presumably inclides fftw_malloc. */
	fftw_init_threads();
	fftw_plan_with_nthreads(N_THREADS);
	

	/* Number of points per side of box in real space */
	size_t X = 256;

	/* Numbers of points in k-space */
	size_t KX = X;

	/* Number of modes, relevant for normalisation of FFTW output */
	size_t N = X * X * X;
	size_t KN = KX * KX * KX;

	/* k-space and real-space buffers */
	complex double *field_ksp_buf = (complex double *) fftw_malloc(sizeof(complex double) * KN);
	complex double *field_rsp_buf = (complex double *) fftw_malloc(sizeof(complex double) * N);


	/* create plan for Fourier transforms */
	eprintf("Planning...");

	/* Check for FFTW wisdom */
	int wisdom_success = fftw_import_wisdom_from_filename("FFTW_wisdom");
	if (wisdom_success) {
		eprintf("Found FFTW wisdom...");
	}
	
	fftw_plan plan_r_to_k = fftw_plan_dft_3d(X, X, X, field_rsp_buf,
			field_ksp_buf, FFTW_FORWARD, FFTW_MEASURE);


	/* create plan for inverse Fourier transforms */
	eprintf("Planning inverse...");
	/* you would think the third argument would be KX_3 but it isn't */
	fftw_plan plan_k_to_r = fftw_plan_dft_3d(KX, KX, KX, field_ksp_buf,
			field_rsp_buf, FFTW_BACKWARD, FFTW_MEASURE);
	
	/* Always save the wisdom; it may be useful */
	fftw_export_wisdom_to_filename("FFTW_wisdom");
	eprintf("Saved FFTW wisdom...");

	eprintf("Done!\n");


	/* ************ *
	 * Set up Field * 
	 * ************ */

	/* Physical size of the box in units of Mpc/h */
	/* TODO: work out a sensible value for this */
	double L = 2500.0;
	double mode_spacing = 2.0 * M_PI / L;
	double real_spacing = L / X;
	// 1/(this * N) is also (deltak/2pi)^3
	double real_dV = real_spacing * real_spacing * real_spacing;
	
	/* Generate the linear field into the Fourier-space buffer */
	double (*spec_fn) (double) = &spec_bbks;
	eprintf("Generating spectrum...");
	/* might want to parallelise this, but managing the RNG is a little subtle */
	gen_field(field_ksp_buf, KX, mode_spacing, spec_fn);


#if (PARAM_SMOOTH)
	eprintf("Smooth...");
	smooth(field_ksp_buf, KX, mode_spacing, PARAM_SMOOTH_LEN);
#endif


	eprintf("Done!\n");

	fftw_execute(plan_k_to_r);
	normalise_k_to_r_t(field_rsp_buf, real_dV, N, N_THREADS);




	/* output generated linear field in the losfile and taufile formats */
	/* no careful error checking here... */

	/* ********************* *
	 * Clean up to be polite *
	 * ********************* */
	eprintf("Destroying plans...");
	fftw_destroy_plan(plan_r_to_k);
	fftw_destroy_plan(plan_k_to_r);

	eprintf("Freeing buffers...");
	fftw_free(field_ksp_buf);
	fftw_free(field_rsp_buf);

	eprintf("Done!\n");
	fftw_cleanup_threads();
	return 0;
}
