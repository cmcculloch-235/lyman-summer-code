#include <stdlib.h>
#include <complex.h>
#include <fftw3.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>

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
	//double L = 2500.0;
	//80Mpc box
	double L = 80;
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
	transform_smooth(field_ksp_buf, KX, mode_spacing, PARAM_SMOOTH_LEN);
#endif


	eprintf("Done!\n");

	fftw_execute(plan_k_to_r);
	normalise_k_to_r_t(field_rsp_buf, real_dV, N, N_THREADS);




	/* output generated linear field in the losfile and taufile formats */
	/* no careful error checking here... */

	{

		/* LoS file format:
		 * ****************
		 * Header: all doubles except the last two, which are uint32_t
		 * ztime omegam omegal omegab h box_length XH nbins nlos
		 * lengths in ckpc/h
		 * ****************
		 * Body:
		 * ilos: LoS axis direction: nlos * uint32_t
		 * xlos: LoS x-position: nlos * double
		 * ylos
		 * zlos
		 * posaxis: distances along LoS for each bin: nbins * double
		 * velaxis: same but in different units? (recession velocity?)
		 * rhoker_H: gas overdensity: nlos * nbins * double
		 * rhoker_H1: neutral hydrogen
		 * tempker
		 * velker
		 * rhoker_DM: CDM overdensity
		 * rhoker_total: total matter (baryons+CDM) overdensity
		 * */
		eprintf("Generate output...");

		FILE *losfile_fp = fopen("losfile_test.dat", "wb");
		/* write header. nonsense except where it describes field dimensions */
		/* ztime */
		eprintf("LoS header...");
		double out_d = 0.0;
		fwrite(&out_d, sizeof(double), 1, losfile_fp);
		/* omegam */
		fwrite(&out_d, sizeof(double), 1, losfile_fp);
		/* omegal */
		fwrite(&out_d, sizeof(double), 1, losfile_fp);
		/* omegab */
		fwrite(&out_d, sizeof(double), 1, losfile_fp);
		/* h */
		out_d = 1.0;
		fwrite(&out_d, sizeof(double), 1, losfile_fp);

		/* box_length in ckpc/h */
		out_d = L * 1000;
		fwrite(&out_d, sizeof(double), 1, losfile_fp);

		/* XH */
		out_d = 0.0;
		fwrite(&out_d, sizeof(double), 1, losfile_fp);

		/* nbins */
		uint32_t nbins = X;
		fwrite(&nbins, sizeof(uint32_t), 1, losfile_fp);

		/* nlos */
		uint32_t nlos = X * X;
		fwrite(&nlos, sizeof(uint32_t), 1, losfile_fp);

		/* data: nonsense in unused fields; all the ones we care about
		 * are just the generated Gaussian random fields. */


		/* add one to the field so it gives overdensity and not contrast */
		/* also we should be writing doubles, not complex doubles! */
		double *field_out_buf = calloc(N, sizeof(double));
		for (size_t i = 0; i < N; ++i) {
			field_out_buf[i] = 1.0 + field_rsp_buf[i];
		}


		/* ilos */
		eprintf("LoS body...");
		fwrite(field_out_buf, sizeof(uint32_t), nlos, losfile_fp);
		/* x,y,z los */
		fwrite(field_out_buf, sizeof(double), nlos, losfile_fp);
		fwrite(field_out_buf, sizeof(double), nlos, losfile_fp);
		fwrite(field_out_buf, sizeof(double), nlos, losfile_fp);
		/* pos,vel axis */
		fwrite(field_out_buf, sizeof(double), nbins, losfile_fp);
		fwrite(field_out_buf, sizeof(double), nbins, losfile_fp);
		


		/* rhoker_H */
		fwrite(field_out_buf, sizeof(double), N, losfile_fp);
		/* rhoker_H1 */
		fwrite(field_out_buf, sizeof(double), N, losfile_fp);
		/* tempker, velker */
		fwrite(field_out_buf, sizeof(double), N, losfile_fp);
		fwrite(field_out_buf, sizeof(double), N, losfile_fp);
		/* rhoker_CDM */
		fwrite(field_out_buf, sizeof(double), N, losfile_fp);
		/* rhoker_total */
		fwrite(field_out_buf, sizeof(double), N, losfile_fp);

		fclose(losfile_fp);

		eprintf("tau file...");

		FILE *taufile_fp = fopen("taufile_test.dat", "wb");
		/* overdensity field also works as total tau field, having nonzero mean */
		fwrite(field_out_buf, sizeof(double), N, taufile_fp);
		fclose(taufile_fp);
		free(field_out_buf);
		eprintf("Done!\n");


	}

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
