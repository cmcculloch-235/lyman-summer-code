#include <stdlib.h>
#include <complex.h>
#include <fftw3.h>
#include <math.h>
#include <string.h>
#include <stdint.h>

#include "util.h"
#include "estimator.h"
#include "thread_pool.h"
#include "config.h"
#include "file_io.h"

#include "transformations.h"



int main(int argc, char *argv[])
{

	/* *************************** *
	 * Read headers from data file *
	 * *************************** */
	

	eprintf("Reading LoS header...");

	FILE *losfile_fp = fopen(LOSFILE, "rb");
	FILE *taufile_fp = fopen(LOSFILE, "rb");
	if (!losfile_fp) {
		eprintf("Error opening file %s\n", LOSFILE);
		return 1;
	}
	if (!taufile_fp) {
		eprintf("Error opening file %s\n", TAUFILE);
		return 1;
	}
	
	double ztime, omegam, omegal, omegab, h100, box_length, XH;
	uint32_t nbins, nlos;

	{
		int file_read_err = 0;
		file_read_err = read_losfile_header(losfile_fp, &ztime, &omegam, &omegal, &omegab,
				&h100, &box_length, &XH, &nbins, &nlos);
		if (file_read_err) {
			eprintf("Error reading header in LoS file (path: %s)\n", LOSFILE);
			return 1;
		}

	}

	eprintf("Done!\n");

	/* Number of points per side of box in real space */
	size_t X = nbins;

	/* Number of modes, relevant for normalisation of FFTW output */
	/* Note that this is guaranteed to hold because the read stage required
	 * a cubic grid */
	size_t N = X * X * X;

	double real_spacing = box_length / X;
	double real_dV = real_spacing * real_spacing * real_spacing;
	double mode_spacing = 2 * M_PI / real_spacing;



	/* *********** *
	 * Set up FFTW * 
	 * *********** */
	/* Based on `info fftw` and Jeong's thesis */

	const int N_THREADS = 8;

	/* Documentation says to enable multi-threading before calling ANY fftw
	 * functions. Presumably includes fftw_malloc. */
	fftw_init_threads();
	fftw_plan_with_nthreads(N_THREADS);
	
	/* k-space and real-space FFTW buffers */
	complex double *fft_ksp_buf = (complex double *) fftw_malloc(sizeof(complex double) * N);
	complex double *fft_rsp_buf = (complex double *) fftw_malloc(sizeof(complex double) * N);

	/* create plan for Fourier transforms */
	eprintf("Planning FFTs...");

	/* Check for FFTW wisdom */
	int wisdom_success = fftw_import_wisdom_from_filename("FFTW_wisdom");
	if (wisdom_success) {
		eprintf("Found FFTW wisdom...");
	}
	
	fftw_plan plan_r_to_k = fftw_plan_dft_3d(X, X, X, fft_rsp_buf,
			fft_ksp_buf, FFTW_FORWARD, FFTW_MEASURE);

	/* create plan for inverse Fourier transforms */
	eprintf("Planning inverse...");
	fftw_plan plan_k_to_r = fftw_plan_dft_3d(X, X, X, fft_ksp_buf,
			fft_rsp_buf, FFTW_BACKWARD, FFTW_MEASURE);
	
	/* Always save the wisdom; it may be useful */
	fftw_export_wisdom_to_filename("FFTW_wisdom");
	eprintf("Saved FFTW wisdom...");

	eprintf("Done!\n");


	/* ************** *
	 * Read in fields * 
	 * ************** */
	eprintf("Load data...");



	complex double *delta_H, *delta_H1, *delta_DM, *delta_matter, *tau_field, *delta_tau;

	{
		int file_read_err = 0;
		eprintf("LoS file...");
		file_read_err =  read_losfile_data(losfile_fp, N, nlos, &delta_H, &delta_H1,
				&delta_DM, &delta_matter);
		if (file_read_err) {
			eprintf("Error reading fields from LoS file (path: %s)\n", LOSFILE);
			return 1;
		}

		eprintf("tau file...");
		file_read_err =  read_taufile_data(taufile_fp, N, &tau_field, &delta_tau);
		if (file_read_err) {
			eprintf("Error reading fields from tau file (path: %s)\n", TAUFILE);
			return 1;
		}
	}
	

	fclose(losfile_fp);
	fclose(taufile_fp);

	eprintf("Done!\n");


	/* ************************ *
	 * Time to process the data *
	 * ************************ */

	/* NOTE: if there's going to be configuration-space field processing,
	 * do it here. */

#define N_FIELDS 5
	complex double *field_list[N_FIELDS] = {delta_H, delta_H1, delta_DM,
		delta_matter, delta_tau};
	const char *field_names[N_FIELDS] = {"delta_H", "delta_H1", "delta_DM",
		"delta_matter", "delta_tau"};

	/* need Fourier-space fields */
	eprintf("FFT fields...");
	for (size_t i = 0; i < N_FIELDS; ++i) {
		eprintf("%s...", field_names[i]);
		memcpy(fft_rsp_buf, field_list[i], N * sizeof(complex double));
		fftw_execute(plan_r_to_k);
		normalise_r_to_k_t(fft_ksp_buf, real_dV, N, N_THREADS);
		memcpy(fft_ksp_buf, field_list[i], N * sizeof(complex double));
	}
	eprintf("Done!\n");

	/* do all sorts of wonderful auto- and cross-correlations */
	eprintf("Cross-correlations...");


	double *xcorr_output_buffer = calloc(xcorr_bin_count, sizeof(double));
	double *xcorr_k_buffer = calloc(xcorr_bin_count, sizeof(double));
	size_t *xcorr_count_buffer = calloc(xcorr_bin_count, sizeof(size_t));
	if (!xcorr_output_buffer || !xcorr_k_buffer || !xcorr_count_buffer) {
		eprintf("Error allocating memory for cross-correlation output.\n");
		return 1;
	}

	for (int i = 0; i < N_FIELDS; ++i) {
		for (int j = i; j < N_FIELDS; ++j) {
			eprintf("%s-%s...", field_names[i], field_names[j]);
			/* k range and bin count specified in config.h */
			correlator(field_list[i], field_list[j], X, mode_spacing, xcorr_k_buffer,
					xcorr_output_buffer, xcorr_count_buffer, xcorr_k_min,
					xcorr_k_max, xcorr_bin_count);

			write_xcorr_data(field_names[i], field_names[j], xcorr_k_buffer,
					xcorr_output_buffer, xcorr_count_buffer, xcorr_bin_count);
		}
		eprintf("\n");
	}
	free(xcorr_output_buffer);
	free(xcorr_k_buffer);
	free(xcorr_count_buffer);

	eprintf("Done!\n");

/* this is mainly here for reference */

	/* Extract the power spectrum at first and second order and print it to 
	 * standard output */
	//double *k_buffer = calloc(pspec_bins, sizeof(double));
	//double *power_buffer = calloc(pspec_bins, sizeof(double));
	//size_t *n_buffer = calloc(pspec_bins, sizeof(size_t));

	/* will want to replace these with a call to correlate(...) or whatever
	 * it's called. can modify that to call thread_map and make very fast
	 * estimates */

	//power_spectrum(linear_ksp, X, mode_spacing, k_buffer, power_buffer_lin,
			// 5 mode_spaceing  up to 1.3X mode_spacing / 2
	//		 K_MIN, K_MAX, pspec_bins);

	//free(k_buffer);
	//free(power_buffer);
	//free(n_buffer);





	/* ********************* *
	 * Clean up to be polite *
	 * ********************* */

	free(delta_H);
	free(delta_H1);
	free(delta_DM);
	free(delta_matter);
	free(tau_field);
	free(delta_tau);


	eprintf("Destroying plans...");
	fftw_destroy_plan(plan_r_to_k);
	fftw_destroy_plan(plan_k_to_r);

	eprintf("Freeing buffers...");
	fftw_free(fft_ksp_buf);
	fftw_free(fft_rsp_buf);

	eprintf("Done!\n");
	fftw_cleanup_threads();
	return 0;
}
