#include <stdlib.h>
#include <complex.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "config.h"
#include "util.h"
#include "file_io.h"


/* this reads parameters from losfile's header into the pointed variables.
 * return value is nonzero if something went wrong. */
int read_losfile_header(FILE *losfile_fp, double *ztime, double *omegam, double *omegal,
		double *omegab, double *h100, double *box_length, double *XH,
		uint32_t *nbins, uint32_t *nlos)
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
	 * */


	/* for error detection */
	size_t n_items_read = 0;
	size_t n_items_expected = 0;

	/* expect doubles to be 64-bit */
	/* this should not be a problem on 64-bit systems, but standards don't
	 * guarantee it */
	if (sizeof(double) != 8) {
		eprintf("sizeof(double) = %zd; expected 8. Run on a different system.\n",
				sizeof(double));
		return 1;
	}

	n_items_expected = 7 + 2;

	n_items_read += fread(ztime, sizeof(double), 1, losfile_fp);
	n_items_read += fread(omegam, sizeof(double), 1, losfile_fp);
	n_items_read += fread(omegal, sizeof(double), 1, losfile_fp);
	n_items_read += fread(omegab, sizeof(double), 1, losfile_fp);
	n_items_read += fread(h100, sizeof(double), 1, losfile_fp);
	n_items_read += fread(box_length, sizeof(double), 1, losfile_fp);
	/* convert to cMpc/h */
	*box_length /= 1000.0;

#ifndef USE_MPC_H
	/* convert to cMpc */
	box_length *= h100;
#endif

	n_items_read += fread(XH, sizeof(double), 1, losfile_fp);
	
	n_items_read += fread(nbins, sizeof(uint32_t), 1, losfile_fp);
	n_items_read += fread(nlos, sizeof(uint32_t), 1, losfile_fp);

	if (n_items_read != n_items_expected) {
		eprintf("LoS file ended before header was read.\n");
		return 1;
	}
	
	/* simplifying assumption: equal lattice resolution in each direction. */
	if (*nbins * *nbins != *nlos) {
		eprintf("Non-cubic lattice in LoS file.\n");
		return 1;
	}

	return 0;
}

static int read_contrast(complex double **field, double *tmp_buffer,
		size_t N, FILE *losfile_fp);


int read_losfile_data(FILE *losfile_fp, size_t N, size_t X, complex double **delta_H,
		complex double **delta_H1, complex double **delta_DM, complex double **delta_matter,
		double omegam, double omegab)
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
	 * rhoker_H: gas overdensity: nlos * nbins * double. this is also the total baryon overdensity
	 * rhoker_H1: neutral hydrogen
	 * tempker
	 * velker
	 * rhoker_DM: CDM overdensity
	 * */

	int fseek_err = 0;

	/* skip past header */
	fseek_err = fseek(losfile_fp, 7 * sizeof(double) + 2 * sizeof(uint32_t), SEEK_SET);

	/* since we assume a sensible grid, ilos, xlos, ylos, zlos, posaxis, velaxis
	 * are not interesting to us. */
	int nlos = X * X;
	int nbins = X;
	fseek_err = fseek(losfile_fp, nlos * sizeof(uint32_t) + 3 * nlos * sizeof(double)
			+ 2 * nbins * sizeof(double), SEEK_CUR);

	if (fseek_err) {
		eprintf("Error skipping irrelevant data in LoS file\n");
		return 1;
	}

	/* fields are real doubles, but need to convert them to complex doubles
	 * for FFTW */
	double *tmp_buffer = calloc(N, sizeof(double));
	if (!tmp_buffer) {
		eprintf("Error allocating memory to read LoS file.\n");
		return 1;
	}

	/* this is a little perverse, but it will work. */
	/* it could have just been a function... but it would complicate error
	 * checking fractionally. */

	int read_contrast_err = 0;
	eprintf("delta_H...");
	read_contrast_err |= read_contrast(delta_H, tmp_buffer, N, losfile_fp);
	eprintf("delta_H1...");
	read_contrast_err |= read_contrast(delta_H1, tmp_buffer, N, losfile_fp);
	if (read_contrast_err) {
		return 1;
	}

	fseek_err = fseek(losfile_fp, 2 * N * sizeof(double),
			SEEK_CUR);
	if (fseek_err) {
		eprintf("Error skipping irrelevant data in LoS file\n");
		return 1;
	}


	eprintf("delta_DM...");
	read_contrast_err |= read_contrast(delta_DM, tmp_buffer, N, losfile_fp);

	if (read_contrast_err) {
		return 1;
	}


	eprintf("delta_matter...");
	/* no longer read from the file, since the format changed. calculate it instead */
	(*delta_matter) = calloc(N, sizeof(complex double));
	if (!(*delta_matter)) {
		eprintf("Failed to allocate memory.\n");
		return 1;
	}

	for (size_t i = 0; i < N; ++i) {
		(*delta_matter)[i] = (omegab * (*delta_H)[i] + (omegam - omegab) * (*delta_DM)[i]) / omegam;
	}




	/* dark matter normalisation issue */
	/* output is n_particles / particles per snapfile, but want
	 * n_particles / average particles = n_particles / (particles per snapfile 
	 * * n_snaps) * n_cells */
	/* and is particles per snapfile even constant? */
	/* IN FACT, all the fields seem to be coming in with mean 1, after all,
	 * so maybe it's fine */
	/* solve this by just contrasting fields in read_contrast */
	/*
	for (size_t i = 0; i < N; ++i) {
		(*delta_DM)[i] = (((*delta_DM)[i] + 1.0) * N / N_SNAPS) - 1.0;
	}
	*/
	
	/*
	double ev_DM = 0.0;
	double ev_H = 0.0;
	for (size_t i = 0; i < N; ++i) {
		ev_DM += (*delta_matter)[i] + 1.0;
		ev_H += (*delta_H)[i] + 1.0;
	}
	ev_DM /= N;
	ev_H /= N;
	eprintf("\nDM: %f | H: %f\n", ev_DM, ev_H);
	*/


	free(tmp_buffer);
	return 0;

}

static void contrast(complex double **field, double *tmp_buffer, size_t N);
static int read_contrast(complex double **field, double *tmp_buffer,
		size_t N, FILE *losfile_fp) {
	bzero(tmp_buffer, N * sizeof(double));
	size_t n_items_read = 0;
	n_items_read = fread(tmp_buffer, sizeof(double), N, losfile_fp);
	if (n_items_read != N) {
		eprintf("LoS file ended while reading data.\n");
		return 1;
	}
	*field = calloc(N, sizeof(complex double));
	if (!*field) {
		eprintf("Error allocating memory while reading LoS file.\n");
		return 1;
	}
	contrast(field, tmp_buffer, N);
	/* if we had very fine grids, this would be a good place to use
	 * thread_map; as it stands, it's not. */
	/*
	for (size_t i = 0; i < N; ++i) {
		(*field)[i] = (complex double) (tmp_buffer[i] - 1.0);
//		(*field)[i] = (complex double) (tmp_buffer[i] );
	}
	*/

	return 0;
}

static void contrast(complex double **field, double *tmp_buffer, size_t N)
{

	double mean = 0.0;
	for (size_t i = 0; i < N; ++i) {
		mean += tmp_buffer[i];
	}
	mean /= N;
	eprintf("\n %f \n", mean);
	/* some kind of sensible cut-off */
	if (fabs(mean) < 1e-8) {
		for (size_t i = 0; i < N; ++i) {
			(*field)[i] = tmp_buffer[i];
		}
		return;
	}

	for (size_t i = 0; i < N; ++i) {
		(*field)[i] = tmp_buffer[i] / mean - 1.0;
	}

}


int read_taufile_data(FILE *taufile_fp, size_t N, complex double **tau_field,
		complex double **delta_tau, complex double **delta_flux)
{
	/* and now read data from tau file */
	/* this just contains the tau LoSs; no headers nor additional data */
	
	/* fields are real doubles, but need to convert them to complex doubles
	 * for FFTW */
	double *tmp_buffer = calloc(N, sizeof(double));
	if (!tmp_buffer) {
		eprintf("Error allocating memory to read tau file.\n");
		return 1;
	}

	size_t n_items_read = 0;
	n_items_read = fread(tmp_buffer, sizeof(double), N, taufile_fp);
	if (n_items_read != N) {
		eprintf("LoS file ended while reading tau_field\n");
		return 1;
	}


	*tau_field = calloc(N, sizeof(complex double));
	if (!tau_field) {
		eprintf("Error allocating memory for tau_field.\n");
		return 1;
	}
	*delta_tau = calloc(N, sizeof(complex double));
	if (!delta_tau) {
		eprintf("Error allocating memory for delta_tau.\n");
		return 1;
	}
	*delta_flux = calloc(N, sizeof(complex double));
	if (!delta_flux) {
		eprintf("Error allocating memory for delta_flux.\n");
		return 1;
	}

	/* while we're doing this, work on normalisation */
	double tau_mean = 0.0;
	for (size_t i = 0; i < N; ++i) {
		(*tau_field)[i] = (complex double) (tmp_buffer[i]);
		tau_mean += tmp_buffer[i];
	}
	tau_mean /= N;
	eprintf("\n%f\n", tau_mean);
	/* generate also optical depth contrast */
	for (size_t i = 0; i < N; ++i) {
		(*delta_tau)[i] = (complex double) (tmp_buffer[i] / tau_mean - 1.0);
	}

	/* and flux contrast */

	/* first, the flux field (with scaling. related to its cumulants, because
	 * it's <ln F>) */
	complex double mean_flux = 0.0;
	for (size_t i = 0; i < N; ++i) {
		complex double dflux_i = 0.0;
		/* subtract off the mean optical depth to get hopefully-sensible exponents.
		 * this represents an overall scaling of the flux field, which doesn't
		 * matter when we take the contrast */
		/*
		if ((double) (*tau_field)[i] > 100) {
			eprintf("%f  ", (double) (*tau_field)[i]);
		}
		*/
		dflux_i = ((*tau_field)[i]); //- tau_mean);
		dflux_i = cexp(-dflux_i);
		(*delta_flux)[i] = dflux_i;
		mean_flux += dflux_i;
	}
	mean_flux /= N;
	eprintf("mean flux :%e\n", (double) mean_flux);

	for (size_t i = 0; i < N; ++i) {
		/* and now contrast the flux field */
		(*delta_flux)[i] /= mean_flux;
		(*delta_flux)[i] -= 1.0;
	}


	free(tmp_buffer);

	return 0;

}



int write_xcorr_data(char *field_name_1, char *field_name_2, double *xcorr_k_buffer,
		double *xcorr_output_buffer, size_t *xcorr_count_buffer, size_t xcorr_bin_count)
{
	/* generate output path */
	const size_t path_bsize = 1024;
	char *output_path_buffer = calloc(path_bsize, sizeof(char));
	if (!output_path_buffer) {
		eprintf("Unable to allocate memory to generate output path\n");
		return 1;
	}

	size_t path_length = 0;
	path_length = snprintf(output_path_buffer, path_bsize, "%s/xcorr_%s_%s.dat",
			OUT_DIR, field_name_1, field_name_2);
	if (path_length == path_bsize) {
		eprintf("Over-long output path\n");
		return 1;
	}

	/* open output file */
	FILE *output_fp = fopen(output_path_buffer, "w");
	if (!output_fp) {
		eprintf("Error opening output file %s\n", output_path_buffer);
		return 1;
	}

	/* write output line by line */
	/* format: k xcorr count */
	for (size_t i = 0; i < xcorr_bin_count; ++i) {
		fprintf(output_fp, "%e %e %zd\n", xcorr_k_buffer[i], xcorr_output_buffer[i],
				xcorr_count_buffer[i]);
	}

	/* clean up */
	fclose(output_fp);
	free(output_path_buffer);

	return 0;
}
