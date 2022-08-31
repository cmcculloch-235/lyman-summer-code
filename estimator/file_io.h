#ifndef FILE_IO_H_INC
#define FILE_IO_H_INC

int read_losfile_header(FILE *losfile_fp, double *ztime, double *omegam,
		double *omegal, double *omegab, double *h100, double *box_length,
		double *XH, uint32_t *nbins, uint32_t *nlos);

/* N: no of points on lattice */
/* these functions handle memory allocation, which is why they have **s */
int read_losfile_data(FILE *losfile_fp, size_t N, int nlos, complex double **delta_H,
		complex double **delta_H1, complex double **delta_DM, complex double **delta_matter);

int read_taufile_data(FILE *taufile_fp, size_t N, complex double **tau_field, complex double **delta_tau);

/* yes, the interface is inconsistent with the read routines. I should fix that
 * at some point */
int write_xcorr_data(char *field_name_1, char *field_name_2, double *xcorr_k_buffer,
		double *xcorr_output_buffer, size_t *xcorr_count_buffer, size_t xcorr_bin_count);

#endif
