#ifndef POWER_SPECTRUM_H_INC
#define POWER_SPECTRUM_H_INC


/* bin_buffer: mean power in the nth bin; k_buffer: the mean k for that bin */
void power_spectrum(complex double *field, size_t KX, double mode_spacing,
		double *k_buffer, double *bin_buffer, size_t *n_buffer, double k_min,
		double k_max, size_t n_bins);


void correlator(complex double *field_1, complex double *field_2, size_t KX,
		double mode_spacing, double h, double *k_buffer, double *bin_buffer,
		size_t *n_buffer, double k_min, double k_max, size_t n_bins);


#endif
