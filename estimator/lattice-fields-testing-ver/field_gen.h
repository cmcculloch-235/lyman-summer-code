#ifndef FIELD_GEN_H_INC
#define FIELD_GEN_H_INC


double spec_flat(double k);
double spec_delta(double k);
double spec_gaussian(double k);
double spec_linear(double k);
double spec_bbks(double k);


/* grid spacing in k-space is mode_spacing */
void gen_field(fftw_complex *field_buffer, size_t KX, double mode_spacing,
		double (*power_spectrum) (double));

#endif
