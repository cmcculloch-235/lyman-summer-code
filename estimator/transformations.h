#ifndef TRANS_H_INC
#define TRANS_H_INC

void transform_identity(complex double *field, size_t N, size_t N_THREADS);

void transform_log_normal_contrast(complex double *field, size_t N, size_t N_THREADS);

void transform_power_law_contrast(complex double *field, size_t N, size_t N_THREADS);


void transform_smooth(complex double *field_ksp, size_t KX, double mode_spacing, double scale);

#endif
