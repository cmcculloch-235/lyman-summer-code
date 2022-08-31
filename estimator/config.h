#ifndef CONFIG_H_INC
#define CONFIG_H_INC

const char* LOSFILE = "../../src/SpecExtract_mpi_with_total/out_80_snap18/los256_n65536_z2.000.dat";
const char* TAUFILE = "../../src/SpecExtract_mpi_with_total/out_80_snap18/tau256_n65536_z2.000.dat";


/* Specify parameters for cross-correlation estimation */
const size_t xcorr_bin_count = 100;
const double xcorr_k_min = 0.01;
const double xcorr_k_max = 1.0;

#endif
