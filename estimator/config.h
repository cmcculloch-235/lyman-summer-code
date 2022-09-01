#ifndef CONFIG_H_INC
#define CONFIG_H_INC

/* config parameters are static so they can be included into multiple files
 * without issues of multiple definitions of global variables */

//static const char* LOSFILE = "../../src/SpecExtract_mpi_with_total/out_80_snap18/los256_n65536_z2.000.dat";
//static const char* TAUFILE = "../../src/SpecExtract_mpi_with_total/out_80_snap18/tau256_n65536_z2.000.dat";
static const char* LOSFILE = "lattice-fields-testing-ver/losfile_test.dat";
static const char* TAUFILE = "lattice-fields-testing-ver/taufile_test.dat";

static const char* OUT_DIR = "out";


/* Specify parameters for cross-correlation estimation */
static const size_t xcorr_bin_count = 100;
static const double xcorr_k_min = 0.1;
static const double xcorr_k_max = 1;

#endif
