#ifndef CONFIG_H_INC
#define CONFIG_H_INC
#include <math.h>

#define TRUE 1
#define FALSE 0

/* output in Mpc/h etc.? */
#define USE_MPC_H TRUE


/* config parameters are static so they can be included into multiple files
 * without issues of multiple definitions of global variables */

//static const char* LOSFILE = "../../src/SpecExtract_mpi_with_total/out_80_snap18/los256_n65536_z2.000.dat";
//static const char* TAUFILE = "../../src/SpecExtract_mpi_with_total/out_80_snap18/tau256_n65536_z2.000.dat";
//static const char* LOSFILE = "lattice-fields-testing-ver/losfile_test.dat";
//static const char* TAUFILE = "lattice-fields-testing-ver/taufile_test.dat";
//static const char* LOSFILE = "out_80_snap0/los256_n65536_z10.000.dat";
//static const char* TAUFILE = "out_80_snap0/tau256_n65536_z10.000.dat";
//static const char* LOSFILE = "out_80_snap18_new/los256_n65536_z2.000.dat";
//static const char* TAUFILE = "out_80_snap18_new/tau_r.dat";
static const char* LOSFILE = "out_80_ics/los256_n65536_z99.000.dat";
static const char* TAUFILE = "out_80_ics/tau256_n65536_z99.000.dat";

//static const char* OUT_DIR = "out/80_snap18_new/";
//static const char* OUT_DIR = "out/80_snap0";
static const char* OUT_DIR = "out/80_ics/";

/* relevant for dark matter normalisation */
// not anymore; just fix the normalisation by calculating it at runtime from the field.
//static const size_t N_SNAPS = 32;

/* Specify parameters for cross-correlation estimation */
/* Obviously, things get interesting when the number of bins is ~ or > the number
 * of LoSs */
static const size_t xcorr_bin_count = 221;
//static const double xcorr_k_min = 2 * M_PI / 80;
//static const double xcorr_k_max = 128 * 2 * M_PI / 80;
//Pylians binning
static const double xcorr_k_min = 2 * M_PI / 80;
//static const double xcorr_k_max = 1.7379e1;
// 128 * 2pi/80 * sqrt(3) = max possible k
static const double xcorr_k_max =  17.412473;
static const size_t N_PARTICLES = 2048ll * 2048ll * 2048ll;

#endif
