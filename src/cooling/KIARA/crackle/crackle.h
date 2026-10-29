
#include "crackle_namespace.h"
#include "grackle_types.h"
#include "grackle_chemistry_data.h"

#define TCMB0	2.726

#define EDOT_EXT_FACTOR 1000.0

/* Crackle default floating point precision */
#define cr_float double

typedef struct
{

  /* For Crackle's integration scheme */
  double edot;
  double edot_ext;
  double edot_end;  // edot near the end of the last iteration, for the timestep
  double dedot;
  double HIdot;
  double H2Idot;
  double fSShHI;
  double fSShHeI;
  double fSShHeII;
  double rhoH;
  double rhoHe;
  double rhoH2;
  double nH;
  double mmw;
  double tgas;
  double logtem;
  double tstart;
  double tdust;
  double metallicity;
  double dust2gas;
  double u_cmb;
  double delta_HI;
  double delta_HII;
  double delta_HeI;
  double delta_HeII;
  double delta_HeIII;
  double delta_H2I;
  double delta_HM;
  double delta_H2II;
  double delta_e;
  int verbose;
  int freeze_shielding;  // if set, compute_edot reuses the self-shielding factors and shielded rates already computed
  int id;

  /* Scalar copy of grackle_field_data */

  int grid_rank;
  int grid_dimension;
  int grid_start;
  int grid_end;

  cr_float grid_dx;

  cr_float density;
  cr_float HI_density;
  cr_float HII_density;
  cr_float HM_density;
  cr_float HeI_density;
  cr_float HeII_density;
  cr_float HeIII_density;
  cr_float H2I_density;
  cr_float H2II_density;
  cr_float DI_density;
  cr_float DII_density;
  cr_float HDI_density;
  cr_float e_density;
  cr_float metal_density;
  cr_float dust_density;

  cr_float internal_energy;
  cr_float x_velocity;
  cr_float y_velocity;
  cr_float z_velocity;

  cr_float volumetric_heating_rate;
  cr_float specific_heating_rate;

  cr_float temperature_floor;

  cr_float RT_heating_rate;
  cr_float RT_HI_ionization_rate;
  cr_float RT_HeI_ionization_rate;
  cr_float RT_HeII_ionization_rate;
  cr_float RT_H2_dissociation_rate;
  cr_float H2_self_shielding_length;
  cr_float H2_custom_shielding_factor;

  cr_float isrf_habing;
  cr_float f_shield;

  /* Individual metal+dust species stored in arrays 
   * Species order [0-9]: He, C, N, O, Ne, Mg, Si, D, Ca, Fe */
  cr_float gas_metalDensity[NUM_METAL_SPECIES_GRACKLE];
  cr_float dust_metalDensity[NUM_METAL_SPECIES_GRACKLE];

  /* densities of individual metal species in gas phase
  cr_float He_gas_metalDensity;
  cr_float C_gas_metalDensity;
  cr_float N_gas_metalDensity;
  cr_float O_gas_metalDensity;
  cr_float Ne_gas_metalDensity;
  cr_float Mg_gas_metalDensity;
  cr_float Si_gas_metalDensity;
  cr_float S_gas_metalDensity;
  cr_float Ca_gas_metalDensity;
  cr_float Fe_gas_metalDensity;

  // densities of individual metal species in dust grains
  cr_float He_dust_metalDensity;
  cr_float C_dust_metalDensity;
  cr_float N_dust_metalDensity;
  cr_float O_dust_metalDensity;
  cr_float Ne_dust_metalDensity;
  cr_float Mg_dust_metalDensity;
  cr_float Si_dust_metalDensity;
  cr_float S_dust_metalDensity;
  cr_float Ca_dust_metalDensity;
  cr_float Fe_dust_metalDensity;*/

  cr_float SNe_density; // Crackle needs the number of SNe per unit volume (in code units)

} grackle_part_data;

/**
 * @brief Initialise the evolved (non-input) quantities of a particle before
 * a crackle call; the inputs are set by the caller.
 */
static inline void crackle_init_part(grackle_part_data *gp)
{
	gp->edot = gp->edot_ext = gp->edot_end = gp->dedot = gp->HIdot = gp->H2Idot = gp->fSShHI = gp->fSShHeI = gp->fSShHeII = gp->rhoH = gp->rhoHe = gp->rhoH2 = gp->mmw = gp->tgas = gp->logtem = gp->delta_HI = gp->delta_HII = gp->delta_HeI = gp->delta_HeII = gp->delta_HeIII = gp->delta_HM = gp->delta_H2I = gp->delta_H2II = gp->delta_e = 0.;
	gp->tdust = -1.;  // initialize to negative value so it makes an initial guess
	gp->verbose = 0;
	gp->freeze_shielding = 0;
}

typedef struct {
  /* Parameters for the chemistry rates table */
  double logT;
  double logT_lo;
  double logT_hi;
  double dlogT;
  double dlogT_inv;
  int index;
  double binfrac;
  /* Parameters for the cloudy cooling table */
  double cloudy_par[3];
  int cloudy_index[3];
  double cloudy_dxinv[3];
  //double cloudy_delta[3];
  /* Per-call constants, set in init_temperature_interpolation() */
  double cloudy_w[2];  // interpolation weights in density and redshift (fixed for the call)
  double dust_logT_lo; // dust temperature table
  double dust_dlogT;
  double dust_dlogT_inv;
  double ssh_nfac_HI;  // UVB-dependent part of the HI, HeI self-shielding densities (0 if no UVB)
  double ssh_nfac_HeI;
} interp_struct;

typedef struct {
  double dom;
  double dom_inv;
  double time_to_cgs;
  double length_to_cgs;
  double density_to_cgs;
  double temp_to_K;
  double coolunit;
  double coolunit_inv;
  double velunit;
  double chunit;
  double redshift;
  double zr1;  // redshift plus 1
  double compton1;
  double compton2;  // TCMB at this redshift
  double log10_TCMB;
  double c_ljeans; // constant for H2 self-shielding
} crackle_units;

/******************************************
 *** Chemistry and cooling data storage ***
 ******************************************/
typedef struct
{

  /**********************************
   * primordial chemistry rate data *
   **********************************/

  /* 6 species rates */
  double k1;
  double k2;
  double k3;
  double k4;
  double k5;
  double k6;

  /* 9 species rates (including H2) */
  double k7;
  double k8;
  double k9;
  double k10;
  double k11;
  double k12;
  double k13;
  double k14;
  double k15;
  double k16;
  double k17;
  double k18;
  double k19;
  double k20;  /* currently not used */
  double k21;  /* currently not used */
  double k22;  /* 3-body H2 formation */
  double k23;  /* H2-H2 dissociation */
  double k13dd;  /* density dependent version of k13 (collisional H2
                    dissociation); actually 7 functions instead of 1. */

  /* Radiative rates for 6-species (for external field). */
  double k24;
  double k25;
  double k26;

  /* Radiative rates for 9-species (for external field). */
  double k27;
  double k28;
  double k29;
  double k30;
  double k31;

  /* 12 species rates (with Deuterium). */
  double k50;
  double k51;
  double k52;
  double k53;
  double k54;
  double k55;
  double k56;

  /* New H-ionizing reactions, used for 6, 9 & 12 species chemistry */
  double k57;
  double k58;

  /* Self-shielded rates */
  double k24shield;
  double k25shield;
  double k26shield;
  double k27shield;
  double k28shield;
  double k29shield;
  double k30shield;
  double k31shield;

  /* H2 formation on dust grains */
  double h2dust;

  /* Chemical heating from H2 formation. */
  /* numerator and denominator of Eq 23 of Omukai ea. 2000. */
  double n_cr_n;
  double n_cr_d1;
  double n_cr_d2;

  /********************************
   * primordial cooling rate data *
   ********************************/

  /* 6 species rates */
  double ceHI;                   // collisional excitation rates
  double ceHeI;
  double ceHeII;
  double ciHI;                   // collisional ionization
  double ciHeI;
  double ciHeIS;
  double ciHeII;
  double reHII;                  // recombination
  double reHeII1;
  double reHeII2;
  double reHeIII;
  double brem;                   // free-free (Bremsstrahlung)
  double comp;                    // Compton cooling
  double comp_xray;               // X-ray compton heating coefficient
  double temp_xray;               // X-ray compton heating temperature (K)

  /* radiative rates (external field). */
  double piHI;                    // photo-ionization cooling
  double piHeI;                   //    (no temperature dependance)
  double piHeII;

  // spectrum averaged absorption cross sections
  double crsHI;
  double crsHeI;
  double crsHeII;

  /* 9 species rates (including H2) 
       The first five are for the Lepp & Shull rates.
       The next two are for the (better) Galli & Palla 1999 rates. 
       The selection is controlled by a flag in cool1d_multi_g.F. */
  double hyd01k;
  double h2k01;
  double vibh;
  double roth;
  double rotl;
  double GP99LowDensityLimit;
  double GP99HighDensityLimit;

  /* Revised H2 cooling rates from Glover & Abel 2008 */
  double GAHI;
  double GAH2;
  double GAHe;
  double GAHp;
  double GAel;

  /* Updated H2 LTE rate from Glover (2015, MNRAS, 451, 2082) */
  double H2LTE;

  /* 12 species rates (including HD) */
  double HDlte;
  double HDlow;

  /* CIE cooling */
  double cieco;

  /*******************************
   * dust chemistry/cooling data *
   *******************************/

  // Photo-electric heating (code units)
  double gammah;

  // Electron recombination onto dust grains
  double regr;

  // Heating of dust by interstellar radiation field
  double gamma_isrf;

  /* Gas/grain energy transfer. */
  double gas_grain;

} chemistry_rate_storage;

/* New stuff for crackle */

int crackle_solve_chemistry(grackle_field_data *p, chemistry_data *chemistry, chemistry_data_storage grackle_rates, photo_rate_storage uvb_rates, code_units *units, double dt);

/* Direct per-particle API: the caller fills a grackle_part_data (after
 * crackle_init_part) and prepares the UV background rates, units and crackle
 * units once per step with crackle_prepare_step. */
int crackle_prepare_step(chemistry_data *chemistry, const chemistry_data_storage *grackle_rates, code_units *units, photo_rate_storage *uvb_rates, crackle_units *cunits);
int crackle_get_uvb_rates(chemistry_data *chemistry, chemistry_data_storage *grackle_rates, photo_rate_storage *uvb_rates, code_units *units);
int crackle_solve_part(grackle_part_data *gp, chemistry_data *chemistry, const chemistry_data_storage *grackle_rates, const photo_rate_storage *uvb_rates, code_units *units, const crackle_units *cunits, double dt);
void crackle_cooling_time_part(grackle_part_data *gp, chemistry_data *chemistry, const chemistry_data_storage *grackle_rates, const photo_rate_storage *uvb_rates, code_units *units, const crackle_units *cunits, gr_float *tcool);
void crackle_temperature_part(grackle_part_data *gp, chemistry_data *chemistry, code_units *units, gr_float *temp);

void evolve_internal_energy(grackle_part_data *gp, chemistry_data *chemistry, double dtit);

void evolve_energy_ros2(grackle_part_data *gp, double u0, double dtit, chemistry_data *chemistry, const chemistry_data_storage *grackle_rates, chemistry_rate_storage *my_rates, const photo_rate_storage *uvb_rates, interp_struct *interpolation, code_units *units, const crackle_units *cunits, int ism_flag);

void compute_edot(grackle_part_data *gp, chemistry_data *chemistry, const chemistry_data_storage *rates_table, chemistry_rate_storage *my_rates, const photo_rate_storage *uvb_rates, interp_struct *interpolation, code_units *units, const crackle_units *cunits, int ism_flag);

double compute_dedot(int chemistry_flag, const grackle_part_data *gp, chemistry_data *chemistry, const chemistry_rate_storage *my_rates, code_units *units);

double compute_HIdot(int chemistry_flag, const grackle_part_data *gp, chemistry_data *chemistry, const chemistry_rate_storage *my_rates, code_units *units);

double compute_iteration_dt(grackle_part_data *gp, double edot_prev, chemistry_data *chemistry, double dt, double dtcool, float accuracy, double *dtsuppress);

void evolve_helium(grackle_part_data *gp, grackle_part_data *gp_old, chemistry_data *chemistry, const chemistry_rate_storage *my_rates, double dtit);

void evolve_hydrogen(grackle_part_data *gp, grackle_part_data *gp_old, chemistry_data *chemistry, const chemistry_rate_storage *my_rates, double dtit);

void evolve_H2(grackle_part_data *gp, int ism_flag, chemistry_data *chemistry, const chemistry_rate_storage *my_rates, const crackle_units *cunits, double dtit);

void evolve_elements(grackle_part_data *gp, grackle_part_data *gp_old, chemistry_data *chemistry);

void evolve_pred_corr(grackle_part_data *gp, grackle_part_data *gp_old, chemistry_data *chemistry);

void save_species_data(grackle_part_data *gp_orig, grackle_part_data *gp_target);

//void evolve_dust(grackle_part_data *gp, chemistry_data *chemistry, code_units *units, int ism_flag, double dtit);

//double dust_thermal_balance(double tdust, double tgas, double trad, double trad4, double gamma_isrf, double gasgr, double nh);

//double calculate_dust_temp(double tgas, double nh, double gasgr, double gamma_isrf, double trad, double tdust);

void crackle_cooling_time(grackle_field_data *p, chemistry_data *chemistry, chemistry_data_storage grackle_rates, photo_rate_storage uvb_rates, code_units *units, gr_float *tcool);

void crackle_temperature(grackle_field_data *p, chemistry_data *chemistry, code_units *units, gr_float *temp);

//double cloudy_metal_cooling(double tgas, double rhoH, double zr, double dom, double metallicity, chemistry_data_storage gr);

//double interpolate_3d(double *par, int *ind, double *delta, double *cloudy_data, long long int *griddim, double *gridpar[3]);
