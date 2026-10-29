/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
/**
 * @file src/cooling/KIARA/cooling.c
 * @brief Cooling using the GRACKLE 3.x library.
 */

#include <config.h>

/* Some standard headers. */
#include <float.h>
#include <hdf5.h>
#include <math.h>
#include <stdint.h>
#include <string.h>
#include <time.h>

/* Include header */
#include "cooling.h"

/* Some standard headers. */
#include <assert.h>
#include <fenv.h>
#include <float.h>
#include <math.h>

/* The grackle library itself */
#include "crackle/grackle.h"
#include "crackle_validate.h"
extern chemistry_data *grackle_data;

/* Local includes. */
#include "chemistry.h"
#include "cooling_io.h"
#include "entropy_floor.h"
#include "error.h"
#include "fof.h"
#include "hydro.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "star_formation.h"
#include "units.h"

/* crackle's (Grackle's) element arrays (gas_metalDensity, dust_metalDensity,
 * SolarAbundances) hold He, C, N, O, Ne, Mg, Si, S, Ca, Fe: the slot of each
 * of chem5's elements in them, -1 if crackle does not know it */
#define CRACKLE_SLOT_H -1
#define CRACKLE_SLOT_He 0
#define CRACKLE_SLOT_Li -1
#define CRACKLE_SLOT_Be -1
#define CRACKLE_SLOT_B -1
#define CRACKLE_SLOT_C 1
#define CRACKLE_SLOT_N 2
#define CRACKLE_SLOT_O 3
#define CRACKLE_SLOT_F -1
#define CRACKLE_SLOT_Ne 4
#define CRACKLE_SLOT_Na -1
#define CRACKLE_SLOT_Mg 5
#define CRACKLE_SLOT_Al -1
#define CRACKLE_SLOT_Si 6
#define CRACKLE_SLOT_P -1
#define CRACKLE_SLOT_S 7
#define CRACKLE_SLOT_Cl -1
#define CRACKLE_SLOT_Ar -1
#define CRACKLE_SLOT_K -1
#define CRACKLE_SLOT_Ca 8
#define CRACKLE_SLOT_Sc -1
#define CRACKLE_SLOT_Ti -1
#define CRACKLE_SLOT_V -1
#define CRACKLE_SLOT_Cr -1
#define CRACKLE_SLOT_Mn -1
#define CRACKLE_SLOT_Fe 9
#define CRACKLE_SLOT_Co -1
#define CRACKLE_SLOT_Ni -1
#define CRACKLE_SLOT_Cu -1
#define CRACKLE_SLOT_Zn -1
#define CRACKLE_SLOT_Ga -1
#define CRACKLE_SLOT_Ge -1

/* crackle's dust model only knows its own elements */
#define CRACKLE_CHECK_DUST_SLOT(e)                    \
  _Static_assert(CRACKLE_SLOT_##e >= 0, "dust element " #e \
                                        " is not one of crackle's elements");
KIARA_DUST_ELEMENTS(CRACKLE_CHECK_DUST_SLOT)
#undef CRACKLE_CHECK_DUST_SLOT

/**
 * @brief The slot of a tracked element in crackle's element arrays, or -1 if
 * crackle does not know it (it is then left unchanged by crackle).
 *
 * @param elem The #chemistry_element.
 */
__attribute__((always_inline)) INLINE static int cooling_crackle_slot(
    const int elem) {
  static const int slot[chemistry_element_count] = {
#define CRACKLE_SLOT_OF(e) CRACKLE_SLOT_##e,
      KIARA_ELEMENTS(CRACKLE_SLOT_OF)
#undef CRACKLE_SLOT_OF
  };
  return slot[elem];
}

/**
 * @brief Common operations performed on the cooling function at a
 * given time-step or redshift. Predominantly used to read cooling tables
 * above and below the current redshift, if not already read in.
 *
 * Also calls the additional H reionisation energy injection if need be.
 *
 * @param phys_const The physical constants in internal units.
 * @param cosmo The current cosmological model.
 * @param pressure_floor Properties of the pressure floor.
 * @param cooling The #cooling_function_data used in the run.
 * @param s The space data, including a pointer to array of particles
 * @param time The current system time
 */
void cooling_update(const struct phys_const *phys_const,
                    const struct cosmology *cosmo,
                    const struct pressure_floor_props *pressure_floor,
                    struct cooling_function_data *cooling, struct space *s,
                    const double time) {

  /* set current time */
  if (cooling->redshift == -1) {
    cooling->units.a_value = cosmo->a;
  } else {
    cooling->units.a_value = 1. / (1. + cooling->redshift);
  }

  /* Everything crackle needs that depends only on the redshift */
  cooling_prepare_crackle_step(cooling);
}

/**
 * @brief Print the chemical network
 *
 * @param xp The #xpart to print
 */
void cooling_print_fractions(const struct xpart *restrict xp) {

  const struct cooling_xpart_data tmp = xp->cooling_data;
#if COOLING_GRACKLE_MODE > 0
  message("HI %g, HII %g, HeI %g, HeII %g, HeIII %g, e %g", tmp.HI_frac,
          tmp.HII_frac, tmp.HeI_frac, tmp.HeII_frac, tmp.HeIII_frac,
          tmp.e_frac);
#endif

#if COOLING_GRACKLE_MODE > 1
  message("HM %g, H2I %g, H2II %g", tmp.HM_frac, tmp.H2I_frac, tmp.H2II_frac);
#endif

#if COOLING_GRACKLE_MODE > 2
  message("DI %g, DII %g, HDI %g", tmp.DI_frac, tmp.DII_frac, tmp.HDI_frac);
#endif
  message("Metal: %g", tmp.metal_frac);
}

/**
 * @brief Initializes grackle particle quantities
 * assuming purely neutral gas
 *
 * H and He start fully neutral with the particle's element mass fractions,
 * all other species at zero (deuterium at the primordial D/H ratio).
 *
 * @param cooling The properties of the cooling function.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 */
__attribute__((always_inline)) INLINE void cooling_grackle_init_part(
    const struct cooling_function_data *cooling, struct part *restrict p,
    struct xpart *restrict xp) {

#if COOLING_GRACKLE_MODE >= 1
  gr_float zero = 1.e-20;
  zero = 0.f;

  /* primordial chemistry >= 1: Start with everything neutral (as in dark ages)
   */
  xp->cooling_data.HI_frac =
      p->chemistry_data.metal_mass_fraction[chemistry_element_H];
  xp->cooling_data.HII_frac = zero;
  xp->cooling_data.HeI_frac =
      p->chemistry_data.metal_mass_fraction[chemistry_element_He];
  xp->cooling_data.HeII_frac = zero;
  xp->cooling_data.HeIII_frac = zero;
  xp->cooling_data.e_frac = xp->cooling_data.HII_frac +
                            0.25 * xp->cooling_data.HeII_frac +
                            0.5 * xp->cooling_data.HeIII_frac;
#endif  // MODE >= 1

#if COOLING_GRACKLE_MODE >= 2
  /* primordial chemistry >= 2 */
  xp->cooling_data.HM_frac = zero;
  xp->cooling_data.H2I_frac = zero;
  xp->cooling_data.H2II_frac = zero;
#endif  // MODE >= 2

#if COOLING_GRACKLE_MODE >= 3
  /* primordial chemistry >= 3 */
  xp->cooling_data.DI_frac = grackle_data->DeuteriumToHydrogenRatio *
                             grackle_data->HydrogenFractionByMass;
  xp->cooling_data.DII_frac = zero;
  xp->cooling_data.HDI_frac = zero;
#endif  // MODE >= 3
}

/**
 * @brief Sets the cooling properties of the (x-)particles to a valid start
 * state.
 *
 * @param phys_const The #phys_const.
 * @param us The #unit_system.
 * @param hydro_props The #hydro_props.
 * @param cosmo The #cosmology.
 * @param cooling The properties of the cooling function.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 */
void cooling_first_init_part(const struct phys_const *restrict phys_const,
                             const struct unit_system *restrict us,
                             const struct hydro_props *hydro_props,
                             const struct cosmology *restrict cosmo,
                             const struct cooling_function_data *cooling,
                             struct part *restrict p,
                             struct xpart *restrict xp) {

  xp->cooling_data.radiated_energy = 0.f;
  xp->cooling_data.time_last_event = -cooling->thermal_time;

  /* Initialize grackle ionization fractions */
  cooling_grackle_init_part(cooling, p, xp);

  p->cooling_data.subgrid_temp = 0.f;
  p->cooling_data.subgrid_dens = 0.f;
  p->cooling_data.subgrid_fcold = 0.f;
  p->cooling_data.G0 = 0.f;
  p->cooling_data.mixing_layer_cool_time = 0.f;

  /* Initialize dust properties */
#if COOLING_GRACKLE_MODE >= 2
  p->cooling_data.dust_mass = 0.f;
  for (int d = 0; d < dust_element_count; d++) {
    p->cooling_data.dust_mass_fraction[d] = 0.f;
  }

  p->cooling_data.dust_temperature = 0.f;
#endif

  /* Initialize debug properties */
  p->cooling_data.du_dt_cooling = 0;
  p->cooling_data.u_cooling = 0;
  p->cooling_data.u_chemistry = 0;
  p->cooling_data.u_cooling_before = 0;
  p->cooling_data.u_wind_after = 0;
  p->cooling_data.u_wind_before = 0;
  p->cooling_data.u_sn_after = 0;
  p->cooling_data.u_sn_before = 0;
  p->cooling_data.du_wind_this_step = 0;
  p->cooling_data.du_sn_this_step = 0;
}

/**
 * @brief Perform additional init on the cooling properties of the
 * (x-)particles that requires the density to be known.
 *
 * Nothing to do here.
 *
 * @param phys_const The physical constant in internal units.
 * @param us The unit system.
 * @param hydro_props The properties of the hydro scheme.
 * @param cosmo The current cosmological model.
 * @param cooling The properties of the cooling function.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 */
__attribute__((always_inline)) INLINE void cooling_post_init_part(
    const struct phys_const *restrict phys_const,
    const struct unit_system *restrict us,
    const struct hydro_props *hydro_props,
    const struct cosmology *restrict cosmo,
    const struct cooling_function_data *cooling, const struct part *restrict p,
    struct xpart *restrict xp) {}

/**
 * @brief Returns the total radiated energy by this particle.
 *
 * @param xp The extended particle data
 * @return The total radiated energy (internal units).
 */
float cooling_get_radiated_energy(const struct xpart *restrict xp) {

  return xp->cooling_data.radiated_energy;
}

/**
 * @brief Prints the properties of the cooling model to stdout.
 *
 * @param cooling The properties of the cooling function.
 */
void cooling_print_backend(const struct cooling_function_data *cooling) {

  if (cooling->chemistry.use_grackle == 1) {
    message("Cooling function is 'Grackle', mode = %i",
            cooling->chemistry.primordial_chemistry);
  } else if (cooling->chemistry.use_grackle == 2) {
    message("Cooling function is 'Crackle', mode = %i",
            cooling->chemistry.primordial_chemistry);
  }

  message("CloudyTable = %s", cooling->cloudy_table);
  message("Redshift = %g", cooling->redshift);
  message("UV background flag = %d", cooling->with_uv_background);
  message("Metal cooling flag = %i", cooling->chemistry.metal_cooling);
  message("Self Shielding flag = %i", cooling->self_shielding_method);
  message("Max subgrid density (internal units) = %g",
          cooling->max_subgrid_density);
  message("Thermal time = %g", cooling->thermal_time);
  message("Specific Heating Rates flag = %i",
          cooling->provide_specific_heating_rates);
  message("Volumetric Heating Rates flag = %i",
          cooling->provide_volumetric_heating_rates);
  message("Units:");
  message("\tComoving = %i", cooling->units.comoving_coordinates);
  message("\tLength = %g", cooling->units.length_units);
  message("\tNumber Density = %g", cooling->units.density_units);
  message("\tTime = %g", cooling->units.time_units);
  message("\tScale Factor = %g (units: %g)", cooling->units.a_value,
          cooling->units.a_units);
}

#if COOLING_GRACKLE_MODE >= 1
/**
 * @brief Copy the primordial (H, He, e) species of a #xpart to the grackle
 * data.
 *
 * The species densities are the stored mass fractions times the density, and
 * are kept in species_densities, which the grackle data point to.
 *
 * @param data The grackle_field_data structure from grackle.
 * @param p The #part.
 * @param xp The #xpart.
 * @param rho The particle density (grackle units).
 * @param species_densities Buffer holding the densities passed to grackle.
 */
void cooling_copy_to_grackle1(grackle_field_data *data, const struct part *p,
                              const struct xpart *xp, gr_float rho,
                              gr_float species_densities[N_SPECIES]) {
  /* HI */
  species_densities[0] = xp->cooling_data.HI_frac * rho;
  data->HI_density = &species_densities[0];

  /* HII */
  species_densities[1] = xp->cooling_data.HII_frac * rho;
  data->HII_density = &species_densities[1];

  /* HeI */
  species_densities[2] = xp->cooling_data.HeI_frac * rho;
  data->HeI_density = &species_densities[2];

  /* HeII */
  species_densities[3] = xp->cooling_data.HeII_frac * rho;
  data->HeII_density = &species_densities[3];

  /* HeIII */
  species_densities[4] = xp->cooling_data.HeIII_frac * rho;
  data->HeIII_density = &species_densities[4];

  /* electrons */
  species_densities[5] = xp->cooling_data.e_frac * rho;
  data->e_density = &species_densities[5];
}
#else
/**
 * @brief Copy the primordial species to the grackle data (none in mode 0).
 *
 * @param data The grackle_field_data structure from grackle.
 * @param p The #part.
 * @param xp The #xpart.
 * @param rho The particle density (grackle units).
 * @param species_densities Buffer holding the densities passed to grackle.
 */
void cooling_copy_to_grackle1(grackle_field_data *data, const struct part *p,
                              const struct xpart *xp, gr_float rho,
                              gr_float species_densities[N_SPECIES]) {
  data->HI_density = NULL;
  data->HII_density = NULL;
  data->HeI_density = NULL;
  data->HeII_density = NULL;
  data->HeIII_density = NULL;
  data->e_density = NULL;
}
#endif

#if COOLING_GRACKLE_MODE >= 2
/**
 * @brief Copy the H2 species and the dust-model inputs of a particle to the
 * grackle data.
 *
 * With the dust model on, this also passes the dust density, the number of
 * SNe in this step, the interstellar radiation field G0 (-1 outside the
 * subgrid ISM), the H2 self-shielding length (rho/|grad rho|, limited to
 * [MIN_SHIELD_H_FRAC h, kernel_gamma h]) and the gas and dust element
 * densities, and sets the slot in which crackle returns its dust temperature.
 * Requires the metal density (species_densities[12]) to be set.
 *
 * @param data The grackle_field_data structure from grackle.
 * @param p The #part.
 * @param xp The #xpart.
 * @param cooling The #cooling_function_data used in the run.
 * @param dt The time-step of this particle.
 * @param rho The particle density (grackle units).
 * @param species_densities Buffer holding the densities passed to grackle.
 */
void cooling_copy_to_grackle2(
    grackle_field_data *data, const struct part *p, const struct xpart *xp,
    const struct cooling_function_data *restrict cooling, const double dt,
    gr_float rho, gr_float species_densities[N_SPECIES]) {
  /* HM */
  species_densities[6] = xp->cooling_data.HM_frac * rho;
  data->HM_density = &species_densities[6];

  /* H2I */
  species_densities[7] = xp->cooling_data.H2I_frac * rho;
  data->H2I_density = &species_densities[7];

  /* H2II */
  species_densities[8] = xp->cooling_data.H2II_frac * rho;
  data->H2II_density = &species_densities[8];

  /* Dust model */
  if (cooling->use_grackle_dust_evol == 1) {
    /* Load dust and metal info */
    species_densities[20] =
        p->cooling_data.dust_mass / p->mass * species_densities[12];
    data->dust_density = &species_densities[20];
    species_densities[21] = p->feedback_data.SNe_ThisTimeStep * dt / p->mass *
                            species_densities[12];
    /* need to pass the number of SNe per volume of particle;
     * recall SNe_ThisTimeStep is the SNe rate */

    // if (species_densities[21] > 1.e15) message("SNe_density: %g %g %g
    // %g\n",p->mass, species_densities[12], p->feedback_data.SNe_ThisTimeStep,
    // species_densities[21]);
    data->SNe_ThisTimeStep = &species_densities[21];
    // if( chemistry_get_total_metal_mass_fraction_for_cooling(p)>0.f)
    // message("Zsm= %g Zp= %g Z= %g Zd=
    // %g",chemistry_get_total_metal_mass_fraction_for_cooling(p),
    // p->chemistry_data.metal_mass_fraction_total, species_densities[19],
    // species_densities[20]);

    /* Determine ISRF in Habing units based on chosen method, -1 == non-ISM */
    if (p->cooling_data.subgrid_temp == 0.f) {
      species_densities[22] = -1.f;
    } else {
      species_densities[22] = p->cooling_data.G0;
    }

    data->isrf_habing = &species_densities[22];

    const double rho_grad_norm2 = p->rho_gradient[0] * p->rho_gradient[0] +
                                  p->rho_gradient[1] * p->rho_gradient[1] +
                                  p->rho_gradient[2] * p->rho_gradient[2];
    const double rho_grad_norm_inv =
        (rho_grad_norm2 > 0.) ? 1. / sqrt(rho_grad_norm2) : 0.;
    const double rho_com = hydro_get_comoving_density(p);
    double L_eff_com = rho_com * rho_grad_norm_inv;
    const double L_eff_com_max = kernel_gamma * p->h;
    const double L_eff_com_min = MIN_SHIELD_H_FRAC * p->h;
    L_eff_com = fmin(L_eff_com, L_eff_com_max);
    L_eff_com = fmax(L_eff_com, L_eff_com_min);

    species_densities[23] = L_eff_com * cooling->units.a_value;
    data->H2_self_shielding_length = &species_densities[23];

    /* Load gas metallicities */
    /* (by crackle slot: gas He..Fe in 25..34, dust in 36..45) */
    for (int k = 24; k < 46; k++) species_densities[k] = 0.;
    for (int i = 0; i < chemistry_element_count; i++) {
      const int slot = cooling_crackle_slot(i);
      if (slot < 0) continue;
      species_densities[25 + slot] = max(
          p->chemistry_data.metal_mass_fraction[i] * species_densities[12], 0.);
    }
    for (int d = 0; d < dust_element_count; d++) {
      const int slot = cooling_crackle_slot(chemistry_element_of_dust(d));
      species_densities[36 + slot] =
          max(p->cooling_data.dust_mass_fraction[d] * species_densities[20], 0);
      // if (i>0) printf("dust densities: %d %g %g
      // %g\n",i,species_densities[23+chemistry_element_count+i],species_densities[23+chemistry_element_count+i]/data->dust_density[0],p->cooling_data.dust_mass_fraction[i]*p->mass
      // / p->cooling_data.dust_mass);
    }

    data->He_gas_metalDensity = &species_densities[25];
    data->C_gas_metalDensity = &species_densities[26];
    data->N_gas_metalDensity = &species_densities[27];
    data->O_gas_metalDensity = &species_densities[28];
    data->Ne_gas_metalDensity = &species_densities[29];
    data->Mg_gas_metalDensity = &species_densities[30];
    data->Si_gas_metalDensity = &species_densities[31];
    data->S_gas_metalDensity = &species_densities[32];
    data->Ca_gas_metalDensity = &species_densities[33];
    data->Fe_gas_metalDensity = &species_densities[34];
    /* Load dust metallicities */
    data->He_dust_metalDensity = &species_densities[36];
    data->C_dust_metalDensity = &species_densities[37];
    data->N_dust_metalDensity = &species_densities[38];
    data->O_dust_metalDensity = &species_densities[39];
    data->Ne_dust_metalDensity = &species_densities[40];
    data->Mg_dust_metalDensity = &species_densities[41];
    data->Si_dust_metalDensity = &species_densities[42];
    data->S_dust_metalDensity = &species_densities[43];
    data->Ca_dust_metalDensity = &species_densities[44];
    data->Fe_dust_metalDensity = &species_densities[45];

    /* Dust temperature returned by crackle (unchanged if not computed) */
    species_densities[46] = p->cooling_data.dust_temperature;
    data->dust_temperature = &species_densities[46];
  } else {
    data->dust_density = NULL;
    data->SNe_ThisTimeStep = NULL;
    data->isrf_habing = NULL;
    data->He_gas_metalDensity = NULL;
    data->C_gas_metalDensity = NULL;
    data->N_gas_metalDensity = NULL;
    data->O_gas_metalDensity = NULL;
    data->Ne_gas_metalDensity = NULL;
    data->Mg_gas_metalDensity = NULL;
    data->Si_gas_metalDensity = NULL;
    data->S_gas_metalDensity = NULL;
    data->Ca_gas_metalDensity = NULL;
    data->Fe_gas_metalDensity = NULL;
    data->He_dust_metalDensity = NULL;
    data->C_dust_metalDensity = NULL;
    data->N_dust_metalDensity = NULL;
    data->O_dust_metalDensity = NULL;
    data->Ne_dust_metalDensity = NULL;
    data->Mg_dust_metalDensity = NULL;
    data->Si_dust_metalDensity = NULL;
    data->S_dust_metalDensity = NULL;
    data->Ca_dust_metalDensity = NULL;
    data->Fe_dust_metalDensity = NULL;
  }
}
#else
/**
 * @brief Copy the H2 species to the grackle data (none below mode 2).
 *
 * @param data The grackle_field_data structure from grackle.
 * @param p The #part.
 * @param xp The #xpart.
 * @param cooling The #cooling_function_data used in the run.
 * @param dt The time-step of this particle.
 * @param rho The particle density (grackle units).
 * @param species_densities Buffer holding the densities passed to grackle.
 */
void cooling_copy_to_grackle2(
    grackle_field_data *data, const struct part *p, const struct xpart *xp,
    const struct cooling_function_data *restrict cooling, const double dt,
    gr_float rho, gr_float species_densities[N_SPECIES]) {
  data->HM_density = NULL;
  data->H2I_density = NULL;
  data->H2II_density = NULL;
}
#endif

#if COOLING_GRACKLE_MODE >= 3
/**
 * @brief Copy the deuterium species of a #xpart to the grackle data.
 *
 * @param data The grackle_field_data structure from grackle.
 * @param p The #part.
 * @param xp The #xpart.
 * @param rho The particle density (grackle units).
 * @param species_densities Buffer holding the densities passed to grackle.
 */
void cooling_copy_to_grackle3(grackle_field_data *data, const struct part *p,
                              const struct xpart *xp, gr_float rho,
                              gr_float species_densities[N_SPECIES]) {
  /* DI */
  species_densities[9] = xp->cooling_data.DI_frac * rho;
  data->DI_density = &species_densities[9];

  /* DII */
  species_densities[10] = xp->cooling_data.DII_frac * rho;
  data->DII_density = &species_densities[10];

  /* HDI */
  species_densities[11] = xp->cooling_data.HDI_frac * rho;
  data->HDI_density = &species_densities[11];
}
#else
/**
 * @brief Copy the deuterium species to the grackle data (none below mode 3).
 *
 * @param data The grackle_field_data structure from grackle.
 * @param p The #part.
 * @param xp The #xpart.
 * @param rho The particle density (grackle units).
 * @param species_densities Buffer holding the densities passed to grackle.
 */
void cooling_copy_to_grackle3(grackle_field_data *data, const struct part *p,
                              const struct xpart *xp, gr_float rho,
                              gr_float species_densities[N_SPECIES]) {
  data->DI_density = NULL;
  data->DII_density = NULL;
  data->HDI_density = NULL;
}
#endif

#if COOLING_GRACKLE_MODE >= 1
/**
 * @brief Copy the primordial (H, He, e) species from the grackle data back
 * to the mass fractions of a #xpart.
 *
 * @param data The grackle_field_data structure from grackle.
 * @param p The #part.
 * @param xp The #xpart.
 * @param rho The particle density (grackle units).
 */
void cooling_copy_from_grackle1(grackle_field_data *data, const struct part *p,
                                struct xpart *xp, gr_float rho) {

  const double rhoinv = 1.f / rho;
  /* HI */
  xp->cooling_data.HI_frac = *data->HI_density * rhoinv;

  /* HII */
  xp->cooling_data.HII_frac = *data->HII_density * rhoinv;

  /* HeI */
  xp->cooling_data.HeI_frac = *data->HeI_density * rhoinv;

  /* HeII */
  xp->cooling_data.HeII_frac = *data->HeII_density * rhoinv;

  /* HeIII */
  xp->cooling_data.HeIII_frac = *data->HeIII_density * rhoinv;

  /* e */
  xp->cooling_data.e_frac = *data->e_density * rhoinv;
}
#else
/**
 * @brief Copy the primordial species from the grackle data (none in mode 0).
 *
 * @param data The grackle_field_data structure from grackle.
 * @param p The #part.
 * @param xp The #xpart.
 * @param rho The particle density (grackle units).
 */
void cooling_copy_from_grackle1(grackle_field_data *data, const struct part *p,
                                struct xpart *xp, gr_float rho) {}
#endif

#if COOLING_GRACKLE_MODE >= 2
/**
 * @brief Copy the H2 species and, with the dust model on, the evolved gas
 * metal and dust masses from the grackle data back to the particle.
 *
 * The gas element mass fractions and total metallicity are overwritten with
 * grackle's values (which include dust growth and destruction), and the dust
 * mass and per-element dust mass fractions are updated (index 0 holds the
 * total).
 *
 * @param data The grackle_field_data structure from grackle.
 * @param p The #part.
 * @param xp The #xpart.
 * @param cooling The #cooling_function_data used in the run.
 * @param rho The particle density (grackle units).
 */
void cooling_copy_from_grackle2(
    grackle_field_data *data, struct part *p, struct xpart *xp,
    const struct cooling_function_data *restrict cooling, gr_float rho) {
  double rhoinv = 1.f / rho;

  /* HM */
  xp->cooling_data.HM_frac = *data->HM_density * rhoinv;
  /* H2I */
  xp->cooling_data.H2I_frac = *data->H2I_density * rhoinv;
  /* H2II */
  xp->cooling_data.H2II_frac = *data->H2II_density * rhoinv;

  /* Dust model */
  if (cooling->use_grackle_dust_evol == 1) {
    /* Load gas metallicities */
    p->chemistry_data.metal_mass_fraction_total = *data->metal_density * rhoinv;
    const gr_float *rhogas[NUM_METAL_SPECIES_GRACKLE] = {
        data->He_gas_metalDensity, data->C_gas_metalDensity,
        data->N_gas_metalDensity,  data->O_gas_metalDensity,
        data->Ne_gas_metalDensity, data->Mg_gas_metalDensity,
        data->Si_gas_metalDensity, data->S_gas_metalDensity,
        data->Ca_gas_metalDensity, data->Fe_gas_metalDensity};
    for (int i = 0; i < chemistry_element_count; i++) {
      const int slot = cooling_crackle_slot(i);
      if (slot >= 0)
        p->chemistry_data.metal_mass_fraction[i] = *rhogas[slot] * rhoinv;
    }

    /* Dust temperature used by crackle */
    p->cooling_data.dust_temperature = *data->dust_temperature;

    /* Load dust metallicities */
    p->cooling_data.dust_mass = *data->dust_density * p->mass * rhoinv;

    if (p->cooling_data.dust_mass > 0.5 * p->mass) {
      warning("DUST > METALS Mg=%g Zg=%g mdust=%g mmet=%g\n", p->mass,
              chemistry_get_total_metal_mass_fraction_for_cooling(p),
              p->cooling_data.dust_mass,
              p->mass * chemistry_get_total_metal_mass_fraction_for_cooling(p));
    }

    if (*data->dust_density > 0.f) {
      rhoinv = 1. / *data->dust_density;

      const double rhodust[chemistry_element_count - 1] = {
          *data->He_dust_metalDensity, *data->C_dust_metalDensity,
          *data->N_dust_metalDensity,  *data->O_dust_metalDensity,
          *data->Ne_dust_metalDensity, *data->Mg_dust_metalDensity,
          *data->Si_dust_metalDensity, *data->S_dust_metalDensity,
          *data->Ca_dust_metalDensity, *data->Fe_dust_metalDensity};

      /* no Helium metal density */
      assert(rhodust[0] == 0.);

      for (int d = 0; d < dust_element_count; d++) {
        p->cooling_data.dust_mass_fraction[d] =
            rhodust[cooling_crackle_slot(chemistry_element_of_dust(d))] *
            rhoinv;
      }
    } else {
      for (int d = 0; d < dust_element_count; d++) {
        p->cooling_data.dust_mass_fraction[d] = 0.f;
      }
    }
  }
}
#else
/**
 * @brief Copy the H2 species from the grackle data (none below mode 2).
 *
 * @param data The grackle_field_data structure from grackle.
 * @param p The #part.
 * @param xp The #xpart.
 * @param cooling The #cooling_function_data used in the run.
 * @param rho The particle density (grackle units).
 */
void cooling_copy_from_grackle2(
    grackle_field_data *data, struct part *p, struct xpart *xp,
    const struct cooling_function_data *restrict cooling, gr_float rho) {}
#endif

#if COOLING_GRACKLE_MODE > 2
/**
 * @brief Copy the deuterium species from the grackle data back to a #xpart.
 *
 * @param data The grackle_field_data structure from grackle.
 * @param p The #part.
 * @param xp The #xpart.
 * @param rho The particle density (grackle units).
 */
void cooling_copy_from_grackle3(grackle_field_data *data, const struct part *p,
                                struct xpart *xp, gr_float rho) {

  const double rhoinv = 1.f / rho;
  /* DI */
  xp->cooling_data.DI_frac = *data->DI_density * rhoinv;

  /* DII */
  xp->cooling_data.DII_frac = *data->DII_density * rhoinv;

  /* HDI */
  xp->cooling_data.HDI_frac = *data->HDI_density * rhoinv;
}
#else
/**
 * @brief Copy the deuterium species from the grackle data (none below mode
 * 3).
 *
 * @param data The grackle_field_data structure from grackle.
 * @param p The #part.
 * @param xp The #xpart.
 * @param rho The particle density (grackle units).
 */
void cooling_copy_from_grackle3(grackle_field_data *data, const struct part *p,
                                struct xpart *xp, gr_float rho) {}
#endif

/**
 * @brief Check that a value is finite (not NaN or inf) from its bits.
 *
 * fpclassify()/isfinite() are folded to "always finite" under -ffast-math
 * (-ffinite-math-only), so test the exponent bits directly instead.
 *
 * @param x The value to check.
 * @return 1 if x is finite, 0 if it is NaN or infinite.
 */
static INLINE int cooling_value_is_finite(const gr_float x) {
#ifdef GRACKLE_FLOAT_4
  uint32_t bits;
  memcpy(&bits, &x, sizeof(bits));
  return ((bits >> 23) & 0xff) != 0xff;
#else
  uint64_t bits;
  memcpy(&bits, &x, sizeof(bits));
  return ((bits >> 52) & 0x7ff) != 0x7ff;
#endif
}

#ifdef CRACKLE_VALIDATE_DIRECT
/*! The crackle input state of the last field-path and direct-path calls of
 * this thread */
static __thread grackle_part_data crackle_dbg_gp_fields, crackle_dbg_gp_direct;

/*! If set, the field path runs on these inputs instead of its own */
static __thread int crackle_dbg_override = 0;

/*! Number of crackle inputs that differed at round-off level */
static long long crackle_dbg_n_roundoff = 0;

/**
 * @brief Give the field path the inputs of the direct path.
 *
 * Under -ffast-math the two copies of the input preparation may round
 * differently (e.g. the float sum in cooling_get_inverse_mu()), so the
 * inputs are compared to round-off separately and everything after them
 * bit for bit.
 */
static void cooling_crackle_override_fields(const grackle_part_data *gp,
                                            gr_float *sd, gr_float *iact) {
  for (int i = 0; i < 9; i++) {
    const gr_float v[9] = {gp->HI_density,   gp->HII_density,
                           gp->HeI_density,  gp->HeII_density,
                           gp->HeIII_density, gp->e_density,
                           gp->HM_density,   gp->H2I_density,
                           gp->H2II_density};
    sd[i] = v[i];
  }
  sd[12] = gp->density;
  sd[13] = gp->internal_energy;
  sd[14] = gp->temperature_floor;
  sd[15] = gp->specific_heating_rate;
  sd[16] = gp->x_velocity;
  sd[17] = gp->y_velocity;
  sd[18] = gp->z_velocity;
  sd[19] = gp->metal_density;
  sd[20] = gp->dust_density;
  sd[21] = gp->SNe_density;
  sd[22] = gp->isrf_habing;
  sd[23] = gp->H2_self_shielding_length;
  for (int i = 0; i < NUM_METAL_SPECIES_GRACKLE; i++) {
    sd[25 + i] = gp->gas_metalDensity[i];
    sd[36 + i] = gp->dust_metalDensity[i];
  }
  if (iact != NULL) {
    iact[0] = gp->RT_heating_rate;
    iact[1] = gp->RT_HI_ionization_rate;
    iact[2] = gp->RT_HeI_ionization_rate;
    iact[3] = gp->RT_HeII_ionization_rate;
    iact[4] = gp->RT_H2_dissociation_rate;
  }
}

/**
 * @brief Make the following cooling_copy_to_grackle() calls of this thread
 * use the inputs of its last cooling_copy_to_crackle() call (1) or not (0).
 */
void cooling_crackle_override_inputs(int on) { crackle_dbg_override = on; }

/**
 * @brief Compare the crackle inputs of the field and direct paths: stop if
 * any differs by more than round-off, count those that differ by round-off.
 */
void cooling_crackle_compare_inputs(long long id) {
  const grackle_part_data *g1 = &crackle_dbg_gp_fields;
  const grackle_part_data *g2 = &crackle_dbg_gp_direct;
  if (g1->grid_rank != g2->grid_rank || g1->grid_end != g2->grid_end ||
      g1->id != g2->id)
    error("crackle direct validation: id=%lld grid/id inputs differ", id);
  const double *d1 = &g1->density;
  const double *d2 = &g2->density;
  const size_t n =
      (sizeof(grackle_part_data) - offsetof(grackle_part_data, density)) /
      sizeof(double);
  const size_t k_shield = (offsetof(grackle_part_data, f_shield) -
                           offsetof(grackle_part_data, density)) /
                          sizeof(double);
  /* The RT rates are compared after the call (cooling_copy_to_grackle()
   * builds them for the instance with use_radiative_transfer) */
  const size_t k_rt0 = (offsetof(grackle_part_data, RT_heating_rate) -
                        offsetof(grackle_part_data, density)) /
                       sizeof(double);
  const size_t k_rt1 = (offsetof(grackle_part_data, RT_H2_dissociation_rate) -
                        offsetof(grackle_part_data, density)) /
                       sizeof(double);
  int nroundoff = 0;
  for (size_t k = 0; k < n; k++) {
    if (k == k_shield) continue; /* set by crackle */
    if (k >= k_rt0 && k <= k_rt1) continue;
    if (memcmp(&d1[k], &d2[k], sizeof(double)) == 0) continue;
    const double a = d1[k], b = d2[k];
    if (fabs(a - b) <= 1.e-6 * fmax(fabs(a), fabs(b))) {
      nroundoff++;
      continue;
    }
    error("crackle direct validation: id=%lld input at gp offset %zu "
          "(density %zu, internal_energy %zu, gas_metalDensity %zu) "
          "fields %.17g direct %.17g",
          id, offsetof(grackle_part_data, density) + k * sizeof(double),
          offsetof(grackle_part_data, density),
          offsetof(grackle_part_data, internal_energy),
          offsetof(grackle_part_data, gas_metalDensity), a, b);
  }
  if (nroundoff > 0) {
    const long long n_tot =
        atomic_add(&crackle_dbg_n_roundoff, (long long)nroundoff) + nroundoff;
    /* Report the 1st, 2nd, 4th, 8th... */
    if ((n_tot & (n_tot - 1)) == 0 || nroundoff > 1)
      message("crackle direct validation: inputs differing at round-off "
              "level so far: %lld (id=%lld, %d in this call)",
              n_tot, id, nroundoff);
  }
}
#endif

/**
 * @brief Fill the grackle data for a single particle.
 *
 * The grackle data point into species_densities and iact_rates, so these
 * must stay alive (in the caller's scope) until the grackle call is done.
 *
 * The density, internal energy and temperature floor depend on the case:
 * - mode 1 (cooling time): the particle's own density and energy, CMB
 *   floor, no heating and no UV background;
 * - gas outside the subgrid ISM (or with a negligible cold fraction): the
 *   particle's density and energy, the warm-ISM EoS temperature as the floor
 *   and the hydro du/dt as specific heating rate;
 * - subgrid ISM gas: the cold-cloud density and subgrid temperature, the CMB
 *   floor, and the hydro du/dt only if ism_adiabatic_heating_method == 1.
 * Errors if any value passed to grackle is not finite.
 *
 * @param data The grackle_field_data structure from grackle.
 * @param us The internal system of units.
 * @param cosmo The current cosmological model.
 * @param cooling The #cooling_function_data used in the run.
 * @param p The #part.
 * @param xp The #xpart.
 * @param dt The time-step of this particle.
 * @param T_warm The warm-ISM EoS temperature of the particle (K).
 * @param species_densities Buffer holding the values passed to grackle.
 * @param iact_rates The RT heating and ionisation rates passed to grackle.
 * @param mode The grackle driver mode (1 = cooling time).
 */
void cooling_copy_to_grackle(
    grackle_field_data *data, const struct unit_system *restrict us,
    const struct cosmology *restrict cosmo,
    const struct cooling_function_data *restrict cooling, const struct part *p,
    const struct xpart *xp, const double dt, const double T_warm,
    gr_float species_densities[N_SPECIES], gr_float *iact_rates, int mode) {

  /* Grid description for a single particle. Grackle only reads these, so
   * they are shared read-only (no per-call allocation). grid_end[0] = -1
   * signals crackle to turn off the UVB (used for cooling times). */
  static int grid_dimension[GRACKLE_RANK] = {GRACKLE_NPART, 1, 1};
  static int grid_start[GRACKLE_RANK] = {0, 0, 0};
  static int grid_end[GRACKLE_RANK] = {GRACKLE_NPART - 1, 0, 0};
  static int grid_end_no_uvb[GRACKLE_RANK] = {-1, 0, 0};

  /* Start from a clean struct so any field not set below is NULL rather
   * than garbage (callers' stack structs are uninitialized) */
  memset(data, 0, sizeof(*data));

  int i;
  /* set values */
  /* grid */
  data->grid_dx = 0.f;
  data->grid_rank = GRACKLE_RANK;
  data->grid_dimension = grid_dimension;
  data->grid_start = grid_start;
  data->grid_end = grid_end;

  /* get particle density, internal energy in
     physical coordinates (still code units) */
  double T_subgrid = cooling_get_subgrid_temperature(p, xp);
  /* mode 1 is cooling time, here we want non-subgrid values in all cases */
  if (mode == 1) {
    species_densities[12] = hydro_get_physical_density(p, cosmo);
    species_densities[13] = hydro_get_physical_internal_energy(p, xp, cosmo);
    species_densities[14] = cooling->T_CMB_0 * (1.f + cosmo->z);
    species_densities[15] = 0.f;
    data->grid_end = grid_end_no_uvb;  // signals crackle to turn off UVB
  }
  /* non-subgrid case, here we set the floor temperature
   * by the EoS (if applicable).  Note the cold fraction has a small
   * limit otherwise one can get underflows in crackle. */
  else if (p->cooling_data.subgrid_temp == 0. ||
           p->cooling_data.subgrid_fcold <= 1.e-6) {
    species_densities[12] = hydro_get_physical_density(p, cosmo);
    species_densities[13] = hydro_get_physical_internal_energy(p, xp, cosmo);
    species_densities[14] = T_warm;
    /* specific_heating_rate has to be in cgs units;
       no unit conversion done within grackle */
    species_densities[15] =
        hydro_get_physical_internal_energy_dt(p, cosmo) * cooling->dudt_units;
  }
  /* subgrid ISM case, use subgrid ISM values and set floor by T_CMB */
  else {
    /* Physical sub-grid (cold-cloud) density */
    species_densities[12] = cooling_get_subgrid_density(p, xp);
    /* Physical internal energy */
    species_densities[13] = cooling_convert_temp_to_u(
        T_subgrid, xp->cooling_data.e_frac, cooling, p, xp);
    /* CMB temp is floor*/
    species_densities[14] = cooling->T_CMB_0 * (1.f + cosmo->z);
    /* If tracking H2, turn off specific heating rate in ISM. */
    species_densities[15] = 0.f;
    if (cooling->ism_adiabatic_heating_method == 1) {
      species_densities[15] +=
          hydro_get_physical_internal_energy_dt(p, cosmo) * cooling->dudt_units;
    }
  }
  /* load into grackle structure */
  data->density = &species_densities[12];
  data->internal_energy = &species_densities[13];
  data->temperature_floor = &species_densities[14];
  data->specific_heating_rate = &species_densities[15];

  /* velocity (maybe not needed?) */
  species_densities[16] = xp->v_full[0] * cosmo->a_inv;
  species_densities[16] = (int)p->id;
  species_densities[17] = xp->v_full[1] * cosmo->a_inv;
  species_densities[18] = xp->v_full[2] * cosmo->a_inv;
  data->x_velocity = &species_densities[16];
  data->y_velocity = &species_densities[17];
  data->z_velocity = &species_densities[18];

  cooling_copy_to_grackle1(data, p, xp, species_densities[12],
                           species_densities);
  cooling_copy_to_grackle2(data, p, xp, cooling, dt, species_densities[12],
                           species_densities);
  cooling_copy_to_grackle3(data, p, xp, species_densities[12],
                           species_densities);

  /* In the subgrid ISM, RT photons are absorbed by the cold clouds only
   * (see rt_do_thermochemistry_with_subgrid), and the warm medium is taken
   * to be ionized: no RT heating and ionisation in the subgrid ISM */
  if (p->cooling_data.subgrid_temp > 0.f) {
    for (i = 0; i < 5; i++) iact_rates[i] = 0.f;
  }

  /*const double Ephot_cgs = p->rt_data.radiation[0].energy_density *
  p->geometry.volume * cooling->energy_units; if (Ephot_cgs > 1.e60) {
    message("IACT: z=%g id=%lld nH=%g T=%g tdec=%g %d SFR=%g RTheat=%g fluxE=%g
  Eph=%g Edens=%g V=%g h=%g Eunit=%g", cosmo->z, p->id, nH,
  cooling_convert_u_to_temp(species_densities[13], xp->cooling_data.e_frac,
  cooling, p), p->feedback_data.decoupling_delay_time, p->decoupled,
  p->sf_data.SFR, iact_rates[0], p->rt_data.flux[0].energy, Ephot_cgs,
  p->rt_data.radiation[0].energy_density, p->geometry.volume, p->h,
  cooling->energy_units);
  }*/

  /* RT heating and ionisation rates */
  data->RT_heating_rate = &iact_rates[0];
  data->RT_HI_ionization_rate = &iact_rates[1];
  data->RT_HeI_ionization_rate = &iact_rates[2];
  data->RT_HeII_ionization_rate = &iact_rates[3];
  data->RT_H2_dissociation_rate = &iact_rates[4];

  species_densities[19] =
      chemistry_get_total_metal_mass_fraction_for_cooling(p) *
      species_densities[12];
  data->metal_density = &species_densities[19];

  for (i = 0; i < N_SPECIES; i++) {
    if (!cooling_value_is_finite(species_densities[i])) {
      error(
          "Passing a non-finite value to grackle! "
          "i=%d / %d, species_densities[i]=%g\n",
          i, N_SPECIES, species_densities[i]);
    }
  }
#ifdef CRACKLE_VALIDATE_DIRECT
  /* The inputs crackle builds from these fields (with this run's cooling
   * chemistry flags, as cooling_copy_to_crackle() uses) */
  memset(&crackle_dbg_gp_fields, 0, sizeof(grackle_part_data));
  crackle_fields_to_part(data, &crackle_dbg_gp_fields,
                         (chemistry_data *)&cooling->chemistry);
  if (crackle_dbg_override)
    cooling_crackle_override_fields(&crackle_dbg_gp_direct, species_densities,
                                    iact_rates);
#endif
}

/**
 * @brief Copy the evolved species (and dust/metals) from the grackle data
 * back to the particle.
 *
 * @param data The grackle_field_data structure from grackle.
 * @param p The #part.
 * @param xp The #xpart.
 * @param cooling The #cooling_function_data used in the run.
 * @param rho The particle density passed to grackle (grackle units).
 */
void cooling_copy_from_grackle(
    grackle_field_data *data, struct part *p, struct xpart *xp,
    const struct cooling_function_data *restrict cooling, gr_float rho) {

  cooling_copy_from_grackle1(data, p, xp, rho);
  cooling_copy_from_grackle2(data, p, xp, cooling, rho);
  cooling_copy_from_grackle3(data, p, xp, rho);
}

/**
 * @brief free memory associated with grackle driver
 *
 * Nothing to do: the grid arrays set in cooling_copy_to_grackle are static.
 * Kept so that callers (e.g. KIARA RT) need not change.
 *
 * @param data The grackle_field_data structure from grackle.
 */
void cooling_grackle_free_data(grackle_field_data *data) {}


/*! What crackle needs at the current redshift (UV background rates and
 * crackle units), set by cooling_prepare_crackle_step() in cooling_update()
 * and read-only during the step. Kept out of #cooling_function_data so that
 * the restart files are unchanged; it is recomputed after a restart. */
static struct {
  photo_rate_storage uvb_rates;
  crackle_units cunits;
} cooling_crackle_step;

/**
 * @brief Everything crackle needs that depends only on the parameters, units
 * and redshift (temperature units, UV background rates, crackle units).
 *
 * Called in cooling_update (every step) and after initialising crackle.
 *
 * @param cooling The #cooling_function_data used in the run.
 */
void cooling_prepare_crackle_step(struct cooling_function_data *cooling) {
  if (crackle_prepare_step(&cooling->chemistry, &grackle_rates, &cooling->units,
                           &cooling_crackle_step.uvb_rates,
                           &cooling_crackle_step.cunits) == 0)
    error("Error in crackle_prepare_step.");
}

/**
 * @brief Fill the crackle state of a particle directly from the particle.
 *
 * This gives crackle exactly the inputs that cooling_copy_to_grackle() and
 * crackle's copy from the grackle_field_data give it, without the
 * intermediate arrays.
 *
 * @param gp The crackle particle state to fill.
 * @param cosmo The current cosmological model.
 * @param cooling The #cooling_function_data used in the run.
 * @param p The #part.
 * @param xp The #xpart.
 * @param dt The time-step of this particle.
 * @param T_warm The warm-ISM EoS temperature (used outside the subgrid ISM).
 * @param iact_rates The RT interaction rates (zeroed in the subgrid ISM).
 * @param mode 0 solve, 1 cooling time, 2 temperature.
 */
void cooling_copy_to_crackle(grackle_part_data *gp,
                             const struct cosmology *restrict cosmo,
                             const struct cooling_function_data *restrict
                                 cooling,
                             const struct part *p, const struct xpart *xp,
                             const double dt, const double T_warm,
                             gr_float *iact_rates, int mode) {

  const chemistry_data *chemistry = &cooling->chemistry;

  crackle_init_part(gp);

  /* A single particle; grid_end < 0 turns the UVB off (cooling times) */
  gp->grid_rank = GRACKLE_RANK;
  gp->grid_dimension = GRACKLE_NPART;
  gp->grid_start = 0;
  gp->grid_end = (mode == 1) ? -1 : GRACKLE_NPART - 1;

  /* Density, internal energy (physical), floor temperature and heating rate,
   * as in cooling_copy_to_grackle() */
  double T_subgrid = cooling_get_subgrid_temperature(p, xp);
  gr_float rho, u, T_floor, heating;
  if (mode == 1) {
    rho = hydro_get_physical_density(p, cosmo);
    u = hydro_get_physical_internal_energy(p, xp, cosmo);
    T_floor = cooling->T_CMB_0 * (1.f + cosmo->z);
    heating = 0.f;
  } else if (p->cooling_data.subgrid_temp == 0. ||
             p->cooling_data.subgrid_fcold <= 1.e-6) {
    rho = hydro_get_physical_density(p, cosmo);
    u = hydro_get_physical_internal_energy(p, xp, cosmo);
    T_floor = T_warm;
    heating =
        hydro_get_physical_internal_energy_dt(p, cosmo) * cooling->dudt_units;
  } else {
    rho = cooling_get_subgrid_density(p, xp);
    u = cooling_convert_temp_to_u(T_subgrid, xp->cooling_data.e_frac, cooling,
                                  p, xp);
    T_floor = cooling->T_CMB_0 * (1.f + cosmo->z);
    heating = 0.f;
    if (cooling->ism_adiabatic_heating_method == 1) {
      heating +=
          hydro_get_physical_internal_energy_dt(p, cosmo) * cooling->dudt_units;
    }
  }
  gp->density = rho;
  gp->internal_energy = u;
  gp->temperature_floor = T_floor;
  if (chemistry->use_specific_heating_rate)
    gp->specific_heating_rate = heating;
  else
    gp->specific_heating_rate = 0.;
  if (chemistry->use_volumetric_heating_rate)
    error("crackle: a volumetric heating rate is not supported");
  gp->volumetric_heating_rate = 0.;

  /* The x velocity slot holds the particle ID (for crackle's debugging) */
  const gr_float x_velocity = (int)p->id;
  gp->x_velocity = x_velocity;
  gp->y_velocity = xp->v_full[1] * cosmo->a_inv;
  gp->z_velocity = xp->v_full[2] * cosmo->a_inv;
  gp->id = (int)x_velocity;

  /* Primordial species */
  gp->HI_density = xp->cooling_data.HI_frac * rho;
  gp->HII_density = xp->cooling_data.HII_frac * rho;
  gp->HeI_density = xp->cooling_data.HeI_frac * rho;
  gp->HeII_density = xp->cooling_data.HeII_frac * rho;
  gp->HeIII_density = xp->cooling_data.HeIII_frac * rho;
  gp->e_density = xp->cooling_data.e_frac * rho;
  gp->HM_density = xp->cooling_data.HM_frac * rho;
  if (chemistry->primordial_chemistry >= 2) {
    gp->H2I_density = xp->cooling_data.H2I_frac * rho;
    gp->H2II_density = xp->cooling_data.H2II_frac * rho;
  } else {
    gp->H2I_density = 0.;
    gp->H2II_density = 0.;
  }
  if (chemistry->primordial_chemistry >= 3)
    error("crackle does not work with primordial_chemistry >= 3");
  gp->DI_density = 0.;
  gp->DII_density = 0.;
  gp->HDI_density = 0.;

  /* In the subgrid ISM, RT photons are absorbed by the cold clouds only
   * (see rt_do_thermochemistry_with_subgrid), and the warm medium is taken
   * to be ionized: no RT heating and ionisation in the subgrid ISM */
  if (p->cooling_data.subgrid_temp > 0.f) {
    for (int i = 0; i < 5; i++) iact_rates[i] = 0.f;
  }
  /* Passed whatever this instance's use_radiative_transfer, as in the field
   * path: the KIARA-RT chemistry instance uses them, and the cooling one may
   * have them off */
  gp->RT_heating_rate = iact_rates[0];
  gp->RT_HI_ionization_rate = iact_rates[1];
  gp->RT_HeI_ionization_rate = iact_rates[2];
  gp->RT_HeII_ionization_rate = iact_rates[3];
  gp->RT_H2_dissociation_rate = iact_rates[4];

  gp->metal_density =
      chemistry_get_total_metal_mass_fraction_for_cooling(p) * rho;

  /* Dust model, ISRF and H2 shielding length (as cooling_copy_to_grackle2) */
  gr_float isrf = 0., L_shield = 0., dust_density = 0., SNe_density = 0.;
  /* By crackle slot */
  gr_float gas_metals[NUM_METAL_SPECIES_GRACKLE] = {0.};
  gr_float dust_metals[NUM_METAL_SPECIES_GRACKLE] = {0.};
  if (cooling->use_grackle_dust_evol == 1) {
    dust_density = p->cooling_data.dust_mass / p->mass * rho;
    SNe_density = p->feedback_data.SNe_ThisTimeStep * dt / p->mass * rho;
    if (p->cooling_data.subgrid_temp == 0.f) {
      isrf = -1.f;
    } else {
      isrf = p->cooling_data.G0;
    }
    const double rho_grad_norm2 = p->rho_gradient[0] * p->rho_gradient[0] +
                                  p->rho_gradient[1] * p->rho_gradient[1] +
                                  p->rho_gradient[2] * p->rho_gradient[2];
    const double rho_grad_norm_inv =
        (rho_grad_norm2 > 0.) ? 1. / sqrt(rho_grad_norm2) : 0.;
    const double rho_com = hydro_get_comoving_density(p);
    double L_eff_com = rho_com * rho_grad_norm_inv;
    const double L_eff_com_max = kernel_gamma * p->h;
    const double L_eff_com_min = MIN_SHIELD_H_FRAC * p->h;
    L_eff_com = fmin(L_eff_com, L_eff_com_max);
    L_eff_com = fmax(L_eff_com, L_eff_com_min);
    L_shield = L_eff_com * cooling->units.a_value;
    for (int i = 0; i < chemistry_element_count; i++) {
      const int slot = cooling_crackle_slot(i);
      if (slot >= 0)
        gas_metals[slot] =
            max(p->chemistry_data.metal_mass_fraction[i] * rho, 0.);
    }
    for (int d = 0; d < dust_element_count; d++) {
      dust_metals[cooling_crackle_slot(chemistry_element_of_dust(d))] =
          max(p->cooling_data.dust_mass_fraction[d] * dust_density, 0);
    }
  }
  if (chemistry->H2_custom_shielding == 2) {
    if (cooling->use_grackle_dust_evol != 1)
      error("crackle: H2_custom_shielding = 2 needs the dust model");
    gp->H2_self_shielding_length = L_shield;
    gp->H2_custom_shielding_factor = 0.; /* unused */
  } else if (chemistry->H2_custom_shielding == 1) {
    error("crackle: H2_custom_shielding = 1 is not supported");
  } else {
    gp->H2_self_shielding_length = 0.;
    gp->H2_custom_shielding_factor = 0.;
  }
  if (chemistry->use_isrf_field)
    gp->isrf_habing = isrf;
  else
    gp->isrf_habing = 0.;
  if (chemistry->use_dust_evol) {
    gp->dust_density = dust_density;
    for (int i = 0; i < NUM_METAL_SPECIES_GRACKLE; i++) {
      gp->gas_metalDensity[i] = gas_metals[i];
      gp->dust_metalDensity[i] = dust_metals[i];
    }
    gp->SNe_density = SNe_density;
  } else {
    gp->dust_density = 0.;
    gp->SNe_density = 0.;
    for (int i = 0; i < NUM_METAL_SPECIES_GRACKLE; i++) {
      gp->gas_metalDensity[i] = 0.;
      gp->dust_metalDensity[i] = 0.;
    }
  }

  /* Never pass a non-finite value to crackle */
  const gr_float check[] = {rho,
                            u,
                            T_floor,
                            heating,
                            gp->y_velocity,
                            gp->z_velocity,
                            gp->HI_density,
                            gp->HII_density,
                            gp->HeI_density,
                            gp->HeII_density,
                            gp->HeIII_density,
                            gp->e_density,
                            gp->HM_density,
                            gp->H2I_density,
                            gp->H2II_density,
                            gp->metal_density,
                            dust_density,
                            SNe_density,
                            isrf,
                            L_shield};
  for (size_t i = 0; i < sizeof(check) / sizeof(check[0]); i++) {
    if (!cooling_value_is_finite(check[i]))
      error("Passing a non-finite value to crackle! id=%lld i=%zu value=%g",
            p->id, i, check[i]);
  }
  for (int i = 0; i < NUM_METAL_SPECIES_GRACKLE; i++) {
    if (!cooling_value_is_finite(gas_metals[i]) ||
        !cooling_value_is_finite(dust_metals[i]))
      error("Passing a non-finite metal density to crackle! id=%lld i=%d",
            p->id, i);
  }
#ifdef CRACKLE_VALIDATE_DIRECT
  memcpy(&crackle_dbg_gp_direct, gp, sizeof(grackle_part_data));
#endif
}

/**
 * @brief Copy the evolved species (and dust/metals) from the crackle state
 * back to the particle.
 *
 * Gives the particle exactly what crackle's copy to the grackle_field_data
 * and cooling_copy_from_grackle() give it.
 *
 * @param gp The crackle particle state after the solve.
 * @param p The #part.
 * @param xp The #xpart.
 * @param cooling The #cooling_function_data used in the run.
 */
void cooling_copy_from_crackle(const grackle_part_data *gp, struct part *p,
                               struct xpart *xp,
                               const struct cooling_function_data *restrict
                                   cooling) {

  const gr_float rho = gp->density;
  double rhoinv = 1.f / rho;

  /* Primordial species */
  xp->cooling_data.HI_frac = gp->HI_density * rhoinv;
  xp->cooling_data.HII_frac = gp->HII_density * rhoinv;
  xp->cooling_data.HeI_frac = gp->HeI_density * rhoinv;
  xp->cooling_data.HeII_frac = gp->HeII_density * rhoinv;
  xp->cooling_data.HeIII_frac = gp->HeIII_density * rhoinv;
  xp->cooling_data.e_frac = gp->e_density * rhoinv;
  xp->cooling_data.HM_frac = gp->HM_density * rhoinv;
  xp->cooling_data.H2I_frac = gp->H2I_density * rhoinv;
  xp->cooling_data.H2II_frac = gp->H2II_density * rhoinv;

  /* Dust model */
  if (cooling->use_grackle_dust_evol == 1) {
    p->chemistry_data.metal_mass_fraction_total = gp->metal_density * rhoinv;
    for (int i = 0; i < chemistry_element_count; i++) {
      const int slot = cooling_crackle_slot(i);
      if (slot >= 0)
        p->chemistry_data.metal_mass_fraction[i] =
            gp->gas_metalDensity[slot] * rhoinv;
    }

    /* Dust temperature used by crackle (unchanged if not computed) */
    if (gp->tdust > 0.) p->cooling_data.dust_temperature = gp->tdust;

    p->cooling_data.dust_mass = gp->dust_density * p->mass * rhoinv;
    if (p->cooling_data.dust_mass > 0.5 * p->mass) {
      warning("DUST > METALS Mg=%g Zg=%g mdust=%g mmet=%g\n", p->mass,
              chemistry_get_total_metal_mass_fraction_for_cooling(p),
              p->cooling_data.dust_mass,
              p->mass * chemistry_get_total_metal_mass_fraction_for_cooling(p));
    }
    if (gp->dust_density > 0.f) {
      rhoinv = 1. / gp->dust_density;
      /* no Helium metal density */
      assert(gp->dust_metalDensity[0] == 0.);
      for (int d = 0; d < dust_element_count; d++) {
        p->cooling_data.dust_mass_fraction[d] =
            gp->dust_metalDensity[cooling_crackle_slot(
                chemistry_element_of_dust(d))] *
            rhoinv;
      }
    } else {
      for (int d = 0; d < dust_element_count; d++) {
        p->cooling_data.dust_mass_fraction[d] = 0.f;
      }
    }
  }
}

/**
 * @brief Renomalize individual H and He species to account for
 * changes made to total H,He in other parts of code
 * (assuming all species are scaled equally).
 *
 * @param p Pointer to the particle data.
 * @param xp Pointer to the particle extra data
 *
 */
__attribute__((always_inline)) INLINE void cooling_normalize_primordial_species(
    struct part *restrict p, struct xpart *restrict xp) {

  /* Normalize hydrogen species */
  float X_H = xp->cooling_data.HI_frac + xp->cooling_data.HII_frac;
#if COOLING_GRACKLE_MODE >= 2
  X_H += xp->cooling_data.H2I_frac + xp->cooling_data.H2II_frac +
         xp->cooling_data.HM_frac;
#endif
  const float new_H_frac =
      p->chemistry_data.metal_mass_fraction[chemistry_element_H];
  const float H_ratio = (X_H > 0.f) ? new_H_frac / X_H : 1.f;

  xp->cooling_data.HI_frac *= H_ratio;
  xp->cooling_data.HII_frac *= H_ratio;
#if COOLING_GRACKLE_MODE >= 2
  xp->cooling_data.H2I_frac *= H_ratio;
  xp->cooling_data.H2II_frac *= H_ratio;
  xp->cooling_data.HM_frac *= H_ratio;
#endif

  /* Normalize helium species */
  float X_He = xp->cooling_data.HeI_frac + xp->cooling_data.HeII_frac +
               xp->cooling_data.HeIII_frac;
  const float new_He_frac =
      p->chemistry_data.metal_mass_fraction[chemistry_element_He];
  const float He_ratio = (X_He > 0.f) ? new_He_frac / X_He : 1.f;

  xp->cooling_data.HeI_frac *= He_ratio;
  xp->cooling_data.HeII_frac *= He_ratio;
  xp->cooling_data.HeIII_frac *= He_ratio;

  /* Recompute electron density (grackle convention: n_e * m_H / rho) */
  xp->cooling_data.e_frac = xp->cooling_data.HII_frac +
                            0.25f * xp->cooling_data.HeII_frac +
                            0.5f * xp->cooling_data.HeIII_frac;
#if COOLING_GRACKLE_MODE >= 2
  xp->cooling_data.e_frac +=
      0.5f * xp->cooling_data.H2II_frac - xp->cooling_data.HM_frac;
#endif
#if COOLING_GRACKLE_MODE >= 3
  xp->cooling_data.e_frac += 0.5f * xp->cooling_data.DII_frac;
#endif
}

/**
 * @brief Compute the energy of a particle after dt and update the particle
 * chemistry data
 *
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param cosmo The #cosmology.
 * @param hydro_props The #hydro_props.
 * @param cooling The #cooling_function_data used in the run.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the particle extra data
 * @param iact_rates The RT heating and ionisation rates passed to grackle.
 * @param dt The time-step of this particle.
 * @param T_warm The warm-ISM EoS temperature of the particle (K).
 * @param mode 0=energy, 1=cooling time, 2=temperature, 3=pressure, 4=gamma
 *
 * @return desired quantity based on mode
 */
static gr_float cooling_grackle_driver_fields(
    const struct phys_const *restrict phys_const,
    const struct unit_system *restrict us,
    const struct cosmology *restrict cosmo,
    const struct hydro_props *hydro_props,
    const struct cooling_function_data *restrict cooling,
    struct part *restrict p, struct xpart *restrict xp, gr_float *iact_rates,
    double dt, double T_warm, int mode) {

  /* set current units for conversion to physical quantities */
  code_units units = cooling->units;

  /* initialize data to send to grackle */
  gr_float species_densities[N_SPECIES] = {0};
  grackle_field_data data;
  // cooling_grackle_malloc_fields(&data, 1, cooling->chemistry.use_dust_evol);

  /* Make an initial guess of H2 fractions, to help grackle converge
  const float H2_frac = xp->cooling_data.H2I_frac + xp->cooling_data.H2II_frac;
  if (p->cooling_data.subgrid_temp > 0.f && H2_frac > 0.f && H2_frac < 1.f &&
  p->chemistry_data.metal_mass_fraction_total > 1.e-6) {
    //float H2_frac_KMT = cooling_KMT_H2_fraction(p, cosmo, cooling);
    //if (H2_frac_KMT < 0.f) H2_frac_KMT = 0.f;
    //if (H2_frac_KMT > 1.f) H2_frac_KMT = 1.f;
    const float H_frac = H2_frac + xp->cooling_data.HI_frac +
  xp->cooling_data.HII_frac + xp->cooling_data.HM_frac; float P_over_k_cgs =
  hydro_get_physical_pressure(p, cosmo) * units_cgs_conversion_factor(us,
  UNIT_CONV_PRESSURE) / 1.38e-16; float Rmol =
  pow(p->chemistry_data.metal_mass_fraction_total * P_over_k_cgs / (0.0134
  * 4.3e4), 0.92); float H2_frac_BR = Rmol / (1.f + Rmol) * H_frac; if (p->id %
  10 == 0) message("z=%g T=%g nH=%g P/k=%g Rmol=%g H2old=%g H2BR=%g H=%g",
  cosmo->z, p->cooling_data.subgrid_temp, p->cooling_data.subgrid_dens *
  cooling->units.density_units * 0.75 / 1.673e-24, P_over_k_cgs, Rmol, H2_frac,
  H2_frac_BR, H_frac); const float H2_guess = H2_frac_BR; const float H2_factor
  = H2_guess / H2_frac; assert(H2_factor > 0.f); assert(non_H2_factor > 0.f);
    xp->cooling_data.H2I_frac *= H2_factor;
    xp->cooling_data.H2II_frac *= H2_factor;
    xp->cooling_data.HI_frac *= non_H2_factor;
    xp->cooling_data.HII_frac *= non_H2_factor;
    xp->cooling_data.HM_frac *= non_H2_factor;
  }*/

  /* load particle information from particle to grackle data */
  cooling_copy_to_grackle(&data, us, cosmo, cooling, p, xp, dt, T_warm,
                          species_densities, iact_rates, mode);


  /* Run Grackle in desired mode */
  gr_float return_value = 0.f;

#ifdef CRACKLE_VALIDATE
  /* Save the inputs, to run the external library on them after each call */
  struct crackle_validate_state crv;
  crackle_validate_save(&crv, &units, &data, species_densities, N_SPECIES,
                        iact_rates, 5);
#endif

  switch (mode) {
    case 0:
      /* solve chemistry, advance thermal energy by dt */
      if (solve_chemistry(&units, &data, dt) == 0) {
        error("Error in Grackle solve_chemistry.");
      }
#ifdef CRACKLE_VALIDATE
      crackle_validate_check(&crv, crackle_validate_global, 0, dt,
                             species_densities, iact_rates, 0., p->id);
#endif
      // if (solve_chemistry(&units, &data, -dt) == 0) {
      //   error("Error in Crackle solve_chemistry.");
      // }
      /* copy from grackle data to particle */
      cooling_copy_from_grackle(&data, p, xp, cooling, species_densities[12]);
      return_value = data.internal_energy[0];
      /* The dust temperature is returned by crackle itself (with the dust
       * model on); grackle's calculate_dust_temperature() uses a different
       * dust opacity normalisation, so is not used. */
      break;

    case 1:
      /* compute cooling time */
      if (calculate_cooling_time(&units, &data, &return_value) == 0) {
        error("Error in Grackle calculate_cooling_time.");
      }
#ifdef CRACKLE_VALIDATE
      crackle_validate_check(&crv, crackle_validate_global, 1, dt,
                             species_densities, iact_rates, return_value,
                             p->id);
#endif
      break;

    case 2:
      /* compute temperature */
      if (calculate_temperature(&units, &data, &return_value) == 0) {
        error("Error in Grackle calculate_temperature.");
      }
#ifdef CRACKLE_VALIDATE
      crackle_validate_check(&crv, crackle_validate_global, 2, dt,
                             species_densities, iact_rates, return_value,
                             p->id);
#endif
      break;

    case 3:
      /* compute pressure */
      if (calculate_pressure(&units, &data, &return_value) == 0) {
        error("Error in Grackle calculate_pressure.");
      }
#ifdef CRACKLE_VALIDATE
      crackle_validate_check(&crv, crackle_validate_global, 3, dt,
                             species_densities, iact_rates, return_value,
                             p->id);
#endif
      break;
    case 4:
      /* compute gamma */
      if (calculate_gamma(&units, &data, &return_value) == 0) {
        error("Error in Grackle calculate_gamma.");
      }
#ifdef CRACKLE_VALIDATE
      crackle_validate_check(&crv, crackle_validate_global, 4, dt,
                             species_densities, iact_rates, return_value,
                             p->id);
#endif
      break;
  }

  cooling_grackle_free_data(&data);

  return return_value;
}

/**
 * @brief Run crackle directly on the particle (modes 0-2), without the
 * grackle_field_data arrays.
 *
 * @param cosmo The #cosmology.
 * @param cooling The #cooling_function_data used in the run.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the particle extra data
 * @param iact_rates The RT heating and ionisation rates passed to crackle.
 * @param dt The time-step of this particle.
 * @param T_warm The warm-ISM EoS temperature of the particle (K).
 * @param mode 0=energy, 1=cooling time, 2=temperature
 *
 * @return desired quantity based on mode
 */
static gr_float cooling_crackle_driver(
    const struct cosmology *restrict cosmo,
    const struct cooling_function_data *restrict cooling,
    struct part *restrict p, struct xpart *restrict xp, gr_float *iact_rates,
    double dt, double T_warm, int mode) {

  /* The units with the temperature units of the current step */
  code_units units = cooling->units;
  chemistry_data *chemistry = (chemistry_data *)&cooling->chemistry;

  grackle_part_data gp;
  cooling_copy_to_crackle(&gp, cosmo, cooling, p, xp, dt, T_warm, iact_rates,
                          mode);

  gr_float return_value = 0.f;
  switch (mode) {
    case 0:
      /* solve chemistry, advance thermal energy by dt */
      if (crackle_solve_part(&gp, chemistry, &grackle_rates,
                             &cooling_crackle_step.uvb_rates, &units,
                             &cooling_crackle_step.cunits, fabs(dt)) == 0) {
        error("Error in crackle_solve_part.");
      }
      cooling_copy_from_crackle(&gp, p, xp, cooling);
      return_value = gp.internal_energy;
      break;

    case 1:
      /* compute cooling time */
      crackle_cooling_time_part(&gp, chemistry, &grackle_rates,
                                &cooling_crackle_step.uvb_rates, &units,
                                &cooling_crackle_step.cunits, &return_value);
      break;

    case 2:
      /* compute temperature */
      crackle_temperature_part(&gp, chemistry, &units, &return_value);
      break;

    default:
      error("crackle driver: unknown mode %d", mode);
  }

  return return_value;
}

#ifdef CRACKLE_VALIDATE_DIRECT
/**
 * @brief Stop if the direct crackle path did not reproduce the field-array
 * path bit for bit.
 */
void cooling_crackle_validate_direct(
    const struct part *p, const struct part *p_ref, const struct xpart *xp,
    const struct xpart *xp_ref, const gr_float *iact, const gr_float *iact_ref,
    const gr_float result, const gr_float result_ref, const int mode) {

  const int diff_p = memcmp(p, p_ref, sizeof(struct part));
  const int diff_xp = memcmp(xp, xp_ref, sizeof(struct xpart));
  const int diff_iact =
      (iact != NULL) ? memcmp(iact, iact_ref, 5 * sizeof(gr_float)) : 0;
  const int diff_r = memcmp(&result, &result_ref, sizeof(gr_float));
  if (!diff_p && !diff_xp && !diff_iact && !diff_r) return;

  const unsigned char *a = (const unsigned char *)p;
  const unsigned char *b = (const unsigned char *)p_ref;
  for (size_t k = 0; k < sizeof(struct part); k++)
    if (a[k] != b[k]) {
      message("crackle direct validation: part byte %zu differs "
              "(cooling_data at %zu, chemistry_data at %zu)",
              k, offsetof(struct part, cooling_data),
              offsetof(struct part, chemistry_data));
      break;
    }
  a = (const unsigned char *)xp;
  b = (const unsigned char *)xp_ref;
  for (size_t k = 0; k < sizeof(struct xpart); k++)
    if (a[k] != b[k]) {
      message("crackle direct validation: xpart byte %zu differs "
              "(cooling_data at %zu)",
              k, offsetof(struct xpart, cooling_data));
      break;
    }
  message("crackle direct validation: id=%lld mode=%d result %.17g ref %.17g; "
          "HI %.17g/%.17g e %.17g/%.17g H2I %.17g/%.17g u_dt %.17g/%.17g "
          "Tdust %.17g/%.17g mdust %.17g/%.17g",
          p->id, mode, result, result_ref, xp->cooling_data.HI_frac,
          xp_ref->cooling_data.HI_frac, xp->cooling_data.e_frac,
          xp_ref->cooling_data.e_frac, xp->cooling_data.H2I_frac,
          xp_ref->cooling_data.H2I_frac, p->u_dt, p_ref->u_dt,
          p->cooling_data.dust_temperature,
          p_ref->cooling_data.dust_temperature, p->cooling_data.dust_mass,
          p_ref->cooling_data.dust_mass);
  error("crackle direct validation: the direct path differs from the field "
        "path (part %d, xpart %d, iact %d, result %d)",
        diff_p != 0, diff_xp != 0, diff_iact != 0, diff_r != 0);
}
#endif

/**
 * @brief Compute the energy of a particle after dt and update the particle
 * chemistry data, or another quantity depending on the mode.
 *
 * Modes 0-2 run crackle directly on the particle; modes 3 and 4 (pressure,
 * gamma; output only) go through Grackle's field arrays. With
 * CRACKLE_VALIDATE (comparison with an external library) everything goes
 * through the field arrays; with CRACKLE_VALIDATE_DIRECT both paths are run
 * and must agree bit for bit.
 *
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param cosmo The #cosmology.
 * @param hydro_props The #hydro_props.
 * @param cooling The #cooling_function_data used in the run.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the particle extra data
 * @param iact_rates The RT heating and ionisation rates passed to grackle.
 * @param dt The time-step of this particle.
 * @param T_warm The warm-ISM EoS temperature of the particle (K).
 * @param mode 0=energy, 1=cooling time, 2=temperature, 3=pressure, 4=gamma
 *
 * @return desired quantity based on mode
 */
gr_float cooling_grackle_driver(
    const struct phys_const *restrict phys_const,
    const struct unit_system *restrict us,
    const struct cosmology *restrict cosmo,
    const struct hydro_props *hydro_props,
    const struct cooling_function_data *restrict cooling,
    struct part *restrict p, struct xpart *restrict xp, gr_float *iact_rates,
    double dt, double T_warm, int mode) {

  /* Renomalize H,He species to account for any changes due to chemistry,
   * stellar evol, etc */
  cooling_normalize_primordial_species(p, xp);

#if defined(CRACKLE_VALIDATE)
  return cooling_grackle_driver_fields(phys_const, us, cosmo, hydro_props,
                                       cooling, p, xp, iact_rates, dt, T_warm,
                                       mode);
#else
  if (mode > 2)
    return cooling_grackle_driver_fields(phys_const, us, cosmo, hydro_props,
                                         cooling, p, xp, iact_rates, dt,
                                         T_warm, mode);
#ifdef CRACKLE_VALIDATE_DIRECT
  /* Both paths start from a snapshot of the particle, as other tasks may
   * write parts of it (e.g. black_holes_data) while this one runs; the real
   * update is then done by the direct path on the particle itself */
  const struct part p_in = *p;
  const struct xpart xp_in = *xp;
  gr_float iact_in[5] = {0., 0., 0., 0., 0.};
  if (iact_rates != NULL) memcpy(iact_in, iact_rates, sizeof(iact_in));

  struct part p_dir = p_in;
  struct xpart xp_dir = xp_in;
  gr_float iact_dir[5];
  memcpy(iact_dir, iact_in, sizeof(iact_dir));
  const gr_float result_dir = cooling_crackle_driver(
      cosmo, cooling, &p_dir, &xp_dir, (iact_rates != NULL) ? iact_dir : NULL,
      dt, T_warm, mode);

  /* The field path, on the inputs of the direct path */
  struct part p_ref = p_in;
  struct xpart xp_ref = xp_in;
  gr_float iact_ref[5];
  memcpy(iact_ref, iact_in, sizeof(iact_ref));
  cooling_crackle_override_inputs(1);
  const gr_float result_ref = cooling_grackle_driver_fields(
      phys_const, us, cosmo, hydro_props, cooling, &p_ref, &xp_ref,
      (iact_rates != NULL) ? iact_ref : NULL, dt, T_warm, mode);
  cooling_crackle_override_inputs(0);

  cooling_crackle_compare_inputs(p->id);
  cooling_crackle_validate_direct(&p_dir, &p_ref, &xp_dir, &xp_ref,
                                  (iact_rates != NULL) ? iact_dir : NULL,
                                  iact_ref, result_dir, result_ref, mode);

  const gr_float result = cooling_crackle_driver(cosmo, cooling, p, xp,
                                                 iact_rates, dt, T_warm, mode);
  return result;
#else
  return cooling_crackle_driver(cosmo, cooling, p, xp, iact_rates, dt, T_warm,
                                mode);
#endif
#endif
}

/**
 * @brief Compute the cooling time. Optionally uses input values
 * for rho and u, but leaves all other properties of p the same.
 *
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param hydro_properties The #hydro_props.
 * @param cosmo The #cosmology.
 * @param cooling The #cooling_function_data used in the run.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the particle extra data.
 * @param rhocool Comoving density used in the tcool calculation (if > 0).
 * @param ucool Comoving internal energy used in the tcool calculation (if >
 * 0).
 *
 * @return The cooling time (internal units; positive if heating, negative if
 * cooling).
 */
gr_float cooling_time(const struct phys_const *restrict phys_const,
                      const struct unit_system *restrict us,
                      const struct hydro_props *hydro_properties,
                      const struct cosmology *restrict cosmo,
                      const struct cooling_function_data *restrict cooling,
                      const struct part *restrict p, struct xpart *restrict xp,
                      const float rhocool, const float ucool) {

  /* Removes const in declaration*/
  struct part p_temp = *p;

  if (rhocool > 0.f) p_temp.rho = rhocool;
  if (ucool > 0.f) p_temp.u = ucool;

  gr_float iact_rates[5] = {0., 0., 0., 0., 0.};

  gr_float cooling_time =
      cooling_grackle_driver(phys_const, us, cosmo, hydro_properties, cooling,
                             &p_temp, xp, iact_rates, 0., 0., 1);

  return cooling_time;
}

/**
 * @brief Compute the temperature of a #part based on the cooling function.
 *
 * @param phys_const #phys_const data structure.
 * @param hydro_properties The properties of the hydro scheme.
 * @param us The internal system of units.
 * @param cosmo #cosmology data structure.
 * @param cooling #cooling_function_data struct.
 * @param p #part data.
 * @param xp Pointer to the #xpart data.
 * @return The temperature of the gas (internal units).
 */
float cooling_get_temperature(
    const struct phys_const *restrict phys_const,
    const struct hydro_props *restrict hydro_properties,
    const struct unit_system *restrict us,
    const struct cosmology *restrict cosmo,
    const struct cooling_function_data *restrict cooling,
    const struct part *restrict p, const struct xpart *restrict xp) {

  /* Remove const in declaration*/
  struct part p_temp = *p;
  struct xpart xp_temp = *xp;

  gr_float iact_rates[5] = {0., 0., 0., 0., 0.};

  const float temperature =
      cooling_grackle_driver(phys_const, us, cosmo, hydro_properties, cooling,
                             &p_temp, &xp_temp, iact_rates, 0., 0., 2);

  return temperature;
}

/**
 * @brief Do dust sputtering and return metals back to gas phase
 *
 * Above 1e4 K, dust is destroyed on the thermal sputtering time-scale of
 * Tsai & Mathews (1995) and its element masses are returned to the gas.
 *
 * @param us The internal system of units.
 * @param cosmo The current cosmological model.
 * @param cooling The #cooling_function_data used in the run.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the #xpart data.
 * @param dt The time-step of this particle.
 */
__attribute__((always_inline)) INLINE void cooling_sputter_dust(
    const struct unit_system *restrict us,
    const struct cosmology *restrict cosmo,
    const struct cooling_function_data *restrict cooling,
    struct part *restrict p, struct xpart *restrict xp, const double dt) {

  /* Do dust destruction in stream particle */
  if (cooling->use_grackle_dust_evol && p->cooling_data.dust_mass > 0.f) {
    const float u_phys = hydro_get_physical_internal_energy(p, xp, cosmo);
    const float Tstream = cooling_convert_u_to_temp(
        u_phys, xp->cooling_data.e_frac, cooling, p, xp);
    const double Tstream_K =
        Tstream * units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE);

    /* Sputtering negligible at low-T */
    if (Tstream_K > 1.e4) {
      const double rho_cgs = hydro_get_physical_density(p, cosmo) *
                             units_cgs_conversion_factor(us, UNIT_CONV_DENSITY);

      /* sputtering timescale, Tsai & Mathews (1995) */
      const double tsp = 1.7e8 * 3.15569251e7 /
                         units_cgs_conversion_factor(us, UNIT_CONV_TIME) *
                         (cooling->dust_grainsize / 0.1) * (1.e-27 / rho_cgs) *
                         (pow(2.e6 / Tstream, 2.5) + 1.0);

      const float dust_mass_old = p->cooling_data.dust_mass;

      /* Update dust mass but limit the destruction */
      p->cooling_data.dust_mass -=
          dust_mass_old * (1.f - exp(-3.f * min(dt / tsp, 5.f)));
      if (p->cooling_data.dust_mass < 0.5f * dust_mass_old) {
        p->cooling_data.dust_mass = 0.5f * dust_mass_old;
      }

#ifdef FIREHOSE_DEBUG_CHECKS
      message(
          "FIREHOSE_SPUT: id=%lld mdust=%g mdustnew=%g T=%g rho=%g "
          "tsp_dt_ratio=%g",
          p->id, dust_mass_old, p->cooling_data.dust_mass, Tstream,
          rho_cgs / 1.673e-24, /* units of mp */
          tsp / dt);
#endif
      /* Dust mass returned to the gas phase */
      const float dust_mass_lost = dust_mass_old - p->cooling_data.dust_mass;

      /* All of the dust is metals */
      p->chemistry_data.metal_mass_fraction_total +=
          dust_mass_lost / hydro_get_mass(p);

      for (int elem = 0; elem < chemistry_element_count; ++elem) {
        if (elem == chemistry_element_H || elem == chemistry_element_He)
          continue;

        /* dust_mass_fraction is the fraction of the dust mass in each
         * element; sputtering does not change the grain composition, so
         * it is left unchanged and each element loses the same fraction. */
        const float Z_elem_old = p->chemistry_data.metal_mass_fraction[elem];
        const float elem_mass_old = Z_elem_old * hydro_get_mass(p);

        /* This is the positive amount of metal mass to add since we
         * are losing dust mass when sputtering. */
        const int d = dust_element_of_chemistry(elem);
        const float delta_metal_mass_elem =
            (d >= 0) ? p->cooling_data.dust_mass_fraction[d] * dust_mass_lost
                     : 0.f;
        const float elem_mass_new = elem_mass_old + delta_metal_mass_elem;
        const float Z_elem_new = elem_mass_new / hydro_get_mass(p);

        p->chemistry_data.metal_mass_fraction[elem] = Z_elem_new;
      }

      /* Make sure that X + Y + Z + D = 1 (Z is gas-phase only) */
      const float Y_He =
          p->chemistry_data.metal_mass_fraction[chemistry_element_He];
      p->chemistry_data.metal_mass_fraction[chemistry_element_H] =
          1.f - Y_He - p->chemistry_data.metal_mass_fraction_total -
          p->cooling_data.dust_mass / hydro_get_mass(p);

      /* Make sure H fraction does not go out of bounds */
      if (p->chemistry_data.metal_mass_fraction[chemistry_element_H] > 1.f ||
          p->chemistry_data.metal_mass_fraction[chemistry_element_H] < 0.f) {
        for (int i = chemistry_element_H; i < chemistry_element_count; i++) {
          warning("\telem[%d] is %g", i,
                  p->chemistry_data.metal_mass_fraction[i]);
        }

        error(
            "Hydrogen fraction exeeds unity or is negative for"
            " particle id=%lld due to dust sputtering",
            p->id);
      }
    }
  }
}

/**
 * @brief Compute particle quantities for the firehose model
 *
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param cosmo The current cosmological model.
 * @param hydro_props The #hydro_props.
 * @param cooling The #cooling_function_data used in the run.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the #xpart data.
 * @param dt The time-step of this particle.
 */
__attribute__((always_inline)) INLINE void firehose_cooling_and_dust(
    const struct phys_const *restrict phys_const,
    const struct unit_system *restrict us,
    const struct cosmology *restrict cosmo,
    const struct hydro_props *restrict hydro_props,
    const struct cooling_function_data *restrict cooling,
    struct part *restrict p, struct xpart *restrict xp, const double dt) {

  /* Initialize cooling time */
  p->cooling_data.mixing_layer_cool_time = 0.f;

  /* Only firehose stream particles need a cooling time during the run (the
   * mixing-layer time, used in chemistry_iact.h). For all other particles it
   * is only an output quantity, computed at snapshot time (cooling_io.h), so
   * skip the extra grackle call here. */
  if (p->chemistry_data.radius_stream <= 0.f ||
      p->chemistry_data.rho_ambient <= 0.f) {
    return;
  }

  const float u = hydro_get_comoving_internal_energy(p, xp);

  /* It's a firehose particles, so compute the cooling rate
   * in the mixing layer. The ambient density is already capped at the
   * firehose maximum; cap the stream's own (SPH) density too, which inside
   * the ISM reflects the dense gas around it, otherwise the cap has no
   * effect there and streams could never escape the ISM. Comoving. */
  const float rho_max =
      cooling->firehose_rho_max * cosmo->a * cosmo->a * cosmo->a;
  const float rhocool = 0.5f * (min(p->chemistry_data.rho_ambient, rho_max) +
                                min(p->rho, rho_max));
  const float ucool = 0.5f * (p->chemistry_data.u_ambient + u);

  /* +ive if heating -ive if cooling*/
  p->cooling_data.mixing_layer_cool_time = cooling_time(
      phys_const, us, hydro_props, cosmo, cooling, p, xp, rhocool, ucool);
#ifdef FIREHOSE_DEBUG_CHECKS
  message("FIREHOSE_COOLING: id=%lld nH=%g T=%g rhoamb=%g Tamb=%g tcool=%g",
          p->id,
          hydro_get_physical_density(p, cosmo) * cooling->units.density_units *
              0.75 / 1.673e-24,
          u * cosmo->a_factor_internal_energy / cooling->temp_to_u_factor,
          p->chemistry_data.rho_ambient,
          p->chemistry_data.u_ambient * cosmo->a_factor_internal_energy /
              cooling->temp_to_u_factor,
          p->cooling_data.mixing_layer_cool_time);
#endif

  cooling_sputter_dust(us, cosmo, cooling, p, xp, dt);
}

/**
 * @brief Set subgrid ISM properties and initialise chemistry
 *
 * @param cooling The #cooling_function_data used in the run.
 * @param p Pointer to the particle data.
 */
void cooling_init_chemistry(
    const struct cooling_function_data *restrict cooling,
    struct part *restrict p) {

  /* If it's eligible for SF and metal-free, crudely self-enrich to
     very small level; needed to kick-start Grackle dust. Arbitrary
     amount to kick-start the model */
  const float init_dust_to_gas = 0.2f;
  const float total_Z = chemistry_get_total_metal_mass_fraction_for_cooling(p);
  const float self_Z =
      (1.f - init_dust_to_gas) * cooling->self_enrichment_metallicity;
  /* Small margin so round-off in the enriched total doesn't re-trigger this
   * every step (which would reset any dust grown since) */
  if (p->cooling_data.subgrid_temp > 0.f && total_Z < 0.99f * self_Z) {
    /* Solar metal mass fraction, total and in the dust elements, from
     * crackle's SolarAbundances (elements crackle does not know are not
     * self-enriched) */
    float Z_sun = 0.f, Z_sun_dust = 0.f;
    for (int i = chemistry_element_C; i < chemistry_element_count; i++) {
      const int slot = cooling_crackle_slot(i);
      if (slot < 0) continue;
      Z_sun += cooling->chemistry.SolarAbundances[slot];
      if (dust_element_of_chemistry(i) >= 0)
        Z_sun_dust += cooling->chemistry.SolarAbundances[slot];
    }

    /* Distribute the self-enrichment metallicity among elements
       assuming solar abundance ratios; init_dust_to_gas of the dust
       elements is in dust */
    p->chemistry_data.metal_mass_fraction_total = 0.f;
    for (int i = chemistry_element_C; i < chemistry_element_count; i++) {
      const int slot = cooling_crackle_slot(i);

      /* Elements crackle does not know keep what they have */
      if (slot < 0) {
        p->chemistry_data.metal_mass_fraction_total +=
            p->chemistry_data.metal_mass_fraction[i];
        continue;
      }
      const int d = dust_element_of_chemistry(i);

      /* fraction of gas mass in each element (gas phase) */
      const float f_gas = (d >= 0) ? 1.f - init_dust_to_gas : 1.f;
      p->chemistry_data.metal_mass_fraction[i] =
          f_gas * cooling->chemistry.SolarAbundances[slot] *
          cooling->self_enrichment_metallicity;
      p->chemistry_data.metal_mass_fraction[i] /= Z_sun;

      /* Update to the new metal mass fraction */
      p->chemistry_data.metal_mass_fraction_total +=
          p->chemistry_data.metal_mass_fraction[i];

      /* fraction of dust mass in each element (sums to 1) */
      if (d >= 0)
        p->cooling_data.dust_mass_fraction[d] =
            cooling->chemistry.SolarAbundances[slot] / Z_sun_dust;
    }

    /* Dust holds the remaining init_dust_to_gas of the dust elements */
    p->cooling_data.dust_mass = init_dust_to_gas *
                                cooling->self_enrichment_metallicity *
                                (Z_sun_dust / Z_sun) * p->mass;

    p->chemistry_data.metal_mass_fraction[chemistry_element_He] =
        cooling->chemistry.SolarAbundances[0];
    /* Since He is fixed at SolarAbundances[0], make sure the hydrogen
       fraction makes sense, i.e. X_H + Y_He + Z + D = 1. */
    p->chemistry_data.metal_mass_fraction[chemistry_element_H] =
        1.f - p->chemistry_data.metal_mass_fraction[chemistry_element_He] -
        p->chemistry_data.metal_mass_fraction_total -
        p->cooling_data.dust_mass / p->mass;

    if (p->chemistry_data.metal_mass_fraction[chemistry_element_H] > 1.f ||
        p->chemistry_data.metal_mass_fraction[chemistry_element_H] < 0.f) {
      for (int i = chemistry_element_H; i < chemistry_element_count; i++) {
        warning("\telem[%d] is %g", i,
                p->chemistry_data.metal_mass_fraction[i]);
      }

      error(
          "Hydrogen fraction exeeds unity or is negative for"
          " particle id=%lld",
          p->id);
    }
  }
}

/**
 * @brief Apply the cooling function to a particle.
 *
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param cosmo The current cosmological model.
 * @param hydro_props The #hydro_props.
 * @param floor_props Properties of the entropy floor.
 * @param cooling The #cooling_function_data used in the run.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the particle' extended data.
 * @param iact_rates Interaction rates for radiative transfer (if used)
 * @param dt The time-step of this particle.
 * @param dt_therm The time-step operator used for thermal quantities.
 * @param set_u_directly If 1 (RT), set the new internal energy directly and
 * zero du/dt; if 0, set du/dt to the Grackle energy change over dt_therm.
 */
void cooling_do_grackle_cooling(
    const struct phys_const *restrict phys_const,
    const struct unit_system *restrict us,
    const struct cosmology *restrict cosmo,
    const struct hydro_props *hydro_props,
    const struct entropy_floor_properties *floor_props,
    const struct cooling_function_data *restrict cooling,
    struct part *restrict p, struct xpart *restrict xp, gr_float *iact_rates,
    const double dt, const double dt_therm, const int set_u_directly) {

  /* Self-enrich gas if very low metallicity */
  cooling_init_chemistry(cooling, p);

  /* Collect information about galaxy that the particle belongs to */
  const float galaxy_mstar = p->galaxy_data.stellar_mass;
  const float galaxy_ssfr = p->galaxy_data.specific_sfr;

  /* Compute the entropy floor */
  // const double T_warm = entropy_floor_temperature(p, cosmo, floor_props);
  const double T_warm = warm_ISM_temperature(p, cooling, phys_const, cosmo);
  const double u_warm = cooling_convert_warm_ISM_temp_to_u(T_warm, cooling, p);

  /* Compute the ISRF */
  p->cooling_data.G0 =
      fmax(cooling_compute_G0(p, xp, p->cooling_data.subgrid_dens, T_warm,
                              cooling, galaxy_mstar, galaxy_ssfr),
           0.);

  /* Do grackle cooling */
  const float u_old = hydro_get_physical_internal_energy(p, xp, cosmo);
  // const float T_old = cooling_get_temperature( phys_const, hydro_props, us,
  // cosmo, cooling, p, xp); // for debugging only
  gr_float u_new = u_old;
  u_new = cooling_grackle_driver(phys_const, us, cosmo, hydro_props, cooling, p,
                                 xp, iact_rates, dt, T_warm, 0);

  /* Apply simulation-wide minimum temperature */
  u_new = max(u_new, hydro_props->minimal_internal_energy);

  /* Assign new thermal energy to particle */
  /* Calculate the cooling rate */
  float cool_du_dt = (u_new - u_old) / dt_therm;

  if (p->cooling_data.subgrid_temp == 0.) {
    /* Normal cooling; check that we are not going to go
     * below any of the limits */
    if (u_new > GRACKLE_HEATLIM * u_old) u_new = GRACKLE_HEATLIM * u_old;
    if (u_new < GRACKLE_COOLLIM * u_old) u_new = GRACKLE_COOLLIM * u_old;
    // u_new = max(u_new, u_warm);
    /* Get physical internal energy du_dt*/
    // const float hydro_du_dt = hydro_get_physical_internal_energy_dt(p,
    // cosmo);

    /* Rennehan: Recompute the actual thermal evolution after setting min/max */
    cool_du_dt = (u_new - u_old) / dt_therm;

    if (set_u_directly) {
      /* RT: directly set the new internal energy of the particle, and set
       * du/dt = 0 */
      hydro_set_physical_internal_energy(p, xp, cosmo, u_new);
      hydro_set_physical_internal_energy_dt(p, cosmo, 0.);
    } else {
      /* Update the internal energy time derivative,
       * which will be evolved later */
      hydro_set_physical_internal_energy_dt(p, cosmo, cool_du_dt);
    }

    /* If there is any dust outside of the ISM, sputter it
     * back into gas phase metals */
    cooling_sputter_dust(us, cosmo, cooling, p, xp, dt);
  } else {
    /* Particle is in subgrid mode; result is stored in subgrid_temp */
    p->cooling_data.subgrid_temp = cooling_convert_u_to_temp(
        u_new, xp->cooling_data.e_frac, cooling, p, xp);

    /* Set the subgrid cold ISM fraction for particle */
    /* Get H number density */
    const double rho = hydro_get_physical_density(p, cosmo);
    const float X_H =
        chemistry_get_metal_mass_fraction_for_cooling(p)[chemistry_element_H];
    const double n_H = rho * X_H / phys_const->const_proton_mass;

    const double fcold_max = cooling_compute_cold_ISM_fraction(n_H, cooling);

    /* Compute cooling time in warm ISM component */
    float rhocool = p->rho * (1.f - p->cooling_data.subgrid_fcold);
    rhocool = fmax(rhocool, 0.001f * p->rho);
    const float u = hydro_get_comoving_internal_energy(p, xp);
    float tcool = cooling_time(phys_const, us, hydro_props, cosmo, cooling, p,
                               xp, rhocool, u);

    /* Evolve fcold upwards by cooling from warm ISM on
     * relevant cooling timescale */
    if (tcool < 0.f) {
      p->cooling_data.subgrid_fcold +=
          (fcold_max - p->cooling_data.subgrid_fcold) * (1.f - exp(dt / tcool));
    }

    /* Compare the adiabatic heating to the heating required to heat the
     * gas back up to the equation of state line, and assume this is the
     * fraction of cold cloud mass destroyed. */
    if (cooling->ism_adiabatic_heating_method == 2) {
      /* Use adiabatic du/dt to evaporate cold gas clouds, into warm phase */
      const double f_evap =
          hydro_get_physical_internal_energy_dt(p, cosmo) * dt_therm /
          (hydro_get_physical_internal_energy(p, xp, cosmo) - u_new);

      /* If it's in the ISM of a galaxy, suppress cold fraction */
      if (f_evap > 0.f && galaxy_mstar > 0.f) {
        p->cooling_data.subgrid_fcold *= max(1. - f_evap, 0.f);
      }
    }

    /* Set internal energy time derivative to 0 for overall particle */
    hydro_set_physical_internal_energy_dt(p, cosmo, 0.f);

    /* No cooling in warm phase since it is fixed on the EoS */
    cool_du_dt = 0.f;

    /* Force the overall particle to lie on the equation of state
    hydro_set_physical_internal_energy(p, xp, cosmo, u_warm);*/

    /* set subgrid properties for use in SF routine */
    cooling_set_particle_subgrid_properties(phys_const, us, cosmo, hydro_props,
                                            floor_props, cooling, p, xp);

    /* Overall particle u is combination of EOS u and subgrid u */
    const float u_part = p->cooling_data.subgrid_fcold * u_new +
                         (1. - p->cooling_data.subgrid_fcold) * u_warm;
    hydro_set_physical_internal_energy(p, xp, cosmo, u_part);
  } /* subgrid mode */

  /* Store the radiated energy */
  xp->cooling_data.radiated_energy -= hydro_get_mass(p) * cool_du_dt * dt_therm;
}

/**
 * @brief Decay the SN rate deposited on a gas particle by stellar feedback,
 * once it has been used over a step of length dt. With no smoothing time,
 * it is dropped, so that each deposit is used for one step.
 *
 * @param cooling The #cooling_function_data used in the run.
 * @param p Pointer to the particle data.
 * @param dt The time-step of this particle.
 */
static INLINE void cooling_decay_SNe_rate(
    const struct cooling_function_data *restrict cooling,
    struct part *restrict p, const double dt) {
#if COOLING_GRACKLE_MODE >= 2
  if (cooling->SNe_smoothing_time > 0.) {
    p->feedback_data.SNe_ThisTimeStep *= exp(-dt / cooling->SNe_smoothing_time);
  } else {
    p->feedback_data.SNe_ThisTimeStep = 0.f;
  }
#endif
}

/**
 * @brief Apply the cooling function to a particle.
 *
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param cosmo The current cosmological model.
 * @param hydro_props The #hydro_props.
 * @param floor_props Properties of the entropy floor.
 * @param pressure_floor_props Properties of the pressure floor.
 * @param cooling The #cooling_function_data used in the run.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the particle' extended data.
 * @param dt The time-step of this particle.
 * @param dt_therm The time-step operator used for thermal quantities.
 * @param time The current time (since the Big Bang or start of the run) in
 * internal units.
 */
void cooling_cool_part(const struct phys_const *restrict phys_const,
                       const struct unit_system *restrict us,
                       const struct cosmology *restrict cosmo,
                       const struct hydro_props *hydro_props,
                       const struct entropy_floor_properties *floor_props,
                       const struct pressure_floor_props *pressure_floor_props,
                       const struct cooling_function_data *restrict cooling,
                       struct part *restrict p, struct xpart *restrict xp,
                       const double dt, const double dt_therm,
                       const double time) {

  /* Compute cooling time and other quantities needed for firehose */
  firehose_cooling_and_dust(phys_const, us, cosmo, hydro_props, cooling, p, xp,
                            dt);

  /* Update the subgrid properties */
  cooling_set_particle_subgrid_properties(phys_const, us, cosmo, hydro_props,
                                          floor_props, cooling, p, xp);

  /* No cooling happens over zero time */
  if (dt == 0.f || dt_therm == 0.f) return;

  /* In KIARA_RT we do cooling in rt_thermochemistry, so we don't do it here.
   * This task runs before stellar feedback and RT in each step, so the SN
   * rate deposited last step has already been used by RT: decay it here. */
  if (cooling->do_cooling_in_rt) {
    cooling_decay_SNe_rate(cooling, p, dt);
    return;
  }

  /* Cool unless the particle is decoupled or its cooling is shut off */
  if (!p->decoupled && !(p->feedback_data.cooling_shutoff_delay_time > 0.f)) {

    /* Interaction rates for RT; not used here */
    gr_float iact_rates[5] = {0., 0., 0., 0., 0.};

    /* Output u for debug. */
    p->cooling_data.u_cooling_before =
        hydro_get_physical_internal_energy(p, xp, cosmo);

    /* Do the cooling and chemistry */
    cooling_do_grackle_cooling(phys_const, us, cosmo, hydro_props, floor_props,
                               cooling, p, xp, iact_rates, dt, dt_therm,
                               /*set_u_directly=*/0);

    /* Output du_dt and u for debug. */
    p->cooling_data.du_dt_cooling =
        hydro_get_physical_internal_energy_dt(p, cosmo);
    p->cooling_data.u_cooling =
        hydro_get_physical_internal_energy(p, xp, cosmo);

    /* Update cooling times this step. */
    p->cooling_data.ncool_this_step += 1;

    /* Record this cooling event */
    xp->cooling_data.time_last_event = time;
  }

  /* The SN rate deposited by stellar feedback has now been used (for the
   * dust and G0) over this step */
  cooling_decay_SNe_rate(cooling, p, dt);
}

/**
 * @brief Set the subgrid properties (rho, T) of the gas particle for use
 *        in SF routine
 *
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param cosmo The current cosmological model.
 * @param hydro_props the hydro_props struct
 * @param floor_props Properties of the entropy floor.
 * @param cooling The #cooling_function_data used in the run.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 */
void cooling_set_particle_subgrid_properties(
    const struct phys_const *phys_const, const struct unit_system *us,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct entropy_floor_properties *floor_props,
    const struct cooling_function_data *cooling, struct part *p,
    struct xpart *xp) {

  /* No subgrid ISM if particle is decoupled */
  if (p->decoupled || p->feedback_data.cooling_shutoff_delay_time > 0.f) {
    /* Make sure these are always set for the wind particles */
    cooling_subgrid_exit(&p->cooling_data, &xp->cooling_data,
                         hydro_get_physical_density(p, cosmo));

    return;
  }

  /* Get temperature of overall particle */
  const double u = hydro_get_physical_internal_energy(p, xp, cosmo);
  const float temperature =
      cooling_convert_u_to_particle_temp(u, cooling, p, xp);

  /* Get density */
  const double rho = hydro_get_physical_density(p, cosmo);

  /* Subgrid model is on if particle is in the Jeans EOS regime */
  const double T_warm = warm_ISM_temperature(p, cooling, phys_const, cosmo);
  // entropy_floor_gas_temperature( rho, rho_com, cosmo, floor_props);
  const double u_warm = cooling_convert_warm_ISM_temp_to_u(T_warm, cooling, p);

  /* Check if it is in subgrid mode: Must be in Jeans EoS regime
   * and have nonzero cold gas */
  if (T_warm > 0 && u < u_warm * cooling->entropy_floor_margin) {
    /* YES: If first time in subgrid, set temperature to particle T,
       otherwise limit to particle T */
    if (p->cooling_data.subgrid_temp == 0.) {
      p->cooling_data.subgrid_temp = temperature;
#if COOLING_GRACKLE_MODE >= 2
      /* H2 the particle still carries (crackle evolves it outside the ISM) */
      const float H2_frac_part =
          xp->cooling_data.H2I_frac + xp->cooling_data.H2II_frac;
#endif
      /* Reset grackle subgrid quantities assuming neutral gas */
      cooling_grackle_init_part(cooling, p, xp);
      /* Initialize ISM cold fraction */
      p->cooling_data.subgrid_fcold = cooling->cold_ISM_frac;
#if COOLING_GRACKLE_MODE >= 2
      /* Assign the H2 to the cold phase (species are cold-phase values in
       * the subgrid ISM), limited to all of its hydrogen being molecular.
       * Like the other species it is reset to neutral, so H2+ becomes H2. */
      if (H2_frac_part > 0.f && p->cooling_data.subgrid_fcold > 0.f) {
        const float H2_frac_cold =
            min(H2_frac_part / p->cooling_data.subgrid_fcold,
                xp->cooling_data.HI_frac);
        xp->cooling_data.H2I_frac = H2_frac_cold;
        xp->cooling_data.HI_frac -= H2_frac_cold;
      }
#endif
    } else {
      /* Subgrid temperature should be no higher
       * than overall particle temperature */
      p->cooling_data.subgrid_temp =
          min(p->cooling_data.subgrid_temp, temperature);
    }

    /* Compute the cold-cloud density assuming pressure equilibrium. The
     * cold-phase energy uses the cold-phase species (reset to neutral above
     * on entry), so the pressure balance accounts for both phases' mu. */
    const double u_subgrid = cooling_convert_temp_to_u(
        p->cooling_data.subgrid_temp, xp->cooling_data.e_frac, cooling, p, xp);
    p->cooling_data.subgrid_dens =
        cooling_compute_subgrid_density(rho, u, u_subgrid, cooling);
  } else {
    /* NO: not in subgrid mode (subgrid temperature 0, density the actual
     * particle's physical density, no cold gas or ISRF) */
    cooling_subgrid_exit(&p->cooling_data, &xp->cooling_data, rho);
  }
}

/**
 * @brief Returns the subgrid temperature of a particle.
 *
 * @param p The particle.
 * @param xp The extended particle data.
 * @return The subgrid temperature in internal units.
 */
float cooling_get_subgrid_temperature(const struct part *p,
                                      const struct xpart *xp) {
  return p->cooling_data.subgrid_temp;
}

/**
 * @brief Returns the subgrid density of a particle.
 *
 * @param p The particle.
 * @param xp The extended particle data.
 * @return The subgrid density in physical internal units.
 */
float cooling_get_subgrid_density(const struct part *p,
                                  const struct xpart *xp) {
  return p->cooling_data.subgrid_dens;
}

/**
 * @brief Compute the y-Compton contribution of a #part based on the cooling
 * function.
 *
 * Does not exist in this model. We return 0.
 *
 * @param phys_const #phys_const data structure.
 * @param hydro_props The properties of the hydro scheme.
 * @param us The internal system of units.
 * @param cosmo #cosmology data structure.
 * @param cooling #cooling_function_data struct.
 * @param p #part data.
 * @param xp Pointer to the #xpart data.
 * @return The y-Compton contribution (always 0 here).
 */
double cooling_get_ycompton(const struct phys_const *phys_const,
                            const struct hydro_props *hydro_props,
                            const struct unit_system *us,
                            const struct cosmology *cosmo,
                            const struct cooling_function_data *cooling,
                            const struct part *p, const struct xpart *xp) {

  return 0.;
}

/**
 * @brief Computes the cooling time-step.
 *
 * We return FLT_MAX so as to impose no limit on the time-step,
 * since Grackle sub-cycles the cooling as needed.
 *
 * @param cooling The #cooling_function_data used in the run.
 * @param phys_const The physical constants in internal units.
 * @param cosmo The #cosmology.
 * @param us The internal system of units.
 * @param hydro_props The #hydro_props.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the particle extra data
 * @return The cooling time-step limit (FLT_MAX: no limit).
 */
float cooling_timestep(const struct cooling_function_data *restrict cooling,
                       const struct phys_const *restrict phys_const,
                       const struct cosmology *restrict cosmo,
                       const struct unit_system *restrict us,
                       const struct hydro_props *hydro_props,
                       const struct part *restrict p,
                       const struct xpart *restrict xp) {

  return FLT_MAX;
}

/**
 * @brief Split the coolong content of a particle into n pieces
 *
 * @param p The #part.
 * @param xp The #xpart.
 * @param n The number of pieces to split into.
 */
void cooling_split_part(struct part *p, struct xpart *xp, double n) {

  xp->cooling_data.radiated_energy /= n;
#if COOLING_GRACKLE_MODE >= 2
  p->cooling_data.dust_mass /= n;
#endif
}

/**
 * @brief Initialises the cooling unit system.
 *
 * @param us The current internal system of units.
 * @param phys_const The #phys_const.
 * @param cooling The cooling properties to initialize
 */
void cooling_init_units(const struct unit_system *us,
                        const struct phys_const *phys_const,
                        struct cooling_function_data *cooling) {

  /* These are conversions from code units to cgs. */

  /* first cosmo:
     units for the expansion factor,
     use cosmological redshift!
     a_value is updated in cooling_update() */
  cooling->units.a_units = 1.0;
  if (cooling->redshift == -1) {
    cooling->units.a_value = 0.01;
  } else {
    cooling->units.a_value = 1.0;
  }

  /* We assume here all quantities to be in physical coordinate */
  cooling->units.comoving_coordinates = 0;

  /* CMB temperature */
  cooling->T_CMB_0 = phys_const->const_T_CMB_0 *
                     units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE);

  /* then units */
  /* converts physical density to cgs number density for H */
  cooling->units.density_units =
      units_cgs_conversion_factor(us, UNIT_CONV_DENSITY);
  cooling->units.length_units =
      units_cgs_conversion_factor(us, UNIT_CONV_LENGTH);
  cooling->units.time_units = units_cgs_conversion_factor(us, UNIT_CONV_TIME);
  cooling->units.velocity_units =
      units_cgs_conversion_factor(us, UNIT_CONV_VELOCITY);

  cooling->temp_to_u_factor =
      phys_const->const_boltzmann_k /
      (hydro_gamma_minus_one * phys_const->const_proton_mass *
       units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE));
  cooling->dudt_units =
      units_cgs_conversion_factor(us, UNIT_CONV_ENERGY_PER_UNIT_MASS) /
      units_cgs_conversion_factor(us, UNIT_CONV_TIME);
  cooling->energy_units = units_cgs_conversion_factor(us, UNIT_CONV_ENERGY);

  /* converts galaxy sSFR into G0 by scaling to MW values */
  const double time_to_yr = units_cgs_conversion_factor(us, UNIT_CONV_TIME) /
                            (365.25f * 24.f * 60.f * 60.f);
  const double mass_to_solar_mass = 1.f / phys_const->const_solar_mass;
  const double length_to_pc =
      units_cgs_conversion_factor(us, UNIT_CONV_LENGTH) / 3.08567758e18f;

  cooling->time_to_Myr = time_to_yr * 1.e-6;
  cooling->mass_to_solar_mass = mass_to_solar_mass;

  const double vel_to_km_s =
      units_cgs_conversion_factor(us, UNIT_CONV_VELOCITY) * 1.e-5;
  cooling->potential_to_kms2 = vel_to_km_s * vel_to_km_s;
  cooling->ff_const = sqrt(3. * M_PI / (32. * phys_const->const_newton_G));

  /* Get the surface density unit Msun / pc^2 in internal units */
  const double Msun_per_pc2 =
      phys_const->const_solar_mass /
      (phys_const->const_parsec * phys_const->const_parsec);

  cooling->surface_rho_to_Msun_per_parsec2 = 1. / Msun_per_pc2;

  /* G0 for MW=1.6 (Parravano etal 2003).  */
  /* Scaled to SFR density in solar neighborhood =0.002 Mo/Gyr/pc^3
     (J. Isern 2019) */
  // const float sfr_density_solar = 0.002f;
  /* SFR within solar circle is 1.5 Mo/yr */
  const float sfr_density_solar = 0.066f;
  cooling->G0_factor1 = 1.6f * mass_to_solar_mass /
                        (sfr_density_solar * time_to_yr * 1.e-9) /
                        (length_to_pc * length_to_pc * length_to_pc);

  /* Calibrated to sSFR for MW=2.71e-11 (Licquia etal 2015) */
  cooling->G0_factor2 = 1.6f / (2.71e-11f * time_to_yr);

  /* Calibrated to an SN rate of 0.015 per year for MW (BC Reed 2005);
   * multiplies the SN rate in internal units */
  cooling->G0_factorSNe = 1.6f / (0.015f * time_to_yr);
}

/**
 * @brief Initialises Grackle.
 *
 * @param cooling The cooling properties to initialize
 */
void cooling_init_grackle(struct cooling_function_data *cooling) {

#ifdef SWIFT_DEBUG_CHECKS
  /* enable verbose for grackle */
  grackle_verbose = 1;
#endif

  chemistry_data *chemistry = &cooling->chemistry;

  /* Create a chemistry object for parameters and rate data. */
  if (set_default_chemistry_parameters(chemistry) == 0) {
    error("Error in set_default_chemistry_parameters.");
  }

  // Set parameter values for chemistry & cooling

  // Flag to activate the grackle machinery:
  chemistry->use_grackle = 2;  // 1=original grackle, 2=crackle
  /* Flag to include radiative cooling and actually update the thermal energy
   * during the chemistry solver. If off, the chemistry species will still be
   * updated. The most common reason to set this to off is to iterate the
   * chemistry network to an equilibrium state. */
  chemistry->with_radiative_cooling = 1;

  /* Flag to control which primordial chemistry network is used (set by Config
   * file) */
  chemistry->primordial_chemistry = COOLING_GRACKLE_MODE;

  /** Flag to enable H2 formation on dust grains, dust cooling, and dust-gas
   * heat transfer follow Omukai (2000). This assumes that the dust to gas ratio
   * scales with the metallicity. Default: 0. */
  chemistry->h2_on_dust = 0;

  /* Flag to enable metal cooling using the Cloudy tables. If enabled, the
   * cooling table to be used must be specified with the grackle_data_file
   * parameter. Default: 0. */
  chemistry->metal_cooling = cooling->with_metal_cooling;

  /* Flag to enable an effective CMB temperature floor. This is implemented by
   * subtracting the value of the cooling rate at TCMB from the total cooling
   * rate. Default: 1. */
  chemistry->cmb_temperature_floor = 1;

  /* Flag to enable a UV background. If enabled, the cooling table to be used
   * must be specified with the grackle_data_file parameter. Default: 0. */
  chemistry->UVbackground = cooling->with_uv_background;

  /* Path to the data file containing the metal cooling and UV background
   * tables: */
  chemistry->grackle_data_file = cooling->cloudy_table;

  /* The ratio of specific heats for an ideal gas. A direct calculation for the
   * molecular component is used if primordial_chemistry > 1. Default: 5/3. */
  chemistry->Gamma = hydro_gamma;

  /* Flag to control which three-body H2 formation rate is used. */
  chemistry->three_body_rate = 0;

  /* Flag to enable H2 collision-induced emission cooling from Ripamonti & Abel
   * (2004). Default: 0. */
  chemistry->cie_cooling = 0;

  /* Flag to enable H2 cooling attenuation from Ripamonti & Abel (2004).
   * Default: 0 */
  chemistry->h2_optical_depth_approximation = 0;
  /* Flag to enable a spatially uniform heating term approximating
   * photo-electric heating from dust from Tasker & Bryan (2008). Default: 0. */
  chemistry->photoelectric_heating = 0;
  /* photo-electric on [but not adjusted to local background, beware!] */
  chemistry->photoelectric_heating_rate = 8.5e-26;

  /* Flag to enable Compton heating from an X-ray background following Madau &
   * Efstathiou (1999). Default: 0. */
  chemistry->Compton_xray_heating = 0;

  /* Intensity of a constant Lyman-Werner H2 photo-dissociating radiation field,
   * in units of 10-21 erg s-1 cm-2 Hz-1 sr-1. Default: 0. */
  chemistry->LWbackground_intensity = 0;

  /* Flag to enable suppression of Lyman-Werner flux due to Lyman-series
   * absorption
   *    (giving a sawtooth pattern), taken from Haiman & Abel, & Rees (2000).
   *    Default: 0. */
  chemistry->LWbackground_sawtooth_suppression = 0;

  /* volumetric heating rates is being provided in the volumetric_heating_rate
   * field of grackle_field_data */
  chemistry->use_volumetric_heating_rate =
      cooling->provide_volumetric_heating_rates;

  /* specific heating rates is being provided in the specific_heating_rate field
   * of grackle_field_data */
  chemistry->use_specific_heating_rate =
      cooling->provide_specific_heating_rates;

  /* Set parameters of temperature floor: 0=none, 1=provide scalar,
   * 2=provide array */
  chemistry->use_temperature_floor = 2;

  /* arrays of ionization and heating rates from radiative transfer solutions
   * are being provided */
  chemistry->use_radiative_transfer = cooling->with_radiative_transfer;

  /* must be enabled to couple the passed radiative transfer fields to the
   * chemistry solver */
  chemistry->radiative_transfer_coupled_rate_solver = 0;

  /* enable intermediate stepping in applying radiative transfer fields to
   * chemistry solver. */
  chemistry->radiative_transfer_intermediate_step = 0;

  /* only use hydrogen ionization and heating rates from the radiative transfer
   * solutions. */
  chemistry->radiative_transfer_hydrogen_only = 0;

  /* Use self-shielding; 0=none, 1=HI only, 2=HI+HeI, 3=HI+HeI but
   * set HeII rates to 0 */
  chemistry->self_shielding_method = cooling->self_shielding_method;

  /* control behaviour of Grackle sub-step integrator */
  chemistry->max_iterations = cooling->max_step;
  chemistry->exit_after_iterations_exceeded = 0;
  chemistry->accuracy = cooling->timestep_accuracy;

  /* control behaviour of Grackle sub-step integration damping */
  if (cooling->grackle_damping_interval > 0) {
    chemistry->use_subcycle_timestep_damping = 1;
    chemistry->subcycle_timestep_damping_interval =
        cooling->grackle_damping_interval;
  } else {
    chemistry->use_subcycle_timestep_damping = 0;
    chemistry->subcycle_timestep_damping_interval = 0;
  }
  /* run on a single thread since Swift sends each particle to a single thread
   *chemistry->omp_nthreads = 1; */

  /* Turn on Li+ 2019 dust evolution model */
  chemistry->use_dust_evol = cooling->use_grackle_dust_evol;
  chemistry->use_dust_density_field = cooling->use_grackle_dust_evol;

  /* Load dust evolution parameters */
  if (cooling->use_grackle_dust_evol == 1) {
    chemistry->dust_destruction_eff = cooling->dust_destruction_eff;
    chemistry->sne_coeff = cooling->dust_sne_coeff;
    chemistry->sne_shockspeed = cooling->dust_sne_shockspeed;
    chemistry->dust_grainsize = cooling->dust_grainsize;
    chemistry->dust_growth_densref = cooling->dust_growth_densref;
    chemistry->dust_growth_tauref = cooling->dust_growth_tauref;
    /* Enable dust temperature calculation using ISRF */
    chemistry->metal_cooling = 1;
    chemistry->dust_chemistry = 1;
    chemistry->h2_on_dust = 1;
    chemistry->use_isrf_field = 1;
    chemistry->H2_self_shielding = 3;
    /* 2 means we specify the H2 shielding
       length ourselves (the gas smoothing length) */
    chemistry->H2_custom_shielding = 2;
    /* Solar abundances to pass to Grackle:
       He  (10.93 in units where log[H]=12, so photospheric mass fraction
          -> Y=0.2485 [Hydrogen X=0.7381]; Anders+Grevesse Y=0.2485, X=0.7314)
       C   (8.43 -> 2.38e-3, AG=3.18e-3)
       N   (7.83 -> 0.70e-3, AG=1.15e-3)
       O   (8.69 -> 5.79e-3, AG=9.97e-3)
       Ne  (7.93 -> 1.26e-3, AG=1.72e-3)
       Mg  (7.60 -> 7.14e-4, AG=6.75e-4)
       Si  (7.51 -> 6.71e-4, AG=7.30e-4)
       S   (7.12 -> 3.12e-4, AG=3.80e-4)
       Ca  (6.34 -> 0.65e-4, AG=0.67e-4)
       Fe (7.50 -> 1.31e-3, AG=1.92e-3)
       */
    chemistry->SolarAbundances[0] = 0.2485;
    chemistry->SolarAbundances[1] = 2.38e-3;
    chemistry->SolarAbundances[2] = 0.70e-3;
    chemistry->SolarAbundances[3] = 5.79e-3;
    chemistry->SolarAbundances[4] = 1.26e-3;
    chemistry->SolarAbundances[5] = 7.14e-4;
    chemistry->SolarAbundances[6] = 6.71e-4;
    chemistry->SolarAbundances[7] = 3.12e-4;
    chemistry->SolarAbundances[8] = 0.65e-4;
    chemistry->SolarAbundances[9] = 1.31e-3;
  } else {
    chemistry->use_dust_evol = 0;
  }

  cooling->use_grackle_h2_form =
      cooling->use_grackle_dust_evol && COOLING_GRACKLE_MODE >= 2;

#ifdef CRACKLE_VALIDATE
  /* Initialise the external library with the same parameters */
  crackle_validate_init(crackle_validate_global, chemistry, &cooling->units,
                        sizeof(chemistry_data), sizeof(code_units),
                        sizeof(grackle_field_data));
#endif

  /* Initialize the chemistry object. */
  if (initialize_chemistry_data(&cooling->units) == 0) {
    error("Error in initialize_chemistry_data.");
  }

  /* What crackle needs at the current redshift (updated every step) */
  cooling_prepare_crackle_step(cooling);
}

/**
 * @brief Initialises the cooling properties.
 *
 * @param parameter_file The parsed parameter file.
 * @param us The current internal system of units.
 * @param phys_const The physical constants in internal units.
 * @param hydro_props The properties of the hydro scheme.
 * @param cooling The cooling properties to initialize
 */
void cooling_init_backend(struct swift_params *parameter_file,
                          const struct unit_system *us,
                          const struct phys_const *phys_const,
                          const struct hydro_props *hydro_props,
                          struct cooling_function_data *cooling) {

  if (GRACKLE_NPART != 1)
    error("Grackle with multiple particles not implemented");

  /* read parameters */
  cooling_read_parameters(parameter_file, cooling, phys_const, us);

  /* Maximal density of the firehose mixing layer: the same cap as the
   * ambient density in the chemistry
   * (KIARAChemistry:firehose_nh_ambient_max_cgs, a physical H number density),
   * converted to a physical mass density */
  cooling->firehose_rho_max =
      parser_get_opt_param_double(
          parameter_file, "KIARAChemistry:firehose_nh_ambient_max_cgs", 0.1) *
      phys_const->const_proton_mass /
      (hydro_props->hydrogen_mass_fraction *
       units_cgs_conversion_factor(us, UNIT_CONV_NUMBER_DENSITY));

  /* Set up the units system. */
  cooling_init_units(us, phys_const, cooling);

  /* Set up grackle */
  cooling_init_grackle(cooling);
}

/**
 * @brief Clean-up the memory allocated for the cooling routines
 *
 * @param cooling the cooling data structure.
 */
void cooling_clean(struct cooling_function_data *cooling) {
  //_free_chemistry_data(&cooling->chemistry, &grackle_rates);
}

/**
 * @brief Write a cooling struct to the given FILE as a stream of bytes.
 *
 * Nothing to do beyond writing the structure from the stream.
 *
 * @param cooling the struct
 * @param stream the file stream
 */
void cooling_struct_dump(const struct cooling_function_data *cooling,
                         FILE *stream) {
  restart_write_blocks((void *)cooling, sizeof(struct cooling_function_data), 1,
                       stream, "cooling", "cooling function");
}

/**
 * @brief Restore a hydro_props struct from the given FILE as a stream of
 * bytes.
 *
 * Nothing to do beyond reading the structure from the stream.
 *
 * @param cooling the struct
 * @param stream the file stream
 * @param cosmo #cosmology structure
 */
void cooling_struct_restore(struct cooling_function_data *cooling, FILE *stream,
                            const struct cosmology *cosmo) {
  restart_read_blocks((void *)cooling, sizeof(struct cooling_function_data), 1,
                      stream, NULL, "cooling function");

  /* Set up grackle */
  cooling_init_grackle(cooling);
}
