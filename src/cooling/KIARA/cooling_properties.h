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
#ifndef SWIFT_COOLING_PROPERTIES_KIARA_H
#define SWIFT_COOLING_PROPERTIES_KIARA_H

/* include grackle */
#include "crackle/grackle.h"

/**
 * @file src/cooling/KIARA/cooling_properties.h
 * @brief Properties of the KIARA (grackle/crackle) cooling function
 */

/**
 * @brief Properties of the cooling function.
 */
struct cooling_function_data {

  /*! Filename of the Cloudy Table */
  char cloudy_table[200];

  /*! Enable/Disable UV backgroud */
  int with_uv_background;

  /*! Enable/Disable UV backgroud */
  int with_radiative_transfer;

  /*! Redshift to use for the UV backgroud (-1 to use cosmological one) */
  double redshift;

  /*! unit system */
  code_units units;

  /*! Temperature of the CMB at present day (for quick access) */
  double T_CMB_0;

  /*! k_Boltz/m_p plus conversion factor for converting u<->T */
  double temp_to_u_factor;

  /*! Convert time to Myr */
  double time_to_Myr;

  /*! Convert mass to Msun */
  double mass_to_solar_mass;

  /*! Convert potential to (km/s)^2 */
  double potential_to_kms2;

  /*! Converts to Msun/pc^2 */
  double surface_rho_to_Msun_per_parsec2;

  /*! Constant needed for computing free-fall time */
  double ff_const;

  /*! conversion unit factor for rate of change of thermal energy */
  double dudt_units;

  /*! conversion unit factor for energy into ergs */
  double energy_units;

  /*! grackle chemistry data */
  chemistry_data chemistry;

  /*! Enable/Disable metal cooling */
  int with_metal_cooling;

  /*! User provide volumetric heating rates */
  int provide_volumetric_heating_rates;

  /*! User provide specific heating rates */
  int provide_specific_heating_rates;

  /*! Self shielding method (1 -> 3 for grackle's ones, 0 for none and -1 for
   * GEAR) */
  int self_shielding_method;

  /*! What to do with adiabatic du/dt when in ISM mode: 0=no adiabatic heating,
   * 1=evolve in grackle, 2=use to evaporate cold ism */
  int ism_adiabatic_heating_method;

  /*! number of step max for first init */
  int max_step;

  /*! max fractional change in quantities in single grackle substep */
  double timestep_accuracy;

  /*! parameter to control how fast grackle damps oscillatory behaviour
   * (lower=more aggressive) */
  int grackle_damping_interval;

  /*! Duration for switching off cooling after an event (e.g. supernovae) */
  double thermal_time;

  /*! track dust growth and destruction (only available in KIARA) */
  int use_grackle_dust_evol;

  /*! Maximal (physical) density of the gas in the firehose mixing layer, for
   * its cooling time (internal units) */
  double firehose_rho_max;

  /*! Timescale (internal units) over which the SN rate deposited on gas
   * particles by stellar feedback decays (0: used for one cooling step) */
  double SNe_smoothing_time;

  /*! track H2 formation; this is set within the code based on selection options
   */
  int use_grackle_h2_form;

  /* G0 conversion factors, scaling to the MW value (G0 = 1.6) based on
   * local/global galaxy properties */

  /*! G0 per unit local SFR density (internal units), normalised to the
   * solar-circle SFR density */
  double G0_factor1;

  /*! G0 per unit galaxy sSFR (internal units), normalised to the MW sSFR */
  double G0_factor2;

  /*! G0 per unit SN rate (internal units), normalised to the MW SN rate */
  double G0_factorSNe;

  /*! G0 per unit galaxy SFR / distance^2 (internal units), normalised to the
   * MW at the solar circle (method 7) */
  double G0_factor_galaxy;

  /*! Method 7: minimum (physical) softening of the galaxy-scale field */
  double G0_softening_min;

  /*! Method 7: fraction of G0 from the local SFR density; the rest is from
   * the galaxy-scale model */
  double G0_local_fraction;

  /* Dust parameters passed to grackle; see sample yml file */

  /*! Dust destruction efficiency by SN shocks (destroyed/shocked mass) */
  double dust_destruction_eff;

  /*! SN energy in units of 1e51 erg, for dust destruction */
  double dust_sne_coeff;

  /*! SN shock speed for dust destruction (km/s) */
  double dust_sne_shockspeed;

  /*! Dust grain size (micron) */
  double dust_grainsize;

  /*! Limit the oxygen in dust to the oxides of the refractory dust elements
   * (MgO, SiO2, CaO, FeO), as grains hold O only bound to them */
  int dust_oxygen_stoichiometric_cap;

  /*! H2 photodissociation by the ISRF G0, self- and dust-shielded (0 = off:
   * only the UVB dissociates H2) */
  int H2_isrf_dissociation;

  /*! Unshielded ISRF H2 photodissociation rate per Habing unit of G0 (s^-1) */
  double H2_isrf_dissociation_rate;

  /*! LW dust absorption cross-section per H nucleus at the local dust-to-gas
   * ratio (cm^2), for the dust shielding of the ISRF */
  double H2_dust_shielding_cross_section;

  /*! H2 shielding column from the particle's mean density, rather than the
   * subgrid cold-cloud density, times the shielding length */
  int H2_shield_mean_density;

  /*! Limit the subgrid ISM molecular fraction to the part of the particle's
   * column shielded by the Sternberg+14 atomic layer */
  int H2_sternberg_limit;

  /*! Density at which the dust growth time-scale is 1 Gyr (g/cm^3) */
  double dust_growth_densref;

  /*! Always 1 Gyr (the growth time-scale at dust_growth_densref); kept so the
   * struct, and so the restart files, are unchanged */
  double dust_growth_tauref;

  /*! For dust model, need self-enrichment up to a small metallicity to
   * kick-start dust */
  double self_enrichment_metallicity;

  /*! For subgrid model (eg KIARA) need a subgrid ISM fraction */
  double cold_ISM_frac;

  /*! For Grackle subgrid model, choose way to determine G0: 0=none; 1=Local
   * SFR density; 2=Global sSFR; 3=2 if sSFR != 0, else 1; -3: vice versa;
   * 4=SN rate; 5=sSFR + SN rate; 6=FIRE-based fit */
  int G0_computation_method;

  /*! For Grackle subgrid model, arbitrary multiplier for G0 */
  double G0_multiplier;

  /*! For Grackle subgrid model, set max density to avoid pointlessly
   * over-iterating in Grackle */
  double max_subgrid_density;

  /*! For Grackle subgrid model, minimum (physical) cold-cloud density, which
   * overrides pressure equilibrium when the ISM pressure is low */
  double min_subgrid_density;

  /*! Minimum self-shielding column length as a fraction of h:
   * MIN_SHIELD_SPACING_FRAC mean interparticle spacings */
  double min_shield_h_frac;

  /*! For Grackle subgrid model, inverse of threshold nH above which multi-phase
   * ISM model kicks in */
  double subgrid_threshold_n_H_inv;

  /*! For Grackle subgrid model, temperature at threshold nH */
  double subgrid_threshold_T;

  /*! For Grackle subgrid model, Power-law eqn of state for warm ISM component
   * above threshold n_H */
  double subgrid_warm_ism_EOS;

  /*! For Grackle subgrid model, factor above the warm-ISM EoS energy below
   * which gas is in subgrid mode */
  double entropy_floor_margin;

  /*! Option to use Cloudy lookup tables when outside ISM */
  int use_tables_outside_ism;

  /*! When using radiative transfer, set this on if you want thermochemistry
   * done within RT modules (eg in KIARART) */
  int do_cooling_in_rt;
};

#endif /* SWIFT_COOLING_PROPERTIES_KIARA_H */
