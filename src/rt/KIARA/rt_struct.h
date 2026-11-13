/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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
#ifndef SWIFT_RT_STRUCT_KIARA_H
#define SWIFT_RT_STRUCT_KIARA_H

/**
 * @file src/rt/KIARA/rt_struct.h
 * @brief Main header file for the KIARA M1 Closure radiative transfer struct.
 */

/* Additional RT data in hydro particle struct */
struct rt_part_data {

#ifdef SWIFT_KIARA_RT_DEBUG
  /*! DEBUG: cumulative photon energy changes (all groups, internal energy
   * units): 0 injected, 1 absorbed, 2 transport exchange, 3 negative-energy
   * clamps after transport, 4 redshifting, 5 drift correction */
  float debug_dE[6];

  /*! DEBUG: the part of debug_dE[2] from exchanges in which only this
   * particle was updated, its active partner being updated in its own call
   * (non-symmetric interactions); sums to zero over all the particles if
   * every such exchange has its counterpart */
  float debug_dE_onesided;

  /*! DEBUG: every transport increment of flux[].energy, added atomically
   * when it is made (as energy); equals debug_dE[2] plus the pending energy
   * unless increments to the accumulators were lost */
  float debug_dE_acc;

  /*! DEBUG: cumulative RT photoheating (internal energy units), filled only
   * with KIARART:debug_grackle_heat: 0 delivered by Grackle (m * (u with RT
   * heating - u without)), 1 intended (heating rate * start-of-step N_HI *
   * dt) */
  float debug_heat[2];
#endif

  /* Radiation state vector. */
  struct {
    /*! Radiation energy density */
    float energy_density;
    /*! Radiation flux */
    float flux[3];
  } radiation[RT_NGROUPS];

  /* Fluxes in the conservation law sense */
  struct {
    /*! Radiation energy flux */
    float energy;
    /*! Flux of the radiation flux (momentum flux) */
    float flux[3];
  } flux[RT_NGROUPS];

  /*! Interstellar radiation field from the FUV photon group, in Habing units
   * (0 below the density threshold, -1 if there is no FUV group). Used by
   * the cooling as G0 with G0_computation_method 8. */
  float G0;

  /* Particle RT time step. */
  float flux_dt;

  /*! Stellar emission received and not yet added to the radiation (total
   * energy): stars spread their emission over their time step, see
   * rt_inject_pending_emission() */
  float emission_pending[RT_NGROUPS];

  /*! Rate (energy per physical time) at which emission_pending is added to
   * the radiation, the sum of the rates of the stars it came from */
  float emission_rate[RT_NGROUPS];

  /*! Sum of the kernel weights of the neighbours eligible for transport,
   * for a particle that is not (decoupled or with its cooling shut off) and
   * gives its radiation to them (see rt_shed_collect()) */
  float shed_wsum;

  /*! Factor (<= 1) applied to the particle's reduced speed of light so that
   * its RT time step satisfies the Courant condition (1 unless
   * KIARART:local_cfl_reduced_c is set and the time step was floored) */
  float c_scale;

  /*! Particle (comoving) volume when the radiation state was last updated
   * for a change of volume (0 until the first volume computation) */
  float volume_previous;

  /*! Sum of the (comoving) interface areas |A_ij| with the neighbours,
   * accumulated during the transport since the last finalise */
  float area_sum;

  /*! Interface area sum of the last completed transport step, used for the
   * RT time-step (0 until the first transport) */
  float area_sum_last;

  /*! Sum of the (comoving) interface areas |A_ij| with the neighbours of the
   * current RT step, from the gradient loop (so it is complete before the
   * transport and the same on all nodes) */
  float area_sum_grad;

  /* gradients of the radiation state. */
  /* for the flux[3][3] quantity:
   *    first index: x, y, z coordinate of the flux.
   *    Second index: gradient along x, y, z direction. */
  struct {
    /*! Gradient of c times the energy density (the transported quantity, see
     * rt_part_get_transport_state_vector) */
    float energy_density[3];
    /*! Gradient of the flux components */
    float flux[3][3];
  } gradient[RT_NGROUPS];

  /* cell slope limiter quantities */
  /* array of length two: store min among all neighbours
   * at first index, store max among all neighbours at
   * second index */
  /* the Gizmo-style slope limiting doesn't help for RT as is,
   * so we're skipping it for now. */
  /* struct { */
  /*   float energy_density[2]; */
  /*   float flux[3][2]; */
  /*   [> float maxr; [> just use the hydro one <] <] */
  /* } limiter[RT_NGROUPS]; */

  /* Data for thermochemistry */
  struct {
    float mass_fraction_HI;         /* mass fraction taken by HI */
    float mass_fraction_HII;        /* mass fraction taken by HII */
    float mass_fraction_HeI;        /* mass fraction taken by HeI */
    float mass_fraction_HeII;       /* mass fraction taken by HeII */
    float mass_fraction_HeIII;      /* mass fraction taken by HeIII */
    float number_density_electrons; /* number density of electrons */
  } tchem;

#ifdef GIZMO_MFV_SPH
  /* Keep track of the actual mass fluxes of the gas species */
  struct {
    float HI;    /* mass fraction taken by HI */
    float HII;   /* mass fraction taken by HII */
    float HeI;   /* mass fraction taken by HeI */
    float HeII;  /* mass fraction taken by HeII */
    float HeIII; /* mass fraction taken by HeIII */
  } mass_flux;
#endif

#ifdef SWIFT_RT_DEBUG_CHECKS
  /* debugging data to store during entire run */

  /*! how much radiation this part received from stars during total lifetime */
  unsigned long long debug_radiation_absorbed_tot;

  /* data to store during one time step */

  /*! how many stars this part interacted with during injection*/
  /* Note: It's useless to write this in outputs, as it gets reset
   * at the end of every step. */
  int debug_iact_stars_inject;

  /*! calls from gradient interaction loop in actual function */
  int debug_calls_iact_gradient_interaction;

  /*! calls from transport interaction loop in actual function */
  int debug_calls_iact_transport_interaction;

  /* Task completion flags */

  /*! part got kicked? */
  int debug_kicked;

  /*! calls from ghost1 tasks */
  int debug_injection_done;

  /*! finalised computing gradients? */
  int debug_gradients_done;

  /*! transport step done? */
  int debug_transport_done;

  /*! thermochemistry done? */
  int debug_thermochem_done;

  /* Subcycling flags */

  /*! Current subcycle wrt (last) hydro step */
  int debug_nsubcycles;

#endif
};

/* Additional RT data in star particle struct */
struct rt_spart_data {

  /* Stellar energy emission that will be injected in to gas.
   * Total energy, not density, not rate! */
  float emission_this_step[RT_NGROUPS];

  /*! Physical time over which emission_this_step is emitted (the star's
   * time step): the gas neighbours add it to their radiation at the rate
   * emission_this_step / emission_dt (0: at once) */
  float emission_dt;

  /*! Neighbour weigths in each octant surrounding the star */
  float octant_weights[8];

  /*! Energy to add to the next emission (total energy): the emission of
   * steps without gas neighbours, and the radiation of the gas particle the
   * star formed from */
  float emission_carry[RT_NGROUPS];

#ifdef SWIFT_KIARA_RT_DEBUG
  /*! DEBUG: photon energy budget (rt_part_data.debug_dE) of the gas particle
   * the star formed from, frozen at the conversion, with the radiation
   * handed to the star moved out of the injected column (it is injected
   * again, and counted there, by the star) */
  float debug_dE[6];

  /*! DEBUG: rt_part_data.debug_dE_onesided of that gas particle */
  float debug_dE_onesided;

  /*! DEBUG: rt_part_data.debug_dE_acc of that gas particle */
  float debug_dE_acc;
#endif

#ifdef SWIFT_RT_DEBUG_CHECKS
  /* data to store during entire run */

  /*! how much radiation this star emitted during total lifetime */
  unsigned long long debug_radiation_emitted_tot;

  /* data to store during one time step */

  /*! how many hydro particles this particle interacted with
   * during injection */
  int debug_iact_hydro_inject;

  /*! how many hydro particles this particle interacted with
   * during injection prep*/
  int debug_iact_hydro_inject_prep;

  /*! stellar photon emisison rate computed? */
  int debug_emission_rate_set;

  /*! how much energy this star particle actually has injected into the gas */
  float debug_injected_energy[RT_NGROUPS];

  /*! how much energy this star particle actually has injected into the gas over
   * the entire run*/
  float debug_injected_energy_tot[RT_NGROUPS];

  /*! sum up total weights used during injection to compare consistency */
  float debug_psi_sum;
#endif
};

#endif /* SWIFT_RT_STRUCT_KIARA_H */
