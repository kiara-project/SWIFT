/*******************************************************************************
 * Li et al. (2021) grain coagulation on KIARA's existing size grid.
 ******************************************************************************/
#ifndef SWIFT_KIARA_DUST_COAGULATION_H
#define SWIFT_KIARA_DUST_COAGULATION_H

#include "dust_collisions.h"

/**
 * @brief Sticking threshold in cm/s: Hirashita & Yan (2009), equation 8.
 *
 * Li et al. (2021), section 2.6, adopts this threshold, including F_stick=10.
 * Use graphite/quartz surface energy, Young's modulus and Poisson ratio from
 * Hirashita & Yan table 2 for carbonaceous/silicate material. Their equation
 * for the reduced elastic modulus is 1/E = 2(1-nu)^2/E_material for identical
 * species. KIARA retains Li's common configurable material density (default 2.4
 * g/cm^3). Radii are cm, surface energies erg/cm^2, and elastic moduli
 * dyn/cm^2. This material-dependent threshold extends Li's shared size
 * distribution; collisions between different materials are not represented.
 */
static inline double kiara_dust_coagulation_threshold(
    const double radius_i, const double radius_j, const double grain_density,
    const enum kiara_dust_material material) {

  assert(radius_i > 0. && radius_j > 0. && grain_density > 0.);
  assert(material == kiara_dust_carbonaceous ||
         material == kiara_dust_silicate);
  const int carbon = material == kiara_dust_carbonaceous;
  const double surface_energy = carbon ? 12. : 25.;
  const double young = carbon ? 3.4e10 : 5.4e11;
  const double poisson = carbon ? 0.5 : 0.17;
  const double reduced_modulus = young / (2. * (1. - poisson) * (1. - poisson));
  const double sum = radius_i + radius_j;
  const double x = radius_i / sum, y = radius_j / sum;
  const double reduced_radius = radius_i * y;
  return 21.4 * sqrt(x * x * x + y * y * y) * pow(surface_energy, 5. / 6.) /
         (cbrt(reduced_modulus) * pow(reduced_radius, 5. / 6.) *
          sqrt(grain_density));
}

/**
 * @brief Build the coagulation operator for one material and fixed velocities.
 *
 * The caller supplies Li's Gaussian-drawn bin velocities, in cm/s. Following
 * Li, grains within one bin share a velocity, so same-bin collisions have zero
 * rate. No RMS-speed substitution or density/temperature switch overrides the
 * sticking threshold. Below threshold, the sticking probability is unity.
 * Representative-radius cross-sections and conservative interpolation of the
 * product mass adapt Li's equations 16, 18 and 21 to the configured mass
 * pivots.
 */
static inline void kiara_dust_coagulation_operator_init(
    struct kiara_dust_collision_operator *op,
    const struct kiara_dust_collision_grid *grid, const double grain_density,
    const enum kiara_dust_material material,
    const double velocity[KIARA_DUST_N_BINS][3]) {

  const double pi = 3.14159265358979323846;
  for (int i = 0; i < KIARA_DUST_N_BINS; ++i) {
    for (int j = i; j < KIARA_DUST_N_BINS; ++j) {
      op->kernel[i][j] = 0.;
      for (int k = 0; k < KIARA_DUST_N_BINS; ++k) op->product[i][j][k] = 0.;
      if (i == j) continue;
      double speed_squared = 0.;
      for (int axis = 0; axis < 3; ++axis) {
        const double dv = velocity[i][axis] - velocity[j][axis];
        speed_squared += dv * dv;
      }
      const double speed = sqrt(speed_squared);
      const double threshold = kiara_dust_coagulation_threshold(
          grid->radius[i], grid->radius[j], grain_density, material);
      if (speed == 0. || speed >= threshold) continue;
      const double radius_sum = grid->radius[i] + grid->radius[j];
      op->kernel[i][j] = pi * radius_sum * radius_sum * speed;
      const double product_mass = grid->mass[i] + grid->mass[j];
      /* Li Eq. 21: merge masses, interpolating onto the configured pivots. */
      kiara_dust_collision_deposit(grid, product_mass, product_mass,
                                   op->product[i][j]);
    }
  }
}

#endif /* SWIFT_KIARA_DUST_COAGULATION_H */
