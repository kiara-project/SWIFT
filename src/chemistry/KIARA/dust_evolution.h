/* Physical radius transport on arbitrary KIARA dust grids. */
#ifndef SWIFT_KIARA_DUST_EVOLUTION_H
#define SWIFT_KIARA_DUST_EVOLUTION_H

#include "dust_grid.h"

#if KIARA_DUST_N_BINS > 0

/** Remap a' = scale*a + shift, with optional fixed-number mass evolution
 * proportional to a'^3. */
static inline double kiara_dust_radius_remap(
    const double radius[KIARA_DUST_N_BINS],
    const double input[KIARA_DUST_N_BINS], const double shift,
    const double scale, const int evolve_mass,
    double output[KIARA_DUST_N_BINS]) {
  if (!isfinite(shift) || !isfinite(scale) || scale < 0.) return NAN;
  for (int i = 0; i < KIARA_DUST_N_BINS; ++i) output[i] = 0.;
  for (int i = 0; i < KIARA_DUST_N_BINS; ++i) {
    if (!isfinite(input[i]) || input[i] < 0.) return NAN;
    if (input[i] == 0.) continue;
    const double evolved_radius = fmax(0., scale * radius[i] + shift);
    if (!isfinite(evolved_radius)) return NAN;
    if (evolve_mass && evolved_radius == 0.) continue;
    const double ratio = evolved_radius / radius[i];
    const double mass =
        evolve_mass ? input[i] * ratio * ratio * ratio : input[i];
    if (!isfinite(mass)) return NAN;
    /* Only physical erosion to zero destroys mass; positive unresolved grains
     * stay in endpoint bins. */
    kiara_dust_grid_deposit(radius, fmax(evolved_radius, radius[0]), mass,
                            output);
  }
  double sum = 0.;
  for (int i = 0; i < KIARA_DUST_N_BINS; ++i) sum += output[i];
  return sum;
}

/**
 * @brief Grain surface area per dust mass of a material-by-size distribution,
 * relative to grains of radius a_ref: a_ref * sum(m / a) / sum(m).
 *
 * Every bin is taken to hold grains of its representative radius, whose area
 * per mass is proportional to 1 / a. Returns 1 for an empty distribution.
 *
 * @param radius The bin radii.
 * @param distribution The masses (or mass fractions) in each material and
 * bin.
 * @param a_ref The reference radius (same units as radius).
 */
static inline double kiara_dust_area_factor(
    const double radius[KIARA_DUST_N_BINS],
    const double distribution[KIARA_DUST_N_MATERIALS][KIARA_DUST_N_BINS],
    const double a_ref) {
  double mass = 0., moment = 0.;
  for (int material = 0; material < KIARA_DUST_N_MATERIALS; ++material)
    for (int i = 0; i < KIARA_DUST_N_BINS; ++i) {
      mass += distribution[material][i];
      moment += distribution[material][i] / radius[i];
    }
  return (mass > 0. && moment > 0.) ? a_ref * moment / mass : 1.;
}

/**
 * @brief The common radius increase da >= 0 of all grains that adds the mass
 * dm to a size distribution: sum_i m_i ((1 + da / a_i)^3 - 1) = dm.
 *
 * Accretion growth in Li et al. (2021) has da/dt independent of the radius,
 * so all grains grow by the same da over a step. The left-hand side is convex
 * and increasing in da, so Newton's method from the linear estimate
 * converges.
 *
 * @param radius The bin radii.
 * @param mass The mass in each bin.
 * @param dm The mass to add (>= 0).
 * @return da (in the units of radius), 0 if there is nothing to grow.
 */
static inline double kiara_dust_growth_shift(
    const double radius[KIARA_DUST_N_BINS],
    const double mass[KIARA_DUST_N_BINS], const double dm) {
  double moment = 0.;
  for (int i = 0; i < KIARA_DUST_N_BINS; ++i) moment += mass[i] / radius[i];
  if (!(dm > 0.) || !(moment > 0.)) return 0.;
  double da = dm / (3. * moment);
  for (int iter = 0; iter < 100; ++iter) {
    double f = -dm, df = 0.;
    for (int i = 0; i < KIARA_DUST_N_BINS; ++i) {
      const double x = 1. + da / radius[i];
      f += mass[i] * (x * x * x - 1.);
      df += 3. * mass[i] * x * x / radius[i];
    }
    const double step = f / df;
    da -= step;
    if (da < 0.) da = 0.;
    if (fabs(step) <= 1.e-12 * da) break;
  }
  return da;
}

#endif /* KIARA_DUST_N_BINS > 0 */

#endif
