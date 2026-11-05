/* Physical radius transport on arbitrary KIARA dust grids. */
#ifndef SWIFT_KIARA_DUST_EVOLUTION_H
#define SWIFT_KIARA_DUST_EVOLUTION_H

#include "dust_grid.h"

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

#endif
