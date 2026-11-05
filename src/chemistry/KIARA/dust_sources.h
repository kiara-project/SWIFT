/* Physical stellar grain-size spectra on the configured KIARA radius grid. */
#ifndef SWIFT_KIARA_DUST_SOURCES_H
#define SWIFT_KIARA_DUST_SOURCES_H

#include "dust_grid.h"

/** Integrate Li et al. (2021), Eq. 5, as mass in geometric cells, folding both
 * tails into endpoint bins. */
static inline int kiara_dust_source_lognormal(
    const double radius[KIARA_DUST_N_BINS], const double a0, const double sigma,
    const double p, double mass[KIARA_DUST_N_BINS]) {
  if (!kiara_dust_grid_valid(radius) || !isfinite(a0) || a0 <= 0. ||
      !isfinite(sigma) || sigma <= 0. || !isfinite(p))
    return 0;
  const double mean = log(a0) + (3. - p) * sigma * sigma;
  const double width = sqrt(2.) * sigma;
  if (!isfinite(mean) || !isfinite(width)) return 0;
  double sum = 0.;
  for (int i = 0; i < KIARA_DUST_N_BINS; ++i) {
    const double lo =
        i == 0 ? -INFINITY
               : (0.5 * (log(radius[i - 1]) + log(radius[i])) - mean) / width;
    const double hi =
        i == KIARA_DUST_N_BINS - 1
            ? INFINITY
            : (0.5 * (log(radius[i]) + log(radius[i + 1])) - mean) / width;
    mass[i] = 0.5 * (lo >= 0. ? erfc(lo) - erfc(hi) : erfc(-hi) - erfc(-lo));
    sum += mass[i];
  }
  if (!isfinite(sum) || sum <= 0.) return 0;
  for (int i = 0; i < KIARA_DUST_N_BINS; ++i) mass[i] /= sum;
  return 1;
}

/** Seed each material from a normalized physical source spectrum, retaining its
 * own dust mass. */
static inline void kiara_dust_distribution_seed(
    const double shape[KIARA_DUST_N_BINS], const float carbon_fraction,
    float distribution[KIARA_DUST_N_MATERIALS][KIARA_DUST_N_BINS]) {
  for (int material = 0; material < KIARA_DUST_N_MATERIALS; ++material) {
    const double fraction = material == kiara_dust_carbonaceous
                                ? carbon_fraction
                                : 1. - carbon_fraction;
    for (int i = 0; i < KIARA_DUST_N_BINS; ++i)
      distribution[material][i] = fraction * shape[i];
  }
}

#endif
