/* Grid geometry and conservative deposition for KIARA dust. */
#ifndef SWIFT_KIARA_DUST_GRID_H
#define SWIFT_KIARA_DUST_GRID_H

#include "dust.h"

#include <assert.h>
#include <math.h>

/** Validate representative grain radii in any consistent length unit. */
static inline int kiara_dust_grid_valid(
    const double radius[KIARA_DUST_N_BINS]) {
  for (int i = 0; i < KIARA_DUST_N_BINS; ++i)
    if (!isfinite(radius[i]) || radius[i] <= 0. ||
        (i > 0 && radius[i] <= radius[i - 1]))
      return 0;
  return 1;
}

/** Construct logarithmic radius pivots, including both supplied endpoints. */
static inline int kiara_dust_grid_logarithmic(
    const double minimum, const double maximum,
    double radius[KIARA_DUST_N_BINS]) {
  if (!isfinite(minimum) || !isfinite(maximum) || minimum <= 0. ||
      maximum <= minimum)
    return 0;
  for (int i = 0; i < KIARA_DUST_N_BINS; ++i) {
    const double x = (double)i / (KIARA_DUST_N_BINS - 1);
    radius[i] = exp((1. - x) * log(minimum) + x * log(maximum));
  }
  radius[0] = minimum;
  radius[KIARA_DUST_N_BINS - 1] = maximum;
  return kiara_dust_grid_valid(radius);
}

/** Interpolate grain number between mass pivots; clamp out-of-range mass to the
 * endpoint bins. */
static inline void kiara_dust_grid_deposit(
    const double radius[KIARA_DUST_N_BINS], const double grain_radius,
    const double total_mass, double bins[KIARA_DUST_N_BINS]) {
  assert(isfinite(grain_radius) && grain_radius > 0.);
  assert(isfinite(total_mass) && total_mass >= 0.);
  if (grain_radius <= radius[0]) {
    bins[0] += total_mass;
    return;
  }
  const int last = KIARA_DUST_N_BINS - 1;
  if (grain_radius >= radius[last]) {
    bins[last] += total_mass;
    return;
  }
  for (int hi = 1; hi < KIARA_DUST_N_BINS; ++hi) {
    if (grain_radius <= radius[hi]) {
      const double x = radius[hi - 1] / radius[hi];
      const double y = radius[hi - 1] / grain_radius;
      const double z = grain_radius / radius[hi];
      const double lo_fraction =
          y * y * y * (-expm1(3. * log(z))) / (-expm1(3. * log(x)));
      const double lo_mass = total_mass * fmin(fmax(lo_fraction, 0.), 1.);
      bins[hi - 1] += lo_mass;
      bins[hi] += total_mass - lo_mass;
      return;
    }
  }
}

#endif
