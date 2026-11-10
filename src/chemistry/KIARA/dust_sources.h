/* Physical stellar grain-size spectra on the configured KIARA radius grid. */
#ifndef SWIFT_KIARA_DUST_SOURCES_H
#define SWIFT_KIARA_DUST_SOURCES_H

#include "dust_grid.h"

#if KIARA_DUST_N_BINS > 0

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

/**
 * @brief The dust mass in each material and size bin of a particle.
 *
 * Each material's row of the distribution (fractions of the dust mass) is
 * rescaled to the material's share of the dust, from the element composition,
 * so that the masses agree with it even if the rows are out of date. An empty
 * row of a material with dust gets the shape seed (if given; else it stays
 * empty).
 *
 * @param distribution The distribution (fractions of the dust mass).
 * @param dust The dust mass.
 * @param f_carbon The carbonaceous share of the dust (the rest is silicate).
 * @param seed A normalised shape for empty rows, or NULL.
 * @param mass (return) The masses.
 */
static inline void kiara_dust_bin_masses(
    const float distribution[KIARA_DUST_N_MATERIALS][KIARA_DUST_N_BINS],
    const double dust, const double f_carbon,
    const double seed[KIARA_DUST_N_BINS],
    double mass[KIARA_DUST_N_MATERIALS][KIARA_DUST_N_BINS]) {
  const double f_C = fmin(fmax(f_carbon, 0.), 1.);
  for (int material = 0; material < KIARA_DUST_N_MATERIALS; ++material) {
    const double m_material =
        (dust > 0.) * dust *
        (material == kiara_dust_carbonaceous ? f_C : 1. - f_C);
    double row = 0.;
    for (int i = 0; i < KIARA_DUST_N_BINS; ++i)
      row += distribution[material][i];
    for (int i = 0; i < KIARA_DUST_N_BINS; ++i)
      mass[material][i] = (row > 0.)
                              ? distribution[material][i] / row * m_material
                              : ((seed != NULL) ? seed[i] * m_material : 0.);
  }
}

/**
 * @brief Set a particle's size distribution from dust masses in its bins.
 *
 * The bin masses give each material's shape; the rows are rescaled to the
 * materials' dust masses m_material (from the element composition) and
 * divided by their sum, so that the distribution always agrees with the
 * element masses. An empty row of a material with dust gets the shape seed
 * (if given; else it stays empty, to be seeded later).
 *
 * @param distribution (return) The distribution (fractions of the dust mass).
 * @param mass The masses in each material and bin (only their shapes count).
 * @param m_material The dust mass of each material.
 * @param seed A normalised shape for empty rows, or NULL.
 */
static inline void kiara_dust_set_distribution(
    float distribution[KIARA_DUST_N_MATERIALS][KIARA_DUST_N_BINS],
    const double mass[KIARA_DUST_N_MATERIALS][KIARA_DUST_N_BINS],
    const double m_material[KIARA_DUST_N_MATERIALS],
    const double seed[KIARA_DUST_N_BINS]) {
  double total = 0.;
  for (int material = 0; material < KIARA_DUST_N_MATERIALS; ++material)
    total += fmax(m_material[material], 0.);
  for (int material = 0; material < KIARA_DUST_N_MATERIALS; ++material) {
    const double m = fmax(m_material[material], 0.);
    double row = 0.;
    for (int i = 0; i < KIARA_DUST_N_BINS; ++i)
      row += fmax(mass[material][i], 0.);
    for (int i = 0; i < KIARA_DUST_N_BINS; ++i) {
      double f = 0.;
      if (total > 0. && m > 0.) {
        if (row > 0.)
          f = fmax(mass[material][i], 0.) / row * m / total;
        else if (seed != NULL)
          f = seed[i] * m / total;
      }
      distribution[material][i] = (float)f;
    }
  }
}

#endif /* KIARA_DUST_N_BINS > 0 */

#endif
