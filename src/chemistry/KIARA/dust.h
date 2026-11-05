/*******************************************************************************
 * Compile-time layout and small helpers for the KIARA binned dust model.
 ******************************************************************************/
#ifndef SWIFT_KIARA_DUST_H
#define SWIFT_KIARA_DUST_H

/* config.h normally provides this.  Keep headers usable by stand-alone tools. */
#ifndef KIARA_DUST_N_BINS
#define KIARA_DUST_N_BINS 6
#endif
#if KIARA_DUST_N_BINS < 2
#error "KIARA dust requires at least two size bins per material"
#endif

#define KIARA_DUST_N_MATERIALS 2
#define KIARA_DUST_N_DISTRIBUTION_VALUES \
  (KIARA_DUST_N_MATERIALS * KIARA_DUST_N_BINS)

enum kiara_dust_material {
  kiara_dust_carbonaceous = 0,
  kiara_dust_silicate = 1,
  kiara_dust_material_count = KIARA_DUST_N_MATERIALS
};

/** Zero a material-by-size distribution. */
static inline void
kiara_dust_distribution_zero(
    float distribution[KIARA_DUST_N_MATERIALS][KIARA_DUST_N_BINS]) {
  for (int material = 0; material < KIARA_DUST_N_MATERIALS; ++material)
    for (int bin = 0; bin < KIARA_DUST_N_BINS; ++bin)
      distribution[material][bin] = 0.f;
}

/** Return the sum of all material-by-size mass fractions. */
static inline float
kiara_dust_distribution_sum(
    const float distribution[KIARA_DUST_N_MATERIALS][KIARA_DUST_N_BINS]) {
  float sum = 0.f;
  for (int material = 0; material < KIARA_DUST_N_MATERIALS; ++material)
    for (int bin = 0; bin < KIARA_DUST_N_BINS; ++bin)
      sum += distribution[material][bin];
  return sum;
}

/** Rescale a material's size shape, using the supplied physical seed for an empty row. */
static inline void kiara_dust_distribution_rescale_material(
    float distribution[KIARA_DUST_N_BINS], const float material_fraction,
    const double seed[KIARA_DUST_N_BINS]) {
  double sum = 0.;
  for (int bin = 0; bin < KIARA_DUST_N_BINS; ++bin) sum += distribution[bin];
  const double scale = sum > 0. ? material_fraction / sum : 0.;
  for (int bin = 0; bin < KIARA_DUST_N_BINS; ++bin)
    distribution[bin] *= scale;
  if (sum == 0.)
    for (int bin = 0; bin < KIARA_DUST_N_BINS; ++bin)
      distribution[bin] = material_fraction * seed[bin];
}

/** Clamp negative entries and normalize a non-empty distribution to unity. */
static inline void
kiara_dust_distribution_normalize(
    float distribution[KIARA_DUST_N_MATERIALS][KIARA_DUST_N_BINS]) {
  float sum = 0.f;
  for (int material = 0; material < KIARA_DUST_N_MATERIALS; ++material) {
    for (int bin = 0; bin < KIARA_DUST_N_BINS; ++bin) {
      if (distribution[material][bin] < 0.f)
        distribution[material][bin] = 0.f;
      sum += distribution[material][bin];
    }
  }
  if (sum > 0.f) {
    const float sum_inv = 1.f / sum;
    for (int material = 0; material < KIARA_DUST_N_MATERIALS; ++material)
      for (int bin = 0; bin < KIARA_DUST_N_BINS; ++bin)
        distribution[material][bin] *= sum_inv;
  }
}

#endif /* SWIFT_KIARA_DUST_H */
