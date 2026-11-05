/* Shared runtime grain-radius parameters for KIARA cooling and feedback. */
#ifndef SWIFT_KIARA_DUST_PARAMETERS_H
#define SWIFT_KIARA_DUST_PARAMETERS_H

#include "dust_grid.h"
#include "dust_sources.h"
#include "error.h"
#include "parser.h"

#include <stdio.h>
#include <string.h>

/** Read one common grid, retaining the legacy endpoint parameter aliases. */
static inline void kiara_dust_read_grid(struct swift_params *params,
                                        double radius[KIARA_DUST_N_BINS]) {
  const double legacy_min = parser_get_opt_param_double(
      params, "KIARACooling:dust_small_grainsize", 0.005);
  const double legacy_max = parser_get_opt_param_double(
      params, "KIARACooling:dust_large_grainsize", 0.1);
  const double minimum = parser_get_opt_param_double(
      params, "KIARACooling:dust_grain_size_min_micron", legacy_min);
  const double maximum = parser_get_opt_param_double(
      params, "KIARACooling:dust_grain_size_max_micron", legacy_max);
  if (parser_does_param_exist(params, "KIARACooling:dust_grain_sizes_micron")) {
    parser_get_opt_param_double_array(params,
                                      "KIARACooling:dust_grain_sizes_micron",
                                      KIARA_DUST_N_BINS, radius);
  } else if (!kiara_dust_grid_logarithmic(minimum, maximum, radius)) {
    error("KIARA dust grid requires finite radii with 0 < minimum < maximum");
  }
  if (!kiara_dust_grid_valid(radius))
    error(
        "KIARACooling:dust_grain_sizes_micron must be finite, positive and "
        "strictly increasing");
}

/** Read Li Eq. 5 source parameters once per consumer, using the same defaults
 * for injection and initialization. */
static inline void kiara_dust_read_source(
    struct swift_params *params, const int agb,
    const double radius[KIARA_DUST_N_BINS], double shape[KIARA_DUST_N_BINS]) {
  char key[80];
  const char *channel = agb ? "AGB" : "SNII";
  snprintf(key, sizeof(key), "KIARAFeedback:dust_%s_a0_micron", channel);
  const double a0 = parser_get_opt_param_double(params, key, 0.1);
  snprintf(key, sizeof(key), "KIARAFeedback:dust_%s_sigma", channel);
  const double sigma =
      parser_get_opt_param_double(params, key, agb ? 0.47 : 0.6);
  snprintf(key, sizeof(key), "KIARAFeedback:dust_%s_p", channel);
  const double p = parser_get_opt_param_double(params, key, agb ? 4. : 0.);
  if (!kiara_dust_source_lognormal(radius, a0, sigma, p, shape))
    error(
        "KIARA dust source spectra require finite a0 > 0, sigma > 0 and finite "
        "p");
}

/** Unknown-history and self-enrichment dust default to the SNII spectrum, with
 * an explicit AGB override. */
static inline int kiara_dust_read_initial_source(
    struct swift_params *params, const double radius[KIARA_DUST_N_BINS],
    double shape[KIARA_DUST_N_BINS]) {
  char channel[PARSER_MAX_LINE_SIZE];
  parser_get_opt_param_string(params, "KIARACooling:dust_initial_source",
                              channel, "SNII");
  if (strcmp(channel, "SNII") && strcmp(channel, "AGB"))
    error("KIARACooling:dust_initial_source must be SNII or AGB");
  const int agb = strcmp(channel, "AGB") == 0;
  kiara_dust_read_source(params, agb, radius, shape);
  return agb;
}

#endif
