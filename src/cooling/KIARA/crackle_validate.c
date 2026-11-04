/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 Romeel Dave (romeel.dave@ed.ac.uk)
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

/**
 * @file cooling/KIARA/crackle_validate.c
 * @brief Run the external Grackle/crackle library on the same inputs as the
 * built-in crackle and check the outputs agree bit for bit
 * (configure --enable-crackle-validate).
 *
 * This file includes the EXTERNAL grackle.h (found through --with-grackle)
 * and none of the built-in crackle headers, so the Grackle names used here
 * are those of the external library; the built-in ones are all crk_-prefixed.
 */

/* Config parameters. */
#include <config.h>

#ifdef CRACKLE_VALIDATE

/* Some standard headers. */
#include <stdint.h>
#include <string.h>

/* The external Grackle library */
#include <grackle.h>

/* Local includes. */
#include "crackle_validate.h"
#include "error.h"

/*! External chemistry instance used through the global API (KIARA cooling) */
static chemistry_data crackle_validate_chem_global;

/*! External chemistry instance used through the local_ API (KIARA-RT) */
static chemistry_data crackle_validate_chem_local;
static chemistry_data_storage crackle_validate_rates_local;

/**
 * @brief Initialise the external library with the parameters of the built-in
 * crackle.
 *
 * Call right before the built-in initialise, with the parameters set.
 *
 * @param api Which external instance (global or local_ API).
 * @param chemistry The built-in chemistry_data (parameters only).
 * @param units The built-in code_units.
 * @param chemistry_size sizeof(chemistry_data) of the built-in crackle.
 * @param units_size sizeof(code_units) of the built-in crackle.
 * @param field_size sizeof(grackle_field_data) of the built-in crackle.
 */
void crackle_validate_init(enum crackle_validate_api api,
                           const void *chemistry, const void *units,
                           size_t chemistry_size, size_t units_size,
                           size_t field_size) {

  /* The two libraries must share the struct layouts */
  if (chemistry_size != sizeof(chemistry_data) ||
      units_size != sizeof(code_units) ||
      field_size != sizeof(grackle_field_data))
    error(
        "crackle validation: struct sizes differ between the built-in crackle "
        "(chemistry_data %zu, code_units %zu, grackle_field_data %zu) and the "
        "external library (%zu, %zu, %zu)",
        chemistry_size, units_size, field_size, sizeof(chemistry_data),
        sizeof(code_units), sizeof(grackle_field_data));
  if (sizeof(code_units) > CRACKLE_VALIDATE_MAX_UNITS_BYTES ||
      sizeof(grackle_field_data) > CRACKLE_VALIDATE_MAX_FIELD_BYTES)
    error("crackle validation: increase CRACKLE_VALIDATE_MAX_*_BYTES");
  if (sizeof(gr_float) != sizeof(double))
    error("crackle validation: the external library must use double gr_float");

  code_units u;
  memcpy(&u, units, sizeof(code_units));

  if (api == crackle_validate_global) {
    /* Points the external global grackle_data at our instance */
    if (set_default_chemistry_parameters(&crackle_validate_chem_global) == 0)
      error("crackle validation: external set_default_chemistry_parameters");
    memcpy(&crackle_validate_chem_global, chemistry, sizeof(chemistry_data));
    if (initialize_chemistry_data(&u) == 0)
      error("crackle validation: external initialize_chemistry_data");
  } else {
    if (local_initialize_chemistry_parameters(&crackle_validate_chem_local) ==
        0)
      error("crackle validation: external local_initialize_chemistry_params");
    memcpy(&crackle_validate_chem_local, chemistry, sizeof(chemistry_data));
    if (local_initialize_chemistry_data(&crackle_validate_chem_local,
                                        &crackle_validate_rates_local,
                                        &u) == 0)
      error("crackle validation: external local_initialize_chemistry_data");
  }

  message("crackle validation: external library initialised (%s API)",
          api == crackle_validate_global ? "global" : "local");
}

/**
 * @brief Save the inputs of a built-in crackle call.
 *
 * @param s The state to fill.
 * @param units The code_units passed to the call.
 * @param field The grackle_field_data passed to the call.
 * @param species The array the field pointers point into.
 * @param nspecies Length of species.
 * @param iact The RT rate array the field pointers point into (or NULL).
 * @param niact Length of iact.
 */
void crackle_validate_save(struct crackle_validate_state *s, const void *units,
                           const void *field, const double *species,
                           int nspecies, const double *iact, int niact) {

  if (nspecies > CRACKLE_VALIDATE_MAX_SPECIES ||
      niact > CRACKLE_VALIDATE_MAX_IACT)
    error("crackle validation: increase CRACKLE_VALIDATE_MAX_SPECIES/IACT");

  memcpy(s->units, units, sizeof(code_units));
  memcpy(s->field, field, sizeof(grackle_field_data));
  memcpy(s->species, species, nspecies * sizeof(double));
  if (iact != NULL) memcpy(s->iact, iact, niact * sizeof(double));
  s->species_orig = species;
  s->iact_orig = iact;
  s->nspecies = nspecies;
  s->niact = (iact != NULL) ? niact : 0;
}

/**
 * @brief Redirect a pointer into an original array to the same element of
 * its copy.
 */
static void crackle_validate_rebase(void **ptr, const double *orig, int n,
                                    double *copy) {
  if (orig == NULL || n == 0) return;
  const uintptr_t v = (uintptr_t)(*ptr);
  const uintptr_t lo = (uintptr_t)orig;
  const uintptr_t hi = (uintptr_t)(orig + n);
  if (v >= lo && v < hi) *ptr = (void *)(copy + (v - lo) / sizeof(double));
}

/**
 * @brief Report the differences between two arrays, if any.
 *
 * @return The number of elements that differ (bitwise).
 */
static int crackle_validate_compare(const char *what, const double *ref,
                                    const double *ext, int n, long long id,
                                    int mode) {
  int ndiff = 0;
  for (int i = 0; i < n; i++) {
    if (memcmp(&ref[i], &ext[i], sizeof(double)) != 0) {
      message(
          "crackle validation: id=%lld mode=%d %s[%d]: built-in %.17g "
          "external %.17g",
          id, mode, what, i, ref[i], ext[i]);
      ndiff++;
    }
  }
  return ndiff;
}

/**
 * @brief Run the external library on the saved inputs and compare its
 * outputs with those of the built-in crackle.
 *
 * Stops the run at the first difference.
 *
 * @param s The saved inputs.
 * @param api Which external instance (global or local_ API).
 * @param mode 0 solve_chemistry, 1 cooling time, 2 temperature, 3 pressure,
 * 4 gamma.
 * @param dt The time-step passed to solve_chemistry.
 * @param species_after The species array after the built-in call.
 * @param iact_after The RT rate array after the built-in call (or NULL).
 * @param result The value returned by the built-in call (modes 1-4).
 * @param id The particle ID (for the report).
 */
void crackle_validate_check(struct crackle_validate_state *s,
                            enum crackle_validate_api api, int mode, double dt,
                            const double *species_after,
                            const double *iact_after, double result,
                            long long id) {

  /* Point the saved field struct at the saved arrays */
  grackle_field_data *f = (grackle_field_data *)s->field;
  void **words = (void **)s->field;
  const size_t nwords = sizeof(grackle_field_data) / sizeof(void *);
  for (size_t k = 0; k < nwords; k++) {
    crackle_validate_rebase(&words[k], s->species_orig, s->nspecies,
                            s->species);
    crackle_validate_rebase(&words[k], s->iact_orig, s->niact, s->iact);
  }
  code_units *u = (code_units *)s->units;

  /* Run the external library */
  gr_float ext_result = 0.;
  int ok = 1;
  if (api == crackle_validate_global) {
    switch (mode) {
      case 0:
        ok = solve_chemistry(u, f, dt);
        break;
      case 1:
        ok = calculate_cooling_time(u, f, &ext_result);
        break;
      case 2:
        ok = calculate_temperature(u, f, &ext_result);
        break;
      case 3:
        ok = calculate_pressure(u, f, &ext_result);
        break;
      case 4:
        ok = calculate_gamma(u, f, &ext_result);
        break;
      default:
        error("crackle validation: unknown mode %d", mode);
    }
  } else {
    switch (mode) {
      case 0:
        ok = local_solve_chemistry(&crackle_validate_chem_local,
                                   &crackle_validate_rates_local, u, f, dt);
        break;
      case 1:
        ok = local_calculate_cooling_time(&crackle_validate_chem_local,
                                          &crackle_validate_rates_local, u, f,
                                          &ext_result);
        break;
      default:
        error("crackle validation: unknown local mode %d", mode);
    }
  }
  if (ok == 0)
    error("crackle validation: external library failed (id=%lld mode=%d)", id,
          mode);

  /* Compare everything the call can write */
  int ndiff = crackle_validate_compare("species", species_after, s->species,
                                       s->nspecies, id, mode);
  if (s->niact > 0)
    ndiff += crackle_validate_compare("iact", iact_after, s->iact, s->niact,
                                      id, mode);
  if (mode != 0) {
    const double ext = ext_result;
    ndiff += crackle_validate_compare("result", &result, &ext, 1, id, mode);
  }
  if (ndiff > 0)
    error(
        "crackle validation: the built-in crackle differs from the external "
        "library in %d value(s) (id=%lld mode=%d api=%d dt=%.17g)",
        ndiff, id, mode, api, dt);
}

#endif /* CRACKLE_VALIDATE */
