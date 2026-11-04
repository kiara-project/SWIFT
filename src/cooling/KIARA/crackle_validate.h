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
#ifndef SWIFT_COOLING_KIARA_CRACKLE_VALIDATE_H
#define SWIFT_COOLING_KIARA_CRACKLE_VALIDATE_H

/**
 * @file cooling/KIARA/crackle_validate.h
 * @brief Bit-for-bit comparison of the built-in crackle with an external
 * Grackle/crackle library (configure --enable-crackle-validate).
 *
 * Before each call to the built-in crackle, the inputs (units, field struct
 * and the arrays it points into) are saved; after it, the external library
 * is run on the saved copy and every output is compared bit for bit with the
 * built-in one. The run stops at the first difference.
 *
 * This header uses no Grackle types, as it is included both by SWIFT code
 * (with the built-in crackle headers) and by crackle_validate.c (with the
 * external Grackle headers).
 */

#include <stddef.h>

#ifdef CRACKLE_VALIDATE

/*! Maximal sizes of the saved copies (checked at initialisation) */
#define CRACKLE_VALIDATE_MAX_FIELD_BYTES 4096
#define CRACKLE_VALIDATE_MAX_UNITS_BYTES 512
#define CRACKLE_VALIDATE_MAX_SPECIES 64
#define CRACKLE_VALIDATE_MAX_IACT 8

/*! Which external chemistry instance to run */
enum crackle_validate_api {
  crackle_validate_global = 0, /*!< Grackle's global API (KIARA cooling) */
  crackle_validate_local = 1   /*!< Grackle's local_ API (KIARA-RT) */
};

/*! The inputs of one crackle call, saved for the external library */
struct crackle_validate_state {

  /*! Copy of the grackle_field_data struct */
  unsigned char field[CRACKLE_VALIDATE_MAX_FIELD_BYTES];

  /*! Copy of the code_units struct */
  unsigned char units[CRACKLE_VALIDATE_MAX_UNITS_BYTES];

  /*! Copies of the arrays the field struct points into */
  double species[CRACKLE_VALIDATE_MAX_SPECIES];
  double iact[CRACKLE_VALIDATE_MAX_IACT];

  /*! The original arrays (to redirect the field pointers into the copies) */
  const double *species_orig;
  const double *iact_orig;
  int nspecies;
  int niact;
};

void crackle_validate_init(enum crackle_validate_api api,
                           const void *chemistry, const void *units,
                           size_t chemistry_size, size_t units_size,
                           size_t field_size);

void crackle_validate_save(struct crackle_validate_state *s, const void *units,
                           const void *field, const double *species,
                           int nspecies, const double *iact, int niact);

void crackle_validate_check(struct crackle_validate_state *s,
                            enum crackle_validate_api api, int mode, double dt,
                            const double *species_after,
                            const double *iact_after, double result,
                            long long id);

#endif /* CRACKLE_VALIDATE */

#endif /* SWIFT_COOLING_KIARA_CRACKLE_VALIDATE_H */
