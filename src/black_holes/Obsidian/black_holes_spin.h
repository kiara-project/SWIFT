/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2022 Filip Husko (filip.husko@durham.ac.uk)
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

#ifndef SWIFT_SPIN_JET_BLACK_HOLES_SPIN_H
#define SWIFT_SPIN_JET_BLACK_HOLES_SPIN_H

/* Standard headers */
#include <float.h>

/* Local includes */
#include "black_holes_properties.h"
#include "black_holes_struct.h"
#include "hydro_properties.h"
#include "inline.h"
#include "physical_constants.h"

/**
 * @brief Compute the gravitational radius of a black hole.
 *
 * @param a Black hole mass.
 * @param constants Physical constants (in internal units).
 */
__attribute__((always_inline)) INLINE static float
black_hole_gravitational_radius(float mass,
                                const struct phys_const* constants) {

  const float r_G =
      mass * constants->const_newton_G /
      (constants->const_speed_light_c * constants->const_speed_light_c);

#ifdef SWIFT_DEBUG_CHECKS
  if (r_G <= 0.f) {
    error(
        "Something went wrong with calculation of R_G of black holes. "
        " R_G is %f instead of R_G > 0.",
        r_G);
  }
#endif

  return r_G;
}

/**
 * @brief Compute the radius of the horizon of a BH particle in gravitational
 * units.
 *
 * @param a Black hole spin, -1 < a < 1.
 */
__attribute__((always_inline)) INLINE static float black_hole_horizon_radius(
    float a) {
  return 1.f + sqrtf((1.f - a) * (1.f + a));
}

/**
 * @brief Compute the radius of the innermost stable circular orbit of a
 * BH particle in gravitational units.
 *
 * The expression is given in Appendix B of Fiacconi et al. (2018) or eqn. 4 in
 * Griffin et al. (2019).
 *
 * @param a Black hole spin, -1 < a < 1.
 */
__attribute__((always_inline)) INLINE static float black_hole_isco_radius(
    float a) {
  const float Z1 = 1.f + (cbrtf((1.f + fabsf(a)) * (1.f - a * a)) +
                          cbrtf((1.f - fabsf(a)) * (1.f - a * a)));
  const float Z2 = sqrtf(3.f * a * a + Z1 * Z1);

  const float R_ISCO =
      3. + Z2 - a / fabsf(a) * sqrtf((3.f - Z1) * (3.f + Z1 + 2.f * Z2));

#ifdef SWIFT_DEBUG_CHECKS
  if (Z1 > 3.f) {
    error(
        "Something went wrong with calculation of Z1 factor for r_isco of"
        " black holes. Z1 is %f instead of Z1 > 3.",
        Z1);
  }

  if ((3.f + Z1 + 2.f * Z2) < 0.f) {
    error(
        "Something went wrong with calculation of (3. + Z1 + 2. * Z2 ) "
        "factor for r_isco of black holes. (3. + Z1 + 2. * Z2 ) is %f instead "
        "of"
        " (3. + Z1 + 2. * Z2 ) > 0.",
        3.f + Z1 + 2.f * Z2);
  }

  if (R_ISCO < 1.f) {
    error(
        "Something went wrong with calculation of R_ISCO of black holes. "
        "R_ISCO is %f instead >= 1.",
        R_ISCO);
  }
#endif

  return R_ISCO;
}

/**
 * @brief Compute the magnitude of the angular momentum of the black hole
 * given its spin.
 *
 * @param a Black hole spin magnitude, 0 < a < 1.
 * @param constants Physical constants (in internal units).
 */
__attribute__((always_inline)) INLINE static float
black_hole_angular_momentum_magnitude(struct bpart* bp,
                                      const struct phys_const* constants) {

  const float J_BH =
      fabs(bp->subgrid_mass * bp->subgrid_mass * bp->spin *
           constants->const_newton_G / constants->const_speed_light_c);

#ifdef SWIFT_DEBUG_CHECKS
  if (J_BH <= 0.f) {
    error(
        "Something went wrong with calculation of j_BH of black holes. "
        " J_BH is %f instead of J_BH > 0.",
        J_BH);
  }
#endif
  bp->spin=0.03f;
  return J_BH;
}

/* ------------test spin, NLT ----------------*/

/* spin direction */
/*const float spin_vec[3] = {bp->spin * bp->angular_momentum_direction[0],
			   bp->spin * bp->angular_momentum_direction[1],
			   bp->spin * bp->angular_momentum_direction[2]}; */



/*-----------------end test---------------------*/


/* jet kicks skipped */
#endif /* SWIFT_SPIN_JET_BLACK_HOLES_SPIN_H */
