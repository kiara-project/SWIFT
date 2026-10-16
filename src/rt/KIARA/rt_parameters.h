/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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

#ifndef SWIFT_KIARA_RT_PARAMETERS_H
#define SWIFT_KIARA_RT_PARAMETERS_H

/**
 * @file src/rt/KIARA/rt_parameters.h
 * @brief Global RT parameters.
 */

extern struct rt_parameters rt_params;

/**
 * @brief Some global RT related parameters.
 */
struct rt_parameters {

  /*! Reduced speed of light (internal units) */
  float reduced_speed_of_light;

  /*! Inverse of the reduced speed of light (internal units) */
  float reduced_speed_of_light_inverse;

  /*! Variable reduced speed of light: a particle with comoving smoothing
   * length h gets c_red * (h / varc_h_ref)^varc_alpha, limited to the range
   * [reduced_speed_of_light, reduced_speed_of_light_max]. Disabled (uniform
   * c_red) if varc_h_ref_inv == 0. */
  float reduced_speed_of_light_max;

  /*! Inverse of the comoving smoothing length above which c increases */
  float varc_h_ref_inv;

  /*! Power of h / varc_h_ref (1: same RT time step for all particles) */
  float varc_alpha;
};

#endif /* SWIFT_KIARA_RT_PARAMETERS_H */
