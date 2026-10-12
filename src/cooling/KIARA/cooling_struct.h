/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_COOLING_STRUCT_KIARA_H
#define SWIFT_COOLING_STRUCT_KIARA_H

/**
 * @brief Properties of the cooling stored in the #part data.
 */
struct cooling_part_data {

  /*! Subgrid temperature */
  float subgrid_temp;

  /*! Subgrid density (internal units, physical frame) */
  float subgrid_dens;

  /*! Subgrid fraction of cold mass */
  float subgrid_fcold;

#if COOLING_GRACKLE_MODE >= 2
  /* Dust properties when use_grackle_dust_evol=1 */

  /*! Total mass in dust */
  float dust_mass;

  /*! Fraction of the dust mass in each element (index 0 holds the sum) */
  float dust_mass_fraction[chemistry_element_count];

  /*! Dust temperature computed by grackle (K) */
  float dust_temperature;
#endif

  /*! Cooling time in mixing layer between stream and ambient gas */
  float mixing_layer_cool_time;

  /*! The strength of the interstellar radiation field in Habing units */
  float G0;
};

/**
 * @brief Properties of the cooling stored in the extra particle data
 */
struct cooling_xpart_data {

  /*! Energy radiated away by this particle since the start of the run */
  float radiated_energy;

  /*! Last time the cooling was switch off */
  double time_last_event;

/* Here all fractions are mass fractions (the electron fraction is
 * n_e m_H / rho, the grackle convention) */
#if COOLING_GRACKLE_MODE >= 1
  /*! Neutral hydrogen mass fraction */
  float HI_frac;

  /*! Ionised hydrogen mass fraction */
  float HII_frac;

  /*! Neutral helium mass fraction */
  float HeI_frac;

  /*! Singly ionised helium mass fraction */
  float HeII_frac;

  /*! Doubly ionised helium mass fraction */
  float HeIII_frac;

  /*! Electron fraction */
  float e_frac;

#if COOLING_GRACKLE_MODE >= 2
  /*! H- mass fraction */
  float HM_frac;

  /*! Molecular hydrogen mass fraction */
  float H2I_frac;

  /*! Ionised molecular hydrogen (H2+) mass fraction */
  float H2II_frac;

#if COOLING_GRACKLE_MODE >= 3
  /*! Neutral deuterium mass fraction */
  float DI_frac;

  /*! Ionised deuterium mass fraction */
  float DII_frac;

  /*! HD mass fraction */
  float HDI_frac;
#endif  // MODE >= 3

#endif  // MODE >= 2

#endif  // MODE >= 1

  /*! Metal mass fraction (only printed in debugging output) */
  float metal_frac;
};

#endif /* SWIFT_COOLING_STRUCT_KIARA_H */
