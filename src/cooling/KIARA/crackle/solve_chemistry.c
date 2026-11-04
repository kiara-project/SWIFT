/***********************************************************************
/
/ Solve the chemistry and cooling
/
/
/ Copyright (c) 2013, Enzo/Grackle Development Team.
/
/ Distributed under the terms of the Enzo Public Licence.
/
/ The full license is in the file LICENSE, distributed with this 
/ software.
************************************************************************/

#include <stdio.h>
#include <math.h>
#include "grackle_macros.h"
#include "grackle_types.h"
#include "grackle_chemistry_data.h"
#include "phys_constants.h"
#include "crackle_chemistry.c"

extern chemistry_data *grackle_data;
extern chemistry_data_storage grackle_rates;

/* function prototypes */

double get_temperature_units(code_units *my_units);

int update_UVbackground_rates(chemistry_data *my_chemistry,
                              chemistry_data_storage *my_rates,
                              photo_rate_storage *my_uvb_rates,
                              code_units *my_units);

/**
 * @brief The UV background rates for the current redshift (my_units->a_value),
 * or the constant photo rates if the UV background is off.
 */
int crackle_get_uvb_rates(chemistry_data *my_chemistry,
                          chemistry_data_storage *my_rates,
                          photo_rate_storage *my_uvb_rates,
                          code_units *my_units)
{

  my_uvb_rates->k24 = my_uvb_rates->k25 = my_uvb_rates->k26 =
    my_uvb_rates->k27 = my_uvb_rates->k28 = my_uvb_rates->k29 =
    my_uvb_rates->k30 = my_uvb_rates->k31 = my_uvb_rates->piHI =
    my_uvb_rates->piHeI = my_uvb_rates->piHeII = my_uvb_rates->crsHI =
    my_uvb_rates->crsHeI = my_uvb_rates->crsHeII =
    my_uvb_rates->comp_xray = my_uvb_rates->temp_xray = 0.;

  if (my_chemistry->UVbackground == 1) {
    if (update_UVbackground_rates(my_chemistry, my_rates,
                                  my_uvb_rates, my_units) == FAIL) {
      fprintf(stderr, "Error in update_UVbackground_rates.\n");
      return FAIL;
    }
  }
  else {
    my_uvb_rates->k24       = my_rates->k24;
    my_uvb_rates->k25       = my_rates->k25;
    my_uvb_rates->k26       = my_rates->k26;
    my_uvb_rates->k27       = my_rates->k27;
    my_uvb_rates->k28       = my_rates->k28;
    my_uvb_rates->k29       = my_rates->k29;
    my_uvb_rates->k30       = my_rates->k30;
    my_uvb_rates->k31       = my_rates->k31;
    my_uvb_rates->piHI      = my_rates->piHI;
    my_uvb_rates->piHeI     = my_rates->piHeI;
    my_uvb_rates->piHeII    = my_rates->piHeII;
    my_uvb_rates->crsHI     = my_rates->crsHI;
    my_uvb_rates->crsHeI    = my_rates->crsHeI;
    my_uvb_rates->crsHeII   = my_rates->crsHeII;
    my_uvb_rates->comp_xray = my_rates->comp_xray;
    my_uvb_rates->temp_xray = my_rates->temp_xray;
  }

  return SUCCESS;
}

/**
 * @brief Everything the per-particle crackle calls need that depends only on
 * the parameters, units and redshift: call once per step, then use the
 * crackle_*_part functions.
 *
 * Sets my_units->temperature_units, the UV background rates and the crackle
 * units, exactly as solve_chemistry() does for every call.
 */
int crackle_prepare_step(chemistry_data *my_chemistry,
                         const chemistry_data_storage *my_rates,
                         code_units *my_units,
                         photo_rate_storage *my_uvb_rates,
                         crackle_units *cunits)
{
  if (crackle_get_uvb_rates(my_chemistry, (chemistry_data_storage *)my_rates,
                            my_uvb_rates, my_units) == FAIL)
    return FAIL;
  my_units->temperature_units = get_temperature_units(my_units);
  set_crackle_units(my_units, my_rates, my_chemistry->Gamma, cunits);
  return SUCCESS;
}

int local_solve_chemistry(chemistry_data *my_chemistry,
                          chemistry_data_storage *my_rates,
                          code_units *my_units,
                          grackle_field_data *my_fields,
                          double dt_value)
{

  /* Return if this doesn't concern us. */

  if (!my_chemistry->use_grackle)
    return SUCCESS;

  /* Update UV background rates. */
  photo_rate_storage my_uvb_rates;
  if (crackle_get_uvb_rates(my_chemistry, my_rates, &my_uvb_rates, my_units) == FAIL)
    return FAIL;

  /* Check for a metal field. */

  int metal_field_present = TRUE;
  if (my_fields->metal_density == NULL)
    metal_field_present = FALSE;

  double co_length_units, co_density_units;
  if (my_units->comoving_coordinates == TRUE) {
    co_length_units = my_units->length_units;
    co_density_units = my_units->density_units;
  }
  else {
    co_length_units = my_units->length_units *
      my_units->a_value * my_units->a_units;
    co_density_units = my_units->density_units /
      POW(my_units->a_value * my_units->a_units, 3);
  }

  /* Error checking for H2 shielding approximation */
  if (my_chemistry->H2_self_shielding == 1 && my_fields->grid_rank != 3){
    fprintf(stderr, "Error in solve_chemistry: H2 self-shielding option 1 "
                    "will only work for 3D Cartesian grids. Use option 2 "
                    "to provide an array of shielding lengths with "
                    "H2_self_shielding_length or option 3 to use the "
                    "local Jeans length.");
    return FAIL;
  }

  /* Calculate temperature units. */

  double temperature_units = get_temperature_units(my_units);
  my_units->temperature_units = temperature_units;

  int ierr = 0;

  /* Only crackle (use_grackle = 2) is available; Grackle's Fortran solver
   * is not part of KIARA */
  if (my_chemistry->use_grackle == 2 || dt_value < 0.) {
    ierr = crackle_solve_chemistry(my_fields, my_chemistry, *my_rates, my_uvb_rates, my_units, fabs(dt_value));
  }
  else {
    fprintf(stderr, "use_grackle=%d not available: only crackle (use_grackle=2)\n",my_chemistry->use_grackle);
    ierr = FAIL;
  }

  if (ierr == FAIL) {
    fprintf(stderr, "Error in crackle_solve_chemistry.\n");
  }

  return ierr;

}

int _solve_chemistry(chemistry_data *my_chemistry,
                     chemistry_data_storage *my_rates,
                     code_units *my_units, double dt_value, double dx_value,
                     int grid_rank, int *grid_dimension,
                     int *grid_start, int *grid_end,
                     gr_float *density, gr_float *internal_energy,
                     gr_float *x_velocity, gr_float *y_velocity, gr_float *z_velocity,
                     gr_float *HI_density, gr_float *HII_density, gr_float *HM_density,
                     gr_float *HeI_density, gr_float *HeII_density, gr_float *HeIII_density,
                     gr_float *H2I_density, gr_float *H2II_density,
                     gr_float *DI_density, gr_float *DII_density, gr_float *HDI_density,
                     gr_float *e_density, gr_float *metal_density,
                     gr_float *volumetric_heating_rate, gr_float *specific_heating_rate,
                     gr_float *RT_heating_rate, gr_float *RT_HI_ionization_rate, gr_float *RT_HeI_ionization_rate,
                     gr_float *RT_HeII_ionization_rate, gr_float *RT_H2_dissociation_rate,
                     gr_float *H2_self_shielding_length)
{

  grackle_field_data my_fields;
  my_fields.grid_dx                  = dx_value;
  my_fields.grid_rank                = grid_rank;
  my_fields.grid_dimension           = grid_dimension;
  my_fields.grid_start               = grid_start;
  my_fields.grid_end                 = grid_end;
  my_fields.density                  = density;
  my_fields.internal_energy          = internal_energy;
  my_fields.x_velocity               = x_velocity;
  my_fields.y_velocity               = y_velocity;
  my_fields.z_velocity               = z_velocity;
  my_fields.HI_density               = HI_density;
  my_fields.HII_density              = HII_density;
  my_fields.HM_density               = HM_density;
  my_fields.HeI_density              = HeI_density;
  my_fields.HeII_density             = HeII_density;
  my_fields.HeIII_density            = HeIII_density;
  my_fields.H2I_density              = H2I_density;
  my_fields.H2II_density             = H2II_density;
  my_fields.DI_density               = DI_density;
  my_fields.DII_density              = DII_density;
  my_fields.HDI_density              = HDI_density;
  my_fields.e_density                = e_density;
  my_fields.metal_density            = metal_density;
  my_fields.volumetric_heating_rate  = volumetric_heating_rate;
  my_fields.specific_heating_rate    = specific_heating_rate;
  my_fields.RT_heating_rate          = RT_heating_rate;
  my_fields.RT_HI_ionization_rate    = RT_HI_ionization_rate;
  my_fields.RT_HeI_ionization_rate   = RT_HeI_ionization_rate;
  my_fields.RT_HeII_ionization_rate  = RT_HeII_ionization_rate;
  my_fields.RT_H2_dissociation_rate  = RT_H2_dissociation_rate;
  my_fields.H2_self_shielding_length = H2_self_shielding_length;

  if (local_solve_chemistry(my_chemistry, my_rates,
                            my_units, &my_fields, dt_value) == FAIL) {
    fprintf(stderr, "Error in local_solve_chemistry.\n");
    return FAIL;
  }
  return SUCCESS;
}

int solve_chemistry(code_units *my_units,
                    grackle_field_data *my_fields,
                    double dt_value)
{
  if (local_solve_chemistry(grackle_data, &grackle_rates,
                            my_units, my_fields, dt_value) == FAIL) {
    fprintf(stderr, "Error in local_solve_chemistry.\n");
    return FAIL;
  }
  return SUCCESS;
}
