/*******************************************************************************
 * Li et al. (2021) grain shattering on KIARA's existing size grid.
 ******************************************************************************/
#ifndef SWIFT_KIARA_DUST_SHATTERING_H
#define SWIFT_KIARA_DUST_SHATTERING_H

#include "dust_collisions.h"

/** Deposit a fragment spectrum dn/da proportional to a^-3.3.
 * Integrate MASS analytically (dM/da proportional to a^-0.3) between geometric
 * bin boundaries. The end bins include unresolved fragments outside the grid:
 * no fragment mass is discarded or returned to gas. Grain number below the
 * minimum represented radius is unresolved. This is the mass-grid approximation
 * to Hirashita & Yan (2009), equations 12--14, adopted by Li (2021).
 */
static inline void kiara_dust_shattering_deposit_fragments(
    const struct kiara_dust_collision_grid *grid, const double radius_min,
    const double radius_max, const double fragment_mass,
    double product[KIARA_DUST_N_BINS]) {

  assert(radius_min > 0. && radius_max >= radius_min && fragment_mass >= 0.);
  if (fragment_mass == 0.) return;
  if (radius_max == radius_min) {
    const double ratio = radius_max / grid->radius[0];
    kiara_dust_collision_deposit(grid, grid->mass[0] * ratio * ratio * ratio,
                                 fragment_mass, product);
    return;
  }
  double weight[KIARA_DUST_N_BINS] = {0}, sum = 0.;
  int largest = 0;
  for (int k = 0; k < KIARA_DUST_N_BINS; ++k) {
    const double lo =
        k == 0 ? radius_min
               : fmax(radius_min, sqrt(grid->radius[k - 1] * grid->radius[k]));
    const double hi =
        k == KIARA_DUST_N_BINS - 1
            ? radius_max
            : fmin(radius_max, sqrt(grid->radius[k] * grid->radius[k + 1]));
    if (hi <= lo) continue;
    /* Divide radii by a_max and use expm1 for narrow intervals. The spectrum's
     * common normalization cancels when assigning its mass to the bins. */
    weight[k] = pow(hi / radius_max, 0.7) * (-expm1(0.7 * log(lo / hi)));
    sum += weight[k];
    if (weight[k] > weight[largest]) largest = k;
  }
  assert(sum > 0. && isfinite(sum));
  double remaining = fragment_mass;
  for (int k = 0; k < KIARA_DUST_N_BINS; ++k) {
    if (k == largest) continue;
    const double mass = fragment_mass * weight[k] / sum;
    product[k] += mass;
    remaining -= mass;
  }
  /* Put the round-off residual in the largest contribution, where it cannot
   * make a very small fragment population negative. */
  product[largest] += remaining;
}

/** Shock parameter sigma(M): Hirashita & Yan (2009), equation 11. */
static inline double kiara_dust_shattering_sigma(const double mach,
                                                 const double shock_slope) {
  const double x = shock_slope + 1. / mach;
  return 0.30 * pow(x - 0.11, 0.13) / (x - 1.);
}

/**
 * @brief Fragment one pair; return the mass fragmented from the larger grain.
 *
 * i < j: i is the projectile, j the target. Both have the same material.
 * Speed must exceed Li's common threshold of 2 km/s. Material sound speeds,
 * shock slopes, cratering and vaporization pressures are from Hirashita & Yan
 * table 2; the configured grain density replaces that table's densities.
 *
 * Implement their section 2.3, equations 9--19, as prescribed by Li section
 * 2.6: disrupt the entire projectile; crater 40% of the shocked target mass,
 * or disrupt the entire target if more than half its mass is shocked. Deposit
 * the surviving target remnant separately. Vaporization pressure sets fragment
 * sizes; it does not remove mass in this collision model.
 */
static inline double kiara_dust_shattering_products(
    const struct kiara_dust_collision_grid *grid, const int i, const int j,
    const double speed, const double grain_density,
    const enum kiara_dust_material material,
    double product[KIARA_DUST_N_BINS]) {

  assert(i >= 0 && i < j && j < KIARA_DUST_N_BINS);
  assert(speed > 2.e5 && grain_density > 0.);
  assert(material == kiara_dust_carbonaceous ||
         material == kiara_dust_silicate);
  const int carbon = material == kiara_dust_carbonaceous;
  const double c0 = carbon ? 1.8e5 : 5.e5;
  const double shock_slope = carbon ? 1.9 : 1.2;
  const double pressure_critical = carbon ? 4.e10 : 3.e11;
  const double pressure_vaporization = carbon ? 5.8e12 : 5.4e12;
  const double phi = pressure_critical / (grain_density * c0 * c0);
  const double mach_critical =
      2. * phi / (1. + sqrt(1. + 4. * shock_slope * phi));
  const double mach = speed / c0;
  const double sigma_critical =
      kiara_dust_shattering_sigma(mach_critical, shock_slope);
  const double sigma_impact =
      kiara_dust_shattering_sigma(mach / 2., shock_slope);
  const double shocked_mass =
      grid->mass[i] * (3. / pow(2., 9. / 16.)) * pow(sigma_impact, -1. / 9.) *
      pow(mach * mach / (sigma_critical * mach_critical * mach_critical),
          8. / 9.);
  const int catastrophic = shocked_mass > 0.5 * grid->mass[j];
  const double target_fragments =
      catastrophic ? grid->mass[j] : 0.4 * shocked_mass;

  for (int parent = i; parent <= j; parent += j - i) {
    const int target = parent == j;
    const double mass = target ? target_fragments : grid->mass[i];
    double a_max, a_min;
    if (target && !catastrophic) {
      /* Equation 18 with z=3.4. */
      const double z = 3.4;
      const double volume_factor = (16. / 3.) * 3.14159265358979323846 *
                                   (z * z * z * (z - 2.) / (z + 1.));
      a_max = cbrt(mass / (grain_density * volume_factor));
      a_min = a_max * pow(pressure_critical / pressure_vaporization, 1.47);
    } else {
      const int other = target ? i : j;
      const double critical_speed =
          c0 * pow(grid->mass[other] / (3. * grid->mass[parent]), 9. / 16.) *
          sqrt(sigma_critical) * pow(sigma_impact, 1. / 16.) * 2. *
          mach_critical;
      /* On a coarse grid, do not assign a fragment a larger radius than its
       * parent when the critical-spalling estimate exceeds that radius. */
      a_max = fmin(grid->radius[parent],
                   0.22 * grid->radius[parent] * critical_speed / speed);
      a_min = 0.03 * a_max;
    }
    kiara_dust_shattering_deposit_fragments(grid, a_min, a_max, mass, product);
  }

  const double remnant = grid->mass[j] - target_fragments;
  if (remnant > 0.)
    kiara_dust_collision_deposit(grid, remnant, remnant, product);
  return target_fragments;
}

/** Build the shattering operator for one material and fixed bin velocities.
 * Li (2021) uses 2 km/s for both materials. Same-bin velocities coincide;
 * different-material collisions and changes to elemental depletion are absent.
 */
static inline void kiara_dust_shattering_operator_init(
    struct kiara_dust_collision_operator *op,
    const struct kiara_dust_collision_grid *grid, const double grain_density,
    const enum kiara_dust_material material,
    const double velocity[KIARA_DUST_N_BINS][3]) {

  const double pi = 3.14159265358979323846;
  for (int i = 0; i < KIARA_DUST_N_BINS; ++i) {
    for (int j = i; j < KIARA_DUST_N_BINS; ++j) {
      op->kernel[i][j] = 0.;
      for (int k = 0; k < KIARA_DUST_N_BINS; ++k) op->product[i][j][k] = 0.;
      if (i == j) continue;
      double speed_squared = 0.;
      for (int axis = 0; axis < 3; ++axis) {
        const double dv = velocity[i][axis] - velocity[j][axis];
        speed_squared += dv * dv;
      }
      const double speed = sqrt(speed_squared);
      /* Li section 2.6: common 2 km/s shattering threshold. */
      if (speed <= 2.e5) continue;
      const double radius_sum = grid->radius[i] + grid->radius[j];
      op->kernel[i][j] = pi * radius_sum * radius_sum * speed;
      kiara_dust_shattering_products(grid, i, j, speed, grain_density, material,
                                     op->product[i][j]);
    }
  }
}

#endif /* SWIFT_KIARA_DUST_SHATTERING_H */
