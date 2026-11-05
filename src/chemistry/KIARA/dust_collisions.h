/*******************************************************************************
 * Shared numerical helpers for KIARA grain collisions.
 ******************************************************************************/
#ifndef SWIFT_KIARA_DUST_COLLISIONS_H
#define SWIFT_KIARA_DUST_COLLISIONS_H

#include "dust.h"

#include <assert.h>
#include <float.h>
#include <math.h>

/** Representative radii and single-grain masses, in cgs units.
 *
 * The existing KIARA radii are mass pivots, not bin edges. Collision products
 * are interpolated between these pivots to avoid suppressing coagulation on
 * a coarse radius grid. No extra particle fields or size bins are needed.
 */
struct kiara_dust_collision_grid {
  double radius[KIARA_DUST_N_BINS];
  double mass[KIARA_DUST_N_BINS];
};

/** Fixed collision kernels (cm^3/s) and products (grams per collision).
 * Only i <= j is stored/used. Products must sum to m_i + m_j, including
 * surviving remnants. Keeping products separate from populations ensures
 * newly deposited mass is not processed again in the same Euler stage.
 */
struct kiara_dust_collision_operator {
  double kernel[KIARA_DUST_N_BINS][KIARA_DUST_N_BINS];
  double product[KIARA_DUST_N_BINS][KIARA_DUST_N_BINS][KIARA_DUST_N_BINS];
};

/**
 * @brief Build a collision grid from radii in micron and density in g/cm^3.
 * @return One on success, zero for non-finite, non-positive or unordered data.
 * The output is only changed on success.
 */
static inline int kiara_dust_collision_grid_init(
    struct kiara_dust_collision_grid *grid,
    const double radius_micron[KIARA_DUST_N_BINS], const double grain_density) {

  if (!isfinite(grain_density) || grain_density <= 0.) return 0;

  struct kiara_dust_collision_grid candidate;
  const double four_pi_over_three = 4.1887902047863909846;
  for (int i = 0; i < KIARA_DUST_N_BINS; ++i) {
    candidate.radius[i] = radius_micron[i] * 1.e-4;
    const double a = candidate.radius[i];
    candidate.mass[i] = four_pi_over_three * grain_density * a * a * a;
    if (!isfinite(a) || a <= 0. || !isfinite(candidate.mass[i]) ||
        candidate.mass[i] <= 0.)
      return 0;
    if (i > 0 && (a <= candidate.radius[i - 1] ||
                  candidate.mass[i] <= candidate.mass[i - 1]))
      return 0;
  }
  *grid = candidate;
  return 1;
}

/**
 * @brief Deposit collision products into one material's mass array.
 *
 * @param grid Valid collision grid.
 * @param grain_mass Mass of one product grain, in grams (positive and finite).
 * @param total_mass Total mass of all such products (non-negative and finite).
 * May also be a mass fraction: the destination uses the same units.
 * @param bin_mass Destination material row, updated additively.
 *
 * Between pivots, interpolate product NUMBER linearly in grain mass, then
 * convert to bin masses. This conserves both deposited mass and grain number:
 * sum(delta_M_i) = total_mass and sum(delta_M_i / m_i) = total_mass/grain_mass.
 * Interpolating mass fractions directly would violate the number constraint.
 *
 * Below/above the represented range, retain all mass in the nearest endpoint.
 * This closed-boundary approximation conserves material mass, but cannot
 * conserve grain number outside the grid. It must not be interpreted as
 * vaporization or as resolved growth beyond the largest radius.
 *
 * This is a coarse-grid adaptation of Li et al. (2021), section 2.6;
 * their 41-bin, piecewise-constant discretization is not used here.
 */
static inline void kiara_dust_collision_deposit(
    const struct kiara_dust_collision_grid *grid, const double grain_mass,
    const double total_mass, double bin_mass[KIARA_DUST_N_BINS]) {

  assert(isfinite(grain_mass) && grain_mass > 0.);
  assert(isfinite(total_mass) && total_mass >= 0.);
  if (total_mass == 0.) return;

  if (grain_mass <= grid->mass[0]) {
    bin_mass[0] += total_mass;
    return;
  }
  const int last = KIARA_DUST_N_BINS - 1;
  if (grain_mass >= grid->mass[last]) {
    bin_mass[last] += total_mass;
    return;
  }

  for (int hi = 1; hi < KIARA_DUST_N_BINS; ++hi) {
    if (grain_mass <= grid->mass[hi]) {
      const double lo_mass = grid->mass[hi - 1];
      const double hi_mass = grid->mass[hi];
      const double lo_fraction = (lo_mass / grain_mass) *
                                 ((hi_mass - grain_mass) / (hi_mass - lo_mass));
      const double deposited_lo = total_mass * lo_fraction;
      bin_mass[hi - 1] += deposited_lo;
      bin_mass[hi] += total_mass - deposited_lo;
      return;
    }
  }
}

/**
 * @brief Li et al. (2021), equation 20, grain velocity dispersion in cm/s.
 *
 * All inputs are cgs: grain radius in cm, gas/grain densities in g/cm^3,
 * sound speed in cm/s. This returns the THREE-dimensional dispersion sigma_gr;
 * each Gaussian velocity component has standard deviation sigma_gr/sqrt(3).
 * The caller chooses the physical gas phase and performs code-unit conversion.
 *
 * L_J = sqrt(pi c_g^2 / (G rho_g))/2 and v_J = 0.7 sqrt(L_J/pc) km/s.
 * Reference: https://doi.org/10.1093/mnras/stab2196, section 2.6.
 */
static inline double kiara_dust_collision_velocity_dispersion(
    const double radius, const double gas_density, const double sound_speed,
    const double grain_density) {

  assert(isfinite(radius) && radius > 0.);
  assert(isfinite(gas_density) && gas_density > 0.);
  assert(isfinite(sound_speed) && sound_speed > 0.);
  assert(isfinite(grain_density) && grain_density > 0.);

  const double pi = 3.14159265358979323846;
  const double gravity_cgs = 6.67430e-8;
  const double hydrogen_mass_cgs = 1.6735575e-24;
  const double parsec_cgs = 3.0856775814913673e18;
  const double jeans_length =
      0.5 * sound_speed * sqrt(pi / (gravity_cgs * gas_density));
  const double mach_squared =
      (7.e4 / sound_speed) * (7.e4 / sound_speed) * (jeans_length / parsec_cgs);

  return 6.e3 * mach_squared * (radius / 1.e-5) *
         sqrt(hydrogen_mass_cgs / gas_density) * (grain_density / 2.4);
}

/** Compute one material's derivative and maximum net depletion frequency.
 * Fractions are relative to TOTAL dust mass, not normalized within a material;
 * dust_density is the total dust density in g/cm^3. No material conversion or
 * mass renormalization occurs. The factor 1/2 counts identical pairs once.
 */
static inline double kiara_dust_collision_derivative(
    const struct kiara_dust_collision_grid *grid,
    const struct kiara_dust_collision_operator *op, const double dust_density,
    const double fraction[KIARA_DUST_N_BINS],
    double derivative[KIARA_DUST_N_BINS]) {

  for (int i = 0; i < KIARA_DUST_N_BINS; ++i) derivative[i] = 0.;
  for (int i = 0; i < KIARA_DUST_N_BINS; ++i) {
    for (int j = i; j < KIARA_DUST_N_BINS; ++j) {
      const double rate = dust_density * op->kernel[i][j];
      if (rate == 0.) continue;
      const double events = rate * (fraction[i] / grid->mass[i]) *
                            (fraction[j] / grid->mass[j]) * (i == j ? 0.5 : 1.);
      for (int k = 0; k < KIARA_DUST_N_BINS; ++k) {
        /* Subtract parent masses before multiplying by the rate: unresolved
         * fragments returned to the same bin cause no net depletion. */
        const double net_mass = op->product[i][j][k] -
                                (k == i ? grid->mass[i] : 0.) -
                                (k == j ? grid->mass[j] : 0.);
        derivative[k] += events * net_mass;
      }
    }
  }
  double max_loss = 0., material_fraction = 0.;
  for (int i = 0; i < KIARA_DUST_N_BINS; ++i) material_fraction += fraction[i];
  for (int i = 0; i < KIARA_DUST_N_BINS; ++i) {
    if (!isfinite(derivative[i])) return INFINITY;
    if (derivative[i] < 0. && fraction[i] > 0.)
      max_loss =
          fmax(max_loss,
               -derivative[i] / fmax(fraction[i], 1.e-12 * material_fraction));
  }
  return max_loss;
}

/** Conservative Euler stage with a bounded collision flux.
 * Limit each pair's WHOLE event count to its available reactants, including
 * fragments and remnants in that same event count. This allows a trace bin
 * to be exhausted without clipping a negative population after the update.
 */
static inline int kiara_dust_collision_forward(
    const struct kiara_dust_collision_grid *grid,
    const struct kiara_dust_collision_operator *op, const double dust_density,
    const double step, const double fraction[KIARA_DUST_N_BINS],
    double forward[KIARA_DUST_N_BINS]) {

  double events[KIARA_DUST_N_BINS][KIARA_DUST_N_BINS];
  double sink[KIARA_DUST_N_BINS] = {0}, limiter[KIARA_DUST_N_BINS];
  for (int i = 0; i < KIARA_DUST_N_BINS; ++i) {
    for (int j = i; j < KIARA_DUST_N_BINS; ++j) {
      const double count =
          dust_density * op->kernel[i][j] * (fraction[i] / grid->mass[i]) *
          (fraction[j] / grid->mass[j]) * (i == j ? 0.5 : 1.) * step;
      if (!isfinite(count)) return 0;
      events[i][j] = count;
      for (int k = 0; k < KIARA_DUST_N_BINS; ++k) {
        const double net = op->product[i][j][k] -
                           (k == i ? grid->mass[i] : 0.) -
                           (k == j ? grid->mass[j] : 0.);
        if (net < 0.) sink[k] -= count * net;
      }
    }
  }
  for (int k = 0; k < KIARA_DUST_N_BINS; ++k) {
    if (!isfinite(sink[k])) return 0;
    limiter[k] = sink[k] > 0. && sink[k] >= fraction[k]
                     ? (1. - 4. * KIARA_DUST_N_BINS * DBL_EPSILON) *
                           fraction[k] / sink[k]
                     : 1.;
    forward[k] = fraction[k];
  }
  for (int i = 0; i < KIARA_DUST_N_BINS; ++i) {
    for (int j = i; j < KIARA_DUST_N_BINS; ++j) {
      double scale = 1.;
      for (int k = 0; k < KIARA_DUST_N_BINS; ++k) {
        const double net = op->product[i][j][k] -
                           (k == i ? grid->mass[i] : 0.) -
                           (k == j ? grid->mass[j] : 0.);
        if (net < 0.) scale = fmin(scale, limiter[k]);
      }
      for (int k = 0; k < KIARA_DUST_N_BINS; ++k) {
        const double net = op->product[i][j][k] -
                           (k == i ? grid->mass[i] : 0.) -
                           (k == j ? grid->mass[j] : 0.);
        forward[k] += events[i][j] * scale * net;
      }
    }
  }
  for (int k = 0; k < KIARA_DUST_N_BINS; ++k)
    if (!isfinite(forward[k]) || forward[k] < 0.) return 0;
  return 1;
}

/**
 * @brief Integrate a fixed collision operator with positive SSP-RK2 substeps.
 *
 * @param dust_density Total physical dust density in g/cm^3.
 * @param dt Duration in seconds. The caller holds velocities fixed over dt.
 * @param fraction One material's fractions of total dust mass, updated on
 * success.
 * @return One if the entire interval was integrated, zero on invalid input,
 * numerical failure or exceeding the substep budget. On failure the input is
 * unchanged, so the caller can report the failure or retry with shorter steps.
 *
 * Limit initial net depletion to 5% per substep, with an absolute scale of
 * 10^-12 of this material's mass for trace bins. Bound event counts when trace
 * reactants would be exhausted. Retry on numerical failure. Do not clip bins
 * or globally renormalize: both would hide conservation errors. A finite budget
 * prevents unbounded work on excessively stiff input; no time is skipped.
 */
static inline int kiara_dust_collision_evolve(
    const struct kiara_dust_collision_grid *grid,
    const struct kiara_dust_collision_operator *op, const double dust_density,
    const double dt, double fraction[KIARA_DUST_N_BINS]) {

  if (!isfinite(dt) || dt < 0. || !isfinite(dust_density) || dust_density < 0.)
    return 0;
  double current[KIARA_DUST_N_BINS];
  for (int i = 0; i < KIARA_DUST_N_BINS; ++i) {
    if (!isfinite(fraction[i]) || fraction[i] < 0.) return 0;
    current[i] = fraction[i];
  }
  if (dt == 0. || dust_density == 0.) return 1;

  double remaining = dt;
  int steps = 0;
  while (remaining > 0.) {
    if (++steps > 8192) return 0;
    double slope[KIARA_DUST_N_BINS];
    const double loss =
        kiara_dust_collision_derivative(grid, op, dust_density, current, slope);
    if (!isfinite(loss)) return 0;
    if (loss == 0.) break;
    double step = fmin(remaining, 0.05 / loss);
    double next[KIARA_DUST_N_BINS];
    int accepted = 0;
    for (int attempt = 0; attempt < 32 && !accepted; ++attempt) {
      if (!(step > 0.) || remaining - step == remaining) return 0;
      double stage[KIARA_DUST_N_BINS], forward[KIARA_DUST_N_BINS];
      if (kiara_dust_collision_forward(grid, op, dust_density, step, current,
                                       stage) &&
          kiara_dust_collision_forward(grid, op, dust_density, step, stage,
                                       forward)) {
        for (int i = 0; i < KIARA_DUST_N_BINS; ++i) {
          next[i] = 0.5 * current[i] + 0.5 * forward[i];
        }
        accepted = 1;
      }
      if (!accepted) step *= 0.5;
    }
    if (!accepted) return 0;
    for (int i = 0; i < KIARA_DUST_N_BINS; ++i) current[i] = next[i];
    remaining -= step;
  }
  for (int i = 0; i < KIARA_DUST_N_BINS; ++i) fraction[i] = current[i];
  return 1;
}

#endif /* SWIFT_KIARA_DUST_COLLISIONS_H */
