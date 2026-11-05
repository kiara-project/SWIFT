/* Analytic radius evolution, mass loss and arbitrary-grid boundaries. */
#include "chemistry/KIARA/dust_coagulation.h"
#include "chemistry/KIARA/dust_evolution.h"
#include "chemistry/KIARA/dust_shattering.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>

static double sum(const double *x) {
  double total = 0.;
  for (int i = 0; i < KIARA_DUST_N_BINS; ++i) {
    assert(isfinite(x[i]) && x[i] >= 0.);
    total += x[i];
  }
  return total;
}

static void check_transport(const double *radius, int bin, double shift,
                            double scale, int evolve_mass) {
  double input[KIARA_DUST_N_BINS] = {0}, output[KIARA_DUST_N_BINS];
  input[bin] = .7;
  const double a = fmax(0., radius[bin] * scale + shift);
  const double expected = evolve_mass ? .7 * pow(a / radius[bin], 3.) : .7;
  const double total =
      kiara_dust_radius_remap(radius, input, shift, scale, evolve_mass, output);
  assert(fabs(total - expected) <= 1.e-12 * fmax(expected, 1.e-30));
  if (a <= radius[0]) assert(output[0] == total);
  if (a >= radius[KIARA_DUST_N_BINS - 1])
    assert(output[KIARA_DUST_N_BINS - 1] == total);
  if (a > radius[0] && a < radius[KIARA_DUST_N_BINS - 1]) {
    double grains = 0.;
    for (int i = 0; i < KIARA_DUST_N_BINS; ++i)
      grains += output[i] / pow(radius[i], 3.);
    assert(fabs(grains / (total / pow(a, 3.)) - 1.) < 1.e-10);
  }
}

int main(void) {
  double radius[KIARA_DUST_N_BINS];
  for (int irregular = 0; irregular < 2; ++irregular) {
    assert(kiara_dust_grid_logarithmic(1.e-4, 1., radius));
    if (irregular)
      for (int i = 0; i < KIARA_DUST_N_BINS; ++i) {
        const double x = (double)i / (KIARA_DUST_N_BINS - 1);
        radius[i] = exp(log(1.e-4) * (1. - x * x));
      }
    for (int i = 0; i < KIARA_DUST_N_BINS; ++i) {
      check_transport(radius, i, 0., 1., 1);
      check_transport(radius, i, .03, 1., 1);
      check_transport(radius, i, -.003, 1., 1);
      check_transport(radius, i, -2., 1., 1);
      check_transport(radius, i, 0., .9, 0);
      check_transport(radius, i, 0., 0., 0);
      check_transport(radius, i, 100., 1., 0);
    }
    struct kiara_dust_collision_grid grid;
    assert(kiara_dust_collision_grid_init(&grid, radius, 2.4));
    double product[KIARA_DUST_N_BINS] = {0};
    kiara_dust_collision_deposit(&grid, grid.mass[0] * .01, 2., product);
    assert(product[0] == 2.);
    kiara_dust_collision_deposit(&grid, grid.mass[KIARA_DUST_N_BINS - 1] * 10.,
                                 3., product);
    assert(fabs(sum(product) - 5.) < 1.e-13);
    for (int i = 0; i < KIARA_DUST_N_BINS; ++i) product[i] = 0.;
    kiara_dust_shattering_deposit_fragments(&grid, grid.radius[0] * .001,
                                            grid.radius[0] * .01, 1., product);
    assert(product[0] == 1.);
    for (int i = 0; i < KIARA_DUST_N_BINS; ++i) product[i] = 0.;
    kiara_dust_shattering_deposit_fragments(
        &grid, grid.radius[KIARA_DUST_N_BINS - 1] * 2.,
        grid.radius[KIARA_DUST_N_BINS - 1] * 4., 1., product);
    assert(product[KIARA_DUST_N_BINS - 1] == 1.);
    for (int i = 0; i < KIARA_DUST_N_BINS; ++i) product[i] = 0.;
    kiara_dust_shattering_deposit_fragments(
        &grid, grid.radius[0] * .001, grid.radius[KIARA_DUST_N_BINS - 1] * 4.,
        1., product);
    assert(fabs(sum(product) - 1.) < 1.e-13);
  }
  printf("PASS dust evolution and collision boundaries: %d bins\n",
         KIARA_DUST_N_BINS);
  return 0;
}
