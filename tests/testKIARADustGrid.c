/* Standalone grid invariants, independent of SWIFT's configured bin count. */
#include "chemistry/KIARA/dust_grid.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>

static void check_deposit(const double *radius, const double a) {
  double bins[KIARA_DUST_N_BINS] = {0}, mass = 0., number = 0.;
  kiara_dust_grid_deposit(radius, a, 2., bins);
  for (int i = 0; i < KIARA_DUST_N_BINS; ++i) {
    assert(bins[i] >= 0.);
    mass += bins[i];
    number += bins[i] / pow(radius[i], 3.);
  }
  assert(fabs(mass - 2.) < 1.e-13);
  if (a < radius[0])
    assert(bins[0] == 2.);
  else if (a > radius[KIARA_DUST_N_BINS - 1])
    assert(bins[KIARA_DUST_N_BINS - 1] == 2.);
  else
    assert(fabs(number / (2. / pow(a, 3.)) - 1.) < 1.e-11);
}

int main(void) {
  double radius[KIARA_DUST_N_BINS];
  assert(!kiara_dust_grid_logarithmic(NAN, 1., radius));
  assert(!kiara_dust_grid_logarithmic(0., 1., radius));
  assert(!kiara_dust_grid_logarithmic(1., 1., radius));
  assert(!kiara_dust_grid_logarithmic(1., INFINITY, radius));
  for (int irregular = 0; irregular < 2; ++irregular) {
    assert(kiara_dust_grid_logarithmic(1.e-4, 1., radius));
    if (irregular)
      for (int i = 0; i < KIARA_DUST_N_BINS; ++i) {
        const double x = (double)i / (KIARA_DUST_N_BINS - 1);
        radius[i] = exp(log(1.e-4) * (1. - x * x));
      }
    assert(kiara_dust_grid_valid(radius));
    check_deposit(radius, radius[0] / 100.);
    check_deposit(radius, radius[KIARA_DUST_N_BINS - 1] * 100.);
    for (int i = 0; i < KIARA_DUST_N_BINS; ++i)
      check_deposit(radius, radius[i]);
    for (int i = 1; i < KIARA_DUST_N_BINS; ++i)
      check_deposit(radius, sqrt(radius[i - 1] * radius[i]));
  }
  radius[1] = NAN;
  assert(!kiara_dust_grid_valid(radius));
  printf("PASS dust grid: %d bins\n", KIARA_DUST_N_BINS);
  return 0;
}
