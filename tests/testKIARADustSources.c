/* Physical source distributions and endpoint mass retention. */
#include "chemistry/KIARA/dust_sources.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>

int main(void) {
  double radius[KIARA_DUST_N_BINS], mass[KIARA_DUST_N_BINS];
  assert(kiara_dust_grid_logarithmic(1.e-4, 1., radius));
  for (int source = 0; source < 2; ++source) {
    const double sigma = source ? .47 : .6, p = source ? 4. : 0.;
    assert(kiara_dust_source_lognormal(radius, .1, sigma, p, mass));
    double sum = 0., cumulative = 0.;
    for (int i = 0; i < KIARA_DUST_N_BINS; ++i) {
      assert(isfinite(mass[i]) && mass[i] >= 0.);
      sum += mass[i];
      cumulative += mass[i];
      if (i < KIARA_DUST_N_BINS - 1) {
        const double boundary = .5 * (log(radius[i]) + log(radius[i + 1]));
        const double expected =
            .5 * erfc(-(boundary - log(.1) - (3. - p) * sigma * sigma) /
                      (sqrt(2.) * sigma));
        assert(fabs(cumulative - expected) < 1.e-13);
      }
    }
    assert(fabs(sum - 1.) < 1.e-13);
  }
  assert(!kiara_dust_source_lognormal(radius, .1, 0., 0., mass));
  assert(!kiara_dust_source_lognormal(radius, NAN, .6, 0., mass));
  assert(kiara_dust_source_lognormal(radius, 1.e-20, .01, 3., mass));
  assert(mass[0] == 1.);
  assert(kiara_dust_source_lognormal(radius, 1.e20, .01, 3., mass));
  assert(mass[KIARA_DUST_N_BINS - 1] == 1.);
  float distribution[KIARA_DUST_N_MATERIALS][KIARA_DUST_N_BINS];
  assert(kiara_dust_source_lognormal(radius, .1, .6, 0., mass));
  kiara_dust_distribution_seed(mass, .3f, distribution);
  double carbon = 0., silicate = 0.;
  for (int i = 0; i < KIARA_DUST_N_BINS; ++i) {
    carbon += distribution[0][i];
    silicate += distribution[1][i];
    assert(fabs(distribution[0][i] - .3 * mass[i]) < 1.e-7);
    assert(fabs(distribution[1][i] - .7 * mass[i]) < 1.e-7);
  }
  assert(fabs(carbon - .3) < 1.e-7 && fabs(silicate - .7) < 1.e-7);
  float empty[KIARA_DUST_N_BINS] = {0};
  kiara_dust_distribution_rescale_material(empty, .4f, mass);
  for (int i = 0; i < KIARA_DUST_N_BINS; ++i)
    assert(fabs(empty[i] - .4 * mass[i]) < 1.e-7);
  assert(kiara_dust_source_lognormal(radius, 1.e-20, .01, 3., mass));
  kiara_dust_distribution_seed(mass, .3f, distribution);
  assert(fabs(distribution[0][0] - .3) < 1.e-7);
  assert(kiara_dust_source_lognormal(radius, 1.e20, .01, 3., mass));
  kiara_dust_distribution_seed(mass, .3f, distribution);
  assert(fabs(distribution[1][KIARA_DUST_N_BINS - 1] - .7) < 1.e-7);
  printf("PASS dust sources: %d bins\n", KIARA_DUST_N_BINS);
  return 0;
}
