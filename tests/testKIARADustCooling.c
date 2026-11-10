/* The size-resolved dust steps of KIARA's cooling: the bins kept consistent
 * with the element dust, the SN shock shrink, and the grain collisions
 * (Li et al. 2021, section 2.5.3; Hirashita & Yan 2009), against
 * independent evaluations of the published formulae. */
#include "chemistry/KIARA/dust_coagulation.h"
#include "chemistry/KIARA/dust_evolution.h"
#include "chemistry/KIARA/dust_shattering.h"
#include "chemistry/KIARA/dust_sources.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#define N KIARA_DUST_N_BINS
#define approx(x, y, tol) (fabs((x) - (y)) <= (tol) * fmax(fabs(y), 1.e-300))

static double row_sum(const double *x) {
  double total = 0.;
  for (int i = 0; i < N; ++i) total += x[i];
  return total;
}

/* Mass-weighted mean log radius of a row */
static double mean_log_radius(const double *radius, const double *x) {
  double m = 0., ml = 0.;
  for (int i = 0; i < N; ++i) {
    m += x[i];
    ml += x[i] * log(radius[i]);
  }
  return ml / m;
}

/* The bins follow the materials' dust masses from the element composition,
 * whatever the (possibly stale) rows hold; their shapes are kept, and empty
 * rows of a material with dust get the seed */
static void check_composition(void) {
  float dist[KIARA_DUST_N_MATERIALS][N];
  double seed[N], mass[KIARA_DUST_N_MATERIALS][N];
  for (int i = 0; i < N; ++i) {
    dist[kiara_dust_carbonaceous][i] = .1f * (i + 1);
    dist[kiara_dust_silicate][i] = (i == N - 1) ? .5f : 0.f;
    seed[i] = 1. / N;
  }
  const double dust = 2., f_C = .3;
  kiara_dust_bin_masses((const float(*)[N])dist, dust, f_C, seed, mass);
  assert(approx(row_sum(mass[kiara_dust_carbonaceous]), dust * f_C, 1.e-6));
  assert(approx(row_sum(mass[kiara_dust_silicate]), dust * (1. - f_C), 1.e-6));
  for (int i = 1; i < N; ++i)
    assert(approx(
        mass[kiara_dust_carbonaceous][i] / mass[kiara_dust_carbonaceous][0],
        i + 1., 1.e-6));
  assert(mass[kiara_dust_silicate][N - 1] == dust * (1. - f_C));

  /* An empty row takes the seed, or stays empty without one */
  for (int i = 0; i < N; ++i) dist[kiara_dust_silicate][i] = 0.f;
  kiara_dust_bin_masses((const float(*)[N])dist, dust, f_C, seed, mass);
  for (int i = 0; i < N; ++i)
    assert(approx(mass[kiara_dust_silicate][i], dust * (1. - f_C) / N, 1.e-12));
  kiara_dust_bin_masses((const float(*)[N])dist, dust, f_C, NULL, mass);
  assert(row_sum(mass[kiara_dust_silicate]) == 0.);

  /* Back to fractions of the dust mass: rows follow the material masses,
   * negative masses are dropped, the shapes are kept */
  double m_material[KIARA_DUST_N_MATERIALS] = {.6, 1.4};
  for (int i = 0; i < N; ++i) {
    mass[kiara_dust_carbonaceous][i] = 3. * (i + 1);
    mass[kiara_dust_silicate][i] = (i == 0) ? -1. : 1.;
  }
  kiara_dust_set_distribution(dist, (const double(*)[N])mass, m_material, seed);
  double f[KIARA_DUST_N_MATERIALS] = {0., 0.};
  for (int m = 0; m < KIARA_DUST_N_MATERIALS; ++m)
    for (int i = 0; i < N; ++i) {
      assert(dist[m][i] >= 0.f);
      f[m] += dist[m][i];
    }
  assert(approx(f[kiara_dust_carbonaceous], .3, 1.e-6));
  assert(approx(f[kiara_dust_silicate], .7, 1.e-6));
  assert(dist[kiara_dust_silicate][0] == 0.f);
  assert(approx(
      dist[kiara_dust_carbonaceous][N - 1] / dist[kiara_dust_carbonaceous][0],
      (double)N, 1.e-6));

  /* No dust, no distribution */
  m_material[0] = m_material[1] = 0.;
  kiara_dust_set_distribution(dist, (const double(*)[N])mass, m_material, seed);
  for (int m = 0; m < KIARA_DUST_N_MATERIALS; ++m)
    for (int i = 0; i < N; ++i) assert(dist[m][i] == 0.f);
}

/* SN shocks shrink every grain by the same factor cbrt(1 - f), which removes
 * the fraction f of the mass of any distribution */
static void check_sn_shrink(const double *radius) {
  double in[N], out[N];
  for (int i = 0; i < N; ++i) in[i] = 1. + (i % 3);
  for (double f = .01; f < 1.; f += .2) {
    const double left =
        kiara_dust_radius_remap(radius, in, 0., cbrt(1. - f), 1, out);
    assert(approx(left, (1. - f) * row_sum(in), 1.e-12));
    assert(approx(row_sum(out), left, 1.e-12));
  }
}

/* Li et al. (2021), Eq. 15, and Hirashita & Yan (2009), Eq. 8, evaluated
 * here from their published forms */
static void check_velocities(void) {
  const double pi = 3.14159265358979323846, G = 6.67430e-8, m_H = 1.6735575e-24,
               pc = 3.0856775814913673e18;
  const double cases[3][4] = {/* a (cm), n (cm^-3), c_s (cm/s), rho_gr */
                              {1.e-5, 1., 1.e5, 2.4},
                              {3.e-6, 100., 3.e4, 3.3},
                              {1.e-4, .1, 1.5e6, 2.2}};
  for (int k = 0; k < 3; ++k) {
    const double a = cases[k][0], n = cases[k][1], c = cases[k][2],
                 rho_gr = cases[k][3];
    const double rho = n * m_H;
    const double L_J = .5 * sqrt(pi * c * c / (G * rho));
    const double v_J = .7e5 * sqrt(L_J / pc);
    const double sigma = .06e5 * (v_J / c) * (v_J / c) * (a / 1.e-5) *
                         pow(n, -.5) * (rho_gr / 2.4);
    assert(approx(kiara_dust_collision_velocity_dispersion(a, rho, c, rho_gr),
                  sigma, 1.e-12));
  }

  /* Sticking threshold: graphite and quartz of Hirashita & Yan's table 2,
   * with their reduced modulus 1/E = (1 - nu)^2/E_1 + (1 - nu)^2/E_2 */
  const double gamma[2] = {12., 25.}, E[2] = {3.4e10, 5.4e11},
               nu[2] = {.5, .17};
  for (int m = 0; m < KIARA_DUST_N_MATERIALS; ++m) {
    const double a1 = 2.e-7, a2 = 3.e-6, rho_gr = 2.4;
    const double E_red = 1. / (2. * (1. - nu[m]) * (1. - nu[m]) / E[m]);
    const double R = a1 * a2 / (a1 + a2);
    const double v =
        2.14 * 10. * sqrt((a1 * a1 * a1 + a2 * a2 * a2) / pow(a1 + a2, 3.)) *
        pow(gamma[m], 5. / 6.) / (cbrt(E_red) * pow(R, 5. / 6.) * sqrt(rho_gr));
    assert(approx(kiara_dust_coagulation_threshold(a1, a2, rho_gr,
                                                   (enum kiara_dust_material)m),
                  v, 1.e-12));
    /* Smaller grains stick at higher speeds */
    assert(kiara_dust_coagulation_threshold(a1, a1, rho_gr,
                                            (enum kiara_dust_material)m) >
           kiara_dust_coagulation_threshold(a2, a2, rho_gr,
                                            (enum kiara_dust_material)m));
  }
}

/* Each pair's products hold the mass of the colliding grains */
static void check_products(const struct kiara_dust_collision_grid *grid,
                           const struct kiara_dust_collision_operator *op) {
  for (int i = 0; i < N; ++i)
    for (int j = i; j < N; ++j) {
      if (op->kernel[i][j] == 0.) continue;
      assert(op->kernel[i][j] > 0. && i != j);
      double products = 0.;
      for (int k = 0; k < N; ++k) {
        assert(op->product[i][j][k] >= 0.);
        products += op->product[i][j][k];
      }
      assert(approx(products, grid->mass[i] + grid->mass[j], 1.e-12));
    }
}

/* Coagulation moves mass to larger grains, shattering to smaller ones; both
 * conserve the material's mass and keep the bins non-negative */
static void check_collisions(const double *radius) {
  struct kiara_dust_collision_grid grid;
  assert(kiara_dust_collision_grid_init(&grid, radius, 2.4));
  static struct kiara_dust_collision_operator op;
  double velocity[N][3];
  const double rho_dust = 1.e-25; /* g/cm^3 */

  for (int m = 0; m < KIARA_DUST_N_MATERIALS; ++m) {
    const enum kiara_dust_material material = (enum kiara_dust_material)m;
    for (int process = 0; process < 2; ++process) {
      /* Slow grains coagulate (a few cm/s, below every threshold); fast ones
       * shatter (tens of km/s, above 2 km/s for all pairs) */
      for (int i = 0; i < N; ++i)
        for (int k = 0; k < 3; ++k)
          velocity[i][k] =
              (k == 0) * (process ? 1.e6 * (1. + 4. * i) : .1 * (1. + i));
      if (process)
        kiara_dust_shattering_operator_init(&op, &grid, 2.4, material,
                                            (const double(*)[3])velocity);
      else
        kiara_dust_coagulation_operator_init(&op, &grid, 2.4, material,
                                             (const double(*)[3])velocity);
      check_products(&grid, &op);

      /* Equal-mass bins: lots of small grains to coagulate or large ones to
       * shatter; the row is a fraction of all the dust (.4 of it here) */
      double row[N], before[N];
      for (int i = 0; i < N; ++i) row[i] = before[i] = .4 / N;
      double rate = 0.;
      for (int i = 0; i < N; ++i)
        for (int j = i + 1; j < N; ++j) rate = fmax(rate, op.kernel[i][j]);
      assert(rate > 0.);

      /* A few depletion times of the fastest-depleted bin */
      double slope[N];
      const double loss =
          kiara_dust_collision_derivative(&grid, &op, rho_dust, row, slope);
      assert(isfinite(loss) && loss > 0.);
      assert(kiara_dust_collision_evolve(&grid, &op, rho_dust, 3. / loss, row));
      assert(fabs(mean_log_radius(radius, row) -
                  mean_log_radius(radius, before)) > .01);
      for (int i = 0; i < N; ++i) assert(row[i] >= 0.);
      /* Conserved up to the round-off of m_i + m_j - m_i - m_j for grain
       * masses 1e9 apart (the bins are floats, and rescaled to the element
       * composition before each crackle call) */
      assert(approx(row_sum(row), .4, 1.e-7));
      if (process)
        assert(mean_log_radius(radius, row) < mean_log_radius(radius, before));
      else
        assert(mean_log_radius(radius, row) > mean_log_radius(radius, before));
    }

    /* No collisions between grains moving together, or without dust */
    for (int i = 0; i < N; ++i)
      for (int k = 0; k < 3; ++k) velocity[i][k] = 3.e5;
    kiara_dust_shattering_operator_init(&op, &grid, 2.4, material,
                                        (const double(*)[3])velocity);
    double row[N];
    for (int i = 0; i < N; ++i) row[i] = 1. / N;
    assert(kiara_dust_collision_evolve(&grid, &op, rho_dust, 1.e16, row));
    for (int i = 0; i < N; ++i) assert(row[i] == 1. / N);
    for (int i = 0; i < N; ++i) velocity[i][0] = 1.e6 * (1. + i);
    kiara_dust_shattering_operator_init(&op, &grid, 2.4, material,
                                        (const double(*)[3])velocity);
    assert(kiara_dust_collision_evolve(&grid, &op, 0., 1.e16, row));
    for (int i = 0; i < N; ++i) assert(row[i] == 1. / N);
  }

  /* The rate of a single pair: grains of bins i < j stick at the rate
   * rho_d K (f_i / m_i) (f_j / m_j), K = pi (a_i + a_j)^2 |v_i - v_j| */
  const int i = 0, j = N - 1;
  memset(velocity, 0, sizeof(velocity));
  velocity[j][1] = .5;
  kiara_dust_coagulation_operator_init(&op, &grid, 2.4, kiara_dust_silicate,
                                       (const double(*)[3])velocity);
  double f[N] = {0}, dfdt[N];
  f[i] = .2;
  f[j] = .5;
  kiara_dust_collision_derivative(&grid, &op, rho_dust, f, dfdt);
  const double K = 3.14159265358979323846 * (grid.radius[i] + grid.radius[j]) *
                   (grid.radius[i] + grid.radius[j]) * .5;
  const double events =
      rho_dust * K * (f[i] / grid.mass[i]) * (f[j] / grid.mass[j]);
  assert(approx(dfdt[i], -events * grid.mass[i], 1.e-12));
  assert(fabs(row_sum(dfdt)) <= 1.e-12 * events * grid.mass[j]);
}

int main(void) {
  double radius[N];
  check_composition();
  check_velocities();
  for (int irregular = 0; irregular < 2; ++irregular) {
    assert(kiara_dust_grid_logarithmic(1.e-3, 1., radius));
    if (irregular)
      for (int i = 0; i < N; ++i) {
        const double x = (double)i / (N - 1);
        radius[i] = exp(log(1.e-3) * (1. - x * x));
      }
    check_sn_shrink(radius);
    check_collisions(radius);
  }
  printf("PASS dust composition, SN shrink and grain collisions: %d bins\n", N);
  return 0;
}
