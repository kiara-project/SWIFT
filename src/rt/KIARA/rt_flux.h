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

#ifndef SWIFT_KIARA_RT_FLUX_H
#define SWIFT_KIARA_RT_FLUX_H

#if defined(RT_RIEMANN_SOLVER_GLF)
#include "rt_riemann_GLF.h"
#elif defined(RT_RIEMANN_SOLVER_HLL)
#error "KIARA RT transports W = (c E, F) (variable reduced speed of light); \
the HLL solver has not been converted, use GLF"
#include "rt_riemann_HLL.h"
#else
#error "No valid choice of RT Riemann solver has been selected"
#endif

#include "rt_unphysical.h"

#include <float.h>
#ifdef __AVX__
#include <immintrin.h>
#endif

/**
 * @file src/rt/KIARA/rt_flux.h
 * @brief Functions related to compute the interparticle flux term of the
 * conservation law. This is its own file so we can switch between Riemann
 * solvers more easily.
 */

/**
 * @brief Reset the fluxes for the given particle.
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void rt_part_reset_fluxes(
    struct part *restrict p) {

  for (int g = 0; g < RT_NGROUPS; g++) {
    p->rt_data.flux[g].energy = 0.f;
    p->rt_data.flux[g].flux[0] = 0.f;
    p->rt_data.flux[g].flux[1] = 0.f;
    p->rt_data.flux[g].flux[2] = 0.f;
  }
}

/**
 * @brief Compute the flux between a left state UL and a right
 * state UR along the direction of the unit vector n_unit
 * through a surface of size Anorm.
 *
 * @param UL left state (energy density, fluxes)
 * @param UR right state (energy density, fluxes)
 * @param n_unit unit vector of the direction of the surface
 * @param Anorm size of the surface through which the flux goes
 * @param fluxes (return) the resulting flux
 */
__attribute__((always_inline)) INLINE static void rt_compute_flux(
    float UL[4], float UR[4], const float n_unit[3], const float Anorm,
    float fluxes[4]) {

  /* Unphysical check not necessary here.
   * It's done in gradients_predict as well. */

  const float FLnorm = sqrtf(UL[1] * UL[1] + UL[2] * UL[2] + UL[3] * UL[3]);
  const float FRnorm = sqrtf(UR[1] * UR[1] + UR[2] * UR[2] + UR[3] * UR[3]);

  /* Get the fluxes in the hyperbolic conservation law sense. */
  float hyperFluxL[4][3];
  rt_get_hyperbolic_flux(UL, FLnorm, hyperFluxL);
  float hyperFluxR[4][3];
  rt_get_hyperbolic_flux(UR, FRnorm, hyperFluxR);

#ifdef SWIFT_RT_DEBUG_CHECKS
  rt_check_unphysical_hyperbolic_flux(hyperFluxL);
  rt_check_unphysical_hyperbolic_flux(hyperFluxR);
#endif

  rt_riemann_solve_for_flux(UL, UR, FLnorm, FRnorm, hyperFluxL, hyperFluxR,
                            n_unit, fluxes);

  /* get the actual flux */
  fluxes[0] *= Anorm;
  fluxes[1] *= Anorm;
  fluxes[2] *= Anorm;
  fluxes[3] *= Anorm;
}

/**
 * @brief Correct an extrapolated face state W = (c E, F) as
 * rt_check_unphysical_state() does (E <= 0: no radiation; |F| <= c E), and
 * return its flux projected on the face normal n, (F.n, P.n), with the M1
 * closure P = c E [(1 - chi)/2 I + (3 chi - 1)/2 F F / |F|^2]. Written with
 * the normalised flux f = F / (c E), so that it is free of under/overflows,
 * and (unless the flux has to be limited) with one reciprocal, one square
 * root and two divisions.
 *
 * @param U (in/out) The face state.
 * @param n_unit Unit vector normal to the face.
 * @param fluxn (return) The flux through the face (per unit area).
 */
__attribute__((always_inline)) INLINE static void rt_face_state_flux(
    float U[4], const float n_unit[3], float fluxn[4]) {

  if (isinf(U[0]) || isnan(U[0]))
    error("Got inf/nan radiation energy at a face | %.6e | %.6e %.6e %.6e",
          U[0], U[1], U[2], U[3]);

  if (U[0] <= 0.f) {
    U[0] = U[1] = U[2] = U[3] = 0.f;
    fluxn[0] = fluxn[1] = fluxn[2] = fluxn[3] = 0.f;
    return;
  }

  /* |F| > c E: limit the flux to c E. Done before normalising, as c E
   * may be tiny at a face (cancellation of the extrapolation) while F is not,
   * so that F / (c E) would overflow. */
  const float F2 = U[1] * U[1] + U[2] * U[2] + U[3] * U[3];
  if (F2 > U[0] * U[0]) {
    const float corr = U[0] / sqrtf(F2);
    U[1] *= corr;
    U[2] *= corr;
    U[3] *= corr;
  }

  const float U0_inv = 1.f / U[0];
  float fv[3] = {U[1] * U0_inv, U[2] * U0_inv, U[3] * U0_inv};
  float f2 = fv[0] * fv[0] + fv[1] * fv[1] + fv[2] * fv[2];
  if (f2 > 1.f) {
    /* Round-off (or F2 and (c E)^2 below the float range): normalise */
    const float corr = 1.f / sqrtf(f2);
    fv[0] *= corr;
    fv[1] *= corr;
    fv[2] *= corr;
    U[1] = fv[0] * U[0];
    U[2] = fv[1] * U[0];
    U[3] = fv[2] * U[0];
    f2 = 1.f;
  }

  const float fn = fv[0] * n_unit[0] + fv[1] * n_unit[1] + fv[2] * n_unit[2];
  fluxn[0] = fn * U[0];

  if (f2 == 0.f) {
    /* Isotropic: P = c E / 3 I */
    const float p = U[0] * (1.f / 3.f);
    fluxn[1] = p * n_unit[0];
    fluxn[2] = p * n_unit[1];
    fluxn[3] = p * n_unit[2];
    return;
  }

  const float chi = (3.f + 4.f * f2) / (5.f + 2.f * sqrtf(4.f - 3.f * f2));
  const float a_iso = 0.5f * (1.f - chi) * U[0];
  const float a_dir = 0.5f * (3.f * chi - 1.f) * U[0] * fn / f2;
  fluxn[1] = a_iso * n_unit[0] + a_dir * fv[0];
  fluxn[2] = a_iso * n_unit[1] + a_dir * fv[1];
  fluxn[3] = a_iso * n_unit[2] + a_dir * fv[2];
}

/*! Number of face states per interface: both sides for every group */
#define RT_NFACE (2 * RT_NGROUPS)

/**
 * @brief rt_face_state_flux() for all the face states of an interface at
 * once: the states of both particles in all groups are independent, so the
 * square roots and divisions of the M1 closure are done as one vector
 * operation (the loop has no branches, so the compiler vectorises it).
 * The isotropic case needs no branch: f = 0 gives chi = 1/3.
 *
 * @param W (in/out) Face states (c E, Fx, Fy, Fz) [component][state].
 * @param n_unit Unit vector normal to the face.
 * @param Fn (return) Fluxes through the face (per unit area).
 */
__attribute__((always_inline)) INLINE static void rt_face_state_fluxes(
    float W[4][RT_NFACE], const float n_unit[3], float Fn[4][RT_NFACE]) {

  for (int s = 0; s < RT_NFACE; s++)
    if (isinf(W[0][s]) || isnan(W[0][s]))
      error("Got inf/nan radiation energy at a face | %.6e | %.6e %.6e %.6e",
            W[0][s], W[1][s], W[2][s], W[3][s]);

#if defined(__AVX__) && (RT_NFACE == 8)
  /* The same as the loop below, written with AVX intrinsics: the compiler
   * does not vectorise the loop once it is inlined in the transport loops */
  const __m256 zero = _mm256_setzero_ps();
  const __m256 one = _mm256_set1_ps(1.f);
  const __m256 half = _mm256_set1_ps(0.5f);
  const __m256 three = _mm256_set1_ps(3.f);
  const __m256 four = _mm256_set1_ps(4.f);
  const __m256 tiny = _mm256_set1_ps(FLT_MIN);

  /* E <= 0: no radiation */
  const __m256 E_in = _mm256_loadu_ps(W[0]);
  const __m256 pos = _mm256_cmp_ps(E_in, zero, _CMP_GT_OQ);
  const __m256 E = _mm256_and_ps(E_in, pos);
  __m256 Fx = _mm256_and_ps(_mm256_loadu_ps(W[1]), pos);
  __m256 Fy = _mm256_and_ps(_mm256_loadu_ps(W[2]), pos);
  __m256 Fz = _mm256_and_ps(_mm256_loadu_ps(W[3]), pos);

  /* |F| > c E: limit the flux to c E, before normalising */
  const __m256 F2 = _mm256_add_ps(
      _mm256_add_ps(_mm256_mul_ps(Fx, Fx), _mm256_mul_ps(Fy, Fy)),
      _mm256_mul_ps(Fz, Fz));
  const __m256 lim1 = _mm256_cmp_ps(F2, _mm256_mul_ps(E, E), _CMP_GT_OQ);
  const __m256 corr1 = _mm256_blendv_ps(
      one, _mm256_div_ps(E, _mm256_sqrt_ps(_mm256_max_ps(F2, tiny))), lim1);

  const __m256 s1 = _mm256_div_ps(corr1, _mm256_max_ps(E, tiny));
  __m256 fx = _mm256_mul_ps(Fx, s1);
  __m256 fy = _mm256_mul_ps(Fy, s1);
  __m256 fz = _mm256_mul_ps(Fz, s1);
  const __m256 f2_raw = _mm256_add_ps(
      _mm256_add_ps(_mm256_mul_ps(fx, fx), _mm256_mul_ps(fy, fy)),
      _mm256_mul_ps(fz, fz));
  /* Round-off (or F2 and E^2 below the float range): normalise */
  const __m256 lim2 = _mm256_cmp_ps(f2_raw, one, _CMP_GT_OQ);
  const __m256 corr2 =
      _mm256_blendv_ps(one, _mm256_div_ps(one, _mm256_sqrt_ps(f2_raw)), lim2);
  fx = _mm256_mul_ps(fx, corr2);
  fy = _mm256_mul_ps(fy, corr2);
  fz = _mm256_mul_ps(fz, corr2);
  const __m256 f2 = _mm256_min_ps(f2_raw, one);
  const __m256 corr = _mm256_mul_ps(corr1, corr2);
  Fx = _mm256_mul_ps(Fx, corr);
  Fy = _mm256_mul_ps(Fy, corr);
  Fz = _mm256_mul_ps(Fz, corr);

  const __m256 n0 = _mm256_set1_ps(n_unit[0]);
  const __m256 n1 = _mm256_set1_ps(n_unit[1]);
  const __m256 n2 = _mm256_set1_ps(n_unit[2]);
  const __m256 fn = _mm256_add_ps(
      _mm256_add_ps(_mm256_mul_ps(fx, n0), _mm256_mul_ps(fy, n1)),
      _mm256_mul_ps(fz, n2));

  /* chi = (3 + 4 f^2) / (5 + 2 sqrt(4 - 3 f^2)) */
  const __m256 chi = _mm256_div_ps(
      _mm256_add_ps(three, _mm256_mul_ps(four, f2)),
      _mm256_add_ps(
          _mm256_set1_ps(5.f),
          _mm256_mul_ps(_mm256_set1_ps(2.f),
                        _mm256_sqrt_ps(_mm256_sub_ps(
                            four, _mm256_mul_ps(three, f2))))));
  const __m256 a_iso = _mm256_mul_ps(_mm256_mul_ps(half, _mm256_sub_ps(one, chi)), E);
  const __m256 a_dir = _mm256_div_ps(
      _mm256_mul_ps(
          _mm256_mul_ps(half, _mm256_sub_ps(_mm256_mul_ps(three, chi), one)),
          _mm256_mul_ps(E, fn)),
      _mm256_max_ps(f2, _mm256_set1_ps(1.e-30f)));

  _mm256_storeu_ps(Fn[0], _mm256_mul_ps(fn, E));
  _mm256_storeu_ps(Fn[1], _mm256_add_ps(_mm256_mul_ps(a_iso, n0),
                                        _mm256_mul_ps(a_dir, fx)));
  _mm256_storeu_ps(Fn[2], _mm256_add_ps(_mm256_mul_ps(a_iso, n1),
                                        _mm256_mul_ps(a_dir, fy)));
  _mm256_storeu_ps(Fn[3], _mm256_add_ps(_mm256_mul_ps(a_iso, n2),
                                        _mm256_mul_ps(a_dir, fz)));
  _mm256_storeu_ps(W[0], E);
  _mm256_storeu_ps(W[1], Fx);
  _mm256_storeu_ps(W[2], Fy);
  _mm256_storeu_ps(W[3], Fz);
#else
  for (int s = 0; s < RT_NFACE; s++) {
    /* E <= 0: no radiation */
    const int pos = (W[0][s] > 0.f);
    const float E = pos ? W[0][s] : 0.f;
    float Fx = pos ? W[1][s] : 0.f;
    float Fy = pos ? W[2][s] : 0.f;
    float Fz = pos ? W[3][s] : 0.f;

    /* |F| > c E: limit the flux to c E, before normalising (c E may be tiny
     * at a face while F is not, see rt_face_state_flux()) */
    const float F2 = Fx * Fx + Fy * Fy + Fz * Fz;
    const float corr1 = (F2 > E * E) ? E / sqrtf(fmaxf(F2, FLT_MIN)) : 1.f;

    const float E_inv = 1.f / fmaxf(E, FLT_MIN);
    float fx = Fx * corr1 * E_inv;
    float fy = Fy * corr1 * E_inv;
    float fz = Fz * corr1 * E_inv;
    const float f2_raw = fx * fx + fy * fy + fz * fz;
    /* Round-off (or F2 and E^2 below the float range): normalise */
    const float corr2 = (f2_raw > 1.f) ? 1.f / sqrtf(f2_raw) : 1.f;
    fx *= corr2;
    fy *= corr2;
    fz *= corr2;
    const float f2 = fminf(f2_raw, 1.f);
    const float corr = corr1 * corr2;
    Fx *= corr;
    Fy *= corr;
    Fz *= corr;

    const float fn = fx * n_unit[0] + fy * n_unit[1] + fz * n_unit[2];
    const float chi = (3.f + 4.f * f2) / (5.f + 2.f * sqrtf(4.f - 3.f * f2));
    const float a_iso = 0.5f * (1.f - chi) * E;
    const float a_dir = 0.5f * (3.f * chi - 1.f) * E * fn / fmaxf(f2, 1.e-30f);

    Fn[0][s] = fn * E;
    Fn[1][s] = a_iso * n_unit[0] + a_dir * fx;
    Fn[2][s] = a_iso * n_unit[1] + a_dir * fy;
    Fn[3][s] = a_iso * n_unit[2] + a_dir * fz;

    W[0][s] = E;
    W[1][s] = Fx;
    W[2][s] = Fy;
    W[3][s] = Fz;
  }
#endif
}

/**
 * @brief reset the mass fluxes of constituent species for a particle.
 *
 * @param p particle to work on.
 **/
__attribute__((always_inline)) INLINE static void rt_part_reset_mass_fluxes(
    struct part *restrict p) {
#ifdef GIZMO_MFV_SPH
  p->rt_data.mass_flux.HI = 0.f;
  p->rt_data.mass_flux.HII = 0.f;
  p->rt_data.mass_flux.HeI = 0.f;
  p->rt_data.mass_flux.HeII = 0.f;
  p->rt_data.mass_flux.HeIII = 0.f;
#endif
}

#endif /* SWIFT_KIARA_RT_FLUX_H */
