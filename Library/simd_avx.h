#ifndef LIBRARY_SIMD_AVX_H_
#define LIBRARY_SIMD_AVX_H_

#ifdef USE_SIMD

#include <math.h>
#include <immintrin.h>

/****************************************************
 * Define vector sizes and alignment
 ****************************************************/

#define VFLOAT 8
#define VDOUBLE 4
#define VALIGN 64

/****************************************************
 * Define vector types
 ****************************************************/

#define vecf __m256
#define vecd __m256d

/****************************************************
 * Define arithmetic operations
 ****************************************************/

/****************************************************
 * Arithmetic operations for float
 ****************************************************/

#define vadd_ps _mm256_add_ps
#define vsub_ps _mm256_sub_ps
#define vmul_ps _mm256_mul_ps
#define vdiv_ps _mm256_div_ps
#define vsqrt_ps _mm256_sqrt_ps
#define vrsqrt_ps _mm256_rsqrt_ps_
#define vfmadd_ps _mm256_fmadd_ps
#define vfmsub_ps _mm256_fmsub_ps
#define vfnmadd_ps _mm256_fnmadd_ps
#define vfnmsub_ps _mm256_fnmsub_ps

/****************************************************
 * Arithmetic operations for double
 ****************************************************/

#define vadd_pd _mm256_add_pd
#define vsub_pd _mm256_sub_pd
#define vmul_pd _mm256_mul_pd
#define vdiv_pd _mm256_div_pd
#define vsqrt_pd _mm256_sqrt_pd
#define vrsqrt_pd _mm256_rsqrt_pd_
#define vfmadd_pd _mm256_fmadd_pd
#define vfmsub_pd _mm256_fmsub_pd
#define vfnmadd_pd _mm256_fnmadd_pd
#define vfnmsub_pd _mm256_fnmsub_pd

/****************************************************
 * Define advanced arithmetic operations
 ****************************************************/

/****************************************************
 * Advanced arithmetic operations for float
 ****************************************************/

#define vsin_ps _mm256_sin_ps_
#define vcos_ps _mm256_cos_ps_
#define vsincos_ps _mm256_sincos_ps_
#define vexp_ps _mm256_exp_ps_

/****************************************************
 * Advanced arithmetic operations for double
 ****************************************************/

#define vsin_pd _mm256_sin_pd_
#define vcos_pd _mm256_cos_pd_
#define vsincos_pd _mm256_sincos_pd_
#define vexp_pd _mm256_exp_pd_

/****************************************************
 * Define load/store operations
 ****************************************************/

/****************************************************
 * Load/store operations for float
 ****************************************************/

#define vload_ps _mm256_load_ps
#define vload1_ps _mm256_broadcast_ss
#define vloadu_ps _mm256_loadu_ps
#define vset1_ps _mm256_set1_ps
#define vsetzero_ps _mm256_setzero_ps
#define vstore_ps _mm256_store_ps
#define vstoreu_ps _mm256_storeu_ps

/****************************************************
 * Load/store operations for double
 ****************************************************/

#define vload_pd _mm256_load_pd
#define vload1_pd _mm256_broadcast_sd
#define vloadu_pd _mm256_loadu_pd
#define vset1_pd _mm256_set1_pd
#define vsetzero_pd _mm256_setzero_pd
#define vstore_pd _mm256_store_pd
#define vstoreu_pd _mm256_storeu_pd

/****************************************************
 * Define compare operations
 ****************************************************/

/****************************************************
 * Define compare operations for float
 ****************************************************/

#define vcmpeq_ps(a,b) _mm256_cmp_ps(a,b, _CMP_EQ_OQ)
#define vcmpneq_ps(a,b) _mm256_cmp_ps(a,b, _CMP_NEQ_OQ)
#define vcmpge_ps(a,b) _mm256_cmp_ps(a,b, _CMP_GE_OQ)
#define vcmpgt_ps(a,b) _mm256_cmp_ps(a,b, _CMP_GT_OQ)
#define vcmpnge_ps(a,b) _mm256_cmp_ps(a,b, _CMP_NGE_OQ)
#define vcmpngt_ps(a,b) _mm256_cmp_ps(a,b, _CMP_NGT_OQ)
#define vcmple_ps(a,b) _mm256_cmp_ps(a,b, _CMP_LE_OQ)
#define vcmplt_ps(a,b) _mm256_cmp_ps(a,b, _CMP_LT_OQ)
#define vcmpnle_ps(a,b) _mm256_cmp_ps(a,b, _CMP_NLE_OQ)
#define vcmpnlt_ps(a,b) _mm256_cmp_ps(a,b, _CMP_NLT_OQ)

/****************************************************
 * Define compare operations for float
 ****************************************************/

#define vcmpeq_pd(a,b) _mm256_cmp_pd(a,b, _CMP_EQ_OQ)
#define vcmpneq_pd(a,b) _mm256_cmp_pd(a,b, _CMP_NEQ_OQ)
#define vcmpge_pd(a,b) _mm256_cmp_pd(a,b, _CMP_GE_OQ)
#define vcmpgt_pd(a,b) _mm256_cmp_pd(a,b, _CMP_GT_OQ)
#define vcmpnge_pd(a,b) _mm256_cmp_pd(a,b, _CMP_NGE_OQ)
#define vcmpngt_pd(a,b) _mm256_cmp_pd(a,b, _CMP_NGT_OQ)
#define vcmple_pd(a,b) _mm256_cmp_pd(a,b, _CMP_LE_OQ)
#define vcmplt_pd(a,b) _mm256_cmp_pd(a,b, _CMP_LT_OQ)
#define vcmpnle_pd(a,b) _mm256_cmp_pd(a,b, _CMP_NLE_OQ)
#define vcmpnlt_pd(a,b) _mm256_cmp_pd(a,b, _CMP_NLT_OQ)

/****************************************************
 * Definitions of bit operations
 ****************************************************/

/****************************************************
 * Definitions of bit operations for float
 ****************************************************/

#define vand_ps _mm256_and_ps
#define vandnot_ps _mm256_andnot_ps
#define vor_ps _mm256_or_ps
#define vxor_ps _mm256_xor_ps

/****************************************************
 * Definitions of bit operations for double
 ****************************************************/

#define vand_pd _mm256_and_pd
#define vandnot_pd _mm256_andnot_pd
#define vor_pd _mm256_or_pd
#define vxor_pd _mm256_xor_pd

/****************************************************
 * Define reductions of vector registers
 ****************************************************/

/****************************************************
 * Define reductions of vector registers for floats
 ****************************************************/

#define vreduce_ps _mm256_reduce_ps_

/****************************************************
 * Define reductions of vector registers for floats
 ****************************************************/

#define vreduce_pd _mm256_reduce_pd_

/****************************************************
 * Definition of little helper functions
 ****************************************************/

/****************************************************
 * Definition of little helfer functions for float
 ****************************************************/

#define vdot3_ps _mm256_dot3_ps_

/****************************************************
 * Definition of little helfer functions for double
 ****************************************************/

#define vdot3_pd _mm256_dot3_pd_

/****************************************************
 * Implementation of FMA function, if needed
 ****************************************************/

/****************************************************
 * Implementation of FMA function, if needed for float
 ****************************************************/

#ifndef __FMA__
static inline vecf _mm256_fmadd_ps(vecf a, vecf b, vecf c) {
  return vadd_ps(vmul_ps(a, b), c);
}

static inline vecf _mm256_fmsub_ps(vecf a, vecf b, vecf c) {
  return vsub_ps(vmul_ps(a, b), c);
}

static inline vecf _mm256_fnmadd_ps(vecf a, vecf b, vecf c) {
  return vsub_ps(c, vmul_ps(a, b));
}

static inline vecf _mm256_fnmsub_ps(vecf a, vecf b, vecf c) {
  return vmul_ps(vset1_ps(-1.0f), vadd_ps(vmul_ps(a, b), c));
}
#endif

/****************************************************
 * Implementation of FMA function, if needed for double
 ****************************************************/

#ifndef __FMA__
static inline vecd _mm256_fmadd_pd(vecd a, vecd b, vecd c) {
  return vadd_pd(vmul_pd(a, b), c);
}

static inline vecd _mm256_fmsub_pd(vecd a, vecd b, vecd c) {
  return vsub_pd(vmul_pd(a, b), c);
}

static inline vecd _mm256_fnmadd_pd(vecd a, vecd b, vecd c) {
  return vsub_pd(c, vmul_pd(a, b));
}

static inline vecd _mm256_fnmsub_pd(vecd a, vecd b, vecd c) {
  return vmul_pd(vset1_pd(-1.0), vadd_pd(vmul_pd(a, b), c));
}
#endif

/****************************************************
 * Implementations of basic vector functions
 ****************************************************/

static inline vecf _mm256_rsqrt_ps_(vecf x) {
  vecf x2, t1, t2, t3;

  x2 = vmul_ps(x, vset1_ps(0.5f));
  x = _mm256_rsqrt_ps(x);

  t1 = vmul_ps(x, x);
  t2 = vmul_ps(x2, x);
  t3 = vmul_ps(vset1_ps(1.5f), x);
  x = vfnmadd_ps(t1, t2, t3);

  return x;
}

static inline vecd _mm256_rsqrt_pd_(vecd x) {
  vecd x2, t1, t2, t3;
  vecd cd_15;

  cd_15 = vset1_pd(1.5);

  x2 = vmul_pd(x, vset1_pd(0.5));
  x = _mm256_cvtps_pd(_mm_rsqrt_ps(_mm256_cvtpd_ps(x)));

  t1 = vmul_pd(x, x);
  t2 = vmul_pd(x2, x);
  t3 = vmul_pd(cd_15, x);
  x = vfnmadd_pd(t1, t2, t3);

  t1 = vmul_pd(x, x);
  t2 = vmul_pd(x2, x);
  t3 = vmul_pd(cd_15, x);
  x = vfnmadd_pd(t1, t2, t3);

  return x;
}

/****************************************************
 * Implementations of basic vector functions for float
 ****************************************************/

static inline float _mm256_reduce_ps_(vecf x) {
  return x[0] + x[1] + x[2] + x[3] + x[4] + x[5] + x[6] + x[7];
}

/****************************************************
 * Implementations of basic vector functions for double
 ****************************************************/

static inline double _mm256_reduce_pd_(vecd x) {
  return x[0] + x[1] + x[2] + x[3];
}

/****************************************************
 * Implementation of advanced arithmetic operations
 ****************************************************/

/****************************************************
 * Implementation of advanced arithmetic operations for float
 ****************************************************/

static const float coeff_fcos[] = {
    1.0 / (1.0 * 2.0), 1.0 / (3.0 * 4.0), 1.0 / (5.0 * 6.0), 1.0 / (7.0 * 8.0),
    1.0 / (9.0 * 10.0), 1.0 / (11.0 * 12.0), 1.0 / (13.0 * 14.0), 1.0
        / (15.0 * 16.0), 1.0 / (17.0 * 18.0), 1.0 / (19.0 * 20.0), 1.0
        / (21.0 * 22.0), 1.0 / (23.0 * 24.0), 1.0 / (25.0 * 26.0), 1.0
        / (27.0 * 28.0) };
static const float coeff_fsin[] = {
    1.0 / (2.0 * 3.0), 1.0 / (4.0 * 5.0), 1.0 / (6.0 * 7.0), 1.0 / (8.0 * 9.0),
    1.0 / (10.0 * 11.0), 1.0 / (12.0 * 13.0), 1.0 / (14.0 * 15.0), 1.0
        / (16.0 * 17.0), 1.0 / (18.0 * 19.0), 1.0 / (20.0 * 21.0), 1.0
        / (22.0 * 23.0), 1.0 / (24.0 * 25.0), 1.0 / (26.0 * 27.0), 1.0
        / (28.0 * 29.0) };
static const int f_terms = 5;

static const float fpi = 3.14159265358979323846264338327950288L;
static const float fone_2pi = 1.0L / 6.28318530717958647692528676655900577L;

static inline vecf _mm256_cos_ps_(vecf xi) {
  vecf xi2, n, v_one, c, c_const, c_const2;
  int i;

  v_one = vset1_ps(1.0f);

  n = _mm256_floor_ps(vmul_ps(xi, vset1_ps(fone_2pi)));
  n = vfmadd_ps(vset1_ps(2.0f), n, v_one);

  xi = vfnmadd_ps(n, vset1_ps(fpi), xi);
  xi2 = vmul_ps(xi, xi);

  c_const2 = vmul_ps(xi2, vset1_ps(coeff_fcos[f_terms + 1]));
  c_const = vmul_ps(xi2, vset1_ps(coeff_fcos[f_terms]));
  c = vfnmadd_ps(c_const2, c_const, c_const);

  for (i = f_terms - 1; i >= 0; i--) {
    c_const = vmul_ps(xi2, vset1_ps(coeff_fcos[i]));
    c = vfnmadd_ps(c, c_const, c_const);
    ;
  }

  c = vsub_ps(c, v_one);

  return c;
}

static inline vecf _mm256_sin_ps_(vecf xi) {
  vecf xi2, n, v_one, s, s_const;
  int i;

  v_one = vset1_ps(1.0f);

  n = _mm256_floor_ps(vmul_ps(xi, vset1_ps(fone_2pi)));
  n = vfmadd_ps(vset1_ps(2.0f), n, v_one);

  xi = vfnmadd_ps(n, vset1_ps(fpi), xi);
  xi2 = vmul_ps(xi, xi);

  s = vmul_ps(xi2, vset1_ps(coeff_fsin[f_terms]));

  for (i = f_terms - 1; i >= 0; i--) {
    s_const = vmul_ps(xi2, vset1_ps(coeff_fsin[i]));
    s = vfnmadd_ps(s, s_const, s_const);
  }

  s = vfmsub_ps(xi, s, xi);

  return s;
}

static inline void _mm256_sincos_ps_(vecf xi, vecf *cp, vecf *sp) {
  vecf xi2, n, v_one, c, c_const, c_const2, s, s_const;
  int i;

  v_one = vset1_ps(1.0f);

  n = _mm256_floor_ps(vmul_ps(xi, vset1_ps(fone_2pi)));
  n = vfmadd_ps(vset1_ps(2.0f), n, v_one);

  xi = vfnmadd_ps(n, vset1_ps(fpi), xi);
  xi2 = vmul_ps(xi, xi);

  s = vmul_ps(xi2, vset1_ps(coeff_fsin[f_terms]));

  c_const2 = vmul_ps(xi2, vset1_ps(coeff_fcos[f_terms + 1]));
  c_const = vmul_ps(xi2, vset1_ps(coeff_fcos[f_terms]));
  c = vfnmadd_ps(c_const2, c_const, c_const);

  for (i = f_terms - 1; i >= 0; i--) {
    c_const = vmul_ps(xi2, vset1_ps(coeff_fcos[i]));
    s_const = vmul_ps(xi2, vset1_ps(coeff_fsin[i]));

    c = vfnmadd_ps(c, c_const, c_const);
    s = vfnmadd_ps(s, s_const, s_const);
  }

  *cp = vsub_ps(c, v_one);
  *sp = vfmsub_ps(xi, s, xi);
}

static const float coeff_fexp[] = {
    1.0 / 1.0, 1.0 / 2.0, 1.0 / 3.0, 1.0 / 4.0, 1.0 / 5.0, 1.0 / 6.0, 1.0 / 7.0,
    1.0 / 8.0, 1.0 / 9.0, 1.0 / 10.0, 1.0 / 11.0, 1.0 / 12.0, 1.0 / 13.0, 1.0
        / 14.0, 1.0 / 15.0 };

static const float flde = 1.4426950408889634e+00f;
static const float fln2 = 6.9314718055994529e-01f;

#ifdef __AVX2__
static const int f_expterms = 6;

static inline vecf _mm256_exp_ps_(vecf x) {
  int i;
  __m256 x1, y1, y, c_one, c;
  __m256i by1, pow2;

  c_one = vset1_ps(1.0f);

  pow2 = _mm256_cvtps_epi32(vmul_ps(vset1_ps(flde), x));
  x1 = vfnmadd_ps(vset1_ps(fln2), _mm256_cvtepi32_ps(pow2), x);

  c = vset1_ps(coeff_fexp[f_expterms - 1]);
  y1 = vfmadd_ps(c, x1, c_one);
  for (i = f_expterms - 2; i >= 0; i--) {
    c = vset1_ps(coeff_fexp[i]);
    y1 = vfmadd_ps(vmul_ps(c, x1), y1, c_one);
  }

  /* Multiply by 2^pow2 by adding to exponent in binary representation */
  by1 = _mm256_castps_si256(y1);
  by1 = _mm256_add_epi32(by1, _mm256_slli_epi32(pow2, 23));
  y = _mm256_castsi256_ps(by1);

  return y;
}

#else
static const int f_expterms = 7;

static inline vecf _mm256_exp_ps_(vecf x) {
  int i;
  __m256 x1, y1, y, c_one, c;
  __m128i pow2_1, pow2_2, by1_1, by1_2;
  __m256i by1, pow2;

  c_one = vset1_ps(1.0f);

  x1 = vmul_ps(vset1_ps(flde), x);

  pow2 = _mm256_cvttps_epi32(x1);
  pow2_1 = _mm256_extractf128_si256(pow2, 0);
  pow2_2 = _mm256_extractf128_si256(pow2, 1);

  x1 = vfnmadd_ps(vset1_ps(fln2), _mm256_cvtepi32_ps(pow2), x);

  c = vset1_ps(coeff_fexp[f_expterms - 1]);
  y1 = vfmadd_ps(c, x1, c_one);
  for (i = f_expterms - 2; i >= 0; i--) {
    c = vset1_ps(coeff_fexp[i]);
    y1 = vfmadd_ps(vmul_ps(c, x1), y1, c_one);
  }

  /* Multiply by 2^pow2 by adding to exponent in binary representation */
  by1_1 = _mm_castps_si128(_mm256_extractf128_ps(y1, 0));
  by1_1 = _mm_add_epi32(by1_1, _mm_slli_epi32(pow2_1, 23));
  by1_2 = _mm_castps_si128(_mm256_extractf128_ps(y1, 1));
  by1_2 = _mm_add_epi32(by1_2, _mm_slli_epi32(pow2_2, 23));
  by1 = _mm256_setzero_si256();
  by1 = _mm256_insertf128_si256(by1, by1_1, 0);
  by1 = _mm256_insertf128_si256(by1, by1_2, 1);
  y = _mm256_castsi256_ps(by1);

  return y;
}
#endif

/****************************************************
 * Implementation of advanced arithmetic operations for double
 ****************************************************/

static const double coeff_dcos[] = {
    1.0 / (1.0 * 2.0), 1.0 / (3.0 * 4.0), 1.0 / (5.0 * 6.0), 1.0 / (7.0 * 8.0),
    1.0 / (9.0 * 10.0), 1.0 / (11.0 * 12.0), 1.0 / (13.0 * 14.0), 1.0
        / (15.0 * 16.0), 1.0 / (17.0 * 18.0), 1.0 / (19.0 * 20.0), 1.0
        / (21.0 * 22.0), 1.0 / (23.0 * 24.0), 1.0 / (25.0 * 26.0), 1.0
        / (27.0 * 28.0) };
static const double coeff_dsin[] = {
    1.0 / (2.0 * 3.0), 1.0 / (4.0 * 5.0), 1.0 / (6.0 * 7.0), 1.0 / (8.0 * 9.0),
    1.0 / (10.0 * 11.0), 1.0 / (12.0 * 13.0), 1.0 / (14.0 * 15.0), 1.0
        / (16.0 * 17.0), 1.0 / (18.0 * 19.0), 1.0 / (20.0 * 21.0), 1.0
        / (22.0 * 23.0), 1.0 / (24.0 * 25.0), 1.0 / (26.0 * 27.0), 1.0
        / (28.0 * 29.0) };
static const int d_terms = 10;

static const double dpi = 3.14159265358979323846264338327950288L;
static const double done_2pi = 1.0L / 6.28318530717958647692528676655900577L;

static inline vecd _mm256_cos_pd_(vecd xi) {
  vecd xi2, n, v_one, c, c_const, c_const2;
  int i;

  v_one = vset1_pd(1.0);

  n = _mm256_floor_pd(vmul_pd(xi, vset1_pd(done_2pi)));
  n = vfmadd_pd(vset1_pd(2.0), n, v_one);

  xi = vfnmadd_pd(n, vset1_pd(dpi), xi);
  xi2 = vmul_pd(xi, xi);

  c_const2 = vmul_pd(xi2, vset1_pd(coeff_dcos[d_terms + 1]));
  c_const = vmul_pd(xi2, vset1_pd(coeff_dcos[d_terms]));
  c = vfnmadd_pd(c_const2, c_const, c_const);

  for (i = d_terms - 1; i >= 0; i--) {
    c_const = vmul_pd(xi2, vset1_pd(coeff_dcos[i]));
    c = vfnmadd_pd(c, c_const, c_const);
  }

  c = vsub_pd(c, v_one);

  return c;
}

static inline vecd _mm256_sin_pd_(vecd xi) {
  vecd xi2, n, v_one, s, s_const;
  int i;

  v_one = vset1_pd(1.0);

  n = _mm256_floor_pd(vmul_pd(xi, vset1_pd(done_2pi)));
  n = vfmadd_pd(vset1_pd(2.0), n, v_one);

  xi = vfnmadd_pd(n, vset1_pd(dpi), xi);
  xi2 = vmul_pd(xi, xi);

  s = vmul_pd(xi2, vset1_pd(coeff_dsin[d_terms]));

  for (i = d_terms - 1; i >= 0; i--) {
    s_const = vmul_pd(xi2, vset1_pd(coeff_dsin[i]));
    s = vfnmadd_pd(s, s_const, s_const);
  }

  s = vfmsub_pd(xi, s, xi);

  return s;
}

static inline void _mm256_sincos_pd_(vecd xi, vecd *cp, vecd *sp) {
  vecd xi2, n, v_one, c, c_const, c_const2, s, s_const;
  int i;

  v_one = vset1_pd(1.0);

  n = _mm256_floor_pd(vmul_pd(xi, vset1_pd(done_2pi)));
  n = vfmadd_pd(vset1_pd(2.0), n, v_one);

  xi = vfnmadd_pd(n, vset1_pd(dpi), xi);
  xi2 = vmul_pd(xi, xi);

  s = vmul_pd(xi2, vset1_pd(coeff_dsin[d_terms]));

  c_const2 = vmul_pd(xi2, vset1_pd(coeff_dcos[d_terms + 1]));
  c_const = vmul_pd(xi2, vset1_pd(coeff_dcos[d_terms]));
  c = vfnmadd_pd(c_const2, c_const, c_const);

  for (i = d_terms - 1; i >= 0; i--) {
    c_const = vmul_pd(xi2, vset1_pd(coeff_dcos[i]));
    s_const = vmul_pd(xi2, vset1_pd(coeff_dsin[i]));

    c = vfnmadd_pd(c, c_const, c_const);
    s = vfnmadd_pd(s, s_const, s_const);
  }

  *cp = vsub_pd(c, v_one);
  *sp = vfmsub_pd(xi, s, xi);
}

static const double coeff_dexp[] = {
    1.0 / 1.0, 1.0 / 2.0, 1.0 / 3.0, 1.0 / 4.0, 1.0 / 5.0, 1.0 / 6.0, 1.0 / 7.0,
    1.0 / 8.0, 1.0 / 9.0, 1.0 / 10.0, 1.0 / 11.0, 1.0 / 12.0, 1.0 / 13.0, 1.0
        / 14.0, 1.0 / 15.0 };

static const double dlde = 1.4426950408889634e+00;
static const double dln2 = 6.9314718055994529e-01;

#ifdef __AVX2__
static const int d_expterms = 11;

static inline vecd _mm256_exp_pd_(vecd x) {
  int i;
  __m256d x1, y1, y, c_one, c;
  __m128i pow2;
  __m256i by1;

  c_one = vset1_pd(1.0);

  pow2 = _mm256_cvtpd_epi32(vmul_pd(vset1_pd(dlde), x));
  x1 = vfnmadd_pd(vset1_pd(dln2), _mm256_cvtepi32_pd(pow2), x);

  c = vset1_pd(coeff_dexp[d_expterms - 1]);
  y1 = vfmadd_pd(c, x1, c_one);
  for (i = d_expterms - 2; i >= 0; i--) {
    c = vset1_pd(coeff_dexp[i]);
    y1 = vfmadd_pd(vmul_pd(c, x1), y1, c_one);
  }

  /* Multiply by 2^pow2 by adding to exponent in binary representation */
  by1 = _mm256_castpd_si256(y1);
  by1 = _mm256_add_epi64(by1,
      _mm256_slli_epi64(_mm256_cvtepi32_epi64(pow2), 52));
  y = _mm256_castsi256_pd(by1);

  return y;
}

#else
static const int d_expterms = 14;

static inline vecd _mm256_exp_pd_(vecd x) {
  int i;
  __m256d x1, y1, y, c_one, c;
  __m128i pow2, pow2_1, pow2_2, by1_1, by1_2;
  __m256i by1;

  c_one = vset1_pd(1.0);

  x1 = vmul_pd(vset1_pd(dlde), x);
  pow2_1 = _mm_cvttpd_epi32(_mm256_extractf128_pd(x1, 0));
  pow2_2 = _mm_cvttpd_epi32(_mm256_extractf128_pd(x1, 1));
  pow2 = _mm_setzero_si128();
  pow2 = _mm_unpacklo_epi32(pow2_1, pow2_2);
  pow2 = _mm_shuffle_epi32(pow2, 216);

  x1 = vfnmadd_pd(vset1_pd(dln2), _mm256_cvtepi32_pd(pow2), x);

  c = vset1_pd(coeff_dexp[d_expterms - 1]);
  y1 = vfmadd_pd(c, x1, c_one);
  for (i = d_expterms - 2; i >= 0; i--) {
    c = vset1_pd(coeff_dexp[i]);
    y1 = vfmadd_pd(vmul_pd(c, x1), y1, c_one);
  }

  /* Multiply by 2^pow2 by adding to exponent in binary representation */
  by1_1 = _mm_castpd_si128(_mm256_extractf128_pd(y1, 0));
  by1_1 = _mm_add_epi64(by1_1, _mm_slli_epi64(_mm_cvtepi32_epi64(pow2_1), 52));
  by1_2 = _mm_castpd_si128(_mm256_extractf128_pd(y1, 1));
  by1_2 = _mm_add_epi64(by1_2, _mm_slli_epi64(_mm_cvtepi32_epi64(pow2_2), 52));
  by1 = _mm256_setzero_si256();
  by1 = _mm256_insertf128_si256(by1, by1_1, 0);
  by1 = _mm256_insertf128_si256(by1, by1_2, 1);
  y = _mm256_castsi256_pd(by1);

  return y;
}
#endif

/****************************************************
 * Implementation of little helper functions
 ****************************************************/

/****************************************************
 * Implementation of little helper functions for float
 ****************************************************/

static inline vecf _mm256_dot3_ps_(vecf x[3], vecf y[3]) {
  vecf res;

  res = vmul_ps(x[0], y[0]);
  res = vfmadd_ps(x[1], y[1], res);
  res = vfmadd_ps(x[2], y[2], res);

  return res;
}

/****************************************************
 * Implementation of little helper functions for double
 ****************************************************/

static inline vecd _mm256_dot3_pd_(vecd x[3], vecd y[3]) {
  vecd res;

  res = vmul_pd(x[0], y[0]);
  res = vfmadd_pd(x[1], y[1], res);
  res = vfmadd_pd(x[2], y[2], res);

  return res;
}

#endif

#endif /* LIBRARY_SIMD_AVX_H_ */
