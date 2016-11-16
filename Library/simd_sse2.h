#ifndef LIBRARY_SIMD_SSE2_H_
#define LIBRARY_SIMD_SSE2_H_

#ifdef USE_SIMD

#include <math.h>
#include <immintrin.h>

/****************************************************
 * Define vector sizes and alignment
 ****************************************************/

#define VFLOAT 4
#define VDOUBLE 2
#define VALIGN 64

/****************************************************
 * Define vector types
 ****************************************************/

#define vecf __m128
#define vecd __m128d

/****************************************************
 * Define arithmetic operations
 ****************************************************/

/****************************************************
 * Arithmetic operations for float
 ****************************************************/

#define vadd_ps _mm_add_ps
#define vsub_ps _mm_sub_ps
#define vmul_ps _mm_mul_ps
#define vdiv_ps _mm_div_ps
#define vsqrt_ps _mm_sqrt_ps
#define vrsqrt_ps _mm_rsqrt_ps_
#define vfmadd_ps _mm_fmadd_ps
#define vfmsub_ps _mm_fmsub_ps
#define vfnmadd_ps _mm_fnmadd_ps
#define vfnmsub_ps _mm_fnmsub_ps

/****************************************************
 * Arithmetic operations for double
 ****************************************************/

#define vadd_pd _mm_add_pd
#define vsub_pd _mm_sub_pd
#define vmul_pd _mm_mul_pd
#define vdiv_pd _mm_div_pd
#define vsqrt_pd _mm_sqrt_pd
#define vrsqrt_pd _mm_rsqrt_pd_
#define vfmadd_pd _mm_fmadd_pd
#define vfmsub_pd _mm_fmsub_pd
#define vfnmadd_pd _mm_fnmadd_pd
#define vfnmsub_pd _mm_fnmsub_pd

/****************************************************
 * Define advanced arithmetic operations
 ****************************************************/

/****************************************************
 * Advanced arithmetic operations for float
 ****************************************************/

#define vsin_ps _mm_sin_ps_
#define vcos_ps _mm_cos_ps_
#define vsincos_ps _mm_sincos_ps_
#define vexp_ps _mm_exp_ps_

/****************************************************
 * Advanced arithmetic operations for double
 ****************************************************/

#define vsin_pd _mm_sin_pd_
#define vcos_pd _mm_cos_pd_
#define vsincos_pd _mm_sincos_pd_
#define vexp_pd _mm_exp_pd_

/****************************************************
 * Define load/store operations
 ****************************************************/

/****************************************************
 * Load/store operations for float
 ****************************************************/

#define vload_ps _mm_load_ps
#define vload1_ps _mm_load1_ps
#define vloadu_ps _mm_loadu_ps
#define vset1_ps _mm_set1_ps
#define vsetzero_ps _mm_setzero_ps
#define vstore_ps _mm_store_ps
#define vstoreu_ps _mm_storeu_ps

/****************************************************
 * Load/store operations for double
 ****************************************************/

#define vload_pd _mm_load_pd
#define vload1_pd _mm_load1_pd
#define vloadu_pd _mm_loadu_pd
#define vset1_pd _mm_set1_pd
#define vsetzero_pd _mm_setzero_pd
#define vstore_pd _mm_store_pd
#define vstoreu_pd _mm_storeu_pd

/****************************************************
 * Define compare operations
 ****************************************************/

/****************************************************
 * Define compare operations for float
 ****************************************************/

#define vcmpeq_ps(a,b) _mm_cmpeq_ps(a,b)
#define vcmpneq_ps(a,b) _mm_cmpneq_ps(a,b)
#define vcmpge_ps(a,b) _mm_cmpge_ps(a,b)
#define vcmpgt_ps(a,b) _mm_cmpgt_ps(a,b)
#define vcmpnge_ps(a,b) _mm_cmpnge_ps(a,b)
#define vcmpngt_ps(a,b) _mm_cmpngt_ps(a,b)
#define vcmple_ps(a,b) _mm_cmple_ps(a,b)
#define vcmplt_ps(a,b) _mm_cmplt_ps(a,b)
#define vcmpnle_ps(a,b) _mm_cmpnle_ps(a,b)
#define vcmpnlt_ps(a,b) _mm_cmpnlt_ps(a,b)

/****************************************************
 * Define compare operations for float
 ****************************************************/

#define vcmpeq_pd(a,b) _mm_cmpeq_pd(a,b)
#define vcmpneq_pd(a,b) _mm_cmpneq_pd(a,b)
#define vcmpge_pd(a,b) _mm_cmpge_pd(a,b)
#define vcmpgt_pd(a,b) _mm_cmpgt_pd(a,b)
#define vcmpnge_pd(a,b) _mm_cmpnge_pd(a,b)
#define vcmpngt_pd(a,b) _mm_cmpngt_pd(a,b)
#define vcmple_pd(a,b) _mm_cmple_pd(a,b)
#define vcmplt_pd(a,b) _mm_cmplt_pd(a,b)
#define vcmpnle_pd(a,b) _mm_cmpnle_pd(a,b)
#define vcmpnlt_pd(a,b) _mm_cmpnlt_pd(a,b)

/****************************************************
 * Definitions of bit operations
 ****************************************************/

/****************************************************
 * Definitions of bit operations for float
 ****************************************************/

#define vand_ps _mm_and_ps
#define vandnot_ps _mm_andnot_ps
#define vor_ps _mm_or_ps
#define vxor_ps _mm_xor_ps

/****************************************************
 * Definitions of bit operations for double
 ****************************************************/

#define vand_pd _mm_and_pd
#define vandnot_pd _mm_andnot_pd
#define vor_pd _mm_or_pd
#define vxor_pd _mm_xor_pd

/****************************************************
 * Define reductions of vector registers
 ****************************************************/

/****************************************************
 * Define reductions of vector registers for floats
 ****************************************************/

#define vreduce_ps _mm_reduce_ps_

/****************************************************
 * Define reductions of vector registers for floats
 ****************************************************/

#define vreduce_pd _mm_reduce_pd_

/****************************************************
 * Definition of little helper functions
 ****************************************************/

/****************************************************
 * Definition of little helfer functions for float
 ****************************************************/

#define vdot3_ps _mm_dot3_ps_

/****************************************************
 * Definition of little helfer functions for double
 ****************************************************/

#define vdot3_pd _mm_dot3_pd_

/****************************************************
 * Implementation of FMA function, if needed
 ****************************************************/

/****************************************************
 * Implementation of FMA function, if needed for float
 ****************************************************/

#ifndef __FMA__
static inline vecf _mm_fmadd_ps(vecf a, vecf b, vecf c) {
  return vadd_ps(vmul_ps(a, b), c);
}

static inline vecf _mm_fmsub_ps(vecf a, vecf b, vecf c) {
  return vsub_ps(vmul_ps(a, b), c);
}

static inline vecf _mm_fnmadd_ps(vecf a, vecf b, vecf c) {
  return vsub_ps(c, vmul_ps(a, b));
}

static inline vecf _mm_fnmsub_ps(vecf a, vecf b, vecf c) {
  return vmul_ps(vset1_ps(-1.0f), vadd_ps(vmul_ps(a, b), c));
}
#endif

/****************************************************
 * Implementation of FMA function, if needed for double
 ****************************************************/

#ifndef __FMA__
static inline vecd _mm_fmadd_pd(vecd a, vecd b, vecd c) {
  return vadd_pd(vmul_pd(a, b), c);
}

static inline vecd _mm_fmsub_pd(vecd a, vecd b, vecd c) {
  return vsub_pd(vmul_pd(a, b), c);
}

static inline vecd _mm_fnmadd_pd(vecd a, vecd b, vecd c) {
  return vsub_pd(c, vmul_pd(a, b));
}

static inline vecd _mm_fnmsub_pd(vecd a, vecd b, vecd c) {
  return vmul_pd(vset1_pd(-1.0), vadd_pd(vmul_pd(a, b), c));
}
#endif

#ifndef __SSE4_1__
static inline vecf _mm_floor_ps_(vecf a) {
  return _mm_cvtepi32_ps(_mm_cvttps_epi32(a));
}

static inline vecd _mm_floor_pd_(vecd a) {
  return _mm_cvtepi32_pd(_mm_cvttpd_epi32(a));
}
#else
static inline vecf _mm_floor_ps_(vecf a) {
  return _mm_floor_ps(a);
}

static inline vecd _mm_floor_pd_(vecd a) {
  return _mm_floor_pd(a);
}
#endif

#ifndef __SSE4_1__
static inline __m128i _mm_cvtepi32_epi64_(__m128i a) {
  __m128i sel = _mm_set1_epi32(0x80000000);
  __m128i mask;
  __m128i a1;

  __m128i b, c;

  mask = _mm_cmpeq_epi32(_mm_and_si128(a, sel), sel);
  a1 = _mm_andnot_si128(sel, a);

  b = _mm_or_si128(_mm_and_si128(mask, a1), _mm_andnot_si128(mask, a));
  c = _mm_or_si128(_mm_and_si128(mask, sel),
      _mm_andnot_si128(mask, _mm_setzero_si128()));
  b = _mm_unpacklo_epi32(b,c);

  return b;
}
#else
static inline __m128i _mm_cvtepi32_epi64_(__m128i a) {
  return _mm_cvtepi32_epi64(a);
}
#endif

/****************************************************
 * Implementations of basic vector functions
 ****************************************************/

static inline vecf _mm_rsqrt_ps_(vecf x) {
  vecf x2, t1, t2, t3;

  x2 = vmul_ps(x, vset1_ps(0.5f));
  x = _mm_rsqrt_ps(x);

  t1 = vmul_ps(x, x);
  t2 = vmul_ps(x2, x);
  t3 = vmul_ps(vset1_ps(1.5f), x);
  x = vfnmadd_ps(t1, t2, t3);

  return x;
}

static inline vecd _mm_rsqrt_pd_(vecd x) {
  vecd x2, t1, t2, t3;
  vecd cd_15;

  cd_15 = vset1_pd(1.5);

  x2 = vmul_pd(x, vset1_pd(0.5));
  x = _mm_cvtps_pd(_mm_rsqrt_ps(_mm_cvtpd_ps(x)));

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

static inline float _mm_reduce_ps_(vecf x) {
  return x[0] + x[1] + x[2] + x[3];
}

/****************************************************
 * Implementations of basic vector functions for double
 ****************************************************/

static inline double _mm_reduce_pd_(vecd x) {
  return x[0] + x[1];
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

static inline vecf _mm_cos_ps_(vecf xi) {
  vecf xi2, n, v_one, c, c_const, c_const2;
  int i;

  v_one = vset1_ps(1.0f);

  n = _mm_floor_ps_(vmul_ps(xi, vset1_ps(fone_2pi)));
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

static inline vecf _mm_sin_ps_(vecf xi) {
  vecf xi2, n, v_one, s, s_const;
  int i;

  v_one = vset1_ps(1.0f);

  n = _mm_floor_ps_(vmul_ps(xi, vset1_ps(fone_2pi)));
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

static inline void _mm_sincos_ps_(vecf xi, vecf *cp, vecf *sp) {
  vecf xi2, n, v_one, c, c_const, c_const2, s, s_const;
  int i;

  v_one = vset1_ps(1.0f);

  n = _mm_floor_ps_(vmul_ps(xi, vset1_ps(fone_2pi)));
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

static const int f_expterms = 7;

static inline vecf _mm_exp_ps_(vecf x) {
  int i;
  __m128 x1, y1, y, c_one, c;
  __m128i by1, pow2;

  c_one = vset1_ps(1.0f);

  x1 = vmul_ps(vset1_ps(flde), x);

  pow2 = _mm_cvttps_epi32(x1);

  x1 = vfnmadd_ps(vset1_ps(fln2), _mm_cvtepi32_ps(pow2), x);

  c = vset1_ps(coeff_fexp[f_expterms - 1]);
  y1 = vfmadd_ps(c, x1, c_one);
  for (i = f_expterms - 2; i >= 0; i--) {
    c = vset1_ps(coeff_fexp[i]);
    y1 = vfmadd_ps(vmul_ps(c, x1), y1, c_one);
  }

  /* Multiply by 2^pow2 by adding to exponent in binary representation */
  by1 = _mm_castps_si128(y1);
  by1 = _mm_add_epi32(by1, _mm_slli_epi32(pow2, 23));
  ;
  y = _mm_castsi128_ps(by1);

  return y;
}

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

static inline vecd _mm_cos_pd_(vecd xi) {
  vecd xi2, n, v_one, c, c_const, c_const2;
  int i;

  v_one = vset1_pd(1.0);

  n = _mm_floor_pd_(vmul_pd(xi, vset1_pd(done_2pi)));
  n = vfmadd_pd(vset1_pd(2.0), n, v_one);

  xi = vfnmadd_pd(n, vset1_pd(dpi), xi);
  xi2 = vmul_pd(xi, xi);

  c_const2 = vmul_pd(xi2, vset1_pd(coeff_dcos[d_terms + 1]));
  c_const = vmul_pd(xi2, vset1_pd(coeff_dcos[d_terms]));
  c = vfnmadd_pd(c_const2, c_const, c_const);

  for (i = d_terms - 1; i >= 0; i--) {
    c_const = vmul_pd(xi2, vset1_pd(coeff_dcos[i]));
    c = vfnmadd_pd(c, c_const, c_const);
    ;
  }

  c = vsub_pd(c, v_one);

  return c;
}

static inline vecd _mm_sin_pd_(vecd xi) {
  vecd xi2, n, v_one, s, s_const;
  int i;

  v_one = vset1_pd(1.0);

  n = _mm_floor_pd_(vmul_pd(xi, vset1_pd(done_2pi)));
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

static inline void _mm_sincos_pd_(vecd xi, vecd *cp, vecd *sp) {
  vecd xi2, n, v_one, c, c_const, c_const2, s, s_const;
  int i;

  v_one = vset1_pd(1.0);

  n = _mm_floor_pd_(vmul_pd(xi, vset1_pd(done_2pi)));
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

static const int d_expterms = 14;

static inline vecd _mm_exp_pd_(vecd x) {
  int i;
  __m128d x1, y1, y, c_one, c;
  __m128i pow2, by1;

  c_one = vset1_pd(1.0);

  x1 = vmul_pd(vset1_pd(dlde), x);
  pow2 = _mm_cvttpd_epi32(x1);

  x1 = vfnmadd_pd(vset1_pd(dln2), _mm_cvtepi32_pd(pow2), x);

  c = vset1_pd(coeff_dexp[d_expterms - 1]);
  y1 = vfmadd_pd(c, x1, c_one);
  for (i = d_expterms - 2; i >= 0; i--) {
    c = vset1_pd(coeff_dexp[i]);
    y1 = vfmadd_pd(vmul_pd(c, x1), y1, c_one);
  }

  /* Multiply by 2^pow2 by adding to exponent in binary representation */
  by1 = _mm_castpd_si128(y1);
  by1 = _mm_add_epi64(by1, _mm_slli_epi64(_mm_cvtepi32_epi64_(pow2), 52));
  y = _mm_castsi128_pd(by1);

  return y;
}

/****************************************************
 * Implementation of little helper functions
 ****************************************************/

/****************************************************
 * Implementation of little helper functions for float
 ****************************************************/

static inline vecf _mm_dot3_ps_(vecf x[3], vecf y[3]) {
  vecf res;

  res = vmul_ps(x[0], y[0]);
  res = vfmadd_ps(x[1], y[1], res);
  res = vfmadd_ps(x[2], y[2], res);

  return res;
}

/****************************************************
 * Implementation of little helper functions for double
 ****************************************************/

static inline vecd _mm_dot3_pd_(vecd x[3], vecd y[3]) {
  vecd res;

  res = vmul_pd(x[0], y[0]);
  res = vfmadd_pd(x[1], y[1], res);
  res = vfmadd_pd(x[2], y[2], res);

  return res;
}

#endif

#endif /* LIBRARY_SIMD_AVX_H_ */
