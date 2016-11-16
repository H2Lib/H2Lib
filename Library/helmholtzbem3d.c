/* ------------------------------------------------------------
 This is the file "helmholtzbem3d.c" of the H2Lib package.
 All rights reserved, Sven Christophersen 2015
 ------------------------------------------------------------ */

/**
 * @file helmholtzbem3d.c
 * @author Sven Christophersen
 * @date 2015
 */

#ifdef USE_COMPLEX

/* C STD LIBRARY */
/* CORE 0 */
/* CORE 1 */
/* CORE 2 */
/* CORE 3 */
/* SIMPLE */
/* PARTICLES */
/* BEM */
#include "helmholtzbem3d.h"

#define KERNEL_CONST_HELMHOLTZBEM3D 0.0795774715459476679

static inline field
slp_kernel_helmholtzbem3d(const real * x, const real * y,
			  const real * nx, const real * ny, void *data)
{
  pcbem3d   bem = (pcbem3d) data;
  real      k_real = REAL(bem->k);
  real      k_imag = -IMAG(bem->k);
  real      dist[3];
  real      norm, norm2, rnorm;

  field     res;

  (void) nx;
  (void) ny;

  dist[0] = x[0] - y[0];
  dist[1] = x[1] - y[1];
  dist[2] = x[2] - y[2];
  norm2 = REAL_NORMSQR3(dist[0], dist[1], dist[2]);
  rnorm = REAL_RSQRT(norm2);

  norm = norm2 * rnorm;
  if (k_imag != 0.0) {
    rnorm *= REAL_EXP(k_imag * norm);
  }
  norm = k_real * norm;
  res = (rnorm * REAL_COS(norm)) + (REAL_SIN(norm) * rnorm) * I;

  return res;
}

#ifdef USE_SIMD
static inline void
slp_kernel_simd_helmholtzbem3d(const vreal * x,
			       const vreal * y, const vreal * nx,
			       const vreal * ny, void *data, vreal * res_re,
			       vreal * res_im)
{
  pcbem3d   bem = (pcbem3d) data;
  vreal     k_real = vset1(REAL(bem->k));
  vreal     k_imag = vset1(-IMAG(bem->k));
  vreal     dist[3];
  vreal     norm, norm2, rnorm, c, s;
  (void) nx;
  (void) ny;

  dist[0] = vsub(x[0], y[0]);
  dist[1] = vsub(x[1], y[1]);
  dist[2] = vsub(x[2], y[2]);

  norm2 = vdot3(dist, dist);
  rnorm = vrsqrt(norm2);

  norm = vmul(norm2, rnorm);
  if (IMAG(bem->k) != 0.0) {
    rnorm = vmul(rnorm, vexp(vmul(k_imag, norm)));
  }

  norm = vmul(k_real, norm);

  vsincos(norm, &c, &s);

  *res_re = vmul(rnorm, c);
  *res_im = vmul(rnorm, s);
}
#endif

static field
slp_kernel_directional_helmholtzbem3d(const real * x,
				      const real * y, const real * nx,
				      const real * ny, pcreal dir, void *data)
{
  real      d[3], angle, kr, ki, inorm, norm, product;
  pcbem3d   bem = (pcbem3d) data;
  real      k_real = REAL(bem->k);
  real      k_imag = -IMAG(bem->k);

  (void) nx;
  (void) ny;

  d[0] = x[0] - y[0];
  d[1] = x[1] - y[1];
  d[2] = x[2] - y[2];

  /* Norm |x-y| */
  norm = REAL_NORM3(d[0], d[1], d[2]);

  /* Product <x-y, c> */
  product = REAL_DOT3(d, dir);

  /* Quotient 1 / (|x-y|) */
  inorm = 1.0 / norm;

  if (k_imag != 0.0) {
    /* Angle k_imag |x-y| */
    angle = k_imag * (norm - product);

    /* exp(angle) */
    inorm *= REAL_EXP(angle);
  }

  /* Angle k_real |x-y| */
  angle = k_real * (norm - product);

  /* Real part cos(k |x-y|) / (4 pi |x-y|) */
  kr = inorm * REAL_COS(angle);

  /* Imaginary part sin(k |x-y|) / (4 pi |x-y|) */
  ki = inorm * REAL_SIN(angle);

  return kr + I * ki;
}

static inline field
dlp_kernel_helmholtzbem3d(const real * x, const real * y,
			  const real * nx, const real * ny, void *data)
{
  pcbem3d   bem = (pcbem3d) data;
  real      k_real = REAL(bem->k);
  real      k_imag = -IMAG(bem->k);
  real      dist[3];
  real      norm, norm2, rnorm, s, c, re;

  field     res;

  (void) nx;

  dist[0] = x[0] - y[0];
  dist[1] = x[1] - y[1];
  dist[2] = x[2] - y[2];
  norm2 = REAL_NORMSQR3(dist[0], dist[1], dist[2]);
  rnorm = REAL_RSQRT(norm2);

  re = (rnorm * rnorm * rnorm) * DOT3(dist, ny);
  norm = norm2 * rnorm;
  if (k_imag != 0.0) {
    re *= REAL_EXP(k_imag * norm);
  }
  norm = k_real * norm;
  s = REAL_SIN(norm);
  c = REAL_COS(norm);
  res = ((c + norm * s) * re + (s - c * norm) * re * I);

  return res;
}

#ifdef USE_SIMD
static inline void
dlp_kernel_simd_helmholtzbem3d(const vreal * x,
			       const vreal * y, const vreal * nx,
			       const vreal * ny, void *data, vreal * res_re,
			       vreal * res_im)
{
  pcbem3d   bem = (pcbem3d) data;
  vreal     k_real = vset1(REAL(bem->k));
  vreal     k_imag = vset1(-IMAG(bem->k));
  vreal     dist[3];
  vreal     norm, norm2, rnorm, s, c;

  (void) nx;

  dist[0] = vsub(x[0], y[0]);
  dist[1] = vsub(x[1], y[1]);
  dist[2] = vsub(x[2], y[2]);
  norm2 = vdot3(dist, dist);
  rnorm = vrsqrt(norm2);

  norm = vmul(norm2, rnorm);
  rnorm = vmul(vmul(rnorm, vmul(rnorm, rnorm)), vdot3(dist, (vreal *) ny));
  if (IMAG(bem->k) != 0.0) {
    rnorm = vmul(rnorm, vexp(vmul(k_imag, norm)));
  }
  norm = vmul(k_real, norm);

  vsincos(norm, &c, &s);

  *res_re = vmul(vfmadd(norm, s, c), rnorm);
  *res_im = vmul(vfnmadd(norm, c, s), rnorm);

}
#endif

static inline field
adlp_kernel_helmholtzbem3d(const real * x, const real * y,
			   const real * nx, const real * ny, void *data)
{
  pcbem3d   bem = (pcbem3d) data;
  real      k_real = REAL(bem->k);
  real      k_imag = -IMAG(bem->k);
  real      dist[3];
  real      norm, norm2, rnorm, s, c, re;

  field     res;

  (void) ny;

  dist[0] = y[0] - x[0];
  dist[1] = y[1] - x[1];
  dist[2] = y[2] - x[2];
  norm2 = REAL_NORMSQR3(dist[0], dist[1], dist[2]);
  rnorm = REAL_RSQRT(norm2);

  re = (rnorm * rnorm * rnorm) * DOT3(dist, nx);
  norm = norm2 * rnorm;
  if (k_imag != 0.0) {
    re *= REAL_EXP(k_imag * norm);
  }
  norm = k_real * norm;
  s = REAL_SIN(norm);
  c = REAL_COS(norm);
  res = ((c + norm * s) * re + (s - c * norm) * re * I);

  return res;
}

#ifdef USE_SIMD
static inline void
adlp_kernel_simd_helmholtzbem3d(const vreal * x,
				const vreal * y, const vreal * nx,
				const vreal * ny, void *data, vreal * res_re,
				vreal * res_im)
{
  pcbem3d   bem = (pcbem3d) data;
  vreal     k_real = vset1(REAL(bem->k));
  vreal     k_imag = vset1(-IMAG(bem->k));
  vreal     dist[3];
  vreal     norm, norm2, rnorm, s, c;

  (void) ny;

  dist[0] = vsub(y[0], x[0]);
  dist[1] = vsub(y[1], x[1]);
  dist[2] = vsub(y[2], x[2]);
  norm2 = vdot3(dist, dist);
  rnorm = vrsqrt(norm2);

  norm = vmul(norm2, rnorm);
  rnorm = vmul(vmul(rnorm, vmul(rnorm, rnorm)), vdot3(dist, (vreal *) nx));
  if (IMAG(bem->k) != 0.0) {
    rnorm = vmul(rnorm, vexp(vmul(k_imag, norm)));
  }
  norm = vmul(k_real, norm);

  vsincos(norm, &c, &s);

  *res_re = vmul(vfmadd(norm, s, c), rnorm);
  *res_im = vmul(vfnmadd(norm, c, s), rnorm);

}
#endif

static inline field
hs_kernel_helmholtzbem3d(const real * x, const real * y,
			 const real * nx, const real * ny, void *data)
{
  pcbem3d   bem = (pcbem3d) data;
  real      k_real = REAL(bem->k);
  real      k_imag = -IMAG(bem->k);
  real      dist[3];
  real      norm, norm2, rnorm, s, c, dot, dotxy, hr, hi;

  field     res;

  dist[0] = x[0] - y[0];
  dist[1] = x[1] - y[1];
  dist[2] = x[2] - y[2];
  norm2 = REAL_NORMSQR3(dist[0], dist[1], dist[2]);
  rnorm = REAL_RSQRT(norm2);

  norm = norm2 * rnorm;
  rnorm = (rnorm * rnorm) * (rnorm * rnorm) * rnorm;
  if (k_imag != 0.0) {
    rnorm *= REAL_EXP(k_imag + norm);
  }
  norm = k_real * norm;
  s = REAL_SIN(norm);
  c = REAL_COS(norm);

  dot = REAL_DOT3(dist, nx) * REAL_DOT3(dist, ny);
  dotxy = REAL_DOT3(nx, ny);

  hr = (norm * norm - 3.0) * dot + norm2 * dotxy;
  hi = 3.0 * norm * dot - norm * norm2 * dotxy;

  res = rnorm * (c * hr - s * hi + (c * hi + s * hr) * I);

  return res;
}

#ifdef USE_SIMD
static inline void
hs_kernel_simd_helmholtzbem3d(const vreal * x,
			      const vreal * y, const vreal * nx,
			      const vreal * ny, void *data, vreal * res_re,
			      vreal * res_im)
{
  pcbem3d   bem = (pcbem3d) data;
  vreal     k_real = vset1(REAL(bem->k));
  vreal     k_imag = vset1(-IMAG(bem->k));
  vreal     dist[3];
  vreal     norm, norm2, rnorm, s, c, dot, dotxy, hr, hi, c3;

  c3 = vset1(3.0);

  dist[0] = vsub(x[0], y[0]);
  dist[1] = vsub(x[1], y[1]);
  dist[2] = vsub(x[2], y[2]);
  norm2 = vdot3(dist, dist);
  rnorm = vrsqrt(norm2);

  norm = vmul(norm2, rnorm);
  rnorm = vmul(vmul(rnorm, rnorm), vmul(vmul(rnorm, rnorm), rnorm));
  if (IMAG(bem->k) != 0.0) {
    rnorm = vmul(rnorm, vexp(vmul(k_imag, norm)));
  }
  norm = vmul(k_real, norm);

  vsincos(norm, &c, &s);

  dot = vdot3(dist, (vreal *) nx) * vdot3(dist, (vreal *) ny);
  dotxy = vdot3((vreal *) nx, (vreal *) ny);

  hr = vadd(vmul(vsub(vmul(norm, norm), c3), dot), vmul(norm2, dotxy));
  hi = vsub(vmul(c3, vmul(norm, dot)), vmul(norm, vmul(norm2, dotxy)));

  *res_re = vmul(rnorm, vsub(vmul(c, hr), vmul(s, hi)));
  *res_im = vmul(rnorm, vadd(vmul(c, hi), vmul(s, hr)));
}
#endif

#ifdef USE_SIMD
static void
fill_slp_cc_near_helmholtzbem3d(const uint * ridx,
				const uint * cidx, pcbem3d bem, bool ntrans,
				pamatrix N)
{
  assemble_cc_simd_near_bem3d(ridx, cidx, bem, ntrans, N,
			      (kernel_simd_func3d)
			      slp_kernel_simd_helmholtzbem3d);
}
#else
static void
fill_slp_cc_near_helmholtzbem3d(const uint * ridx, const uint * cidx,
				pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_cc_near_bem3d(ridx, cidx, bem, ntrans, N,
			 (kernel_func3d) slp_kernel_helmholtzbem3d);
}
#endif

#ifdef USE_SIMD
static void
fill_slp_cc_far_helmholtzbem3d(const uint * ridx, const uint * cidx,
			       pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_cc_simd_far_bem3d(ridx, cidx, bem, ntrans, N,
			     (kernel_simd_func3d)
			     slp_kernel_simd_helmholtzbem3d);
}
#else
static void
fill_slp_cc_far_helmholtzbem3d(const uint * ridx, const uint * cidx,
			       pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_cc_far_bem3d(ridx, cidx, bem, ntrans, N,
			(kernel_func3d) slp_kernel_helmholtzbem3d);
}
#endif

#ifdef USE_SIMD
static void
fill_dlp_cc_near_helmholtzbem3d(const uint * ridx,
				const uint * cidx, pcbem3d bem, bool ntrans,
				pamatrix N)
{
  assemble_cc_simd_near_bem3d(ridx, cidx, bem, ntrans, N,
			      (kernel_simd_func3d)
			      dlp_kernel_simd_helmholtzbem3d);
}
#else
static void
fill_dlp_cc_near_helmholtzbem3d(const uint * ridx, const uint * cidx,
				pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_cc_near_bem3d(ridx, cidx, bem, ntrans, N,
			 (kernel_func3d) dlp_kernel_helmholtzbem3d);
}
#endif

#ifdef USE_SIMD
static void
fill_dlp_cc_far_helmholtzbem3d(const uint * ridx, const uint * cidx,
			       pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_cc_simd_far_bem3d(ridx, cidx, bem, ntrans, N,
			     (kernel_simd_func3d)
			     dlp_kernel_simd_helmholtzbem3d);
}
#else
static void
fill_dlp_cc_far_helmholtzbem3d(const uint * ridx, const uint * cidx,
			       pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_cc_far_bem3d(ridx, cidx, bem, ntrans, N,
			(kernel_func3d) dlp_kernel_helmholtzbem3d);
}
#endif

#ifdef USE_SIMD
static void
fill_adlp_cc_near_helmholtzbem3d(const uint * ridx,
				 const uint * cidx, pcbem3d bem, bool ntrans,
				 pamatrix N)
{
  assemble_cc_simd_near_bem3d(ridx, cidx, bem, ntrans, N,
			      (kernel_simd_func3d)
			      adlp_kernel_simd_helmholtzbem3d);
}
#else
static void
fill_adlp_cc_near_helmholtzbem3d(const uint * ridx, const uint * cidx,
				 pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_cc_near_bem3d(ridx, cidx, bem, ntrans, N,
			 (kernel_func3d) adlp_kernel_helmholtzbem3d);
}
#endif

#ifdef USE_SIMD
static void
fill_adlp_cc_far_helmholtzbem3d(const uint * ridx,
				const uint * cidx, pcbem3d bem, bool ntrans,
				pamatrix N)
{
  assemble_cc_simd_far_bem3d(ridx, cidx, bem, ntrans, N,
			     (kernel_simd_func3d)
			     adlp_kernel_simd_helmholtzbem3d);
}
#else
static void
fill_adlp_cc_far_helmholtzbem3d(const uint * ridx, const uint * cidx,
				pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_cc_far_bem3d(ridx, cidx, bem, ntrans, N,
			(kernel_func3d) adlp_kernel_helmholtzbem3d);
}
#endif

#ifdef USE_SIMD
static void
fill_slp_cl_near_helmholtzbem3d(const uint * ridx,
				const uint * cidx, pcbem3d bem, bool ntrans,
				pamatrix N)
{
  assemble_cl_simd_near_bem3d(ridx, cidx, bem, ntrans, N,
			      (kernel_simd_func3d)
			      slp_kernel_simd_helmholtzbem3d);
}
#else
static void
fill_slp_cl_near_helmholtzbem3d(const uint * ridx, const uint * cidx,
				pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_cl_near_bem3d(ridx, cidx, bem, ntrans, N,
			 (kernel_func3d) slp_kernel_helmholtzbem3d);
}
#endif

#ifdef USE_SIMD
static void
fill_slp_cl_far_helmholtzbem3d(const uint * ridx, const uint * cidx,
			       pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_cl_simd_far_bem3d(ridx, cidx, bem, ntrans, N,
			     (kernel_simd_func3d)
			     slp_kernel_simd_helmholtzbem3d);
}
#else
static void
fill_slp_cl_far_helmholtzbem3d(const uint * ridx, const uint * cidx,
			       pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_cl_far_bem3d(ridx, cidx, bem, ntrans, N,
			(kernel_func3d) slp_kernel_helmholtzbem3d);
}
#endif

#ifdef USE_SIMD
static void
fill_dlp_cl_near_helmholtzbem3d(const uint * ridx,
				const uint * cidx, pcbem3d bem, bool ntrans,
				pamatrix N)
{
  assemble_cl_simd_near_bem3d(ridx, cidx, bem, ntrans, N,
			      (kernel_simd_func3d)
			      dlp_kernel_simd_helmholtzbem3d);
}
#else
static void
fill_dlp_cl_near_helmholtzbem3d(const uint * ridx, const uint * cidx,
				pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_cl_near_bem3d(ridx, cidx, bem, ntrans, N,
			 (kernel_func3d) dlp_kernel_helmholtzbem3d);
}
#endif

#ifdef USE_SIMD
static void
fill_dlp_cl_far_helmholtzbem3d(const uint * ridx, const uint * cidx,
			       pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_cl_simd_far_bem3d(ridx, cidx, bem, ntrans, N,
			     (kernel_simd_func3d)
			     dlp_kernel_simd_helmholtzbem3d);
}
#else
static void
fill_dlp_cl_far_helmholtzbem3d(const uint * ridx, const uint * cidx,
			       pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_cl_far_bem3d(ridx, cidx, bem, ntrans, N,
			(kernel_func3d) dlp_kernel_helmholtzbem3d);
}
#endif

#ifdef USE_SIMD
static void
fill_adlp_cl_near_helmholtzbem3d(const uint * ridx,
				 const uint * cidx, pcbem3d bem, bool ntrans,
				 pamatrix N)
{
  assemble_cl_simd_near_bem3d(ridx, cidx, bem, ntrans, N,
			      (kernel_simd_func3d)
			      adlp_kernel_simd_helmholtzbem3d);
}
#else
static void
fill_adlp_cl_near_helmholtzbem3d(const uint * ridx, const uint * cidx,
				 pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_cl_near_bem3d(ridx, cidx, bem, ntrans, N,
			 (kernel_func3d) adlp_kernel_helmholtzbem3d);
}
#endif

#ifdef USE_SIMD
static void
fill_adlp_cl_far_helmholtzbem3d(const uint * ridx,
				const uint * cidx, pcbem3d bem, bool ntrans,
				pamatrix N)
{
  assemble_cl_simd_far_bem3d(ridx, cidx, bem, ntrans, N,
			     (kernel_simd_func3d)
			     adlp_kernel_simd_helmholtzbem3d);
}
#else
static void
fill_adlp_cl_far_helmholtzbem3d(const uint * ridx, const uint * cidx,
				pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_cl_far_bem3d(ridx, cidx, bem, ntrans, N,
			(kernel_func3d) adlp_kernel_helmholtzbem3d);
}
#endif

#ifdef USE_SIMD
static void
fill_slp_lc_near_helmholtzbem3d(const uint * ridx,
				const uint * cidx, pcbem3d bem, bool ntrans,
				pamatrix N)
{
  assemble_lc_simd_near_bem3d(ridx, cidx, bem, ntrans, N,
			      (kernel_simd_func3d)
			      slp_kernel_simd_helmholtzbem3d);
}
#else
static void
fill_slp_lc_near_helmholtzbem3d(const uint * ridx, const uint * cidx,
				pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_lc_near_bem3d(ridx, cidx, bem, ntrans, N,
			 (kernel_func3d) slp_kernel_helmholtzbem3d);
}
#endif

#ifdef USE_SIMD
static void
fill_slp_lc_far_helmholtzbem3d(const uint * ridx, const uint * cidx,
			       pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_lc_simd_far_bem3d(ridx, cidx, bem, ntrans, N,
			     (kernel_simd_func3d)
			     slp_kernel_simd_helmholtzbem3d);
}
#else
static void
fill_slp_lc_far_helmholtzbem3d(const uint * ridx, const uint * cidx,
			       pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_lc_far_bem3d(ridx, cidx, bem, ntrans, N,
			(kernel_func3d) slp_kernel_helmholtzbem3d);
}
#endif

#ifdef USE_SIMD
static void
fill_dlp_lc_near_helmholtzbem3d(const uint * ridx,
				const uint * cidx, pcbem3d bem, bool ntrans,
				pamatrix N)
{
  assemble_lc_simd_near_bem3d(ridx, cidx, bem, ntrans, N,
			      (kernel_simd_func3d)
			      dlp_kernel_simd_helmholtzbem3d);
}
#else
static void
fill_dlp_lc_near_helmholtzbem3d(const uint * ridx, const uint * cidx,
				pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_lc_near_bem3d(ridx, cidx, bem, ntrans, N,
			 (kernel_func3d) dlp_kernel_helmholtzbem3d);
}
#endif

#ifdef USE_SIMD
static void
fill_dlp_lc_far_helmholtzbem3d(const uint * ridx, const uint * cidx,
			       pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_lc_simd_far_bem3d(ridx, cidx, bem, ntrans, N,
			     (kernel_simd_func3d)
			     dlp_kernel_simd_helmholtzbem3d);
}
#else
static void
fill_dlp_lc_far_helmholtzbem3d(const uint * ridx, const uint * cidx,
			       pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_lc_far_bem3d(ridx, cidx, bem, ntrans, N,
			(kernel_func3d) dlp_kernel_helmholtzbem3d);
}
#endif

#ifdef USE_SIMD
static void
fill_adlp_lc_near_helmholtzbem3d(const uint * ridx,
				 const uint * cidx, pcbem3d bem, bool ntrans,
				 pamatrix N)
{
  assemble_lc_simd_near_bem3d(ridx, cidx, bem, ntrans, N,
			      (kernel_simd_func3d)
			      adlp_kernel_simd_helmholtzbem3d);
}
#else
static void
fill_adlp_lc_near_helmholtzbem3d(const uint * ridx, const uint * cidx,
				 pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_lc_near_bem3d(ridx, cidx, bem, ntrans, N,
			 (kernel_func3d) adlp_kernel_helmholtzbem3d);
}
#endif

#ifdef USE_SIMD
static void
fill_adlp_lc_far_helmholtzbem3d(const uint * ridx,
				const uint * cidx, pcbem3d bem, bool ntrans,
				pamatrix N)
{
  assemble_lc_simd_far_bem3d(ridx, cidx, bem, ntrans, N,
			     (kernel_simd_func3d)
			     adlp_kernel_simd_helmholtzbem3d);
}
#else
static void
fill_adlp_lc_far_helmholtzbem3d(const uint * ridx, const uint * cidx,
				pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_lc_far_bem3d(ridx, cidx, bem, ntrans, N,
			(kernel_func3d) adlp_kernel_helmholtzbem3d);
}
#endif

#ifdef USE_SIMD
static void
fill_slp_ll_near_helmholtzbem3d(const uint * ridx,
				const uint * cidx, pcbem3d bem, bool ntrans,
				pamatrix N)
{
  assemble_ll_simd_near_bem3d(ridx, cidx, bem, ntrans, N,
			      (kernel_simd_func3d)
			      slp_kernel_simd_helmholtzbem3d);
}
#else
static void
fill_slp_ll_near_helmholtzbem3d(const uint * ridx, const uint * cidx,
				pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_ll_near_bem3d(ridx, cidx, bem, ntrans, N,
			 (kernel_func3d) slp_kernel_helmholtzbem3d);
}
#endif

#ifdef USE_SIMD
static void
fill_slp_ll_far_helmholtzbem3d(const uint * ridx, const uint * cidx,
			       pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_ll_simd_far_bem3d(ridx, cidx, bem, ntrans, N,
			     (kernel_simd_func3d)
			     slp_kernel_simd_helmholtzbem3d);
}
#else
static void
fill_slp_ll_far_helmholtzbem3d(const uint * ridx, const uint * cidx,
			       pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_ll_far_bem3d(ridx, cidx, bem, ntrans, N,
			(kernel_func3d) slp_kernel_helmholtzbem3d);
}
#endif

#ifdef USE_SIMD
static void
fill_dlp_ll_near_helmholtzbem3d(const uint * ridx,
				const uint * cidx, pcbem3d bem, bool ntrans,
				pamatrix N)
{
  assemble_ll_simd_near_bem3d(ridx, cidx, bem, ntrans, N,
			      (kernel_simd_func3d)
			      dlp_kernel_simd_helmholtzbem3d);
}
#else
static void
fill_dlp_ll_near_helmholtzbem3d(const uint * ridx, const uint * cidx,
				pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_ll_near_bem3d(ridx, cidx, bem, ntrans, N,
			 (kernel_func3d) dlp_kernel_helmholtzbem3d);
}
#endif

#ifdef USE_SIMD
static void
fill_dlp_ll_far_helmholtzbem3d(const uint * ridx, const uint * cidx,
			       pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_ll_simd_far_bem3d(ridx, cidx, bem, ntrans, N,
			     (kernel_simd_func3d)
			     dlp_kernel_simd_helmholtzbem3d);
}
#else
static void
fill_dlp_ll_far_helmholtzbem3d(const uint * ridx, const uint * cidx,
			       pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_ll_far_bem3d(ridx, cidx, bem, ntrans, N,
			(kernel_func3d) dlp_kernel_helmholtzbem3d);
}
#endif

#ifdef USE_SIMD
static void
fill_adlp_ll_near_helmholtzbem3d(const uint * ridx,
				 const uint * cidx, pcbem3d bem, bool ntrans,
				 pamatrix N)
{
  assemble_ll_simd_near_bem3d(ridx, cidx, bem, ntrans, N,
			      (kernel_simd_func3d)
			      adlp_kernel_simd_helmholtzbem3d);
}
#else
static void
fill_adlp_ll_near_helmholtzbem3d(const uint * ridx, const uint * cidx,
				 pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_ll_near_bem3d(ridx, cidx, bem, ntrans, N,
			 (kernel_func3d) adlp_kernel_helmholtzbem3d);
}
#endif

#ifdef USE_SIMD
static void
fill_adlp_ll_far_helmholtzbem3d(const uint * ridx,
				const uint * cidx, pcbem3d bem, bool ntrans,
				pamatrix N)
{
  assemble_ll_simd_far_bem3d(ridx, cidx, bem, ntrans, N,
			     (kernel_simd_func3d)
			     adlp_kernel_simd_helmholtzbem3d);
}
#else
static void
fill_adlp_ll_far_helmholtzbem3d(const uint * ridx, const uint * cidx,
				pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_ll_far_bem3d(ridx, cidx, bem, ntrans, N,
			(kernel_func3d) adlp_kernel_helmholtzbem3d);
}
#endif
static void
fill_kernel_helmholtzbem3d(pcbem3d bem, const real(*X)[3],
			   const real(*Y)[3], pamatrix V)
{
  fill_bem3d(bem, X, Y, NULL, NULL, V, slp_kernel_helmholtzbem3d);
}

static void
fill_kernel_wave_helmholtzbem3d(pcbem3d bem, const real(*X)[3],
				const real(*Y)[3], pcreal dir, pamatrix V)
{
  fill_wave_bem3d(bem, X, Y, NULL, NULL, V, dir,
		  slp_kernel_directional_helmholtzbem3d);
}

static void
fill_dny_kernel_helmholtzbem3d(pcbem3d bem, const real(*X)[3],
			       const real(*Y)[3], const real(*NY)[3],
			       pamatrix V)
{
  fill_bem3d(bem, X, Y, NULL, NY, V, dlp_kernel_helmholtzbem3d);
}

static void
fill_dnx_dny_kernel_helmholtzbem3d(pcbem3d bem, const real(*X)[3],
				   const real(*NX)[3], const real(*Y)[3],
				   const real(*NY)[3], pamatrix V)
{
  fill_bem3d(bem, X, Y, NX, NY, V, hs_kernel_helmholtzbem3d);
}

#ifdef USE_SIMD
static void
fill_kernel_c_helmholtzbem3d(const uint * idx, const real(*Z)[3],
			     pcbem3d bem, pamatrix V)
{
  fill_row_simd_c_bem3d(idx, Z, bem, V, slp_kernel_simd_helmholtzbem3d);
}
#else
static void
fill_kernel_c_helmholtzbem3d(const uint * idx, const real(*Z)[3],
			     pcbem3d bem, pamatrix V)
{
  fill_row_c_bem3d(idx, Z, bem, V, slp_kernel_helmholtzbem3d);
}
#endif

static void
fill_kernel_l_helmholtzbem3d(const uint * idx, const real(*Z)[3],
			     pcbem3d bem, pamatrix V)
{
  fill_row_l_bem3d(idx, Z, bem, V, slp_kernel_helmholtzbem3d);
}

#ifdef USE_SIMD
static void
fill_dnz_kernel_c_helmholtzbem3d(const uint * idx,
				 const real(*Z)[3], const real(*N)[3],
				 pcbem3d bem, pamatrix V)
{
  fill_dnz_row_simd_c_bem3d(idx, Z, N, bem, V,
			    dlp_kernel_simd_helmholtzbem3d);
}
#else
static void
fill_dnz_kernel_c_helmholtzbem3d(const uint * idx,
				 const real(*Z)[3], const real(*N)[3],
				 pcbem3d bem, pamatrix V)
{
  fill_dnz_row_c_bem3d(idx, Z, N, bem, V, dlp_kernel_helmholtzbem3d);
}
#endif

static void
fill_dnz_kernel_l_helmholtzbem3d(const uint * idx,
				 const real(*Z)[3], const real(*N)[3],
				 pcbem3d bem, pamatrix V)
{
  fill_dnz_row_l_bem3d(idx, Z, N, bem, V, dlp_kernel_helmholtzbem3d);
}

#ifdef USE_SIMD
static void
fill_dnzdrow_kernel_c_helmholtzbem3d(const uint * idx,
				     const real(*Z)[3], const real(*N)[3],
				     pcbem3d bem, pamatrix V)
{
  fill_dnz_row_simd_c_bem3d(idx, Z, N, bem, V, hs_kernel_simd_helmholtzbem3d);
}
#else
static void
fill_dnzdrow_kernel_c_helmholtzbem3d(const uint * idx,
				     const real(*Z)[3], const real(*N)[3],
				     pcbem3d bem, pamatrix V)
{
  fill_dnz_row_c_bem3d(idx, Z, N, bem, V, hs_kernel_helmholtzbem3d);
}
#endif

static void
fill_dnzdrow_kernel_l_helmholtzbem3d(const uint * idx,
				     const real(*Z)[3], const real(*N)[3],
				     pcbem3d bem, pamatrix V)
{
  fill_dnz_row_l_bem3d(idx, Z, N, bem, V, hs_kernel_helmholtzbem3d);
}

#ifdef USE_SIMD
static void
fill_dnzdcol_kernel_c_helmholtzbem3d(const uint * idx,
				     const real(*Z)[3], const real(*N)[3],
				     pcbem3d bem, pamatrix V)
{
  fill_dnz_col_simd_c_bem3d(idx, Z, N, bem, V, hs_kernel_simd_helmholtzbem3d);
}
#else
static void
fill_dnzdcol_kernel_c_helmholtzbem3d(const uint * idx,
				     const real(*Z)[3], const real(*N)[3],
				     pcbem3d bem, pamatrix V)
{
  fill_dnz_col_c_bem3d(idx, Z, N, bem, V, hs_kernel_helmholtzbem3d);
}
#endif

static void
fill_dnzdcol_kernel_l_helmholtzbem3d(const uint * idx,
				     const real(*Z)[3], const real(*N)[3],
				     pcbem3d bem, pamatrix V)
{
  fill_dnz_col_l_bem3d(idx, Z, N, bem, V, hs_kernel_helmholtzbem3d);
}

#ifdef USE_SIMD
static void
fill_drow_kernel_c_helmholtzbem3d(const uint * idx,
				  const real(*Z)[3], pcbem3d bem, pamatrix V)
{
  fill_row_simd_c_bem3d(idx, Z, bem, V, adlp_kernel_simd_helmholtzbem3d);
}
#else
static void
fill_drow_kernel_c_helmholtzbem3d(const uint * idx,
				  const real(*Z)[3], pcbem3d bem, pamatrix V)
{
  fill_row_c_bem3d(idx, Z, bem, V, adlp_kernel_helmholtzbem3d);
}
#endif

static void
fill_drow_kernel_l_helmholtzbem3d(const uint * idx,
				  const real(*Z)[3], pcbem3d bem, pamatrix V)
{
  fill_row_l_bem3d(idx, Z, bem, V, adlp_kernel_helmholtzbem3d);
}

#ifdef USE_SIMD
static void
fill_dcol_kernel_c_helmholtzbem3d(const uint * idx,
				  const real(*Z)[3], pcbem3d bem, pamatrix V)
{
  fill_col_simd_c_bem3d(idx, Z, bem, V, dlp_kernel_simd_helmholtzbem3d);
}

#else
static void
fill_dcol_kernel_c_helmholtzbem3d(const uint * idx,
				  const real(*Z)[3], pcbem3d bem, pamatrix V)
{
  fill_col_c_bem3d(idx, Z, bem, V, dlp_kernel_helmholtzbem3d);
}
#endif

static void
fill_dcol_kernel_l_helmholtzbem3d(const uint * idx,
				  const real(*Z)[3], pcbem3d bem, pamatrix V)
{
  fill_col_l_bem3d(idx, Z, bem, V, dlp_kernel_helmholtzbem3d);
}

pbem3d
new_slp_helmholtz_bem3d(field k, pcsurface3d gr, uint q_regular,
			uint q_singular, basisfunctionbem3d row_basis,
			basisfunctionbem3d col_basis)
{
  pkernelbem3d kernels;

  pbem3d    bem;

  bem = new_bem3d(gr, row_basis, col_basis);
  kernels = bem->kernels;

  bem->sq = build_singquad2d(gr, q_regular, q_singular);

  bem->kernel_const = KERNEL_CONST_HELMHOLTZBEM3D;
  bem->k = k;
  bem->alpha = 0.0;

  kernels->fundamental = fill_kernel_helmholtzbem3d;
  kernels->fundamental_wave = fill_kernel_wave_helmholtzbem3d;
  kernels->dny_fundamental = fill_dny_kernel_helmholtzbem3d;
  kernels->dnx_dny_fundamental = fill_dnx_dny_kernel_helmholtzbem3d;

  switch (row_basis) {
  case BASIS_CONSTANT_BEM3D:
    kernels->lagrange_row = assemble_bem3d_lagrange_c_amatrix;
    kernels->lagrange_wave_row = assemble_bem3d_lagrange_wave_c_amatrix;

    kernels->fundamental_row = fill_kernel_c_helmholtzbem3d;
    kernels->dnz_fundamental_row = fill_dnz_kernel_c_helmholtzbem3d;
    kernels->kernel_row = fill_kernel_c_helmholtzbem3d;
    kernels->dnz_kernel_row = fill_dnz_kernel_c_helmholtzbem3d;
    break;
  case BASIS_LINEAR_BEM3D:
    kernels->lagrange_row = assemble_bem3d_lagrange_l_amatrix;
    kernels->lagrange_wave_row = assemble_bem3d_lagrange_wave_l_amatrix;

    kernels->fundamental_row = fill_kernel_l_helmholtzbem3d;
    kernels->dnz_fundamental_row = fill_dnz_kernel_l_helmholtzbem3d;
    kernels->kernel_row = fill_kernel_l_helmholtzbem3d;
    kernels->dnz_kernel_row = fill_dnz_kernel_l_helmholtzbem3d;
    break;
  default:
    fprintf(stderr, "Unknown basis type detected!\n");
    abort();
    break;
  }

  switch (col_basis) {
  case BASIS_CONSTANT_BEM3D:
    kernels->lagrange_col = assemble_bem3d_lagrange_c_amatrix;
    kernels->lagrange_wave_col = assemble_bem3d_lagrange_wave_c_amatrix;

    kernels->fundamental_col = fill_kernel_c_helmholtzbem3d;
    kernels->dnz_fundamental_col = fill_dnz_kernel_c_helmholtzbem3d;
    kernels->kernel_col = fill_kernel_c_helmholtzbem3d;
    kernels->dnz_kernel_col = fill_dnz_kernel_c_helmholtzbem3d;
    break;
  case BASIS_LINEAR_BEM3D:
    kernels->lagrange_col = assemble_bem3d_lagrange_l_amatrix;
    kernels->lagrange_wave_col = assemble_bem3d_lagrange_wave_l_amatrix;

    kernels->fundamental_col = fill_kernel_l_helmholtzbem3d;
    kernels->dnz_fundamental_col = fill_dnz_kernel_l_helmholtzbem3d;
    kernels->kernel_col = fill_kernel_l_helmholtzbem3d;
    kernels->dnz_kernel_col = fill_dnz_kernel_l_helmholtzbem3d;
    break;
  default:
    fprintf(stderr, "Unknown basis type detected!\n");
    abort();
    break;
  }

  switch (row_basis) {
  case BASIS_CONSTANT_BEM3D:
    switch (col_basis) {
    case BASIS_CONSTANT_BEM3D:
      bem->nearfield = fill_slp_cc_near_helmholtzbem3d;
      bem->nearfield_far = fill_slp_cc_far_helmholtzbem3d;
      break;
    case BASIS_LINEAR_BEM3D:
      bem->nearfield = fill_slp_cl_near_helmholtzbem3d;
      bem->nearfield_far = fill_slp_cl_far_helmholtzbem3d;

      weight_basisfunc_cl_singquad2d(bem->sq->x_id, bem->sq->y_id,
				     bem->sq->w_id, bem->sq->n_id);
      weight_basisfunc_cl_singquad2d(bem->sq->x_edge, bem->sq->y_edge,
				     bem->sq->w_edge, bem->sq->n_edge);
      weight_basisfunc_cl_singquad2d(bem->sq->x_vert, bem->sq->y_vert,
				     bem->sq->w_vert, bem->sq->n_vert);
      weight_basisfunc_cl_singquad2d(bem->sq->x_dist, bem->sq->y_dist,
				     bem->sq->w_dist, bem->sq->n_dist);
      weight_basisfunc_l_singquad2d(bem->sq->x_single, bem->sq->y_single,
				    bem->sq->w_single, bem->sq->n_single);

      bem->mass = allocreal(3);
      bem->mass[0] = 1.0 / 6.0;
      bem->mass[1] = 1.0 / 6.0;
      bem->mass[2] = 1.0 / 6.0;
      break;
    default:
      fprintf(stderr, "Unknown basis type detected!\n");
      abort();
      break;
    }
    break;
  case BASIS_LINEAR_BEM3D:
    switch (col_basis) {
    case BASIS_CONSTANT_BEM3D:
      bem->nearfield = fill_slp_lc_near_helmholtzbem3d;
      bem->nearfield_far = fill_slp_lc_far_helmholtzbem3d;

      weight_basisfunc_lc_singquad2d(bem->sq->x_id, bem->sq->y_id,
				     bem->sq->w_id, bem->sq->n_id);
      weight_basisfunc_lc_singquad2d(bem->sq->x_edge, bem->sq->y_edge,
				     bem->sq->w_edge, bem->sq->n_edge);
      weight_basisfunc_lc_singquad2d(bem->sq->x_vert, bem->sq->y_vert,
				     bem->sq->w_vert, bem->sq->n_vert);
      weight_basisfunc_lc_singquad2d(bem->sq->x_dist, bem->sq->y_dist,
				     bem->sq->w_dist, bem->sq->n_dist);
      weight_basisfunc_l_singquad2d(bem->sq->x_single, bem->sq->y_single,
				    bem->sq->w_single, bem->sq->n_single);

      bem->mass = allocreal(3);
      bem->mass[0] = 1.0 / 6.0;
      bem->mass[1] = 1.0 / 6.0;
      bem->mass[2] = 1.0 / 6.0;
      break;
    case BASIS_LINEAR_BEM3D:
      bem->nearfield = fill_slp_ll_near_helmholtzbem3d;
      bem->nearfield_far = fill_slp_ll_far_helmholtzbem3d;

      weight_basisfunc_ll_singquad2d(bem->sq->x_id, bem->sq->y_id,
				     bem->sq->w_id, bem->sq->n_id);
      weight_basisfunc_ll_singquad2d(bem->sq->x_edge, bem->sq->y_edge,
				     bem->sq->w_edge, bem->sq->n_edge);
      weight_basisfunc_ll_singquad2d(bem->sq->x_vert, bem->sq->y_vert,
				     bem->sq->w_vert, bem->sq->n_vert);
      weight_basisfunc_ll_singquad2d(bem->sq->x_dist, bem->sq->y_dist,
				     bem->sq->w_dist, bem->sq->n_dist);
      weight_basisfunc_l_singquad2d(bem->sq->x_single, bem->sq->y_single,
				    bem->sq->w_single, bem->sq->n_single);

      bem->mass = allocreal(9);
      bem->mass[0] = 1.0 / 12.0;
      bem->mass[1] = 1.0 / 24.0;
      bem->mass[2] = 1.0 / 24.0;
      bem->mass[3] = 1.0 / 24.0;
      bem->mass[4] = 1.0 / 12.0;
      bem->mass[5] = 1.0 / 24.0;
      bem->mass[6] = 1.0 / 24.0;
      bem->mass[7] = 1.0 / 24.0;
      bem->mass[8] = 1.0 / 12.0;
      break;
    default:
      fprintf(stderr, "Unknown basis type detected!\n");
      abort();
      break;
    }
    break;
  default:
    fprintf(stderr, "Unknown basis type detected!\n");
    abort();
    break;
  }

  return bem;
}

pbem3d
new_dlp_helmholtz_bem3d(field k, pcsurface3d gr, uint q_regular,
			uint q_singular, basisfunctionbem3d row_basis,
			basisfunctionbem3d col_basis, field alpha)
{
  pkernelbem3d kernels;

  pbem3d    bem;

  bem = new_bem3d(gr, row_basis, col_basis);
  kernels = bem->kernels;

  bem->sq = build_singquad2d(gr, q_regular, q_singular);

  bem->kernel_const = KERNEL_CONST_HELMHOLTZBEM3D;
  bem->k = k;
  bem->alpha = alpha;

  kernels->fundamental = fill_kernel_helmholtzbem3d;
  kernels->fundamental_wave = fill_kernel_wave_helmholtzbem3d;
  kernels->dny_fundamental = fill_dny_kernel_helmholtzbem3d;
  kernels->dnx_dny_fundamental = fill_dnx_dny_kernel_helmholtzbem3d;

  switch (row_basis) {
  case BASIS_CONSTANT_BEM3D:
    kernels->lagrange_row = assemble_bem3d_lagrange_c_amatrix;
    kernels->lagrange_wave_row = assemble_bem3d_lagrange_wave_c_amatrix;

    kernels->fundamental_row = fill_kernel_c_helmholtzbem3d;
    kernels->dnz_fundamental_row = fill_dnz_kernel_c_helmholtzbem3d;
    kernels->kernel_row = fill_kernel_c_helmholtzbem3d;
    kernels->dnz_kernel_row = fill_dnz_kernel_c_helmholtzbem3d;
    break;
  case BASIS_LINEAR_BEM3D:
    kernels->lagrange_row = assemble_bem3d_lagrange_l_amatrix;
    kernels->lagrange_wave_row = assemble_bem3d_lagrange_wave_l_amatrix;

    kernels->fundamental_row = fill_kernel_l_helmholtzbem3d;
    kernels->dnz_fundamental_row = fill_dnz_kernel_l_helmholtzbem3d;
    kernels->kernel_row = fill_kernel_l_helmholtzbem3d;
    kernels->dnz_kernel_row = fill_dnz_kernel_l_helmholtzbem3d;
    break;
  default:
    fprintf(stderr, "Unknown basis type detected!\n");
    abort();
    break;
  }

  switch (col_basis) {
  case BASIS_CONSTANT_BEM3D:
    kernels->lagrange_col = assemble_bem3d_dn_lagrange_c_amatrix;
    kernels->lagrange_wave_col = assemble_bem3d_dn_lagrange_wave_c_amatrix;

    kernels->fundamental_col = fill_kernel_c_helmholtzbem3d;
    kernels->dnz_fundamental_col = fill_dnz_kernel_c_helmholtzbem3d;
    kernels->kernel_col = fill_dcol_kernel_c_helmholtzbem3d;
    kernels->dnz_kernel_col = fill_dnzdcol_kernel_c_helmholtzbem3d;
    break;
  case BASIS_LINEAR_BEM3D:
    kernels->lagrange_col = assemble_bem3d_dn_lagrange_l_amatrix;
    kernels->lagrange_wave_col = assemble_bem3d_dn_lagrange_wave_l_amatrix;

    kernels->fundamental_col = fill_kernel_l_helmholtzbem3d;
    kernels->dnz_fundamental_col = fill_dnz_kernel_l_helmholtzbem3d;
    kernels->kernel_col = fill_dcol_kernel_l_helmholtzbem3d;
    kernels->dnz_kernel_col = fill_dnzdcol_kernel_l_helmholtzbem3d;
    break;
  default:
    fprintf(stderr, "Unknown basis type detected!\n");
    abort();
    break;
  }

  switch (row_basis) {
  case BASIS_CONSTANT_BEM3D:
    switch (col_basis) {
    case BASIS_CONSTANT_BEM3D:
      bem->nearfield = fill_dlp_cc_near_helmholtzbem3d;
      bem->nearfield_far = fill_dlp_cc_far_helmholtzbem3d;
      break;
    case BASIS_LINEAR_BEM3D:
      bem->nearfield = fill_dlp_cl_near_helmholtzbem3d;
      bem->nearfield_far = fill_dlp_cl_far_helmholtzbem3d;

      weight_basisfunc_cl_singquad2d(bem->sq->x_id, bem->sq->y_id,
				     bem->sq->w_id, bem->sq->n_id);
      weight_basisfunc_cl_singquad2d(bem->sq->x_edge, bem->sq->y_edge,
				     bem->sq->w_edge, bem->sq->n_edge);
      weight_basisfunc_cl_singquad2d(bem->sq->x_vert, bem->sq->y_vert,
				     bem->sq->w_vert, bem->sq->n_vert);
      weight_basisfunc_cl_singquad2d(bem->sq->x_dist, bem->sq->y_dist,
				     bem->sq->w_dist, bem->sq->n_dist);
      weight_basisfunc_l_singquad2d(bem->sq->x_single, bem->sq->y_single,
				    bem->sq->w_single, bem->sq->n_single);

      bem->mass = allocreal(3);
      bem->mass[0] = 1.0 / 6.0;
      bem->mass[1] = 1.0 / 6.0;
      bem->mass[2] = 1.0 / 6.0;
      break;
    default:
      fprintf(stderr, "Unknown basis type detected!\n");
      abort();
      break;
    }
    break;
  case BASIS_LINEAR_BEM3D:
    switch (col_basis) {
    case BASIS_CONSTANT_BEM3D:
      bem->nearfield = fill_dlp_lc_near_helmholtzbem3d;
      bem->nearfield_far = fill_dlp_lc_far_helmholtzbem3d;

      weight_basisfunc_lc_singquad2d(bem->sq->x_id, bem->sq->y_id,
				     bem->sq->w_id, bem->sq->n_id);
      weight_basisfunc_lc_singquad2d(bem->sq->x_edge, bem->sq->y_edge,
				     bem->sq->w_edge, bem->sq->n_edge);
      weight_basisfunc_lc_singquad2d(bem->sq->x_vert, bem->sq->y_vert,
				     bem->sq->w_vert, bem->sq->n_vert);
      weight_basisfunc_lc_singquad2d(bem->sq->x_dist, bem->sq->y_dist,
				     bem->sq->w_dist, bem->sq->n_dist);
      weight_basisfunc_l_singquad2d(bem->sq->x_single, bem->sq->y_single,
				    bem->sq->w_single, bem->sq->n_single);

      bem->mass = allocreal(3);
      bem->mass[0] = 1.0 / 6.0;
      bem->mass[1] = 1.0 / 6.0;
      bem->mass[2] = 1.0 / 6.0;
      break;
    case BASIS_LINEAR_BEM3D:
      bem->nearfield = fill_dlp_ll_near_helmholtzbem3d;
      bem->nearfield_far = fill_dlp_ll_far_helmholtzbem3d;

      weight_basisfunc_ll_singquad2d(bem->sq->x_id, bem->sq->y_id,
				     bem->sq->w_id, bem->sq->n_id);
      weight_basisfunc_ll_singquad2d(bem->sq->x_edge, bem->sq->y_edge,
				     bem->sq->w_edge, bem->sq->n_edge);
      weight_basisfunc_ll_singquad2d(bem->sq->x_vert, bem->sq->y_vert,
				     bem->sq->w_vert, bem->sq->n_vert);
      weight_basisfunc_ll_singquad2d(bem->sq->x_dist, bem->sq->y_dist,
				     bem->sq->w_dist, bem->sq->n_dist);
      weight_basisfunc_l_singquad2d(bem->sq->x_single, bem->sq->y_single,
				    bem->sq->w_single, bem->sq->n_single);

      bem->mass = allocreal(9);
      bem->mass[0] = 1.0 / 12.0;
      bem->mass[1] = 1.0 / 24.0;
      bem->mass[2] = 1.0 / 24.0;
      bem->mass[3] = 1.0 / 24.0;
      bem->mass[4] = 1.0 / 12.0;
      bem->mass[5] = 1.0 / 24.0;
      bem->mass[6] = 1.0 / 24.0;
      bem->mass[7] = 1.0 / 24.0;
      bem->mass[8] = 1.0 / 12.0;
      break;
    default:
      fprintf(stderr, "Unknown basis type detected!\n");
      abort();
      break;
    }
    break;
  default:
    fprintf(stderr, "Unknown basis type detected!\n");
    abort();
    break;
  }

  return bem;
}

pbem3d
new_adlp_helmholtz_bem3d(field k, pcsurface3d gr, uint q_regular,
			 uint q_singular, basisfunctionbem3d row_basis,
			 basisfunctionbem3d col_basis, field alpha)
{
  pkernelbem3d kernels;

  pbem3d    bem;

  bem = new_bem3d(gr, row_basis, col_basis);
  kernels = bem->kernels;

  bem->sq = build_singquad2d(gr, q_regular, q_singular);

  bem->kernel_const = KERNEL_CONST_HELMHOLTZBEM3D;
  bem->k = k;
  bem->alpha = alpha;

  kernels->fundamental = fill_kernel_helmholtzbem3d;
  kernels->fundamental_wave = fill_kernel_wave_helmholtzbem3d;
  kernels->dny_fundamental = fill_dny_kernel_helmholtzbem3d;
  kernels->dnx_dny_fundamental = fill_dnx_dny_kernel_helmholtzbem3d;

  switch (row_basis) {
  case BASIS_CONSTANT_BEM3D:
    kernels->lagrange_row = assemble_bem3d_dn_lagrange_c_amatrix;
    kernels->lagrange_wave_row = assemble_bem3d_dn_lagrange_wave_c_amatrix;

    kernels->fundamental_row = fill_kernel_c_helmholtzbem3d;
    kernels->dnz_fundamental_row = fill_dnz_kernel_c_helmholtzbem3d;
    kernels->kernel_row = fill_drow_kernel_c_helmholtzbem3d;
    kernels->dnz_kernel_row = fill_dnzdrow_kernel_c_helmholtzbem3d;
    break;
  case BASIS_LINEAR_BEM3D:
    kernels->lagrange_row = assemble_bem3d_dn_lagrange_l_amatrix;
    kernels->lagrange_wave_row = assemble_bem3d_dn_lagrange_wave_l_amatrix;

    kernels->fundamental_row = fill_kernel_l_helmholtzbem3d;
    kernels->dnz_fundamental_row = fill_dnz_kernel_l_helmholtzbem3d;
    kernels->kernel_row = fill_drow_kernel_l_helmholtzbem3d;
    kernels->dnz_kernel_row = fill_dnzdrow_kernel_l_helmholtzbem3d;
    break;
  default:
    fprintf(stderr, "Unknown basis type detected!\n");
    abort();
    break;
  }

  switch (col_basis) {
  case BASIS_CONSTANT_BEM3D:
    kernels->lagrange_col = assemble_bem3d_lagrange_c_amatrix;
    kernels->lagrange_wave_col = assemble_bem3d_lagrange_wave_c_amatrix;

    kernels->fundamental_col = fill_kernel_c_helmholtzbem3d;
    kernels->dnz_fundamental_col = fill_dnz_kernel_c_helmholtzbem3d;
    kernels->kernel_col = fill_kernel_c_helmholtzbem3d;
    kernels->dnz_kernel_col = fill_dnz_kernel_c_helmholtzbem3d;
    break;
  case BASIS_LINEAR_BEM3D:
    kernels->lagrange_col = assemble_bem3d_lagrange_l_amatrix;
    kernels->lagrange_wave_col = assemble_bem3d_lagrange_wave_l_amatrix;

    kernels->fundamental_col = fill_kernel_l_helmholtzbem3d;
    kernels->dnz_fundamental_col = fill_dnz_kernel_l_helmholtzbem3d;
    kernels->kernel_col = fill_kernel_l_helmholtzbem3d;
    kernels->dnz_kernel_col = fill_dnz_kernel_l_helmholtzbem3d;
    break;
  default:
    fprintf(stderr, "Unknown basis type detected!\n");
    abort();
    break;
  }

  switch (row_basis) {
  case BASIS_CONSTANT_BEM3D:
    switch (col_basis) {
    case BASIS_CONSTANT_BEM3D:
      bem->nearfield = fill_adlp_cc_near_helmholtzbem3d;
      bem->nearfield_far = fill_adlp_cc_far_helmholtzbem3d;
      break;
    case BASIS_LINEAR_BEM3D:
      bem->nearfield = fill_adlp_cl_near_helmholtzbem3d;
      bem->nearfield_far = fill_adlp_cl_far_helmholtzbem3d;

      weight_basisfunc_cl_singquad2d(bem->sq->x_id, bem->sq->y_id,
				     bem->sq->w_id, bem->sq->n_id);
      weight_basisfunc_cl_singquad2d(bem->sq->x_edge, bem->sq->y_edge,
				     bem->sq->w_edge, bem->sq->n_edge);
      weight_basisfunc_cl_singquad2d(bem->sq->x_vert, bem->sq->y_vert,
				     bem->sq->w_vert, bem->sq->n_vert);
      weight_basisfunc_cl_singquad2d(bem->sq->x_dist, bem->sq->y_dist,
				     bem->sq->w_dist, bem->sq->n_dist);
      weight_basisfunc_l_singquad2d(bem->sq->x_single, bem->sq->y_single,
				    bem->sq->w_single, bem->sq->n_single);

      bem->mass = allocreal(3);
      bem->mass[0] = 1.0 / 6.0;
      bem->mass[1] = 1.0 / 6.0;
      bem->mass[2] = 1.0 / 6.0;
      break;
    default:
      fprintf(stderr, "Unknown basis type detected!\n");
      abort();
      break;
    }
    break;
  case BASIS_LINEAR_BEM3D:
    switch (col_basis) {
    case BASIS_CONSTANT_BEM3D:
      bem->nearfield = fill_adlp_lc_near_helmholtzbem3d;
      bem->nearfield_far = fill_adlp_lc_far_helmholtzbem3d;

      weight_basisfunc_lc_singquad2d(bem->sq->x_id, bem->sq->y_id,
				     bem->sq->w_id, bem->sq->n_id);
      weight_basisfunc_lc_singquad2d(bem->sq->x_edge, bem->sq->y_edge,
				     bem->sq->w_edge, bem->sq->n_edge);
      weight_basisfunc_lc_singquad2d(bem->sq->x_vert, bem->sq->y_vert,
				     bem->sq->w_vert, bem->sq->n_vert);
      weight_basisfunc_lc_singquad2d(bem->sq->x_dist, bem->sq->y_dist,
				     bem->sq->w_dist, bem->sq->n_dist);
      weight_basisfunc_l_singquad2d(bem->sq->x_single, bem->sq->y_single,
				    bem->sq->w_single, bem->sq->n_single);

      bem->mass = allocreal(3);
      bem->mass[0] = 1.0 / 6.0;
      bem->mass[1] = 1.0 / 6.0;
      bem->mass[2] = 1.0 / 6.0;
      break;
    case BASIS_LINEAR_BEM3D:
      bem->nearfield = fill_adlp_ll_near_helmholtzbem3d;
      bem->nearfield_far = fill_adlp_ll_far_helmholtzbem3d;

      weight_basisfunc_ll_singquad2d(bem->sq->x_id, bem->sq->y_id,
				     bem->sq->w_id, bem->sq->n_id);
      weight_basisfunc_ll_singquad2d(bem->sq->x_edge, bem->sq->y_edge,
				     bem->sq->w_edge, bem->sq->n_edge);
      weight_basisfunc_ll_singquad2d(bem->sq->x_vert, bem->sq->y_vert,
				     bem->sq->w_vert, bem->sq->n_vert);
      weight_basisfunc_ll_singquad2d(bem->sq->x_dist, bem->sq->y_dist,
				     bem->sq->w_dist, bem->sq->n_dist);
      weight_basisfunc_l_singquad2d(bem->sq->x_single, bem->sq->y_single,
				    bem->sq->w_single, bem->sq->n_single);

      bem->mass = allocreal(9);
      bem->mass[0] = 1.0 / 12.0;
      bem->mass[1] = 1.0 / 24.0;
      bem->mass[2] = 1.0 / 24.0;
      bem->mass[3] = 1.0 / 24.0;
      bem->mass[4] = 1.0 / 12.0;
      bem->mass[5] = 1.0 / 24.0;
      bem->mass[6] = 1.0 / 24.0;
      bem->mass[7] = 1.0 / 24.0;
      bem->mass[8] = 1.0 / 12.0;
      break;
    default:
      fprintf(stderr, "Unknown basis type detected!\n");
      abort();
      break;
    }
    break;
  default:
    fprintf(stderr, "Unknown basis type detected!\n");
    abort();
    break;
  }

  return bem;
}

void
del_helmholtz_bem3d(pbem3d bem)
{
  del_bem3d(bem);
}

field
rhs_dirichlet_point_helmholtzbem3d(const real * x, const real * n,
				   const void *data)
{
  helmholtz_data *hdata = (helmholtz_data *) data;
  real     *kvec = hdata->kvec;
  field     k =
    REAL_SQRT(REAL_SQR(kvec[0]) + REAL_SQR(kvec[1]) + REAL_SQR(kvec[2]));
  real     *s = hdata->source;

  real      r[3];
  real      norm, norm2, rnorm;

  field     res;

  (void) n;

  r[0] = x[0] - s[0];
  r[1] = x[1] - s[1];
  r[2] = x[2] - s[2];

  norm2 = REAL_DOT3(r, r);
  rnorm = REAL_RSQRT(norm2);
  norm = norm2 * rnorm;

  res = cexp(I * k * norm) * rnorm;

  return res;
}

field
rhs_neumann_point_helmholtzbem3d(const real * x, const real * n,
				 const void *data)
{
  helmholtz_data *hdata = (helmholtz_data *) data;
  real     *kvec = hdata->kvec;
  field     k =
    REAL_SQRT(REAL_SQR(kvec[0]) + REAL_SQR(kvec[1]) + REAL_SQR(kvec[2]));
  real     *s = hdata->source;

  real      r[3];
  real      norm, norm2, rnorm;

  field     res;

  r[0] = x[0] - s[0];
  r[1] = x[1] - s[1];
  r[2] = x[2] - s[2];

  norm2 = REAL_DOT3(r, r);
  rnorm = REAL_RSQRT(norm2);
  norm = norm2 * rnorm;

  res = cexp(I * k * norm) * rnorm;
  res = res * rnorm * rnorm * (I * k * norm - 1.0) * REAL_DOT3(n, r);

  return res;
}

field
rhs_robin_point_helmholtzbem3d(const real * x, const real * n,
			       const void *data)
{
  helmholtz_data *hdata = (helmholtz_data *) data;
  real     *kvec = hdata->kvec;
  field     k =
    REAL_SQRT(REAL_SQR(kvec[0]) + REAL_SQR(kvec[1]) + REAL_SQR(kvec[2]));
  real     *s = hdata->source;

  real      r[3];
  real      norm, norm2, rnorm;

  field     res;

  r[0] = x[0] - s[0];
  r[1] = x[1] - s[1];
  r[2] = x[2] - s[2];

  norm2 = REAL_DOT3(r, r);
  rnorm = REAL_RSQRT(norm2);
  norm = norm2 * rnorm;

  res = cexp(I * k * norm) * rnorm;

  res *= (rnorm * rnorm * (I * k * norm - 1.0) * REAL_DOT3(n, r)
	  - I * hdata->eta);

  return res;
}

field
rhs_dirichlet_plane_helmholtzbem3d(const real * x, const real * n,
				   const void *data)
{
  helmholtz_data *hdata = (helmholtz_data *) data;
  real     *kvec = hdata->kvec;
  real     *s = hdata->source;

  real      r[3];

  field     res;

  (void) n;

  r[0] = x[0] - s[0];
  r[1] = x[1] - s[1];
  r[2] = x[2] - s[2];

  res = cexp(I * REAL_DOT3(kvec, r));

  return res;
}

field
rhs_neumann_plane_helmholtzbem3d(const real * x, const real * n,
				 const void *data)
{
  helmholtz_data *hdata = (helmholtz_data *) data;
  real     *kvec = hdata->kvec;
  real     *s = hdata->source;

  real      r[3];

  field     res;

  r[0] = x[0] - s[0];
  r[1] = x[1] - s[1];
  r[2] = x[2] - s[2];

  res = I * REAL_DOT3(kvec, n) * cexp(I * REAL_DOT3(kvec, r));

  return res;
}

field
rhs_robin_plane_helmholtzbem3d(const real * x, const real * n,
			       const void *data)
{
  helmholtz_data *hdata = (helmholtz_data *) data;
  real     *kvec = hdata->kvec;
  real     *s = hdata->source;

  real      r[3];

  field     res;

  r[0] = x[0] - s[0];
  r[1] = x[1] - s[1];
  r[2] = x[2] - s[2];

  res = I * (REAL_DOT3(kvec, n) - hdata->eta) * cexp(I * REAL_DOT3(kvec, r));

  return res;
}

#endif
