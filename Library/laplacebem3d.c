/* ------------------------------------------------------------
 This is the file "laplacebem3d.c" of the H2Lib package.
 All rights reserved, Sven Christophersen 2011
 ------------------------------------------------------------ */

/**
 * @file laplacebem3d.c
 * @author Sven Christophersen
 * @date 2011
 */

/* C STD LIBRARY */
#include <string.h>
/* CORE 0 */
#include "basic.h"
#include "parameters.h"
#ifdef USE_OPENCL
#include "opencl.h"
#endif
/* CORE 1 */
/* CORE 2 */
/* CORE 3 */
/* SIMPLE */
/* PARTICLES */
/* BEM */
#include "laplacebem3d.h"
#include "singquad2d.h"

/* This constant comes from the fundamental solution of the Laplace-equation,
 * which is defined as: @f$ g(x,y) := \frac{1}{4 \pi} \frac{1}{ \left\lVert x - y \right\rVert}  @f$ .
 * Therefore the constant takes the value @f$ \frac{1}{4 \pi}@f$ .
 * */
#define KERNEL_CONST_LAPLACEBEM3D 0.0795774715459476679

static inline field
slp_kernel_laplacebem3d(const real * x, const real * y,
			const real * nx, const real * ny, void *data)
{
  real      dist[3];
  real      norm;

  field     res;

  (void) nx;
  (void) ny;
  (void) data;

  dist[0] = x[0] - y[0];
  dist[1] = x[1] - y[1];
  dist[2] = x[2] - y[2];
  norm = REAL_NORMSQR3(dist[0], dist[1], dist[2]);

  res = REAL_RSQRT(norm);

  return res;
}

#ifdef USE_SIMD
static inline void
slp_kernel_simd_laplacebem3d(const vreal * x,
			     const vreal * y, const vreal * nx,
			     const vreal * ny, void *data, vreal * res_re,
			     vreal * res_im)
{
  vreal     dist[3];
  vreal     rnorm;

  (void) nx;
  (void) ny;
  (void) data;

  dist[0] = vsub(x[0], y[0]);
  dist[1] = vsub(x[1], y[1]);
  dist[2] = vsub(x[2], y[2]);

  rnorm = vrsqrt(vdot3(dist, dist));

  *res_re = rnorm;
  *res_im = vset1(0.0);
}
#endif

static inline field
dlp_kernel_laplacebem3d(const real * x, const real * y,
			const real * nx, const real * ny, void *data)
{
  real      dist[3];
  real      norm;

  field     res;

  (void) nx;
  (void) data;

  dist[0] = x[0] - y[0];
  dist[1] = x[1] - y[1];
  dist[2] = x[2] - y[2];
  norm = REAL_NORMSQR3(dist[0], dist[1], dist[2]);

  norm = REAL_RSQRT(norm);
  norm *= norm * norm;

  res = REAL_DOT3(dist, ny) * norm;

  return res;
}

#ifdef USE_SIMD
static inline void
dlp_kernel_simd_laplacebem3d(const vreal * x,
			     const vreal * y, const vreal * nx,
			     const vreal * ny, void *data, vreal * res_re,
			     vreal * res_im)
{
  vreal     dist[3];
  vreal     rnorm, res;

  (void) nx;
  (void) ny;
  (void) data;

  dist[0] = vsub(x[0], y[0]);
  dist[1] = vsub(x[1], y[1]);
  dist[2] = vsub(x[2], y[2]);

  rnorm = vrsqrt(vdot3(dist, dist));
  rnorm = vmul(vmul(rnorm, rnorm), rnorm);
  res = vmul(vdot3(dist, (vreal *) ny), rnorm);

  *res_re = res;
  *res_im = vset1(0.0);
}
#endif

static inline field
adlp_kernel_laplacebem3d(const real * x, const real * y,
			 const real * nx, const real * ny, void *data)
{
  real      dist[3];
  real      norm;

  field     res;

  (void) ny;
  (void) data;

  dist[0] = y[0] - x[0];
  dist[1] = y[1] - x[1];
  dist[2] = y[2] - x[2];
  norm = REAL_NORMSQR3(dist[0], dist[1], dist[2]);

  norm = REAL_RSQRT(norm);
  norm *= norm * norm;

  res = REAL_DOT3(dist, nx) * norm;

  return res;
}

#ifdef USE_SIMD
static inline void
adlp_kernel_simd_laplacebem3d(const vreal * x,
			      const vreal * y, const vreal * nx,
			      const vreal * ny, void *data, vreal * res_re,
			      vreal * res_im)
{
  vreal     dist[3];
  vreal     rnorm, res;

  (void) ny;
  (void) data;

  dist[0] = vsub(y[0], x[0]);
  dist[1] = vsub(y[1], x[1]);
  dist[2] = vsub(y[2], x[2]);

  rnorm = vrsqrt(vdot3(dist, dist));
  rnorm = vmul(vmul(rnorm, rnorm), rnorm);
  res = vmul(vdot3(dist, (vreal *) nx), rnorm);

  *res_re = res;
  *res_im = vset1(0.0);
}
#endif

static inline field
hs_kernel_laplacebem3d(const real * x, const real * y,
		       const real * nx, const real * ny, void *data)
{
  real      dist[3];
  real      norm, norm2, dot, dotxy;

  field     res;

  (void) data;

  dist[0] = x[0] - y[0];
  dist[1] = x[1] - y[1];
  dist[2] = x[2] - y[2];
  norm2 = REAL_NORMSQR3(dist[0], dist[1], dist[2]);
  norm = REAL_RSQRT(norm2);
  norm = (norm * norm) * (norm * norm) * norm;

  dot = REAL_DOT3(dist, nx) * REAL_DOT3(dist, ny);
  dotxy = REAL_DOT3(nx, ny);

  res = norm * (-3.0 * dot + norm2 * dotxy);

  return res;
}

#ifdef USE_SIMD
static inline void
hs_kernel_simd_laplacebem3d(const vreal * x, const vreal * y,
			    const vreal * nx, const vreal * ny, void *data,
			    vreal * res_re, vreal * res_im)
{
  vreal     dist[3];
  vreal     norm2, rnorm, res, dot, c3;

  (void) data;

  c3 = vset1(-3.0);

  dist[0] = vsub(x[0], y[0]);
  dist[1] = vsub(x[1], y[1]);
  dist[2] = vsub(x[2], y[2]);

  norm2 = vdot3(dist, dist);
  rnorm = vrsqrt(norm2);
  rnorm = vmul(vmul(rnorm, rnorm), vmul(vmul(rnorm, rnorm), rnorm));
  dot = vmul(vdot3(dist, (vreal *) nx), vdot3(dist, (vreal *) ny));

  res = vmul(norm2, vdot3((vreal *) nx, (vreal *) ny));
  res = vmul(rnorm, vadd(vmul(c3, dot), res));

  *res_re = res;
  *res_im = vset1(0.0);
}
#endif

#ifdef USE_SIMD
static void
fill_slp_cc_near_laplacebem3d(const uint * ridx, const uint * cidx,
			      pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_cc_simd_near_bem3d(ridx, cidx, bem, ntrans, N, (kernel_simd_func3d)
			      slp_kernel_simd_laplacebem3d);
}
#else
static void
fill_slp_cc_near_laplacebem3d(const uint * ridx, const uint * cidx,
			      pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_cc_near_bem3d(ridx, cidx, bem, ntrans, N,
			 (kernel_func3d) slp_kernel_laplacebem3d);
}
#endif

#ifdef USE_SIMD
static void
fill_slp_cc_far_laplacebem3d(const uint * ridx, const uint * cidx,
			     pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_cc_simd_far_bem3d(ridx, cidx, bem, ntrans, N, (kernel_simd_func3d)
			     slp_kernel_simd_laplacebem3d);
}
#else
static void
fill_slp_cc_far_laplacebem3d(const uint * ridx, const uint * cidx,
			     pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_cc_far_bem3d(ridx, cidx, bem, ntrans, N,
			(kernel_func3d) slp_kernel_laplacebem3d);
}
#endif

#ifdef USE_SIMD
static void
fill_dlp_cc_near_laplacebem3d(const uint * ridx, const uint * cidx,
			      pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_cc_simd_near_bem3d(ridx, cidx, bem, ntrans, N, (kernel_simd_func3d)
			      dlp_kernel_simd_laplacebem3d);
}
#else
static void
fill_dlp_cc_near_laplacebem3d(const uint * ridx, const uint * cidx,
			      pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_cc_near_bem3d(ridx, cidx, bem, ntrans, N,
			 (kernel_func3d) dlp_kernel_laplacebem3d);
}
#endif

#ifdef USE_SIMD
static void
fill_dlp_cc_far_laplacebem3d(const uint * ridx, const uint * cidx,
			     pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_cc_simd_far_bem3d(ridx, cidx, bem, ntrans, N, (kernel_simd_func3d)
			     dlp_kernel_simd_laplacebem3d);
}
#else
static void
fill_dlp_cc_far_laplacebem3d(const uint * ridx, const uint * cidx,
			     pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_cc_far_bem3d(ridx, cidx, bem, ntrans, N,
			(kernel_func3d) dlp_kernel_laplacebem3d);
}
#endif

#ifdef USE_SIMD
static void
fill_adlp_cc_near_laplacebem3d(const uint * ridx, const uint * cidx,
			       pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_cc_simd_near_bem3d(ridx, cidx, bem, ntrans, N, (kernel_simd_func3d)
			      adlp_kernel_simd_laplacebem3d);
}
#else
static void
fill_adlp_cc_near_laplacebem3d(const uint * ridx, const uint * cidx,
			       pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_cc_near_bem3d(ridx, cidx, bem, ntrans, N,
			 (kernel_func3d) adlp_kernel_laplacebem3d);
}
#endif

#ifdef USE_SIMD
static void
fill_adlp_cc_far_laplacebem3d(const uint * ridx, const uint * cidx,
			      pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_cc_simd_far_bem3d(ridx, cidx, bem, ntrans, N, (kernel_simd_func3d)
			     adlp_kernel_simd_laplacebem3d);
}
#else
static void
fill_adlp_cc_far_laplacebem3d(const uint * ridx, const uint * cidx,
			      pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_cc_far_bem3d(ridx, cidx, bem, ntrans, N,
			(kernel_func3d) adlp_kernel_laplacebem3d);
}
#endif

#ifdef USE_SIMD
static void
fill_slp_cl_near_laplacebem3d(const uint * ridx, const uint * cidx,
			      pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_cl_simd_near_bem3d(ridx, cidx, bem, ntrans, N, (kernel_simd_func3d)
			      slp_kernel_simd_laplacebem3d);
}
#else
static void
fill_slp_cl_near_laplacebem3d(const uint * ridx, const uint * cidx,
			      pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_cl_near_bem3d(ridx, cidx, bem, ntrans, N,
			 (kernel_func3d) slp_kernel_laplacebem3d);
}
#endif

#ifdef USE_SIMD
static void
fill_slp_cl_far_laplacebem3d(const uint * ridx, const uint * cidx,
			     pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_cl_simd_far_bem3d(ridx, cidx, bem, ntrans, N, (kernel_simd_func3d)
			     slp_kernel_simd_laplacebem3d);
}
#else
static void
fill_slp_cl_far_laplacebem3d(const uint * ridx, const uint * cidx,
			     pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_cl_far_bem3d(ridx, cidx, bem, ntrans, N,
			(kernel_func3d) slp_kernel_laplacebem3d);
}
#endif

#ifdef USE_SIMD
static void
fill_dlp_cl_near_laplacebem3d(const uint * ridx, const uint * cidx,
			      pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_cl_simd_near_bem3d(ridx, cidx, bem, ntrans, N, (kernel_simd_func3d)
			      dlp_kernel_simd_laplacebem3d);
}
#else
static void
fill_dlp_cl_near_laplacebem3d(const uint * ridx, const uint * cidx,
			      pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_cl_near_bem3d(ridx, cidx, bem, ntrans, N,
			 (kernel_func3d) dlp_kernel_laplacebem3d);
}
#endif

#ifdef USE_SIMD
static void
fill_dlp_cl_far_laplacebem3d(const uint * ridx, const uint * cidx,
			     pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_cl_simd_far_bem3d(ridx, cidx, bem, ntrans, N, (kernel_simd_func3d)
			     dlp_kernel_simd_laplacebem3d);
}
#else
static void
fill_dlp_cl_far_laplacebem3d(const uint * ridx, const uint * cidx,
			     pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_cl_far_bem3d(ridx, cidx, bem, ntrans, N,
			(kernel_func3d) dlp_kernel_laplacebem3d);
}
#endif

#ifdef USE_SIMD
static void
fill_adlp_cl_near_laplacebem3d(const uint * ridx, const uint * cidx,
			       pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_cl_simd_near_bem3d(ridx, cidx, bem, ntrans, N, (kernel_simd_func3d)
			      adlp_kernel_simd_laplacebem3d);
}
#else
static void
fill_adlp_cl_near_laplacebem3d(const uint * ridx, const uint * cidx,
			       pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_cl_near_bem3d(ridx, cidx, bem, ntrans, N,
			 (kernel_func3d) adlp_kernel_laplacebem3d);
}
#endif

#ifdef USE_SIMD
static void
fill_adlp_cl_far_laplacebem3d(const uint * ridx, const uint * cidx,
			      pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_cl_simd_far_bem3d(ridx, cidx, bem, ntrans, N, (kernel_simd_func3d)
			     adlp_kernel_simd_laplacebem3d);
}
#else
static void
fill_adlp_cl_far_laplacebem3d(const uint * ridx, const uint * cidx,
			      pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_cl_far_bem3d(ridx, cidx, bem, ntrans, N,
			(kernel_func3d) adlp_kernel_laplacebem3d);
}
#endif

#ifdef USE_SIMD
static void
fill_slp_lc_near_laplacebem3d(const uint * ridx, const uint * cidx,
			      pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_lc_simd_near_bem3d(ridx, cidx, bem, ntrans, N, (kernel_simd_func3d)
			      slp_kernel_simd_laplacebem3d);
}
#else
static void
fill_slp_lc_near_laplacebem3d(const uint * ridx, const uint * cidx,
			      pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_lc_near_bem3d(ridx, cidx, bem, ntrans, N,
			 (kernel_func3d) slp_kernel_laplacebem3d);
}
#endif

#ifdef USE_SIMD
static void
fill_slp_lc_far_laplacebem3d(const uint * ridx, const uint * cidx,
			     pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_lc_simd_far_bem3d(ridx, cidx, bem, ntrans, N, (kernel_simd_func3d)
			     slp_kernel_simd_laplacebem3d);
}
#else
static void
fill_slp_lc_far_laplacebem3d(const uint * ridx, const uint * cidx,
			     pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_lc_far_bem3d(ridx, cidx, bem, ntrans, N,
			(kernel_func3d) slp_kernel_laplacebem3d);
}
#endif

#ifdef USE_SIMD
static void
fill_dlp_lc_near_laplacebem3d(const uint * ridx, const uint * cidx,
			      pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_lc_simd_near_bem3d(ridx, cidx, bem, ntrans, N, (kernel_simd_func3d)
			      dlp_kernel_simd_laplacebem3d);
}
#else
static void
fill_dlp_lc_near_laplacebem3d(const uint * ridx, const uint * cidx,
			      pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_lc_near_bem3d(ridx, cidx, bem, ntrans, N,
			 (kernel_func3d) dlp_kernel_laplacebem3d);
}
#endif

#ifdef USE_SIMD
static void
fill_dlp_lc_far_laplacebem3d(const uint * ridx, const uint * cidx,
			     pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_lc_simd_far_bem3d(ridx, cidx, bem, ntrans, N, (kernel_simd_func3d)
			     dlp_kernel_simd_laplacebem3d);
}
#else
static void
fill_dlp_lc_far_laplacebem3d(const uint * ridx, const uint * cidx,
			     pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_lc_far_bem3d(ridx, cidx, bem, ntrans, N,
			(kernel_func3d) dlp_kernel_laplacebem3d);
}
#endif

#ifdef USE_SIMD
static void
fill_adlp_lc_near_laplacebem3d(const uint * ridx, const uint * cidx,
			       pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_lc_simd_near_bem3d(ridx, cidx, bem, ntrans, N, (kernel_simd_func3d)
			      adlp_kernel_simd_laplacebem3d);
}
#else
static void
fill_adlp_lc_near_laplacebem3d(const uint * ridx, const uint * cidx,
			       pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_lc_near_bem3d(ridx, cidx, bem, ntrans, N,
			 (kernel_func3d) adlp_kernel_laplacebem3d);
}
#endif

#ifdef USE_SIMD
static void
fill_adlp_lc_far_laplacebem3d(const uint * ridx, const uint * cidx,
			      pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_lc_simd_far_bem3d(ridx, cidx, bem, ntrans, N, (kernel_simd_func3d)
			     adlp_kernel_simd_laplacebem3d);
}
#else
static void
fill_adlp_lc_far_laplacebem3d(const uint * ridx, const uint * cidx,
			      pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_lc_far_bem3d(ridx, cidx, bem, ntrans, N,
			(kernel_func3d) adlp_kernel_laplacebem3d);
}
#endif

#ifdef USE_SIMD
static void
fill_slp_ll_near_laplacebem3d(const uint * ridx, const uint * cidx,
			      pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_ll_simd_near_bem3d(ridx, cidx, bem, ntrans, N, (kernel_simd_func3d)
			      slp_kernel_simd_laplacebem3d);
}
#else
static void
fill_slp_ll_near_laplacebem3d(const uint * ridx, const uint * cidx,
			      pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_ll_near_bem3d(ridx, cidx, bem, ntrans, N,
			 (kernel_func3d) slp_kernel_laplacebem3d);
}
#endif

#ifdef USE_SIMD
static void
fill_slp_ll_far_laplacebem3d(const uint * ridx, const uint * cidx,
			     pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_ll_simd_far_bem3d(ridx, cidx, bem, ntrans, N, (kernel_simd_func3d)
			     slp_kernel_simd_laplacebem3d);
}
#else
static void
fill_slp_ll_far_laplacebem3d(const uint * ridx, const uint * cidx,
			     pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_ll_far_bem3d(ridx, cidx, bem, ntrans, N,
			(kernel_func3d) slp_kernel_laplacebem3d);
}
#endif

#ifdef USE_SIMD
static void
fill_dlp_ll_near_laplacebem3d(const uint * ridx, const uint * cidx,
			      pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_ll_simd_near_bem3d(ridx, cidx, bem, ntrans, N, (kernel_simd_func3d)
			      dlp_kernel_simd_laplacebem3d);
}
#else
static void
fill_dlp_ll_near_laplacebem3d(const uint * ridx, const uint * cidx,
			      pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_ll_near_bem3d(ridx, cidx, bem, ntrans, N,
			 (kernel_func3d) dlp_kernel_laplacebem3d);
}
#endif

#ifdef USE_SIMD
static void
fill_dlp_ll_far_laplacebem3d(const uint * ridx, const uint * cidx,
			     pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_ll_simd_far_bem3d(ridx, cidx, bem, ntrans, N, (kernel_simd_func3d)
			     dlp_kernel_simd_laplacebem3d);
}
#else
static void
fill_dlp_ll_far_laplacebem3d(const uint * ridx, const uint * cidx,
			     pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_ll_far_bem3d(ridx, cidx, bem, ntrans, N,
			(kernel_func3d) dlp_kernel_laplacebem3d);
}
#endif

#ifdef USE_SIMD
static void
fill_adlp_ll_near_laplacebem3d(const uint * ridx, const uint * cidx,
			       pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_ll_simd_near_bem3d(ridx, cidx, bem, ntrans, N, (kernel_simd_func3d)
			      adlp_kernel_simd_laplacebem3d);
}
#else
static void
fill_adlp_ll_near_laplacebem3d(const uint * ridx, const uint * cidx,
			       pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_ll_near_bem3d(ridx, cidx, bem, ntrans, N,
			 (kernel_func3d) adlp_kernel_laplacebem3d);
}
#endif

#ifdef USE_SIMD
static void
fill_adlp_ll_far_laplacebem3d(const uint * ridx, const uint * cidx,
			      pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_ll_simd_far_bem3d(ridx, cidx, bem, ntrans, N, (kernel_simd_func3d)
			     adlp_kernel_simd_laplacebem3d);
}
#else
static void
fill_adlp_ll_far_laplacebem3d(const uint * ridx, const uint * cidx,
			      pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_ll_far_bem3d(ridx, cidx, bem, ntrans, N,
			(kernel_func3d) adlp_kernel_laplacebem3d);
}
#endif

static void
fill_kernel_laplacebem3d(pcbem3d bem, const real(*X)[3],
			 const real(*Y)[3], pamatrix V)
{
  fill_bem3d(bem, X, Y, NULL, NULL, V, slp_kernel_laplacebem3d);
}

static void
fill_dny_kernel_laplacebem3d(pcbem3d bem, const real(*X)[3],
			     const real(*Y)[3], const real(*NY)[3],
			     pamatrix V)
{
  fill_bem3d(bem, X, Y, NULL, NY, V, dlp_kernel_laplacebem3d);
}

static void
fill_dnx_dny_kernel_laplacebem3d(pcbem3d bem, const real(*X)[3],
				 const real(*NX)[3], const real(*Y)[3],
				 const real(*NY)[3], pamatrix V)
{
  fill_bem3d(bem, X, Y, NX, NY, V, hs_kernel_laplacebem3d);
}

#ifdef USE_SIMD
static void
fill_kernel_c_laplacebem3d(const uint * idx, const real(*Z)[3],
			   pcbem3d bem, pamatrix V)
{
  fill_row_simd_c_bem3d(idx, Z, bem, V, slp_kernel_simd_laplacebem3d);
}
#else
static void
fill_kernel_c_laplacebem3d(const uint * idx, const real(*Z)[3],
			   pcbem3d bem, pamatrix V)
{
  fill_row_c_bem3d(idx, Z, bem, V, slp_kernel_laplacebem3d);
}
#endif

static void
fill_kernel_l_laplacebem3d(const uint * idx, const real(*Z)[3],
			   pcbem3d bem, pamatrix V)
{
  fill_row_l_bem3d(idx, Z, bem, V, slp_kernel_laplacebem3d);
}

#ifdef USE_SIMD
static void
fill_dnz_kernel_c_laplacebem3d(const uint * idx, const real(*Z)[3],
			       const real(*N)[3], pcbem3d bem, pamatrix V)
{
  fill_dnz_row_simd_c_bem3d(idx, Z, N, bem, V, dlp_kernel_simd_laplacebem3d);
}
#else
static void
fill_dnz_kernel_c_laplacebem3d(const uint * idx,
			       const real(*Z)[3], const real(*N)[3],
			       pcbem3d bem, pamatrix V)
{
  fill_dnz_row_c_bem3d(idx, Z, N, bem, V, dlp_kernel_laplacebem3d);
}
#endif

static void
fill_dnz_kernel_l_laplacebem3d(const uint * idx, const real(*Z)[3],
			       const real(*N)[3], pcbem3d bem, pamatrix V)
{
  fill_dnz_row_l_bem3d(idx, Z, N, bem, V, dlp_kernel_laplacebem3d);
}

#ifdef USE_SIMD
static void
fill_dnzdrow_kernel_c_laplacebem3d(const uint * idx,
				   const real(*Z)[3], const real(*N)[3],
				   pcbem3d bem, pamatrix V)
{
  fill_dnz_row_simd_c_bem3d(idx, Z, N, bem, V, hs_kernel_simd_laplacebem3d);
}
#else
static void
fill_dnzdrow_kernel_c_laplacebem3d(const uint * idx,
				   const real(*Z)[3], const real(*N)[3],
				   pcbem3d bem, pamatrix V)
{
  fill_dnz_row_c_bem3d(idx, Z, N, bem, V, hs_kernel_laplacebem3d);
}
#endif

static void
fill_dnzdrow_kernel_l_laplacebem3d(const uint * idx,
				   const real(*Z)[3], const real(*N)[3],
				   pcbem3d bem, pamatrix V)
{
  fill_dnz_row_l_bem3d(idx, Z, N, bem, V, hs_kernel_laplacebem3d);
}

#ifdef USE_SIMD
static void
fill_dnzdcol_kernel_c_laplacebem3d(const uint * idx,
				   const real(*Z)[3], const real(*N)[3],
				   pcbem3d bem, pamatrix V)
{
  fill_dnz_col_simd_c_bem3d(idx, Z, N, bem, V, hs_kernel_simd_laplacebem3d);
}
#else
static void
fill_dnzdcol_kernel_c_laplacebem3d(const uint * idx,
				   const real(*Z)[3], const real(*N)[3],
				   pcbem3d bem, pamatrix V)
{
  fill_dnz_col_c_bem3d(idx, Z, N, bem, V, hs_kernel_laplacebem3d);
}
#endif

static void
fill_dnzdcol_kernel_l_laplacebem3d(const uint * idx,
				   const real(*Z)[3], const real(*N)[3],
				   pcbem3d bem, pamatrix V)
{
  fill_dnz_col_l_bem3d(idx, Z, N, bem, V, hs_kernel_laplacebem3d);
}

#ifdef USE_SIMD
static void
fill_drow_kernel_c_laplacebem3d(const uint * idx,
				const real(*Z)[3], pcbem3d bem, pamatrix V)
{
  fill_row_simd_c_bem3d(idx, Z, bem, V, adlp_kernel_simd_laplacebem3d);
}
#else
static void
fill_drow_kernel_c_laplacebem3d(const uint * idx,
				const real(*Z)[3], pcbem3d bem, pamatrix V)
{
  fill_row_c_bem3d(idx, Z, bem, V, adlp_kernel_laplacebem3d);
}
#endif

static void
fill_drow_kernel_l_laplacebem3d(const uint * idx,
				const real(*Z)[3], pcbem3d bem, pamatrix V)
{
  fill_row_l_bem3d(idx, Z, bem, V, adlp_kernel_laplacebem3d);
}

#ifdef USE_SIMD
static void
fill_dcol_kernel_c_laplacebem3d(const uint * idx,
				const real(*Z)[3], pcbem3d bem, pamatrix V)
{
  fill_col_simd_c_bem3d(idx, Z, bem, V, dlp_kernel_simd_laplacebem3d);
}

#else
static void
fill_dcol_kernel_c_laplacebem3d(const uint * idx,
				const real(*Z)[3], pcbem3d bem, pamatrix V)
{
  fill_col_c_bem3d(idx, Z, bem, V, dlp_kernel_laplacebem3d);
}
#endif

static void
fill_dcol_kernel_l_laplacebem3d(const uint * idx,
				const real(*Z)[3], pcbem3d bem, pamatrix V)
{
  fill_col_l_bem3d(idx, Z, bem, V, dlp_kernel_laplacebem3d);
}

pbem3d
new_slp_laplace_bem3d(pcsurface3d gr, uint q_regular, uint q_singular,
		      basisfunctionbem3d row_basis,
		      basisfunctionbem3d col_basis)
{
  pkernelbem3d kernels;

  pbem3d    bem;

  bem = new_bem3d(gr, row_basis, col_basis);
  kernels = bem->kernels;

  bem->sq = build_singquad2d(gr, q_regular, q_singular);

  bem->kernel_const = KERNEL_CONST_LAPLACEBEM3D;
  bem->alpha = 0.0;

  kernels->fundamental = fill_kernel_laplacebem3d;
  kernels->fundamental_wave = NULL;
  kernels->dny_fundamental = fill_dny_kernel_laplacebem3d;
  kernels->dnx_dny_fundamental = fill_dnx_dny_kernel_laplacebem3d;

  switch (row_basis) {
  case BASIS_CONSTANT_BEM3D:
    kernels->lagrange_row = assemble_bem3d_lagrange_c_amatrix;
    kernels->lagrange_wave_row = NULL;

    kernels->fundamental_row = fill_kernel_c_laplacebem3d;
    kernels->dnz_fundamental_row = fill_dnz_kernel_c_laplacebem3d;
    kernels->kernel_row = fill_kernel_c_laplacebem3d;
    kernels->dnz_kernel_row = fill_dnz_kernel_c_laplacebem3d;
    break;
  case BASIS_LINEAR_BEM3D:
    kernels->lagrange_row = assemble_bem3d_lagrange_l_amatrix;
    kernels->lagrange_wave_row = NULL;

    kernels->fundamental_row = fill_kernel_l_laplacebem3d;
    kernels->dnz_fundamental_row = fill_dnz_kernel_l_laplacebem3d;
    kernels->kernel_row = fill_kernel_l_laplacebem3d;
    kernels->dnz_kernel_row = fill_dnz_kernel_l_laplacebem3d;
    break;
  default:
    fprintf(stderr, "Unknown basis type detected!\n");
    abort();
    break;
  }

  switch (col_basis) {
  case BASIS_CONSTANT_BEM3D:
    kernels->lagrange_col = assemble_bem3d_lagrange_c_amatrix;
    kernels->lagrange_wave_col = NULL;

    kernels->fundamental_col = fill_kernel_c_laplacebem3d;
    kernels->dnz_fundamental_col = fill_dnz_kernel_c_laplacebem3d;
    kernels->kernel_col = fill_kernel_c_laplacebem3d;
    kernels->dnz_kernel_col = fill_dnz_kernel_c_laplacebem3d;
    break;
  case BASIS_LINEAR_BEM3D:
    kernels->lagrange_col = assemble_bem3d_lagrange_l_amatrix;
    kernels->lagrange_wave_col = NULL;

    kernels->fundamental_col = fill_kernel_l_laplacebem3d;
    kernels->dnz_fundamental_col = fill_dnz_kernel_l_laplacebem3d;
    kernels->kernel_col = fill_kernel_l_laplacebem3d;
    kernels->dnz_kernel_col = fill_dnz_kernel_l_laplacebem3d;
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
      bem->nearfield = fill_slp_cc_near_laplacebem3d;
      bem->nearfield_far = fill_slp_cc_far_laplacebem3d;
      break;
    case BASIS_LINEAR_BEM3D:
      bem->nearfield = fill_slp_cl_near_laplacebem3d;
      bem->nearfield_far = fill_slp_cl_far_laplacebem3d;

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
      bem->nearfield = fill_slp_lc_near_laplacebem3d;
      bem->nearfield_far = fill_slp_lc_far_laplacebem3d;

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
      bem->nearfield = fill_slp_ll_near_laplacebem3d;
      bem->nearfield_far = fill_slp_ll_far_laplacebem3d;

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
new_dlp_laplace_bem3d(pcsurface3d gr, uint q_regular, uint q_singular,
		      basisfunctionbem3d row_basis,
		      basisfunctionbem3d col_basis, field alpha)
{
  pkernelbem3d kernels;

  pbem3d    bem;

  bem = new_bem3d(gr, row_basis, col_basis);
  kernels = bem->kernels;

  bem->sq = build_singquad2d(gr, q_regular, q_singular);

  bem->kernel_const = KERNEL_CONST_LAPLACEBEM3D;
  bem->alpha = alpha;

  kernels->fundamental = fill_kernel_laplacebem3d;
  kernels->fundamental_wave = NULL;
  kernels->dny_fundamental = fill_dny_kernel_laplacebem3d;
  kernels->dnx_dny_fundamental = fill_dnx_dny_kernel_laplacebem3d;

  switch (row_basis) {
  case BASIS_CONSTANT_BEM3D:
    kernels->lagrange_row = assemble_bem3d_lagrange_c_amatrix;
    kernels->lagrange_wave_row = NULL;

    kernels->fundamental_row = fill_kernel_c_laplacebem3d;
    kernels->dnz_fundamental_row = fill_dnz_kernel_c_laplacebem3d;
    kernels->kernel_row = fill_kernel_c_laplacebem3d;
    kernels->dnz_kernel_row = fill_dnz_kernel_c_laplacebem3d;
    break;
  case BASIS_LINEAR_BEM3D:
    kernels->lagrange_row = assemble_bem3d_lagrange_l_amatrix;
    kernels->lagrange_wave_row = NULL;

    kernels->fundamental_row = fill_kernel_l_laplacebem3d;
    kernels->dnz_fundamental_row = fill_dnz_kernel_l_laplacebem3d;
    kernels->kernel_row = fill_kernel_l_laplacebem3d;
    kernels->dnz_kernel_row = fill_dnz_kernel_l_laplacebem3d;
    break;
  default:
    fprintf(stderr, "Unknown basis type detected!\n");
    abort();
    break;
  }

  switch (col_basis) {
  case BASIS_CONSTANT_BEM3D:
    kernels->lagrange_col = assemble_bem3d_dn_lagrange_c_amatrix;
    kernels->lagrange_wave_col = NULL;

    kernels->fundamental_col = fill_kernel_c_laplacebem3d;
    kernels->dnz_fundamental_col = fill_dnz_kernel_c_laplacebem3d;
    kernels->kernel_col = fill_dcol_kernel_c_laplacebem3d;
    kernels->dnz_kernel_col = fill_dnzdcol_kernel_c_laplacebem3d;
    break;
  case BASIS_LINEAR_BEM3D:
    kernels->lagrange_col = assemble_bem3d_dn_lagrange_l_amatrix;
    kernels->lagrange_wave_col = NULL;

    kernels->fundamental_col = fill_kernel_l_laplacebem3d;
    kernels->dnz_fundamental_col = fill_dnz_kernel_l_laplacebem3d;
    kernels->kernel_col = fill_dcol_kernel_l_laplacebem3d;
    kernels->dnz_kernel_col = fill_dnzdcol_kernel_l_laplacebem3d;
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
      bem->nearfield = fill_dlp_cc_near_laplacebem3d;
      bem->nearfield_far = fill_dlp_cc_far_laplacebem3d;
      break;
    case BASIS_LINEAR_BEM3D:
      bem->nearfield = fill_dlp_cl_near_laplacebem3d;
      bem->nearfield_far = fill_dlp_cl_far_laplacebem3d;

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
      bem->nearfield = fill_dlp_lc_near_laplacebem3d;
      bem->nearfield_far = fill_dlp_lc_far_laplacebem3d;

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
      bem->nearfield = fill_dlp_ll_near_laplacebem3d;
      bem->nearfield_far = fill_dlp_ll_far_laplacebem3d;

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
new_adlp_laplace_bem3d(pcsurface3d gr, uint q_regular, uint q_singular,
		       basisfunctionbem3d row_basis,
		       basisfunctionbem3d col_basis, field alpha)
{
  pkernelbem3d kernels;

  pbem3d    bem;

  bem = new_bem3d(gr, row_basis, col_basis);
  kernels = bem->kernels;

  bem->sq = build_singquad2d(gr, q_regular, q_singular);

  bem->kernel_const = KERNEL_CONST_LAPLACEBEM3D;
  bem->alpha = alpha;

  kernels->fundamental = fill_kernel_laplacebem3d;
  kernels->fundamental_wave = NULL;
  kernels->dny_fundamental = fill_dny_kernel_laplacebem3d;
  kernels->dnx_dny_fundamental = fill_dnx_dny_kernel_laplacebem3d;

  switch (row_basis) {
  case BASIS_CONSTANT_BEM3D:
    kernels->lagrange_row = assemble_bem3d_dn_lagrange_c_amatrix;
    kernels->lagrange_wave_row = NULL;

    kernels->fundamental_row = fill_kernel_c_laplacebem3d;
    kernels->dnz_fundamental_row = fill_dnz_kernel_c_laplacebem3d;
    kernels->kernel_row = fill_drow_kernel_c_laplacebem3d;
    kernels->dnz_kernel_row = fill_dnzdrow_kernel_c_laplacebem3d;
    break;
  case BASIS_LINEAR_BEM3D:
    kernels->lagrange_row = assemble_bem3d_dn_lagrange_l_amatrix;
    kernels->lagrange_wave_row = NULL;

    kernels->fundamental_row = fill_kernel_l_laplacebem3d;
    kernels->dnz_fundamental_row = fill_dnz_kernel_l_laplacebem3d;
    kernels->kernel_row = fill_drow_kernel_l_laplacebem3d;
    kernels->dnz_kernel_row = fill_dnzdrow_kernel_l_laplacebem3d;
    break;
  default:
    fprintf(stderr, "Unknown basis type detected!\n");
    abort();
    break;
  }

  switch (col_basis) {
  case BASIS_CONSTANT_BEM3D:
    kernels->lagrange_col = assemble_bem3d_lagrange_c_amatrix;
    kernels->lagrange_wave_col = NULL;

    kernels->fundamental_col = fill_kernel_c_laplacebem3d;
    kernels->dnz_fundamental_col = fill_dnz_kernel_c_laplacebem3d;
    kernels->kernel_col = fill_kernel_c_laplacebem3d;
    kernels->dnz_kernel_col = fill_dnz_kernel_c_laplacebem3d;
    break;
  case BASIS_LINEAR_BEM3D:
    kernels->lagrange_col = assemble_bem3d_lagrange_l_amatrix;
    kernels->lagrange_wave_col = NULL;

    kernels->fundamental_col = fill_kernel_l_laplacebem3d;
    kernels->dnz_fundamental_col = fill_dnz_kernel_l_laplacebem3d;
    kernels->kernel_col = fill_kernel_l_laplacebem3d;
    kernels->dnz_kernel_col = fill_dnz_kernel_l_laplacebem3d;
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
      bem->nearfield = fill_adlp_cc_near_laplacebem3d;
      bem->nearfield_far = fill_adlp_cc_far_laplacebem3d;
      break;
    case BASIS_LINEAR_BEM3D:
      bem->nearfield = fill_adlp_cl_near_laplacebem3d;
      bem->nearfield_far = fill_adlp_cl_far_laplacebem3d;

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
      bem->nearfield = fill_adlp_lc_near_laplacebem3d;
      bem->nearfield_far = fill_adlp_lc_far_laplacebem3d;

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
      bem->nearfield = fill_adlp_ll_near_laplacebem3d;
      bem->nearfield_far = fill_adlp_ll_far_laplacebem3d;

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
del_laplace_bem3d(pbem3d bem)
{
  del_bem3d(bem);
}

field
eval_dirichlet_linear_laplacebem3d(const real * x, const real * n, void *data)
{
  (void) n;
  (void) data;

  return (x[0] + x[1] + x[2]);
}

field
eval_neumann_linear_laplacebem3d(const real * x, const real * n, void *data)
{
  (void) x;
  (void) data;

  return (n[0] + n[1] + n[2]);
}

field
eval_dirichlet_quadratic_laplacebem3d(const real * x, const real * n,
				      void *data)
{
  (void) n;
  (void) data;

  return (x[0] * x[0] - x[2] * x[2]);
}

field
eval_neumann_quadratic_laplacebem3d(const real * x, const real * n,
				    void *data)
{
  (void) data;

  return 2.0 * (n[0] * x[0] - n[2] * x[2]);
}

field
eval_dirichlet_fundamental_laplacebem3d(const real * x, const real * n,
					void *data)
{
  real      d[3];

  field     kernel;

  (void) n;
  (void) data;

  d[0] = x[0] - (1.2);
  d[1] = x[1] - (1.2);
  d[2] = x[2] - (1.2);

  kernel = REAL_SQR(d[0]) + REAL_SQR(d[1]) + REAL_SQR(d[2]);
  kernel = REAL_RSQRT(kernel);

  return kernel;
}

field
eval_neumann_fundamental_laplacebem3d(const real * x, const real * n,
				      void *data)
{
  real      d[3];

  field     kernel;

  (void) data;

  d[0] = x[0] - (1.2);
  d[1] = x[1] - (1.2);
  d[2] = x[2] - (1.2);

  kernel = REAL_SQR(d[0]) + REAL_SQR(d[1]) + REAL_SQR(d[2]);
  kernel = -(d[0] * n[0] + d[1] * n[1] + d[2] * n[2])
    / (kernel * REAL_SQRT(kernel));

  return kernel;
}

field
eval_dirichlet_fundamental2_laplacebem3d(const real * x, const real * n,
					 void *data)
{
  real      d[3];

  field     kernel;

  (void) n;
  (void) data;

  d[0] = x[0] - (0.25);
  d[1] = x[1] - (0.25);
  d[2] = x[2] - (0.25);

  kernel = REAL_SQR(d[0]) + REAL_SQR(d[1]) + REAL_SQR(d[2]);
  kernel = REAL_RSQRT(kernel);

  return kernel;
}

field
eval_neumann_fundamental2_laplacebem3d(const real * x, const real * n,
				       void *data)
{
  real      d[3];

  field     kernel;

  (void) data;

  d[0] = x[0] - (0.25);
  d[1] = x[1] - (0.25);
  d[2] = x[2] - (0.25);

  kernel = REAL_RSQRT(REAL_SQR(d[0]) + REAL_SQR(d[1]) + REAL_SQR(d[2]));
  kernel = -(d[0] * n[0] + d[1] * n[1] + d[2] * n[2])
    * (kernel * kernel * kernel);

  return kernel;
}
