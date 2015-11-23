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

static void
fill_slp_cc_near_laplacebem3d(const uint * ridx, const uint * cidx,
			      pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_cc_near_bem3d(ridx, cidx, bem, ntrans, N,
			 (kernel_func3d) slp_kernel_laplacebem3d);
}

static void
fill_slp_cc_far_laplacebem3d(const uint * ridx, const uint * cidx,
			     pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_cc_far_bem3d(ridx, cidx, bem, ntrans, N,
			(kernel_func3d) slp_kernel_laplacebem3d);
}

static void
fill_dlp_cc_near_laplacebem3d(const uint * ridx, const uint * cidx,
			      pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_cc_near_bem3d(ridx, cidx, bem, ntrans, N,
			 (kernel_func3d) dlp_kernel_laplacebem3d);
}

static void
fill_dlp_cc_far_laplacebem3d(const uint * ridx, const uint * cidx,
			     pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_cc_far_bem3d(ridx, cidx, bem, ntrans, N,
			(kernel_func3d) dlp_kernel_laplacebem3d);
}

static void
fill_dlp_cl_near_laplacebem3d(const uint * ridx, const uint * cidx,
			      pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_cl_near_bem3d(ridx, cidx, bem, ntrans, N, dlp_kernel_laplacebem3d);
}

static void
fill_dlp_cl_far_laplacebem3d(const uint * ridx, const uint * cidx,
			     pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_cl_far_bem3d(ridx, cidx, bem, ntrans, N, dlp_kernel_laplacebem3d);
}

static void
fill_slp_ll_near_laplacebem3d(const uint * ridx, const uint * cidx,
			      pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_ll_near_bem3d(ridx, cidx, bem, ntrans, N, slp_kernel_laplacebem3d);
}

static void
fill_slp_ll_far_laplacebem3d(const uint * ridx, const uint * cidx,
			     pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_ll_far_bem3d(ridx, cidx, bem, ntrans, N, slp_kernel_laplacebem3d);
}

static void
fill_dlp_ll_near_laplacebem3d(const uint * ridx, const uint * cidx,
			      pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_ll_near_bem3d(ridx, cidx, bem, ntrans, N, dlp_kernel_laplacebem3d);
}

static void
fill_dlp_ll_far_laplacebem3d(const uint * ridx, const uint * cidx,
			     pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_ll_far_bem3d(ridx, cidx, bem, ntrans, N, dlp_kernel_laplacebem3d);
}

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

static void
fill_kernel_c_laplacebem3d(const uint * idx, const real(*Z)[3],
			   pcbem3d bem, pamatrix V)
{
  fill_row_c_bem3d(idx, Z, bem, V, slp_kernel_laplacebem3d);
}

static void
fill_kernel_l_laplacebem3d(const uint * idx, const real(*Z)[3],
			   pcbem3d bem, pamatrix V)
{
  fill_row_l_bem3d(idx, Z, bem, V, slp_kernel_laplacebem3d);
}

static void
fill_dnz_kernel_c_laplacebem3d(const uint * idx, const real(*Z)[3],
			       const real(*N)[3], pcbem3d bem, pamatrix V)
{
  fill_dnz_row_c_bem3d(idx, Z, N, bem, V, dlp_kernel_laplacebem3d);
}

static void
fill_dnz_kernel_l_laplacebem3d(const uint * idx, const real(*Z)[3],
			       const real(*N)[3], pcbem3d bem, pamatrix V)
{
  fill_dnz_row_l_bem3d(idx, Z, N, bem, V, dlp_kernel_laplacebem3d);
}

static void
fill_dnzdcol_kernel_c_laplacebem3d(const uint * idx,
				   const real(*Z)[3], const real(*N)[3],
				   pcbem3d bem, pamatrix V)
{
  fill_dnz_col_c_bem3d(idx, Z, N, bem, V, hs_kernel_laplacebem3d);
}

static void
fill_dnzdcol_kernel_l_laplacebem3d(const uint * idx,
				   const real(*Z)[3], const real(*N)[3],
				   pcbem3d bem, pamatrix V)
{
  fill_dnz_col_l_bem3d(idx, Z, N, bem, V, hs_kernel_laplacebem3d);
}

static void
fill_dcol_kernel_col_c_laplacebem3d(const uint * idx,
				    const real(*Z)[3], pcbem3d bem,
				    pamatrix V)
{
  fill_col_c_bem3d(idx, Z, bem, V, dlp_kernel_laplacebem3d);
}

static void
fill_dcol_kernel_col_l_laplacebem3d(const uint * idx,
				    const real(*Z)[3], pcbem3d bem,
				    pamatrix V)
{
  fill_col_l_bem3d(idx, Z, bem, V, dlp_kernel_laplacebem3d);
}

pbem3d
new_slp_laplace_bem3d(pcsurface3d gr, uint q_regular, uint q_singular,
		      basisfunctionbem3d basis)
{
  pkernelbem3d kernels;

  pbem3d    bem;

  bem = new_bem3d(gr);
  kernels = bem->kernels;

  bem->sq = build_singquad2d(gr, q_regular, q_singular);

  if (basis == BASIS_LINEAR_BEM3D) {
    setup_vertex_to_triangle_map_bem3d(bem);
  }

  bem->basis_neumann = basis;

  bem->kernel_const = KERNEL_CONST_LAPLACEBEM3D;
  bem->alpha = 0.0;

  kernels->fundamental = fill_kernel_laplacebem3d;
  kernels->dny_fundamental = fill_dny_kernel_laplacebem3d;
  kernels->dnx_dny_fundamental = fill_dnx_dny_kernel_laplacebem3d;

  if (basis == BASIS_CONSTANT_BEM3D) {
    bem->N_neumann = gr->triangles;

    bem->nearfield = fill_slp_cc_near_laplacebem3d;
    bem->nearfield_far = fill_slp_cc_far_laplacebem3d;

    kernels->lagrange_row = assemble_bem3d_lagrange_const_amatrix;
    kernels->lagrange_col = assemble_bem3d_lagrange_const_amatrix;

    kernels->fundamental_row = fill_kernel_c_laplacebem3d;
    kernels->fundamental_col = fill_kernel_c_laplacebem3d;
    kernels->dnz_fundamental_row = fill_dnz_kernel_c_laplacebem3d;
    kernels->dnz_fundamental_col = fill_dnz_kernel_c_laplacebem3d;

    kernels->kernel_row = fill_kernel_c_laplacebem3d;
    kernels->kernel_col = fill_kernel_c_laplacebem3d;
    kernels->dnz_kernel_row = fill_dnz_kernel_c_laplacebem3d;
    kernels->dnz_kernel_col = fill_dnz_kernel_c_laplacebem3d;
  }
  else {
    assert(basis == BASIS_LINEAR_BEM3D);

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

    bem->N_neumann = gr->vertices;

    bem->nearfield = fill_slp_ll_near_laplacebem3d;
    bem->nearfield_far = fill_slp_ll_far_laplacebem3d;

    kernels->lagrange_row = assemble_bem3d_lagrange_linear_amatrix;
    kernels->lagrange_col = assemble_bem3d_lagrange_linear_amatrix;

    kernels->fundamental_row = fill_kernel_l_laplacebem3d;
    kernels->fundamental_col = fill_kernel_l_laplacebem3d;
    kernels->dnz_fundamental_row = fill_dnz_kernel_l_laplacebem3d;
    kernels->dnz_fundamental_col = fill_dnz_kernel_l_laplacebem3d;

    kernels->kernel_row = fill_kernel_l_laplacebem3d;
    kernels->kernel_col = fill_kernel_l_laplacebem3d;
    kernels->dnz_kernel_row = fill_dnz_kernel_l_laplacebem3d;
    kernels->dnz_kernel_col = fill_dnz_kernel_l_laplacebem3d;
  }

  return bem;
}

pbem3d
new_dlp_laplace_bem3d(pcsurface3d gr, uint q_regular, uint q_singular,
		      basisfunctionbem3d basis_neumann,
		      basisfunctionbem3d basis_dirichlet, field alpha)
{
  pkernelbem3d kernels;

  pbem3d    bem;

  bem = new_bem3d(gr);
  kernels = bem->kernels;

  bem->sq = build_singquad2d(gr, q_regular, q_singular);

  if (basis_neumann == BASIS_LINEAR_BEM3D || basis_dirichlet
      == BASIS_LINEAR_BEM3D) {
    setup_vertex_to_triangle_map_bem3d(bem);
  }

  bem->basis_neumann = basis_neumann;
  bem->basis_dirichlet = basis_dirichlet;

  bem->kernel_const = KERNEL_CONST_LAPLACEBEM3D;
  bem->alpha = alpha;

  kernels->fundamental = fill_kernel_laplacebem3d;
  kernels->dny_fundamental = fill_dny_kernel_laplacebem3d;
  kernels->dnx_dny_fundamental = fill_dnx_dny_kernel_laplacebem3d;

  if (basis_neumann == BASIS_CONSTANT_BEM3D) {
    bem->N_neumann = gr->triangles;
  }
  else {
    assert(basis_neumann == BASIS_LINEAR_BEM3D);
    bem->N_neumann = gr->vertices;
  }
  if (basis_dirichlet == BASIS_CONSTANT_BEM3D) {
    bem->N_dirichlet = gr->triangles;
  }
  else {
    assert(basis_dirichlet == BASIS_LINEAR_BEM3D);
    bem->N_dirichlet = gr->vertices;
  }

  if (basis_neumann == BASIS_CONSTANT_BEM3D && basis_dirichlet
      == BASIS_CONSTANT_BEM3D) {
    bem->nearfield = fill_dlp_cc_near_laplacebem3d;
    bem->nearfield_far = fill_dlp_cc_far_laplacebem3d;

    kernels->lagrange_row = assemble_bem3d_lagrange_const_amatrix;
    kernels->lagrange_col = assemble_bem3d_dn_lagrange_const_amatrix;

    kernels->fundamental_row = fill_kernel_c_laplacebem3d;
    kernels->fundamental_col = fill_kernel_c_laplacebem3d;
    kernels->dnz_fundamental_row = fill_dnz_kernel_c_laplacebem3d;
    kernels->dnz_fundamental_col = fill_dnz_kernel_c_laplacebem3d;

    kernels->kernel_row = fill_kernel_c_laplacebem3d;
    kernels->kernel_col = fill_dcol_kernel_col_c_laplacebem3d;
    kernels->dnz_kernel_row = fill_dnz_kernel_c_laplacebem3d;
    kernels->dnz_kernel_col = fill_dnzdcol_kernel_c_laplacebem3d;
  }
  else if (basis_neumann == BASIS_LINEAR_BEM3D
	   && basis_dirichlet == BASIS_CONSTANT_BEM3D) {
    bem->mass = allocreal(3);
    /* TODO MASS-MATRIX */

  }
  else if (basis_neumann == BASIS_CONSTANT_BEM3D
	   && basis_dirichlet == BASIS_LINEAR_BEM3D) {

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

    bem->nearfield = fill_dlp_cl_near_laplacebem3d;
    bem->nearfield_far = fill_dlp_cl_far_laplacebem3d;

    kernels->lagrange_row = assemble_bem3d_lagrange_const_amatrix;
    kernels->lagrange_col = assemble_bem3d_dn_lagrange_linear_amatrix;

    kernels->fundamental_row = fill_kernel_c_laplacebem3d;
    kernels->fundamental_col = fill_kernel_l_laplacebem3d;
    kernels->dnz_fundamental_row = fill_dnz_kernel_c_laplacebem3d;
    kernels->dnz_fundamental_col = fill_dnz_kernel_l_laplacebem3d;

    kernels->kernel_row = fill_kernel_c_laplacebem3d;
    kernels->kernel_col = fill_dcol_kernel_col_l_laplacebem3d;
    kernels->dnz_kernel_row = fill_dnz_kernel_c_laplacebem3d;
    kernels->dnz_kernel_col = fill_dnzdcol_kernel_l_laplacebem3d;

    bem->mass = allocreal(3);
    bem->mass[0] = 1.0 / 6.0;
    bem->mass[1] = 1.0 / 6.0;
    bem->mass[2] = 1.0 / 6.0;

  }
  else {
    assert(basis_neumann == BASIS_LINEAR_BEM3D && basis_dirichlet
	   == BASIS_LINEAR_BEM3D);

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

    bem->nearfield = fill_dlp_ll_near_laplacebem3d;
    bem->nearfield_far = fill_dlp_ll_far_laplacebem3d;

    kernels->lagrange_row = assemble_bem3d_lagrange_linear_amatrix;
    kernels->lagrange_col = assemble_bem3d_dn_lagrange_linear_amatrix;

    kernels->fundamental_row = fill_kernel_l_laplacebem3d;
    kernels->fundamental_col = fill_kernel_l_laplacebem3d;
    kernels->dnz_fundamental_row = fill_dnz_kernel_l_laplacebem3d;
    kernels->dnz_fundamental_col = fill_dnz_kernel_l_laplacebem3d;

    kernels->kernel_row = fill_kernel_l_laplacebem3d;
    kernels->kernel_col = fill_dcol_kernel_col_l_laplacebem3d;
    kernels->dnz_kernel_row = NULL;
    kernels->dnz_kernel_col = fill_dnzdcol_kernel_l_laplacebem3d;

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
