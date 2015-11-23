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
  real      k = bem->k;
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

  norm = k * norm2 * rnorm;
  res = (cos(norm) + I * sin(norm)) * rnorm;

  return res;
}

static inline field
dlp_kernel_helmholtzbem3d(const real * x, const real * y,
			  const real * nx, const real * ny, void *data)
{
  pcbem3d   bem = (pcbem3d) data;
  real      k = bem->k;
  real      dist[3];
  real      norm, norm2, rnorm, s, c;

  field     res;

  (void) nx;

  dist[0] = x[0] - y[0];
  dist[1] = x[1] - y[1];
  dist[2] = x[2] - y[2];
  norm2 = REAL_NORMSQR3(dist[0], dist[1], dist[2]);
  rnorm = REAL_RSQRT(norm2);

  norm = k * norm2 * rnorm;
  s = sin(norm);
  c = cos(norm);
  res = (c + norm * s + (s - c * norm) * I);
  res *= (rnorm * rnorm * rnorm) * DOT3(dist, ny);

  return res;
}

static inline field
hs_kernel_helmholtzbem3d(const real * x, const real * y,
			 const real * nx, const real * ny, void *data)
{
  pcbem3d   bem = (pcbem3d) data;
  real      k = bem->k;
  real      dist[3];
  real      norm, norm2, rnorm, s, c, dot, dotxy, hr, hi;

  field     res;

  dist[0] = x[0] - y[0];
  dist[1] = x[1] - y[1];
  dist[2] = x[2] - y[2];
  norm2 = REAL_NORMSQR3(dist[0], dist[1], dist[2]);
  rnorm = REAL_RSQRT(norm2);

  norm = k * norm2 * rnorm;
  rnorm = (rnorm * rnorm) * (rnorm * rnorm) * rnorm;
  s = sin(norm);
  c = cos(norm);

  dot = REAL_DOT3(dist, nx) * REAL_DOT3(dist, ny);
  dotxy = REAL_DOT3(nx, ny);

  hr = (norm * norm - 3.0) * dot + norm2 * dotxy;
  hi = 3.0 * norm * dot - norm * norm2 * dotxy;

  res = rnorm * (c * hr - s * hi + (c * hi + s * hr) * I);

  return res;
}

static void
fill_slp_cc_near_helmholtzbem3d(const uint * ridx,
				const uint * cidx, pcbem3d bem, bool ntrans,
				pamatrix N)
{
  assemble_cc_near_bem3d(ridx, cidx, bem, ntrans, N,
			 (kernel_func3d) slp_kernel_helmholtzbem3d);
}

static void
fill_slp_cc_far_helmholtzbem3d(const uint * ridx, const uint * cidx,
			       pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_cc_far_bem3d(ridx, cidx, bem, ntrans, N,
			(kernel_func3d) slp_kernel_helmholtzbem3d);
}

static void
fill_dlp_cc_near_helmholtzbem3d(const uint * ridx,
				const uint * cidx, pcbem3d bem, bool ntrans,
				pamatrix N)
{
  assemble_cc_near_bem3d(ridx, cidx, bem, ntrans, N,
			 (kernel_func3d) dlp_kernel_helmholtzbem3d);
}

static void
fill_dlp_cc_far_helmholtzbem3d(const uint * ridx, const uint * cidx,
			       pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_cc_far_bem3d(ridx, cidx, bem, ntrans, N,
			(kernel_func3d) dlp_kernel_helmholtzbem3d);
}

static void
fill_dlp_cl_near_helmholtzbem3d(const uint * ridx,
				const uint * cidx, pcbem3d bem, bool ntrans,
				pamatrix N)
{
  assemble_cl_near_bem3d(ridx, cidx, bem, ntrans, N,
			 dlp_kernel_helmholtzbem3d);
}

static void
fill_dlp_cl_far_helmholtzbem3d(const uint * ridx, const uint * cidx,
			       pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_cl_far_bem3d(ridx, cidx, bem, ntrans, N,
			dlp_kernel_helmholtzbem3d);
}

static void
fill_slp_ll_near_helmholtzbem3d(const uint * ridx,
				const uint * cidx, pcbem3d bem, bool ntrans,
				pamatrix N)
{
  assemble_ll_near_bem3d(ridx, cidx, bem, ntrans, N,
			 slp_kernel_helmholtzbem3d);
}

static void
fill_slp_ll_far_helmholtzbem3d(const uint * ridx, const uint * cidx,
			       pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_ll_far_bem3d(ridx, cidx, bem, ntrans, N,
			slp_kernel_helmholtzbem3d);
}

static void
fill_dlp_ll_near_helmholtzbem3d(const uint * ridx,
				const uint * cidx, pcbem3d bem, bool ntrans,
				pamatrix N)
{
  assemble_ll_near_bem3d(ridx, cidx, bem, ntrans, N,
			 dlp_kernel_helmholtzbem3d);
}

static void
fill_dlp_ll_far_helmholtzbem3d(const uint * ridx, const uint * cidx,
			       pcbem3d bem, bool ntrans, pamatrix N)
{
  assemble_ll_far_bem3d(ridx, cidx, bem, ntrans, N,
			dlp_kernel_helmholtzbem3d);
}

static void
fill_kernel_helmholtzbem3d(pcbem3d bem, const real(*X)[3],
			   const real(*Y)[3], pamatrix V)
{
  fill_bem3d(bem, X, Y, NULL, NULL, V, slp_kernel_helmholtzbem3d);
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

static void
fill_kernel_c_helmholtzbem3d(const uint * idx, const real(*Z)[3],
			     pcbem3d bem, pamatrix V)
{
  fill_row_c_bem3d(idx, Z, bem, V, slp_kernel_helmholtzbem3d);
}

static void
fill_kernel_l_helmholtzbem3d(const uint * idx, const real(*Z)[3],
			     pcbem3d bem, pamatrix V)
{
  fill_row_l_bem3d(idx, Z, bem, V, slp_kernel_helmholtzbem3d);
}

static void
fill_dnz_kernel_c_helmholtzbem3d(const uint * idx,
				 const real(*Z)[3], const real(*N)[3],
				 pcbem3d bem, pamatrix V)
{
  fill_dnz_row_c_bem3d(idx, Z, N, bem, V, dlp_kernel_helmholtzbem3d);
}

static void
fill_dnz_kernel_l_helmholtzbem3d(const uint * idx,
				 const real(*Z)[3], const real(*N)[3],
				 pcbem3d bem, pamatrix V)
{
  fill_dnz_row_l_bem3d(idx, Z, N, bem, V, dlp_kernel_helmholtzbem3d);
}

static void
fill_dnzdcol_kernel_c_helmholtzbem3d(const uint * idx,
				     const real(*Z)[3], const real(*N)[3],
				     pcbem3d bem, pamatrix V)
{
  fill_dnz_col_c_bem3d(idx, Z, N, bem, V, hs_kernel_helmholtzbem3d);
}

static void
fill_dnzdcol_kernel_l_helmholtzbem3d(const uint * idx,
				     const real(*Z)[3], const real(*N)[3],
				     pcbem3d bem, pamatrix V)
{
  fill_dnz_col_l_bem3d(idx, Z, N, bem, V, hs_kernel_helmholtzbem3d);
}

static void
fill_dcol_kernel_col_c_helmholtzbem3d(const uint * idx,
				      const real(*Z)[3], pcbem3d bem,
				      pamatrix V)
{
  fill_col_c_bem3d(idx, Z, bem, V, dlp_kernel_helmholtzbem3d);
}

static void
fill_dcol_kernel_col_l_helmholtzbem3d(const uint * idx,
				      const real(*Z)[3], pcbem3d bem,
				      pamatrix V)
{
  fill_col_l_bem3d(idx, Z, bem, V, dlp_kernel_helmholtzbem3d);
}

pbem3d
new_slp_helmholtz_bem3d(field * kvec, pcsurface3d gr, uint q_regular,
			uint q_singular, basisfunctionbem3d basis)
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

  bem->kernel_const = KERNEL_CONST_HELMHOLTZBEM3D;
  bem->kvec = kvec;
  bem->k = NORM3(kvec[0], kvec[1], kvec[2]);
  bem->alpha = 0.0;

  kernels->fundamental = fill_kernel_helmholtzbem3d;
  kernels->dny_fundamental = fill_dny_kernel_helmholtzbem3d;
  kernels->dnx_dny_fundamental = fill_dnx_dny_kernel_helmholtzbem3d;

  if (basis == BASIS_CONSTANT_BEM3D) {
    bem->N_neumann = gr->triangles;

    bem->nearfield = fill_slp_cc_near_helmholtzbem3d;
    bem->nearfield_far = fill_slp_cc_far_helmholtzbem3d;

    kernels->lagrange_row = assemble_bem3d_lagrange_const_amatrix;
    kernels->lagrange_col = assemble_bem3d_lagrange_const_amatrix;

    kernels->fundamental_row = fill_kernel_c_helmholtzbem3d;
    kernels->fundamental_col = fill_kernel_c_helmholtzbem3d;
    kernels->dnz_fundamental_row = fill_dnz_kernel_c_helmholtzbem3d;
    kernels->dnz_fundamental_col = fill_dnz_kernel_c_helmholtzbem3d;

    kernels->kernel_row = fill_kernel_c_helmholtzbem3d;
    kernels->kernel_col = fill_kernel_c_helmholtzbem3d;
    kernels->dnz_kernel_row = fill_dnz_kernel_c_helmholtzbem3d;
    kernels->dnz_kernel_col = fill_dnz_kernel_c_helmholtzbem3d;
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

    bem->nearfield = fill_slp_ll_near_helmholtzbem3d;
    bem->nearfield_far = fill_slp_ll_far_helmholtzbem3d;

    kernels->lagrange_row = assemble_bem3d_lagrange_linear_amatrix;
    kernels->lagrange_col = assemble_bem3d_lagrange_linear_amatrix;

    kernels->fundamental_row = fill_kernel_l_helmholtzbem3d;
    kernels->fundamental_col = fill_kernel_l_helmholtzbem3d;
    kernels->dnz_fundamental_row = fill_dnz_kernel_l_helmholtzbem3d;
    kernels->dnz_fundamental_col = fill_dnz_kernel_l_helmholtzbem3d;

    kernels->kernel_row = fill_kernel_l_helmholtzbem3d;
    kernels->kernel_col = fill_kernel_l_helmholtzbem3d;
    kernels->dnz_kernel_row = fill_dnz_kernel_l_helmholtzbem3d;
    kernels->dnz_kernel_col = fill_dnz_kernel_l_helmholtzbem3d;
  }

  return bem;
}

pbem3d
new_dlp_helmholtz_bem3d(field * kvec, pcsurface3d gr, uint q_regular,
			uint q_singular, basisfunctionbem3d basis_neumann,
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

  bem->kernel_const = KERNEL_CONST_HELMHOLTZBEM3D;
  bem->kvec = kvec;
  bem->k = NORM3(kvec[0], kvec[1], kvec[2]);
  bem->alpha = alpha;

  kernels->fundamental = fill_kernel_helmholtzbem3d;
  kernels->dny_fundamental = fill_dny_kernel_helmholtzbem3d;
  kernels->dnx_dny_fundamental = fill_dnx_dny_kernel_helmholtzbem3d;

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
    bem->nearfield = fill_dlp_cc_near_helmholtzbem3d;
    bem->nearfield_far = fill_dlp_cc_far_helmholtzbem3d;

    kernels->lagrange_row = assemble_bem3d_lagrange_const_amatrix;
    kernels->lagrange_col = assemble_bem3d_dn_lagrange_const_amatrix;

    kernels->fundamental_row = fill_kernel_c_helmholtzbem3d;
    kernels->fundamental_col = fill_kernel_c_helmholtzbem3d;
    kernels->dnz_fundamental_row = fill_dnz_kernel_c_helmholtzbem3d;
    kernels->dnz_fundamental_col = fill_dnz_kernel_c_helmholtzbem3d;

    kernels->kernel_row = fill_kernel_c_helmholtzbem3d;
    kernels->kernel_col = fill_dcol_kernel_col_c_helmholtzbem3d;
    kernels->dnz_kernel_row = fill_dnz_kernel_c_helmholtzbem3d;
    kernels->dnz_kernel_col = fill_dnzdcol_kernel_c_helmholtzbem3d;
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

    bem->nearfield = fill_dlp_cl_near_helmholtzbem3d;
    bem->nearfield_far = fill_dlp_cl_far_helmholtzbem3d;

    kernels->lagrange_row = assemble_bem3d_lagrange_const_amatrix;
    kernels->lagrange_col = assemble_bem3d_dn_lagrange_linear_amatrix;

    kernels->fundamental_row = fill_kernel_c_helmholtzbem3d;
    kernels->fundamental_col = fill_kernel_l_helmholtzbem3d;
    kernels->dnz_fundamental_row = fill_dnz_kernel_c_helmholtzbem3d;
    kernels->dnz_fundamental_col = fill_dnz_kernel_l_helmholtzbem3d;

    kernels->kernel_row = fill_kernel_c_helmholtzbem3d;
    kernels->kernel_col = fill_dcol_kernel_col_l_helmholtzbem3d;
    kernels->dnz_kernel_row = fill_dnz_kernel_c_helmholtzbem3d;
    kernels->dnz_kernel_col = fill_dnzdcol_kernel_l_helmholtzbem3d;

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

    bem->nearfield = fill_dlp_ll_near_helmholtzbem3d;
    bem->nearfield_far = fill_dlp_ll_far_helmholtzbem3d;

    kernels->lagrange_row = assemble_bem3d_lagrange_linear_amatrix;
    kernels->lagrange_col = assemble_bem3d_dn_lagrange_linear_amatrix;

    kernels->fundamental_row = fill_kernel_l_helmholtzbem3d;
    kernels->fundamental_col = fill_kernel_l_helmholtzbem3d;
    kernels->dnz_fundamental_row = fill_dnz_kernel_l_helmholtzbem3d;
    kernels->dnz_fundamental_col = fill_dnz_kernel_l_helmholtzbem3d;

    kernels->kernel_row = fill_kernel_l_helmholtzbem3d;
    kernels->kernel_col = fill_dcol_kernel_col_l_helmholtzbem3d;
    kernels->dnz_kernel_row = NULL;
    kernels->dnz_kernel_col = fill_dnzdcol_kernel_l_helmholtzbem3d;
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
del_helmholtz_bem3d(pbem3d bem)
{
  del_bem3d(bem);
}

field
rhs_dirichlet_point_helmholtzbem3d(const real * x, const real * n,
				   const void *data)
{

  helmholtz_data *hdata = (helmholtz_data *) data;
  field    *kvec = hdata->kvec;
  field     k =
    REAL_SQRT(ABSSQR(kvec[0]) + ABSSQR(kvec[1]) + ABSSQR(kvec[2]));
  real     *s = hdata->source;

  real      r =
    REAL_SQRT(REAL_SQR(x[0] - s[0]) + REAL_SQR(x[1] - s[1]) +
	      REAL_SQR(x[2] - s[2]));

  (void) n;

  return cexp(I * k * r) / r;
}

field
rhs_neumann_point_helmholtzbem3d(const real * x, const real * n,
				 const void *data)
{
  helmholtz_data *hdata = (helmholtz_data *) data;
  field    *kvec = hdata->kvec;
  field     k =
    REAL_SQRT(ABSSQR(kvec[0]) + ABSSQR(kvec[1]) + ABSSQR(kvec[2]));
  real     *s = hdata->source;

  real      r =
    REAL_SQRT(REAL_SQR(x[0] - s[0]) + REAL_SQR(x[1] - s[1]) +
	      REAL_SQR(x[2] - s[2]));

  field     res;

  res = cexp(I * k * r) / (r * r * r) * (I * k * r - 1.0);
  res = res
    * (n[0] * (x[0] - s[0]) + n[1] * (x[1] - s[1]) + n[2] * (x[2] - s[2]));

  return res;
}

#endif
