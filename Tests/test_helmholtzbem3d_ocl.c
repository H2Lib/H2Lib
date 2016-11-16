#ifdef USE_COMPLEX
#ifdef USE_OPENCL
#ifdef USE_OPENMP

#include "basic.h"
#include "krylov.h"
#include "helmholtzoclbem3d.h"
#include "matrixnorms.h"

#ifdef USE_CAIRO
#include <cairo.h>
#endif

static uint problems = 0;

#define IS_IN_RANGE(a, b, c) (((a) <= (b)) && ((b) <= (c)))

struct _eval_A {
  matrixtype Vtype;
  void     *V;
  matrixtype KMtype;
  void     *KM;
  field     eta;
};

void
addeval_A(field alpha, void *matrix, pcavector x, pavector y)
{
  struct _eval_A *eval = (struct _eval_A *) matrix;

  field     beta;

  beta = -I * alpha * eval->eta;

  switch (eval->KMtype) {
  case AMATRIX:
    addeval_amatrix_avector(alpha, (pamatrix) eval->KM, x, y);
    break;
  case HMATRIX:
    addeval_hmatrix_avector(alpha, (phmatrix) eval->KM, x, y);
    break;
  case H2MATRIX:
    addeval_h2matrix_avector(alpha, (ph2matrix) eval->KM, x, y);
    break;
  default:
    printf("ERROR: unknown matrix type!\n");
    abort();
    break;
  }

  switch (eval->Vtype) {
  case AMATRIX:
    addeval_amatrix_avector(beta, (pamatrix) eval->V, x, y);
    break;
  case HMATRIX:
    addeval_hmatrix_avector(beta, (phmatrix) eval->V, x, y);
    break;
  case H2MATRIX:
    addeval_h2matrix_avector(beta, (ph2matrix) eval->V, x, y);
    break;
  default:
    printf("ERROR: unknown matrix type!\n");
    abort();
    break;
  }
}

field
eval_brakhage_werner_c(pcbem3d bem, pcavector w, field eta, real * x)
{
  pcsurface3d gr = bem->gr;
  const     real(*gr_x)[3] = (const real(*)[3]) gr->x;
  const     uint(*gr_t)[3] = (const uint(*)[3]) gr->t;
  const     real(*gr_n)[3] = (const real(*)[3]) gr->n;
  const preal gr_g = (const preal) gr->g;
  const uint rows = gr->triangles;

  uint      nq = bem->sq->n_single;
  uint      vnq = ROUNDUP(nq, VREAL);
  real     *xx = bem->sq->x_single;
  real     *yy = bem->sq->y_single;
  real     *ww = bem->sq->w_single + 3 * vnq;

  const real *A, *B, *C, *ns;
  uint      s, ss, q;
  real      gs_fac, dx, dy, dz, tx, sx, Ax, Bx, Cx, norm, rnorm, norm2;
  field     k, sum;

  field     res;

  k = bem->k;

  res = 0.0;

  for (s = 0; s < rows; ++s) {
    ss = s;
    gs_fac = gr_g[ss] * bem->kernel_const;
    ns = gr_n[ss];
    A = gr_x[gr_t[ss][0]];
    B = gr_x[gr_t[ss][1]];
    C = gr_x[gr_t[ss][2]];

    sum = 0.0;

    for (q = 0; q < nq; ++q) {
      tx = xx[q];
      sx = yy[q];
      Ax = 1.0 - tx;
      Bx = tx - sx;
      Cx = sx;

      dx = x[0] - (A[0] * Ax + B[0] * Bx + C[0] * Cx);
      dy = x[1] - (A[1] * Ax + B[1] * Bx + C[1] * Cx);
      dz = x[2] - (A[2] * Ax + B[2] * Bx + C[2] * Cx);

      norm2 = dx * dx + dy * dy + dz * dz;
      rnorm = REAL_RSQRT(norm2);
      norm = norm2 * rnorm;
      rnorm *= rnorm * rnorm;

      sum += ww[q]
	* cexp(I * k * norm)
	* rnorm
	* ((1.0 - I * k * norm) * (dx * ns[0] + dy * ns[1] + dz * ns[2]) - I
	   * eta * norm2);
    }

    res += sum * gs_fac * w->v[ss];
  }

  return res;
}

real
max_rel_outer_error(pcbem3d bem, helmholtz_data * hdata, pcavector x,
		    boundary_func3d rhs)
{
  uint      nx, nz, npoints;
  real(*xdata)[3];
  field    *ydata;
  uint      i, j;
  real      error, maxerror;
  real      eta_bw = ABS(bem->k);

  nx = 20;
  nz = 20;
  npoints = nx * nz;

  xdata = (real(*)[3]) allocreal(3 * npoints);
  npoints = 0;
  for (j = 0; j < nz; ++j) {
    for (i = 0; i < nx; ++i) {
      xdata[npoints][0] = -10.0 + (20.0 / (nx - 1)) * i;
      xdata[npoints][1] = 0.0;
      xdata[npoints][2] = -10.0 + (20.0 / (nz - 1)) * j;
      if (REAL_SQR(xdata[npoints][0]) + REAL_SQR(xdata[npoints][2]) > 1) {
	npoints++;
      }
    }
  }

  ydata = allocfield(npoints);

#pragma omp parallel for
  for (j = 0; j < npoints; ++j) {
    ydata[j] = eval_brakhage_werner_c(bem, x, eta_bw, xdata[j]);
  }

  j = 0;
  maxerror =
    ABS(ydata[j] -
	rhs(xdata[j], NULL, hdata)) / ABS(rhs(xdata[j], NULL, hdata));
  for (j = 1; j < npoints; ++j) {
    error =
      ABS(ydata[j] -
	  rhs(xdata[j], NULL, hdata)) / ABS(rhs(xdata[j], NULL, hdata));

    maxerror = error > maxerror ? error : maxerror;
  }

  freemem(ydata);
  freemem(xdata);

  return maxerror;
}

real
max_rel_inner_error(pcbem3d bem, helmholtz_data * hdata, pcavector x,
		    boundary_func3d rhs)
{
  uint      nx, nz, npoints;
  real(*xdata)[3];
  field    *ydata;
  uint      i, j;
  real      error, maxerror;
  real      eta_bw = ABS(bem->k);

  nx = 20;
  nz = 20;
  npoints = nx * nz;

  xdata = (real(*)[3]) allocreal(3 * npoints);
  npoints = 0;
  for (j = 0; j < nz; ++j) {
    for (i = 0; i < nx; ++i) {
      xdata[npoints][0] = -1.0 + (2.0 / (nx - 1)) * i;
      xdata[npoints][1] = 0.0;
      xdata[npoints][2] = -1.0 + (2.0 / (nz - 1)) * j;
      if (REAL_SQR(xdata[npoints][0]) + REAL_SQR(xdata[npoints][2]) < 1) {
	npoints++;
      }
    }
  }

  ydata = allocfield(npoints);

#pragma omp parallel for
  for (j = 0; j < npoints; ++j) {
    ydata[j] = eval_brakhage_werner_c(bem, x, eta_bw, xdata[j]);
  }

  j = 0;
  maxerror =
    ABS(ydata[j] -
	rhs(xdata[j], NULL, hdata)) / ABS(rhs(xdata[j], NULL, hdata));
  for (j = 1; j < npoints; ++j) {
    error =
      ABS(ydata[j] -
	  rhs(xdata[j], NULL, hdata)) / ABS(rhs(xdata[j], NULL, hdata));

    maxerror = error > maxerror ? error : maxerror;
  }

  freemem(ydata);
  freemem(xdata);

  return maxerror;
}

static void
test_hmatrix_system(const char *apprxtype, pcamatrix Vfull,
		    pcamatrix KMfull, pblock block, pbem3d bem_slp,
		    phmatrix V, pbem3d bem_dlp, phmatrix KM, bool linear,
		    bool exterior, real low, real high)
{
  struct _eval_A eval;
  helmholtz_data hdata;
  pavector  x, b;
  real      errorV, errorKM, error_solve, eps_solve;
  uint      steps, iter;
  boundary_func3d rhs = (boundary_func3d) rhs_dirichlet_point_helmholtzbem3d;

  eps_solve = 1.0e-12;
  steps = 1000;

  printf("Testing: %s Hmatrix %s\n"
	 "====================================\n\n",
	 (exterior == true ? "exterior" : "interior"), apprxtype);

  SCHEDULE_OPENCL(0, 1, assemble_bem3d_hmatrix, bem_slp, block, V);
  SCHEDULE_OPENCL(0, 1, assemble_bem3d_hmatrix, bem_dlp, block, KM);

  errorV = norm2diff_amatrix_hmatrix(V, Vfull) / norm2_amatrix(Vfull);
  printf("rel. error V       : %.5e\n", errorV);
  errorKM = norm2diff_amatrix_hmatrix(KM, KMfull) / norm2_amatrix(KMfull);
  printf("rel. error K%c0.5*M : %.5e\n", (exterior == true ? '-' : '+'),
	 errorKM);

  eval.V = V;
  eval.Vtype = HMATRIX;
  eval.KM = KM;
  eval.KMtype = HMATRIX;
  eval.eta = ABS(bem_slp->k);

  hdata.kvec = allocreal(3);
  hdata.kvec[0] = ABS(bem_slp->k);
  hdata.kvec[1] = 0.0;
  hdata.kvec[2] = 0.0;
  hdata.source = allocreal(3);
  if (exterior) {
    hdata.source[0] = 0.0, hdata.source[1] = 0.0, hdata.source[2] = 0.2;
  }
  else {
    hdata.source[0] = 0.0, hdata.source[1] = 0.0, hdata.source[2] = 5.0;
  }

  x = new_avector(Vfull->rows);
  b = new_avector(KMfull->cols);

  printf("Solving Dirichlet problem:\n");

  integrate_bem3d_c_avector(bem_dlp, rhs, b, (void *) &hdata);

  iter = solve_gmres_avector(&eval, addeval_A, b, x, eps_solve, steps, 50);
  printf("GMRES iterations:\n");
  printf("  %d\n", iter);

  if (exterior) {
    error_solve = max_rel_outer_error(bem_slp, &hdata, x, rhs);
  }
  else {
    error_solve = max_rel_inner_error(bem_slp, &hdata, x, rhs);
  }

  printf("max. rel. error : %.5e       %s\n", error_solve,
	 (IS_IN_RANGE(low, error_solve, high) ? "    okay" : "NOT okay"));

  if (!IS_IN_RANGE(low, error_solve, high))
    problems++;

  printf("\n");

  del_avector(x);
  del_avector(b);
  freemem(hdata.source);
  freemem(hdata.kvec);
}

static void
test_h2matrix_system(const char *apprxtype, pcamatrix Vfull,
		     pcamatrix KMfull, pblock block, pbem3d bem_slp,
		     ph2matrix V, pbem3d bem_dlp, ph2matrix KM, bool linear,
		     bool exterior, real low, real high)
{
  struct _eval_A eval;
  helmholtz_data hdata;
  pavector  x, b;
  real      errorV, errorKM, error_solve, eps_solve;
  uint      steps, iter;
  boundary_func3d rhs = (boundary_func3d) rhs_dirichlet_point_helmholtzbem3d;

  eps_solve = 1.0e-12;
  steps = 1000;

  printf("Testing: %s H2matrix %s\n"
	 "====================================\n\n",
	 (exterior == true ? "exterior" : "interior"), apprxtype);

  assemble_bem3d_h2matrix_row_clusterbasis(bem_slp, V->rb);
  assemble_bem3d_h2matrix_col_clusterbasis(bem_slp, V->cb);
  SCHEDULE_OPENCL(0, 1, assemble_bem3d_h2matrix, bem_slp, V);

  assemble_bem3d_h2matrix_row_clusterbasis(bem_dlp, KM->rb);
  assemble_bem3d_h2matrix_col_clusterbasis(bem_dlp, KM->cb);
  SCHEDULE_OPENCL(0, 1, assemble_bem3d_h2matrix, bem_dlp, KM);

  eval.V = V;
  eval.Vtype = H2MATRIX;
  eval.KM = KM;
  eval.KMtype = H2MATRIX;
  eval.eta = ABS(bem_slp->k);

  hdata.kvec = allocreal(3);
  hdata.kvec[0] = ABS(bem_slp->k);
  hdata.kvec[1] = 0.0;
  hdata.kvec[2] = 0.0;
  hdata.source = allocreal(3);
  if (exterior) {
    hdata.source[0] = 0.0, hdata.source[1] = 0.0, hdata.source[2] = 0.2;
  }
  else {
    hdata.source[0] = 0.0, hdata.source[1] = 0.0, hdata.source[2] = 5.0;
  }

  errorV = norm2diff_amatrix_h2matrix(V, Vfull) / norm2_amatrix(Vfull);
  printf("rel. error V       : %.5e\n", errorV);
  errorKM = norm2diff_amatrix_h2matrix(KM, KMfull) / norm2_amatrix(KMfull);
  printf("rel. error K%c0.5*M : %.5e\n", (exterior == true ? '-' : '+'),
	 errorKM);

  x = new_avector(Vfull->rows);
  b = new_avector(KMfull->cols);

  printf("Solving Dirichlet problem:\n");

  integrate_bem3d_c_avector(bem_dlp, rhs, b, (void *) &hdata);

  iter = solve_gmres_avector(&eval, addeval_A, b, x, eps_solve, steps, 50);
  printf("GMRES iterations:\n");
  printf("  %d\n", iter);

  error_solve = max_rel_outer_error(bem_slp, &hdata, x, rhs);

  printf("max. rel. error : %.5e       %s\n", error_solve,
	 (IS_IN_RANGE(low, error_solve, high) ? "    okay" : "NOT okay"));

  if (!IS_IN_RANGE(low, error_solve, high))
    problems++;

  printf("\n");

  del_avector(x);
  del_avector(b);
  freemem(hdata.source);
}

int
main(int argc, char **argv)
{
  pmacrosurface3d mg;
  psurface3d gr;
  pamatrix  Vfull, KMfull;
  pbem3d    bem_slp, bem_dlp;
  pcluster  root;
  pblock block;
  phmatrix  V, KM;
  pclusterbasis Vrb, Vcb, KMrb, KMcb;
  ph2matrix V2, KM2;
  uint      n, q, clf, m, l;
  real      eta, delta, eps_aca;
  field     k;
  cl_device_id *devices;
  cl_uint   ndevices;
  uint     *idx;
  uint      i;

  init_h2lib(&argc, &argv);

  get_opencl_devices(&devices, &ndevices);
  ndevices = 1;
  set_opencl_devices(devices, ndevices, 2);

  k = 2.0;
  n = 512;
  q = 2;
  clf = 16;
  eta = 1.0;

  mg = new_sphere_macrosurface3d();
  gr = build_from_macrosurface3d_surface3d(mg, REAL_SQRT(n * 0.125));
  n = gr->triangles;

  printf("Testing unit sphere with %d triangles\n", n);

  bem_slp = new_slp_helmholtz_ocl_bem3d(k, gr, q, q + 2, BASIS_CONSTANT_BEM3D,
					BASIS_CONSTANT_BEM3D);
  bem_dlp = new_dlp_helmholtz_ocl_bem3d(k, gr, q, q + 2, BASIS_CONSTANT_BEM3D,
					BASIS_CONSTANT_BEM3D, 0.5);
  root = build_bem3d_cluster(bem_slp, clf, BASIS_CONSTANT_BEM3D);
  block = build_nonstrict_block(root, root, &eta, admissible_max_cluster);

  max_pardepth = 0;

  Vfull = new_amatrix(n, n);
  KMfull = new_amatrix(n, n);
  idx = allocuint(n);
  for (i = 0; i < n; ++i) {
    idx[i] = i;
  }
  SCHEDULE_OPENCL(0, 1, bem_slp->nearfield, idx, idx, bem_slp, false, Vfull);
  SCHEDULE_OPENCL(0, 1, bem_dlp->nearfield, idx, idx, bem_dlp, false, KMfull);

  V = build_from_block_hmatrix(block, 0);
  KM = build_from_block_hmatrix(block, 0);

  printf("----------------------------------------\n");
  printf("Testing outer Boundary integral equations:\n");
  printf("----------------------------------------\n\n");

  /*
   * Test Interpolation
   */

  m = 4;

  setup_hmatrix_aprx_inter_row_bem3d(bem_slp, root, root, block, m);
  setup_hmatrix_aprx_inter_row_bem3d(bem_dlp, root, root, block, m);
  test_hmatrix_system("Interpolation row", Vfull, KMfull, block, bem_slp, V,
		      bem_dlp, KM, false, true, 1.0e-3, 2.0e-3);

  setup_hmatrix_aprx_inter_col_bem3d(bem_slp, root, root, block, m);
  setup_hmatrix_aprx_inter_col_bem3d(bem_dlp, root, root, block, m);
  test_hmatrix_system("Interpolation column", Vfull, KMfull, block, bem_slp,
		      V, bem_dlp, KM, false, true, 1.0e-3, 2.0e-3);

  setup_hmatrix_aprx_inter_mixed_bem3d(bem_slp, root, root, block, m);
  setup_hmatrix_aprx_inter_mixed_bem3d(bem_dlp, root, root, block, m);
  test_hmatrix_system("Interpolation mixed", Vfull, KMfull, block, bem_slp, V,
		      bem_dlp, KM, false, true, 1.0e-3, 2.0e-3);

  /*
   * Test Green
   */

  m = 5;
  l = 1;
  delta = 0.5;

  setup_hmatrix_aprx_green_row_bem3d(bem_slp, root, root, block, m, l, delta,
				     build_bem3d_cube_quadpoints);
  setup_hmatrix_aprx_green_row_bem3d(bem_dlp, root, root, block, m, l, delta,
				     build_bem3d_cube_quadpoints);
  test_hmatrix_system("Green row", Vfull, KMfull, block, bem_slp, V, bem_dlp,
		      KM, false, true, 1.0e-3, 2.0e-3);

  setup_hmatrix_aprx_green_col_bem3d(bem_slp, root, root, block, m, l, delta,
				     build_bem3d_cube_quadpoints);
  setup_hmatrix_aprx_green_col_bem3d(bem_dlp, root, root, block, m, l, delta,
				     build_bem3d_cube_quadpoints);
  test_hmatrix_system("Green column", Vfull, KMfull, block, bem_slp, V,
		      bem_dlp, KM, false, true, 1.0e-3, 2.0e-3);

  setup_hmatrix_aprx_green_mixed_bem3d(bem_slp, root, root, block, m, l,
				       delta, build_bem3d_cube_quadpoints);
  setup_hmatrix_aprx_green_mixed_bem3d(bem_dlp, root, root, block, m, l,
				       delta, build_bem3d_cube_quadpoints);
  test_hmatrix_system("Green mixed", Vfull, KMfull, block, bem_slp, V,
		      bem_dlp, KM, false, true, 1.0e-3, 2.0e-3);

  /*
   * Test Greenhybrid
   */

  m = 2;
  l = 1;
  delta = 1.0;
  eps_aca = 2.0e-2;

  setup_hmatrix_aprx_greenhybrid_row_bem3d(bem_slp, root, root, block, m, l,
					   delta, eps_aca,
					   build_bem3d_cube_quadpoints);
  setup_hmatrix_aprx_greenhybrid_row_bem3d(bem_dlp, root, root, block, m, l,
					   delta, eps_aca,
					   build_bem3d_cube_quadpoints);
  test_hmatrix_system("Greenhybrid row", Vfull, KMfull, block, bem_slp, V,
		      bem_dlp, KM, false, true, 1.0e-3, 2.0e-3);

  setup_hmatrix_aprx_greenhybrid_col_bem3d(bem_slp, root, root, block, m, l,
					   delta, eps_aca,
					   build_bem3d_cube_quadpoints);
  setup_hmatrix_aprx_greenhybrid_col_bem3d(bem_dlp, root, root, block, m, l,
					   delta, eps_aca,
					   build_bem3d_cube_quadpoints);
  test_hmatrix_system("Greenhybrid column", Vfull, KMfull, block, bem_slp, V,
		      bem_dlp, KM, false, true, 1.0e-3, 2.0e-3);

  setup_hmatrix_aprx_greenhybrid_mixed_bem3d(bem_slp, root, root, block, m, l,
					     delta, eps_aca,
					     build_bem3d_cube_quadpoints);
  setup_hmatrix_aprx_greenhybrid_mixed_bem3d(bem_dlp, root, root, block, m, l,
					     delta, eps_aca,
					     build_bem3d_cube_quadpoints);
  test_hmatrix_system("Greenhybrid mixed", Vfull, KMfull, block, bem_slp, V,
		      bem_dlp, KM, false, true, 1.0e-3, 2.0e-3);

  /*
   * Test ACA / PACA / HCA
   */

  m = 3;
  eps_aca = 1.0e-2;

  /* Nearfield computation on GPU not applicable here yet! */

//  setup_hmatrix_aprx_aca_bem3d(bem_slp, root, root, block, eps_aca);
//  setup_hmatrix_aprx_aca_bem3d(bem_dlp, root, root, block, eps_aca);
//  test_hmatrix_system("ACA full pivoting", Vfull, KMfull, block, bem_slp, V,
//      bem_dlp, KM, false, true, 1.0e-3, 2.0e-3);
//
//  setup_hmatrix_aprx_paca_bem3d(bem_slp, root, root, block, eps_aca);
//  setup_hmatrix_aprx_paca_bem3d(bem_dlp, root, root, block, eps_aca);
//  test_hmatrix_system("ACA partial pivoting", Vfull, KMfull, block, bem_slp, V,
//      bem_dlp, KM, false, true, 1.0e-3, 2.0e-3);
  setup_hmatrix_aprx_hca_bem3d(bem_slp, root, root, block, m, eps_aca);
  setup_hmatrix_aprx_hca_bem3d(bem_dlp, root, root, block, m, eps_aca);
  test_hmatrix_system("HCA2", Vfull, KMfull, block, bem_slp, V, bem_dlp, KM,
		      false, true, 1.0e-3, 2.0e-3);

  del_hmatrix(V);
  del_hmatrix(KM);
  del_block(block);

  /*
   * H2-matrix
   */

  block = build_strict_block(root, root, &eta, admissible_max_cluster);

  Vrb = build_from_cluster_clusterbasis(root);
  Vcb = build_from_cluster_clusterbasis(root);
  KMrb = build_from_cluster_clusterbasis(root);
  KMcb = build_from_cluster_clusterbasis(root);

  V2 = build_from_block_h2matrix(block, Vrb, Vcb);
  KM2 = build_from_block_h2matrix(block, KMrb, KMcb);

  /*
   * Test Interpolation
   */

  m = 3;

  setup_h2matrix_aprx_inter_bem3d(bem_slp, Vrb, Vcb, block, m);
  setup_h2matrix_aprx_inter_bem3d(bem_dlp, KMrb, KMcb, block, m);
  test_h2matrix_system("Interpolation", Vfull, KMfull, block, bem_slp, V2,
		       bem_dlp, KM2, false, true, 1.0e-3, 2.0e-3);

  /*
   * Test Greenhybrid
   */

  m = 2;
  l = 1;
  delta = 1.0;
  eps_aca = 2.0e-2;

  setup_h2matrix_aprx_greenhybrid_bem3d(bem_slp, Vrb, Vcb, block, m, l, delta,
					eps_aca, build_bem3d_cube_quadpoints);
  setup_h2matrix_aprx_greenhybrid_bem3d(bem_dlp, KMrb, KMcb, block, m, l,
					delta, eps_aca,
					build_bem3d_cube_quadpoints);
  test_h2matrix_system("Greenhybrid", Vfull, KMfull, block, bem_slp, V2,
		       bem_dlp, KM2, false, true, 1.0e-3, 2.0e-3);

  /* Nearfield computation on GPU not applicable here yet! */

//  setup_h2matrix_aprx_greenhybrid_ortho_bem3d(bem_slp, Vrb, Vcb, block, m, l,
//      delta, eps_aca, build_bem3d_cube_quadpoints);
//  setup_h2matrix_aprx_greenhybrid_ortho_bem3d(bem_dlp, KMrb, KMcb, block, m, l,
//      delta, eps_aca, build_bem3d_cube_quadpoints);
//  test_h2matrix_system("Greenhybrid ortho", Vfull, KMfull, block, bem_slp, V2,
//      bem_dlp, KM2, false, true, 1.0e-3, 2.0e-3);
  del_h2matrix(V2);
  del_h2matrix(KM2);
  del_block(block);
  freemem(root->idx);
  del_cluster(root);
  del_helmholtz_ocl_bem3d(bem_slp);
  del_helmholtz_ocl_bem3d(bem_dlp);

  del_amatrix(Vfull);
  del_amatrix(KMfull);
  del_surface3d(gr);
  del_macrosurface3d(mg);
  freemem(idx);

  (void) printf("----------------------------------------\n"
		"  %u matrices and\n"
		"  %u vectors still active\n"
		"  %u errors found\n", getactives_amatrix(),
		getactives_avector(), problems);

  uninit_h2lib();

  return problems;
}

#else
#include <stdio.h>

int
main(int argc, char **argv)
{
  (void) argc;
  (void) argv;

  fprintf(stderr,
	  "Need to activate OpenMP by \'-DUSE_OPENMP\' to use this test\n");

  return 0;
}
#endif

#else
#include <stdio.h>

int
main(int argc, char **argv)
{
  (void) argc;
  (void) argv;

  fprintf(stderr,
	  "Need to activate OpenCL by \'-DUSE_OPENCL\' to use this test\n");

  return 0;
}
#endif

#else
#include <stdio.h>

int
main(int argc, char **argv)
{
  (void) argc;
  (void) argv;

  fprintf(stderr,
	  "Need to activate complex numbers by \'-DUSE_COMPLEX\' to use this test\n");

  return 0;
}
#endif
