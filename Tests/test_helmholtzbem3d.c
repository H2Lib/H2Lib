#ifdef USE_COMPLEX

#include "basic.h"
#include "krylov.h"
#include "helmholtzbem3d.h"
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

/* eval Brakhage-Werner system matrix */
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

field
eval_brakhage_werner_l(pcbem3d bem, pcavector w, field eta, real * x)
{
  pcsurface3d gr = bem->gr;
  const     real(*gr_x)[3] = (const real(*)[3]) gr->x;
  const     uint(*gr_t)[3] = (const uint(*)[3]) gr->t;
  const     real(*gr_n)[3] = (const real(*)[3]) gr->n;
  const preal gr_g = (const preal) gr->g;
  const uint triangles = gr->triangles;
  const uint vertices = gr->vertices;
  real      k = bem->k;

  uint      nq = bem->sq->n_single;
  uint      vnq = ROUNDUP(nq, VREAL);
  real     *xx = bem->sq->x_single;
  real     *yy = bem->sq->y_single;
  real     *ww = bem->sq->w_single;
  real      base = bem->sq->base_single;

  pavector  v;
  const real *A, *B, *C, *N;
  uint      s, i, ii, q;
  real      gs_fac, tx, sx, Ax, Bx, Cx, norm, rnorm, norm2, dx[3];
  field     sum;
  field    *quad;

  field     res;

  assert(gr->vertices == w->dim);

  quad = allocfield(vnq);
  v = new_zero_avector(vertices);

  for (s = 0; s < triangles; s++) {
    gs_fac = gr_g[s] * bem->kernel_const;
    A = gr_x[gr_t[s][0]];
    B = gr_x[gr_t[s][1]];
    C = gr_x[gr_t[s][2]];
    N = gr_n[s];

    for (q = 0; q < nq; ++q) {
      tx = xx[q];
      sx = yy[q];
      Ax = 1.0 - tx;
      Bx = tx - sx;
      Cx = sx;

      dx[0] = x[0] - (A[0] * Ax + B[0] * Bx + C[0] * Cx);
      dx[1] = x[1] - (A[1] * Ax + B[1] * Bx + C[1] * Cx);
      dx[2] = x[2] - (A[2] * Ax + B[2] * Bx + C[2] * Cx);

      norm2 = REAL_NORMSQR3(dx[0], dx[1], dx[2]);
      rnorm = REAL_RSQRT(norm2);
      norm = norm2 * rnorm;
      rnorm *= rnorm * rnorm;

      quad[q] =
	cexp(I * k * norm) * rnorm
	* ((1.0 - I * k * norm) *
	   (dx[0] * N[0] + dx[1] * N[1] + dx[2] * N[2]) - I * eta * norm2);
    }

    ww = bem->sq->w_single;

    for (i = 0; i < 3; ++i) {
      ii = gr_t[s][i];
      sum = base;

      for (q = 0; q < nq; ++q) {
	sum += ww[q] * quad[q];
      }

      assert(ii < vertices);
      v->v[ii] += sum * gs_fac;

      ww += vnq;
    }
  }

  res = 0.0;
  for (i = 0; i < vertices; ++i) {
    res += w->v[i] * v->v[i];
  }

  del_avector(v);
  freemem(quad);

  return res;
}

real
max_rel_outer_error(pcbem3d bem, helmholtz_data * hdata, pcavector x,
		    boundary_func3d rhs, basisfunctionbem3d basis)
{
  uint      nx, nz, npoints;
  real(*xdata)[3];
  field    *ydata;
  uint      i, j;
  real      error, maxerror;
  real      eta_bw = bem->k;

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

  if (basis == BASIS_CONSTANT_BEM3D) {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for (j = 0; j < npoints; ++j) {
      ydata[j] = eval_brakhage_werner_c(bem, x, eta_bw, xdata[j]);
    }
  }
  else {
    assert(basis == BASIS_LINEAR_BEM3D);
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for (j = 0; j < npoints; ++j) {
      ydata[j] = eval_brakhage_werner_l(bem, x, eta_bw, xdata[j]);
    }
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
		    boundary_func3d rhs, basisfunctionbem3d basis)
{
  uint      nx, nz, npoints;
  real(*xdata)[3];
  field    *ydata;
  uint      i, j;
  real      error, maxerror;
  real      eta_bw = bem->k;

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

  if (basis == BASIS_CONSTANT_BEM3D) {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for (j = 0; j < npoints; ++j) {
      ydata[j] = eval_brakhage_werner_c(bem, x, eta_bw, xdata[j]);
    }
  }
  else {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for (j = 0; j < npoints; ++j) {
      ydata[j] = eval_brakhage_werner_l(bem, x, eta_bw, xdata[j]);
    }
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
test_system(matrixtype mattype, const char *apprxtype,
	    pcamatrix Vfull, pcamatrix KMfull, pblock brootV, pbem3d bem_slp,
	    void *V, pblock brootKM, pbem3d bem_dlp, void *KM,
	    basisfunctionbem3d basis_neumann,
	    basisfunctionbem3d basis_dirichlet, bool exterior, real low,
	    real high)
{
  pavector  x, b;
  struct _eval_A eval;
  helmholtz_data hdata;
  real      errorV, errorKM, error_solve, eps_solve;
  uint      steps, iter;
  boundary_func3d rhs = (boundary_func3d) rhs_dirichlet_point_helmholtzbem3d;

  eps_solve = 1.0e-12;
  steps = 1000;

  printf("Testing: %c%c %s %s %s\n"
	 "====================================\n\n",
	 basis_neumann == BASIS_LINEAR_BEM3D ? 'l' : 'c',
	 basis_dirichlet == BASIS_LINEAR_BEM3D ? 'l' : 'c',
	 (exterior == true ? "exterior" : "interior"),
	 (mattype == HMATRIX ? "Hmatrix" : "H2matrix"), apprxtype);

  eval.V = V;
  eval.KM = KM;
  eval.eta = bem_slp->k;

  if (mattype == HMATRIX) {
    assemble_bem3d_hmatrix(bem_slp, brootV, (phmatrix) V);
    assemble_bem3d_hmatrix(bem_dlp, brootKM, KM);
    errorV = norm2diff_amatrix_hmatrix((phmatrix) V, Vfull)
      / norm2_amatrix(Vfull);
    printf("rel. error V       : %.5e\n", errorV);
    errorKM = norm2diff_amatrix_hmatrix(KM, KMfull) / norm2_amatrix(KMfull);
    printf("rel. error K%c0.5*M : %.5e\n", (exterior == true ? '+' : '-'),
	   errorKM);

    eval.Vtype = HMATRIX;
    eval.KMtype = HMATRIX;
  }
  else {
    assert(mattype == H2MATRIX);
    assemble_bem3d_h2matrix_row_clusterbasis(bem_slp, ((ph2matrix) V)->rb);
    assemble_bem3d_h2matrix_col_clusterbasis(bem_slp, ((ph2matrix) V)->cb);
    assemble_bem3d_h2matrix(bem_slp, (ph2matrix) V);

    assemble_bem3d_h2matrix_row_clusterbasis(bem_dlp, ((ph2matrix) KM)->rb);
    assemble_bem3d_h2matrix_col_clusterbasis(bem_dlp, ((ph2matrix) KM)->cb);
    assemble_bem3d_h2matrix(bem_dlp, (ph2matrix) KM);

    errorV = norm2diff_amatrix_h2matrix((ph2matrix) V, Vfull)
      / norm2_amatrix(Vfull);
    printf("rel. error V       : %.5e\n", errorV);
    errorKM = norm2diff_amatrix_h2matrix((ph2matrix) KM, KMfull)
      / norm2_amatrix(KMfull);
    printf("rel. error K%c0.5*M : %.5e\n", (exterior == true ? '-' : '+'),
	   errorKM);

    eval.Vtype = H2MATRIX;
    eval.KMtype = H2MATRIX;
  }

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

  x = new_avector(Vfull->cols);
  b = new_avector(KMfull->rows);

  printf("Solving Dirichlet problem:\n");

  if (basis_dirichlet == BASIS_LINEAR_BEM3D) {
    integrate_bem3d_l_avector(bem_dlp, rhs, b, (void *) &hdata);
  }
  else {
    integrate_bem3d_c_avector(bem_dlp, rhs, b, (void *) &hdata);
  }

  iter = solve_gmres_avector(&eval, addeval_A, b, x, eps_solve, steps, 50);
  printf("GMRES iterations:\n");
  printf("  %d\n", iter);

  error_solve = max_rel_outer_error(bem_slp, &hdata, x, rhs, basis_neumann);

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

void
test_suite(pcsurface3d gr, field k, uint q, uint clf, real eta,
	   basisfunctionbem3d row_basis, basisfunctionbem3d col_basis,
	   bool exterior, real error_min, real error_max)
{
  pbem3d    bem_slp, bem_dlp;
  pcluster  rootn, rootd;
  pblock    brootV, brootKM;
  pamatrix  Vfull, KMfull;
  phmatrix  V, KM;
  pclusterbasis Vrb, Vcb, KMrb, KMcb;
  ph2matrix V2, KM2;
  uint      nn, nd;

  uint      m;
  uint      l;
  real      delta;
  real      eps_aca;

  nn = row_basis == BASIS_LINEAR_BEM3D ? gr->vertices : gr->triangles;
  nd = col_basis == BASIS_LINEAR_BEM3D ? gr->vertices : gr->triangles;

  bem_slp = new_slp_helmholtz_bem3d(k, gr, q, q + 2, row_basis, row_basis);
  bem_dlp = new_dlp_helmholtz_bem3d(k, gr, q, q + 2, row_basis, col_basis,
				    exterior ? 0.5 : -0.5);

  rootn = build_bem3d_cluster(bem_slp, clf, row_basis);
  rootd = build_bem3d_cluster(bem_dlp, clf, col_basis);

  brootV = build_nonstrict_block(rootn, rootn, &eta, admissible_max_cluster);
  brootKM = build_nonstrict_block(rootn, rootd, &eta, admissible_max_cluster);

  V = build_from_block_hmatrix(brootV, 0);
  KM = build_from_block_hmatrix(brootKM, 0);

  Vfull = new_amatrix(nn, nn);
  KMfull = new_amatrix(nn, nd);
  bem_slp->nearfield(NULL, NULL, bem_slp, false, Vfull);
  bem_dlp->nearfield(NULL, NULL, bem_dlp, false, KMfull);

  /*
   * Test Interpolation
   */

  m = 4;

  setup_hmatrix_aprx_inter_row_bem3d(bem_slp, rootn, rootn, brootV, m);
  setup_hmatrix_aprx_inter_row_bem3d(bem_dlp, rootn, rootd, brootKM, m);
  test_system(HMATRIX, "Interpolation row", Vfull, KMfull, brootV, bem_slp, V,
	      brootKM, bem_dlp, KM, row_basis, col_basis, exterior, error_min,
	      error_max);

  setup_hmatrix_aprx_inter_col_bem3d(bem_slp, rootn, rootn, brootV, m);
  setup_hmatrix_aprx_inter_col_bem3d(bem_dlp, rootn, rootd, brootKM, m);
  test_system(HMATRIX, "Interpolation column", Vfull, KMfull, brootV, bem_slp,
	      V, brootKM, bem_dlp, KM, row_basis, col_basis, exterior,
	      error_min, error_max);

  setup_hmatrix_aprx_inter_mixed_bem3d(bem_slp, rootn, rootn, brootV, m);
  setup_hmatrix_aprx_inter_mixed_bem3d(bem_dlp, rootn, rootd, brootKM, m);
  test_system(HMATRIX, "Interpolation mixed", Vfull, KMfull, brootV, bem_slp,
	      V, brootKM, bem_dlp, KM, row_basis, col_basis, exterior,
	      error_min, error_max);

  /*
   * Test Green
   */

  m = 5;
  l = 1;
  delta = 0.5;

  setup_hmatrix_aprx_green_row_bem3d(bem_slp, rootn, rootn, brootV, m, l,
				     delta, build_bem3d_cube_quadpoints);
  setup_hmatrix_aprx_green_row_bem3d(bem_dlp, rootn, rootd, brootKM, m, l,
				     delta, build_bem3d_cube_quadpoints);
  test_system(HMATRIX, "Green row", Vfull, KMfull, brootV, bem_slp, V,
	      brootKM, bem_dlp, KM, row_basis, col_basis, exterior, error_min,
	      error_max);
  setup_hmatrix_aprx_green_col_bem3d(bem_slp, rootn, rootn, brootV, m, l,
				     delta, build_bem3d_cube_quadpoints);
  setup_hmatrix_aprx_green_col_bem3d(bem_dlp, rootn, rootd, brootKM, m, l,
				     delta, build_bem3d_cube_quadpoints);
  test_system(HMATRIX, "Green column", Vfull, KMfull, brootV, bem_slp, V,
	      brootKM, bem_dlp, KM, row_basis, col_basis, exterior, error_min,
	      error_max);

  setup_hmatrix_aprx_green_mixed_bem3d(bem_slp, rootn, rootn, brootV, m, l,
				       delta, build_bem3d_cube_quadpoints);
  setup_hmatrix_aprx_green_mixed_bem3d(bem_dlp, rootn, rootd, brootKM, m, l,
				       delta, build_bem3d_cube_quadpoints);
  test_system(HMATRIX, "Green mixed", Vfull, KMfull, brootV, bem_slp, V,
	      brootKM, bem_dlp, KM, row_basis, col_basis, exterior, error_min,
	      error_max);

  /*
   * Test Greenhybrid
   */

  m = 2;
  l = 1;
  delta = 1.0;
  eps_aca = 5.0e-3;

  setup_hmatrix_aprx_greenhybrid_row_bem3d(bem_slp, rootn, rootn, brootV, m,
					   l, delta, eps_aca,
					   build_bem3d_cube_quadpoints);
  setup_hmatrix_aprx_greenhybrid_row_bem3d(bem_dlp, rootn, rootd, brootKM, m,
					   l, delta, eps_aca,
					   build_bem3d_cube_quadpoints);
  test_system(HMATRIX, "Greenhybrid row", Vfull, KMfull, brootV, bem_slp, V,
	      brootKM, bem_dlp, KM, row_basis, col_basis, exterior, error_min,
	      error_max);

  setup_hmatrix_aprx_greenhybrid_col_bem3d(bem_slp, rootn, rootn, brootV, m,
					   l, delta, eps_aca,
					   build_bem3d_cube_quadpoints);
  setup_hmatrix_aprx_greenhybrid_col_bem3d(bem_dlp, rootn, rootd, brootKM, m,
					   l, delta, eps_aca,
					   build_bem3d_cube_quadpoints);
  test_system(HMATRIX, "Greenhybrid column", Vfull, KMfull, brootV, bem_slp,
	      V, brootKM, bem_dlp, KM, row_basis, col_basis, exterior,
	      error_min, error_max);

  setup_hmatrix_aprx_greenhybrid_mixed_bem3d(bem_slp, rootn, rootn, brootV, m,
					     l, delta, eps_aca,
					     build_bem3d_cube_quadpoints);
  setup_hmatrix_aprx_greenhybrid_mixed_bem3d(bem_dlp, rootn, rootd, brootKM,
					     m, l, delta, eps_aca,
					     build_bem3d_cube_quadpoints);
  test_system(HMATRIX, "Greenhybrid mixed", Vfull, KMfull, brootV, bem_slp, V,
	      brootKM, bem_dlp, KM, row_basis, col_basis, exterior, error_min,
	      error_max);

  /*
   * Test ACA / PACA / HCA
   */

  m = 3;
  eps_aca = 1.0e-2;

  setup_hmatrix_aprx_aca_bem3d(bem_slp, rootn, rootn, brootV, eps_aca);
  setup_hmatrix_aprx_aca_bem3d(bem_dlp, rootn, rootd, brootKM, eps_aca);
  test_system(HMATRIX, "ACA full pivoting", Vfull, KMfull, brootV, bem_slp, V,
	      brootKM, bem_dlp, KM, row_basis, col_basis, exterior, error_min,
	      error_max);
  setup_hmatrix_aprx_paca_bem3d(bem_slp, rootn, rootn, brootV, eps_aca);
  setup_hmatrix_aprx_paca_bem3d(bem_dlp, rootn, rootd, brootKM, eps_aca);
  test_system(HMATRIX, "ACA partial pivoting", Vfull, KMfull, brootV, bem_slp,
	      V, brootKM, bem_dlp, KM, row_basis, col_basis, exterior,
	      error_min, error_max);

  setup_hmatrix_aprx_hca_bem3d(bem_slp, rootn, rootn, brootV, m, eps_aca);
  setup_hmatrix_aprx_hca_bem3d(bem_dlp, rootn, rootd, brootKM, m, eps_aca);
  test_system(HMATRIX, "HCA2", Vfull, KMfull, brootV, bem_slp, V, brootKM,
	      bem_dlp, KM, row_basis, col_basis, exterior, error_min,
	      error_max);
  /*
   * H2-matrix
   */

  del_hmatrix(V);
  del_hmatrix(KM);
  del_block(brootV);
  del_block(brootKM);

  brootV = build_strict_block(rootn, rootn, &eta, admissible_max_cluster);
  brootKM = build_strict_block(rootn, rootd, &eta, admissible_max_cluster);

  Vrb = build_from_cluster_clusterbasis(rootn);
  Vcb = build_from_cluster_clusterbasis(rootn);
  KMrb = build_from_cluster_clusterbasis(rootn);
  KMcb = build_from_cluster_clusterbasis(rootd);

  V2 = build_from_block_h2matrix(brootV, Vrb, Vcb);
  KM2 = build_from_block_h2matrix(brootKM, KMrb, KMcb);

  /*
   * Test Interpolation
   */

  m = 5;

  setup_h2matrix_aprx_inter_bem3d(bem_slp, Vrb, Vcb, brootV, m);
  setup_h2matrix_aprx_inter_bem3d(bem_dlp, KMrb, KMcb, brootKM, m);
  test_system(H2MATRIX, "Interpolation", Vfull, KMfull, brootV, bem_slp, V2,
	      brootKM, bem_dlp, KM2, row_basis, col_basis, exterior,
	      error_min, error_max);

  /*
   * Test Greenhybrid
   */

  m = 2;
  l = 1;
  delta = 1.0;
  eps_aca = 5.0e-3;

  setup_h2matrix_aprx_greenhybrid_bem3d(bem_slp, Vrb, Vcb, brootV, m, l,
					delta, eps_aca,
					build_bem3d_cube_quadpoints);
  setup_h2matrix_aprx_greenhybrid_bem3d(bem_dlp, KMrb, KMcb, brootKM, m, l,
					delta, eps_aca,
					build_bem3d_cube_quadpoints);
  test_system(H2MATRIX, "Greenhybrid", Vfull, KMfull, brootV, bem_slp, V2,
	      brootKM, bem_dlp, KM2, row_basis, col_basis, exterior,
	      error_min, error_max);

  setup_h2matrix_aprx_greenhybrid_ortho_bem3d(bem_slp, Vrb, Vcb, brootV, m, l,
					      delta, eps_aca,
					      build_bem3d_cube_quadpoints);
  setup_h2matrix_aprx_greenhybrid_ortho_bem3d(bem_dlp, KMrb, KMcb, brootKM, m,
					      l, delta, eps_aca,
					      build_bem3d_cube_quadpoints);
  test_system(H2MATRIX, "Greenhybrid ortho", Vfull, KMfull, brootV, bem_slp,
	      V2, brootKM, bem_dlp, KM2, row_basis, col_basis, exterior,
	      error_min, error_max);

  del_h2matrix(V2);
  del_h2matrix(KM2);
  del_block(brootV);
  del_block(brootKM);
  freemem(rootn->idx);
  freemem(rootd->idx);
  del_cluster(rootn);
  del_cluster(rootd);
  del_amatrix(Vfull);
  del_amatrix(KMfull);
  del_helmholtz_bem3d(bem_slp);
  del_helmholtz_bem3d(bem_dlp);
}

int
main(int argc, char **argv)
{
  pmacrosurface3d mg;
  psurface3d gr;
  uint      n, q, clf;
  real      eta;
  field     k;

  init_h2lib(&argc, &argv);

  n = 512;
  k = 2.0;
  q = 2;
  eta = 1.0;

  mg = new_sphere_macrosurface3d();
  gr = build_from_macrosurface3d_surface3d(mg, REAL_SQRT(n * 0.125));
  n = gr->triangles;

  write_surface3d(gr, "../H2Lib_trunk/sphere.tri");

  printf("Testing unit sphere with %d triangles and %d vertices\n",
	 gr->triangles, gr->vertices);

  clf = 32;

  printf("----------------------------------------\n");
  printf("Testing outer Boundary integral equations:\n");
  printf("----------------------------------------\n\n");

  test_suite(gr, k, q, clf, eta, BASIS_CONSTANT_BEM3D, BASIS_CONSTANT_BEM3D,
	     true, 1.0e-3, 2.0e-3);

  clf = 16;

  printf("----------------------------------------\n");
  printf("Testing outer Boundary integral equations:\n");
  printf("----------------------------------------\n\n");

  test_suite(gr, k, q, clf, eta, BASIS_LINEAR_BEM3D, BASIS_LINEAR_BEM3D,
	     true, 2.5e-4, 3.5e-4);

  del_surface3d(gr);
  del_macrosurface3d(mg);

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
	  "Need to activate complex numbers by \'-DUSE_COMPLEX\' to use this test\n");

  return 0;
}
#endif
