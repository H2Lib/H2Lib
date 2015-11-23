#include "basic.h"
#include "krylov.h"
#include "laplacebem3d.h"

#ifdef USE_CAIRO
#include <cairo.h>
#endif

static uint problems = 0;

#define IS_IN_RANGE(a, b, c) (((a) <= (b)) && ((b) <= (c)))

/* Compute the L_2-error for constant basis functions */
real
L2gamma_c_diff_norm2(pbem3d bem, pavector x, boundary_func3d rhs, void *data)
{
  uint      n = x->dim;
  pcsurface3d gr = bem->gr;
  const     real(*gr_x)[3] = (const real(*)[3]) gr->x;
  const     uint(*gr_t)[3] = (const uint(*)[3]) gr->t;
  const     real(*gr_n)[3] = (const real(*)[3]) gr->n;
  const real *gr_g = gr->g;
  real     *xq = bem->sq->x_single;
  real     *yq = bem->sq->y_single;
  uint      nq = bem->sq->n_single;
  real     *wq = bem->sq->w_single + 3 * nq;
  field    *xv = x->v;

  const real *A, *B, *C, *N;
  real      X[3];
  real      sum, tx, sx, Ax, Bx, Cx;
  uint      t, q;

  real      norm;

  norm = 0.0;
  for (t = 0; t < n; ++t) {
    A = gr_x[gr_t[t][0]];
    B = gr_x[gr_t[t][1]];
    C = gr_x[gr_t[t][2]];
    N = gr_n[t];

    sum = 0.0;
    for (q = 0; q < nq; ++q) {
      tx = xq[q];
      sx = yq[q];
      Ax = 1.0 - tx;
      Bx = tx - sx;
      Cx = sx;

      X[0] = A[0] * Ax + B[0] * Bx + C[0] * Cx;
      X[1] = A[1] * Ax + B[1] * Bx + C[1] * Cx;
      X[2] = A[2] * Ax + B[2] * Bx + C[2] * Cx;

      sum += wq[q] * REAL_SQR(rhs(X, N, data) - xv[t]);
    }
    norm += gr_g[t] * sum;
  }
  norm = REAL_SQRT(norm);

  return norm;
}

/* Compute the L_2-error for linear basis functions */
real
L2gamma_l_diff_norm2(pbem3d bem, pavector x, boundary_func3d rhs, void *data)
{
  uint      n = x->dim;
  pcsurface3d gr = bem->gr;
  const     real(*gr_x)[3] = (const real(*)[3]) gr->x;
  const     uint(*gr_t)[3] = (const uint(*)[3]) gr->t;
  const     real(*gr_n)[3] = (const real(*)[3]) gr->n;
  const real *gr_g = gr->g;
  real     *xq = bem->sq->x_single;
  real     *yq = bem->sq->y_single;
  uint      nq = bem->sq->n_single;
  real     *wq = bem->sq->w_single + 3 * nq;
  field    *xv = x->v;

  const real *A, *B, *C, *N;
  real      X[3];
  real      sum, bf, tx, sx, Ax, Bx, Cx;
  uint      t, q;

  real      norm;

  norm = 0.0;
  for (t = 0; t < n; ++t) {
    A = gr_x[gr_t[t][0]];
    B = gr_x[gr_t[t][1]];
    C = gr_x[gr_t[t][2]];
    N = gr_n[t];

    sum = 0.0;
    for (q = 0; q < nq; ++q) {
      tx = xq[q];
      sx = yq[q];
      Ax = 1.0 - tx;
      Bx = tx - sx;
      Cx = sx;

      X[0] = A[0] * Ax + B[0] * Bx + C[0] * Cx;
      X[1] = A[1] * Ax + B[1] * Bx + C[1] * Cx;
      X[2] = A[2] * Ax + B[2] * Bx + C[2] * Cx;

      bf = xv[gr_t[t][0]] * Ax + xv[gr_t[t][1]] * Bx + xv[gr_t[t][2]] * Cx;

      sum += wq[q] * REAL_SQR(rhs(X, N, data) - bf);
    }
    norm += gr_g[t] * sum;
  }
  norm = REAL_SQRT(norm);

  return norm;
}

/* Simple convenience wrapper for conjugate gradient solver */
static void
solve_cg_bem3d(matrixtype type, void *A, pavector b, pavector x,
	       real accuracy, uint steps)
{
  addeval_t addevalA;
  pavector  r, p, a;
  uint      i, n;
  real      norm;

  switch (type) {
  case AMATRIX:
    addevalA = (addeval_t) addeval_amatrix_avector;
    break;
  case HMATRIX:
    addevalA = (addeval_t) addeval_hmatrix_avector;
    break;
  case H2MATRIX:
    addevalA = (addeval_t) addeval_h2matrix_avector;
    break;
  default:
    printf("ERROR: unknown matrix type!\n");
    abort();
    break;
  }

  n = b->dim;
  assert(x->dim == n);

  r = new_avector(n);
  p = new_avector(n);
  a = new_avector(n);
  random_avector(x);

  init_cg(addevalA, A, b, x, r, p, a);

  for (i = 1; i < steps; i++) {
    step_cg(addevalA, A, b, x, r, p, a);
    norm = norm2_avector(r);
#ifndef NDEBUG
    printf("  Residual: %.5e\t Iterations: %u\r", norm, i);
    fflush(stdout);
#endif
    if (norm <= accuracy) {
      break;
    }
  }
#ifndef NDEBUG
  printf("  Residual: %.5e\t Iterations: %u\n", norm2_avector(r), i);
#endif

  del_avector(r);
  del_avector(p);
  del_avector(a);

}

static void
test_system(matrixtype mattype, const char *apprxtype,
	    pcamatrix Vfull, pcamatrix KMfull, pblock brootV, pbem3d bem_slp,
	    void *V, pblock brootKM, pbem3d bem_dlp, void *KM,
	    basisfunctionbem3d basis_neumann,
	    basisfunctionbem3d basis_dirichlet, bool exterior, real low,
	    real high)
{
  pavector  x, f, b;
  real      errorV, errorKM, error_solve, eps_solve;
  uint      steps;

  eps_solve = 1.0e-12;
  steps = 1000;

  printf("Testing: %c%c %s %s %s\n"
	 "====================================\n\n",
	 basis_neumann == BASIS_LINEAR_BEM3D ? 'l' : 'c',
	 basis_dirichlet == BASIS_LINEAR_BEM3D ? 'l' : 'c',
	 (exterior == true ? "exterior" : "interior"),
	 (mattype == HMATRIX ? "Hmatrix" : "H2matrix"), apprxtype);

  if (mattype == HMATRIX) {
    assemble_bem3d_hmatrix(bem_slp, brootV, (phmatrix) V);
    assemble_bem3d_hmatrix(bem_dlp, brootKM, KM);
    errorV = norm2diff_amatrix_hmatrix((phmatrix) V, Vfull)
      / norm2_amatrix(Vfull);
    printf("rel. error V       : %.5e\n", errorV);
    errorKM = norm2diff_amatrix_hmatrix(KM, KMfull) / norm2_amatrix(KMfull);
    printf("rel. error K%c0.5*M : %.5e\n", (exterior == true ? '-' : '+'),
	   errorKM);
  }
  else {
    assert(mattype == H2MATRIX);
    assemble_bem3d_h2matrix_row_clusterbasis(bem_slp, ((ph2matrix) V)->rb);
    assemble_bem3d_h2matrix_col_clusterbasis(bem_slp, ((ph2matrix) V)->cb);
    assemble_bem3d_h2matrix(bem_slp, brootV, (ph2matrix) V);

    assemble_bem3d_h2matrix_row_clusterbasis(bem_dlp, ((ph2matrix) KM)->rb);
    assemble_bem3d_h2matrix_col_clusterbasis(bem_dlp, ((ph2matrix) KM)->cb);
    assemble_bem3d_h2matrix(bem_dlp, brootKM, (ph2matrix) KM);

    errorV = norm2diff_amatrix_h2matrix((ph2matrix) V, Vfull)
      / norm2_amatrix(Vfull);
    printf("rel. error V       : %.5e\n", errorV);
    errorKM = norm2diff_amatrix_h2matrix((ph2matrix) KM, KMfull)
      / norm2_amatrix(KMfull);
    printf("rel. error K%c0.5*M : %.5e\n", (exterior == true ? '-' : '+'),
	   errorKM);
  }

  x = new_avector(Vfull->cols);
  f = new_avector(KMfull->cols);
  b = new_avector(KMfull->rows);

  printf("Solving Dirichlet problem:\n");

  if (basis_dirichlet == BASIS_LINEAR_BEM3D) {
    projectl2_bem3d_linear_avector(bem_dlp,
				   (exterior ==
				    true ?
				    eval_dirichlet_fundamental2_laplacebem3d :
				    eval_dirichlet_fundamental_laplacebem3d),
				   f, NULL);
  }
  else {
    projectl2_bem3d_const_avector(bem_dlp,
				  (exterior ==
				   true ?
				   eval_dirichlet_fundamental2_laplacebem3d :
				   eval_dirichlet_fundamental_laplacebem3d),
				  f, NULL);
  }
  clear_avector(b);
  if (mattype == HMATRIX) {
    addeval_hmatrix_avector(1.0, (phmatrix) KM, f, b);
  }
  else {
    assert(mattype == H2MATRIX);
    addeval_h2matrix_avector(1.0, (ph2matrix) KM, f, b);
  }

  solve_cg_bem3d(mattype, V, b, x, eps_solve, steps);

  if (basis_neumann == BASIS_LINEAR_BEM3D) {
    error_solve = L2gamma_l_diff_norm2(bem_slp, x,
				       (exterior ==
					true ?
					eval_neumann_fundamental2_laplacebem3d
					:
					eval_neumann_fundamental_laplacebem3d),
				       NULL);
  }
  else {
    error_solve = L2gamma_c_diff_norm2(bem_slp, x,
				       (exterior ==
					true ?
					eval_neumann_fundamental2_laplacebem3d
					:
					eval_neumann_fundamental_laplacebem3d),
				       NULL);
  }

  clear_avector(x);
  if (basis_neumann == BASIS_LINEAR_BEM3D) {
    error_solve = error_solve
      / L2gamma_l_diff_norm2(bem_slp, x,
			     (exterior ==
			      true ? eval_neumann_fundamental2_laplacebem3d :
			      eval_neumann_fundamental_laplacebem3d), NULL);
  }
  else {
    error_solve = error_solve
      / L2gamma_c_diff_norm2(bem_slp, x,
			     (exterior ==
			      true ? eval_neumann_fundamental2_laplacebem3d :
			      eval_neumann_fundamental_laplacebem3d), NULL);
  }

  printf("rel. error neumann : %.5e       %s\n", error_solve,
	 (IS_IN_RANGE(low, error_solve, high) ? "    okay" : "NOT okay"));

  if (!IS_IN_RANGE(low, error_solve, high))
    problems++;

  printf("\n");

  del_avector(x);
  del_avector(f);
  del_avector(b);
}

void
test_suite(pcsurface3d gr, uint q, uint clf, real eta,
	   basisfunctionbem3d basis_neumann,
	   basisfunctionbem3d basis_dirichlet, bool exterior, real error_min,
	   real error_max)
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

  nn = basis_neumann == BASIS_LINEAR_BEM3D ? gr->vertices : gr->triangles;
  nd = basis_dirichlet == BASIS_LINEAR_BEM3D ? gr->vertices : gr->triangles;

  bem_slp = new_slp_laplace_bem3d(gr, q, q + 2, basis_neumann);
  bem_dlp =
    new_dlp_laplace_bem3d(gr, q, q + 2, basis_neumann, basis_dirichlet,
			  exterior ? -0.5 : 0.5);

  rootn = build_bem3d_cluster(bem_slp, clf, basis_neumann);
  rootd = build_bem3d_cluster(bem_dlp, clf, basis_dirichlet);

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

  m = 3;

  setup_hmatrix_aprx_inter_row_bem3d(bem_slp, rootn, rootn, brootV, m);
  setup_hmatrix_aprx_inter_row_bem3d(bem_dlp, rootn, rootd, brootKM, m);
  test_system(HMATRIX, "Interpolation row", Vfull, KMfull, brootV, bem_slp, V,
	      brootKM, bem_dlp, KM, basis_neumann, basis_dirichlet, exterior,
	      error_min, error_max);

  setup_hmatrix_aprx_inter_col_bem3d(bem_slp, rootn, rootn, brootV, m);
  setup_hmatrix_aprx_inter_col_bem3d(bem_dlp, rootn, rootd, brootKM, m);
  test_system(HMATRIX, "Interpolation column", Vfull, KMfull, brootV, bem_slp,
	      V, brootKM, bem_dlp, KM, basis_neumann, basis_dirichlet,
	      exterior, error_min, error_max);

  setup_hmatrix_aprx_inter_mixed_bem3d(bem_slp, rootn, rootn, brootV, m);
  setup_hmatrix_aprx_inter_mixed_bem3d(bem_dlp, rootn, rootd, brootKM, m);
  test_system(HMATRIX, "Interpolation mixed", Vfull, KMfull, brootV, bem_slp,
	      V, brootKM, bem_dlp, KM, basis_neumann, basis_dirichlet,
	      exterior, error_min, error_max);

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
	      brootKM, bem_dlp, KM, basis_neumann, basis_dirichlet, exterior,
	      error_min, error_max);
  setup_hmatrix_aprx_green_col_bem3d(bem_slp, rootn, rootn, brootV, m, l,
				     delta, build_bem3d_cube_quadpoints);
  setup_hmatrix_aprx_green_col_bem3d(bem_dlp, rootn, rootd, brootKM, m, l,
				     delta, build_bem3d_cube_quadpoints);
  test_system(HMATRIX, "Green column", Vfull, KMfull, brootV, bem_slp, V,
	      brootKM, bem_dlp, KM, basis_neumann, basis_dirichlet, exterior,
	      error_min, error_max);

  setup_hmatrix_aprx_green_mixed_bem3d(bem_slp, rootn, rootn, brootV, m, l,
				       delta, build_bem3d_cube_quadpoints);
  setup_hmatrix_aprx_green_mixed_bem3d(bem_dlp, rootn, rootd, brootKM, m, l,
				       delta, build_bem3d_cube_quadpoints);
  test_system(HMATRIX, "Green mixed", Vfull, KMfull, brootV, bem_slp, V,
	      brootKM, bem_dlp, KM, basis_neumann, basis_dirichlet, exterior,
	      error_min, error_max);

  /*
   * Test Greenhybrid
   */

  m = 2;
  l = 1;
  delta = 1.0;
  eps_aca = 1.0e-2;

  setup_hmatrix_aprx_greenhybrid_row_bem3d(bem_slp, rootn, rootn, brootV, m,
					   l, delta, eps_aca,
					   build_bem3d_cube_quadpoints);
  setup_hmatrix_aprx_greenhybrid_row_bem3d(bem_dlp, rootn, rootd, brootKM, m,
					   l, delta, eps_aca,
					   build_bem3d_cube_quadpoints);
  test_system(HMATRIX, "Greenhybrid row", Vfull, KMfull, brootV, bem_slp, V,
	      brootKM, bem_dlp, KM, basis_neumann, basis_dirichlet, exterior,
	      error_min, error_max);

  setup_hmatrix_aprx_greenhybrid_col_bem3d(bem_slp, rootn, rootn, brootV, m,
					   l, delta, eps_aca,
					   build_bem3d_cube_quadpoints);
  setup_hmatrix_aprx_greenhybrid_col_bem3d(bem_dlp, rootn, rootd, brootKM, m,
					   l, delta, eps_aca,
					   build_bem3d_cube_quadpoints);
  test_system(HMATRIX, "Greenhybrid column", Vfull, KMfull, brootV, bem_slp,
	      V, brootKM, bem_dlp, KM, basis_neumann, basis_dirichlet,
	      exterior, error_min, error_max);

  setup_hmatrix_aprx_greenhybrid_mixed_bem3d(bem_slp, rootn, rootn, brootV, m,
					     l, delta, eps_aca,
					     build_bem3d_cube_quadpoints);
  setup_hmatrix_aprx_greenhybrid_mixed_bem3d(bem_dlp, rootn, rootd, brootKM,
					     m, l, delta, eps_aca,
					     build_bem3d_cube_quadpoints);
  test_system(HMATRIX, "Greenhybrid mixed", Vfull, KMfull, brootV, bem_slp, V,
	      brootKM, bem_dlp, KM, basis_neumann, basis_dirichlet, exterior,
	      error_min, error_max);

  /*
   * Test ACA / PACA / HCA
   */

  m = 2;
  eps_aca = 1.0e-2;

  setup_hmatrix_aprx_aca_bem3d(bem_slp, rootn, rootn, brootV, eps_aca);
  setup_hmatrix_aprx_aca_bem3d(bem_dlp, rootn, rootd, brootKM, eps_aca);
  test_system(HMATRIX, "ACA full pivoting", Vfull, KMfull, brootV, bem_slp, V,
	      brootKM, bem_dlp, KM, basis_neumann, basis_dirichlet, exterior,
	      error_min, error_max);
  setup_hmatrix_aprx_paca_bem3d(bem_slp, rootn, rootn, brootV, eps_aca);
  setup_hmatrix_aprx_paca_bem3d(bem_slp, rootn, rootd, brootKM, eps_aca);
  test_system(HMATRIX, "ACA partial pivoting", Vfull, KMfull, brootV, bem_slp,
	      V, brootKM, bem_dlp, KM, basis_neumann, basis_dirichlet,
	      exterior, error_min, error_max);

  setup_hmatrix_aprx_hca_bem3d(bem_slp, rootn, rootn, brootV, m, eps_aca);
  setup_hmatrix_aprx_hca_bem3d(bem_slp, rootn, rootd, brootKM, m, eps_aca);
  test_system(HMATRIX, "HCA2", Vfull, KMfull, brootV, bem_slp, V, brootKM,
	      bem_dlp, KM, basis_neumann, basis_dirichlet, exterior,
	      error_min, error_max);
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

  m = 3;

  setup_h2matrix_aprx_inter_bem3d(bem_slp, Vrb, Vcb, brootV, m);
  setup_h2matrix_aprx_inter_bem3d(bem_dlp, KMrb, KMcb, brootKM, m);
  test_system(H2MATRIX, "Interpolation", Vfull, KMfull, brootV, bem_slp, V2,
	      brootKM, bem_dlp, KM2, basis_neumann, basis_dirichlet, exterior,
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
	      brootKM, bem_dlp, KM2, basis_neumann, basis_dirichlet, exterior,
	      error_min, error_max);

  setup_h2matrix_aprx_greenhybrid_ortho_bem3d(bem_slp, Vrb, Vcb, brootV, m, l,
					      delta, eps_aca,
					      build_bem3d_cube_quadpoints);
  setup_h2matrix_aprx_greenhybrid_ortho_bem3d(bem_dlp, KMrb, KMcb, brootKM, m,
					      l, delta, eps_aca,
					      build_bem3d_cube_quadpoints);
  test_system(H2MATRIX, "Greenhybrid ortho", Vfull, KMfull, brootV, bem_slp,
	      V2, brootKM, bem_dlp, KM2, basis_neumann, basis_dirichlet,
	      exterior, error_min, error_max);

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
  del_laplace_bem3d(bem_slp);
  del_laplace_bem3d(bem_dlp);
}

int
main(int argc, char **argv)
{
  pmacrosurface3d mg;
  psurface3d gr;
  uint      n, q, clf;
  real      eta;

  init_h2lib(&argc, &argv);

  n = 512;
  q = 2;
  eta = 1.0;

  mg = new_sphere_macrosurface3d();
  gr = build_from_macrosurface3d_surface3d(mg, REAL_SQRT(n * 0.125));

  printf("Testing unit sphere with %d triangles and %d vertices\n",
	 gr->triangles, gr->vertices);

  /****************************************************
   * Neumann: constant, Dirichlet: constant
   ****************************************************/

  clf = 32;

  printf("----------------------------------------\n");
  printf("Testing inner Boundary integral equations:\n");
  printf("----------------------------------------\n\n");

  test_suite(gr, q, clf, eta, BASIS_CONSTANT_BEM3D, BASIS_CONSTANT_BEM3D,
	     false, 7.0e-2, 7.5e-2);

  printf("----------------------------------------\n");
  printf("Testing outer Boundary integral equations:\n");
  printf("----------------------------------------\n\n");

  test_suite(gr, q, clf, eta, BASIS_CONSTANT_BEM3D, BASIS_CONSTANT_BEM3D,
	     true, 6.0e-2, 7.0e-2);

  /****************************************************
   * Neumann: constant, Dirichlet: linear
   ****************************************************/

  clf = 16;

  printf("----------------------------------------\n");
  printf("Testing inner Boundary integral equations:\n");
  printf("----------------------------------------\n\n");

  test_suite(gr, q, clf, eta, BASIS_CONSTANT_BEM3D, BASIS_LINEAR_BEM3D, false,
	     7.0e-2, 8.0e-2);

  printf("----------------------------------------\n");
  printf("Testing outer Boundary integral equations:\n");
  printf("----------------------------------------\n\n");

  test_suite(gr, q, clf, eta, BASIS_CONSTANT_BEM3D, BASIS_LINEAR_BEM3D, true,
	     6.0e-2, 7.0e-2);

  /****************************************************
   * Neumann: linear, Dirichlet: linear
   ****************************************************/

  clf = 16;

  printf("----------------------------------------\n");
  printf("Testing inner Boundary integral equations:\n");
  printf("----------------------------------------\n\n");

  test_suite(gr, q, clf, eta, BASIS_LINEAR_BEM3D, BASIS_LINEAR_BEM3D, false,
	     6.0e-2, 7.0e-2);

  printf("----------------------------------------\n");
  printf("Testing outer Boundary integral equations:\n");
  printf("----------------------------------------\n\n");

  test_suite(gr, q, clf, eta, BASIS_LINEAR_BEM3D, BASIS_LINEAR_BEM3D, true,
	     1.5e-2, 2.5e-2);

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
