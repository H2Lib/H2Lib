#include <stdio.h>
#include "settings.h"
#include "hmatrix.h"
#include "harith.h"

#include "laplacebem2d.h"

#ifdef USE_CAIRO
#include "cairo/cairo.h"
#endif

static uint problems = 0;

#define IS_IN_RANGE(a, b, c) (((a) <= (b)) && ((b) <= (c)))

static void
check_addhmatrix(phmatrix a, real tol)
{
  phmatrix  acopy;
  real      error;

  acopy = clone_hmatrix(a);

  add_hmatrix(-1.0, a, 0, tol, acopy);
  error = norm2_hmatrix(acopy);
  (void) printf("Checking add_hmatrix\n"
		"  Accuracy %g, %sokay\n", error,
		(IS_IN_RANGE(0.0, error, 1.0e-14) ? "" : "    NOT "));
  if (!IS_IN_RANGE(0.0, error, 1.0e-14))
    problems++;

  del_hmatrix(acopy);
}

static void
check_triangularsolve(bool lower, bool unit, bool atrans,
		      pchmatrix a, bool xtrans, real tol)
{
  uint      n = a->rc->size;
  amatrix   xtmp, btmp;
  pamatrix  x, b;
  avector   xvtmp, bvtmp;
  pavector  xv, bv;
  pblock block;
  phmatrix  xh, bh;
  real      error, norm, eta;

  assert(n == a->cc->size);

  /*
   * amatrix
   */

  x = init_amatrix(&xtmp, n, n);
  random_amatrix(x);

  b = init_zero_amatrix(&btmp, n, n);

  addmul_hmatrix_amatrix_amatrix(1.0, atrans, a, xtrans, x, xtrans, b);

  triangularinvmul_hmatrix_amatrix(lower, unit, atrans, a, xtrans, b);

  add_amatrix(-1.0, false, x, b);
  norm = norm2_amatrix(x);
  error = norm2_amatrix(b) / norm;

  (void) printf("Checking amatrix triangularsolve\n"
		"  (lower=%s, unit=%s, atrans=%s, xtrans=%s)\n"
		"  Accuracy %g, %sokay\n", (lower ? "tr" : "fl"),
		(unit ? "tr" : "fl"), (atrans ? "tr" : "fl"),
		(xtrans ? "tr" : "fl"), error,
		(IS_IN_RANGE(0.0, error, 1.0e-13) ? "" : "    NOT "));
  if (!IS_IN_RANGE(0.0, error, 1.0e-13))
    problems++;

  copy_amatrix(false, x, b);
  triangularmul_hmatrix_amatrix(lower, unit, atrans, a, xtrans, b);

  triangularinvmul_hmatrix_amatrix(lower, unit, atrans, a, xtrans, b);

  add_amatrix(-1.0, false, x, b);
  error = norm2_amatrix(b) / norm;

  (void) printf("Checking amatrix triangulareval/triangularsolve\n"
		"  (lower=%s, unit=%s, atrans=%s, xtrans=%s):\n"
		"  Accuracy %g, %sokay\n", (lower ? "tr" : "fl"),
		(unit ? "tr" : "fl"), (atrans ? "tr" : "fl"),
		(xtrans ? "tr" : "fl"), error,
		(IS_IN_RANGE(0.0, error, 1.0e-13) ? "" : "    NOT "));
  if (!IS_IN_RANGE(0.0, error, 1.0e-13))
    problems++;

  /*
   * avector
   */

  xv = init_avector(&xvtmp, n);
  random_avector(xv);

  bv = init_avector(&bvtmp, n);
  clear_avector(bv);

  if (atrans) {
    addevaltrans_hmatrix_avector(1.0, a, xv, bv);
  }
  else {
    addeval_hmatrix_avector(1.0, a, xv, bv);
  }

  triangularsolve_hmatrix_avector(lower, unit, atrans, a, bv);

  add_avector(-1.0, xv, bv);
  norm = norm2_avector(xv);
  error = norm2_avector(bv) / norm;

  (void) printf("Checking avector triangularsolve\n"
		"  (lower=%s, unit=%s, atrans=%s)\n"
		"  Accuracy %g, %sokay\n", (lower ? "tr" : "fl"),
		(unit ? "tr" : "fl"), (atrans ? "tr" : "fl"), error,
		(IS_IN_RANGE(0.0, error, 1.0e-13) ? "" : "    NOT "));
  if (!IS_IN_RANGE(0.0, error, 1.0e-13))
    problems++;

  copy_avector(xv, bv);
  triangulareval_hmatrix_avector(lower, unit, atrans, a, bv);

  triangularsolve_hmatrix_avector(lower, unit, atrans, a, bv);

  add_avector(-1.0, xv, bv);
  error = norm2_avector(bv) / norm;

  (void) printf("Checking avector triangulareval/triangularsolve\n"
		"  (lower=%s, unit=%s, atrans=%s):\n"
		"  Accuracy %g, %sokay\n", (lower ? "tr" : "fl"),
		(unit ? "tr" : "fl"), (atrans ? "tr" : "fl"), error,
		(IS_IN_RANGE(0.0, error, 1.0e-13) ? "" : "    NOT "));
  if (!IS_IN_RANGE(0.0, error, 1.0e-13))
    problems++;

  /*
   * hmatrix
   */

  eta = 1.0;

  block = build_nonstrict_block((pcluster) a->rc, (pcluster) a->cc, &eta,
				admissible_max_cluster);

  xh = build_from_block_hmatrix(block, 0);
  clear_hmatrix(xh);
  addmul_hmatrix(1.0, false, a, true, a, 0, tol, xh);

  bh = build_from_block_hmatrix(block, 0);
  clear_hmatrix(bh);

  if (xtrans)
    addmul_hmatrix(1.0, false, xh, !atrans, a, 0, tol, bh);
  else
    addmul_hmatrix(1.0, atrans, a, false, xh, 0, tol, bh);

  triangularinvmul_hmatrix(lower, unit, atrans, a, 0, tol, xtrans, bh);

  add_hmatrix(-1.0, xh, 0, tol, bh);
  norm = norm2_hmatrix(xh);
  error = norm2_hmatrix(bh) / norm;

  (void) printf("Checking hmatrix triangularsolve\n"
		"  (lower=%s, unit=%s, atrans=%s, xtrans=%s)\n"
		"  Accuracy %g, %sokay\n", (lower ? "tr" : "fl"),
		(unit ? "tr" : "fl"), (atrans ? "tr" : "fl"),
		(xtrans ? "tr" : "fl"), error,
		(IS_IN_RANGE(0.0, error, 8.0e-13) ? "" : "    NOT "));
  if (!IS_IN_RANGE(0.0, error, 8.0e-13))
    problems++;

  copy_hmatrix(xh, bh);
  triangularmul_hmatrix(lower, unit, atrans, a, 0, tol, xtrans, bh);

  triangularinvmul_hmatrix(lower, unit, atrans, a, 0, tol, xtrans, bh);

  add_hmatrix(-1.0, xh, 0, tol, bh);
  error = norm2_hmatrix(bh) / norm;

  (void) printf("Checking hmatrix triangulareval/triangularsolve\n"
		"  (lower=%s, unit=%s, atrans=%s, xtrans=%s):\n"
		"  Accuracy %g, %sokay\n", (lower ? "tr" : "fl"),
		(unit ? "tr" : "fl"), (atrans ? "tr" : "fl"),
		(xtrans ? "tr" : "fl"), error,
		(IS_IN_RANGE(0.0, error, 8.0e-13) ? "" : "    NOT "));
  if (!IS_IN_RANGE(0.0, error, 8.0e-13))
    problems++;

  del_block(block);
  del_hmatrix(xh);
  del_hmatrix(bh);
  uninit_amatrix(b);
  uninit_amatrix(x);
  uninit_avector(bv);
  uninit_avector(xv);
}

static void
set_unit(phmatrix R)
{
  const uint sons = R->rsons;
  const uint minsons = UINT_MIN(sons, R->csons);

  uint      i, m;

  if (R->son != NULL) {
    for (i = 0; i < minsons; ++i) {
      set_unit(R->son[i + i * sons]);
    }
  }
  else {
    assert(R->f != NULL);

    m = UINT_MIN(R->f->rows, R->f->cols);
    for (i = 0; i < m; ++i) {
      R->f->a[i + i * R->f->ld] = 1.0;
    }
  }
}

int
main()
{
  phmatrix  a, acopy, L, R;
  pavector  x, b;
  uint      n;
  real      error;
  pcurve2d  gr2;
  pbem2d    bem2;
  pcluster  root2;
  pblock    block2;
#ifdef USE_CAIRO
  cairo_t  *cr;
#endif
  uint      clf, m;
  real      tol, eta, delta, eps_aca;

  n = 579;
  tol = 1.0e-13;

  clf = 16;
  eta = 1.0;
  m = 4;
  delta = 1.0;
  eps_aca = 1.0e-13;

  gr2 = new_circle_curve2d(n, 0.333);
  bem2 = new_slp_laplace_bem2d(gr2, 2, BASIS_CONSTANT_BEM2D);
  root2 = build_bem2d_cluster(bem2, clf, BASIS_CONSTANT_BEM2D);
  block2 = build_nonstrict_block(root2, root2, &eta, admissible_max_cluster);
  setup_hmatrix_aprx_greenhybrid_row_bem2d(bem2, root2, root2, block2, m, 1,
					   delta, eps_aca,
					   build_bem2d_rect_quadpoints);
  setup_hmatrix_recomp_bem2d(bem2, true, eps_aca, false, eps_aca);

  (void) printf("----------------------------------------\n"
		"Check %u x %u H-matrix addition\n", n, n);

  (void) printf("Creating laplacebem2d SLP matrix\n");

  a = build_from_block_hmatrix(block2, 0);
  assemble_bem2d_hmatrix(bem2, block2, a);

  check_addhmatrix(a, tol);

  del_hmatrix(a);

  (void) printf("----------------------------------------\n"
		"Check %u x %u Cholesky factorization\n", n, n);

  (void) printf("Creating laplacebem2d SLP matrix\n");

  a = build_from_block_hmatrix(block2, 0);
  assemble_bem2d_hmatrix(bem2, block2, a);

  (void) printf("Creating random solution and right-hand side\n");
  x = new_avector(n);
  random_avector(x);
  b = new_avector(n);
  clear_avector(b);
  addevalsymm_hmatrix_avector(1.0, a, x, b);

  (void) printf("Copying matrix\n");
  acopy = clone_hmatrix(a);

  (void) printf("Computing Cholesky factorization\n");
  choldecomp_hmatrix(a, 0, tol);

  (void) printf("Solving\n");
  cholsolve_hmatrix_avector(a, b);

  add_avector(-1.0, x, b);
  error = norm2_avector(b) / norm2_avector(x);
  (void) printf("  Accuracy %g, %sokay\n", error,
		IS_IN_RANGE(0.0, error, 4.0e-13) ? "" : "    NOT ");
  if (!IS_IN_RANGE(0.0, error, 4.0e-13))
    problems++;

  (void) printf("Evaluating\n");
  copy_avector(x, b);
  triangulareval_hmatrix_avector(true, false, true, a, b);
  triangulareval_hmatrix_avector(true, false, false, a, b);
  addevalsymm_hmatrix_avector(-1.0, acopy, x, b);
  error = norm2_avector(b) / norm2_avector(x);
  (void) printf("  Accuracy %g, %sokay\n", error,
		IS_IN_RANGE(0.0, error, 4.0e-13) ? "" : "    NOT ");
  if (!IS_IN_RANGE(0.0, error, 4.0e-13))
    problems++;

  (void) printf("Building triangular factor\n");
  L = clone_lower_hmatrix(false, a);

  copy_avector(x, b);
  triangulareval_hmatrix_avector(true, false, false, a, b);
  mvm_hmatrix_avector(-1.0, false, L, x, b);
  error = norm2_avector(b) / norm2_avector(x);
  (void) printf("  Accuracy %g, %sokay\n", error,
		IS_IN_RANGE(0.0, error, 4.0e-13) ? "" : "    NOT ");
  if (!IS_IN_RANGE(0.0, error, 4.0e-13))
    problems++;

  copy_avector(x, b);
  triangulareval_hmatrix_avector(true, false, true, a, b);
  mvm_hmatrix_avector(-1.0, true, L, x, b);
  error = norm2_avector(b) / norm2_avector(x);
  (void) printf("  Accuracy %g, %sokay\n", error,
		IS_IN_RANGE(0.0, error, 4.0e-13) ? "" : "    NOT ");
  if (!IS_IN_RANGE(0.0, error, 4.0e-13))
    problems++;

  (void) printf("Checking factorization\n");
  error = norm2_hmatrix(acopy);
  addmul_hmatrix(-1.0, false, L, true, L, 0, tol, acopy);
  error = norm2_hmatrix(acopy) / error;
  (void) printf("  Accuracy %g, %sokay\n", error,
		IS_IN_RANGE(0.0, error, 1.0e-14) ? "" : "    NOT ");
  if (!IS_IN_RANGE(0.0, error, 1.0e-14))
    problems++;

  del_hmatrix(acopy);
  del_hmatrix(a);
  del_hmatrix(L);
  del_avector(b);
  del_avector(x);

  (void) printf("----------------------------------------\n"
		"Check %u x %u LR factorization\n", n, n);

  (void) printf("Creating laplacebem2d SLP matrix\n");
  a = build_from_block_hmatrix(block2, 0);
  assemble_bem2d_hmatrix(bem2, block2, a);

  (void) printf("Creating random solution and right-hand side\n");
  x = new_avector(n);
  random_avector(x);
  b = new_avector(n);
  clear_avector(b);
  mvm_hmatrix_avector(1.0, false, a, x, b);

  (void) printf("Copying matrix\n");
  acopy = clone_hmatrix(a);

  (void) printf("Computing LR factorization\n");
  lrdecomp_hmatrix(a, 0, tol);

  (void) printf("Solving\n");
  lrsolve_hmatrix_avector(false, a, b);

  add_avector(-1.0, x, b);
  error = norm2_avector(b) / norm2_avector(x);
  (void) printf("  Accuracy %g, %sokay\n", error,
		IS_IN_RANGE(0.0, error, 6.0e-13) ? "" : "    NOT ");
  if (!IS_IN_RANGE(0.0, error, 6.0e-13))
    problems++;

  (void) printf("Evaluating\n");
  copy_avector(x, b);
  triangulareval_hmatrix_avector(false, false, false, a, b);
  triangulareval_hmatrix_avector(true, true, false, a, b);
  mvm_hmatrix_avector(-1.0, false, acopy, x, b);
  error = norm2_avector(b) / norm2_avector(x);
  (void) printf("  Accuracy %g, %sokay\n", error,
		IS_IN_RANGE(0.0, error, 6.0e-13) ? "" : "    NOT ");
  if (!IS_IN_RANGE(0.0, error, 6.0e-13))
    problems++;

  copy_avector(x, b);
  triangulareval_hmatrix_avector(true, true, true, a, b);
  triangulareval_hmatrix_avector(false, false, true, a, b);
  mvm_hmatrix_avector(-1.0, true, acopy, x, b);
  error = norm2_avector(b) / norm2_avector(x);
  (void) printf("  Accuracy %g, %sokay\n", error,
		IS_IN_RANGE(0.0, error, 6.0e-13) ? "" : "    NOT ");
  if (!IS_IN_RANGE(0.0, error, 6.0e-13))
    problems++;

  (void) printf("Building triangular factors\n");
  L = clone_lower_hmatrix(true, a);
  R = clone_upper_hmatrix(false, a);

  copy_avector(x, b);
  triangulareval_hmatrix_avector(false, false, false, a, b);
  mvm_hmatrix_avector(-1.0, false, R, x, b);
  error = norm2_avector(b) / norm2_avector(x);
  (void) printf("  Accuracy %g, %sokay\n", error,
		IS_IN_RANGE(0.0, error, 6.0e-13) ? "" : "    NOT ");
  if (!IS_IN_RANGE(0.0, error, 6.0e-13))
    problems++;

  copy_avector(x, b);
  triangulareval_hmatrix_avector(false, false, true, a, b);
  mvm_hmatrix_avector(-1.0, true, R, x, b);
  error = norm2_avector(b) / norm2_avector(x);
  (void) printf("  Accuracy %g, %sokay\n", error,
		IS_IN_RANGE(0.0, error, 6.0e-13) ? "" : "    NOT ");
  if (!IS_IN_RANGE(0.0, error, 6.0e-13))
    problems++;

  copy_avector(x, b);
  triangulareval_hmatrix_avector(true, true, false, a, b);
  mvm_hmatrix_avector(-1.0, false, L, x, b);
  error = norm2_avector(b) / norm2_avector(x);
  (void) printf("  Accuracy %g, %sokay\n", error,
		IS_IN_RANGE(0.0, error, 6.0e-13) ? "" : "    NOT ");
  if (!IS_IN_RANGE(0.0, error, 6.0e-13))
    problems++;

  copy_avector(x, b);
  triangulareval_hmatrix_avector(true, true, true, a, b);
  mvm_hmatrix_avector(-1.0, true, L, x, b);
  error = norm2_avector(b) / norm2_avector(x);
  (void) printf("  Accuracy %g, %sokay\n", error,
		IS_IN_RANGE(0.0, error, 6.0e-13) ? "" : "    NOT ");
  if (!IS_IN_RANGE(0.0, error, 6.0e-13))
    problems++;

  (void) printf("Checking factorization\n");
  error = norm2_hmatrix(acopy);
  addmul_hmatrix(-1.0, false, L, false, R, 0, tol, acopy);
  error = norm2_hmatrix(acopy) / error;
  (void) printf("  Accuracy %g, %sokay\n", error,
		IS_IN_RANGE(0.0, error, 4.0e-16) ? "" : "    NOT ");
  if (!IS_IN_RANGE(0.0, error, 4.0e-16))
    problems++;

  /* Check forward/backward substitution */
  (void) printf("----------------------------------------\n"
		"Checking lower triangular solve/eval\n");
  check_triangularsolve(true, true, false, L, false, tol);
  check_triangularsolve(true, true, false, L, true, tol);
  check_triangularsolve(true, true, true, L, false, tol);
  check_triangularsolve(true, true, true, L, true, tol);
  check_triangularsolve(true, false, false, L, false, tol);
  check_triangularsolve(true, false, false, L, true, tol);
  check_triangularsolve(true, false, true, L, false, tol);
  check_triangularsolve(true, false, true, L, true, tol);

  (void) printf("----------------------------------------\n"
		"Checking upper triangular solve/eval\n");
  check_triangularsolve(false, false, false, R, false, tol);
  check_triangularsolve(false, false, false, R, true, tol);
  check_triangularsolve(false, false, true, R, false, tol);
  check_triangularsolve(false, false, true, R, true, tol);
  set_unit(R);
  check_triangularsolve(false, true, false, R, false, tol);
  check_triangularsolve(false, true, false, R, true, tol);
  check_triangularsolve(false, true, true, R, false, tol);
  check_triangularsolve(false, true, true, R, true, tol);

#ifdef USE_CAIRO
  /* Check Cairo drawing */
  (void) printf("----------------------------------------\n"
		"Checking Cairo drawing\n" "  Drawing to \"hmatrix.pdf\"\n");
  cr = new_cairopdf("hmatrix.pdf", 600.0, 600.0);
  draw_cairo_hmatrix(cr, a, true, 0);
  cairo_destroy(cr);
#endif

  /* Final clean-up */
  (void) printf("Cleaning up\n");
  del_hmatrix(acopy);
  del_hmatrix(a);
  del_hmatrix(L);
  del_hmatrix(R);
  del_avector(b);
  del_avector(x);

  del_bem2d(bem2);
  del_block(block2);
  freemem(root2->idx);
  del_cluster(root2);
  del_curve2d(gr2);

  (void) printf("----------------------------------------\n"
		"  %u matrices and\n"
		"  %u vectors still active\n"
		"  %u errors found\n", getactives_amatrix(),
		getactives_avector(), problems);

  return problems;
}
