#include <stdio.h>

#include "amatrix.h"
#include "factorizations.h"
#include "settings.h"

static uint problems = 0;
static const real tolerance = 1e-12;

#define IS_IN_RANGE(a, b, c) (((a) <= (b)) && ((b) <= (c)))

static void
check_triangularsolve(bool lower, bool unit, bool atrans,
		      pcamatrix a, bool xtrans)
{
  uint      n = a->rows;
  amatrix   xtmp, btmp;
  pamatrix  x, b;
  avector   xvtmp, bvtmp;
  pavector  xv, bv;
  real      error;

  assert(n == a->cols);

  /*
   * amatrix
   */

  x = init_amatrix(&xtmp, n, n);
  random_amatrix(x);

  b = init_zero_amatrix(&btmp, n, n);

  if (xtrans)
    addmul_amatrix(1.0, false, x, !atrans, a, b);
  else
    addmul_amatrix(1.0, atrans, a, false, x, b);

  triangularsolve_amatrix(lower, unit, atrans, a, xtrans, b);

  add_amatrix(-1.0, false, x, b);
  error = norm2_amatrix(b) / norm2_amatrix(x);

  (void) printf("Checking amatrix triangularsolve\n"
		"  (lower=%s, unit=%s, atrans=%s, xtrans=%s)\n"
		"  Accuracy %g, %sokay\n", (lower ? "tr" : "fl"),
		(unit ? "tr" : "fl"), (atrans ? "tr" : "fl"),
		(xtrans ? "tr" : "fl"), error,
		(IS_IN_RANGE(0.0, error, 1.0e-14) ? "" : "    NOT "));
  if (!IS_IN_RANGE(0.0, error, 1.0e-14))
    problems++;

  copy_amatrix(false, x, b);
  triangulareval_amatrix(lower, unit, atrans, a, xtrans, b);

  triangularsolve_amatrix(lower, unit, atrans, a, xtrans, b);

  add_amatrix(-1.0, false, x, b);
  error = norm2_amatrix(b) / norm2_amatrix(x);

  (void) printf("Checking amatrix triangulareval/triangularsolve\n"
		"  (lower=%s, unit=%s, atrans=%s, xtrans=%s):\n"
		"  Accuracy %g, %sokay\n", (lower ? "tr" : "fl"),
		(unit ? "tr" : "fl"), (atrans ? "tr" : "fl"),
		(xtrans ? "tr" : "fl"), error,
		(IS_IN_RANGE(0.0, error, 1.0e-14) ? "" : "    NOT "));
  if (!IS_IN_RANGE(0.0, error, 1.0e-14))
    problems++;

  /*
   * avector
   */

  xv = init_avector(&xvtmp, n);
  random_avector(xv);

  bv = init_avector(&bvtmp, n);
  clear_avector(bv);

  if (atrans) {
    addevaltrans_amatrix_avector(1.0, a, xv, bv);
  }
  else {
    addeval_amatrix_avector(1.0, a, xv, bv);
  }

  triangularsolve_amatrix_avector(lower, unit, atrans, a, bv);

  add_avector(-1.0, xv, bv);
  error = norm2_avector(bv) / norm2_avector(xv);

  (void) printf("Checking avector triangularsolve\n"
		"  (lower=%s, unit=%s, atrans=%s)\n"
		"  Accuracy %g, %sokay\n", (lower ? "tr" : "fl"),
		(unit ? "tr" : "fl"), (atrans ? "tr" : "fl"), error,
		(IS_IN_RANGE(0.0, error, 1.0e-14) ? "" : "    NOT "));
  if (!IS_IN_RANGE(0.0, error, 1.0e-14))
    problems++;

  copy_avector(xv, bv);
  triangulareval_amatrix_avector(lower, unit, atrans, a, bv);

  triangularsolve_amatrix_avector(lower, unit, atrans, a, bv);

  add_avector(-1.0, xv, bv);
  error = norm2_avector(bv) / norm2_avector(xv);

  (void) printf("Checking avector triangulareval/triangularsolve\n"
		"  (lower=%s, unit=%s, atrans=%s):\n"
		"  Accuracy %g, %sokay\n", (lower ? "tr" : "fl"),
		(unit ? "tr" : "fl"), (atrans ? "tr" : "fl"), error,
		(IS_IN_RANGE(0.0, error, 1.0e-14) ? "" : "    NOT "));
  if (!IS_IN_RANGE(0.0, error, 1.0e-14))
    problems++;

  uninit_amatrix(b);
  uninit_amatrix(x);
  uninit_avector(bv);
  uninit_avector(xv);
}

static void
check_lowereval(bool unit, bool atrans, pcamatrix a, bool xtrans)
{
  pamatrix  a2, x, b, b2;
  amatrix   a2tmp, xtmp, btmp;
  real      error;

  a2 = init_amatrix(&a2tmp, a->rows, a->cols);
  if (atrans)
    copy_upper_amatrix(a, unit, a2);
  else
    copy_lower_amatrix(a, unit, a2);

  x = (atrans ? init_amatrix(&xtmp, a->rows, a->rows) :
       init_amatrix(&xtmp, a->cols, a->cols));
  random_amatrix(x);

  b = (atrans ?
       (xtrans ? init_amatrix(&btmp, a->rows, UINT_MAX(a->cols, a->rows)) :
	init_amatrix(&btmp, UINT_MAX(a->rows, a->cols), a->rows)) :
       (xtrans ? init_amatrix(&btmp, a->cols, UINT_MAX(a->cols, a->rows)) :
	init_amatrix(&btmp, UINT_MAX(a->rows, a->cols), a->cols)));

  copy_sub_amatrix(false, x, b);
  triangulareval_amatrix(!atrans, unit, atrans, a, xtrans, b);

  if (xtrans)
    addmul_amatrix(-1.0, false, x, !atrans, a2, b);
  else
    addmul_amatrix(-1.0, atrans, a2, false, x, b);

  uninit_amatrix(a2);
  b2 = (atrans ?
	(xtrans ? init_sub_amatrix(&a2tmp, b, b->rows, 0, a->cols, 0) :
	 init_sub_amatrix(&a2tmp, b, a->cols, 0, b->cols, 0)) :
	(xtrans ? init_sub_amatrix(&a2tmp, b, b->rows, 0, a->rows, 0) :
	 init_sub_amatrix(&a2tmp, b, a->rows, 0, b->cols, 0)));

  error = norm2_amatrix(b2) / norm2_amatrix(x);

  (void) printf("Checking lowereval(unit=%s, atrans=%s, xtrans=%s)\n"
		"  Accuracy %g, %sokay\n", (unit ? "tr" : "fl"),
		(atrans ? "tr" : "fl"), (xtrans ? "tr" : "fl"), error,
		(error < tolerance ? "" : "    NOT "));
  if (error >= tolerance)
    problems++;

  uninit_amatrix(b2);
  uninit_amatrix(b);
  uninit_amatrix(x);
}

static void
check_uppereval(bool unit, bool atrans, pcamatrix a, bool xtrans)
{
  pamatrix  a2, x, b, b2;
  amatrix   a2tmp, xtmp, btmp;
  real      error;

  a2 = init_amatrix(&a2tmp, a->rows, a->cols);
  if (atrans)
    copy_lower_amatrix(a, unit, a2);
  else
    copy_upper_amatrix(a, unit, a2);

  x = (atrans ? init_amatrix(&xtmp, a->rows, a->rows) :
       init_amatrix(&xtmp, a->cols, a->cols));
  random_amatrix(x);

  b = (atrans ?
       (xtrans ? init_amatrix(&btmp, a->rows, UINT_MAX(a->cols, a->rows)) :
	init_amatrix(&btmp, UINT_MAX(a->rows, a->cols), a->rows)) :
       (xtrans ? init_amatrix(&btmp, a->cols, UINT_MAX(a->cols, a->rows)) :
	init_amatrix(&btmp, UINT_MAX(a->rows, a->cols), a->cols)));

  copy_sub_amatrix(false, x, b);
  triangulareval_amatrix(atrans, unit, atrans, a, xtrans, b);

  if (xtrans)
    addmul_amatrix(-1.0, false, x, !atrans, a2, b);
  else
    addmul_amatrix(-1.0, atrans, a2, false, x, b);

  uninit_amatrix(a2);
  b2 = (atrans ?
	(xtrans ? init_sub_amatrix(&a2tmp, b, b->rows, 0, a->cols, 0) :
	 init_sub_amatrix(&a2tmp, b, a->cols, 0, b->cols, 0)) :
	(xtrans ? init_sub_amatrix(&a2tmp, b, b->rows, 0, a->rows, 0) :
	 init_sub_amatrix(&a2tmp, b, a->rows, 0, b->cols, 0)));

  error = norm2_amatrix(b2) / norm2_amatrix(x);

  (void) printf("Checking uppereval(unit=%s, atrans=%s, xtrans=%s)\n"
		"  Accuracy %g, %sokay\n", (unit ? "tr" : "fl"),
		(atrans ? "tr" : "fl"), (xtrans ? "tr" : "fl"), error,
		(error < tolerance ? "" : "    NOT "));
  if (error >= tolerance)
    problems++;

  uninit_amatrix(b2);
  uninit_amatrix(b);
  uninit_amatrix(x);
}

static void
check_triangularaddmul(bool alower, bool atrans, bool blower, bool btrans)
{
  pamatrix  a, b, at, bt, x;
  amatrix   atmp, btmp, attmp, bttmp, xtmp;
  real      error;
  uint      dim1, dim2, dim3;

  dim1 = 100;
  dim2 = 80;
  dim3 = 90;

  a =
    (atrans ? init_amatrix(&atmp, dim2, dim1) :
     init_amatrix(&atmp, dim1, dim2));
  random_amatrix(a);

  b =
    (btrans ? init_amatrix(&btmp, dim3, dim2) :
     init_amatrix(&btmp, dim2, dim3));
  random_amatrix(b);

  at = init_amatrix(&attmp, a->rows, a->cols);
  if (alower)
    copy_lower_amatrix(a, false, at);
  else
    copy_upper_amatrix(a, false, at);

  bt = init_amatrix(&bttmp, b->rows, b->cols);
  if (blower)
    copy_lower_amatrix(b, false, bt);
  else
    copy_upper_amatrix(b, false, bt);

  x = init_amatrix(&xtmp, dim1, dim3);
  clear_amatrix(x);

  triangularaddmul_amatrix(1.0, alower, atrans, a, blower, btrans, b, x);

  addmul_amatrix(-1.0, atrans, at, btrans, bt, x);

  error = norm2_amatrix(x);

  (void)
    printf
    ("Checking triangularaddmul(alower=%s, atrans=%s, blower=%s, btrans=%s)\n"
     "  Accuracy %g, %sokay\n", (alower ? "tr" : "fl"),
     (atrans ? "tr" : "fl"), (blower ? "tr" : "fl"), (btrans ? "tr" : "fl"),
     error, (error < tolerance ? "" : "    NOT "));
  if (error >= tolerance)
    problems++;

  uninit_amatrix(x);
  uninit_amatrix(bt);
  uninit_amatrix(at);
  uninit_amatrix(b);
  uninit_amatrix(a);

  dim1 = 70;
  dim2 = 80;
  dim3 = 90;

  a =
    (atrans ? init_amatrix(&atmp, dim2, dim1) :
     init_amatrix(&atmp, dim1, dim2));
  random_amatrix(a);

  b =
    (btrans ? init_amatrix(&btmp, dim3, dim2) :
     init_amatrix(&btmp, dim2, dim3));
  random_amatrix(b);

  at = init_amatrix(&attmp, a->rows, a->cols);
  if (alower)
    copy_lower_amatrix(a, false, at);
  else
    copy_upper_amatrix(a, false, at);

  bt = init_amatrix(&bttmp, b->rows, b->cols);
  if (blower)
    copy_lower_amatrix(b, false, bt);
  else
    copy_upper_amatrix(b, false, bt);

  x = init_amatrix(&xtmp, dim1, dim3);
  clear_amatrix(x);

  triangularaddmul_amatrix(1.0, alower, atrans, a, blower, btrans, b, x);

  addmul_amatrix(-1.0, atrans, at, btrans, bt, x);

  error = norm2_amatrix(x);

  (void)
    printf
    ("Checking triangularaddmul(alower=%s, atrans=%s, blower=%s, btrans=%s)\n"
     "  Accuracy %g, %sokay\n", (alower ? "tr" : "fl"),
     (atrans ? "tr" : "fl"), (blower ? "tr" : "fl"), (btrans ? "tr" : "fl"),
     error, (error < tolerance ? "" : "    NOT "));
  if (error >= tolerance)
    problems++;

  uninit_amatrix(x);
  uninit_amatrix(bt);
  uninit_amatrix(at);
  uninit_amatrix(b);
  uninit_amatrix(a);
}

static void
set_unit(pamatrix R)
{
  uint      i, m;

  m = UINT_MIN(R->rows, R->cols);
  for (i = 0; i < m; ++i) {
    R->a[i + i * R->ld] = 1.0;
  }
}

int
main()
{
  pamatrix  a, acopy, l, ld, r, q, qr;
  pavector  x, b, tau;
  uint      rows, cols;
  real      error;

  rows = 8;
  cols = 5;

  (void) printf("----------------------------------------\n"
		"Check random %u x %u Cholesky factorization\n", rows, rows);

  (void) printf("Creating random self-adjoint positive definite matrix\n");
  a = new_amatrix(rows, rows);
  random_spd_amatrix(a, 1.0);

  (void) printf("Creating random solution and right-hand side\n");
  x = new_avector(rows);
  random_avector(x);
  b = new_avector(rows);
  clear_avector(b);
  mvm_amatrix_avector(1.0, false, a, x, b);

  (void) printf("Copying matrix\n");
  acopy = new_amatrix(rows, rows);
  copy_amatrix(false, a, acopy);

  (void) printf("Computing Cholesky factorization\n");
  choldecomp_amatrix(a);

  (void) printf("Solving\n");
  cholsolve_amatrix_avector(a, b);
  add_avector(-1.0, x, b);
  error = norm2_avector(b);
  (void) printf("  Accuracy %g, %sokay\n", error,
		(error < tolerance ? "" : "    NOT "));
  if (error >= tolerance)
    problems++;

  (void) printf("Building triangular factor\n");
  l = new_amatrix(rows, rows);
  copy_lower_amatrix(a, false, l);

  (void) printf("Checking factorization\n");
  copy_amatrix(false, acopy, a);
  addmul_amatrix(-1.0, false, l, true, l, a);
  error = normfrob_amatrix(a);
  (void) printf("  Accuracy %g, %sokay\n", error,
		(error < tolerance ? "" : "    NOT "));
  if (error >= tolerance)
    problems++;

  del_amatrix(l);
  del_amatrix(acopy);
  del_avector(b);
  del_avector(x);
  del_amatrix(a);

  (void) printf("----------------------------------------\n"
		"Check random %u x %u LDL^T factorization\n", rows, rows);

  (void) printf("Creating random self-adjoint positive definite matrix\n");
  a = new_amatrix(rows, rows);
  random_spd_amatrix(a, 1.0);

  (void) printf("Creating random solution and right-hand side\n");
  x = new_avector(rows);
  random_avector(x);
  b = new_avector(rows);
  clear_avector(b);
  mvm_amatrix_avector(1.0, false, a, x, b);

  (void) printf("Copying matrix\n");
  acopy = new_amatrix(rows, rows);
  copy_amatrix(false, a, acopy);

  (void) printf("Computing LDL^T factorization\n");
  ldltdecomp_amatrix(a);

  (void) printf("Solving\n");
  ldltsolve_amatrix_avector(a, b);
  add_avector(-1.0, x, b);
  error = norm2_avector(b);
  (void) printf("  Accuracy %g, %sokay\n", error,
		(error < tolerance ? "" : "    NOT "));
  if (error >= tolerance)
    problems++;

  (void) printf("Building triangular factor\n");
  l = new_amatrix(rows, rows);
  copy_lower_amatrix(a, true, l);

  (void) printf("Multiplying by diagonal\n");
  ld = new_amatrix(rows, rows);
  copy_amatrix(false, l, ld);
  diageval_amatrix(true, a, true, ld);

  (void) printf("Checking factorization\n");
  copy_amatrix(false, acopy, a);
  addmul_amatrix(-1.0, false, ld, true, l, a);
  error = normfrob_amatrix(a);
  (void) printf("  Accuracy %g, %sokay\n", error,
		(error < tolerance ? "" : "    NOT "));
  if (error >= tolerance)
    problems++;

  del_amatrix(ld);
  del_amatrix(l);
  del_amatrix(acopy);
  del_avector(b);
  del_avector(x);
  del_amatrix(a);

  (void) printf("----------------------------------------\n"
		"Check random %u x %u LR factorization\n", rows, rows);

  (void) printf("Creating random invertible matrix\n");
  a = new_amatrix(rows, rows);
  random_invertible_amatrix(a, 1.0);

  (void) printf("Creating random solution and right-hand side\n");
  x = new_avector(rows);
  random_avector(x);
  b = new_avector(rows);
  clear_avector(b);
  mvm_amatrix_avector(1.0, false, a, x, b);

  (void) printf("Copying matrix\n");
  acopy = new_amatrix(rows, rows);
  copy_amatrix(false, a, acopy);

  (void) printf("Computing LR factorization\n");
  lrdecomp_amatrix(a);

  (void) printf("Solving\n");
  lrsolve_amatrix_avector(a, b);
  add_avector(-1.0, x, b);
  error = norm2_avector(b);
  (void) printf("  Accuracy %g, %sokay\n", error,
		(error < tolerance ? "" : "    NOT "));
  if (error >= tolerance)
    problems++;

  (void) printf("Building triangular factors\n");
  l = new_amatrix(rows, rows);
  r = new_amatrix(rows, rows);
  copy_lower_amatrix(a, true, l);
  copy_upper_amatrix(a, false, r);

  (void) printf("Checking factorization\n");
  copy_amatrix(false, acopy, a);
  addmul_amatrix(-1.0, false, l, false, r, a);
  error = normfrob_amatrix(a);
  (void) printf("  Accuracy %g, %sokay\n", error,
		(error < tolerance ? "" : "    NOT "));
  if (error >= tolerance)
    problems++;

  /* Check forward/backward substitution */
  (void) printf("----------------------------------------\n");
  check_triangularsolve(true, true, false, l, false);
  check_triangularsolve(true, true, false, l, true);
  check_triangularsolve(true, true, true, l, false);
  check_triangularsolve(true, true, true, l, true);
  check_triangularsolve(true, false, false, l, false);
  check_triangularsolve(true, false, false, l, true);
  check_triangularsolve(true, false, true, l, false);
  check_triangularsolve(true, false, true, l, true);

  check_triangularsolve(false, false, false, r, false);
  check_triangularsolve(false, false, false, r, true);
  check_triangularsolve(false, false, true, r, false);
  check_triangularsolve(false, false, true, r, true);
  set_unit(r);
  check_triangularsolve(false, true, false, r, false);
  check_triangularsolve(false, true, false, r, true);
  check_triangularsolve(false, true, true, r, false);
  check_triangularsolve(false, true, true, r, true);

  /* Intermediate clean-up */
  (void) printf("Cleaning up\n");
  del_amatrix(r);
  del_amatrix(l);

  /* Check triangular matrix multiplication */
  (void) printf("----------------------------------------\n");
  check_triangularaddmul(false, false, false, false);
  check_triangularaddmul(true, false, false, false);
  check_triangularaddmul(false, true, false, false);
  check_triangularaddmul(true, true, false, false);

  check_triangularaddmul(false, false, true, false);
  check_triangularaddmul(true, false, true, false);
  check_triangularaddmul(false, true, true, false);
  check_triangularaddmul(true, true, true, false);

  check_triangularaddmul(false, false, false, true);
  check_triangularaddmul(true, false, false, true);
  check_triangularaddmul(false, true, false, true);
  check_triangularaddmul(true, true, false, true);

  check_triangularaddmul(false, false, true, true);
  check_triangularaddmul(true, false, true, true);
  check_triangularaddmul(false, true, true, true);
  check_triangularaddmul(true, true, true, true);

  /* Checking QR factorization */
  (void) printf("----------------------------------------\n"
		"Check square QR factorization\n");
  copy_amatrix(false, acopy, a);
  random_avector(x);
  clear_avector(b);
  mvm_amatrix_avector(1.0, false, a, x, b);

  tau = new_avector(rows);
  qrdecomp_amatrix(a, tau);

  (void) printf("Solving\n");
  qrsolve_amatrix_avector(a, tau, b);
  add_avector(-1.0, x, b);
  error = norm2_avector(b);
  (void) printf("  Accuracy %g, %sokay\n", error,
		(error < tolerance ? "" : "    NOT "));
  if (error >= tolerance)
    problems++;

  q = new_amatrix(rows, rows);
  qrexpand_amatrix(a, tau, q);

  (void) printf("Checking factorization\n");
  r = new_amatrix(rows, rows);
  copy_upper_amatrix(a, false, r);
  qr = new_amatrix(rows, rows);
  copy_amatrix(false, acopy, qr);
  addmul_amatrix(-1.0, false, q, false, r, qr);
  error = normfrob_amatrix(qr);
  (void) printf("  Accuracy %g, %sokay\n", error,
		(error < tolerance ? "" : "    NOT "));
  if (error >= tolerance)
    problems++;

  (void) printf("Checking inverse\n");
  copy_amatrix(false, acopy, a);
  qrinvert_amatrix(a);

  identity_amatrix(qr);
  addmul_amatrix(-1.0, false, acopy, false, a, qr);
  error = normfrob_amatrix(qr);
  (void) printf("  Accuracy %g, %sokay\n", error,
		(error < tolerance ? "" : "    NOT "));
  if (error >= tolerance)
    problems++;

  /* Intermediate clean-up */
  (void) printf("Cleaning up\n");

  del_amatrix(qr);
  del_amatrix(r);
  del_amatrix(q);
  del_avector(tau);
  del_avector(b);
  del_avector(x);
  del_amatrix(acopy);
  del_amatrix(a);

  /* Check forward/backward evaluation */
  (void) printf("----------------------------------------\n");
  a = new_amatrix(rows, cols);
  random_amatrix(a);
  check_lowereval(false, false, a, false);
  check_lowereval(false, false, a, true);
  check_lowereval(false, true, a, false);
  check_lowereval(false, true, a, true);
  check_uppereval(false, false, a, false);
  check_uppereval(false, false, a, true);
  check_uppereval(false, true, a, false);
  check_uppereval(false, true, a, true);
  set_unit(a);
  check_lowereval(true, false, a, false);
  check_lowereval(true, false, a, true);
  check_lowereval(true, true, a, false);
  check_lowereval(true, true, a, true);
  check_uppereval(true, false, a, false);
  check_uppereval(true, false, a, true);
  check_uppereval(true, true, a, false);
  check_uppereval(true, true, a, true);
  (void) printf("Cleaning up\n");
  del_amatrix(a);
  a = new_amatrix(cols, rows);
  random_amatrix(a);
  check_lowereval(false, false, a, false);
  check_lowereval(false, false, a, true);
  check_lowereval(false, true, a, false);
  check_lowereval(false, true, a, true);
  check_uppereval(false, false, a, false);
  check_uppereval(false, false, a, true);
  check_uppereval(false, true, a, false);
  check_uppereval(false, true, a, true);
  set_unit(a);
  check_lowereval(true, false, a, false);
  check_lowereval(true, false, a, true);
  check_lowereval(true, true, a, false);
  check_lowereval(true, true, a, true);
  check_uppereval(true, false, a, false);
  check_uppereval(true, false, a, true);
  check_uppereval(true, true, a, false);
  check_uppereval(true, true, a, true);
  (void) printf("Cleaning up\n");
  del_amatrix(a);

  /* Test random QR factorizations */
  (void) printf("----------------------------------------\n"
		"Check full rectangular QR factorization\n");
  a = new_amatrix(rows, cols);
  random_amatrix(a);
  acopy = new_amatrix(rows, cols);
  copy_amatrix(false, a, acopy);

  tau = new_avector(rows);
  qrdecomp_amatrix(a, tau);
  q = new_amatrix(rows, rows);
  r = new_amatrix(rows, cols);
  qrexpand_amatrix(a, tau, q);
  copy_upper_amatrix(a, false, r);
  qr = new_amatrix(rows, cols);
  copy_amatrix(false, acopy, qr);
  addmul_amatrix(-1.0, false, q, false, r, qr);
  error = normfrob_amatrix(qr);
  printf("  Accuracy %g, %sokay\n", error,
	 (error < tolerance ? "" : "    NOT "));
  if (error >= tolerance)
    problems++;

  del_amatrix(r);
  del_amatrix(q);

  (void) printf("Check skinny QR factorization\n");
  q = new_amatrix(rows, cols);
  r = new_amatrix(cols, cols);
  qrexpand_amatrix(a, tau, q);
  copy_upper_amatrix(a, false, r);
  copy_amatrix(false, acopy, qr);
  addmul_amatrix(-1.0, false, q, false, r, qr);
  error = normfrob_amatrix(qr);
  printf("  Accuracy %g, %sokay\n", error,
	 (error < tolerance ? "" : "    NOT "));
  if (error >= tolerance)
    problems++;

  /* Final clean-up */
  (void) printf("Cleaning up\n");
  del_amatrix(qr);
  del_amatrix(r);
  del_amatrix(q);
  del_avector(tau);
  del_amatrix(acopy);
  del_amatrix(a);

  (void) printf("----------------------------------------\n"
		"  %u matrices and\n"
		"  %u vectors still active\n"
		"  %u errors found\n", getactives_amatrix(),
		getactives_avector(), problems);

  return problems;
}
