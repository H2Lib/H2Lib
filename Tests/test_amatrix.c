
#include <stdio.h>

#include "amatrix.h"
#include "factorizations.h"
#include "settings.h"

static uint problems = 0;

#ifdef USE_FLOAT
static const real tolerance = 5.0e-5;
#else
static const real tolerance = 1.0e-12;
#endif

#define IS_IN_RANGE(a, b, c) (((a) <= (b)) && ((b) <= (c)))

#ifdef USE_COMPLEX
static field alpha = 1.0 + 1.0 * I;
#else
static field alpha = 1.0;
#endif

static void
check_diagsolve(bool atrans, pamatrix a, bool xtrans)
{
  uint      n = a->rows;
  amatrix   xtmp, btmp;
  pamatrix  X, B;
  avector   xvtmp, bvtmp;
  pavector  x, b;
  real      error;
  uint      rhs = UINT_MAX(2, n / 2);

  assert(a->cols == n);

  x = init_avector(&xvtmp, n);
  b = init_avector(&bvtmp, n);
  if (xtrans) {
    X = init_amatrix(&xtmp, rhs, n);
    B = init_amatrix(&btmp, rhs, n);
  }
  else {
    X = init_amatrix(&xtmp, n, rhs);
    B = init_amatrix(&btmp, n, rhs);
  }

  random_avector(x);
  clear_avector(b);
  add_avector(alpha, x, b);
  diageval_amatrix_avector(atrans, a, b);

  diagsolve_amatrix_avector(atrans, a, b);
  add_avector(-alpha, x, b);
  error = norm2_avector(b);
  (void) printf("Checking avector diageval/diagsolve\n"
		"  (atrans=%s)\n"
		"  Accuracy %g, %sokay\n", (atrans ? "tr" : "fl"), error,
		(IS_IN_RANGE(0.0, error, tolerance) ? "" : "    NOT "));
  if (!IS_IN_RANGE(0.0, error, tolerance))
    problems++;

  random_amatrix(X);
  clear_amatrix(B);
  add_amatrix(alpha, false, X, B);
  diageval_amatrix(atrans, a, xtrans, B);

  diagsolve_amatrix(atrans, a, xtrans, B);
  add_amatrix(-alpha, false, X, B);
  error = norm2_amatrix(B);
  (void) printf("Checking amatrix diageval/diagsolve\n"
		"  (atrans=%s, xtrans=%s)\n"
		"  Accuracy %g, %sokay\n", (atrans ? "tr" : "fl"),
		(xtrans ? "tr" : "fl"), error,
		(IS_IN_RANGE(0.0, error, tolerance) ? "" : "    NOT "));
  if (!IS_IN_RANGE(0.0, error, tolerance))
    problems++;

  uninit_avector(x);
  uninit_avector(b);
  uninit_amatrix(X);
  uninit_amatrix(B);
}

static void
check_triangularsolve(bool lower, bool unit, bool atrans,
		      pcamatrix a, bool xtrans)
{
  uint      n = a->rows;
  amatrix   xtmp, btmp;
  pamatrix  x, b;
  avector   xvtmp;
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
    addmul_amatrix(CONJ(alpha), false, x, !atrans, a, b);
  else
    addmul_amatrix(alpha, atrans, a, false, x, b);

  triangularsolve_amatrix(lower, unit, atrans, a, xtrans, b);

  add_amatrix(xtrans ? -CONJ(alpha) : -alpha, false, x, b);
  error = norm2_amatrix(b) / norm2_amatrix(x);

  clear_amatrix(b);
  add_amatrix(alpha, false, x, b);

  triangulareval_amatrix(lower, unit, atrans, a, xtrans, b);
  triangularsolve_amatrix(lower, unit, atrans, a, xtrans, b);

  add_amatrix(-alpha, false, x, b);

  error = norm2_amatrix(b) / norm2_amatrix(x);

  (void) printf("Checking amatrix triangulareval/triangularsolve\n"
		"  (lower=%s, unit=%s, atrans=%s, xtrans=%s):\n"
		"  Accuracy %g, %sokay\n", (lower ? "tr" : "fl"),
		(unit ? "tr" : "fl"), (atrans ? "tr" : "fl"),
		(xtrans ? "tr" : "fl"), error,
		(IS_IN_RANGE(0.0, error, tolerance) ? "" : "    NOT "));
  if (!IS_IN_RANGE(0.0, error, tolerance))
    problems++;

  /*
   * avector
   */

  xv = init_avector(&xvtmp, n);
  random_avector(xv);

  bv = new_zero_avector(n);

  if (atrans) {
    addevaltrans_amatrix_avector(alpha, a, xv, bv);
  }
  else {
    addeval_amatrix_avector(alpha, a, xv, bv);
  }

  triangularsolve_amatrix_avector(lower, unit, atrans, a, bv);

  add_avector(-alpha, xv, bv);
  error = norm2_avector(bv) / norm2_avector(xv);

  (void) printf("Checking avector triangularsolve\n"
		"  (lower=%s, unit=%s, atrans=%s)\n"
		"  Accuracy %g, %sokay\n", (lower ? "tr" : "fl"),
		(unit ? "tr" : "fl"), (atrans ? "tr" : "fl"), error,
		(IS_IN_RANGE(0.0, error, tolerance) ? "" : "    NOT "));
  if (!IS_IN_RANGE(0.0, error, tolerance))
    problems++;

  clear_avector(bv);
  add_avector(alpha, xv, bv);
  triangulareval_amatrix_avector(lower, unit, atrans, a, bv);
  triangularsolve_amatrix_avector(lower, unit, atrans, a, bv);

  add_avector(-alpha, xv, bv);
  error = norm2_avector(bv) / norm2_avector(xv);

  (void) printf("Checking avector triangulareval/triangularsolve\n"
		"  (lower=%s, unit=%s, atrans=%s):\n"
		"  Accuracy %g, %sokay\n", (lower ? "tr" : "fl"),
		(unit ? "tr" : "fl"), (atrans ? "tr" : "fl"), error,
		(IS_IN_RANGE(0.0, error, tolerance) ? "" : "    NOT "));
  if (!IS_IN_RANGE(0.0, error, tolerance))
    problems++;

  uninit_amatrix(b);
  uninit_amatrix(x);
  uninit_avector(bv);
  uninit_avector(xv);
}

static void
check_lowereval(bool unit, bool atrans, pcamatrix a, bool xtrans)
{
  pamatrix  a2, X, B, B2;
  amatrix   a2tmp, Xtmp, Btmp, B2tmp;
  pavector  x, b, b2;
  avector   xtmp, btmp, b2tmp;
  uint      rows[] = {
    a->rows, a->rows, a->rows / 2
  };
  uint      cols[] = {
    a->cols, a->cols / 2, a->cols
  };
  real      error;
  uint      i;

  for (i = 0; i < 3; ++i) {
    a2 = init_amatrix(&a2tmp, rows[i], cols[i]);

    x = init_avector(&xtmp, atrans ? a2->rows : a2->cols);
    b = init_avector(&btmp, UINT_MAX(a2->rows, a2->cols));
    b2 = init_sub_avector(&b2tmp, b, atrans ? a2->cols : a2->rows, 0);

    X = (atrans ? init_amatrix(&Xtmp, a2->rows, a2->rows) :
	 init_amatrix(&Xtmp, a2->cols, a2->cols));

    B = (atrans ?
	 (xtrans ? init_amatrix(&Btmp, X->rows, UINT_MAX(a->cols, a->rows)) :
	  init_amatrix(&Btmp, UINT_MAX(a->rows, a->cols), X->cols)) :
	 (xtrans ? init_amatrix(&Btmp, X->rows, UINT_MAX(a->cols, a->rows)) :
	  init_amatrix(&Btmp, UINT_MAX(a->rows, a->cols), X->cols)));

    B2 = (atrans ?
	  (xtrans ? init_sub_amatrix(&B2tmp, B, B->rows, 0, a2->cols, 0) :
	   init_sub_amatrix(&B2tmp, B, a2->cols, 0, B->cols, 0)) :
	  (xtrans ? init_sub_amatrix(&B2tmp, B, B->rows, 0, a2->rows, 0) :
	   init_sub_amatrix(&B2tmp, B, a2->rows, 0, B->cols, 0)));

    if (atrans) {
      copy_upper_amatrix(a, unit, a2);
    }
    else {
      copy_lower_amatrix(a, unit, a2);
    }

    /****************************************************
     * avector
     ****************************************************/

    random_avector(x);
    clear_avector(b);
    add_avector(alpha, x, b);
    triangulareval_amatrix_avector(!atrans, unit, atrans, a2, b);
    mvm_amatrix_avector(-alpha, atrans, a2, x, b);

    error = norm2_avector(b2) / norm2_avector(x);

    (void) printf("Checking lowereval (avector)(%d x %d * %d = %d)\n"
		  "  (unit=%s, atrans=%s)\n"
		  "  Accuracy %g, %sokay\n", atrans ? a2->cols : a2->rows,
		  atrans ? a2->rows : a2->cols, x->dim, b2->dim,
		  (unit ? "tr" : "fl"), (atrans ? "tr" : "fl"), error,
		  (error < tolerance ? "" : "    NOT "));
    if (error >= tolerance)
      problems++;

    uninit_avector(x);
    uninit_avector(b);
    uninit_avector(b2);

    /****************************************************
     * amatrix
     ****************************************************/

    random_amatrix(X);

    clear_amatrix(B);
    add_amatrix(alpha, false, X, B);
    triangulareval_amatrix(!atrans, unit, atrans, a2, xtrans, B);

    if (xtrans) {
      addmul_amatrix(-alpha, false, X, !atrans, a2, B);
    }
    else {
      addmul_amatrix(-alpha, atrans, a2, false, X, B);
    }

    error = norm2_amatrix(B2) / norm2_amatrix(X);

    (void)
      printf("Checking lowereval (amatrix)(%d x %d * %d x %d = %d x %d)\n"
	     "  (unit=%s, atrans=%s, xtrans=%s)\n" "  Accuracy %g, %sokay\n",
	     atrans ? a2->cols : a2->rows, atrans ? a2->rows : a2->cols,
	     xtrans ? X->cols : X->rows, xtrans ? X->rows : X->cols,
	     xtrans ? B2->cols : B2->rows, xtrans ? B2->rows : B2->cols,
	     (unit ? "tr" : "fl"), (atrans ? "tr" : "fl"),
	     (xtrans ? "tr" : "fl"), error,
	     (error < tolerance ? "" : "    NOT "));
    if (error >= tolerance)
      problems++;

    uninit_amatrix(B2);
    uninit_amatrix(B);
    uninit_amatrix(X);
    uninit_amatrix(a2);
  }
}

static void
check_uppereval(bool unit, bool atrans, pcamatrix a, bool xtrans)
{
  pamatrix  a2, X, B, B2;
  amatrix   a2tmp, Xtmp, Btmp, B2tmp;
  pavector  x, b, b2;
  avector   xtmp, btmp, b2tmp;
  uint      rows[] = {
    a->rows, a->rows, a->rows / 2
  };
  uint      cols[] = {
    a->cols, a->cols / 2, a->cols
  };
  real      error;
  uint      i;

  for (i = 0; i < 3; ++i) {
    a2 = init_amatrix(&a2tmp, rows[i], cols[i]);

    x = init_avector(&xtmp, atrans ? a2->rows : a2->cols);
    b = init_avector(&btmp, UINT_MAX(a2->rows, a2->cols));
    b2 = init_sub_avector(&b2tmp, b, atrans ? a2->cols : a2->rows, 0);

    X = (atrans ? init_amatrix(&Xtmp, a2->rows, a2->rows) :
	 init_amatrix(&Xtmp, a2->cols, a2->cols));

    B = (atrans ?
	 (xtrans ? init_amatrix(&Btmp, X->rows, UINT_MAX(a->cols, a->rows)) :
	  init_amatrix(&Btmp, UINT_MAX(a->rows, a->cols), X->cols)) :
	 (xtrans ? init_amatrix(&Btmp, X->rows, UINT_MAX(a->cols, a->rows)) :
	  init_amatrix(&Btmp, UINT_MAX(a->rows, a->cols), X->cols)));

    B2 = (atrans ?
	  (xtrans ? init_sub_amatrix(&B2tmp, B, B->rows, 0, a2->cols, 0) :
	   init_sub_amatrix(&B2tmp, B, a2->cols, 0, B->cols, 0)) :
	  (xtrans ? init_sub_amatrix(&B2tmp, B, B->rows, 0, a2->rows, 0) :
	   init_sub_amatrix(&B2tmp, B, a2->rows, 0, B->cols, 0)));

    if (atrans) {
      copy_lower_amatrix(a, unit, a2);
    }
    else {
      copy_upper_amatrix(a, unit, a2);
    }

    /****************************************************
     * avector
     ****************************************************/

    random_avector(x);
    clear_avector(b);
    add_avector(alpha, x, b);
    triangulareval_amatrix_avector(atrans, unit, atrans, a2, b);
    mvm_amatrix_avector(-alpha, atrans, a2, x, b);

    error = norm2_avector(b2) / norm2_avector(x);

    (void) printf("Checking uppereval (avector)(%d x %d * %d = %d)\n"
		  "  (unit=%s, atrans=%s)\n"
		  "  Accuracy %g, %sokay\n", atrans ? a2->cols : a2->rows,
		  atrans ? a2->rows : a2->cols, x->dim, b2->dim,
		  (unit ? "tr" : "fl"), (atrans ? "tr" : "fl"), error,
		  (error < tolerance ? "" : "    NOT "));
    if (error >= tolerance)
      problems++;

    uninit_avector(x);
    uninit_avector(b);
    uninit_avector(b2);

    /****************************************************
     * amatrix
     ****************************************************/

    random_amatrix(X);
    clear_amatrix(B);
    add_amatrix(alpha, false, X, B);
    triangulareval_amatrix(atrans, unit, atrans, a2, xtrans, B);

    if (xtrans) {
      addmul_amatrix(-alpha, false, X, !atrans, a2, B);
    }
    else {
      addmul_amatrix(-alpha, atrans, a2, false, X, B);
    }

    error = norm2_amatrix(B2) / norm2_amatrix(X);

    (void)
      printf("Checking uppereval (amatrix)(%d x %d * %d x %d = %d x %d)\n"
	     "  (unit=%s, atrans=%s, xtrans=%s)\n" "  Accuracy %g, %sokay\n",
	     atrans ? a2->cols : a2->rows, atrans ? a2->rows : a2->cols,
	     xtrans ? X->cols : X->rows, xtrans ? X->rows : X->cols,
	     xtrans ? B2->cols : B2->rows, xtrans ? B2->rows : B2->cols,
	     (unit ? "tr" : "fl"), (atrans ? "tr" : "fl"),
	     (xtrans ? "tr" : "fl"), error,
	     (error < tolerance ? "" : "    NOT "));
    if (error >= tolerance)
      problems++;

    uninit_amatrix(B2);
    uninit_amatrix(B);
    uninit_amatrix(X);
    uninit_amatrix(a2);
  }
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

  triangularaddmul_amatrix(alpha, alower, atrans, a, blower, btrans, b, x);

  addmul_amatrix(-alpha, atrans, at, btrans, bt, x);

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

  triangularaddmul_amatrix(alpha, alower, atrans, a, blower, btrans, b, x);

  addmul_amatrix(-alpha, atrans, at, btrans, bt, x);

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
check_clear_copy_lower(pcamatrix a)
{
  amatrix   a2tmp, a3tmp;
  pamatrix  a2, a3;
  real      error;
  uint      i, n;

  n = UINT_MIN(a->rows, a->cols);

  a2 = init_amatrix(&a2tmp, a->rows, a->cols);
  a3 = init_amatrix(&a3tmp, a->rows, a->cols);

  copy_lower_amatrix(a, true, a2);
  copy_amatrix(false, a, a3);
  clear_upper_amatrix(a3, true);

  error = norm2diff_amatrix(a2, a3) / norm2_amatrix(a2);
  (void) printf("Checking clear_upper(strict) / copy_lower\n"
		"  Accuracy %g, %sokay\n", error,
		(error < tolerance ? "" : "    NOT "));
  if (error >= tolerance)
    problems++;

  copy_upper_amatrix(a, true, a2);
  copy_amatrix(false, a, a3);
  clear_lower_amatrix(a3, true);

  error = norm2diff_amatrix(a2, a3) / norm2_amatrix(a2);
  (void) printf("Checking clear_lower(strict) / copy_upper\n"
		"  Accuracy %g, %sokay\n", error,
		(error < tolerance ? "" : "    NOT "));
  if (error >= tolerance)
    problems++;

  for (i = 0; i < n; ++i) {
    a->a[i + i * a->ld] = 0.0;
  }

  copy_lower_amatrix(a, false, a2);
  copy_amatrix(false, a, a3);
  clear_upper_amatrix(a3, false);

  error = norm2diff_amatrix(a2, a3) / norm2_amatrix(a2);
  (void) printf("Checking clear_upper(non strict) / copy_lower\n"
		"  Accuracy %g, %sokay\n", error,
		(error < tolerance ? "" : "    NOT "));
  if (error >= tolerance)
    problems++;

  copy_upper_amatrix(a, false, a2);
  copy_amatrix(false, a, a3);
  clear_lower_amatrix(a3, false);

  error = norm2diff_amatrix(a2, a3) / norm2_amatrix(a2);
  (void) printf("Checking clear_lower(non strict) / copy_upper\n"
		"  Accuracy %g, %sokay\n", error,
		(error < tolerance ? "" : "    NOT "));
  if (error >= tolerance)
    problems++;

  uninit_amatrix(a2);
  uninit_amatrix(a3);
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
  pamatrix  a, acopy, l, ld, dlt, ldltcopy, r, q, qr, X, B;
  pavector  x, b, tau, lvec, dvec;
  field    *data;
  uint      rows, cols, rhs, i, j;
  uint     *colpiv;
  real      error;

  rows = 8;
  cols = 7;
  rhs = UINT_MAX(2, UINT_MIN(rows, cols) / 2);

  (void) printf("----------------------------------------\n"
		"Call utility functions\n");
  a = new_amatrix(rows, cols);
  random_amatrix(a);

  (void) printf("Print matrix to stdout\n");
  print_amatrix(a);
  (void) printf("Print matrix in matlab format\n");
  print_matlab_amatrix(a);

  del_amatrix(a);

  a = new_amatrix(rows, 0);
  random_amatrix(a);

  (void) printf("Print matrix to stdout\n");
  print_amatrix(a);
  (void) printf("Print matrix in matlab format\n");
  print_matlab_amatrix(a);

  del_amatrix(a);

  a = new_amatrix(0, cols);
  random_amatrix(a);

  (void) printf("Print matrix to stdout\n");
  print_amatrix(a);
  (void) printf("Print matrix in matlab format\n");
  print_matlab_amatrix(a);

  del_amatrix(a);

  a = new_amatrix(1, 0);
  random_amatrix(a);

  (void) printf("Print matrix to stdout\n");
  print_amatrix(a);
  (void) printf("Print matrix in matlab format\n");
  print_matlab_amatrix(a);

  del_amatrix(a);

  a = new_amatrix(1, cols);
  random_amatrix(a);

  (void) printf("Print matrix to stdout\n");
  print_amatrix(a);
  (void) printf("Print matrix in matlab format\n");
  print_matlab_amatrix(a);

  del_amatrix(a);

  (void) printf("----------------------------------------\n"
		"Check random %u x %u Cholesky factorization\n", rows, rows);

  (void) printf("Creating random self-adjoint positive definite matrix\n");
  a = new_amatrix(rows, rows);
  random_spd_amatrix(a, 1.0);

  (void) printf("Creating %d random solutions and right-hand sides\n",
		rhs + 1);
  x = new_avector(rows);
  random_avector(x);
  b = new_avector(rows);
  shrink_avector(b, UINT_MAX(rows - 1, 0));
  resize_avector(b, rows);
  clear_avector(b);
  mvm_amatrix_avector(alpha, false, a, x, b);
  printf("  Size of 'b' is %lu Byte\n", getsize_avector(b));
  printf("  Heapsize of 'b' is %lu Byte\n", getsize_heap_avector(b));
  tau = new_sub_avector(b, b->dim, 0);
  printf("  Size of 'tau' is %lu Byte\n", getsize_avector(tau));
  printf("  Heapsize of 'tau' is %lu Byte\n", getsize_heap_avector(tau));
  print_avector(tau);

  X = new_amatrix(rows, rhs);
  random_amatrix(X);
  B = new_amatrix(rows, rhs);
  clear_amatrix(B);
  addmul_amatrix(alpha, false, a, false, X, B);

  (void) printf("Copying matrix\n");
  acopy = new_amatrix(rows, rows);
  copy_amatrix(false, a, acopy);

  (void) printf("Computing Cholesky factorization\n");
  choldecomp_amatrix(a);

  (void) printf("Solving avector\n");
  cholsolve_amatrix_avector(a, b);
  add_avector(-alpha, x, b);
  error = norm2_avector(b);
  (void) printf("  Accuracy %g, %sokay\n", error,
		(error < tolerance ? "" : "    NOT "));
  if (error >= tolerance)
    problems++;

  (void) printf("Solving amatrix\n");
  cholsolve_amatrix(a, B);
  add_amatrix(-alpha, false, X, B);
  error = norm2_amatrix(B);
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

  del_amatrix(X);
  del_amatrix(B);
  del_amatrix(l);
  del_amatrix(acopy);
  del_avector(b);
  del_avector(x);
  del_avector(tau);
  del_amatrix(a);

  (void) printf("----------------------------------------\n"
		"Check random %u x %u LDL^T factorization\n", rows, rows);

  (void) printf("Creating random self-adjoint positive definite matrix\n");
  a = new_amatrix(rows, rows);
  random_spd_amatrix(a, 1.0);

  (void) printf("Creating %d random solutions and right-hand sides\n",
		rhs + 1);
  x = new_avector(rows);
  random_avector(x);
  b = new_avector(rows);
  clear_avector(b);
  mvm_amatrix_avector(alpha, false, a, x, b);

  X = new_amatrix(rows, rhs);
  random_amatrix(X);
  B = new_amatrix(rows, rhs);
  clear_amatrix(B);
  addmul_amatrix(alpha, false, a, false, X, B);

  (void) printf("Copying matrix\n");
  acopy = new_amatrix(rows, rows);
  copy_amatrix(false, a, acopy);

  (void) printf("Computing LDL^T factorization\n");
  ldltdecomp_amatrix(a);

  (void) printf("Solving avector\n");
  ldltsolve_amatrix_avector(a, b);
  add_avector(-alpha, x, b);
  error = norm2_avector(b);
  (void) printf("  Accuracy %g, %sokay\n", error,
		(error < tolerance ? "" : "    NOT "));
  if (error >= tolerance)
    problems++;

  (void) printf("Solving amatrix\n");
  ldltsolve_amatrix(a, B);
  add_amatrix(-alpha, false, X, B);
  error = norm2_amatrix(B);
  (void) printf("  Accuracy %g, %sokay\n", error,
		(error < tolerance ? "" : "    NOT "));
  if (error >= tolerance)
    problems++;

  (void) printf("Building triangular factor\n");
  l = new_amatrix(rows, rows);
  copy_lower_amatrix(a, true, l);
  ldltcopy = new_amatrix(rows, rows);
  copy_amatrix(false, a, ldltcopy);

  (void) printf("Multiplying by diagonal from right\n");
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

  (void) printf("Multiplying by diagonal from left\n");
  copy_amatrix(false, ldltcopy, a);
  dlt = new_amatrix(rows, rows);
  copy_amatrix(true, l, dlt);
  diageval_amatrix(false, a, false, dlt);

  (void) printf("Checking factorization\n");
  copy_amatrix(false, acopy, a);
  addmul_amatrix(-1.0, false, l, false, dlt, a);
  error = normfrob_amatrix(a);
  (void) printf("  Accuracy %g, %sokay\n", error,
		(error < tolerance ? "" : "    NOT "));
  if (error >= tolerance)
    problems++;

  (void) printf("----------------------------------------\n"
		"Check diageval/diagsolve\n");
  copy_amatrix(false, ldltcopy, a);

  check_diagsolve(false, a, false);
  check_diagsolve(false, a, true);
  check_diagsolve(true, a, false);
  check_diagsolve(true, a, true);

  del_amatrix(X);
  del_amatrix(B);
  del_amatrix(ld);
  del_amatrix(dlt);
  del_amatrix(ldltcopy);
  del_amatrix(l);
  del_amatrix(acopy);
  del_avector(b);
  del_avector(x);
  del_amatrix(a);

  (void) printf("----------------------------------------\n"
		"Check random %u x %u Cholesky factorization\n", rows, rows);

  (void) printf("Creating random invertible matrix\n");
  a = new_amatrix(rows, rows);
  random_spd_amatrix(a, 1.0);

  (void) printf("Creating %d random solutions and right-hand sides\n",
		rhs + 1);
  x = new_avector(rows);
  random_avector(x);
  b = new_avector(rows);

  X = new_amatrix(rows, rhs);
  random_amatrix(X);
  B = new_amatrix(rows, rhs);
  clear_amatrix(B);
  addmul_amatrix(alpha, false, a, false, X, B);

  (void) printf("Copying matrix\n");
  acopy = new_amatrix(rows, rows);
  copy_amatrix(false, a, acopy);

  (void) printf("Computing Cholesky factorization\n");
  choldecomp_amatrix(a);

  (void) printf("Solving avector\n");
  clear_avector(b);
  mvm_amatrix_avector(alpha, false, acopy, x, b);
  cholsolve_amatrix_avector(a, b);
  add_avector(-alpha, x, b);
  error = norm2_avector(b);
  (void) printf("  Accuracy %g, %sokay\n", error,
		(error < tolerance ? "" : "    NOT "));
  if (error >= tolerance)
    problems++;

  (void) printf("Computing Cholesky factorization by block algorithm\n");
  copy_amatrix(false, acopy, a);
  choldecomp_blocks_amatrix(a, 4);

  (void) printf("Solving avector\n");
  clear_avector(b);
  mvm_amatrix_avector(alpha, false, acopy, x, b);
  cholsolve_amatrix_avector(a, b);
  add_avector(-alpha, x, b);
  error = norm2_avector(b);
  (void) printf("  Accuracy %g, %sokay\n", error,
		(error < tolerance ? "" : "    NOT "));
  if (error >= tolerance)
    problems++;

  (void)
    printf("Computing Cholesky factorization by parallelized algorithm\n");
  copy_amatrix(false, acopy, a);
  choldecomp_tasks_amatrix(a, 4);

  (void) printf("Solving avector\n");
  clear_avector(b);
  mvm_amatrix_avector(alpha, false, acopy, x, b);
  cholsolve_amatrix_avector(a, b);
  add_avector(-alpha, x, b);
  error = norm2_avector(b);
  (void) printf("  Accuracy %g, %sokay\n", error,
		(error < tolerance ? "" : "    NOT "));
  if (error >= tolerance)
    problems++;

  (void) printf("Solving amatrix\n");
  cholsolve_amatrix(a, B);
  add_amatrix(-alpha, false, X, B);
  error = norm2_amatrix(B);
  (void) printf("  Accuracy %g, %sokay\n", error,
		(error < tolerance ? "" : "    NOT "));
  if (error >= tolerance)
    problems++;

  del_amatrix(B);
  del_amatrix(X);
  del_avector(b);
  del_avector(x);
  del_amatrix(acopy);
  del_amatrix(a);

  (void) printf("----------------------------------------\n"
		"Check random %u x %u LR factorization\n", rows, rows);

  (void) printf("Creating random invertible matrix\n");
  a = new_amatrix(rows, rows);
  random_invertible_amatrix(a, 1.0);

  (void) printf("Creating %d random solutions and right-hand sides\n",
		rhs + 1);
  x = new_avector(rows);
  random_avector(x);
  b = new_avector(rows);

  X = new_amatrix(rows, rhs);
  random_amatrix(X);
  B = new_amatrix(rows, rhs);
  clear_amatrix(B);
  addmul_amatrix(alpha, false, a, false, X, B);

  (void) printf("Copying matrix\n");
  acopy = new_amatrix(rows, rows);
  copy_amatrix(false, a, acopy);

  (void) printf("Computing LR factorization\n");
  lrdecomp_amatrix(a);

  (void) printf("Solving avector\n");
  clear_avector(b);
  mvm_amatrix_avector(alpha, false, acopy, x, b);
  lrsolve_n_amatrix_avector(a, b);
  add_avector(-alpha, x, b);
  error = norm2_avector(b);
  (void) printf("  Accuracy %g, %sokay\n", error,
		(error < tolerance ? "" : "    NOT "));
  if (error >= tolerance)
    problems++;

  (void) printf("Computing LR factorization by block algorithm\n");
  copy_amatrix(false, acopy, a);
  lrdecomp_blocks_amatrix(a, 4);

  (void) printf("Solving avector\n");
  clear_avector(b);
  mvm_amatrix_avector(alpha, false, acopy, x, b);
  lrsolve_n_amatrix_avector(a, b);
  add_avector(-alpha, x, b);
  error = norm2_avector(b);
  (void) printf("  Accuracy %g, %sokay\n", error,
		(error < tolerance ? "" : "    NOT "));
  if (error >= tolerance)
    problems++;

  (void) printf("Computing LR factorization by parallelized algorithm\n");
  copy_amatrix(false, acopy, a);
  lrdecomp_tasks_amatrix(a, 4);

  (void) printf("Solving avector\n");
  clear_avector(b);
  mvm_amatrix_avector(alpha, false, acopy, x, b);
  lrsolve_n_amatrix_avector(a, b);
  add_avector(-alpha, x, b);
  error = norm2_avector(b);
  (void) printf("  Accuracy %g, %sokay\n", error,
		(error < tolerance ? "" : "    NOT "));
  if (error >= tolerance)
    problems++;

  (void) printf("Solving amatrix\n");
  lrsolve_amatrix(a, B);
  add_amatrix(-alpha, false, X, B);
  error = norm2_amatrix(B);
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
  (void) printf("----------------------------------------\n");
  check_triangularsolve(true, false, false, l, false);
  check_triangularsolve(true, false, false, l, true);
  check_triangularsolve(true, false, true, l, false);
  check_triangularsolve(true, false, true, l, true);
  (void) printf("----------------------------------------\n");
  check_triangularsolve(false, false, false, r, false);
  check_triangularsolve(false, false, false, r, true);
  check_triangularsolve(false, false, true, r, false);
  check_triangularsolve(false, false, true, r, true);
  (void) printf("----------------------------------------\n");
  copy_amatrix(true, l, r);
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
  mvm_amatrix_avector(alpha, false, a, x, b);

  tau = new_avector(rows);
  qrdecomp_amatrix(a, tau);

  (void) printf("Solving\n");
  qrsolve_amatrix_avector(a, tau, b);
  add_avector(-alpha, x, b);
  error = norm2_avector(b);
  (void) printf("  Accuracy %g, %sokay\n", error,
		(error < tolerance ? "" : "    NOT "));
  if (error >= tolerance)
    problems++;

  (void) printf("Checking qreval avector\n");
  clear_avector(b);
  add_avector(alpha, x, b);
  triangulareval_amatrix_avector(false, false, false, a, b);
  qreval_amatrix_avector(false, a, tau, b);
  addeval_amatrix_avector(-alpha, acopy, x, b);
  error = norm2_avector(b);
  (void) printf("  Accuracy %g, %sokay\n", error,
		(error < tolerance ? "" : "    NOT "));
  if (error >= tolerance)
    problems++;

  (void) printf("Checking qreval amatrix\n");
  clear_amatrix(B);
  add_amatrix(alpha, false, X, B);
  triangulareval_amatrix(false, false, false, a, false, B);
  qreval_amatrix(false, a, tau, B);
  addmul_amatrix(-alpha, false, acopy, false, X, B);
  error = norm2_amatrix(B);
  (void) printf("  Accuracy %g, %sokay\n", error,
		(error < tolerance ? "" : "    NOT "));
  if (error >= tolerance)
    problems++;

  (void) printf("Checking qreval amatrix trans\n");
  clear_amatrix(B);
  add_amatrix(alpha, false, X, B);
  qreval_amatrix(true, a, tau, B);
  triangulareval_amatrix(false, false, true, a, false, B);
  addmul_amatrix(-alpha, true, acopy, false, X, B);
  error = norm2_amatrix(B);
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

  del_amatrix(X);
  del_amatrix(B);
  del_amatrix(qr);
  del_amatrix(r);
  del_amatrix(q);
  del_avector(tau);
  del_avector(b);
  del_avector(x);
  del_amatrix(acopy);
  del_amatrix(a);

  /* Check column-pivoted QR factorizations */
  (void) printf("----------------------------------------\n"
		"Check column-pivoted QR factorization\n");
  a = new_amatrix(rows, cols);
  random_amatrix(a);

  colpiv = (uint *) allocmem(sizeof(uint) * cols);
  tau = new_avector(cols);
  acopy = clone_amatrix(a);

  qrdecomp_pivot_amatrix(a, tau, colpiv);

  r = new_amatrix(rows, cols);
  copy_colpiv_amatrix(false, acopy, colpiv, r);
  qreval_amatrix(true, a, tau, r);

  clear_lower_amatrix(a, true);

  add_amatrix(-1.0, false, a, r);
  error = normfrob_amatrix(r);

  (void) printf("  Accuracy %g, %sokay\n", error,
		(error < tolerance ? "" : "    NOT "));
  if (error >= tolerance)
    problems++;

  del_amatrix(r);
  del_amatrix(acopy);
  del_avector(tau);
  freemem(colpiv);
  del_amatrix(a);

  /* Check rank-revealing QR factorization */
  (void) printf("----------------------------------------\n"
		"Check rank-revealing QR factorization\n");
  a = new_amatrix(rows, cols);
  for (j = 0; j < cols; j++)
    for (i = 0; i < rows; i++)
      a->a[i + j * a->ld] = j;

  acopy = clone_amatrix(a);

  colpiv = (uint *) allocmem(sizeof(uint) * cols);
  tau = new_avector(cols);
  j = qrdecomp_rank_amatrix(a, tau, 0, tolerance, colpiv);

  (void) printf("  Detected rank %u, %sokay\n", j,
		(j == 1 ? "" : "    NOT "));
  if (j != 1)
    problems++;

  (void) printf("  Chosen pivot %u, %sokay\n", colpiv[0],
		(colpiv[0] == cols - 1 ? "" : "    NOT "));
  if (colpiv[0] != cols - 1)
    problems++;

  r = new_amatrix(j, cols);
  copy_upper_amatrix(a, false, r);

  q = new_amatrix(rows, j);
  qrexpand_amatrix(a, tau, q);

  clear_amatrix(a);
  copy_colpiv_amatrix(false, acopy, colpiv, a);

  addmul_amatrix(-1.0, false, q, false, r, a);
  error = normfrob_amatrix(a);

  (void) printf("  Accuracy %g, %sokay\n", error,
		(error < tolerance ? "" : "    NOT "));
  if (error >= tolerance)
    problems++;

  del_amatrix(q);
  del_amatrix(r);
  del_avector(tau);
  freemem(colpiv);
  del_amatrix(acopy);
  del_amatrix(a);

  /* Check forward/backward evaluation */
  (void) printf("----------------------------------------\n");
  a = new_amatrix(cols, cols);
  random_amatrix(a);
  check_lowereval(false, false, a, false);
  check_lowereval(false, false, a, true);
  check_lowereval(false, true, a, false);
  check_lowereval(false, true, a, true);
  check_uppereval(false, false, a, false);
  check_uppereval(false, false, a, true);
  check_uppereval(false, true, a, false);
  check_uppereval(false, true, a, true);
  (void) printf("----------------------------------------\n");
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
  (void) printf("----------------------------------------\n");
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
  (void) printf("----------------------------------------\n");
  set_unit(a);
  check_lowereval(true, false, a, false);
  check_lowereval(true, false, a, true);
  check_lowereval(true, true, a, false);
  check_lowereval(true, true, a, true);
  check_uppereval(true, false, a, false);
  check_uppereval(true, false, a, true);
  check_uppereval(true, true, a, false);
  check_uppereval(true, true, a, true);

  (void) printf("----------------------------------------\n");
  check_clear_copy_lower(a);

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

  (void) printf("Cleaning up\n");
  del_amatrix(qr);
  del_amatrix(r);
  del_amatrix(q);
  del_avector(tau);
  del_amatrix(acopy);
  del_amatrix(a);

  (void) printf("----------------------------------------\n");
  (void) printf("Testing amatrix / avector conversion\n");
  data = allocfield(rows * rows);

  a = (pamatrix) allocmem(sizeof(amatrix));
  q = new_amatrix(rows, rows);
  r = new_sub_amatrix(q, rows, 0, rows, 0);

  acopy = new_pointer_amatrix(data, rows, rows);
  x = new_pointer_avector(data, rows * rows);
  (void) init_vec_amatrix(a, x, rows, rows);

  random_spd_amatrix(a, 1.0);
  copy_amatrix(true, a, q);
  add_amatrix(-1.0, true, r, acopy);
  error = normfrob_amatrix(acopy);
  printf("  Accuracy %g, %sokay\n", error,
	 (error < tolerance ? "" : "    NOT "));
  if (error >= tolerance)
    problems++;

  del_amatrix(q);
  del_amatrix(r);
  uninit_amatrix(acopy);
  uninit_amatrix(a);
  uninit_avector(x);
  freemem(data);
  freemem(a);
  freemem(acopy);
  freemem(x);

  (void) printf("----------------------------------------\n");
  (void) printf("Testing dotprod_amatrix\n");

  a = new_amatrix(rows, cols);
  random_amatrix(a);
  error = dotprod_amatrix(a, a) - REAL_SQR(normfrob_amatrix(a));
  printf("  Accuracy %g, %sokay\n", error,
	 (error < tolerance ? "" : "    NOT "));
  if (error >= tolerance)
    problems++;

  del_amatrix(a);

  (void) printf("----------------------------------------\n");
  (void) printf("Testing bidiagmul_amatrix (right)\n");
  lvec = new_avector(cols - 1);
  dvec = new_avector(cols);
  random_avector(lvec);
  random_avector(dvec);
  X = new_zero_amatrix(cols, cols);
  for (i = 0; i < cols; ++i) {
    addentry_amatrix(X, i, i, dvec->v[i]);
  }
  for (i = 0; i < cols - 1; ++i) {
    addentry_amatrix(X, i + 1, i, lvec->v[i]);
  }
  a = new_amatrix(rows, cols);
  random_amatrix(a);
  acopy = clone_amatrix(a);

  bidiagmul_amatrix(alpha, false, a, dvec, lvec);
  addmul_amatrix(-alpha, false, acopy, false, X, a);
  error = normfrob_amatrix(a);
  printf("  Accuracy %g, %sokay\n", error,
	 (error < tolerance ? "" : "    NOT "));
  if (error >= tolerance)
    problems++;

  del_amatrix(a);
  del_amatrix(acopy);
  del_amatrix(X);
  del_avector(lvec);
  del_avector(dvec);

  (void) printf("Testing bidiagmul_amatrix (left)\n");
  lvec = new_avector(rows - 1);
  dvec = new_avector(rows);
  random_avector(lvec);
  random_avector(dvec);
  X = new_zero_amatrix(rows, rows);
  for (i = 0; i < rows; ++i) {
    addentry_amatrix(X, i, i, dvec->v[i]);
  }
  for (i = 0; i < rows - 1; ++i) {
    addentry_amatrix(X, i + 1, i, lvec->v[i]);
  }
  a = new_amatrix(rows, cols);
  random_amatrix(a);
  acopy = clone_amatrix(a);

  bidiagmul_amatrix(alpha, true, a, dvec, lvec);
  addmul_amatrix(-CONJ(alpha), true, X, false, acopy, a);
  error = normfrob_amatrix(a);
  printf("  Accuracy %g, %sokay\n", error,
	 (error < tolerance ? "" : "    NOT "));
  if (error >= tolerance)
    problems++;

  /* Final clean-up */
  del_amatrix(a);
  del_amatrix(acopy);
  del_amatrix(X);
  del_avector(lvec);
  del_avector(dvec);

  (void) printf("----------------------------------------\n"
		"  %u matrices and\n"
		"  %u vectors still active\n"
		"  %u errors found\n", getactives_amatrix(),
		getactives_avector(), problems);

  return problems;
}
