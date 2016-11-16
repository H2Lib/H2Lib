/* ------------------------------------------------------------
 * This is the file "factorizations.c" of the H2Lib package.
 * All rights reserved, Steffen Boerm 2010
 * ------------------------------------------------------------ */

#include "factorizations.h"

#include "settings.h"
#include "basic.h"

#include <math.h>
#include <stdio.h>

/* ------------------------------------------------------------
 * Diagonal matrices
 * ------------------------------------------------------------ */

void
diagsolve_amatrix_avector(bool atrans, pcamatrix a, pavector x)
{
  pcfield   aa = a->a;
  longindex lda = a->ld;
  pfield    xv = x->v;
  uint      n = UINT_MIN(a->rows, a->cols);
  uint      i;

  if (atrans) {
    for (i = 0; i < n; i++)
      xv[i] /= CONJ(aa[i + i * lda]);
  }
  else {
    for (i = 0; i < n; i++)
      xv[i] /= aa[i + i * lda];
  }
}

#ifdef USE_BLAS
void
diagsolve_amatrix(bool atrans, pcamatrix a, bool xtrans, pamatrix x)
{
  uint      n = UINT_MIN(a->rows, a->cols);
  longindex lda = a->ld;
  longindex ldx = x->ld;
  field    *aa = a->a;
  field    *xa = x->a;
  field     alpha;
  uint      i;

  if (xtrans) {
    for (i = 0; i < n; i++) {
      alpha = (atrans ? 1.0 / CONJ(aa[i + i * lda]) : 1.0 / aa[i + i * lda]);
      h2_scal(&x->rows, &alpha, xa + i * ldx, &u_one);
    }
  }
  else {
    for (i = 0; i < n; i++) {
      alpha = (atrans ? 1.0 / CONJ(aa[i + i * lda]) : 1.0 / aa[i + i * lda]);
      h2_scal(&x->cols, &alpha, xa + i, &x->ld);
    }
  }
}
#else
void
diagsolve_amatrix(bool atrans, pcamatrix a, bool xtrans, pamatrix x)
{
  uint      n = UINT_MIN(a->rows, a->cols);
  longindex lda = a->ld;
  longindex ldx = x->ld;
  field    *aa = a->a;
  field    *xa = x->a;
  field     alpha;
  uint      i, j;

  if (xtrans) {
    for (i = 0; i < n; i++) {
      alpha = (atrans ? 1.0 / CONJ(aa[i + i * lda]) : 1.0 / aa[i + i * lda]);
      for (j = 0; j < x->rows; j++)
	xa[j + i * ldx] *= alpha;
    }
  }
  else {
    for (i = 0; i < n; i++) {
      alpha = (atrans ? 1.0 / CONJ(aa[i + i * lda]) : 1.0 / aa[i + i * lda]);
      for (j = 0; j < x->cols; j++)
	xa[i + j * ldx] *= alpha;
    }
  }
}
#endif

void
diageval_amatrix_avector(bool atrans, pcamatrix a, pavector x)
{
  uint      n = UINT_MIN(a->rows, a->cols);
  longindex lda = a->ld;
  field    *aa = a->a;
  field    *xv = x->v;
  field     alpha;
  uint      i;

  assert(x->dim >= a->cols);

  for (i = 0; i < n; i++) {
    alpha = (atrans ? CONJ(aa[i + i * lda]) : aa[i + i * lda]);
    xv[i] *= alpha;
  }
}

#ifdef USE_BLAS
void
diageval_amatrix(bool atrans, pcamatrix a, bool xtrans, pamatrix x)
{
  uint      n = UINT_MIN(a->rows, a->cols);
  longindex lda = a->ld;
  longindex ldx = x->ld;
  field    *aa = a->a;
  field    *xa = x->a;
  field     alpha;
  uint      i;

  if (xtrans) {
    for (i = 0; i < n; i++) {
      alpha = (atrans ? CONJ(aa[i + i * lda]) : aa[i + i * lda]);
      h2_scal(&x->rows, &alpha, xa + i * ldx, &u_one);
    }
  }
  else {
    for (i = 0; i < n; i++) {
      alpha = (atrans ? CONJ(aa[i + i * lda]) : aa[i + i * lda]);
      h2_scal(&x->cols, &alpha, xa + i, &x->ld);
    }
  }
}
#else
void
diageval_amatrix(bool atrans, pcamatrix a, bool xtrans, pamatrix x)
{
  uint      n = UINT_MIN(a->rows, a->cols);
  longindex lda = a->ld;
  longindex ldx = x->ld;
  field    *aa = a->a;
  field    *xa = x->a;
  field     alpha;
  uint      i, j;

  if (xtrans) {
    for (i = 0; i < n; i++) {
      alpha = (atrans ? CONJ(aa[i + i * lda]) : aa[i + i * lda]);
      for (j = 0; j < x->rows; j++)
	xa[j + i * ldx] *= alpha;
    }
  }
  else {
    for (i = 0; i < n; i++) {
      alpha = (atrans ? CONJ(aa[i + i * lda]) : aa[i + i * lda]);
      for (j = 0; j < x->cols; j++)
	xa[i + j * ldx] *= alpha;
    }
  }
}
#endif

#ifdef USE_BLAS
void
diageval_realavector_amatrix(field alpha, bool atrans, pcrealavector a,
			     bool xtrans, pamatrix x)
{
  preal     av = a->v;
  pfield    xa = x->a;
  longindex ldx = x->ld;
  field     beta;
  unsigned  j;

  if (xtrans) {
    assert(x->cols <= a->dim);

    if (a->dim < 1 || x->rows < 1) {	/* Quick exit */
      return;
    }

    for (j = 0; j < x->cols; j++) {
      beta = (atrans ? CONJ(alpha) * av[j] : CONJ(alpha * av[j]));
      h2_scal(&x->rows, &beta, xa + j * ldx, &u_one);
    }
  }
  else {
    assert(x->rows <= a->dim);

    if (a->dim < 1 || x->cols < 1) {	/* Quick exit */
      return;
    }

    for (j = 0; j < x->rows; j++) {
      beta = (atrans ? alpha * CONJ(av[j]) : alpha * av[j]);
      h2_scal(&x->cols, &beta, xa + j, &x->ld);
    }
  }
}
#else
void
diageval_realavector_amatrix(field alpha, bool atrans, pcrealavector a,
			     bool xtrans, pamatrix x)
{
  preal     av = a->v;
  pfield    xa = x->a;
  longindex ldx = x->ld;
  field     beta;
  unsigned  i, j;

  if (xtrans) {
    assert(x->cols <= a->dim);

    if (a->dim < 1 || x->rows < 1) {	/* Quick exit */
      return;
    }

    for (i = 0; i < x->cols; i++) {
      beta = (atrans ? CONJ(alpha) * av[i] : CONJ(alpha * av[i]));
      for (j = 0; j < x->rows; j++) {
	xa[j + i * ldx] *= beta;
      }
    }
  }
  else {
    assert(x->rows <= a->dim);

    if (a->dim < 1 || x->cols < 1) {	/* Quick exit */
      return;
    }

    for (i = 0; i < x->rows; i++) {
      beta = (atrans ? alpha * CONJ(av[i]) : alpha * av[i]);
      for (j = 0; j < x->cols; j++) {
	xa[i + j * ldx] *= beta;
      }
    }
  }
}
#endif

/* ------------------------------------------------------------
 * Triangular matrices
 * ------------------------------------------------------------ */

#ifdef USE_BLAS
static void
lowersolve_amatrix_avector(bool aunit, bool atrans, pcamatrix a, pavector x)
{
  uint      n = UINT_MIN(a->rows, a->cols);

  assert(x->dim >= a->rows);
  assert(x->dim >= a->cols);

  if (atrans) {
    h2_trsm(_h2_left, _h2_lower, _h2_adj, (aunit ? _h2_unit : _h2_nonunit),
	    &n, &u_one, &f_one, a->a, &a->ld, x->v, &x->dim);
  }
  else {
    h2_trsm(_h2_left, _h2_lower, _h2_ntrans, (aunit ? _h2_unit : _h2_nonunit),
	    &n, &u_one, &f_one, a->a, &a->ld, x->v, &x->dim);
  }
}

static void
uppersolve_amatrix_avector(bool aunit, bool atrans, pcamatrix a, pavector x)
{
  uint      n = UINT_MIN(a->rows, a->cols);

  assert(x->dim >= a->rows);
  assert(x->dim >= a->cols);

  if (atrans) {
    h2_trsm(_h2_left, _h2_upper, _h2_adj, (aunit ? _h2_unit : _h2_nonunit),
	    &n, &u_one, &f_one, a->a, &a->ld, x->v, &x->dim);
  }
  else {
    h2_trsm(_h2_left, _h2_upper, _h2_ntrans, (aunit ? _h2_unit : _h2_nonunit),
	    &n, &u_one, &f_one, a->a, &a->ld, x->v, &x->dim);
  }
}
#else
static void
lowersolve_amatrix_avector(bool aunit, bool atrans, pcamatrix a, pavector x)
{
  pcfield   aa = a->a;
  longindex lda = a->ld;
  pfield    xv = x->v;
  uint      n = UINT_MIN(a->rows, a->cols);
  field     newval;
  uint      i, j;

  assert(x->dim >= a->rows);
  assert(x->dim >= a->cols);

  if (atrans) {
    for (j = n; j-- > 0;) {
      newval = (aunit ? xv[j] : (xv[j] /= CONJ(aa[j + j * lda])));
      for (i = 0; i < j; i++)
	xv[i] -= CONJ(aa[j + i * lda]) * newval;
    }
  }
  else {
    for (j = 0; j < n; j++) {
      newval = (aunit ? xv[j] : (xv[j] /= aa[j + j * lda]));
      for (i = j + 1; i < n; i++)
	xv[i] -= aa[i + j * lda] * newval;
    }
  }
}

static void
uppersolve_amatrix_avector(bool aunit, bool atrans, pcamatrix a, pavector x)
{
  pcfield   aa = a->a;
  longindex lda = a->ld;
  pfield    xv = x->v;
  uint      n = UINT_MIN(a->rows, a->cols);
  field     newval;
  uint      i, j;

  assert(x->dim >= a->rows);
  assert(x->dim >= a->cols);

  if (atrans) {
    for (j = 0; j < n; j++) {
      newval = (aunit ? xv[j] : (xv[j] /= CONJ(aa[j + j * lda])));
      for (i = j + 1; i < n; i++)
	xv[i] -= CONJ(aa[j + i * lda]) * newval;
    }
  }
  else {
    for (j = n; j-- > 0;) {
      newval = (aunit ? xv[j] : (xv[j] /= aa[j + j * lda]));
      for (i = 0; i < j; i++)
	xv[i] -= aa[i + j * lda] * newval;
    }
  }
}
#endif

void
triangularsolve_amatrix_avector(bool alower, bool aunit, bool atrans,
				pcamatrix a, pavector x)
{
  if (alower)
    lowersolve_amatrix_avector(aunit, atrans, a, x);
  else
    uppersolve_amatrix_avector(aunit, atrans, a, x);
}

#ifdef USE_BLAS
static void
lowersolve_amatrix(bool aunit, bool atrans, pcamatrix a,
		   bool xtrans, pamatrix x)
{
  uint      n = UINT_MIN(a->rows, a->cols);
  field    *aa = a->a;
  field    *xa = x->a;

  if (atrans) {
    if (xtrans) {
      assert(x->cols >= n);

      h2_trsm(_h2_right, _h2_lower, _h2_ntrans,
	      (aunit ? _h2_unit : _h2_nonunit), &x->rows, &n, &f_one, aa,
	      &a->ld, xa, &x->ld);
    }
    else {
      assert(x->rows >= n);

      h2_trsm(_h2_left, _h2_lower, _h2_adj, (aunit ? _h2_unit : _h2_nonunit),
	      &n, &x->cols, &f_one, aa, &a->ld, xa, &x->ld);
    }
  }
  else {
    if (xtrans) {
      assert(x->cols >= n);

      h2_trsm(_h2_right, _h2_lower, _h2_adj, (aunit ? _h2_unit : _h2_nonunit),
	      &x->rows, &n, &f_one, aa, &a->ld, xa, &x->ld);
    }
    else {
      assert(x->rows >= n);

      h2_trsm(_h2_left, _h2_lower, _h2_ntrans,
	      (aunit ? _h2_unit : _h2_nonunit), &n, &x->cols, &f_one, aa,
	      &a->ld, xa, &x->ld);
    }
  }
}

static void
uppersolve_amatrix(bool aunit, bool atrans, pcamatrix a,
		   bool xtrans, pamatrix x)
{
  uint      n = UINT_MIN(a->rows, a->cols);
  field    *aa = a->a;
  field    *xa = x->a;

  if (atrans) {
    if (xtrans) {
      assert(x->cols >= n);

      h2_trsm(_h2_right, _h2_upper, _h2_ntrans,
	      (aunit ? _h2_unit : _h2_nonunit), &x->rows, &n, &f_one, aa,
	      &a->ld, xa, &x->ld);
    }
    else {
      assert(x->rows >= n);

      h2_trsm(_h2_left, _h2_upper, _h2_adj, (aunit ? _h2_unit : _h2_nonunit),
	      &n, &x->cols, &f_one, aa, &a->ld, xa, &x->ld);
    }
  }
  else {
    if (xtrans) {
      assert(x->cols >= n);

      h2_trsm(_h2_right, _h2_upper, _h2_adj, (aunit ? _h2_unit : _h2_nonunit),
	      &x->rows, &n, &f_one, aa, &a->ld, xa, &x->ld);
    }
    else {
      assert(x->rows >= n);

      h2_trsm(_h2_left, _h2_upper, _h2_ntrans,
	      (aunit ? _h2_unit : _h2_nonunit), &n, &x->cols, &f_one, aa,
	      &a->ld, xa, &x->ld);
    }
  }
}
#else
static void
lowersolve_amatrix(bool aunit, bool atrans, pcamatrix a,
		   bool xtrans, pamatrix x)
{
  uint      n = UINT_MIN(a->rows, a->cols);
  longindex lda = a->ld;
  longindex ldx = x->ld;
  pfield    aa = a->a;
  pfield    xa = x->a;
  uint      i, j, k;
  field     alpha;

  if (atrans) {
    if (xtrans) {
      assert(x->cols >= n);

      for (k = n; k-- > 0;) {
	if (!aunit) {
	  alpha = 1.0 / aa[k + k * lda];
	  for (i = 0; i < x->rows; i++)
	    xa[i + k * ldx] *= alpha;
	}
	for (i = 0; i < x->rows; i++)
	  for (j = 0; j < k; j++)
	    xa[i + j * ldx] -= xa[i + k * ldx] * aa[k + j * lda];
      }
    }
    else {
      assert(x->rows >= n);

      for (k = n; k-- > 0;) {
	if (!aunit) {
	  alpha = 1.0 / CONJ(aa[k + k * lda]);
	  for (j = 0; j < x->cols; j++)
	    xa[k + j * ldx] *= alpha;
	}
	for (i = 0; i < k; i++)
	  for (j = 0; j < x->cols; j++)
	    xa[i + j * ldx] -= CONJ(aa[k + i * lda]) * xa[k + j * ldx];
      }
    }
  }
  else {
    if (xtrans) {
      assert(x->cols >= n);

      for (k = 0; k < n; k++) {
	if (!aunit) {
	  alpha = 1.0 / CONJ(aa[k + k * lda]);
	  for (i = 0; i < x->rows; i++)
	    xa[i + k * ldx] *= alpha;
	}
	for (i = 0; i < x->rows; i++)
	  for (j = k + 1; j < n; j++)
	    xa[i + j * ldx] -= xa[i + k * ldx] * CONJ(aa[j + k * lda]);
      }
    }
    else {
      assert(x->rows >= n);

      for (k = 0; k < n; k++) {
	if (!aunit) {
	  alpha = 1.0 / aa[k + k * lda];
	  for (j = 0; j < x->cols; j++)
	    xa[k + j * ldx] *= alpha;
	}
	for (i = k + 1; i < n; i++)
	  for (j = 0; j < x->cols; j++)
	    xa[i + j * ldx] -= aa[i + k * lda] * xa[k + j * ldx];
      }
    }
  }
}

static void
uppersolve_amatrix(bool aunit, bool atrans, pcamatrix a,
		   bool xtrans, pamatrix x)
{
  uint      n = UINT_MIN(a->rows, a->cols);
  longindex lda = a->ld;
  longindex ldx = x->ld;
  pfield    aa = a->a;
  pfield    xa = x->a;
  uint      i, j, k;
  field     alpha;

  if (atrans) {
    if (xtrans) {
      assert(x->cols >= n);

      for (k = 0; k < n; k++) {
	if (!aunit) {
	  alpha = 1.0 / aa[k + k * lda];
	  for (i = 0; i < x->rows; i++)
	    xa[i + k * ldx] *= alpha;
	}
	for (i = 0; i < x->rows; i++)
	  for (j = k + 1; j < n; j++)
	    xa[i + j * ldx] -= xa[i + k * ldx] * aa[k + j * lda];
      }
    }
    else {
      assert(x->rows >= n);

      for (k = 0; k < n; k++) {
	if (!aunit) {
	  alpha = 1.0 / CONJ(aa[k + k * lda]);
	  for (j = 0; j < x->cols; j++)
	    xa[k + j * ldx] *= alpha;
	}
	for (i = k + 1; i < n; i++)
	  for (j = 0; j < x->cols; j++)
	    xa[i + j * ldx] -= CONJ(aa[k + i * lda]) * xa[k + j * ldx];
      }
    }
  }
  else {
    if (xtrans) {
      assert(x->cols >= n);

      for (k = n; k-- > 0;) {
	if (!aunit) {
	  alpha = 1.0 / CONJ(aa[k + k * lda]);
	  for (i = 0; i < x->rows; i++)
	    xa[i + k * ldx] *= alpha;
	}
	for (i = 0; i < x->rows; i++)
	  for (j = 0; j < k; j++)
	    xa[i + j * ldx] -= xa[i + k * ldx] * CONJ(aa[j + k * lda]);
      }
    }
    else {
      assert(x->rows >= n);

      for (k = n; k-- > 0;) {
	if (!aunit) {
	  alpha = 1.0 / aa[k + k * lda];
	  for (j = 0; j < x->cols; j++)
	    xa[k + j * ldx] *= alpha;
	}
	for (i = 0; i < k; i++)
	  for (j = 0; j < x->cols; j++)
	    xa[i + j * ldx] -= aa[i + k * lda] * xa[k + j * ldx];
      }
    }
  }
}
#endif

void
triangularsolve_amatrix(bool alower, bool aunit, bool atrans, pcamatrix a,
			bool xtrans, pamatrix x)
{
  if (alower)
    lowersolve_amatrix(aunit, atrans, a, xtrans, x);
  else
    uppersolve_amatrix(aunit, atrans, a, xtrans, x);
}

#ifdef USE_BLAS
static void
lowereval_amatrix_avector(bool aunit, bool atrans, pcamatrix a, pavector x)
{
  field    *aa = a->a;
  field    *xv = x->v;
  uint      n = UINT_MIN(a->rows, a->cols);
  uint      n1, i;

  assert(x->dim >= a->rows);
  assert(x->dim >= a->cols);

  if (n == 0)			/* Quick exit */
    return;

  if (atrans) {
    /* Left upper part, upper triangular */
    h2_trmv(_h2_lower, _h2_adj, (aunit ? _h2_unit : _h2_nonunit), &n, aa,
	    &a->ld, xv, &u_one);

    /* Right part */
    if (n < a->rows) {
      n1 = a->rows - n;
      h2_gemv(_h2_adj, &n1, &n, &f_one, aa + n, &a->ld, xv + n, &u_one,
	      &f_one, xv, &u_one);
    }

    /* Lower part */
    if (n < a->cols)
      for (i = n; i < a->cols; i++)
	xv[i] = 0.0;
  }
  else {
    /* Lower part */
    if (n < a->rows) {
      for (i = n; i < a->rows; i++)
	xv[i] = 0.0;

      n1 = a->rows - n;
      h2_gemv(_h2_ntrans, &n1, &n, &f_one, aa + n, &a->ld, xv, &u_one, &f_one,
	      xv + n, &u_one);
    }

    /* Top part, lower triangular */
    h2_trmv(_h2_lower, _h2_ntrans, (aunit ? _h2_unit : _h2_nonunit), &n, aa,
	    &a->ld, xv, &u_one);
  }
}

static void
uppereval_amatrix_avector(bool aunit, bool atrans, pcamatrix a, pavector x)
{
  field    *aa = a->a;
  field    *xv = x->v;
  longindex lda = a->ld;
  uint      n = UINT_MIN(a->rows, a->cols);
  uint      n1, i;

  assert(x->dim >= a->rows);
  assert(x->dim >= a->cols);

  if (n == 0)			/* Quick exit */
    return;

  if (atrans) {
    /* Lower part */
    if (n < a->cols) {
      for (i = n; i < a->cols; i++)
	xv[i] = 0.0;

      n1 = a->cols - n;
      h2_gemv(_h2_adj, &n, &n1, &f_one, aa + n * lda, &a->ld, xv, &u_one,
	      &f_one, xv + n, &u_one);
    }

    /* Top part, lower triangular */
    h2_trmv(_h2_upper, _h2_adj, (aunit ? _h2_unit : _h2_nonunit), &n, aa,
	    &a->ld, xv, &u_one);
  }
  else {
    /* Left upper part, upper triangular */
    h2_trmv(_h2_upper, _h2_ntrans, (aunit ? _h2_unit : _h2_nonunit), &n, aa,
	    &a->ld, xv, &u_one);

    /* Right part */
    if (n < a->cols) {
      n1 = a->cols - n;
      h2_gemv(_h2_ntrans, &n, &n1, &f_one, aa + n * lda, &a->ld, xv + n,
	      &u_one, &f_one, xv, &u_one);
    }

    /* Lower part */
    if (n < a->rows)
      for (i = n; i < a->rows; i++)
	xv[i] = 0.0;
  }
}
#else
static void
lowereval_amatrix_avector(bool aunit, bool atrans, pcamatrix a, pavector x)
{
  pcfield   aa = a->a;
  longindex lda = a->ld;
  pfield    xv = x->v;
  uint      n = UINT_MIN(a->rows, a->cols);
  field     newval;
  uint      i, j;

  assert(x->dim >= a->rows);
  assert(x->dim >= a->cols);

  if (atrans) {
    /* Left upper part, upper triangular */
    for (i = 0; i < n; i++) {
      newval = (aunit ? xv[i] : CONJ(aa[i + i * lda]) * xv[i]);
      for (j = i + 1; j < n; j++)
	newval += CONJ(aa[j + i * lda]) * xv[j];
      xv[i] = newval;
    }

    /* Right part */
    if (n < a->rows)
      for (i = 0; i < n; i++) {
	newval = xv[i];
	for (j = n; j < a->rows; j++)
	  newval += CONJ(aa[j + i * lda]) * xv[j];
	xv[i] = newval;
      }

    /* Lower part */
    if (n < a->cols)
      for (i = n; i < a->cols; i++)
	xv[i] = 0.0;
  }
  else {
    /* Lower part */
    if (n < a->rows)
      for (i = n; i < a->rows; i++) {
	newval = 0.0;
	for (j = 0; j < a->cols; j++)
	  newval += aa[i + j * lda] * xv[j];
	xv[i] = newval;
      }

    /* Top part, lower triangular */
    for (i = n; i-- > 0;) {
      newval = (aunit ? xv[i] : aa[i + i * lda] * xv[i]);
      for (j = 0; j < i; j++)
	newval += aa[i + j * lda] * xv[j];
      xv[i] = newval;
    }
  }
}

static void
uppereval_amatrix_avector(bool aunit, bool atrans, pcamatrix a, pavector x)
{
  pcfield   aa = a->a;
  longindex lda = a->ld;
  pfield    xv = x->v;
  uint      n = UINT_MIN(a->rows, a->cols);
  field     newval;
  uint      i, j;

  assert(x->dim >= a->rows);
  assert(x->dim >= a->cols);

  if (atrans) {
    /* Lower part */
    if (n < a->cols)
      for (i = n; i < a->cols; i++) {
	newval = 0.0;
	for (j = 0; j < n; j++)
	  newval += CONJ(aa[j + i * lda]) * xv[j];
	xv[i] = newval;
      }

    /* Top part, lower triangular */
    for (i = n; i-- > 0;) {
      newval = (aunit ? xv[i] : CONJ(aa[i + i * lda]) * xv[i]);
      for (j = 0; j < i; j++)
	newval += CONJ(aa[j + i * lda]) * xv[j];
      xv[i] = newval;
    }
  }
  else {
    /* Left upper part, upper triangular */
    for (i = 0; i < n; i++) {
      newval = (aunit ? xv[i] : xv[i] * aa[i + i * lda]);
      for (j = i + 1; j < n; j++)
	newval += aa[i + j * lda] * xv[j];
      xv[i] = newval;
    }

    /* Right part */
    if (n < a->cols)
      for (i = 0; i < n; i++) {
	newval = xv[i];
	for (j = n; j < a->cols; j++)
	  newval += aa[i + j * lda] * xv[j];
	xv[i] = newval;
      }

    /* Lower part */
    if (n < a->rows)
      for (i = n; i < a->rows; i++)
	xv[i] = 0.0;
  }
}
#endif

void
triangulareval_amatrix_avector(bool alower, bool aunit, bool atrans,
			       pcamatrix a, pavector x)
{
  if (alower)
    lowereval_amatrix_avector(aunit, atrans, a, x);
  else
    uppereval_amatrix_avector(aunit, atrans, a, x);
}

#ifdef USE_BLAS
static void
lowereval_amatrix(bool aunit, bool atrans, pcamatrix a, bool xtrans,
		  pamatrix x)
{
  field    *aa = a->a;
  field    *xa = x->a;
  longindex ldx = x->ld;
  uint      n = UINT_MIN(a->rows, a->cols);
  uint      n1, i, j;

  if (n == 0)			/* Quick exit */
    return;

  if (xtrans) {
    assert(x->cols >= a->rows);
    assert(x->cols >= a->cols);

    if (atrans) {
      /* Left upper part, upper triangular */
      h2_trmm(_h2_right, _h2_lower, _h2_ntrans,
	      (aunit ? _h2_unit : _h2_nonunit), &x->rows, &n, &f_one, aa,
	      &a->ld, xa, &x->ld);

      /* Right part */
      if (n < a->rows) {
	n1 = a->rows - n;
	h2_gemm(_h2_ntrans, _h2_ntrans, &x->rows, &n, &n1, &f_one,
		xa + n * ldx, &x->ld, aa + n, &a->ld, &f_one, xa, &x->ld);
      }

      /* Lower part */
      if (n < a->cols)
	for (j = 0; j < x->rows; j++)
	  for (i = n; i < a->cols; i++)
	    xa[j + i * ldx] = 0.0;
    }
    else {
      /* Lower part */
      if (n < a->rows) {
	for (j = 0; j < x->rows; j++)
	  for (i = n; i < a->rows; i++)
	    xa[j + i * ldx] = 0.0;

	n1 = a->rows - n;
	h2_gemm(_h2_ntrans, _h2_adj, &x->rows, &n1, &n, &f_one, xa, &x->ld,
		aa + n, &a->ld, &f_one, xa + n * ldx, &x->ld);
      }

      /* Top part, lower triangular */
      h2_trmm(_h2_right, _h2_lower, _h2_adj, (aunit ? _h2_unit : _h2_nonunit),
	      &x->rows, &n, &f_one, aa, &a->ld, xa, &x->ld);
    }
  }
  else {
    assert(x->rows >= a->rows);
    assert(x->rows >= a->cols);

    if (atrans) {
      /* Left upper part, upper triangular */
      h2_trmm(_h2_left, _h2_lower, _h2_adj, (aunit ? _h2_unit : _h2_nonunit),
	      &n, &x->cols, &f_one, aa, &a->ld, xa, &x->ld);

      /* Right part */
      if (n < a->rows) {
	n1 = a->rows - n;
	h2_gemm(_h2_adj, _h2_ntrans, &n, &x->cols, &n1, &f_one, aa + n,
		&a->ld, xa + n, &x->ld, &f_one, xa, &x->ld);
      }

      /* Lower part */
      if (n < a->cols)
	for (j = 0; j < x->cols; j++)
	  for (i = n; i < a->cols; i++)
	    xa[i + j * ldx] = 0.0;
    }
    else {
      /* Lower part */
      if (n < a->rows) {
	for (j = 0; j < x->cols; j++)
	  for (i = n; i < a->rows; i++)
	    xa[i + j * ldx] = 0.0;

	n1 = a->rows - n;
	h2_gemm(_h2_ntrans, _h2_ntrans, &n1, &x->cols, &n, &f_one, aa + n,
		&a->ld, xa, &x->ld, &f_one, xa + n, &x->ld);
      }

      /* Top part, lower triangular */
      h2_trmm(_h2_left, _h2_lower, _h2_ntrans,
	      (aunit ? _h2_unit : _h2_nonunit), &n, &x->cols, &f_one, aa,
	      &a->ld, xa, &x->ld);
    }
  }
}

static void
uppereval_amatrix(bool aunit, bool atrans, pcamatrix a, bool xtrans,
		  pamatrix x)
{
  field    *aa = a->a;
  longindex lda = a->ld;
  field    *xa = x->a;
  longindex ldx = x->ld;
  uint      n = UINT_MIN(a->rows, a->cols);
  uint      n1, i, j;

  if (n == 0)			/* Quick exit */
    return;

  if (xtrans) {
    assert(x->cols >= a->rows);
    assert(x->cols >= a->cols);

    if (atrans) {
      /* Lower part */
      if (n < a->cols) {
	for (j = 0; j < x->rows; j++)
	  for (i = n; i < a->cols; i++)
	    xa[j + i * ldx] = 0.0;

	n1 = a->cols - n;
	h2_gemm(_h2_ntrans, _h2_ntrans, &x->rows, &n1, &n, &f_one, xa, &x->ld,
		aa + n * lda, &a->ld, &f_one, xa + n * ldx, &x->ld);
      }

      /* Top part, lower triangular */
      h2_trmm(_h2_right, _h2_upper, _h2_ntrans,
	      (aunit ? _h2_unit : _h2_nonunit), &x->rows, &n, &f_one, aa,
	      &a->ld, xa, &x->ld);
    }
    else {
      /* Left upper part, upper triangular */
      h2_trmm(_h2_right, _h2_upper, _h2_adj, (aunit ? _h2_unit : _h2_nonunit),
	      &x->rows, &n, &f_one, aa, &a->ld, xa, &x->ld);

      /* Right part */
      if (n < a->cols) {
	n1 = a->cols - n;
	h2_gemm(_h2_ntrans, _h2_adj, &x->rows, &n, &n1, &f_one, xa + n * ldx,
		&x->ld, aa + n * lda, &a->ld, &f_one, xa, &x->ld);
      }

      /* Lower part */
      if (n < a->rows)
	for (j = 0; j < x->rows; j++)
	  for (i = n; i < a->rows; i++)
	    xa[j + i * ldx] = 0.0;
    }
  }
  else {
    assert(x->rows >= a->rows);
    assert(x->rows >= a->cols);

    if (atrans) {
      /* Lower part */
      if (n < a->cols) {
	for (j = 0; j < x->cols; j++)
	  for (i = n; i < a->cols; i++)
	    xa[i + j * ldx] = 0.0;

	n1 = a->cols - n;
	h2_gemm(_h2_adj, _h2_ntrans, &n1, &x->cols, &n, &f_one, aa + n * lda,
		&a->ld, xa, &x->ld, &f_one, xa + n, &x->ld);
      }

      /* Top part, lower triangular */
      h2_trmm(_h2_left, _h2_upper, _h2_adj, (aunit ? _h2_unit : _h2_nonunit),
	      &n, &x->cols, &f_one, aa, &a->ld, xa, &x->ld);
    }
    else {
      /* Left upper part, upper triangular */
      h2_trmm(_h2_left, _h2_upper, _h2_ntrans,
	      (aunit ? _h2_unit : _h2_nonunit), &n, &x->cols, &f_one, aa,
	      &a->ld, xa, &x->ld);

      /* Right part */
      if (n < a->cols) {
	n1 = a->cols - n;
	h2_gemm(_h2_ntrans, _h2_ntrans, &n, &x->cols, &n1, &f_one,
		aa + n * lda, &a->ld, xa + n, &x->ld, &f_one, xa, &x->ld);
      }

      /* Lower part */
      if (n < a->rows)
	for (j = 0; j < x->cols; j++)
	  for (i = n; i < a->rows; i++)
	    xa[i + j * ldx] = 0.0;
    }
  }
}
#else
static void
lowereval_amatrix(bool aunit, bool atrans, pcamatrix a, bool xtrans,
		  pamatrix x)
{
  pcfield   aa = a->a;
  longindex lda = a->ld;
  pfield    xa = x->a;
  longindex ldx = x->ld;
  uint      n = UINT_MIN(a->rows, a->cols);
  field     newval;
  uint      i, j, k;

  if (xtrans) {
    assert(x->cols >= a->rows);
    assert(x->cols >= a->cols);

    if (atrans) {
      for (k = 0; k < x->rows; k++) {
	/* Left upper part, upper triangular */
	for (i = 0; i < n; i++) {
	  newval =
	    (aunit ? xa[k + i * ldx] : aa[i + i * lda] * xa[k + i * ldx]);
	  for (j = i + 1; j < n; j++)
	    newval += aa[j + i * lda] * xa[k + j * ldx];
	  xa[k + i * ldx] = newval;
	}

	/* Right part */
	if (n < a->rows)
	  for (i = 0; i < n; i++) {
	    newval = xa[k + i * ldx];
	    for (j = n; j < a->rows; j++)
	      newval += aa[j + i * lda] * xa[k + j * ldx];
	    xa[k + i * ldx] = newval;
	  }

	/* Lower part */
	if (n < a->cols)
	  for (i = n; i < a->cols; i++)
	    xa[k + i * ldx] = 0.0;
      }
    }
    else {
      for (k = 0; k < x->rows; k++) {
	/* Lower part */
	if (n < a->rows)
	  for (i = n; i < a->rows; i++) {
	    newval = 0.0;
	    for (j = 0; j < a->cols; j++)
	      newval += CONJ(aa[i + j * lda]) * xa[k + j * ldx];
	    xa[k + i * ldx] = newval;
	  }

	/* Top part, lower triangular */
	for (i = n; i-- > 0;) {
	  newval =
	    (aunit ? xa[k + i * ldx] : CONJ(aa[i + i * lda]) *
	     xa[k + i * ldx]);
	  for (j = 0; j < i; j++)
	    newval += CONJ(aa[i + j * lda]) * xa[k + j * ldx];
	  xa[k + i * ldx] = newval;
	}
      }
    }
  }
  else {
    assert(x->rows >= a->rows);
    assert(x->rows >= a->cols);

    if (atrans) {
      for (k = 0; k < x->cols; k++) {
	/* Left upper part, upper triangular */
	for (i = 0; i < n; i++) {
	  newval =
	    (aunit ? xa[i + k * ldx] : CONJ(aa[i + i * lda]) *
	     xa[i + k * ldx]);
	  for (j = i + 1; j < n; j++)
	    newval += CONJ(aa[j + i * lda]) * xa[j + k * ldx];
	  xa[i + k * ldx] = newval;
	}

	/* Right part */
	if (n < a->rows)
	  for (i = 0; i < n; i++) {
	    newval = xa[i + k * ldx];
	    for (j = n; j < a->rows; j++)
	      newval += CONJ(aa[j + i * lda]) * xa[j + k * ldx];
	    xa[i + k * ldx] = newval;
	  }

	/* Lower part */
	if (n < a->cols)
	  for (i = n; i < a->cols; i++)
	    xa[i + k * ldx] = 0.0;
      }
    }
    else {
      for (k = 0; k < x->cols; k++) {
	/* Lower part */
	if (n < a->rows)
	  for (i = n; i < a->rows; i++) {
	    newval = 0.0;
	    for (j = 0; j < a->cols; j++)
	      newval += aa[i + j * lda] * xa[j + k * ldx];
	    xa[i + k * ldx] = newval;
	  }

	/* Top part, lower triangular */
	for (i = n; i-- > 0;) {
	  newval =
	    (aunit ? xa[i + k * ldx] : aa[i + i * lda] * xa[i + k * ldx]);
	  for (j = 0; j < i; j++)
	    newval += aa[i + j * lda] * xa[j + k * ldx];
	  xa[i + k * ldx] = newval;
	}
      }
    }
  }
}

static void
uppereval_amatrix(bool aunit, bool atrans, pcamatrix a, bool xtrans,
		  pamatrix x)
{
  pcfield   aa = a->a;
  longindex lda = a->ld;
  pfield    xa = x->a;
  longindex ldx = x->ld;
  uint      n = UINT_MIN(a->rows, a->cols);
  field     newval;
  uint      i, j, k;

  if (xtrans) {
    assert(x->cols >= a->rows);
    assert(x->cols >= a->cols);

    if (atrans) {
      for (k = 0; k < x->rows; k++) {
	/* Lower part */
	if (n < a->cols)
	  for (i = n; i < a->cols; i++) {
	    newval = 0.0;
	    for (j = 0; j < a->rows; j++)
	      newval += aa[j + i * lda] * xa[k + j * ldx];
	    xa[k + i * ldx] = newval;
	  }

	/* Top part, lower triangular */
	for (i = n; i-- > 0;) {
	  newval =
	    (aunit ? xa[k + i * ldx] : aa[i + i * lda] * xa[k + i * ldx]);
	  for (j = 0; j < i; j++)
	    newval += aa[j + i * lda] * xa[k + j * ldx];
	  xa[k + i * ldx] = newval;
	}
      }
    }
    else {
      for (k = 0; k < x->rows; k++) {
	/* Left upper part, upper triangular */
	for (i = 0; i < n; i++) {
	  newval =
	    (aunit ? xa[k + i * ldx] : CONJ(aa[i + i * lda]) *
	     xa[k + i * ldx]);
	  for (j = i + 1; j < n; j++)
	    newval += CONJ(aa[i + j * lda]) * xa[k + j * ldx];
	  xa[k + i * ldx] = newval;
	}

	/* Right part */
	if (n < a->cols)
	  for (i = 0; i < n; i++) {
	    newval = xa[k + i * ldx];
	    for (j = n; j < a->cols; j++)
	      newval += CONJ(aa[i + j * lda]) * xa[k + j * ldx];
	    xa[k + i * ldx] = newval;
	  }

	/* Lower part */
	if (n < a->rows)
	  for (i = n; i < a->rows; i++)
	    xa[k + i * ldx] = 0.0;
      }
    }
  }
  else {
    assert(x->rows >= a->rows);
    assert(x->rows >= a->cols);

    if (atrans) {
      for (k = 0; k < x->cols; k++) {
	/* Lower part */
	if (n < a->cols)
	  for (i = n; i < a->cols; i++) {
	    newval = 0.0;
	    for (j = 0; j < a->rows; j++)
	      newval += CONJ(aa[j + i * lda]) * xa[j + k * ldx];
	    xa[i + k * ldx] = newval;
	  }

	/* Top part, lower triangular */
	for (i = n; i-- > 0;) {
	  newval =
	    (aunit ? xa[i + k * ldx] : CONJ(aa[i + i * lda]) *
	     xa[i + k * ldx]);
	  for (j = 0; j < i; j++)
	    newval += CONJ(aa[j + i * lda]) * xa[j + k * ldx];
	  xa[i + k * ldx] = newval;
	}
      }
    }
    else {
      for (k = 0; k < x->cols; k++) {
	/* Left upper part, upper triangular */
	for (i = 0; i < n; i++) {
	  newval =
	    (aunit ? xa[i + k * ldx] : aa[i + i * lda] * xa[i + k * ldx]);
	  for (j = i + 1; j < n; j++)
	    newval += aa[i + j * lda] * xa[j + k * ldx];
	  xa[i + k * ldx] = newval;
	}

	/* Right part */
	if (n < a->cols)
	  for (i = 0; i < n; i++) {
	    newval = xa[i + k * ldx];
	    for (j = n; j < a->cols; j++)
	      newval += aa[i + j * lda] * xa[j + k * ldx];
	    xa[i + k * ldx] = newval;
	  }

	/* Lower part */
	if (n < a->rows)
	  for (i = n; i < a->rows; i++)
	    xa[i + k * ldx] = 0.0;
      }
    }
  }
}
#endif

void
triangulareval_amatrix(bool alower, bool aunit, bool atrans, pcamatrix a,
		       bool xtrans, pamatrix x)
{
  if (alower)
    lowereval_amatrix(aunit, atrans, a, xtrans, x);
  else
    uppereval_amatrix(aunit, atrans, a, xtrans, x);
}

void
triangularaddmul_amatrix(field alpha, bool alower, bool atrans,
			 pcamatrix a, bool blower, bool btrans, pcamatrix b,
			 pamatrix c)
{
  pcfield   aa = a->a;
  longindex lda = a->ld;
  pcfield   ba = b->a;
  longindex ldb = b->ld;
  pfield    ca = c->a;
  longindex ldc = c->ld;
  uint      aoff, adim, ainc, boff, bdim, binc;
  uint      j;
  uint      i, k;

  if (atrans) {
    assert(c->rows == a->cols);

    ainc = lda;
    lda = 1;

    if (btrans) {		/* C += alpha A^* B^* */
      assert(a->rows == b->cols);
      assert(c->cols == b->rows);

      binc = 1;

      for (j = 0; j < a->rows; j++) {
	if (alower) {		/* A^* upper triangular */
	  aoff = 0;
	  adim = UINT_MIN(j + 1, a->cols);
	}
	else {			/* A^* lower triangular */
	  aoff = j;
	  adim = a->cols - UINT_MIN(j, a->cols);
	}

	if (blower) {		/* B^* upper triangular */
	  boff = j;
	  bdim = b->rows - UINT_MIN(j, b->rows);
	}
	else {			/* B^* lower triangular */
	  boff = 0;
	  bdim = UINT_MIN(j + 1, b->rows);
	}

	//TODO Need a 'ger' equivalent for A \gets A + alpha * CONJ(x) * y ** H
	//#ifdef USE_BLAS
	//        h2_geru(&adim, &bdim, &alpha, aa + aoff * ainc + j * lda, &ainc,
	//            ba + boff * binc + j * ldb, &binc, ca + aoff + boff * ldc, &ldc);
	//#else

	for (k = 0; k < bdim; k++) {
	  for (i = 0; i < adim; i++) {
	    ca[(aoff + i) + (boff + k) * ldc] += alpha
	      * CONJ(aa[(aoff + i) * ainc + j * lda]) *
	      CONJ(ba[(boff + k) * binc + j * ldb]);
	  }
	}
	//#endif
      }
    }
    else {			/* C += alpha A^* B */
      assert(a->rows == b->rows);
      assert(c->cols == b->cols);

      binc = ldb;
      ldb = 1;

      for (j = 0; j < a->rows; j++) {
	if (alower) {		/* A^* upper triangular */
	  aoff = 0;
	  adim = UINT_MIN(j + 1, a->cols);
	}
	else {			/* A^* lower triangular */
	  aoff = j;
	  adim = a->cols - UINT_MIN(j, a->cols);
	}

	if (blower) {		/* B lower triangular */
	  boff = 0;
	  bdim = UINT_MIN(j + 1, b->cols);
	}
	else {			/* B upper triangular */
	  boff = j;
	  bdim = b->cols - UINT_MIN(j, b->cols);
	}

	//TODO Need a 'ger' equivalent for A \gets A + alpha * CONJ(x) * y ** T
	//#ifdef USE_BLAS
	//        h2_gerc(&adim, &bdim, &alpha, aa + aoff * ainc + j * lda, &ainc,
	//            ba + boff * binc + j * ldb, &binc, ca + aoff + boff * ldc, &ldc);
	//#else
	for (k = 0; k < bdim; k++) {
	  for (i = 0; i < adim; i++) {
	    ca[(aoff + i) + (boff + k) * ldc] += alpha
	      * CONJ(aa[(aoff + i) * ainc + j * lda])
	      * ba[(boff + k) * binc + j * ldb];
	  }
	}
	//#endif
      }
    }
  }
  else {
    assert(c->rows == a->rows);

    ainc = 1;

    if (btrans) {		/* C += alpha A B^* */
      assert(a->cols == b->cols);
      assert(c->cols == b->rows);

      binc = 1;

      for (j = 0; j < a->cols; j++) {
	if (alower) {		/* A lower triangular */
	  aoff = j;
	  adim = a->rows - UINT_MIN(j, a->rows);
	}
	else {			/* A upper triangular */
	  aoff = 0;
	  adim = UINT_MIN(j + 1, a->rows);
	}

	if (blower) {		/* B^* upper triangular */
	  boff = j;
	  bdim = b->rows - UINT_MIN(j, b->rows);
	}
	else {			/* B^* lower triangular */
	  boff = 0;
	  bdim = UINT_MIN(j + 1, b->rows);
	}

#ifdef USE_BLAS
	h2_gerc(&adim, &bdim, &alpha, aa + aoff * ainc + j * lda, &ainc,
		ba + boff * binc + j * ldb, &binc, ca + aoff + boff * ldc,
		&c->ld);
#else
	for (k = 0; k < bdim; k++) {
	  for (i = 0; i < adim; i++) {
	    ca[(aoff + i) + (boff + k) * ldc] += alpha
	      * aa[(aoff + i) * ainc + j * lda]
	      * CONJ(ba[(boff + k) * binc + j * ldb]);
	  }
	}
#endif
      }
    }
    else {			/* C += alpha A B */
      assert(a->cols == b->rows);
      assert(c->cols == b->cols);

      binc = ldb;
      ldb = 1;

      for (j = 0; j < a->cols; j++) {
	if (alower) {		/* A lower triangular */
	  aoff = j;
	  adim = a->rows - UINT_MIN(j, a->rows);
	}
	else {			/* A upper triangular */
	  aoff = 0;
	  adim = UINT_MIN(j + 1, a->rows);
	}

	if (blower) {		/* B lower triangular */
	  boff = 0;
	  bdim = UINT_MIN(j + 1, b->cols);
	}
	else {			/* B upper triangular */
	  boff = j;
	  bdim = b->cols - UINT_MIN(j, b->cols);
	}

#ifdef USE_BLAS
	h2_geru(&adim, &bdim, &alpha, aa + aoff * ainc + j * lda, &ainc,
		ba + boff * binc + j * ldb, &binc, ca + aoff + boff * ldc,
		&c->ld);
#else

	for (k = 0; k < bdim; k++) {
	  for (i = 0; i < adim; i++) {
	    ca[(aoff + i) + (boff + k) * ldc] += alpha
	      * aa[(aoff + i) * ainc + j * lda]
	      * ba[(boff + k) * binc + j * ldb];
	  }
	}
#endif
      }
    }
  }
}

void
copy_lower_amatrix(pcamatrix a, bool aunit, pamatrix b)
{
  pfield    aa = a->a;
  longindex lda = a->ld;
  pfield    ba = b->a;
  longindex ldb = b->ld;
  uint      rows;
  uint      cols;
  uint      i, j;

  rows = UINT_MIN(a->rows, b->rows);
  cols = UINT_MIN(a->cols, b->cols);

  for (j = 0; j < cols; j++) {
    if (aunit) {
      for (i = 0; i < rows && i < j; i++) {
	ba[i + j * ldb] = 0.0;
      }

      if (i < rows && i == j) {
	ba[i + j * ldb] = 1.0;
	i++;
      }
    }
    else {
      for (i = 0; i < rows && i < j; i++) {
	ba[i + j * ldb] = 0.0;
      }
    }

    for (; i < rows; i++) {
      ba[i + j * ldb] = aa[i + j * lda];
    }
  }
  for (; j < b->cols; j++) {
    for (i = 0; i < rows; i++) {
      ba[i + j * ldb] = 0.0;
    }
  }
}

void
copy_upper_amatrix(pcamatrix a, bool aunit, pamatrix b)
{
  pfield    aa = a->a;
  longindex lda = a->ld;
  pfield    ba = b->a;
  longindex ldb = b->ld;
  uint      rows;
  uint      cols;
  uint      i, j;

  rows = UINT_MIN(a->rows, b->rows);
  cols = UINT_MIN(a->cols, b->cols);

  for (j = 0; j < cols; j++) {
    if (aunit) {
      for (i = 0; i < rows && i < j; i++) {
	ba[i + j * ldb] = aa[i + j * lda];
      }

      if (i < rows && i == j) {
	ba[i + j * ldb] = 1.0;
	i++;
      }
    }
    else {
      for (i = 0; i < rows && i <= j; i++) {
	ba[i + j * ldb] = aa[i + j * lda];
      }
    }

    for (; i < b->rows; i++) {
      ba[i + j * ldb] = 0.0;
    }
  }
}

/* ------------------------------------------------------------
 * LR decomposition
 * ------------------------------------------------------------ */

#ifdef USE_BLAS
uint
lrdecomp_amatrix(pamatrix a)
{
  field    *aa = a->a;
  longindex lda = a->ld;
  uint      n = a->rows;
  field     alpha;
  uint      i, n1;

  assert(n == a->cols);

  for (i = 0; i < n - 1; i++) {
    if (aa[i + i * lda] == 0.0)
      return i + 1;

    alpha = 1.0 / aa[i + i * lda];

    n1 = n - i - 1;
    h2_scal(&n1, &alpha, aa + (i + 1) + i * lda, &u_one);
    h2_geru(&n1, &n1, &f_minusone, aa + (i + 1) + i * lda, &u_one,
	    aa + i + (i + 1) * lda, &a->ld, aa + (i + 1) + (i + 1) * lda,
	    &a->ld);
  }

  if (aa[i + i * lda] == 0.0)
    return i + 1;

  return 0;
}
#else
uint
lrdecomp_amatrix(pamatrix a)
{
  pfield    aa = a->a;
  longindex lda = a->ld;
  uint      n = a->rows;
  field     alpha;
  uint      i, j, k;

  assert(n == a->cols);

  for (i = 0; i < n - 1; i++) {
    if (aa[i + i * lda] == 0.0)
      return i + 1;

    alpha = 1.0 / aa[i + i * lda];

    for (j = i + 1; j < n; j++) {
      aa[j + i * lda] *= alpha;
      for (k = i + 1; k < n; k++)
	aa[j + k * lda] -= aa[j + i * lda] * aa[i + k * lda];
    }
  }

  if (aa[i + i * lda] == 0.0)
    return i + 1;

  return 0;
}
#endif

uint
lrdecomp_blocks_amatrix(pamatrix a, uint blocksize)
{
  amatrix   tmp1, tmp2, tmp3;
  pamatrix  d, l, r, s;
  uint      blocks;
  uint      n, ioff, isize, joff, jsize, koff, ksize;
  uint      i, j, k;

  assert(a->rows == a->cols);

  n = a->rows;

  /* Use standard function if the matrix is small */
  if (n <= blocksize)
    return lrdecomp_amatrix(a);

  blocks = n / blocksize;

  for (k = 0; k < blocks; k++) {
    /* Initialize k-th diagonal block */
    koff = n * k / blocks;
    ksize = n * (k + 1) / blocks - koff;
    d = init_sub_amatrix(&tmp1, a, ksize, koff, ksize, koff);

    /* Compute LR factorization of diagonal block */
    i = lrdecomp_amatrix(d);
    if (i > 0) {
      uninit_amatrix(d);
      return i + koff;
    }

    /* Solve A_{kj} = L_{kk} R_{kj} */
    for (j = k + 1; j < blocks; j++) {
      joff = n * j / blocks;
      jsize = n * (j + 1) / blocks - joff;
      r = init_sub_amatrix(&tmp2, a, ksize, koff, jsize, joff);

      triangularsolve_amatrix(true, true, false, d, false, r);

      uninit_amatrix(r);
    }

    /* Solve A_{ik} = L_{ik} R_{kk}, i.e., A_{ik}^* = R_{kk}^* L_{ik}^* */
    for (i = k + 1; i < blocks; i++) {
      ioff = n * i / blocks;
      isize = n * (i + 1) / blocks - ioff;
      l = init_sub_amatrix(&tmp3, a, isize, ioff, ksize, koff);

      triangularsolve_amatrix(false, false, true, d, true, l);

      uninit_amatrix(l);
    }

    uninit_amatrix(d);

    /* Update A_{ij} = A_{ij} - L_{ik} R_{kj} */
    for (j = k + 1; j < blocks; j++) {
      joff = n * j / blocks;
      jsize = n * (j + 1) / blocks - joff;
      r = init_sub_amatrix(&tmp2, a, ksize, koff, jsize, joff);

      for (i = k + 1; i < blocks; i++) {
	ioff = n * i / blocks;
	isize = n * (i + 1) / blocks - ioff;
	l = init_sub_amatrix(&tmp3, a, isize, ioff, ksize, koff);

	s = init_sub_amatrix(&tmp1, a, isize, ioff, jsize, joff);

	addmul_amatrix(-1.0, false, l, false, r, s);

	uninit_amatrix(s);
	uninit_amatrix(l);
      }

      uninit_amatrix(r);
    }
  }

  return 0;
}

#if defined(USE_OPENMP) && _OPENMP >= 200805
#define USE_OPENMP_TASKS
#endif

void
lrdecomp_tasks_amatrix(pamatrix a, uint blocksize)
{
#ifndef USE_OPENMP_TASKS
  (void) lrdecomp_blocks_amatrix(a, blocksize);
#else
  pamatrix *sa;
  uint      blocks;
  uint      n, roff, rsize, coff, csize;
  uint      i, j, k;

  assert(a->rows == a->cols);

  n = a->rows;

  /* Use standard function if the matrix is small */
  if (n <= blocksize)
    (void) lrdecomp_amatrix(a);

  blocks = n / blocksize;

  /* Split A into submatrices A_{ij} */
  sa = (pamatrix *) allocmem(sizeof(pamatrix) * blocks * blocks);
  for (j = 0; j < blocks; j++) {
    coff = n * j / blocks;
    csize = n * (j + 1) / blocks - coff;

    for (i = 0; i < blocks; i++) {
      roff = n * i / blocks;
      rsize = n * (i + 1) / blocks - roff;

      sa[i + j * blocks] = new_sub_amatrix(a, rsize, roff, csize, coff);
    }
  }

#pragma omp parallel private(i,j,k)
  {
    for (k = 0; k < blocks; k++) {
      /* Compute decomposition L_{kk} R_{kk} = A_{kk} */
#pragma omp single
      lrdecomp_amatrix(sa[k + k * blocks]);

      /* Solve A_{kj} = L_{kk} R_{kj} */
#pragma omp single
      for (j = k + 1; j < blocks; j++)
#pragma omp task
	triangularsolve_amatrix(true, true,
				false, sa[k + k * blocks],
				false, sa[k + j * blocks]);

      /* Solve A_{ik} = L_{ik} R_{kk}, i.e., A_{ik}^* = R_{kk}^* L_{ik}^* */
#pragma omp single
      for (i = k + 1; i < blocks; i++)
#pragma omp task
	triangularsolve_amatrix(false, false,
				true, sa[k + k * blocks],
				true, sa[i + k * blocks]);

#pragma omp taskwait

      /* Update A_{ij} = A_{ij} - L_{ik} R_{kj} */
#pragma omp single
      for (j = k + 1; j < blocks; j++)
	for (i = k + 1; i < blocks; i++)
#pragma omp task
	  addmul_amatrix(-1.0, false, sa[i + k * blocks], false,
			 sa[k + j * blocks], sa[i + j * blocks]);

#pragma omp taskwait
    }
  }

  /* Clean up */
  for (j = 0; j < blocks; j++)
    for (i = 0; i < blocks; i++)
      del_amatrix(sa[i + j * blocks]);
  freemem(sa);
#endif
}

void
lrsolve_n_amatrix_avector(pcamatrix a, pavector x)
{
  triangularsolve_amatrix_avector(true, true, false, a, x);
  triangularsolve_amatrix_avector(false, false, false, a, x);
}

void
lrsolve_t_amatrix_avector(pcamatrix a, pavector x)
{
  triangularsolve_amatrix_avector(false, false, true, a, x);
  triangularsolve_amatrix_avector(true, true, true, a, x);
}

void
lrsolve_amatrix_avector(bool atrans, pcamatrix a, pavector x)
{
  if (atrans) {
    lrsolve_t_amatrix_avector(a, x);
  }
  else {
    lrsolve_n_amatrix_avector(a, x);
  }
}

void
lrsolve_amatrix(pcamatrix a, pamatrix x)
{
  triangularsolve_amatrix(true, true, false, a, false, x);
  triangularsolve_amatrix(false, false, false, a, false, x);
}

void
lreval_n_amatrix_avector(pcamatrix a, pavector x)
{
  triangulareval_amatrix_avector(false, false, false, a, x);
  triangulareval_amatrix_avector(true, true, false, a, x);
}

void
lreval_t_amatrix_avector(pcamatrix a, pavector x)
{
  triangulareval_amatrix_avector(true, true, true, a, x);
  triangulareval_amatrix_avector(false, false, true, a, x);
}

void
lreval_amatrix_avector(bool atrans, pcamatrix a, pavector x)
{
  if (atrans) {
    lreval_t_amatrix_avector(a, x);
  }
  else {
    lreval_n_amatrix_avector(a, x);
  }
}

/* ------------------------------------------------------------
 * Cholesky decomposition
 * ------------------------------------------------------------ */

#ifdef USE_BLAS
uint
choldecomp_amatrix(pamatrix a)
{
  field    *aa = a->a;
  uint      n = a->rows;

  int       info;

  assert(n == a->cols);

  h2_potrf(_h2_lower, &n, aa, &a->ld, &info);

  return info;
}
#else
uint
choldecomp_amatrix(pamatrix a)
{
  pfield    aa = a->a;
  longindex lda = a->ld;
  uint      n = a->rows;
  real      diag, alpha;
  uint      i, j, k;

  assert(n == a->cols);

  for (i = 0; i < n - 1; i++) {
    diag = REAL(aa[i + i * lda]);

    if (ABS(aa[i + i * lda] - diag) > 1e-12 || diag <= 0.0)
      return i + 1;

    aa[i + i * lda] = REAL_SQRT(diag);
    alpha = 1.0 / aa[i + i * lda];
    for (j = i + 1; j < n; j++)
      aa[j + i * lda] *= alpha;

    for (j = i + 1; j < n; j++)
      for (k = i + 1; k <= j; k++)
	aa[j + k * lda] -= aa[j + i * lda] * CONJ(aa[k + i * lda]);
  }

  diag = REAL(aa[i + i * lda]);
  if (ABS(aa[i + i * lda] - diag) > 1e-12 || diag <= 0.0)
    return i + 1;

  aa[i + i * lda] = REAL_SQRT(diag);

  return 0;
}
#endif

uint
choldecomp_blocks_amatrix(pamatrix a, uint blocksize)
{
  amatrix   tmp1, tmp2, tmp3;
  pamatrix  d, l, r, s;
  uint      blocks;
  uint      n, ioff, isize, joff, jsize, koff, ksize;
  uint      i, j, k;

  assert(a->rows == a->cols);

  n = a->rows;

  /* Use standard function if the matrix is small */
  if (n <= blocksize)
    return choldecomp_amatrix(a);

  blocks = n / blocksize;

  for (k = 0; k < blocks; k++) {
    /* Initialize k-th diagonal block */
    koff = n * k / blocks;
    ksize = n * (k + 1) / blocks - koff;
    d = init_sub_amatrix(&tmp1, a, ksize, koff, ksize, koff);

    /* Compute Cholesky factorization of diagonal block */
    i = choldecomp_amatrix(d);
    if (i > 0) {
      uninit_amatrix(d);
      return i + koff;
    }

    /* Solve A_{ik} = L_{ik} L_{kk}^*, i.e., A_{ik}^* = L_{kk} L_{ik}^* */
    for (i = k + 1; i < blocks; i++) {
      ioff = n * i / blocks;
      isize = n * (i + 1) / blocks - ioff;
      l = init_sub_amatrix(&tmp3, a, isize, ioff, ksize, koff);

      triangularsolve_amatrix(true, false, false, d, true, l);

      uninit_amatrix(l);
    }

    uninit_amatrix(d);

    /* Update A_{ij} = A_{ij} - L_{ik} L_{jk}^* */
    for (j = k + 1; j < blocks; j++) {
      joff = n * j / blocks;
      jsize = n * (j + 1) / blocks - joff;
      r = init_sub_amatrix(&tmp2, a, jsize, joff, ksize, koff);

      for (i = j; i < blocks; i++) {
	ioff = n * i / blocks;
	isize = n * (i + 1) / blocks - ioff;
	l = init_sub_amatrix(&tmp3, a, isize, ioff, ksize, koff);

	s = init_sub_amatrix(&tmp1, a, isize, ioff, jsize, joff);

	addmul_amatrix(-1.0, false, l, true, r, s);

	uninit_amatrix(s);
	uninit_amatrix(l);
      }

      uninit_amatrix(r);
    }
  }

  return 0;
}

void
choldecomp_tasks_amatrix(pamatrix a, uint blocksize)
{
#ifndef USE_OPENMP_TASKS
  (void) choldecomp_blocks_amatrix(a, blocksize);
#else
  pamatrix *sa;
  uint      blocks;
  uint      n, roff, rsize, coff, csize;
  uint      i, j, k;

  assert(a->rows == a->cols);

  n = a->rows;

  /* Use standard function if the matrix is small */
  if (n <= blocksize)
    (void) choldecomp_amatrix(a);

  blocks = n / blocksize;
  if (blocks > 512)
    blocks = 512;

  /* Split A into submatrices A_{ij} */
  sa = (pamatrix *) allocmem(sizeof(pamatrix) * blocks * blocks);
  for (j = 0; j < blocks; j++) {
    coff = n * j / blocks;
    csize = n * (j + 1) / blocks - coff;

    for (i = 0; i < j; i++)
      sa[i + j * blocks] = 0;

    for (; i < blocks; i++) {
      roff = n * i / blocks;
      rsize = n * (i + 1) / blocks - roff;

      sa[i + j * blocks] = new_sub_amatrix(a, rsize, roff, csize, coff);
    }
  }

#pragma omp parallel private(i,j,k)
  {
    for (k = 0; k < blocks; k++) {
      /* Compute decomposition L_{kk} L_{kk}^* = A_{kk} */
#pragma omp single
      choldecomp_amatrix(sa[k + k * blocks]);

      /* Solve A_{ik} = L_{ik} L_{kk}^*, i.e., A_{ik}^* = L_{kk} L_{ik}^* */
#pragma omp single
      for (i = k + 1; i < blocks; i++)
#pragma omp task
	triangularsolve_amatrix(true, false,
				false, sa[k + k * blocks],
				true, sa[i + k * blocks]);

#pragma omp taskwait

      /* Update A_{ij} = A_{ij} - L_{ik} L_{jk} */
#pragma omp single
      for (j = 0; j < (blocks - k) / 2; j++)
#pragma omp task private(i)
      {
	uint      j1 = j + k + 1;
	uint      j2 = blocks - 1 - j;

	for (i = j1; i < blocks; i++)
	  addmul_amatrix(-1.0,
			 false, sa[i + k * blocks],
			 true, sa[j1 + k * blocks], sa[i + j1 * blocks]);

	if (j1 < j2)
	  for (i = j2; i < blocks; i++)
	    addmul_amatrix(-1.0,
			   false, sa[i + k * blocks],
			   true, sa[j2 + k * blocks], sa[i + j2 * blocks]);
      }

#pragma omp taskwait
    }
  }

  /* Clean up */
  for (j = 0; j < blocks; j++)
    for (i = j; i < blocks; i++)
      del_amatrix(sa[i + j * blocks]);
  freemem(sa);
#endif
}

void
cholsolve_amatrix_avector(pcamatrix a, pavector x)
{
  triangularsolve_amatrix_avector(true, false, false, a, x);
  triangularsolve_amatrix_avector(true, false, true, a, x);
}

void
cholsolve_amatrix(pcamatrix a, pamatrix x)
{
  assert(a->cols == x->rows);

  triangularsolve_amatrix(true, false, false, a, false, x);
  triangularsolve_amatrix(true, false, true, a, false, x);
}

void
choleval_amatrix_avector(pcamatrix a, pavector x)
{
  triangulareval_amatrix_avector(true, false, true, a, x);
  triangulareval_amatrix_avector(true, false, false, a, x);
}

/* ------------------------------------------------------------
 * LDL^T decomposition
 * ------------------------------------------------------------ */

#ifdef USE_BLAS
uint
ldltdecomp_amatrix(pamatrix a)
{
  field    *aa = a->a;
  longindex lda = a->ld;
  uint      n = a->rows;
  real      diag, alpha;
  uint      i, n1;

  assert(n == a->cols);

  for (i = 0; i < n - 1; i++) {
    diag = REAL(aa[i + i * lda]);

    if (ABS(aa[i + i * lda] - diag) > 1e-12 || diag == 0.0)
      return i + 1;

    alpha = 1.0 / diag;
    n1 = n - i - 1;
    h2_rscal(&n1, &alpha, aa + (i + 1) + i * lda, &u_one);

    alpha = -diag;
    h2_syr(_h2_lower, &n1, &alpha, aa + (i + 1) + i * lda, &u_one,
	   aa + (i + 1) + (i + 1) * lda, &a->ld);
  }

  diag = REAL(aa[i + i * lda]);
  if (ABS(aa[i + i * lda] - diag) > 1e-12 || diag == 0.0)
    return i + 1;

  return 0;
}
#else
uint
ldltdecomp_amatrix(pamatrix a)
{
  pfield    aa = a->a;
  longindex lda = a->ld;
  uint      n = a->rows;
  real      diag, alpha;
  uint      i, j, k;

  assert(n == a->cols);

  for (i = 0; i < n - 1; i++) {
    diag = REAL(aa[i + i * lda]);

    if (ABS(aa[i + i * lda] - diag) > 1e-12 || diag == 0.0)
      return i + 1;

    alpha = 1.0 / diag;
    for (j = i + 1; j < n; j++)
      aa[j + i * lda] *= alpha;

    for (j = i + 1; j < n; j++)
      for (k = i + 1; k <= j; k++)
	aa[j + k * lda] -= diag * aa[j + i * lda] * CONJ(aa[k + i * lda]);
  }

  diag = REAL(aa[i + i * lda]);
  if (ABS(aa[i + i * lda] - diag) > 1e-12 || diag == 0.0)
    return i + 1;

  return 0;
}
#endif

void
ldltsolve_amatrix_avector(pcamatrix a, pavector x)
{
  triangularsolve_amatrix_avector(true, true, false, a, x);
  diagsolve_amatrix_avector(false, a, x);
  triangularsolve_amatrix_avector(true, true, true, a, x);
}

void
ldltsolve_amatrix(pcamatrix a, pamatrix x)
{
  assert(a->cols == x->rows);

  triangularsolve_amatrix(true, true, false, a, false, x);
  diagsolve_amatrix(false, a, false, x);
  triangularsolve_amatrix(true, true, true, a, false, x);
}

/* ------------------------------------------------------------
 * Orthogonal decompositions
 * ------------------------------------------------------------ */

#ifdef USE_BLAS
void
qrdecomp_amatrix(pamatrix a, pavector tau)
{
  uint      rows = a->rows;
  uint      cols = a->cols;
  uint      refl = UINT_MIN(rows, cols);
  field    *work;
  int       lwork, info;

  assert(a->ld >= rows);
  /* Quick exit if no reflections used */
  if (refl == 0)
    return;

  lwork = 4 * cols;
  work = allocfield(lwork);

  if (tau->dim < refl)
    resize_avector(tau, refl);

  h2_geqrf(&rows, &cols, a->a, &a->ld, tau->v, work, &lwork, &info);
  assert(info == 0);

  freemem(work);
}
#else
void
qrdecomp_amatrix(pamatrix a, pavector tau)
{
  pfield    aa = a->a;
  longindex lda = a->ld;
  pfield    tauv;
  uint      rows = a->rows;
  uint      cols = a->cols;
  uint      refl = UINT_MIN(rows, cols);
  field     alpha, beta, gamma, diag;
  real      norm2, norm, adiag;
  uint      i, j, k;

  /* Provide storage for scaling factors */
  if (tau->dim < refl)
    resize_avector(tau, refl);
  tauv = tau->v;

  for (k = 0; k < refl; k++) {
    /* Compute norm of k-th column */
    norm2 = 0.0;
    for (i = k; i < rows; i++)
      norm2 += ABSSQR(aa[i + k * lda]);
    norm = REAL_SQRT(norm2);

    if (norm2 <= H2_ALMOST_ZERO)
      tauv[k] = 0.0;
    else {
      /* Determine reflection vector v */
      diag = aa[k + k * lda];
      adiag = ABS(diag) + norm;
      alpha = -SIGN1(diag) * norm;

      /* Compute norm of v */
      beta = 1.0 / (norm2 + ABS(diag) * norm);

      /* Rescale to ensure v_1 = 1 */
      beta *= adiag * adiag;
      gamma = CONJ(SIGN1(diag)) / adiag;
      for (i = k + 1; i < rows; i++)
	aa[i + k * lda] *= gamma;
      tauv[k] = beta;

      /* Compute k-th column */
      aa[k + k * lda] = alpha;

      /* Update remaining columns */
      for (j = k + 1; j < cols; j++) {
	gamma = aa[k + j * lda];
	for (i = k + 1; i < rows; i++)
	  gamma += CONJ(aa[i + k * lda]) * aa[i + j * lda];

	gamma *= beta;

	aa[k + j * lda] -= gamma;
	for (i = k + 1; i < rows; i++)
	  aa[i + j * lda] -= gamma * aa[i + k * lda];
      }
    }
  }
}
#endif

#ifdef USE_BLAS
uint
qrdecomp_pivot_amatrix(pamatrix a, pavector tau, uint * colpiv)
{
  pfield    aa = a->a;
  longindex lda = a->ld;
  pfield    tauv;
  pfield    work;
  uint      rows = a->rows;
  uint      cols = a->cols;
  uint      refl = UINT_MIN(rows, cols);
  field     alpha, beta;
  real      norm2, maxnorm2;
  uint      j, k, jmax, rows1, cols1, lwork;

  /* Provide storage for scaling factors */
  if (tau->dim < refl)
    resize_avector(tau, refl);
  tauv = tau->v;

  /* Provide working storage for h2_larf */
  lwork = cols;
  work = allocfield(lwork);

  /* Initialize column pivots */
  if (colpiv) {
    for (j = 0; j < cols; j++)
      colpiv[j] = j;
  }

  for (k = 0; k < refl; k++) {
    /* Compute norm of k-th column */
    rows1 = rows - k;
    norm2 = REAL_SQR(h2_nrm2(&rows1, aa + k + k * lda, &u_one));
    maxnorm2 = norm2;
    jmax = k;

    /* Find maximal norm */
    for (j = k + 1; j < cols; j++) {
      norm2 = REAL_SQR(h2_nrm2(&rows1, aa + k + j * lda, &u_one));
      if (norm2 > maxnorm2) {
	maxnorm2 = norm2;
	jmax = j;
      }
    }

    /* Swap columns */
    if (jmax != k) {
      h2_swap(&rows, aa + k * lda, &u_one, aa + jmax * lda, &u_one);

      if (colpiv) {
	j = colpiv[k];
	colpiv[k] = colpiv[jmax];
	colpiv[jmax] = j;
      }
    }

    /* Determine Householder reflection vector v */
    alpha = aa[k + k * lda];
    h2_larfg(&rows1, &alpha, aa + k + 1 + k * lda, &u_one, &beta);

    /* Store scaling factor */
    tauv[k] = beta;

    /* Update remaining columns */
    cols1 = cols - k - 1;
    aa[k + k * lda] = 1.0;
    beta = CONJ(beta);
    h2_larf(_h2_left, &rows1, &cols1, aa + k + k * lda, &u_one, &beta,
	    aa + k + (k + 1) * lda, &a->ld, work);

    /* Complete k-th column */
    aa[k + k * lda] = alpha;
  }

  freemem(work);

  return refl;
}

uint
qrdecomp_rank_amatrix(pamatrix a, pavector tau, pctruncmode tm, real eps,
		      uint * colpiv)
{
  pfield    aa = a->a;
  longindex lda = a->ld;
  pfield    tauv;
  pfield    work;
  uint      rows = a->rows;
  uint      cols = a->cols;
  uint      refl = UINT_MIN(rows, cols);
  field     alpha, beta;
  real      norm2, maxnorm2, frobnorm2, firstnorm2;
  uint      j, k, jmax, rows1, cols1, lwork;

  /* Provide enough storage for scaling factors */
  if (tau->dim < refl)
    resize_avector(tau, refl);
  tauv = tau->v;

  /* Provide working storage for h2_larf */
  lwork = cols;
  work = allocfield(lwork);

  /* Initialize column pivots */
  if (colpiv) {
    for (j = 0; j < cols; j++)
      colpiv[j] = j;
  }

  /* Initialize norm estimate for relative error criterion */
  firstnorm2 = 0.0;

  for (k = 0; k < refl; k++) {
    /* Compute norm of k-th column */
    rows1 = rows - k;
    norm2 = REAL_SQR(h2_nrm2(&rows1, aa + k + k * lda, &u_one));
    maxnorm2 = norm2;
    frobnorm2 = norm2;
    jmax = k;

    /* Find maximal norm */
    for (j = k + 1; j < cols; j++) {
      norm2 = REAL_SQR(h2_nrm2(&rows1, aa + k + j * lda, &u_one));
      if (norm2 > maxnorm2) {
	maxnorm2 = norm2;
	jmax = j;
      }
      frobnorm2 += norm2;
    }

    /* Swap columns */
    if (jmax != k) {
      h2_swap(&rows, aa + k * lda, &u_one, aa + jmax * lda, &u_one);

      if (colpiv) {
	j = colpiv[k];
	colpiv[k] = colpiv[jmax];
	colpiv[jmax] = j;
      }
    }

    /* Prepare norm */
    norm2 = maxnorm2;

    /* Exit if the norm is small enough */
    if (tm && tm->absolute) {
      /* Check the Frobenius norm */
      if (frobnorm2 <= eps * eps)
	break;
    }
    else {
      if (k == 0)
	/* Compute Frobenius norm or estimate spectral norm */
	firstnorm2 = ((tm && tm->frobenius) ? frobnorm2 : norm2);
      else
	/* Compare Frobenius norm with estimate for the entire matrix */
      if (frobnorm2 <= eps * eps * firstnorm2)
	break;
    }

    /* Determine Householder reflection vector v */
    alpha = aa[k + k * lda];
    h2_larfg(&rows1, &alpha, aa + k + 1 + k * lda, &u_one, &beta);

    /* Store scaling factor */
    tauv[k] = beta;

    /* Update remaining columns */
    cols1 = cols - k - 1;
    aa[k + k * lda] = 1.0;
    beta = CONJ(beta);
    h2_larf(_h2_left, &rows1, &cols1, aa + k + k * lda, &u_one, &beta,
	    aa + k + (k + 1) * lda, &a->ld, work);

    /* Complete k-th column */
    aa[k + k * lda] = alpha;
  }

  freemem(work);

  shrink_avector(tau, k);

  return k;
}
#else
uint
qrdecomp_pivot_amatrix(pamatrix a, pavector tau, uint * colpiv)
{
  pfield    aa = a->a;
  longindex lda = a->ld;
  pfield    tauv;
  uint      rows = a->rows;
  uint      cols = a->cols;
  uint      refl = UINT_MIN(rows, cols);
  field     alpha, beta, gamma, diag;
  real      norm2, norm, maxnorm2, adiag;
  uint      i, j, k, jmax;

  /* Provide enough storage for scaling factors */
  if (tau->dim < refl)
    resize_avector(tau, refl);
  tauv = tau->v;

  /* Initialize column pivots */
  if (colpiv) {
    for (j = 0; j < cols; j++)
      colpiv[j] = j;
  }

  for (k = 0; k < refl; k++) {
    /* Compute norm of k-th column */
    norm2 = 0.0;
    for (i = k; i < rows; i++)
      norm2 += ABSSQR(aa[i + k * lda]);
    maxnorm2 = norm2;
    jmax = k;

    /* Find maximal norm */
    for (j = k + 1; j < cols; j++) {
      norm2 = 0.0;
      for (i = k; i < rows; i++)
	norm2 += ABSSQR(aa[i + j * lda]);
      if (norm2 > maxnorm2) {
	maxnorm2 = norm2;
	jmax = j;
      }
    }

    /* Swap columns */
    if (jmax != k) {
      for (i = 0; i < rows; i++) {
	alpha = aa[i + k * lda];
	aa[i + k * lda] = aa[i + jmax * lda];
	aa[i + jmax * lda] = alpha;
      }

      if (colpiv) {
	j = colpiv[k];
	colpiv[k] = colpiv[jmax];
	colpiv[jmax] = j;
      }
    }

    /* Remember pivot
       if (colpiv)
       colpiv[k] = jmax;
     */

    /* Prepare norm */
    norm2 = maxnorm2;
    norm = REAL_SQRT(norm2);

    if (norm2 <= H2_ALMOST_ZERO)
      tauv[k] = 0.0;
    else {
      /* Determine reflection vector v */
      diag = aa[k + k * lda];
      adiag = ABS(diag) + norm;
      alpha = -SIGN1(diag) * norm;

      /* Compute norm of v */
      beta = 1.0 / (norm2 + ABS(diag) * norm);

      /* Rescale to ensure v_1 = 1 */
      beta *= adiag * adiag;
      gamma = CONJ(SIGN1(diag)) / adiag;
      for (i = k + 1; i < rows; i++)
	aa[i + k * lda] *= gamma;
      tauv[k] = beta;

      /* Compute k-th column */
      aa[k + k * lda] = alpha;

      /* Update remaining columns */
      for (j = k + 1; j < cols; j++) {
	gamma = aa[k + j * lda];
	for (i = k + 1; i < rows; i++)
	  gamma += CONJ(aa[i + k * lda]) * aa[i + j * lda];

	gamma *= beta;

	aa[k + j * lda] -= gamma;
	for (i = k + 1; i < rows; i++)
	  aa[i + j * lda] -= gamma * aa[i + k * lda];
      }
    }
  }

  return refl;
}

uint
qrdecomp_rank_amatrix(pamatrix a, pavector tau, pctruncmode tm, real eps,
		      uint * colpiv)
{
  pfield    aa = a->a;
  longindex lda = a->ld;
  pfield    tauv;
  uint      rows = a->rows;
  uint      cols = a->cols;
  uint      refl = UINT_MIN(rows, cols);
  field     alpha, beta, gamma, diag;
  real      norm2, norm, maxnorm2, frobnorm2, firstnorm2, adiag;
  uint      i, j, k, jmax;

  /* Provide enough storage for scaling factors */
  if (tau->dim < refl)
    resize_avector(tau, refl);
  tauv = tau->v;

  /* Initialize column pivots */
  if (colpiv) {
    for (j = 0; j < cols; j++)
      colpiv[j] = j;
  }

  /* Initialize norm estimate for relative error criterion */
  firstnorm2 = 0.0;

  for (k = 0; k < refl; k++) {
    /* Compute norm of k-th column */
    norm2 = 0.0;
    for (i = k; i < rows; i++)
      norm2 += ABSSQR(aa[i + k * lda]);
    maxnorm2 = norm2;
    frobnorm2 = norm2;
    jmax = k;

    /* Find maximal norm */
    for (j = k + 1; j < cols; j++) {
      norm2 = 0.0;
      for (i = k; i < rows; i++)
	norm2 += ABSSQR(aa[i + j * lda]);
      if (norm2 > maxnorm2) {
	maxnorm2 = norm2;
	jmax = j;
      }
      frobnorm2 += norm2;
    }

    /* Swap columns */
    if (jmax != k) {
      for (i = 0; i < rows; i++) {
	alpha = aa[i + k * lda];
	aa[i + k * lda] = aa[i + jmax * lda];
	aa[i + jmax * lda] = alpha;
      }

      if (colpiv) {
	j = colpiv[k];
	colpiv[k] = jmax;
	colpiv[jmax] = j;
      }
    }

    /* Prepare norm */
    norm2 = maxnorm2;
    norm = REAL_SQRT(norm2);

    /* Exit if the norm is small enough */
    if (tm && tm->absolute) {
      /* Check the Frobenius norm */
      if (frobnorm2 <= eps * eps)
	break;
    }
    else {
      if (k == 0)
	/* Compute Frobenius norm or estimate spectral norm */
	firstnorm2 = ((tm && tm->frobenius) ? frobnorm2 : norm2);
      else
	/* Compare Frobenius norm with estimate for the entire matrix */
      if (frobnorm2 <= eps * eps * firstnorm2)
	break;
    }

    if (norm2 == 0.0)
      tauv[k] = 0.0;
    else {
      /* Determine reflection vector v */
      diag = aa[k + k * lda];
      adiag = ABS(diag) + norm;
      alpha = -SIGN1(diag) * norm;

      /* Compute norm of v */
      beta = 1.0 / (norm2 + ABS(diag) * norm);

      /* Rescale to ensure v_1 = 1 */
      beta *= adiag * adiag;
      gamma = CONJ(SIGN1(diag)) / adiag;
      for (i = k + 1; i < rows; i++)
	aa[i + k * lda] *= gamma;
      tauv[k] = beta;

      /* Compute k-th column */
      aa[k + k * lda] = alpha;

      /* Update remaining columns */
      for (j = k + 1; j < cols; j++) {
	gamma = aa[k + j * lda];
	for (i = k + 1; i < rows; i++)
	  gamma += CONJ(aa[i + k * lda]) * aa[i + j * lda];

	gamma *= beta;

	aa[k + j * lda] -= gamma;
	for (i = k + 1; i < rows; i++)
	  aa[i + j * lda] -= gamma * aa[i + k * lda];
      }
    }
  }

  shrink_avector(tau, k);

  return k;
}
#endif

#ifdef USE_BLAS
/* Remark: if compiled the wrong way, DORMQR, and by extension DORMBR
 * and DGESVD, are currently not thread-safe.
 * gfortran does the right thing if called with "-frecursive", but this
 * appears not to be the standard in, e.g., OpenSUSE Linux. */
#if defined(THREADSAFE_LAPACK) || !defined(USE_OPENMP)
void
qreval_amatrix_avector(bool qtrans, pcamatrix a, pcavector tau, pavector x)
{
  uint      rows = a->rows;
  uint      cols = a->cols;
  uint      refl;
  field     work[4];
  int       lwork, info;

  refl = UINT_MIN3(rows, cols, tau->dim);

  if (refl < 1)
    return;

  assert(x->dim >= rows);

  lwork = 4;

  if (qtrans) {
    h2_ormqr(_h2_left, _h2_adj,
	     &rows, &u_one, &refl,
	     a->a, &a->ld, tau->v, x->v, &x->dim, work, &lwork, &info);
    assert(info == 0);
  }
  else {
    h2_ormqr(_h2_left, _h2_ntrans,
	     &rows, &u_one, &refl,
	     a->a, &a->ld, tau->v, x->v, &x->dim, work, &lwork, &info);
    assert(info == 0);
  }
}

void
qreval_amatrix(bool qtrans, pcamatrix a, pcavector tau, pamatrix x)
{
  uint      rows = a->rows;
  uint      cols = a->cols;
  uint      refl;
  pfield    work;
  int       lwork, info;

  refl = UINT_MIN3(rows, cols, tau->dim);

  if (refl < 1 || x->cols < 1)
    return;

  assert(x->rows >= rows);

  lwork = 4 * x->cols;
  work = allocfield(lwork);

  if (qtrans) {
    h2_ormqr(_h2_left, _h2_adj,
	     &rows, &x->cols, &refl,
	     a->a, &a->ld, tau->v, x->a, &x->ld, work, &lwork, &info);
    assert(info == 0);
  }
  else {
    h2_ormqr(_h2_left, _h2_ntrans,
	     &rows, &x->cols, &refl,
	     a->a, &a->ld, tau->v, x->a, &x->ld, work, &lwork, &info);
    assert(info == 0);
  }

  freemem(work);
}
#else
void
qreval_amatrix_avector(bool qtrans, pcamatrix a, pcavector tau, pavector x)
{
  pfield    aa = a->a;
  field     work[1];
  longindex lda = a->ld;
  pcfield   tauv = tau->v;
  pfield    xv = x->v;
  uint      rows = a->rows;
  uint      cols = a->cols;
  uint      refl, rows1;
  field     beta, gamma;
  uint      k;

  refl = UINT_MIN3(rows, cols, tau->dim);

  assert(x->dim >= rows);

  if (qtrans) {
    for (k = 0; k < refl; k++) {
      rows1 = rows - k;
      beta = CONJ(tauv[k]);
      gamma = aa[k + k * lda];
      aa[k + k * lda] = 1.0;
      h2_larf(_h2_left, &rows1, &u_one, aa + k + k * lda, &u_one, &beta,
	      xv + k, &rows, work);
      aa[k + k * lda] = gamma;
    }
  }
  else {
    for (k = refl; k-- > 0;) {
      rows1 = rows - k;
      beta = tauv[k];
      gamma = aa[k + k * lda];
      aa[k + k * lda] = 1.0;
      h2_larf(_h2_left, &rows1, &u_one, aa + k + k * lda, &u_one, &beta,
	      xv + k, &rows, work);
      aa[k + k * lda] = gamma;
    }
  }
}

void
qreval_amatrix(bool qtrans, pcamatrix a, pcavector tau, pamatrix x)
{
  pfield    aa = a->a;
  pfield    work;
  longindex lda = a->ld;
  pcfield   tauv = tau->v;
  pfield    xa = x->a;
  uint      rows = a->rows;
  uint      cols = a->cols;
  uint      refl, rows1;
  field     beta, gamma;
  uint      k;

  refl = UINT_MIN3(rows, cols, tau->dim);

  assert(x->rows >= rows);

  work = allocfield(UINT_MAX(x->rows, x->cols));

  if (qtrans) {
    for (k = 0; k < refl; k++) {
      rows1 = rows - k;
      beta = CONJ(tauv[k]);
      gamma = aa[k + k * lda];
      aa[k + k * lda] = 1.0;
      h2_larf(_h2_left, &rows1, &x->cols, aa + k + k * lda, &u_one, &beta,
	      xa + k, &x->ld, work);
      aa[k + k * lda] = gamma;
    }
  }
  else {
    for (k = refl; k-- > 0;) {
      rows1 = rows - k;
      beta = tauv[k];
      gamma = aa[k + k * lda];
      aa[k + k * lda] = 1.0;
      h2_larf(_h2_left, &rows1, &x->cols, aa + k + k * lda, &u_one, &beta,
	      xa + k, &x->ld, work);
      aa[k + k * lda] = gamma;
    }
  }

  freemem(work);
}
#endif
#else
void
qreval_amatrix_avector(bool qtrans, pcamatrix a, pcavector tau, pavector x)
{
  pcfield   aa = a->a;
  longindex lda = a->ld;
  pcfield   tauv = tau->v;
  pfield    xv = x->v;
  uint      rows = a->rows;
  uint      cols = a->cols;
  uint      refl;
  field     beta, gamma;
  uint      i, k;

  refl = UINT_MIN3(rows, cols, tau->dim);

  assert(x->dim >= rows);

  if (qtrans) {
    for (k = 0; k < refl; k++) {
      beta = CONJ(tauv[k]);

      if (beta != 0.0) {
	gamma = xv[k];
	for (i = k + 1; i < rows; i++)
	  gamma += CONJ(aa[i + k * lda]) * xv[i];

	gamma *= beta;

	xv[k] -= gamma;
	for (i = k + 1; i < rows; i++)
	  xv[i] -= gamma * aa[i + k * lda];
      }
    }
  }
  else {
    for (k = refl; k-- > 0;) {
      beta = tauv[k];

      if (beta != 0.0) {
	gamma = xv[k];
	for (i = k + 1; i < rows; i++)
	  gamma += CONJ(aa[i + k * lda]) * xv[i];

	gamma *= beta;

	xv[k] -= gamma;
	for (i = k + 1; i < rows; i++)
	  xv[i] -= gamma * aa[i + k * lda];
      }
    }
  }
}

void
qreval_amatrix(bool qtrans, pcamatrix a, pcavector tau, pamatrix x)
{
  pcfield   aa = a->a;
  longindex lda = a->ld;
  pcfield   tauv = tau->v;
  pfield    xa = x->a;
  longindex ldx = x->ld;
  uint      rows = a->rows;
  uint      cols = a->cols;
  uint      refl;
  field     beta, gamma;
  uint      i, j, k;

  refl = UINT_MIN3(rows, cols, tau->dim);

  assert(x->rows >= rows);

  if (qtrans) {
    for (k = 0; k < refl; k++) {
      beta = CONJ(tauv[k]);

      if (beta != 0.0)
	for (j = 0; j < x->cols; j++) {
	  gamma = xa[k + j * ldx];
	  for (i = k + 1; i < rows; i++)
	    gamma += CONJ(aa[i + k * lda]) * xa[i + j * ldx];

	  gamma *= beta;

	  xa[k + j * ldx] -= gamma;
	  for (i = k + 1; i < rows; i++)
	    xa[i + j * ldx] -= gamma * aa[i + k * lda];
	}
    }
  }
  else {
    for (k = refl; k-- > 0;) {
      beta = tauv[k];

      if (beta != 0.0)
	for (j = 0; j < x->cols; j++) {
	  gamma = xa[k + j * ldx];
	  for (i = k + 1; i < rows; i++)
	    gamma += CONJ(aa[i + k * lda]) * xa[i + j * ldx];

	  gamma *= beta;

	  xa[k + j * ldx] -= gamma;
	  for (i = k + 1; i < rows; i++)
	    xa[i + j * ldx] -= gamma * aa[i + k * lda];
	}
    }
  }
}
#endif

void
qrsolve_amatrix_avector(pcamatrix a, pcavector tau, pavector x)
{
  qreval_amatrix_avector(true, a, tau, x);
  triangularsolve_amatrix_avector(false, false, false, a, x);
}

void
qrinvert_amatrix(pamatrix a)
{
  pamatrix  acopy;
  pavector  tau, b;
  amatrix   atmp;
  avector   ttmp, btmp;
  uint      n = a->rows;
  uint      i;

  assert(n == a->cols);

  acopy = init_amatrix(&atmp, n, n);
  tau = init_avector(&ttmp, n);

  copy_amatrix(false, a, acopy);

  qrdecomp_amatrix(acopy, tau);

  identity_amatrix(a);

  for (i = 0; i < n; i++) {
    b = init_column_avector(&btmp, a, i);
    qrsolve_amatrix_avector(acopy, tau, b);
    uninit_avector(b);
  }

  uninit_avector(tau);
  uninit_amatrix(acopy);
}

#ifdef USE_BLAS
void
qrexpand_amatrix(pcamatrix a, pcavector tau, pamatrix q)
{
  field    *work;
  uint      refl, lwork, info;

  refl = UINT_MIN(UINT_MIN(q->cols, tau->dim), UINT_MIN(a->rows, a->cols));

  /* Quick exit if no reflections used */
  if (refl == 0) {
    identity_amatrix(q);
    return;
  }

  copy_sub_amatrix(false, a, q);

  lwork = 4 * a->rows;
  work = allocfield(lwork);

  h2_orgqr(&q->rows, &q->cols, &refl, q->a, &q->ld, tau->v, work, &lwork,
	   &info);
  assert(info == 0);

  freemem(work);
}
#else
void
qrexpand_amatrix(pcamatrix a, pcavector tau, pamatrix q)
{
  pcfield   aa = a->a;
  longindex lda = a->ld;
  pcfield   tauv = tau->v;
  pfield    qa = q->a;
  longindex ldq = q->ld;
  uint      rows = a->rows;
  uint      refl;
  field     beta, gamma;
  uint      i, j, k;

  /* Determine number of relevant elementary reflections */
  refl = UINT_MIN(UINT_MIN(q->cols, tau->dim), UINT_MIN(a->rows, a->cols));

  assert(q->rows >= rows);
  assert(q->cols <= rows);

  /* Create identity matrix */
  for (j = 0; j < q->cols; j++) {
    for (i = 0; i < rows; i++)
      qa[i + j * ldq] = 0.0;
    qa[j + j * ldq] = 1.0;
  }

  /* Apply reflections in reversed order */
  for (k = refl; k-- > 0;) {
    beta = tauv[k];
    if (beta != 0.0)
      for (j = k; j < q->cols; j++) {
	gamma = qa[k + j * ldq];
	for (i = k + 1; i < rows; i++)
	  gamma += CONJ(aa[i + k * lda]) * qa[i + j * ldq];

	gamma *= beta;

	qa[k + j * ldq] -= gamma;
	for (i = k + 1; i < rows; i++)
	  qa[i + j * ldq] -= gamma * aa[i + k * lda];
      }
  }
}
#endif
