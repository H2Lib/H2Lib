/* ------------------------------------------------------------
 This is the file "aca.c" of the H2Lib package.
 All rights reserved, Sven Christophersen 2012
 ------------------------------------------------------------ */

/* C STD LIBRARY */
/* CORE 0 */
#include "basic.h"
/* CORE 1 */
#include "aca.h"
/* CORE 2 */
/* CORE 3 */
/* SIMPLE */
/* PARTICLES */
/* BEM */

void
decomp_fullaca_rkmatrix(pamatrix A, const real accur, uint ** ridx,
			uint ** cidx, prkmatrix R)
{
  uint      rows = A->rows;
  uint      cols = A->cols;

  uint      k, i, j, imax, jmax, ld, ldc, ldd, maxrank;
  uint     *rperm, *cperm, *rpivot, *cpivot;
  field     Aij, y;
  real      absAij, absy, tol;
  pfield    a, c, d;

  maxrank = UINT_MIN(rows, cols);

  rperm = allocmem(maxrank * sizeof(uint));
  cperm = allocmem(maxrank * sizeof(uint));
  rpivot = allocmem(rows * sizeof(uint));
  cpivot = allocmem(cols * sizeof(uint));

  for (i = 0; i < rows; ++i) {
    rpivot[i] = i;
  }

  for (j = 0; j < cols; ++j) {
    cpivot[j] = j;
  }

  k = 0;
  a = A->a;
  ld = A->ld;
  absAij = 0.0;
  imax = 0;
  jmax = 0;

  /* find maximal element */
  for (j = 0; j < cols; j++) {
    for (i = 0; i < rows; i++) {
      absy = ABSSQR(a[i + j * ld]);
      if (absy > absAij) {
	absAij = absy;
	imax = i;
	jmax = j;
      }
    }
  }

  tol = absAij * accur * accur;

  while (k < maxrank && absAij > tol) {

    Aij = 1.0 / a[imax + jmax * ld];

    i = rpivot[k];
    rpivot[k] = rpivot[imax];
    rpivot[imax] = i;

    j = cpivot[k];
    cpivot[k] = cpivot[jmax];
    cpivot[jmax] = j;

    rperm[k] = imax;
    cperm[k] = jmax;

    for (i = 0; i < rows; i++) {
      y = a[i + jmax * ld];
      a[i + jmax * ld] = a[i + k * ld];
      a[i + k * ld] = y;
    }
    for (j = 0; j < cols; j++) {
      y = a[imax + j * ld];
      a[imax + j * ld] = a[k + j * ld];
      a[k + j * ld] = y;
    }

    for (i = k + 1; i < rows; i++) {
      a[i + k * ld] *= Aij;
    }

    for (j = k + 1; j < cols; j++) {
      for (i = k + 1; i < rows; i++) {
	a[i + j * ld] -= a[i + k * ld] * a[k + j * ld];
      }
    }

    k++;
    absAij = 0.0;

    /* find maximal element */
    for (j = k; j < cols; j++) {
      for (i = k; i < rows; i++) {
	absy = ABSSQR(a[i + j * ld]);
	if (absy > absAij) {
	  absAij = absy;
	  imax = i;
	  jmax = j;
	}
      }
    }
  }

  resize_rkmatrix(R, rows, cols, k);

  c = R->A.a;
  d = R->B.a;
  ldc = R->A.ld;
  ldd = R->B.ld;

  /* copy in C and D */
  for (j = 0; j < k; j++) {
    for (i = 0; i < j; i++) {
      c[i + j * ldc] = 0.0;
    }
    c[j + j * ldc] = 1.0;
    for (i = j + 1; i < rows; i++) {
      c[i + j * ldc] = a[i + j * ld];
    }
  }
  for (j = 0; j < k; j++) {
    for (i = 0; i < j; i++) {
      d[i + j * ldd] = 0.0;
    }
    for (i = j; i < cols; i++) {
      d[i + j * ldd] = CONJ(a[j + i * ld]);
    }
  }

  /* reverse pivot permutations */
  for (i = k; i-- > 0;) {
    for (j = 0; j < k; j++) {
      y = c[i + j * ldc];
      c[i + j * ldc] = c[rperm[i] + j * ldc];
      c[rperm[i] + j * ldc] = y;
      y = d[i + j * ldd];
      d[i + j * ldd] = d[cperm[i] + j * ldd];
      d[cperm[i] + j * ldd] = y;
    }
  }

  if (ridx != NULL) {
    *ridx = allocuint(k);
    for (i = 0; i < k; ++i) {
      (*ridx)[i] = rpivot[i];
    }
  }

  if (cidx != NULL) {
    *cidx = allocuint(k);
    for (i = 0; i < k; ++i) {
      (*cidx)[i] = cpivot[i];
    }
  }

  freemem(rperm);
  freemem(cperm);
  freemem(rpivot);
  freemem(cpivot);

}

void
decomp_partialaca_rkmatrix(matrixentry_t entry, void *data,
			   const uint * ridx, const uint rows,
			   const uint * cidx, const uint cols, real accur,
			   uint ** rpivot, uint ** cpivot, prkmatrix R)
{
  pamatrix  A, B;
  amatrix   A_k, B_k;
  uint     *rperm, *cperm, *rpiv, *cpiv;
  uint      i, j, mu, k, i_k, j_k;
  real      error, error2, starterror, M;
  field     Aij;
  pfield    aa, bb;

  k = 0;

  rperm = allocmem(rows * sizeof(uint));
  cperm = allocmem(cols * sizeof(uint));
  rpiv = allocmem(rows * sizeof(uint));
  cpiv = allocmem(cols * sizeof(uint));

  for (i = 0; i < rows; ++i) {
    rpiv[i] = ridx[i];
  }

  for (j = 0; j < cols; ++j) {
    cpiv[j] = cidx[j];
  }

  A = &R->A;
  B = &R->B;

  error = 1.0;
  starterror = 1.0;

  while (error > accur && k < rows && k < cols) {
    resizecopy_amatrix(A, rows, k + 1);
    resizecopy_amatrix(B, cols, k + 1);
    aa = A->a;
    bb = B->a;

    (void) init_sub_amatrix(&A_k, A, rows - k, k, 1, k);
    (void) init_sub_amatrix(&B_k, B, cols - k, k, 1, k);

    if (k > 0) {
      /* Compute next i_k. */
      M = ABSSQR(aa[k + (k - 1) * rows]);
      i_k = k;
      for (i = k + 1; i < rows; ++i) {
	if (ABSSQR(aa[i + (k - 1) * rows]) > M) {
	  M = ABSSQR(aa[i + (k - 1) * rows]);
	  i_k = i;
	}
      }
    }
    else {
      i_k = rows / 2;
    }

    /* Get next column for B. */
    entry(rpiv + i_k, cpiv + k, data, true, &B_k);
    conjugate_amatrix(&B_k);

    /* Subtract current rank-k-approximation. */
    for (j = 0; j < k; ++j) {
      bb[j + k * cols] = 0.0;
    }
    for (mu = 0; mu < k; ++mu) {
      Aij = aa[i_k + mu * rows];
      for (j = k; j < cols; ++j) {
	bb[j + k * cols] = bb[j + k * cols] - bb[j + mu * cols] * Aij;
      }
    }

    /* Compute next j_k. */
    M = ABSSQR(bb[k + k * cols]);
    j_k = k;
    for (j = k + 1; j < cols; ++j) {
      if (ABSSQR(bb[j + k * cols]) > M) {
	M = ABSSQR(bb[j + k * cols]);
	j_k = j;
      }
    }

    /* Get next column for A. */
    entry(rpiv + k, cpiv + j_k, data, false, &A_k);

    /* Subtract current rank-k-approximation. */
    for (i = 0; i < k; ++i) {
      aa[i + k * rows] = 0.0;
    }
    for (mu = 0; mu < k; ++mu) {
      Aij = bb[j_k + mu * cols];
      for (i = k; i < rows; ++i) {
	aa[i + k * rows] = aa[i + k * rows] - aa[i + mu * rows] * Aij;
      }
    }

    Aij = 1.0 / bb[j_k + k * cols];
    for (i = 0; i < rows; ++i) {
      aa[i + k * rows] = aa[i + k * rows] * Aij;
    }

    /* Update permutations. */
    i = rpiv[k];
    rpiv[k] = rpiv[i_k];
    rpiv[i_k] = i;
    rperm[k] = i_k;

    j = cpiv[k];
    cpiv[k] = cpiv[j_k];
    cpiv[j_k] = j;
    cperm[k] = j_k;

    k++;

    /* Apply permutation to current approximation */
    for (i = 0; i < k; ++i) {
      Aij = aa[i_k + i * rows];
      aa[i_k + i * rows] = aa[k - 1 + i * rows];
      aa[k - 1 + i * rows] = Aij;
    }

    for (j = 0; j < k; ++j) {
      Aij = bb[j_k + j * cols];
      bb[j_k + j * cols] = bb[k - 1 + j * cols];
      bb[k - 1 + j * cols] = Aij;
    }

    if (k == 1) {
      /* Computation of starterror. */
      starterror = 0.0;
      for (i = k - 1; i < rows; ++i) {
	starterror += ABSSQR(aa[i]);
      }
      error = 0.0;
      for (j = k - 1; j < cols; ++j) {
	error += ABSSQR(bb[j]);
      }
      starterror = REAL_RSQRT(starterror * error);
    }

    /* Computation of current relative error. */
    error = 0.0;
    for (i = k - 1; i < rows; ++i) {
      error += ABSSQR(aa[i + (k - 1) * rows]);
    }
    error2 = 0.0;
    for (j = k - 1; j < cols; ++j) {
      error2 += ABSSQR(bb[j + (k - 1) * cols]);
    }
    error = REAL_SQRT(error * error2) * starterror;

    uninit_amatrix(&A_k);
    uninit_amatrix(&B_k);
  }

  R->k = k;

  /* Reverse pivot permutations */
  aa = A->a;
  bb = B->a;
  for (i = k; i-- > 0;) {
    for (j = 0; j < k; j++) {
      Aij = aa[i + j * rows];
      aa[i + j * rows] = aa[rperm[i] + j * rows];
      aa[rperm[i] + j * rows] = Aij;

      Aij = bb[i + j * cols];
      bb[i + j * cols] = bb[cperm[i] + j * cols];
      bb[cperm[i] + j * cols] = Aij;
    }
  }

  conjugate_amatrix(B);

  if (rpivot != NULL) {
    *rpivot = allocuint(k);
    for (i = 0; i < k; ++i) {
      (*rpivot)[i] = rpiv[i];
    }
  }

  if (cpivot != NULL) {
    *cpivot = allocuint(k);
    for (i = 0; i < k; ++i) {
      (*cpivot)[i] = cpiv[i];
    }
  }

  freemem(rperm);
  freemem(cperm);
  freemem(rpiv);
  freemem(cpiv);
}

void
copy_lower_aca_amatrix(bool unit, pcamatrix A, uint * xi, pamatrix B)
{
  pcfield   a = A->a;
  pfield    b = B->a;
  const uint lda = A->ld;
  const uint ldb = B->ld;
  const uint rows = B->rows;
  const uint cols = A->cols;

  uint      i, j, ii;

  assert(A->rows >= rows);
  assert(B->cols == cols);

  for (j = 0; j < cols; ++j) {
    if (unit) {
      b[j + j * ldb] = 1.0;
    }
    else {
      ii = xi[j];
      assert(ii < A->rows);
      b[j + j * ldb] = a[ii + j * lda];
    }
    for (i = j + 1; i < rows; ++i) {
      ii = xi[i];
      assert(ii < A->rows);
      b[i + j * ldb] = a[ii + j * lda];
    }
  }
}

void
copy_upper_aca_amatrix(bool unit, pcamatrix A, uint * xi, pamatrix B)
{
  pcfield   a = A->a;
  pfield    b = B->a;
  const uint lda = A->ld;
  const uint ldb = B->ld;
  const uint rows = B->rows;
  const uint cols = A->cols;

  uint      i, j, jj;

  assert(A->rows >= rows);
  assert(B->rows == cols);

  for (j = 0; j < rows; ++j) {
    assert(xi[j] < A->rows);
    jj = xi[j];
    for (i = 0; i < j; ++i) {
      b[i + j * ldb] = a[jj + i * lda];
    }
    if (unit) {
      b[i + j * ldb] = 1.0;
    }
    else {
      b[i + j * ldb] = a[jj + i * lda];
    }
  }
}
