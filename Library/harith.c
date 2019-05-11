/* ------------------------------------------------------------
 * This is the file "harith.c" of the H2Lib package.
 * All rights reserved, Steffen Boerm 2014
 * ------------------------------------------------------------ */

#include "harith.h"
#include "basic.h"
#include "eigensolvers.h"
#include "factorizations.h"

/* Quick exit if a rank zero matrix is added? (handled in options.inc) */
/* #define HARITH_RKMATRIX_QUICK_EXIT */

/* Quick exit if a dense zero matrix is added? (handled in options.inc) */
/* #define HARITH_AMATRIX_QUICK_EXIT */

/* ------------------------------------------------------------
 * Truncation of an rkmatrix
 * ------------------------------------------------------------ */

/* First version: apply SVD directly to A B^*.
 * This approach is only advisable if the rank is high compared to
 * the number of rows and columns. */
static void
trunc_ab_rkmatrix(pctruncmode tm, real eps, prkmatrix r)
{
  amatrix   tmp1, tmp2, tmp3;
  realavector tmp4;
  pamatrix  a, b, c, u, vt;
  prealavector sigma;
  uint      rows, cols;
  uint      k1, knew;

  assert(r->A.cols == r->k);
  assert(r->B.cols == r->k);

  rows = r->A.rows;
  cols = r->B.rows;
  a = &r->A;
  b = &r->B;

  /* Compute C = A B^* */
  c = init_amatrix(&tmp1, rows, cols);
  clear_amatrix(c);
  addmul_amatrix(1.0, false, a, true, b, c);
  k1 = UINT_MIN(rows, cols);

  /* Compute singular value decomposition */
  u = init_amatrix(&tmp2, rows, k1);
  vt = init_amatrix(&tmp3, k1, cols);
  sigma = init_realavector(&tmp4, k1);
  svd_amatrix(c, sigma, u, vt);

  /* Determine rank */
  knew = findrank_truncmode(tm, eps, sigma);

  /* Set new rank */
  setrank_rkmatrix(r, knew);

  /* Copy singular vectors */
  copy_sub_amatrix(false, u, &r->A);
  copy_sub_amatrix(true, vt, &r->B);

  /* Scale singular vectors */
  diageval_realavector_amatrix(1.0, true, sigma, true, &r->A);

  /* Clean up */
  uninit_realavector(sigma);
  uninit_amatrix(vt);
  uninit_amatrix(u);
  uninit_amatrix(c);
}

/* Second version: Turn A into an upper triangular matrix by a QR
 * factorization, then apply SVD to the resulting reduced matrix.
 * Advisable if the number of rows is significantly larger than the
 * rank. */
static void
trunc_qb_rkmatrix(pctruncmode tm, real eps, prkmatrix r)
{
  amatrix   tmp1, tmp2, tmp3, tmp4;
  avector   tmp5;
  realavector tmp6;
  pamatrix  a, b, c, a1, u, vt;
  pavector  tau;
  prealavector sigma;
  uint      rows, cols;
  uint      k, kr, k1, knew;

  assert(r->A.cols == r->k);
  assert(r->B.cols == r->k);

  rows = r->A.rows;
  cols = r->B.rows;
  k = r->k;
  b = &r->B;

  /* Copy factor A */
  a = init_amatrix(&tmp1, rows, k);
  copy_amatrix(false, &r->A, a);

  /* Compute QR factorization of A */
  tau = init_avector(&tmp5, k);
  qrdecomp_amatrix(a, tau);

  /* Overwrite B by C = B A^* (C^* = A B^*) */
  kr = UINT_MIN(rows, k);
  a1 = init_sub_amatrix(&tmp2, a, kr, 0, k, 0);
  triangulareval_amatrix(false, false, false, a1, true, b);
  uninit_amatrix(a1);
  c = init_sub_amatrix(&tmp4, b, cols, 0, kr, 0);

  /* Compute singular value decomposition */
  k1 = UINT_MIN(cols, kr);
  u = init_amatrix(&tmp2, cols, k1);
  vt = init_amatrix(&tmp3, k1, kr);
  sigma = init_realavector(&tmp6, k1);
  svd_amatrix(c, sigma, u, vt);

  /* Determine rank */
  knew = findrank_truncmode(tm, eps, sigma);

  /* Set new rank */
  setrank_rkmatrix(r, knew);

  /* Copy singular vectors */
  copy_sub_amatrix(false, u, &r->B);
  clear_amatrix(&r->A);
  copy_sub_amatrix(true, vt, &r->A);
  qreval_amatrix(false, a, tau, &r->A);
  uninit_avector(tau);

  /* Scale singular vectors */
  diageval_realavector_amatrix(1.0, true, sigma, true, &r->B);

  /* Clean up */
  uninit_realavector(sigma);
  uninit_amatrix(vt);
  uninit_amatrix(u);
  uninit_amatrix(c);
  uninit_amatrix(a);
}

/* Third version: Turn B into an upper triangular matrix by a QR
 * factorization, then apply SVD to the resulting reduced matrix.
 * Advisable if the number of columns is significantly larger than the
 * rank. */
static void
trunc_aq_rkmatrix(pctruncmode tm, real eps, prkmatrix r)
{
  amatrix   tmp1, tmp2, tmp3, tmp4;
  avector   tmp5;
  realavector tmp6;
  pamatrix  a, b, c, b1, u, vt;
  pavector  tau;
  prealavector sigma;
  uint      rows, cols;
  uint      k, kc, k1, knew;

  assert(r->A.cols == r->k);
  assert(r->B.cols == r->k);

  rows = r->A.rows;
  cols = r->B.rows;
  k = r->k;
  a = &r->A;

  /* Copy factor B */
  b = init_amatrix(&tmp1, cols, k);
  copy_amatrix(false, &r->B, b);

  /* Compute QR factorization of B */
  tau = init_avector(&tmp5, k);
  qrdecomp_amatrix(b, tau);

  /* Overwrite A by C = A B^* (C^* = B A^*) */
  kc = UINT_MIN(cols, k);
  b1 = init_sub_amatrix(&tmp2, b, kc, 0, k, 0);
  triangulareval_amatrix(false, false, false, b1, true, a);
  uninit_amatrix(b1);
  c = init_sub_amatrix(&tmp4, a, rows, 0, kc, 0);

  /* Compute singular value decomposition */
  k1 = UINT_MIN(rows, kc);
  u = init_amatrix(&tmp2, rows, k1);
  vt = init_amatrix(&tmp3, k1, kc);
  sigma = init_realavector(&tmp6, k1);
  svd_amatrix(c, sigma, u, vt);

  /* Determine rank */
  knew = findrank_truncmode(tm, eps, sigma);

  /* Set new rank */
  setrank_rkmatrix(r, knew);

  /* Copy singular vectors */
  copy_sub_amatrix(false, u, &r->A);
  clear_amatrix(&r->B);
  copy_sub_amatrix(true, vt, &r->B);
  qreval_amatrix(false, b, tau, &r->B);
  uninit_avector(tau);

  /* Scale singular vectors */
  diageval_realavector_amatrix(1.0, true, sigma, true, &r->A);

  /* Clean up */
  uninit_realavector(sigma);
  uninit_amatrix(vt);
  uninit_amatrix(u);
  uninit_amatrix(c);
  uninit_amatrix(b);
}

/* Fourth version: Reduce both A and B to upper triangular matrices
 * by QR factorizations, then compute SVD of the reduced matrix.
 * Advisable if both the number of rows and the number of columns are
 * larger than the rank. */
static void
trunc_qq_rkmatrix(pctruncmode tm, real eps, prkmatrix r)
{
  amatrix   tmp1, tmp2, tmp3, tmp4, tmp5;
  avector   tmp6, tmp7;
  realavector tmp8;
  pamatrix  a, b, c, a1, b1, u, vt;
  pavector  atau, btau;
  prealavector sigma;
  uint      rows, cols;
  uint      k, ak, bk, k1, knew;

  assert(r->A.cols == r->k);
  assert(r->B.cols == r->k);

  rows = r->A.rows;
  cols = r->B.rows;
  k = r->k;

  /* Copy factor A and B */
  a = init_amatrix(&tmp1, rows, k);
  copy_amatrix(false, &r->A, a);
  b = init_amatrix(&tmp2, cols, k);
  copy_amatrix(false, &r->B, b);

  /* Compute QR factorization Q_A R_A = A */
  atau = init_avector(&tmp6, k);
  qrdecomp_amatrix(a, atau);
  ak = UINT_MIN(k, rows);

  /* Compute QR factorization Q_B R_B = B */
  btau = init_avector(&tmp7, k);
  qrdecomp_amatrix(b, btau);
  bk = UINT_MIN(k, cols);

  /* Compute condensed matrix C = R_A R_B^* */
  c = init_amatrix(&tmp3, ak, bk);
  clear_amatrix(c);
  a1 = init_sub_amatrix(&tmp4, a, ak, 0, k, 0);
  b1 = init_sub_amatrix(&tmp5, b, bk, 0, k, 0);
  triangularaddmul_amatrix(1.0, false, false, a1, false, true, b1, c);
  uninit_amatrix(b1);
  uninit_amatrix(a1);

  /* Find singular value decomposition of Z */
  k1 = UINT_MIN(ak, bk);
  u = init_amatrix(&tmp4, ak, k1);
  vt = init_amatrix(&tmp5, k1, bk);
  sigma = init_realavector(&tmp8, k1);
  svd_amatrix(c, sigma, u, vt);

  /* Determine rank */
  knew = findrank_truncmode(tm, eps, sigma);

  /* Set new rank */
  setrank_rkmatrix(r, knew);

  /* Scale singular vectors */
  diageval_realavector_amatrix(1.0, true, sigma, true, u);

  /* Copy singular vectors */
  clear_amatrix(&r->A);
  copy_sub_amatrix(false, u, &r->A);
  clear_amatrix(&r->B);
  copy_sub_amatrix(true, vt, &r->B);

  /* Multiply singular vectors by Q_A and Q_B, respectively */
  qreval_amatrix(false, a, atau, &r->A);
  qreval_amatrix(false, b, btau, &r->B);

  /* Clean up */
  uninit_realavector(sigma);
  uninit_amatrix(vt);
  uninit_amatrix(u);
  uninit_amatrix(c);
  uninit_avector(btau);
  uninit_avector(atau);
  uninit_amatrix(b);
  uninit_amatrix(a);
}

void
trunc_rkmatrix(pctruncmode tm, real eps, prkmatrix r)
{
  uint      rows, cols, k;

  assert(r->A.cols == r->k);
  assert(r->B.cols == r->k);

#ifdef HARITH_RKMATRIX_QUICK_EXIT
  if (r->k == 0)
    return;
#endif

  rows = r->A.rows;
  cols = r->B.rows;
  k = r->k;

  /* Choose most efficient truncation algorithm */
  if (k < rows) {
    if (k < cols)
      /* rows and cols large, use QR decomposition for A and B */
      trunc_qq_rkmatrix(tm, eps, r);
    else
      /* rows large, cols small, use QR decomposition for A only */
      trunc_qb_rkmatrix(tm, eps, r);
  }
  else {
    if (k < cols)
      /* rows small, cols large, use QR decomposition for B only */
      trunc_aq_rkmatrix(tm, eps, r);
    else
      /* rows and cols small, no QR decomposition required */
      trunc_ab_rkmatrix(tm, eps, r);
  }
}

/* ------------------------------------------------------------
 * Truncated addition of an amatrix to an rkmatrix.
 * ------------------------------------------------------------ */

void
add_amatrix_rkmatrix(field alpha, bool atrans, pcamatrix a, pctruncmode tm,
		     real eps, prkmatrix b)
{
  amatrix   tmp1, tmp2, tmp3, tmp4;
  realavector tmp5;
  pamatrix  z, z1, u, vt, u1, vt1;
  prealavector sigma;
  pavector  tau;
  uint      rows, cols, k, knew;

#ifdef HARITH_AMATRIX_QUICK_EXIT
  if (normfrob2_amatrix(a) == 0.0) {
    return;
  }
#endif

  if (atrans) {
    assert(a->rows == b->B.rows);
    assert(a->cols == b->A.rows);
  }
  else {
    assert(a->rows == b->A.rows);
    assert(a->cols == b->B.rows);
  }

  rows = b->A.rows;
  cols = b->B.rows;

  if (rows > cols) {
    z = init_amatrix(&tmp1, rows, cols);

    /* Compute sum */
    if (atrans)
      copy_amatrix(true, a, z);
    else
      copy_amatrix(false, a, z);

    if (alpha != 1.0)
      scale_amatrix(alpha, z);

    addmul_amatrix(1.0, false, &b->A, true, &b->B, z);

    k = UINT_MIN(rows, cols);
    tau = new_avector(k);
    qrdecomp_amatrix(z, tau);
    z1 = new_amatrix(k, k);
    copy_upper_amatrix(z, false, z1);

    u = init_amatrix(&tmp2, k, k);
    vt = init_amatrix(&tmp3, k, k);
    sigma = init_realavector(&tmp5, k);
    svd_amatrix(z1, sigma, u, vt);

    /* Determine rank */
    knew = findrank_truncmode(tm, eps, sigma);

    /* Set new rank */
    setrank_rkmatrix(b, knew);

    /* Copy singular vectors */
    clear_amatrix(&b->A);
    u1 = init_sub_amatrix(&tmp4, u, k, 0, knew, 0);
    copy_sub_amatrix(false, u1, &b->A);
    qreval_amatrix(false, z, tau, &b->A);
    uninit_amatrix(u1);

    clear_amatrix(&b->B);
    vt1 = init_sub_amatrix(&tmp4, vt, knew, 0, k, 0);
    diageval_realavector_amatrix(1.0, false, sigma, false, vt1);
    copy_amatrix(true, vt1, &b->B);
    uninit_amatrix(vt1);

    del_avector(tau);
    del_amatrix(z1);
  }
  else if (cols > rows) {
    z = init_amatrix(&tmp1, cols, rows);

    /* Compute sum */
    if (atrans)
      copy_amatrix(false, a, z);
    else
      copy_amatrix(true, a, z);

    if (alpha != 1.0)
      scale_amatrix(alpha, z);

    addmul_amatrix(1.0, false, &b->B, true, &b->A, z);

    k = UINT_MIN(rows, cols);
    tau = new_avector(k);
    qrdecomp_amatrix(z, tau);
    z1 = new_amatrix(k, k);
    copy_upper_amatrix(z, false, z1);

    u = init_amatrix(&tmp2, k, k);
    vt = init_amatrix(&tmp3, k, k);
    sigma = init_realavector(&tmp5, k);
    svd_amatrix(z1, sigma, u, vt);

    /* Determine rank */
    knew = findrank_truncmode(tm, eps, sigma);

    /* Set new rank */
    setrank_rkmatrix(b, knew);

    /* Copy singular vectors */
    clear_amatrix(&b->A);
    vt1 = init_sub_amatrix(&tmp4, vt, knew, 0, k, 0);
    copy_amatrix(true, vt1, &b->A);
    uninit_amatrix(vt1);

    clear_amatrix(&b->B);
    u1 = init_sub_amatrix(&tmp4, u, k, 0, knew, 0);
    diageval_realavector_amatrix(1.0, true, sigma, true, u1);
    copy_sub_amatrix(false, u1, &b->B);
    qreval_amatrix(false, z, tau, &b->B);
    uninit_amatrix(u1);

    del_avector(tau);
    del_amatrix(z1);

  }
  else {
    z = init_amatrix(&tmp1, rows, cols);

    /* Compute sum */
    if (atrans)
      copy_amatrix(true, a, z);
    else
      copy_amatrix(false, a, z);

    if (alpha != 1.0)
      scale_amatrix(alpha, z);

    addmul_amatrix(1.0, false, &b->A, true, &b->B, z);

    /* Find singular value decomposition of Z */
    k = UINT_MIN(rows, cols);
    u = init_amatrix(&tmp2, rows, k);
    vt = init_amatrix(&tmp3, k, cols);
    sigma = init_realavector(&tmp5, k);
    svd_amatrix(z, sigma, u, vt);

    /* Determine rank */
    knew = findrank_truncmode(tm, eps, sigma);

    /* Set new rank */
    setrank_rkmatrix(b, knew);

    /* Copy singular vectors */
    clear_amatrix(&b->A);
    u1 = init_sub_amatrix(&tmp4, u, rows, 0, knew, 0);
    diageval_realavector_amatrix(1.0, true, sigma, true, u1);
    copy_amatrix(false, u1, &b->A);
    uninit_amatrix(u1);

    clear_amatrix(&b->B);
    vt1 = init_sub_amatrix(&tmp4, vt, knew, 0, cols, 0);
    copy_amatrix(true, vt1, &b->B);
    uninit_amatrix(vt1);
  }

  /* Clean up */
  uninit_realavector(sigma);
  uninit_amatrix(vt);
  uninit_amatrix(u);
  uninit_amatrix(z);
}

/* ------------------------------------------------------------
 * Add an rkmatrix to an amatrix
 * ------------------------------------------------------------ */

void
add_rkmatrix_amatrix(field alpha, bool atrans, pcrkmatrix a, pamatrix b)
{
#ifdef HARITH_RKMATRIX_QUICK_EXIT
  if (a->k == 0)
    return;
#endif

  if (atrans) {
    assert(a->B.rows == b->rows);
    assert(a->A.rows == b->cols);

    if (a->k > 0)
      addmul_amatrix(alpha, false, &a->B, true, &a->A, b);
  }
  else {
    assert(a->A.rows == b->rows);
    assert(a->B.rows == b->cols);

    if (a->k > 0)
      addmul_amatrix(alpha, false, &a->A, true, &a->B, b);
  }
}

/* ------------------------------------------------------------
 * Add an hmatrix to an amatrix
 * ------------------------------------------------------------ */

void
add_hmatrix_amatrix(field alpha, bool atrans, pchmatrix a, pamatrix b)
{
  amatrix   tmp;
  pamatrix  b1;
  uint      rsons, csons;
  uint      i, j;
  uint      roff, coff;

  if (atrans) {
    assert(a->rc->size == b->cols);
    assert(a->cc->size == b->rows);
  }
  else {
    assert(a->rc->size == b->rows);
    assert(a->cc->size == b->cols);
  }

  if (a->f)
    add_amatrix(alpha, atrans, a->f, b);
  else if (a->r)
    add_rkmatrix_amatrix(alpha, atrans, a->r, b);
  else {
    rsons = a->rsons;
    csons = a->csons;

    coff = 0;
    for (j = 0; j < csons; j++) {
      roff = 0;
      for (i = 0; i < rsons; i++) {
	b1 = (atrans ?
	      init_sub_amatrix(&tmp, b, a->son[i + j * rsons]->cc->size, coff,
			       a->son[i + j * rsons]->rc->size, roff) :
	      init_sub_amatrix(&tmp, b, a->son[i + j * rsons]->rc->size, roff,
			       a->son[i + j * rsons]->cc->size, coff));

	add_hmatrix_amatrix(alpha, atrans, a->son[i + j * rsons], b1);

	uninit_amatrix(b1);

	roff += a->son[i]->rc->size;
      }
      assert(roff == a->rc->size);

      coff += a->son[j * rsons]->cc->size;
    }
    assert(coff == a->cc->size);
  }
}

/* ------------------------------------------------------------
 * Truncated addition of an rkmatrix to another rkmatrix.
 * ------------------------------------------------------------ */

/* First version: Set up A and B and apply SVD directly to A B^*.
 * This approach is only advisable if the rank is high compared to
 * the number of rows and columns. */
static void
add_ab_rkmatrix(field alpha, pcrkmatrix src, pctruncmode tm,
		real eps, prkmatrix trg)
{
  amatrix   tmp1, tmp2, tmp3, tmp4, tmp5;
  realavector tmp6;
  pamatrix  a, b, c, a1, b1, u, vt;
  prealavector sigma;
  uint      rows, cols;
  uint      k, k1, knew;

  rows = trg->A.rows;
  cols = trg->B.rows;
  k = trg->k + src->k;

  assert(src->A.rows == rows);
  assert(src->B.rows == cols);

  /* Create matrices A = (alpha Asrc, Atrg) and B = (Bsrc, Btrg) */
  a = init_amatrix(&tmp1, rows, k);
  b = init_amatrix(&tmp2, cols, k);

  a1 = init_sub_amatrix(&tmp3, a, rows, 0, src->k, 0);
  copy_amatrix(false, &src->A, a1);
  if (alpha != 1.0) {
    scale_amatrix(alpha, a1);
  }
  uninit_amatrix(a1);

  b1 = init_sub_amatrix(&tmp3, b, cols, 0, src->k, 0);
  copy_amatrix(false, &src->B, b1);
  uninit_amatrix(b1);

  a1 = init_sub_amatrix(&tmp3, a, rows, 0, trg->k, src->k);
  copy_amatrix(false, &trg->A, a1);
  uninit_amatrix(a1);

  b1 = init_sub_amatrix(&tmp3, b, cols, 0, trg->k, src->k);
  copy_amatrix(false, &trg->B, b1);
  uninit_amatrix(b1);

  /* Compute C = A B^* */
  c = init_amatrix(&tmp3, rows, cols);
  clear_amatrix(c);
  addmul_amatrix(1.0, false, a, true, b, c);
  k1 = UINT_MIN(rows, cols);

  /* Compute singular value decomposition */
  u = init_amatrix(&tmp4, rows, k1);
  vt = init_amatrix(&tmp5, k1, cols);
  sigma = init_realavector(&tmp6, k1);
  svd_amatrix(c, sigma, u, vt);

  /* Determine rank */
  knew = findrank_truncmode(tm, eps, sigma);

  /* Set new rank */
  setrank_rkmatrix(trg, knew);

  /* Copy singular vectors */
  copy_sub_amatrix(false, u, &trg->A);
  copy_sub_amatrix(true, vt, &trg->B);

  /* Scale singular vectors */
  diageval_realavector_amatrix(1.0, true, sigma, true, &trg->A);

  /* Clean up */
  uninit_realavector(sigma);
  uninit_amatrix(vt);
  uninit_amatrix(u);
  uninit_amatrix(c);
  uninit_amatrix(b);
  uninit_amatrix(a);
}

/* Second version: Turn A into an upper triangular matrix by a QR
 * factorization, then apply SVD to the resulting reduced matrix.
 * Advisable if the number of rows is significantly larger than the
 * rank. */
static void
add_qb_rkmatrix(field alpha, pcrkmatrix src, pctruncmode tm,
		real eps, prkmatrix trg)
{
  amatrix   tmp1, tmp2, tmp3, tmp4, tmp5;
  avector   tmp6;
  realavector tmp7;
  pamatrix  a, b, c, a1, b1, u, vt;
  pavector  tau;
  prealavector sigma;
  uint      rows, cols;
  uint      k, kr, k1, knew;

  rows = trg->A.rows;
  cols = trg->B.rows;
  k = trg->k + src->k;

  assert(src->A.rows == rows);
  assert(src->B.rows == cols);

  /* Create matrices A = (alpha Asrc, Atrg) and B = (Bsrc, Btrg) */
  a = init_amatrix(&tmp1, rows, k);
  b = init_amatrix(&tmp2, cols, k);

  a1 = init_sub_amatrix(&tmp3, a, rows, 0, src->k, 0);
  copy_amatrix(false, &src->A, a1);
  if (alpha != 1.0)
    scale_amatrix(alpha, a1);
  uninit_amatrix(a1);

  b1 = init_sub_amatrix(&tmp3, b, cols, 0, src->k, 0);
  copy_amatrix(false, &src->B, b1);
  uninit_amatrix(b1);

  a1 = init_sub_amatrix(&tmp3, a, rows, 0, trg->k, src->k);
  copy_amatrix(false, &trg->A, a1);
  uninit_amatrix(a1);

  b1 = init_sub_amatrix(&tmp3, b, cols, 0, trg->k, src->k);
  copy_amatrix(false, &trg->B, b1);
  uninit_amatrix(b1);

  /* Compute QR factorization of A */
  tau = init_avector(&tmp6, k);
  qrdecomp_amatrix(a, tau);

  /* Overwrite B by C = B A^* (C^* = A B^*) */
  kr = UINT_MIN(rows, k);
  a1 = init_sub_amatrix(&tmp3, a, kr, 0, k, 0);
  triangulareval_amatrix(false, false, false, a1, true, b);
  uninit_amatrix(a1);
  c = init_sub_amatrix(&tmp5, b, cols, 0, kr, 0);

  /* Compute singular value decomposition */
  k1 = UINT_MIN(cols, kr);
  u = init_amatrix(&tmp3, cols, k1);
  vt = init_amatrix(&tmp4, k1, kr);
  sigma = init_realavector(&tmp7, k1);
  svd_amatrix(c, sigma, u, vt);

  /* Determine rank */
  knew = findrank_truncmode(tm, eps, sigma);

  /* Set new rank */
  setrank_rkmatrix(trg, knew);

  /* Copy singular vectors */
  copy_sub_amatrix(false, u, &trg->B);
  clear_amatrix(&trg->A);
  copy_sub_amatrix(true, vt, &trg->A);
  qreval_amatrix(false, a, tau, &trg->A);
  uninit_avector(tau);

  /* Scale singular vectors */
  diageval_realavector_amatrix(1.0, true, sigma, true, &trg->B);

  /* Clean up */
  uninit_realavector(sigma);
  uninit_amatrix(vt);
  uninit_amatrix(u);
  uninit_amatrix(c);
  uninit_amatrix(b);
  uninit_amatrix(a);
}

/* Third version: Turn B into an upper triangular matrix by a QR
 * factorization, then apply SVD to the resulting reduced matrix.
 * Advisable if the number of columns is significantly larger than the
 * rank. */
static void
add_aq_rkmatrix(field alpha, pcrkmatrix src, pctruncmode tm,
		real eps, prkmatrix trg)
{
  amatrix   tmp1, tmp2, tmp3, tmp4, tmp5;
  avector   tmp6;
  realavector tmp7;
  pamatrix  a, b, c, a1, b1, u, vt;
  pavector  tau;
  prealavector sigma;
  uint      rows, cols;
  uint      k, kc, k1, knew;

  rows = trg->A.rows;
  cols = trg->B.rows;
  k = trg->k + src->k;

  assert(src->A.rows == rows);
  assert(src->B.rows == cols);

  /* Create matrices A = (alpha Asrc, Atrg) and B = (Bsrc, Btrg) */
  a = init_amatrix(&tmp1, rows, k);
  b = init_amatrix(&tmp2, cols, k);

  a1 = init_sub_amatrix(&tmp3, a, rows, 0, src->k, 0);
  copy_amatrix(false, &src->A, a1);
  if (alpha != 1.0)
    scale_amatrix(alpha, a1);
  uninit_amatrix(a1);

  b1 = init_sub_amatrix(&tmp3, b, cols, 0, src->k, 0);
  copy_amatrix(false, &src->B, b1);
  uninit_amatrix(b1);

  a1 = init_sub_amatrix(&tmp3, a, rows, 0, trg->k, src->k);
  copy_amatrix(false, &trg->A, a1);
  uninit_amatrix(a1);

  b1 = init_sub_amatrix(&tmp3, b, cols, 0, trg->k, src->k);
  copy_amatrix(false, &trg->B, b1);
  uninit_amatrix(b1);

  /* Compute QR factorization of B */
  tau = init_avector(&tmp6, k);
  qrdecomp_amatrix(b, tau);

  /* Overwrite A by C = A B^* (C^* = B A^*) */
  kc = UINT_MIN(b->rows, k);
  b1 = init_sub_amatrix(&tmp3, b, kc, 0, k, 0);
  triangulareval_amatrix(false, false, false, b1, true, a);
  uninit_amatrix(b1);
  c = init_sub_amatrix(&tmp5, a, rows, 0, kc, 0);

  /* Compute singular value decomposition */
  k1 = UINT_MIN(rows, kc);
  u = init_amatrix(&tmp3, rows, k1);
  vt = init_amatrix(&tmp4, k1, kc);
  sigma = init_realavector(&tmp7, k1);
  svd_amatrix(c, sigma, u, vt);

  /* Determine rank */
  knew = findrank_truncmode(tm, eps, sigma);

  /* Set new rank */
  setrank_rkmatrix(trg, knew);

  /* Copy singular vectors */
  copy_sub_amatrix(false, u, &trg->A);
  clear_amatrix(&trg->B);
  copy_sub_amatrix(true, vt, &trg->B);
  qreval_amatrix(false, b, tau, &trg->B);
  uninit_avector(tau);

  /* Scale singular vectors */
  diageval_realavector_amatrix(1.0, true, sigma, true, &trg->A);

  /* Clean up */
  uninit_realavector(sigma);
  uninit_amatrix(vt);
  uninit_amatrix(u);
  uninit_amatrix(c);
  uninit_amatrix(b);
  uninit_amatrix(a);
}

/* Fourth version: Reduce both A and B to upper triangular matrices
 * by QR factorizations, then compute SVD of the reduced matrix.
 * Advisable if both the number of rows and the number of columns are
 * larger than the rank. */
static void
add_qq_rkmatrix(field alpha, pcrkmatrix src, pctruncmode tm,
		real eps, prkmatrix trg)
{
  amatrix   tmp1, tmp2, tmp3, tmp4, tmp5;
  avector   tmp6, tmp7;
  realavector tmp8;
  pamatrix  a, b, c, a1, b1, u, vt;
  pavector  atau, btau;
  prealavector sigma;
  uint      rows, cols;
  uint      k, ak, bk, k1, knew;

  rows = trg->A.rows;
  cols = trg->B.rows;
  k = trg->k + src->k;

  assert(src->A.rows == rows);
  assert(src->B.rows == cols);

  /* Create matrices A = (alpha Asrc, Atrg) and B = (Bsrc, Btrg) */
  a = init_amatrix(&tmp1, rows, k);
  b = init_amatrix(&tmp2, cols, k);

  a1 = init_sub_amatrix(&tmp3, a, rows, 0, src->k, 0);
  copy_amatrix(false, &src->A, a1);
  if (alpha != 1.0)
    scale_amatrix(alpha, a1);
  uninit_amatrix(a1);

  b1 = init_sub_amatrix(&tmp3, b, cols, 0, src->k, 0);
  copy_amatrix(false, &src->B, b1);
  uninit_amatrix(b1);

  a1 = init_sub_amatrix(&tmp3, a, rows, 0, trg->k, src->k);
  copy_amatrix(false, &trg->A, a1);
  uninit_amatrix(a1);

  b1 = init_sub_amatrix(&tmp3, b, cols, 0, trg->k, src->k);
  copy_amatrix(false, &trg->B, b1);
  uninit_amatrix(b1);

  /* Compute QR factorization Q_A R_A = A */
  atau = init_avector(&tmp6, k);
  qrdecomp_amatrix(a, atau);
  ak = UINT_MIN(k, a->rows);

  /* Compute QR factorization Q_B R_B = B */
  btau = init_avector(&tmp7, k);
  qrdecomp_amatrix(b, btau);
  bk = UINT_MIN(k, b->rows);

  /* Compute condensed matrix C = R_A R_B^* */
  c = init_amatrix(&tmp3, ak, bk);
  clear_amatrix(c);
  a1 = init_sub_amatrix(&tmp4, a, ak, 0, k, 0);
  b1 = init_sub_amatrix(&tmp5, b, bk, 0, k, 0);
  triangularaddmul_amatrix(1.0, false, false, a1, false, true, b1, c);
  uninit_amatrix(b1);
  uninit_amatrix(a1);

  /* Find singular value decomposition of C */
  k1 = UINT_MIN(ak, bk);
  u = init_amatrix(&tmp4, ak, k1);
  vt = init_amatrix(&tmp5, k1, bk);
  sigma = init_realavector(&tmp8, k1);
  svd_amatrix(c, sigma, u, vt);

  /* Determine rank */
  knew = findrank_truncmode(tm, eps, sigma);

  /* Set new rank */
  setrank_rkmatrix(trg, knew);

  /* Scale singular vectors */
  diageval_realavector_amatrix(1.0, true, sigma, true, u);

  /* Copy singular vectors */
  clear_amatrix(&trg->A);
  copy_sub_amatrix(false, u, &trg->A);
  clear_amatrix(&trg->B);
  copy_sub_amatrix(true, vt, &trg->B);

  /* Multiply singular vectors by Q_A and Q_B, respectively */
  qreval_amatrix(false, a, atau, &trg->A);
  qreval_amatrix(false, b, btau, &trg->B);

  /* Clean up */
  uninit_realavector(sigma);
  uninit_amatrix(vt);
  uninit_amatrix(u);
  uninit_amatrix(c);
  uninit_avector(btau);
  uninit_amatrix(b);
  uninit_avector(atau);
  uninit_amatrix(a);
}

/* User-visible function, chooses appropriate truncation function by
 * considering the rank, number of rows and number of columns. */
void
add_rkmatrix(field alpha, pcrkmatrix src, pctruncmode tm, real eps,
	     prkmatrix trg)
{
  uint      rows, cols, k;

#ifdef HARITH_RKMATRIX_QUICK_EXIT
  if (src->k == 0)
    return;
#endif

  rows = trg->A.rows;
  cols = trg->B.rows;
  k = src->k + trg->k;

  /* Choose most efficient truncation algorithm */
  if (k < rows) {
    if (k < cols) {
      /* rows and cols large, use QR decomposition for A and B */
      add_qq_rkmatrix(alpha, src, tm, eps, trg);
    }
    else {
      /* rows large, cols small, use QR decomposition for A only */
      add_qb_rkmatrix(alpha, src, tm, eps, trg);
    }
  }
  else {
    if (k < cols) {
      /* rows small, cols large, use QR decomposition for B only */
      add_aq_rkmatrix(alpha, src, tm, eps, trg);
    }
    else {
      /* rows and cols small, no QR decomposition required */
      add_ab_rkmatrix(alpha, src, tm, eps, trg);
    }
  }
}

/* ------------------------------------------------------------
 * Add an amatrix to an hmatrix.
 * ------------------------------------------------------------ */

void
add_amatrix_hmatrix(field alpha, bool atrans, pcamatrix a, pctruncmode tm,
		    real eps, phmatrix b)
{
  amatrix   tmp;
  phmatrix  b1;
  pamatrix  a1;
  uint      rsons, csons;
  uint      roff, coff;
  uint      i, j;

#ifdef HARITH_AMATRIX_QUICK_EXIT
  if (normfrob_amatrix(a) == 0.0)
    return;
#endif

  if (b->r)
    add_amatrix_rkmatrix(alpha, atrans, a, tm, eps, b->r);
  else if (b->f)
    add_amatrix(alpha, atrans, a, b->f);
  else {
    rsons = b->rsons;
    csons = b->csons;

    coff = 0;
    for (j = 0; j < csons; j++) {
      roff = 0;
      for (i = 0; i < rsons; i++) {
	b1 = b->son[i + j * rsons];
	a1 = init_sub_amatrix(&tmp, (pamatrix) a, b1->rc->size, roff,
			      b1->cc->size, coff);

	add_amatrix_hmatrix(alpha, atrans, a1, tm, eps, b1);

	uninit_amatrix(a1);

	roff += b1->rc->size;
      }
      assert(roff == b->rc->size);

      coff += b->son[j * rsons]->cc->size;
    }
    assert(coff == b->cc->size);
  }
}

void
add_lower_amatrix_hmatrix(field alpha, bool atrans, pcamatrix a,
			  pctruncmode tm, real eps, phmatrix b)
{
  amatrix   tmp;
  pamatrix  a1;
  phmatrix  b1;
  uint      sons;
  uint      roff, coff;
  uint      i, j;

  assert(b->rc == b->cc);

#ifdef HARITH_AMATRIX_QUICK_EXIT
  if (normfrob_amatrix(a) == 0.0)
    return;
#endif

  if (b->f)
    add_amatrix(alpha, atrans, a, b->f);
  else {
    assert(b->son);
    assert(b->rsons == b->csons);

    sons = b->rsons;

    coff = 0;
    for (j = 0; j < sons; j++) {
      /* Strictly upper triangular part */
      roff = 0;
      for (i = 0; i < j; i++) {
	b1 = b->son[i + j * sons];

	roff += b1->rc->size;
      }

      /* Diagonal */
      b1 = b->son[i + j * sons];
      a1 = (atrans ?
	    init_sub_amatrix(&tmp, (pamatrix) a, b1->cc->size, coff,
			     b1->rc->size, roff) :
	    init_sub_amatrix(&tmp, (pamatrix) a, b1->rc->size, roff,
			     b1->cc->size, coff));

      add_lower_amatrix_hmatrix(alpha, atrans, a1, tm, eps, b1);

      uninit_amatrix(a1);

      roff += b1->rc->size;
      i++;

      /* Strictly lower triangular part */
      for (; i < sons; i++) {
	b1 = b->son[i + j * sons];
	a1 = (atrans ?
	      init_sub_amatrix(&tmp, (pamatrix) a, b1->cc->size, coff,
			       b1->rc->size, roff) :
	      init_sub_amatrix(&tmp, (pamatrix) a, b1->rc->size, roff,
			       b1->cc->size, coff));

	add_amatrix_hmatrix(alpha, atrans, a1, tm, eps, b1);

	uninit_amatrix(a1);

	roff += b1->rc->size;
      }
      assert(roff == b->rc->size);

      coff += b->son[j * sons]->cc->size;
    }
    assert(coff == b->cc->size);
  }
}

/* ------------------------------------------------------------
 * Add an rkmatrix to an hmatrix.
 * ------------------------------------------------------------ */

void
add_rkmatrix_hmatrix(field alpha, pcrkmatrix r, pctruncmode tm, real eps,
		     phmatrix a)
{
  rkmatrix  tmp;
  phmatrix  a1;
  pcrkmatrix r1;
  uint      rsons, csons;
  uint      roff, coff;
  uint      i, j;

#ifdef HARITH_RKMATRIX_QUICK_EXIT
  if (r->k == 0)
    return;
#endif

  if (a->r)
    add_rkmatrix(alpha, r, tm, eps, a->r);
  else if (a->f)
    addmul_amatrix(alpha, false, &r->A, true, &r->B, a->f);
  else {
    rsons = a->rsons;
    csons = a->csons;

    coff = 0;
    for (j = 0; j < csons; j++) {
      roff = 0;
      for (i = 0; i < rsons; i++) {
	a1 = a->son[i + j * rsons];
	r1 =
	  init_sub_rkmatrix(&tmp, r, a1->rc->size, roff, a1->cc->size, coff);

	add_rkmatrix_hmatrix(alpha, r1, tm, eps, a1);

	uninit_rkmatrix((prkmatrix) r1);

	roff += a1->rc->size;
      }
      assert(roff == a->rc->size);

      coff += a->son[j * rsons]->cc->size;
    }
    assert(coff == a->cc->size);
  }
}

void
add_lower_rkmatrix_hmatrix(field alpha, pcrkmatrix a, pctruncmode tm,
			   real eps, phmatrix b)
{
  rkmatrix  tmp;
  pcrkmatrix a1;
  phmatrix  b1;
  uint      sons;
  uint      roff, coff;
  uint      i, j;

  assert(b->rc == b->cc);

#ifdef HARITH_RKMATRIX_QUICK_EXIT
  if (a->k == 0)
    return;
#endif

  if (b->f)
    add_rkmatrix_amatrix(alpha, false, a, b->f);
  else {
    assert(b->son);
    assert(b->rsons == b->csons);

    sons = b->rsons;

    coff = 0;
    for (j = 0; j < sons; j++) {
      /* Strictly upper triangular part */
      roff = 0;
      for (i = 0; i < j; i++) {
	b1 = b->son[i + j * sons];

	roff += b1->rc->size;
      }

      /* Diagonal */
      b1 = b->son[i + j * sons];
      a1 = init_sub_rkmatrix(&tmp, a, b1->rc->size, roff, b1->cc->size, coff);

      add_lower_rkmatrix_hmatrix(alpha, a1, tm, eps, b1);

      uninit_rkmatrix((prkmatrix) a1);

      roff += b1->rc->size;
      i++;

      /* Strictly lower triangular part */
      for (; i < sons; i++) {
	b1 = b->son[i + j * sons];
	a1 =
	  init_sub_rkmatrix(&tmp, a, b1->rc->size, roff, b1->cc->size, coff);

	add_rkmatrix_hmatrix(alpha, a1, tm, eps, b1);

	uninit_rkmatrix((prkmatrix) a1);

	roff += b1->rc->size;
      }
      assert(roff == b->rc->size);

      coff += b->son[j * sons]->cc->size;
    }
    assert(coff == b->cc->size);
  }
}

void
add_hmatrix(field alpha, pchmatrix a, pctruncmode tm, real eps, phmatrix b)
{
  phmatrix  bs;
  prkmatrix rk;
  uint      rsons, csons;
  uint      i, j;

  assert(a->rc == b->rc);
  assert(a->cc == b->cc);

  if (a->r) {
    add_rkmatrix_hmatrix(alpha, a->r, tm, eps, b);
  }
  else if (a->f) {
    add_amatrix_hmatrix(alpha, false, a->f, tm, eps, b);
  }
  else {
    if (b->son != NULL) {
      assert(a->rsons == b->rsons);	/* Should be generalized */
      assert(a->csons == b->csons);

      rsons = a->rsons;
      csons = a->csons;

      for (j = 0; j < csons; j++) {
	for (i = 0; i < rsons; i++) {
	  add_hmatrix(alpha, a->son[i + j * rsons], tm, eps,
		      b->son[i + j * rsons]);
	}
      }
    }
    else if (b->r != NULL) {
      bs = split_rkmatrix(b->r, b->rc, b->cc, a->rsons != b->rsons,
			  a->csons != b->csons, false);

      rsons = a->rsons;
      csons = a->csons;

      for (j = 0; j < csons; j++) {
	for (i = 0; i < rsons; i++) {
	  add_hmatrix(alpha, a->son[i + j * rsons], tm, eps,
		      bs->son[i + j * rsons]);
	}
      }

      rk = merge_hmatrix_rkmatrix(bs, tm, eps);
      add_rkmatrix(1.0, rk, tm, eps, b->r);

      del_rkmatrix(rk);
      del_hmatrix(bs);
    }
    else if (b->f != NULL) {
      bs = split_sub_amatrix(b->f, b->rc, b->cc, a->rsons != b->rsons,
			     a->csons != b->csons);

      rsons = a->rsons;
      csons = a->csons;

      for (j = 0; j < csons; j++) {
	for (i = 0; i < rsons; i++) {
	  add_hmatrix(alpha, a->son[i + j * rsons], tm, eps,
		      bs->son[i + j * rsons]);
	}
      }

      del_hmatrix(bs);
    }
  }
}

/* ------------------------------------------------------------
 * Split a leaf hmatrix
 * ------------------------------------------------------------ */

phmatrix
split_sub_amatrix(pcamatrix f, pccluster rc, pccluster cc, bool rsplit,
		  bool csplit)
{
  phmatrix  s, s1;
  uint      rsons, csons;
  uint      roff, coff;
  uint      i, j;

  s = 0;

  if (rsplit && rc->sons > 0) {
    rsons = rc->sons;

    if (csplit && cc->sons > 0) {
      /* Split full matrix in row and column direction */
      csons = cc->sons;

      s = new_super_hmatrix(rc, cc, rsons, csons);

      coff = 0;
      for (j = 0; j < csons; j++) {
	roff = 0;
	for (i = 0; i < rsons; i++) {
	  s1 = new_hmatrix(rc->son[i], cc->son[j]);
	  s1->f = new_sub_amatrix((pamatrix) f, rc->son[i]->size, roff,
				  cc->son[j]->size, coff);
	  s1->desc = 1;

	  ref_hmatrix(s->son + i + j * rsons, s1);

	  roff += rc->son[i]->size;
	}
	assert(roff == rc->size);

	coff += cc->son[j]->size;
      }
      assert(coff == cc->size);
    }
    else {
      /* Split full matrix in row direction */
      s = new_super_hmatrix(rc, cc, rsons, 1);

      roff = 0;
      for (i = 0; i < rsons; i++) {
	s1 = new_hmatrix(rc->son[i], cc);
	s1->f =
	  new_sub_amatrix((pamatrix) f, rc->son[i]->size, roff, cc->size, 0);
	s1->desc = 1;

	ref_hmatrix(s->son + i, s1);

	roff += rc->son[i]->size;
      }
      assert(roff == rc->size);
    }
  }
  else {
    if (csplit && cc->sons > 0) {
      /* Split full matrix in column direction */
      csons = cc->sons;

      s = new_super_hmatrix(rc, cc, 1, csons);

      coff = 0;
      for (j = 0; j < csons; j++) {
	s1 = new_hmatrix(rc, cc->son[j]);
	s1->f = new_sub_amatrix((pamatrix) f, rc->size, 0, cc->son[j]->size,
				coff);
	s1->desc = 1;

	ref_hmatrix(s->son + j, s1);

	coff += cc->son[j]->size;
      }
      assert(coff == cc->size);
    }
    else {
      /* Do not split full matrix at all */
      s = new_super_hmatrix(rc, cc, 1, 1);

      s1 = new_hmatrix(rc, cc);
      s1->f = new_sub_amatrix((pamatrix) f, rc->size, 0, cc->size, 0);
      s1->desc = 1;

      ref_hmatrix(s->son, s1);
    }
  }

  update_hmatrix(s);

  return s;
}

phmatrix
split_rkmatrix(pcrkmatrix r, pccluster rc, pccluster cc, bool rsplit,
	       bool csplit, bool copy)
{
  phmatrix  s, s1;
  prkmatrix r1;
  rkmatrix  tmp;
  uint      rsons, csons;
  uint      roff, coff;
  uint      i, j;

  s = 0;

  if (rsplit && rc->sons > 0) {
    rsons = rc->sons;

    if (csplit && cc->sons > 0) {
      /* Split rk matrix in row and column direction */
      csons = cc->sons;

      s = new_super_hmatrix(rc, cc, rsons, csons);

      coff = 0;
      for (j = 0; j < csons; j++) {
	roff = 0;
	for (i = 0; i < rsons; i++) {
	  s1 = 0;
	  if (copy) {
	    s1 = new_rk_hmatrix(rc->son[i], cc->son[j], r->k);

	    r1 =
	      (prkmatrix) init_sub_rkmatrix(&tmp, r, rc->son[i]->size, roff,
					    cc->son[j]->size, coff);
	    copy_rkmatrix(false, r1, s1->r);
	    uninit_rkmatrix(r1);
	  }
	  else
	    s1 = new_rk_hmatrix(rc->son[i], cc->son[j], 0);

	  ref_hmatrix(s->son + i + j * rsons, s1);

	  roff += rc->son[i]->size;
	}
	assert(roff == rc->size);

	coff += cc->son[j]->size;
      }
      assert(coff == cc->size);
    }
    else {
      /* Split rk matrix in row direction */
      s = new_super_hmatrix(rc, cc, rsons, 1);

      roff = 0;
      for (i = 0; i < rsons; i++) {
	s1 = 0;
	if (copy) {
	  s1 = new_rk_hmatrix(rc->son[i], cc, r->k);

	  r1 = (prkmatrix) init_sub_rkmatrix(&tmp, r, rc->son[i]->size, roff,
					     cc->size, 0);
	  copy_rkmatrix(false, r1, s1->r);
	  uninit_rkmatrix(r1);
	}
	else
	  s1 = new_rk_hmatrix(rc->son[i], cc, 0);

	ref_hmatrix(s->son + i, s1);

	roff += rc->son[i]->size;
      }
      assert(roff == rc->size);
    }
  }
  else {
    if (csplit && cc->sons > 0) {
      /* Split rk matrix in column direction */
      csons = cc->sons;

      s = new_super_hmatrix(rc, cc, 1, csons);

      coff = 0;
      for (j = 0; j < csons; j++) {
	s1 = 0;
	if (copy) {
	  s1 = new_rk_hmatrix(rc, cc->son[j], r->k);

	  r1 = (prkmatrix) init_sub_rkmatrix(&tmp, r, rc->size, 0,
					     cc->son[j]->size, coff);
	  copy_rkmatrix(false, r1, s1->r);
	  uninit_rkmatrix(r1);
	}
	else
	  s1 = new_rk_hmatrix(rc, cc->son[j], 0);

	ref_hmatrix(s->son + j, s1);

	coff += cc->son[j]->size;
      }
      assert(coff == cc->size);
    }
    else {
      /* Do not split rk matrix at all */
      s = new_super_hmatrix(rc, cc, 1, 1);

      if (copy) {
	s1 = new_rk_hmatrix(rc, cc, r->k);

	copy_rkmatrix(false, r, s1->r);
      }
      else
	s1 = new_rk_hmatrix(rc, cc, 0);

      ref_hmatrix(s->son, s1);
    }
  }

  update_hmatrix(s);

  return s;
}

phmatrix
split_sub_rkmatrix(pcrkmatrix r, pccluster rc, pccluster cc,
		   bool rsplit, bool csplit)
{
  phmatrix  s, s1;
  uint      rsons, csons;
  uint      roff, coff;
  uint      i, j;

  s = 0;

  if (rsplit && rc->sons > 0) {
    rsons = rc->sons;

    if (csplit && cc->sons > 0) {
      /* Split rk matrix in row and column direction */
      csons = cc->sons;

      s = new_super_hmatrix(rc, cc, rsons, csons);

      coff = 0;
      for (j = 0; j < csons; j++) {
	roff = 0;
	for (i = 0; i < rsons; i++) {
	  s1 = new_hmatrix(rc->son[i], cc->son[j]);

	  s1->r = (prkmatrix) new_sub_rkmatrix(r, rc->son[i]->size, roff,
					       cc->son[j]->size, coff);
	  s1->desc = 1;

	  ref_hmatrix(s->son + i + j * rsons, s1);

	  roff += rc->son[i]->size;
	}
	assert(roff == rc->size);

	coff += cc->son[j]->size;
      }
      assert(coff == cc->size);
    }
    else {
      /* Split rk matrix in row direction */
      s = new_super_hmatrix(rc, cc, rsons, 1);

      roff = 0;
      for (i = 0; i < rsons; i++) {
	s1 = new_hmatrix(rc->son[i], cc);

	s1->r = (prkmatrix) new_sub_rkmatrix(r, rc->son[i]->size, roff,
					     cc->size, 0);
	s1->desc = 1;

	ref_hmatrix(s->son + i, s1);

	roff += rc->son[i]->size;
      }
      assert(roff == rc->size);
    }
  }
  else {
    if (csplit && cc->sons > 0) {
      /* Split rk matrix in column direction */
      csons = cc->sons;

      s = new_super_hmatrix(rc, cc, 1, csons);

      coff = 0;
      for (j = 0; j < csons; j++) {
	s1 = new_hmatrix(rc, cc->son[j]);

	s1->r = (prkmatrix) new_sub_rkmatrix(r, rc->size, 0, cc->son[j]->size,
					     coff);
	s1->desc = 1;

	ref_hmatrix(s->son + j, s1);

	coff += cc->son[j]->size;
      }
      assert(coff == cc->size);
    }
    else {
      /* Do not split rk matrix at all */
      s = new_super_hmatrix(rc, cc, 1, 1);

      s1 = new_hmatrix(rc, cc);

      s1->r = (prkmatrix) new_sub_rkmatrix(r, rc->size, 0, cc->size, 0);
      s1->desc = 1;

      ref_hmatrix(s->son, s1);
    }
  }

  update_hmatrix(s);

  return s;
}

/* ------------------------------------------------------------
 * Merge two rk matrices
 * ------------------------------------------------------------ */

static void
merge_aq_rkmatrix(bool colmerge, pcrkmatrix src, pctruncmode tm,
		  real eps, prkmatrix trg)
{
  amatrix   tmp1, tmp2, tmp3, tmp4, tmp5, tmp6;
  avector   tmp7, tmp8;
  realavector tmp9;
  pamatrix  a, a1, a2, a1r, a2r, b, b1, b2, b1r, b2r;
  pamatrix  c, c1, u, u1, vt, v1, vt1;
  pavector  tau1, tau2;
  prealavector sigma;
  uint      rows, cols, rows1, rows2, cols1, cols2;
  uint      k, k1, k2, kk, kmax, knew;

  if (colmerge) {
    assert(src->B.rows == trg->B.rows);

    rows1 = trg->A.rows;
    rows2 = src->A.rows;
    rows = rows1 + rows2;
    cols = trg->B.rows;
    k = trg->k + src->k;

    /* Set up matrix B = (B1 B2) */
    b = init_amatrix(&tmp1, cols, k);

    b1 = init_sub_amatrix(&tmp2, b, cols, 0, trg->k, 0);
    copy_amatrix(false, &trg->B, b1);
    uninit_amatrix(b1);

    b2 = init_sub_amatrix(&tmp2, b, cols, 0, src->k, trg->k);
    copy_amatrix(false, &src->B, b2);
    uninit_amatrix(b2);

    /* Compute factorization A1 = Q1 R1 */
    k1 = UINT_MIN(trg->A.rows, trg->k);
    a1 = init_amatrix(&tmp3, trg->A.rows, trg->k);
    copy_amatrix(false, &trg->A, a1);
    tau1 = init_avector(&tmp7, k1);
    qrdecomp_amatrix(a1, tau1);

    /* Compute factorization A2 = Q2 R2 */
    k2 = UINT_MIN(src->A.rows, src->k);
    a2 = init_amatrix(&tmp4, src->A.rows, src->k);
    copy_amatrix(false, &src->A, a2);
    tau2 = init_avector(&tmp8, k2);
    qrdecomp_amatrix(a2, tau2);

    kk = k1 + k2;
    c = init_amatrix(&tmp5, cols, kk);

    /* Compute B1 R1^* */
    b1 = init_sub_amatrix(&tmp2, b, cols, 0, trg->k, 0);
    a1r = init_sub_amatrix(&tmp6, a1, k1, 0, trg->k, 0);
    triangulareval_amatrix(false, false, false, a1r, true, b1);
    uninit_amatrix(a1r);
    uninit_amatrix(b1);
    b1 = init_sub_amatrix(&tmp2, b, cols, 0, k1, 0);
    c1 = init_sub_amatrix(&tmp6, c, cols, 0, k1, 0);
    copy_amatrix(false, b1, c1);
    uninit_amatrix(c1);
    uninit_amatrix(b1);

    /* Compute B2 R2^* */
    b2 = init_sub_amatrix(&tmp2, b, cols, 0, src->k, trg->k);
    a2r = init_sub_amatrix(&tmp6, a2, k2, 0, src->k, 0);
    triangulareval_amatrix(false, false, false, a2r, true, b2);
    uninit_amatrix(a2r);
    uninit_amatrix(b2);
    b2 = init_sub_amatrix(&tmp2, b, cols, 0, k2, trg->k);
    c1 = init_sub_amatrix(&tmp6, c, cols, 0, k2, k1);
    copy_amatrix(false, b2, c1);
    uninit_amatrix(c1);
    uninit_amatrix(b2);

    /* Intermediate clean-up */
    uninit_amatrix(b);

    /* Compute SVD */
    kmax = UINT_MIN(kk, cols);
    u = init_amatrix(&tmp1, cols, kmax);
    vt = init_amatrix(&tmp2, kmax, kk);
    sigma = init_realavector(&tmp9, kmax);
    svd_amatrix(c, sigma, u, vt);

    /* Determine rank */
    knew = findrank_truncmode(tm, eps, sigma);

    /* Resize target matrix and set new rank */
    resize_rkmatrix(trg, rows, cols, knew);

    /* Copy left singular vectors */
    u1 = init_sub_amatrix(&tmp6, u, cols, 0, knew, 0);
    copy_amatrix(false, u1, &trg->B);
    uninit_amatrix(u1);
    uninit_amatrix(u);

    /* Copy right singular vectors */
    clear_amatrix(&trg->A);

    vt1 = init_sub_amatrix(&tmp6, vt, knew, 0, k1, 0);
    v1 = init_sub_amatrix(&tmp1, &trg->A, k1, 0, knew, 0);
    copy_amatrix(true, vt1, v1);
    uninit_amatrix(v1);
    uninit_amatrix(vt1);

    vt1 = init_sub_amatrix(&tmp6, vt, knew, 0, k2, k1);
    v1 = init_sub_amatrix(&tmp1, &trg->A, k2, rows1, knew, 0);
    copy_amatrix(true, vt1, v1);
    uninit_amatrix(v1);
    uninit_amatrix(vt1);

    uninit_amatrix(vt);

    /* Apply orthogonal transformations */
    v1 = init_sub_amatrix(&tmp1, &trg->A, rows1, 0, knew, 0);
    qreval_amatrix(false, a1, tau1, v1);
    uninit_amatrix(v1);

    v1 = init_sub_amatrix(&tmp1, &trg->A, rows2, rows1, knew, 0);
    qreval_amatrix(false, a2, tau2, v1);
    uninit_amatrix(v1);

    /* Intermediate clean-up */
    uninit_avector(tau2);
    uninit_amatrix(a2);
    uninit_avector(tau1);
    uninit_amatrix(a1);

    /* Scale singular vectors */
    diageval_realavector_amatrix(1.0, true, sigma, true, &trg->B);

    /* Final clean-up */
    uninit_realavector(sigma);
    uninit_amatrix(c);
  }
  else {
    assert(src->A.rows == trg->A.rows);

    rows = src->A.rows;
    cols1 = trg->B.rows;
    cols2 = src->B.rows;
    cols = cols1 + cols2;
    k = trg->k + src->k;

    /* Set up matrix A = (A1 A2) */
    a = init_amatrix(&tmp1, rows, k);

    a1 = init_sub_amatrix(&tmp2, a, rows, 0, trg->k, 0);
    copy_amatrix(false, &trg->A, a1);
    uninit_amatrix(a1);

    a1 = init_sub_amatrix(&tmp2, a, rows, 0, src->k, trg->k);
    copy_amatrix(false, &src->A, a1);
    uninit_amatrix(a1);

    /* Compute factorization B1 = Q1 R1 */
    k1 = UINT_MIN(trg->B.rows, trg->k);
    b1 = init_amatrix(&tmp3, trg->B.rows, trg->k);
    copy_amatrix(false, &trg->B, b1);
    tau1 = init_avector(&tmp7, k1);
    qrdecomp_amatrix(b1, tau1);

    /* Compute factorization B2 = Q2 R2 */
    k2 = UINT_MIN(src->B.rows, src->k);
    b2 = init_amatrix(&tmp4, src->B.rows, src->k);
    copy_amatrix(false, &src->B, b2);
    tau2 = init_avector(&tmp8, k2);
    qrdecomp_amatrix(b2, tau2);

    kk = k1 + k2;
    c = init_amatrix(&tmp5, rows, kk);

    /* Compute A1 R1^* */
    a1 = init_sub_amatrix(&tmp2, a, rows, 0, trg->k, 0);
    b1r = init_sub_amatrix(&tmp6, b1, k1, 0, trg->k, 0);
    triangulareval_amatrix(false, false, false, b1r, true, a1);
    uninit_amatrix(b1r);
    uninit_amatrix(a1);
    a1 = init_sub_amatrix(&tmp2, a, rows, 0, k1, 0);
    c1 = init_sub_amatrix(&tmp6, c, rows, 0, k1, 0);
    copy_amatrix(false, a1, c1);
    uninit_amatrix(c1);
    uninit_amatrix(a1);

    /* Compute A2 R2^* */
    a2 = init_sub_amatrix(&tmp2, a, rows, 0, src->k, trg->k);
    b2r = init_sub_amatrix(&tmp6, b2, k2, 0, src->k, 0);
    triangulareval_amatrix(false, false, false, b2r, true, a2);
    uninit_amatrix(b2r);
    uninit_amatrix(a2);
    a2 = init_sub_amatrix(&tmp2, a, rows, 0, k2, trg->k);
    c1 = init_sub_amatrix(&tmp6, c, rows, 0, k2, k1);
    copy_amatrix(false, a2, c1);
    uninit_amatrix(c1);
    uninit_amatrix(a2);

    /* Intermediate clean-up */
    uninit_amatrix(a);

    /* Compute SVD */
    kmax = UINT_MIN(kk, rows);
    u = init_amatrix(&tmp1, rows, kmax);
    vt = init_amatrix(&tmp2, kmax, kk);
    sigma = init_realavector(&tmp9, kmax);
    svd_amatrix(c, sigma, u, vt);

    /* Determine rank */
    knew = findrank_truncmode(tm, eps, sigma);

    /* Resize target matrix and set new rank */
    resize_rkmatrix(trg, rows, cols, knew);

    /* Copy left singular vectors */
    u1 = init_sub_amatrix(&tmp6, u, rows, 0, knew, 0);
    copy_amatrix(false, u1, &trg->A);
    uninit_amatrix(u1);
    uninit_amatrix(u);

    /* Copy right singular vectors */
    clear_amatrix(&trg->B);

    vt1 = init_sub_amatrix(&tmp6, vt, knew, 0, k1, 0);
    v1 = init_sub_amatrix(&tmp1, &trg->B, k1, 0, knew, 0);
    copy_amatrix(true, vt1, v1);
    uninit_amatrix(v1);
    uninit_amatrix(vt1);

    vt1 = init_sub_amatrix(&tmp6, vt, knew, 0, k2, k1);
    v1 = init_sub_amatrix(&tmp1, &trg->B, k2, cols1, knew, 0);
    copy_amatrix(true, vt1, v1);
    uninit_amatrix(v1);
    uninit_amatrix(vt1);

    uninit_amatrix(vt);

    /* Apply orthogonal transformations */
    v1 = init_sub_amatrix(&tmp1, &trg->B, cols1, 0, knew, 0);
    qreval_amatrix(false, b1, tau1, v1);
    uninit_amatrix(v1);

    v1 = init_sub_amatrix(&tmp1, &trg->B, cols2, cols1, knew, 0);
    qreval_amatrix(false, b2, tau2, v1);
    uninit_amatrix(v1);

    /* Intermediate clean-up */
    uninit_avector(tau2);
    uninit_amatrix(b2);
    uninit_avector(tau1);
    uninit_amatrix(b1);

    /* Scale singular vectors */
    diageval_realavector_amatrix(1.0, true, sigma, true, &trg->A);

    /* Final clean-up */
    uninit_realavector(sigma);
    uninit_amatrix(c);
  }
}

void
merge_rkmatrix(bool colmerge, pcrkmatrix src, pctruncmode tm, real eps,
	       prkmatrix trg)
{
  merge_aq_rkmatrix(colmerge, src, tm, eps, trg);
}

/* ------------------------------------------------------------
 * Merge a depth-one hmatrix
 * ------------------------------------------------------------ */

prkmatrix
merge_hmatrix_rkmatrix(pchmatrix s, pctruncmode tm, real eps)
{
  rkmatrix  tmp1;
  prkmatrix r1, r;
  pcrkmatrix s1;
  uint      rsons, csons;
  uint      i, j;

  rsons = s->rsons;
  csons = s->csons;

  r = 0;
  if (rsons > 1) {
    if (csons > 1) {
      r = new_rkmatrix(s->rc->size, 0, 0);

      for (j = 0; j < csons; j++) {
	assert(s->son[j * rsons]->r != 0);
	s1 = s->son[j * rsons]->r;

	r1 = init_rkmatrix(&tmp1, s1->A.rows, s1->B.rows, s1->k);
	copy_rkmatrix(false, s1, r1);

	for (i = 1; i < rsons; i++) {
	  assert(s->son[i + j * rsons]->r != 0);
	  s1 = s->son[i + j * rsons]->r;
	  merge_rkmatrix(true, s1, tm, eps, r1);
	}

	merge_rkmatrix(false, r1, tm, eps, r);

	uninit_rkmatrix(r1);
      }
    }
    else {
      assert(s->son[0]->r != 0);
      s1 = s->son[0]->r;

      r = new_rkmatrix(s1->A.rows, s1->B.rows, s1->k);
      copy_rkmatrix(false, s1, r);

      for (i = 1; i < rsons; i++) {
	assert(s->son[i]->r != 0);
	s1 = s->son[i]->r;
	merge_rkmatrix(true, s1, tm, eps, r);
      }
    }
  }
  else {
    assert(s->son[0]->r != 0);
    s1 = s->son[0]->r;

    r = new_rkmatrix(s1->A.rows, s1->B.rows, s1->k);
    copy_rkmatrix(false, s1, r);

    for (j = 1; j < csons; j++) {
      assert(s->son[j]->r != 0);
      s1 = s->son[j]->r;
      merge_rkmatrix(false, s1, tm, eps, r);
    }
  }

  return r;
}

/* ------------------------------------------------------------
 * Multiply an rkmatrix by an amatrix.
 * ------------------------------------------------------------ */

void
addmul_rkmatrix_amatrix_amatrix(field alpha, bool xtrans, pcrkmatrix x,
				bool ytrans, pcamatrix y, bool ztrans,
				pamatrix z)
{
  pamatrix  ay, by;
  amatrix   tmp1;
  uint      rows, cols, k;

  rows = x->A.rows;
  cols = x->B.rows;
  k = x->A.cols;

  assert(x->A.cols == k);
  assert(x->B.cols == k);

  if (xtrans) {
    if (ztrans)
      assert(z->cols == cols);
    else
      assert(z->rows == cols);

    if (ytrans) {
      if (ztrans)
	assert(z->rows == y->rows);
      else
	assert(z->cols == y->rows);

      assert(y->cols == rows);
    }
    else {
      if (ztrans)
	assert(z->rows == y->cols);
      else
	assert(z->cols == y->cols);

      assert(y->rows == rows);
    }

    ay = init_amatrix(&tmp1, k, (ytrans ? y->rows : y->cols));
    clear_amatrix(ay);
    addmul_amatrix(1.0, true, &x->A, ytrans, y, ay);

    if (ztrans)
      addmul_amatrix(CONJ(alpha), true, ay, true, &x->B, z);
    else
      addmul_amatrix(alpha, false, &x->B, false, ay, z);

    uninit_amatrix(ay);
  }
  else {
    if (ztrans)
      assert(z->cols == rows);
    else
      assert(z->rows == rows);

    if (ytrans) {
      if (ztrans)
	assert(z->rows == y->rows);
      else
	assert(z->cols == y->rows);

      assert(y->cols == cols);
    }
    else {
      if (ztrans)
	assert(z->rows == y->cols);
      else
	assert(z->cols == y->cols);

      assert(y->rows == cols);
    }

    by = init_amatrix(&tmp1, k, (ytrans ? y->rows : y->cols));
    clear_amatrix(by);
    addmul_amatrix(1.0, true, &x->B, ytrans, y, by);

    if (ztrans)
      addmul_amatrix(CONJ(alpha), true, by, true, &x->A, z);
    else
      addmul_amatrix(alpha, false, &x->A, false, by, z);

    uninit_amatrix(by);
  }
}

/* ------------------------------------------------------------
 * Multiply an hmatrix by an amatrix.
 * ------------------------------------------------------------ */

/* First version: Multiply by A. */
static void
addmul_n_hmatrix_amatrix_amatrix(field alpha, pchmatrix a,
				 bool btrans, pcamatrix bp, bool ctrans,
				 pamatrix cp)
{
  pamatrix  bp1, cp1;
  amatrix   tmp1, tmp2;
  uint      rsons, csons;
  uint      boff, coff, i, j;

  assert((btrans ? bp->cols : bp->rows) == a->cc->size);
  assert((btrans ? bp->rows : bp->cols) == (ctrans ? cp->rows : cp->cols));
  assert((ctrans ? cp->cols : cp->rows) == a->rc->size);

  if (a->r) {
    addmul_rkmatrix_amatrix_amatrix(alpha, false, a->r, btrans, bp, ctrans,
				    cp);
  }
  else if (a->f) {
    if (ctrans) {
      addmul_amatrix(CONJ(alpha), !btrans, bp, true, a->f, cp);
    }
    else {
      addmul_amatrix(alpha, false, a->f, btrans, bp, cp);
    }
  }
  else {
    rsons = a->rsons;
    csons = a->csons;

    boff = 0;
    for (j = 0; j < csons; j++) {
      bp1 = (btrans ?
	     init_sub_amatrix(&tmp1, (pamatrix) bp, bp->rows, 0,
			      a->son[j * rsons]->cc->size, boff) :
	     init_sub_amatrix(&tmp1, (pamatrix) bp,
			      a->son[j * rsons]->cc->size, boff, bp->cols,
			      0));

      coff = 0;
      for (i = 0; i < rsons; i++) {
	cp1 = (ctrans ?
	       init_sub_amatrix(&tmp2, cp, cp->rows, 0, a->son[i]->rc->size,
				coff) :
	       init_sub_amatrix(&tmp2, cp, a->son[i]->rc->size, coff,
				cp->cols, 0));

	addmul_n_hmatrix_amatrix_amatrix(alpha, a->son[i + j * rsons], btrans,
					 bp1, ctrans, cp1);

	uninit_amatrix(cp1);

	coff += a->son[i]->rc->size;
      }
      assert(coff == a->rc->size);

      uninit_amatrix(bp1);

      boff += a->son[j * rsons]->cc->size;
    }
    assert(boff == a->cc->size);
  }
}

/* Second version: Multiply by the adjoint A^*. */
static void
addmul_t_hmatrix_amatrix_amatrix(field alpha, pchmatrix a,
				 bool btrans, pcamatrix bp, bool ctrans,
				 pamatrix cp)
{
  pamatrix  bp1, cp1;
  amatrix   tmp1, tmp2;
  uint      rsons, csons;
  uint      boff, coff, i, j;

  assert((btrans ? bp->cols : bp->rows) == a->rc->size);
  assert((btrans ? bp->rows : bp->cols) == (ctrans ? cp->rows : cp->cols));
  assert((ctrans ? cp->cols : cp->rows) == a->cc->size);

  if (a->r) {
    addmul_rkmatrix_amatrix_amatrix(alpha, true, a->r, btrans, bp, ctrans,
				    cp);
  }
  else if (a->f) {
    if (ctrans)
      addmul_amatrix(CONJ(alpha), !btrans, bp, false, a->f, cp);
    else
      addmul_amatrix(alpha, true, a->f, btrans, bp, cp);
  }
  else {
    rsons = a->rsons;
    csons = a->csons;

    coff = 0;
    for (j = 0; j < csons; j++) {
      cp1 = (ctrans ?
	     init_sub_amatrix(&tmp2, cp, cp->rows, 0,
			      a->son[j * rsons]->cc->size, coff) :
	     init_sub_amatrix(&tmp2, cp, a->son[j * rsons]->cc->size, coff,
			      cp->cols, 0));

      boff = 0;
      for (i = 0; i < rsons; i++) {
	bp1 = (btrans ?
	       init_sub_amatrix(&tmp1, (pamatrix) bp, bp->rows, 0,
				a->son[i]->rc->size, boff) :
	       init_sub_amatrix(&tmp1, (pamatrix) bp, a->son[i]->rc->size,
				boff, bp->cols, 0));

	addmul_t_hmatrix_amatrix_amatrix(alpha, a->son[i + j * rsons], btrans,
					 bp1, ctrans, cp1);

	uninit_amatrix(bp1);

	boff += a->son[i]->rc->size;
      }
      assert(boff == a->rc->size);

      uninit_amatrix(cp1);

      coff += a->son[j * rsons]->cc->size;
    }
    assert(coff == a->cc->size);
  }
}

/* User-visible function, checks parameters and calls appropriate
 * static function. */
void
addmul_hmatrix_amatrix_amatrix(field alpha, bool atrans, pchmatrix a,
			       bool btrans, pcamatrix bp, bool ctrans,
			       pamatrix cp)
{

  if (atrans) {
    if (btrans) {
      assert(bp->cols == a->rc->size);
      if (ctrans)
	assert(bp->rows == cp->rows);
      else
	assert(bp->rows == cp->cols);
    }
    else {
      assert(bp->rows == a->rc->size);
      if (ctrans)
	assert(bp->cols == cp->rows);
      else
	assert(bp->cols == cp->cols);
    }

    if (ctrans)
      assert(cp->cols == a->cc->size);
    else
      assert(cp->rows == a->cc->size);

    addmul_t_hmatrix_amatrix_amatrix(alpha, a, btrans, bp, ctrans, cp);
  }
  else {
    if (btrans) {
      assert(bp->cols == a->cc->size);
      if (ctrans)
	assert(bp->rows == cp->rows);
      else
	assert(bp->rows == cp->cols);
    }
    else {
      assert(bp->rows == a->cc->size);
      if (ctrans)
	assert(bp->cols == cp->rows);
      else
	assert(bp->cols == cp->cols);
    }

    if (ctrans)
      assert(cp->cols == a->rc->size);
    else
      assert(cp->rows == a->rc->size);

    addmul_n_hmatrix_amatrix_amatrix(alpha, a, btrans, bp, ctrans, cp);
  }
}

/* ------------------------------------------------------------
 * Multiply two hmatrices.
 * ------------------------------------------------------------ */

static void
addmul_nn_hmatrix(field alpha, pchmatrix x, pchmatrix y,
		  pctruncmode tm, real eps, phmatrix z)
{
  prkmatrix xy;
  pamatrix  id;
  pcamatrix xf, yf;
  pamatrix  zf, zf1;
  phmatrix  z1, ztmp;
  rkmatrix  tmp1;
  amatrix   tmp2;
  hmatrix   tmp3;
  pccluster rc, cc;
  uint      rsons, msons, csons;
  uint      roff, coff;
  uint      i, j, k;

  assert(z->rc == x->rc);
  assert(z->cc == y->cc);
  assert(x->cc == y->rc);

  rc = z->rc;
  cc = z->cc;

  if (x->f) {
    xf = x->f;
    if (z->f)			/* Z = Z + X Y   <=>   Z^* = Z^* + Y^* X^* */
      addmul_hmatrix_amatrix_amatrix(CONJ(alpha), true, y, true, xf, true,
				     z->f);
    else {
      if (xf->rows > xf->cols) {
	/* Compute rkmatrix X Y = X (Y^* I^*)^* */
	xy = init_rkmatrix(&tmp1, rc->size, cc->size, xf->cols);
	copy_amatrix(false, xf, &xy->A);
	id = init_amatrix(&tmp2, xf->cols, xf->cols);
	identity_amatrix(id);
	clear_amatrix(&xy->B);
	addmul_hmatrix_amatrix_amatrix(CONJ(alpha), true, y, true, id, false,
				       &xy->B);

	/* Add rkmatrix to Z */
	add_rkmatrix_hmatrix(1.0, xy, tm, eps, z);

	/* Clean up */
	uninit_amatrix(id);
	uninit_rkmatrix(xy);
      }
      else {
	/* Compute rkmatrix X Y = I (Y^* X^*)^* */
	xy = init_rkmatrix(&tmp1, rc->size, cc->size, xf->rows);
	identity_amatrix(&xy->A);
	clear_amatrix(&xy->B);
	addmul_hmatrix_amatrix_amatrix(CONJ(alpha), true, y, true, xf, false,
				       &xy->B);

	/* Add rkmatrix to Z */
	add_rkmatrix_hmatrix(1.0, xy, tm, eps, z);

	/* Clean up */
	uninit_rkmatrix(xy);
      }
    }
  }
  else if (x->r) {
    /* Compute rkmatrix X Y = A (Y^* B)^* */
    xy = init_rkmatrix(&tmp1, rc->size, cc->size, x->r->k);
    copy_amatrix(false, &x->r->A, &xy->A);
    clear_amatrix(&xy->B);
    addmul_hmatrix_amatrix_amatrix(CONJ(alpha), true, y, false, &x->r->B,
				   false, &xy->B);

    /* Add rkmatrix to Z */
    add_rkmatrix_hmatrix(1.0, xy, tm, eps, z);

    /* Clean up */
    uninit_rkmatrix(xy);
  }
  else {
    assert(x->son);

    if (y->f) {
      yf = y->f;
      if (z->f)			/* Z = Z + X Y */
	addmul_hmatrix_amatrix_amatrix(alpha, false, x, false, yf, false,
				       z->f);
      else {
	if (yf->cols > yf->rows) {
	  /* Compute rkmatrix X Y = (X I) (Y^*)^* */
	  xy = init_rkmatrix(&tmp1, rc->size, cc->size, yf->rows);
	  copy_amatrix(true, yf, &xy->B);
	  id = init_amatrix(&tmp2, yf->rows, yf->rows);
	  identity_amatrix(id);
	  clear_amatrix(&xy->A);
	  addmul_hmatrix_amatrix_amatrix(alpha, false, x, false, id, false,
					 &xy->A);

	  /* Add rkmatrix to Z */
	  add_rkmatrix_hmatrix(1.0, xy, tm, eps, z);

	  uninit_amatrix(id);
	  uninit_rkmatrix(xy);
	}
	else {
	  /* Compute rkmatrix X Y = (X Y) I^* */
	  xy = init_rkmatrix(&tmp1, rc->size, cc->size, yf->cols);
	  identity_amatrix(&xy->B);
	  clear_amatrix(&xy->A);
	  addmul_hmatrix_amatrix_amatrix(alpha, false, x, false, yf, false,
					 &xy->A);

	  /* Add rkmatrix to Z */
	  add_rkmatrix_hmatrix(1.0, xy, tm, eps, z);

	  /* Clean up */
	  uninit_rkmatrix(xy);
	}
      }
    }
    else if (y->r) {
      /* Compute rkmatrix X Y = (X A) B^* */
      xy = init_rkmatrix(&tmp1, rc->size, cc->size, y->r->k);
      copy_amatrix(false, &y->r->B, &xy->B);
      clear_amatrix(&xy->A);
      addmul_hmatrix_amatrix_amatrix(alpha, false, x, false, &y->r->A, false,
				     &xy->A);

      /* Add rkmatrix to Z */
      add_rkmatrix_hmatrix(1.0, xy, tm, eps, z);

      /* Clean up */
      uninit_rkmatrix(xy);
    }
    else {
      assert(y->son);
      assert(x->csons == y->rsons);

      rsons = x->rsons;
      msons = x->csons;
      csons = y->csons;

      if (z->f) {
	zf = z->f;

	if (y->son[0]->cc == y->cc) {
	  if (x->son[0]->rc == x->rc) {
	    /* Row and column not subdivided */
	    assert(csons == 1);
	    assert(rsons == 1);

	    /* Multiplications of submatrices */
	    for (j = 0; j < msons; j++)
	      addmul_nn_hmatrix(alpha, x->son[j], y->son[j], tm, eps, z);
	  }
	  else {
	    /* Only row subdivided */
	    assert(csons == 1);

	    roff = 0;
	    for (i = 0; i < rsons; i++) {
	      /* Create submatrix of z */
	      zf1 = init_sub_amatrix(&tmp2, zf, rc->son[i]->size, roff,
				     cc->size, 0);
	      z1 = init_hmatrix(&tmp3, rc->son[i], cc);
	      z1->f = zf1;

	      /* Multiplications of submatrices */
	      for (j = 0; j < msons; j++)
		addmul_nn_hmatrix(alpha, x->son[i + j * rsons], y->son[j], tm,
				  eps, z1);

	      /* Clean up */
	      z1->f = 0;
	      uninit_hmatrix(z1);
	      uninit_amatrix(zf1);

	      roff += rc->son[i]->size;
	    }
	    assert(roff == rc->size);
	  }
	}
	else {
	  if (x->son[0]->rc == x->rc) {
	    /* Only column subdivided */
	    assert(rsons == 1);

	    coff = 0;
	    for (k = 0; k < csons; k++) {
	      /* Create submatrix of z */
	      zf1 = init_sub_amatrix(&tmp2, zf, rc->size, 0, cc->son[k]->size,
				     coff);
	      z1 = init_hmatrix(&tmp3, rc, cc->son[k]);
	      z1->f = zf1;

	      /* Multiplications of submatrices */
	      for (j = 0; j < msons; j++)
		addmul_nn_hmatrix(alpha, x->son[j], y->son[j + k * msons], tm,
				  eps, z1);

	      /* Clean up */
	      z1->f = 0;
	      uninit_hmatrix(z1);
	      uninit_amatrix(zf1);

	      coff += cc->son[k]->size;
	    }
	    assert(coff == cc->size);
	  }
	  else {
	    /* Row and column subdivided */
	    coff = 0;
	    for (k = 0; k < csons; k++) {

	      roff = 0;
	      for (i = 0; i < rsons; i++) {
		/* Create submatrix of z */
		zf1 = init_sub_amatrix(&tmp2, zf, rc->son[i]->size, roff,
				       cc->son[k]->size, coff);
		z1 = init_hmatrix(&tmp3, rc->son[i], cc->son[k]);
		z1->f = zf1;

		/* Multiplications of submatrices */
		for (j = 0; j < msons; j++)
		  addmul_nn_hmatrix(alpha, x->son[i + j * rsons],
				    y->son[j + k * msons], tm, eps, z1);

		/* Clean up */
		z1->f = 0;
		uninit_hmatrix(z1);
		uninit_amatrix(zf1);

		roff += rc->son[i]->size;
	      }
	      assert(roff == rc->size);

	      coff += cc->son[k]->size;
	    }
	    assert(coff == cc->size);
	  }
	}
      }
      else if (z->r) {
	ztmp = split_rkmatrix(z->r, z->rc, z->cc, (x->rc != x->son[0]->rc),
			      (y->cc != y->son[0]->cc), false);

	for (k = 0; k < csons; k++)
	  for (i = 0; i < rsons; i++)
	    for (j = 0; j < msons; j++)
	      addmul_nn_hmatrix(alpha, x->son[i + j * rsons],
				y->son[j + k * msons], tm, eps,
				ztmp->son[i + k * rsons]);

	xy = merge_hmatrix_rkmatrix(ztmp, tm, eps);
	add_rkmatrix(1.0, xy, tm, eps, z->r);

	del_rkmatrix(xy);
	del_hmatrix(ztmp);
      }
      else {
	for (k = 0; k < csons; k++)
	  for (i = 0; i < rsons; i++)
	    for (j = 0; j < msons; j++)
	      addmul_nn_hmatrix(alpha, x->son[i + j * rsons],
				y->son[j + k * msons], tm, eps,
				z->son[i + k * rsons]);
      }
    }
  }
}

static void
addmul_nt_hmatrix(field alpha, pchmatrix x, pchmatrix y,
		  pctruncmode tm, real eps, phmatrix z)
{
  prkmatrix xy;
  pamatrix  id;
  pcamatrix xf, yf;
  pamatrix  zf, zf1;
  phmatrix  z1, ztmp;
  rkmatrix  tmp1;
  amatrix   tmp2;
  hmatrix   tmp3;
  pccluster rc, cc;
  uint      rsons, msons, csons;
  uint      roff, coff;
  uint      i, j, k;

  assert(z->rc == x->rc);
  assert(z->cc == y->rc);
  assert(x->cc == y->cc);

  rc = z->rc;
  cc = z->cc;

  if (x->f) {
    xf = x->f;
    if (z->f)			/* Z = Z + X Y^*   <=>   Z^* = Z^* + Y X^* */
      addmul_hmatrix_amatrix_amatrix(CONJ(alpha), false, y, true, xf, true,
				     z->f);
    else {
      if (xf->rows > xf->cols) {
	/* Compute rkmatrix X Y^* = X (Y I^*)^* */
	xy = init_rkmatrix(&tmp1, rc->size, cc->size, xf->cols);
	copy_amatrix(false, xf, &xy->A);
	id = init_amatrix(&tmp2, xf->cols, xf->cols);
	identity_amatrix(id);
	clear_amatrix(&xy->B);
	addmul_hmatrix_amatrix_amatrix(CONJ(alpha), false, y, true, id, false,
				       &xy->B);

	/* Add rkmatrix to Z */
	add_rkmatrix_hmatrix(1.0, xy, tm, eps, z);

	/* Clean up */
	uninit_amatrix(id);
	uninit_rkmatrix(xy);
      }
      else {
	/* Compute rkmatrix X Y^* = I (Y X^*)^* */
	xy = init_rkmatrix(&tmp1, rc->size, cc->size, xf->rows);
	identity_amatrix(&xy->A);
	clear_amatrix(&xy->B);
	addmul_hmatrix_amatrix_amatrix(CONJ(alpha), false, y, true, xf, false,
				       &xy->B);

	/* Add rkmatrix to Z */
	add_rkmatrix_hmatrix(1.0, xy, tm, eps, z);

	/* Clean up */
	uninit_rkmatrix(xy);
      }
    }
  }
  else if (x->r) {
    /* Compute rkmatrix X Y^* = A (Y B)^* */
    xy = init_rkmatrix(&tmp1, rc->size, cc->size, x->r->k);
    copy_amatrix(false, &x->r->A, &xy->A);
    clear_amatrix(&xy->B);
    addmul_hmatrix_amatrix_amatrix(CONJ(alpha), false, y, false, &x->r->B,
				   false, &xy->B);

    /* Add rkmatrix to Z */
    add_rkmatrix_hmatrix(1.0, xy, tm, eps, z);

    /* Clean up */
    uninit_rkmatrix(xy);
  }
  else {
    assert(x->son);

    if (y->f) {
      yf = y->f;
      if (z->f)			/* Z = Z + X Y^* */
	addmul_hmatrix_amatrix_amatrix(alpha, false, x, true, yf, false,
				       z->f);
      else {
	if (yf->cols < yf->rows) {
	  /* Compute rkmatrix X Y^* = (X I) Y^* */
	  xy = init_rkmatrix(&tmp1, rc->size, cc->size, yf->cols);
	  copy_amatrix(false, yf, &xy->B);
	  id = init_amatrix(&tmp2, yf->cols, yf->cols);
	  identity_amatrix(id);
	  clear_amatrix(&xy->A);
	  addmul_hmatrix_amatrix_amatrix(alpha, false, x, false, id, false,
					 &xy->A);

	  /* Add rkmatrix to Z */
	  add_rkmatrix_hmatrix(1.0, xy, tm, eps, z);

	  uninit_amatrix(id);
	  uninit_rkmatrix(xy);
	}
	else {
	  /* Compute rkmatrix X Y^* = (X Y^*) I^* */
	  xy = init_rkmatrix(&tmp1, rc->size, cc->size, yf->rows);
	  identity_amatrix(&xy->B);
	  clear_amatrix(&xy->A);
	  addmul_hmatrix_amatrix_amatrix(alpha, false, x, true, yf, false,
					 &xy->A);

	  /* Add rkmatrix to Z */
	  add_rkmatrix_hmatrix(1.0, xy, tm, eps, z);

	  /* Clean up */
	  uninit_rkmatrix(xy);
	}
      }
    }
    else if (y->r) {
      /* Compute rkmatrix X Y^* = (X B) A^* */
      xy = init_rkmatrix(&tmp1, rc->size, cc->size, y->r->k);
      copy_amatrix(false, &y->r->A, &xy->B);
      clear_amatrix(&xy->A);
      addmul_hmatrix_amatrix_amatrix(alpha, false, x, false, &y->r->B, false,
				     &xy->A);

      /* Add rkmatrix to Z */
      add_rkmatrix_hmatrix(1.0, xy, tm, eps, z);

      /* Clean up */
      uninit_rkmatrix(xy);
    }
    else {
      assert(y->son);
      assert(x->csons == y->csons);

      rsons = x->rsons;
      msons = x->csons;
      csons = y->rsons;

      if (z->f) {
	zf = z->f;

	if (y->son[0]->rc == y->rc) {
	  if (x->son[0]->rc == x->rc) {
	    /* Row and column not subdivided */
	    assert(csons == 1);
	    assert(rsons == 1);

	    /* Multiplications of submatrices */
	    for (j = 0; j < msons; j++)
	      addmul_nt_hmatrix(alpha, x->son[j], y->son[j], tm, eps, z);
	  }
	  else {
	    /* Only row subdivided */
	    assert(csons == 1);

	    roff = 0;
	    for (i = 0; i < rsons; i++) {
	      /* Create submatrix of z */
	      zf1 = init_sub_amatrix(&tmp2, zf, rc->son[i]->size, roff,
				     cc->size, 0);
	      z1 = init_hmatrix(&tmp3, rc->son[i], cc);
	      z1->f = zf1;

	      /* Multiplications of submatrices */
	      for (j = 0; j < msons; j++)
		addmul_nt_hmatrix(alpha, x->son[i + j * rsons], y->son[j], tm,
				  eps, z1);

	      /* Clean up */
	      z1->f = 0;
	      uninit_hmatrix(z1);
	      uninit_amatrix(zf1);

	      roff += rc->son[i]->size;
	    }
	    assert(roff == rc->size);
	  }
	}
	else {
	  if (x->son[0]->rc == x->rc) {
	    /* Only column subdivided */
	    assert(rsons == 1);

	    coff = 0;
	    for (k = 0; k < csons; k++) {
	      /* Create submatrix of z */
	      zf1 = init_sub_amatrix(&tmp2, zf, rc->size, 0, cc->son[k]->size,
				     coff);
	      z1 = init_hmatrix(&tmp3, rc, cc->son[k]);
	      z1->f = zf1;

	      /* Multiplications of submatrices */
	      for (j = 0; j < msons; j++)
		addmul_nt_hmatrix(alpha, x->son[j], y->son[k + j * csons], tm,
				  eps, z1);

	      /* Clean up */
	      z1->f = 0;
	      uninit_hmatrix(z1);
	      uninit_amatrix(zf1);

	      coff += cc->son[k]->size;
	    }
	    assert(coff == cc->size);
	  }
	  else {
	    /* Row and column subdivided */
	    coff = 0;
	    for (k = 0; k < csons; k++) {

	      roff = 0;
	      for (i = 0; i < rsons; i++) {
		/* Create submatrix of z */
		zf1 = init_sub_amatrix(&tmp2, zf, rc->son[i]->size, roff,
				       cc->son[k]->size, coff);
		z1 = init_hmatrix(&tmp3, rc->son[i], cc->son[k]);
		z1->f = zf1;

		/* Multiplications of submatrices */
		for (j = 0; j < msons; j++)
		  addmul_nt_hmatrix(alpha, x->son[i + j * rsons],
				    y->son[k + j * csons], tm, eps, z1);

		/* Clean up */
		z1->f = 0;
		uninit_hmatrix(z1);
		uninit_amatrix(zf1);

		roff += rc->son[i]->size;
	      }
	      assert(roff == rc->size);

	      coff += cc->son[k]->size;
	    }
	    assert(coff == cc->size);
	  }
	}
      }
      else if (z->r) {
	ztmp = split_rkmatrix(z->r, z->rc, z->cc, (x->rc != x->son[0]->rc),
			      (y->rc != y->son[0]->rc), false);

	for (k = 0; k < csons; k++)
	  for (i = 0; i < rsons; i++)
	    for (j = 0; j < msons; j++)
	      addmul_nt_hmatrix(alpha, x->son[i + j * rsons],
				y->son[k + j * csons], tm, eps,
				ztmp->son[i + k * rsons]);

	xy = merge_hmatrix_rkmatrix(ztmp, tm, eps);
	add_rkmatrix(1.0, xy, tm, eps, z->r);

	del_rkmatrix(xy);
	del_hmatrix(ztmp);
      }
      else {
	for (k = 0; k < csons; k++)
	  for (i = 0; i < rsons; i++)
	    for (j = 0; j < msons; j++)
	      addmul_nt_hmatrix(alpha, x->son[i + j * rsons],
				y->son[k + j * csons], tm, eps,
				z->son[i + k * rsons]);
      }
    }
  }
}

static void
addmul_tn_hmatrix(field alpha, pchmatrix x, pchmatrix y,
		  pctruncmode tm, real eps, phmatrix z)
{
  prkmatrix xy;
  pamatrix  id;
  pcamatrix xf, yf;
  pamatrix  zf, zf1;
  phmatrix  z1, ztmp;
  rkmatrix  tmp1;
  amatrix   tmp2;
  hmatrix   tmp3;
  pccluster rc, cc;
  uint      rsons, msons, csons;
  uint      roff, coff;
  uint      i, j, k;

  assert(z->rc == x->cc);
  assert(z->cc == y->cc);
  assert(x->rc == y->rc);

  rc = z->rc;
  cc = z->cc;

  if (x->f) {
    xf = x->f;
    if (z->f)			/* Z = Z + X^* Y   <=>   Z^* = Z^* + Y^* X */
      addmul_hmatrix_amatrix_amatrix(CONJ(alpha), true, y, false, xf, true,
				     z->f);
    else {
      if (xf->rows < xf->cols) {
	/* Compute rkmatrix X^* Y = X^* (Y^* I^*)^* */
	xy = init_rkmatrix(&tmp1, rc->size, cc->size, xf->rows);
	copy_amatrix(true, xf, &xy->A);
	id = init_amatrix(&tmp2, xf->rows, xf->rows);
	identity_amatrix(id);
	clear_amatrix(&xy->B);
	addmul_hmatrix_amatrix_amatrix(CONJ(alpha), true, y, true, id, false,
				       &xy->B);

	/* Add rkmatrix to Z */
	add_rkmatrix_hmatrix(1.0, xy, tm, eps, z);

	/* Clean up */
	uninit_amatrix(id);
	uninit_rkmatrix(xy);
      }
      else {
	/* Compute rkmatrix X^* Y = I (Y^* X)^* */
	xy = init_rkmatrix(&tmp1, rc->size, cc->size, xf->cols);
	identity_amatrix(&xy->A);
	clear_amatrix(&xy->B);
	addmul_hmatrix_amatrix_amatrix(CONJ(alpha), true, y, false, xf, false,
				       &xy->B);

	/* Add rkmatrix to Z */
	add_rkmatrix_hmatrix(1.0, xy, tm, eps, z);

	/* Clean up */
	uninit_rkmatrix(xy);
      }
    }
  }
  else if (x->r) {
    /* Compute rkmatrix X^* Y = B (Y^* A)^* */
    xy = init_rkmatrix(&tmp1, rc->size, cc->size, x->r->k);
    copy_amatrix(false, &x->r->B, &xy->A);
    clear_amatrix(&xy->B);
    addmul_hmatrix_amatrix_amatrix(CONJ(alpha), true, y, false, &x->r->A,
				   false, &xy->B);

    /* Add rkmatrix to Z */
    add_rkmatrix_hmatrix(1.0, xy, tm, eps, z);

    /* Clean up */
    uninit_rkmatrix(xy);
  }
  else {
    assert(x->son);

    if (y->f) {
      yf = y->f;
      if (z->f)			/* Z = Z + X^* Y */
	addmul_hmatrix_amatrix_amatrix(alpha, true, x, false, yf, false,
				       z->f);
      else {
	if (yf->cols > yf->rows) {
	  /* Compute rkmatrix X^* Y = (X^* I) (Y^*)^* */
	  xy = init_rkmatrix(&tmp1, rc->size, cc->size, yf->rows);
	  copy_amatrix(true, yf, &xy->B);
	  id = init_amatrix(&tmp2, yf->rows, yf->rows);
	  identity_amatrix(id);
	  clear_amatrix(&xy->A);
	  addmul_hmatrix_amatrix_amatrix(alpha, true, x, false, id, false,
					 &xy->A);

	  /* Add rkmatrix to Z */
	  add_rkmatrix_hmatrix(1.0, xy, tm, eps, z);

	  uninit_amatrix(id);
	  uninit_rkmatrix(xy);
	}
	else {
	  /* Compute rkmatrix X^* Y = (X^* Y) I^* */
	  xy = init_rkmatrix(&tmp1, rc->size, cc->size, yf->cols);
	  identity_amatrix(&xy->B);
	  clear_amatrix(&xy->A);
	  addmul_hmatrix_amatrix_amatrix(alpha, true, x, false, yf, false,
					 &xy->A);

	  /* Add rkmatrix to Z */
	  add_rkmatrix_hmatrix(1.0, xy, tm, eps, z);

	  /* Clean up */
	  uninit_rkmatrix(xy);
	}
      }
    }
    else if (y->r) {
      /* Compute rkmatrix X^* Y = (X^* A) B^* */
      xy = init_rkmatrix(&tmp1, rc->size, cc->size, y->r->k);
      copy_amatrix(false, &y->r->B, &xy->B);
      clear_amatrix(&xy->A);
      addmul_hmatrix_amatrix_amatrix(alpha, true, x, false, &y->r->A, false,
				     &xy->A);

      /* Add rkmatrix to Z */
      add_rkmatrix_hmatrix(1.0, xy, tm, eps, z);

      /* Clean up */
      uninit_rkmatrix(xy);
    }
    else {
      assert(y->son);
      assert(x->rsons == y->rsons);

      rsons = x->csons;
      msons = x->rsons;
      csons = y->csons;

      if (z->f) {
	zf = z->f;

	if (y->son[0]->cc == y->cc) {
	  if (x->son[0]->cc == x->cc) {
	    /* Row and column not subdivided */
	    assert(csons == 1);
	    assert(rsons == 1);

	    /* Multiplications of submatrices */
	    for (j = 0; j < msons; j++)
	      addmul_tn_hmatrix(alpha, x->son[j], y->son[j], tm, eps, z);
	  }
	  else {
	    /* Only row subdivided */
	    assert(csons == 1);

	    roff = 0;
	    for (i = 0; i < rsons; i++) {
	      /* Create submatrix of z */
	      zf1 = init_sub_amatrix(&tmp2, zf, rc->son[i]->size, roff,
				     cc->size, 0);
	      z1 = init_hmatrix(&tmp3, rc->son[i], cc);
	      z1->f = zf1;

	      /* Multiplications of submatrices */
	      for (j = 0; j < msons; j++)
		addmul_tn_hmatrix(alpha, x->son[j + i * msons], y->son[j], tm,
				  eps, z1);

	      /* Clean up */
	      z1->f = 0;
	      uninit_hmatrix(z1);
	      uninit_amatrix(zf1);

	      roff += rc->son[i]->size;
	    }
	    assert(roff == rc->size);
	  }
	}
	else {
	  if (x->son[0]->cc == x->cc) {
	    /* Only column subdivided */
	    assert(rsons == 1);

	    coff = 0;
	    for (k = 0; k < csons; k++) {
	      /* Create submatrix of z */
	      zf1 = init_sub_amatrix(&tmp2, zf, rc->size, 0, cc->son[k]->size,
				     coff);
	      z1 = init_hmatrix(&tmp3, rc, cc->son[k]);
	      z1->f = zf1;

	      /* Multiplications of submatrices */
	      for (j = 0; j < msons; j++)
		addmul_tn_hmatrix(alpha, x->son[j], y->son[j + k * msons], tm,
				  eps, z1);

	      /* Clean up */
	      z1->f = 0;
	      uninit_hmatrix(z1);
	      uninit_amatrix(zf1);

	      coff += cc->son[k]->size;
	    }
	    assert(coff == cc->size);
	  }
	  else {
	    /* Row and column subdivided */
	    coff = 0;
	    for (k = 0; k < csons; k++) {

	      roff = 0;
	      for (i = 0; i < rsons; i++) {
		/* Create submatrix of z */
		zf1 = init_sub_amatrix(&tmp2, zf, rc->son[i]->size, roff,
				       cc->son[k]->size, coff);
		z1 = init_hmatrix(&tmp3, rc->son[i], cc->son[k]);
		z1->f = zf1;

		/* Multiplications of submatrices */
		for (j = 0; j < msons; j++)
		  addmul_tn_hmatrix(alpha, x->son[j + i * msons],
				    y->son[j + k * msons], tm, eps, z1);

		/* Clean up */
		z1->f = 0;
		uninit_hmatrix(z1);
		uninit_amatrix(zf1);

		roff += rc->son[i]->size;
	      }
	      assert(roff == rc->size);

	      coff += cc->son[k]->size;
	    }
	    assert(coff == cc->size);
	  }
	}
      }
      else if (z->r) {
	ztmp = split_rkmatrix(z->r, z->rc, z->cc, (x->cc != x->son[0]->cc),
			      (y->cc != y->son[0]->cc), false);

	for (k = 0; k < csons; k++)
	  for (i = 0; i < rsons; i++)
	    for (j = 0; j < msons; j++)
	      addmul_tn_hmatrix(alpha, x->son[j + i * msons],
				y->son[j + k * msons], tm, eps,
				ztmp->son[i + k * rsons]);

	xy = merge_hmatrix_rkmatrix(ztmp, tm, eps);
	add_rkmatrix(1.0, xy, tm, eps, z->r);

	del_rkmatrix(xy);
	del_hmatrix(ztmp);
      }
      else {
	for (k = 0; k < csons; k++)
	  for (i = 0; i < rsons; i++)
	    for (j = 0; j < msons; j++)
	      addmul_tn_hmatrix(alpha, x->son[j + i * msons],
				y->son[j + k * msons], tm, eps,
				z->son[i + k * rsons]);
      }
    }
  }
}

static void
addmul_tt_hmatrix(field alpha, pchmatrix x, pchmatrix y,
		  pctruncmode tm, real eps, phmatrix z)
{
  prkmatrix xy;
  pamatrix  id;
  pcamatrix xf, yf;
  pamatrix  zf, zf1;
  phmatrix  z1, ztmp;
  rkmatrix  tmp1;
  amatrix   tmp2;
  hmatrix   tmp3;
  pccluster rc, cc;
  uint      rsons, msons, csons;
  uint      roff, coff;
  uint      i, j, k;

  assert(z->rc == x->cc);
  assert(z->cc == y->rc);
  assert(x->rc == y->cc);

  rc = z->rc;
  cc = z->cc;

  if (x->f) {
    xf = x->f;
    if (z->f)			/* Z = Z + X^* Y^*   <=>   Z^* = Z^* + Y X */
      addmul_hmatrix_amatrix_amatrix(CONJ(alpha), false, y, false, xf, true,
				     z->f);
    else {
      if (xf->rows < xf->cols) {
	/* Compute rkmatrix X^* Y^* = X^* (Y I^*)^* */
	xy = init_rkmatrix(&tmp1, rc->size, cc->size, xf->rows);
	copy_amatrix(true, xf, &xy->A);
	id = init_amatrix(&tmp2, xf->rows, xf->rows);
	identity_amatrix(id);
	clear_amatrix(&xy->B);
	addmul_hmatrix_amatrix_amatrix(CONJ(alpha), false, y, true, id, false,
				       &xy->B);

	/* Add rkmatrix to Z */
	add_rkmatrix_hmatrix(1.0, xy, tm, eps, z);

	/* Clean up */
	uninit_amatrix(id);
	uninit_rkmatrix(xy);
      }
      else {
	/* Compute rkmatrix X^* Y^* = I (Y X)^* */
	xy = init_rkmatrix(&tmp1, rc->size, cc->size, xf->cols);
	identity_amatrix(&xy->A);
	clear_amatrix(&xy->B);
	addmul_hmatrix_amatrix_amatrix(CONJ(alpha), false, y, false, xf,
				       false, &xy->B);

	/* Add rkmatrix to Z */
	add_rkmatrix_hmatrix(1.0, xy, tm, eps, z);

	/* Clean up */
	uninit_rkmatrix(xy);
      }
    }
  }
  else if (x->r) {
    /* Compute rkmatrix X^* Y^* = B (Y A)^* */
    xy = init_rkmatrix(&tmp1, rc->size, cc->size, x->r->k);
    copy_amatrix(false, &x->r->B, &xy->A);
    clear_amatrix(&xy->B);
    addmul_hmatrix_amatrix_amatrix(CONJ(alpha), false, y, false, &x->r->A,
				   false, &xy->B);

    /* Add rkmatrix to Z */
    add_rkmatrix_hmatrix(1.0, xy, tm, eps, z);

    /* Clean up */
    uninit_rkmatrix(xy);
  }
  else {
    assert(x->son);

    if (y->f) {
      yf = y->f;
      if (z->f)			/* Z = Z + X^* Y^* */
	addmul_hmatrix_amatrix_amatrix(alpha, true, x, true, yf, false, z->f);
      else {
	if (yf->cols < yf->rows) {
	  /* Compute rkmatrix X^* Y^* = (X^* I) Y^* */
	  xy = init_rkmatrix(&tmp1, rc->size, cc->size, yf->cols);
	  copy_amatrix(false, yf, &xy->B);
	  id = init_amatrix(&tmp2, yf->cols, yf->cols);
	  identity_amatrix(id);
	  clear_amatrix(&xy->A);
	  addmul_hmatrix_amatrix_amatrix(alpha, true, x, false, id, false,
					 &xy->A);

	  /* Add rkmatrix to Z */
	  add_rkmatrix_hmatrix(1.0, xy, tm, eps, z);

	  uninit_amatrix(id);
	  uninit_rkmatrix(xy);
	}
	else {
	  /* Compute rkmatrix X^* Y^* = (X^* Y^*) I^* */
	  xy = init_rkmatrix(&tmp1, rc->size, cc->size, yf->rows);
	  identity_amatrix(&xy->B);
	  clear_amatrix(&xy->A);
	  addmul_hmatrix_amatrix_amatrix(alpha, true, x, true, yf, false,
					 &xy->A);

	  /* Add rkmatrix to Z */
	  add_rkmatrix_hmatrix(1.0, xy, tm, eps, z);

	  /* Clean up */
	  uninit_rkmatrix(xy);
	}
      }
    }
    else if (y->r) {
      /* Compute rkmatrix X^* Y^* = (X^* B) A^* */
      xy = init_rkmatrix(&tmp1, rc->size, cc->size, y->r->k);
      copy_amatrix(false, &y->r->A, &xy->B);
      clear_amatrix(&xy->A);
      addmul_hmatrix_amatrix_amatrix(alpha, true, x, false, &y->r->B, false,
				     &xy->A);

      /* Add rkmatrix to Z */
      add_rkmatrix_hmatrix(1.0, xy, tm, eps, z);

      /* Clean up */
      uninit_rkmatrix(xy);
    }
    else {
      assert(y->son);
      assert(x->rsons == y->csons);

      rsons = x->csons;
      msons = x->rsons;
      csons = y->rsons;

      if (z->f) {
	zf = z->f;

	if (y->son[0]->rc == y->rc) {
	  if (x->son[0]->cc == x->cc) {
	    /* Row and column not subdivided */
	    assert(csons == 1);
	    assert(rsons == 1);

	    /* Multiplications of submatrices */
	    for (j = 0; j < msons; j++)
	      addmul_tt_hmatrix(alpha, x->son[j], y->son[j], tm, eps, z);
	  }
	  else {
	    /* Only row subdivided */
	    assert(csons == 1);

	    roff = 0;
	    for (i = 0; i < rsons; i++) {
	      /* Create submatrix of z */
	      zf1 = init_sub_amatrix(&tmp2, zf, rc->son[i]->size, roff,
				     cc->size, 0);
	      z1 = init_hmatrix(&tmp3, rc->son[i], cc);
	      z1->f = zf1;

	      /* Multiplications of submatrices */
	      for (j = 0; j < msons; j++)
		addmul_tt_hmatrix(alpha, x->son[j + i * msons], y->son[j], tm,
				  eps, z1);

	      /* Clean up */
	      z1->f = 0;
	      uninit_hmatrix(z1);
	      uninit_amatrix(zf1);

	      roff += rc->son[i]->size;
	    }
	    assert(roff == rc->size);
	  }
	}
	else {
	  if (x->son[0]->cc == x->cc) {
	    /* Only column subdivided */
	    assert(rsons == 1);

	    coff = 0;
	    for (k = 0; k < csons; k++) {
	      /* Create submatrix of z */
	      zf1 = init_sub_amatrix(&tmp2, zf, rc->size, 0, cc->son[k]->size,
				     coff);
	      z1 = init_hmatrix(&tmp3, rc, cc->son[k]);
	      z1->f = zf1;

	      /* Multiplications of submatrices */
	      for (j = 0; j < msons; j++)
		addmul_tt_hmatrix(alpha, x->son[j], y->son[k + j * csons], tm,
				  eps, z1);

	      /* Clean up */
	      z1->f = 0;
	      uninit_hmatrix(z1);
	      uninit_amatrix(zf1);

	      coff += cc->son[k]->size;
	    }
	    assert(coff == cc->size);
	  }
	  else {
	    /* Row and column subdivided */
	    coff = 0;
	    for (k = 0; k < csons; k++) {

	      roff = 0;
	      for (i = 0; i < rsons; i++) {
		/* Create submatrix of z */
		zf1 = init_sub_amatrix(&tmp2, zf, rc->son[i]->size, roff,
				       cc->son[k]->size, coff);
		z1 = init_hmatrix(&tmp3, rc->son[i], cc->son[k]);
		z1->f = zf1;

		/* Multiplications of submatrices */
		for (j = 0; j < msons; j++)
		  addmul_tt_hmatrix(alpha, x->son[j + i * msons],
				    y->son[k + j * csons], tm, eps, z1);

		/* Clean up */
		z1->f = 0;
		uninit_hmatrix(z1);
		uninit_amatrix(zf1);

		roff += rc->son[i]->size;
	      }
	      assert(roff == rc->size);

	      coff += cc->son[k]->size;
	    }
	    assert(coff == cc->size);
	  }
	}
      }
      else if (z->r) {
	ztmp = split_rkmatrix(z->r, z->rc, z->cc, (x->cc != x->son[0]->cc),
			      (y->rc != y->son[0]->rc), false);

	for (k = 0; k < csons; k++)
	  for (i = 0; i < rsons; i++)
	    for (j = 0; j < msons; j++)
	      addmul_tt_hmatrix(alpha, x->son[j + i * msons],
				y->son[k + j * csons], tm, eps,
				ztmp->son[i + k * rsons]);

	xy = merge_hmatrix_rkmatrix(ztmp, tm, eps);
	add_rkmatrix(1.0, xy, tm, eps, z->r);

	del_rkmatrix(xy);
	del_hmatrix(ztmp);
      }
      else {
	for (k = 0; k < csons; k++)
	  for (i = 0; i < rsons; i++)
	    for (j = 0; j < msons; j++)
	      addmul_tt_hmatrix(alpha, x->son[j + i * msons],
				y->son[k + j * csons], tm, eps,
				z->son[i + k * rsons]);
      }
    }
  }
}

void
addmul_hmatrix(field alpha, bool xtrans, pchmatrix x, bool ytrans,
	       pchmatrix y, pctruncmode tm, real eps, phmatrix z)
{
  if (xtrans) {
    if (ytrans)
      addmul_tt_hmatrix(alpha, x, y, tm, eps, z);
    else
      addmul_tn_hmatrix(alpha, x, y, tm, eps, z);
  }
  else {
    if (ytrans)
      addmul_nt_hmatrix(alpha, x, y, tm, eps, z);
    else
      addmul_nn_hmatrix(alpha, x, y, tm, eps, z);
  }
}

void
addmul_lower_hmatrix(field alpha, bool xtrans, pchmatrix x, bool ytrans,
		     pchmatrix y, pctruncmode tm, real eps, phmatrix z)
{
  amatrix   tmp1;
  rkmatrix  tmp2;
  pamatrix  xf, yf, mf;
  prkmatrix xr, yr, mr;
  phmatrix  ztmp;
  uint      sons, msons;
  uint      i, j, k;

  assert(z->rc == z->cc);

  if (x->f) {
    xf = x->f;
    if (xtrans) {
      if (xf->cols <= xf->rows) {
	/* Compute M = X^* Y = (Y^* X)^* */
	mf = init_amatrix(&tmp1, xf->cols,
			  (ytrans ? y->rc->size : y->cc->size));
	clear_amatrix(mf);
	addmul_hmatrix_amatrix_amatrix(1.0, !ytrans, y, !xtrans, xf, true,
				       mf);
	add_lower_amatrix_hmatrix(alpha, false, mf, tm, eps, z);
	uninit_amatrix(mf);
      }
      else {
	/* Compute M = X^* Y = X^* (Y^*)^* */
	mr = init_rkmatrix(&tmp2, xf->cols,
			   (ytrans ? y->rc->size : y->cc->size), xf->rows);
	copy_amatrix(true, xf, &mr->A);
	clear_amatrix(&mr->B);
	add_hmatrix_amatrix(1.0, !ytrans, y, &mr->B);
	add_lower_rkmatrix_hmatrix(alpha, mr, tm, eps, z);
	uninit_rkmatrix(mr);
      }
    }
    else {
      if (xf->rows <= xf->cols) {
	/* Compute M = X Y = (Y^* X^*)^* */
	mf = init_amatrix(&tmp1, xf->rows,
			  (ytrans ? y->rc->size : y->cc->size));
	clear_amatrix(mf);
	addmul_hmatrix_amatrix_amatrix(1.0, !ytrans, y, !xtrans, xf, true,
				       mf);
	add_lower_amatrix_hmatrix(alpha, false, mf, tm, eps, z);
	uninit_amatrix(mf);
      }
      else {
	/* Compute M = X Y = X (Y^*)^* */
	mr = init_rkmatrix(&tmp2, xf->rows,
			   (ytrans ? y->rc->size : y->cc->size), xf->cols);
	copy_amatrix(false, xf, &mr->A);
	clear_amatrix(&mr->B);
	add_hmatrix_amatrix(1.0, !ytrans, y, &mr->B);
	add_lower_rkmatrix_hmatrix(alpha, mr, tm, eps, z);
	uninit_rkmatrix(mr);
      }
    }
  }
  else if (y->f) {
    yf = y->f;
    if (ytrans) {
      if (yf->rows <= yf->cols) {
	/* Compute M = X Y^* */
	mf = init_amatrix(&tmp1, (xtrans ? x->cc->size : x->rc->size),
			  yf->rows);
	clear_amatrix(mf);
	addmul_hmatrix_amatrix_amatrix(1.0, xtrans, x, ytrans, yf, false, mf);
	add_lower_amatrix_hmatrix(alpha, false, mf, tm, eps, z);
	uninit_amatrix(mf);
      }
      else {
	/* Compute M = X Y^* */
	mr = init_rkmatrix(&tmp2, (xtrans ? x->cc->size : x->rc->size),
			   yf->rows, yf->cols);
	clear_amatrix(&mr->A);
	add_hmatrix_amatrix(1.0, xtrans, x, &mr->A);
	copy_amatrix(false, yf, &mr->B);
	add_lower_rkmatrix_hmatrix(alpha, mr, tm, eps, z);
	uninit_rkmatrix(mr);
      }
    }
    else {
      if (yf->cols <= yf->rows) {
	/* Compute M = X Y */
	mf = init_amatrix(&tmp1, (xtrans ? x->cc->size : x->rc->size),
			  yf->cols);
	clear_amatrix(mf);
	addmul_hmatrix_amatrix_amatrix(1.0, xtrans, x, ytrans, yf, false, mf);
	add_lower_amatrix_hmatrix(alpha, false, mf, tm, eps, z);
	uninit_amatrix(mf);
      }
      else {
	/* Compute M = X Y = X (Y^*)^* */
	mr = init_rkmatrix(&tmp2, (xtrans ? x->cc->size : x->rc->size),
			   yf->cols, yf->rows);
	clear_amatrix(&mr->A);
	add_hmatrix_amatrix(1.0, xtrans, x, &mr->A);
	copy_amatrix(true, yf, &mr->B);
	add_lower_rkmatrix_hmatrix(alpha, mr, tm, eps, z);
	uninit_rkmatrix(mr);
      }
    }
  }
  else if (x->r) {
    xr = x->r;
    if (xtrans) {
      /* Compute M = B A^* Y = B (Y^* A)^*, where X = A B^* */
      mr = init_rkmatrix(&tmp2, x->cc->size,
			 (ytrans ? y->rc->size : y->cc->size), xr->k);
      copy_amatrix(false, &xr->B, &mr->A);
      clear_amatrix(&mr->B);
      addmul_hmatrix_amatrix_amatrix(1.0, !ytrans, y, false, &xr->A, false,
				     &mr->B);
      add_lower_rkmatrix_hmatrix(alpha, mr, tm, eps, z);
      uninit_rkmatrix(mr);
    }
    else {
      /* Compute M = A B^* Y = A (Y^* B)^*, where X = A B^* */
      mr = init_rkmatrix(&tmp2, x->rc->size,
			 (ytrans ? y->rc->size : y->cc->size), xr->k);
      copy_amatrix(false, &xr->A, &mr->A);
      clear_amatrix(&mr->B);
      addmul_hmatrix_amatrix_amatrix(1.0, !ytrans, y, false, &xr->B, false,
				     &mr->B);
      add_lower_rkmatrix_hmatrix(alpha, mr, tm, eps, z);
      uninit_rkmatrix(mr);
    }
  }
  else if (y->r) {
    yr = y->r;
    if (ytrans) {
      /* Compute M = X B A^* = (X B) A^*, where Y = A B^* */
      mr = init_rkmatrix(&tmp2, (xtrans ? x->cc->size : x->rc->size),
			 y->rc->size, yr->k);
      clear_amatrix(&mr->A);
      addmul_hmatrix_amatrix_amatrix(1.0, xtrans, x, false, &yr->B, false,
				     &mr->A);
      copy_amatrix(false, &yr->A, &mr->B);
      add_lower_rkmatrix_hmatrix(alpha, mr, tm, eps, z);
      uninit_rkmatrix(mr);
    }
    else {
      /* Compute M = X A B^* = (X A) B^*, where Y = A B^* */
      mr = init_rkmatrix(&tmp2, (xtrans ? x->cc->size : x->rc->size),
			 y->cc->size, yr->k);
      clear_amatrix(&mr->A);
      addmul_hmatrix_amatrix_amatrix(1.0, xtrans, x, false, &yr->A, false,
				     &mr->A);
      copy_amatrix(false, &yr->B, &mr->B);
      add_lower_rkmatrix_hmatrix(alpha, mr, tm, eps, z);
      uninit_rkmatrix(mr);
    }
  }
  else {
    assert(x->son);
    assert(y->son);

    ztmp = z;
    if (ztmp->son == 0) {
      assert(z->f);
      ztmp = split_sub_amatrix(z->f, z->rc, z->cc,
			       (xtrans ? x->son[0]->cc !=
				x->cc : x->son[0]->rc != x->rc),
			       (ytrans ? y->son[0]->rc !=
				y->rc : y->son[0]->cc != y->cc));
    }

    assert(ztmp->rsons == ztmp->csons);
    sons = ztmp->rsons;

    if (xtrans) {
      assert(x->csons == ztmp->rsons);
      if (ytrans) {
	assert(x->rsons == y->csons);
	assert(y->rsons == ztmp->csons);

	msons = x->rsons;

	for (i = 0; i < sons; i++) {
	  for (j = 0; j < i; j++)
	    for (k = 0; k < msons; k++)
	      addmul_hmatrix(alpha, xtrans, x->son[k + i * x->rsons], ytrans,
			     y->son[j + k * y->rsons], tm, eps,
			     ztmp->son[i + j * sons]);

	  for (k = 0; k < msons; k++)
	    addmul_lower_hmatrix(alpha, xtrans, x->son[k + i * x->rsons],
				 ytrans, y->son[j + k * y->rsons], tm, eps,
				 ztmp->son[i + j * sons]);
	}
      }
      else {
	assert(x->rsons == y->rsons);
	assert(y->csons == ztmp->csons);

	msons = x->rsons;

	for (i = 0; i < sons; i++) {
	  for (j = 0; j < i; j++)
	    for (k = 0; k < msons; k++)
	      addmul_hmatrix(alpha, xtrans, x->son[k + i * x->rsons], ytrans,
			     y->son[j + k * y->rsons], tm, eps,
			     ztmp->son[i + j * sons]);

	  for (k = 0; k < msons; k++)
	    addmul_lower_hmatrix(alpha, xtrans, x->son[k + i * x->rsons],
				 ytrans, y->son[k + j * y->rsons], tm, eps,
				 ztmp->son[i + j * sons]);
	}
      }
    }
    else {
      assert(x->rsons == ztmp->rsons);
      if (ytrans) {
	assert(x->csons == y->csons);
	assert(y->rsons == ztmp->csons);

	msons = x->csons;

	for (i = 0; i < sons; i++) {
	  for (j = 0; j < i; j++)
	    for (k = 0; k < msons; k++)
	      addmul_hmatrix(alpha, xtrans, x->son[i + k * x->rsons], ytrans,
			     y->son[j + k * y->rsons], tm, eps,
			     ztmp->son[i + j * sons]);

	  for (k = 0; k < msons; k++)
	    addmul_lower_hmatrix(alpha, xtrans, x->son[i + k * x->rsons],
				 ytrans, y->son[j + k * y->rsons], tm, eps,
				 ztmp->son[i + j * sons]);
	}
      }
      else {
	assert(x->csons == y->rsons);
	assert(y->csons == ztmp->csons);

	msons = x->csons;

	for (i = 0; i < sons; i++) {
	  for (j = 0; j < i; j++)
	    for (k = 0; k < msons; k++)
	      addmul_hmatrix(alpha, xtrans, x->son[i + k * x->rsons], ytrans,
			     y->son[k + j * y->rsons], tm, eps,
			     ztmp->son[i + j * sons]);

	  for (k = 0; k < msons; k++)
	    addmul_lower_hmatrix(alpha, xtrans, x->son[i + k * x->rsons],
				 ytrans, y->son[k + j * y->rsons], tm, eps,
				 ztmp->son[i + j * sons]);
	}
      }
    }

    if (ztmp != z)
      del_hmatrix(ztmp);
  }
}

/* ------------------------------------------------------------
 * Inversion of an H-matrix
 * ------------------------------------------------------------ */

void
invert_hmatrix(phmatrix a, phmatrix work, pctruncmode tm, real eps)
{
  uint      sons;
  uint      i, j, k;

  assert(a->rc == a->cc);
  assert(a->rc == work->rc);
  assert(a->cc == work->cc);
  assert(a->rsons == work->rsons);
  assert(a->csons == work->csons);

  if (a->f)
    qrinvert_amatrix(a->f);
  else {
    assert(a->son != 0);
    assert(a->rsons == a->csons);

    sons = a->rsons;

    for (k = 0; k < sons; k++) {
      /* Invert A_{kk} */
      invert_hmatrix(a->son[k + k * sons], work->son[k + k * sons], tm, eps);

      /* Compute W_{ik} = A_{ik} A_{kk}^{-1} */
      for (i = k + 1; i < sons; i++) {
	clear_hmatrix(work->son[i + k * sons]);
	addmul_hmatrix(1.0, false, a->son[i + k * sons], false,
		       a->son[k + k * sons], tm, eps,
		       work->son[i + k * sons]);
      }

      /* Compute W_{kj} = A_{kk}^{-1} A_{kj} */
      for (j = k + 1; j < sons; j++) {
	clear_hmatrix(work->son[k + j * sons]);
	addmul_hmatrix(1.0, false, a->son[k + k * sons], false,
		       a->son[k + j * sons], tm, eps,
		       work->son[k + j * sons]);
      }

      /* Update Schur complement A_{ij} = A_{ij} - W_{ik} A_{kj} */
      for (j = k + 1; j < sons; j++)
	for (i = k + 1; i < sons; i++)
	  addmul_hmatrix(-1.0, false, work->son[i + k * sons], false,
			 a->son[k + j * sons], tm, eps, a->son[i + j * sons]);
    }

    for (k = sons; k-- > 0;) {
      /* Compute (A^{-1})_{ik} = \sum_{j>k} -S_{ij} W_{jk} */
      for (i = k + 1; i < sons; i++) {
	clear_hmatrix(a->son[i + k * sons]);
	for (j = k + 1; j < sons; j++)
	  addmul_hmatrix(-1.0, false, a->son[i + j * sons], false,
			 work->son[j + k * sons], tm, eps,
			 a->son[i + k * sons]);
      }

      /* Compute (A^{-1})_{kj} = \sum_{i>k} -W_{ki} S_{ij} */
      for (j = k + 1; j < sons; j++) {
	clear_hmatrix(a->son[k + j * sons]);
	for (i = k + 1; i < sons; i++)
	  addmul_hmatrix(-1.0, false, work->son[k + i * sons], false,
			 a->son[i + j * sons], tm, eps, a->son[k + j * sons]);
      }

      /* Compute (A^{-1})_{kk} = A_{kk}^{-1} + \sum_{i>k} (A^{-1})_{ki} W_{ik} */
      for (i = k + 1; i < sons; i++)
	addmul_hmatrix(-1.0, false, a->son[k + i * sons], false,
		       work->son[i + k * sons], tm, eps,
		       a->son[k + k * sons]);
    }
  }
}

/* ------------------------------------------------------------
 * Triangular matrices
 * ------------------------------------------------------------ */

phmatrix
clone_lower_hmatrix(bool aunit, pchmatrix a)
{
  phmatrix  b, b1;
  uint      sons;
  uint      i, j;

  assert(a->rc == a->cc);

  b = 0;
  if (a->f) {
    b = new_full_hmatrix(a->rc, a->cc);
    copy_lower_amatrix(a->f, aunit, b->f);
  }
  else {
    assert(a->son != 0);
    assert(a->rsons == a->csons);

    sons = a->rsons;

    b = new_super_hmatrix(a->rc, a->cc, sons, sons);

    for (j = 0; j < sons; j++) {
      for (i = 0; i < j; i++) {
	b1 = new_rk_hmatrix(a->rc->son[i], a->cc->son[j], 0);
	ref_hmatrix(b->son + i + j * sons, b1);
      }

      b1 = clone_lower_hmatrix(aunit, a->son[j + j * sons]);
      ref_hmatrix(b->son + j + j * sons, b1);

      for (i = j + 1; i < sons; i++) {
	b1 = clone_hmatrix(a->son[i + j * sons]);
	ref_hmatrix(b->son + i + j * sons, b1);
      }
    }

    update_hmatrix(b);
  }

  return b;
}

phmatrix
clone_upper_hmatrix(bool aunit, pchmatrix a)
{
  phmatrix  b, b1;
  uint      sons;
  uint      i, j;

  assert(a->rc == a->cc);

  b = 0;
  if (a->f) {
    b = new_full_hmatrix(a->rc, a->cc);
    copy_upper_amatrix(a->f, aunit, b->f);
  }
  else {
    assert(a->son != 0);
    assert(a->rsons == a->csons);

    sons = a->rsons;

    b = new_super_hmatrix(a->rc, a->cc, sons, sons);

    for (j = 0; j < sons; j++) {
      for (i = 0; i < j; i++) {
	b1 = clone_hmatrix(a->son[i + j * sons]);
	ref_hmatrix(b->son + i + j * sons, b1);
      }

      b1 = clone_upper_hmatrix(aunit, a->son[j + j * sons]);
      ref_hmatrix(b->son + j + j * sons, b1);

      for (i = j + 1; i < sons; i++) {
	b1 = new_rk_hmatrix(a->rc->son[i], a->cc->son[j], 0);
	ref_hmatrix(b->son + i + j * sons, b1);
      }
    }

    update_hmatrix(b);
  }

  return b;
}

/* ------------------------------------------------------------
 * Solve triangular systems
 * ------------------------------------------------------------ */

static void
lowersolve_n_hmatrix_avector(bool aunit, pchmatrix a, pavector xp)
{
  avector   tmp1, tmp2;
  pavector  xp1, xp2;
  uint      sons;
  uint      roff, roff2;
  uint      i, j;

  assert(a->rc == a->cc);
  assert(a->cc->size == xp->dim);

  if (a->f)
    triangularsolve_amatrix_avector(true, aunit, false, a->f, xp);
  else {
    assert(a->son != 0);
    assert(a->rsons == a->csons);

    sons = a->rsons;

    roff = 0;
    for (i = 0; i < sons; i++) {
      xp1 = init_sub_avector(&tmp1, xp, a->son[i]->rc->size, roff);

      lowersolve_n_hmatrix_avector(aunit, a->son[i + i * sons], xp1);

      roff2 = roff + a->son[i]->rc->size;
      for (j = i + 1; j < sons; j++) {
	xp2 = init_sub_avector(&tmp2, xp, a->son[j]->rc->size, roff2);

	fastaddeval_hmatrix_avector(-1.0, a->son[j + i * sons], xp1, xp2);

	uninit_avector(xp2);

	roff2 += a->son[j]->rc->size;
      }
      assert(roff2 == a->rc->size);

      uninit_avector(xp1);

      roff += a->son[i]->rc->size;
    }
    assert(roff == a->rc->size);
  }
}

static void
lowersolve_t_hmatrix_avector(bool aunit, pchmatrix a, pavector xp)
{
  avector   tmp1, tmp2;
  pavector  xp1, xp2;
  uint      sons;
  uint      roff, roff2;
  uint      i, j;

  assert(a->rc == a->cc);
  assert(a->cc->size == xp->dim);

  if (a->f)
    triangularsolve_amatrix_avector(true, aunit, true, a->f, xp);
  else {
    assert(a->son != 0);
    assert(a->rsons == a->csons);

    sons = a->rsons;

    roff = a->rc->size;
    for (i = sons; i-- > 0;) {
      roff -= a->son[i]->rc->size;

      xp1 = init_sub_avector(&tmp1, xp, a->son[i]->rc->size, roff);

      lowersolve_t_hmatrix_avector(aunit, a->son[i + i * sons], xp1);

      roff2 = roff;
      for (j = i; j-- > 0;) {
	roff2 -= a->son[j]->rc->size;

	xp2 = init_sub_avector(&tmp2, xp, a->son[j]->rc->size, roff2);

	fastaddevaltrans_hmatrix_avector(-1.0, a->son[i + j * sons], xp1,
					 xp2);

	uninit_avector(xp2);
      }
      assert(roff2 == 0);

      uninit_avector(xp1);
    }
    assert(roff == 0);
  }
}

static void
lowersolve_hmatrix_avector(bool aunit, bool atrans, pchmatrix a, pavector xp)
{
  if (atrans)
    lowersolve_t_hmatrix_avector(aunit, a, xp);
  else
    lowersolve_n_hmatrix_avector(aunit, a, xp);
}

static void
uppersolve_n_hmatrix_avector(bool aunit, pchmatrix a, pavector xp)
{
  avector   tmp1, tmp2;
  pavector  xp1, xp2;
  uint      sons;
  uint      roff, roff2;
  uint      i, j;

  assert(a->rc == a->cc);
  assert(a->cc->size == xp->dim);

  if (a->f)
    triangularsolve_amatrix_avector(false, aunit, false, a->f, xp);
  else {
    assert(a->son != 0);
    assert(a->rsons == a->csons);

    sons = a->rsons;

    roff = a->rc->size;
    for (i = sons; i-- > 0;) {
      roff -= a->son[i]->rc->size;

      xp1 = init_sub_avector(&tmp1, xp, a->son[i]->rc->size, roff);

      uppersolve_n_hmatrix_avector(aunit, a->son[i + i * sons], xp1);

      roff2 = roff;
      for (j = i; j-- > 0;) {
	roff2 -= a->son[j]->rc->size;

	xp2 = init_sub_avector(&tmp2, xp, a->son[j]->rc->size, roff2);

	fastaddeval_hmatrix_avector(-1.0, a->son[j + i * sons], xp1, xp2);

	uninit_avector(xp2);
      }
      assert(roff2 == 0);

      uninit_avector(xp1);
    }
    assert(roff == 0);
  }
}

static void
uppersolve_t_hmatrix_avector(bool aunit, pchmatrix a, pavector xp)
{
  avector   tmp1, tmp2;
  pavector  xp1, xp2;
  uint      sons;
  uint      roff, roff2;
  uint      i, j;

  assert(a->rc == a->cc);
  assert(a->cc->size == xp->dim);

  if (a->f)
    triangularsolve_amatrix_avector(false, aunit, true, a->f, xp);
  else {
    assert(a->son != 0);
    assert(a->rsons == a->csons);

    sons = a->rsons;

    roff = 0;
    for (i = 0; i < sons; i++) {
      xp1 = init_sub_avector(&tmp1, xp, a->son[i]->rc->size, roff);

      uppersolve_t_hmatrix_avector(aunit, a->son[i + i * sons], xp1);

      roff2 = roff + a->son[i]->rc->size;
      for (j = i + 1; j < sons; j++) {
	xp2 = init_sub_avector(&tmp2, xp, a->son[j]->rc->size, roff2);

	fastaddevaltrans_hmatrix_avector(-1.0, a->son[i + j * sons], xp1,
					 xp2);

	uninit_avector(xp2);

	roff2 += a->son[j]->rc->size;
      }
      assert(roff2 == a->rc->size);

      uninit_avector(xp1);

      roff += a->son[i]->rc->size;
    }
    assert(roff == a->rc->size);
  }
}

static void
uppersolve_hmatrix_avector(bool aunit, bool atrans, pchmatrix a, pavector xp)
{
  if (atrans)
    uppersolve_t_hmatrix_avector(aunit, a, xp);
  else
    uppersolve_n_hmatrix_avector(aunit, a, xp);
}

void
triangularinvmul_hmatrix_avector(bool alower, bool aunit, bool atrans,
				 pchmatrix a, pavector xp)
{
  if (alower)
    lowersolve_hmatrix_avector(aunit, atrans, a, xp);
  else
    uppersolve_hmatrix_avector(aunit, atrans, a, xp);
}

void
triangularsolve_hmatrix_avector(bool alower, bool aunit, bool atrans,
				pchmatrix a, pavector x)
{
  avector   tmp;
  pavector  xp;
  const uint *idx;
  uint      i, n;

  assert(x->dim == a->rc->size);

  n = a->rc->size;
  idx = a->rc->idx;

  xp = init_avector(&tmp, n);
  for (i = 0; i < n; i++)
    xp->v[i] = x->v[idx[i]];

  triangularinvmul_hmatrix_avector(alower, aunit, atrans, a, xp);

  for (i = 0; i < n; i++)
    x->v[idx[i]] = xp->v[i];
  uninit_avector(xp);
}

static void
lowersolve_nn_hmatrix_amatrix(bool aunit, pchmatrix a, pamatrix xp)
{
  amatrix   tmp1, tmp2;
  pamatrix  xp1, xp2;
  uint      sons;
  uint      roff, roff2;
  uint      i, j;

  assert(a->rc == a->cc);
  assert(a->cc->size == xp->rows);

  if (a->f)
    triangularsolve_amatrix(true, aunit, false, a->f, false, xp);
  else {
    assert(a->son != 0);
    assert(a->rsons == a->csons);

    sons = a->rsons;

    roff = 0;
    for (i = 0; i < sons; i++) {
      xp1 =
	init_sub_amatrix(&tmp1, xp, a->son[i]->rc->size, roff, xp->cols, 0);

      lowersolve_nn_hmatrix_amatrix(aunit, a->son[i + i * sons], xp1);

      roff2 = roff + a->son[i]->rc->size;
      for (j = i + 1; j < sons; j++) {
	xp2 =
	  init_sub_amatrix(&tmp2, xp, a->son[j]->rc->size, roff2, xp->cols,
			   0);

	addmul_hmatrix_amatrix_amatrix(-1.0, false, a->son[j + i * sons],
				       false, xp1, false, xp2);

	uninit_amatrix(xp2);

	roff2 += a->son[j]->rc->size;
      }
      assert(roff2 == a->rc->size);

      uninit_amatrix(xp1);

      roff += a->son[i]->rc->size;
    }
    assert(roff == a->rc->size);
  }
}

static void
lowersolve_nt_hmatrix_amatrix(bool aunit, pchmatrix a, pamatrix xp)
{
  amatrix   tmp1, tmp2;
  pamatrix  xp1, xp2;
  uint      sons;
  uint      roff, roff2;
  uint      i, j;

  assert(a->rc == a->cc);
  assert(a->cc->size == xp->cols);

  if (a->f) {
    triangularsolve_amatrix(true, aunit, false, a->f, true, xp);
  }
  else {
    assert(a->son != 0);
    assert(a->rsons == a->csons);

    sons = a->rsons;

    roff = 0;
    for (i = 0; i < sons; i++) {
      xp1 =
	init_sub_amatrix(&tmp1, xp, xp->rows, 0, a->son[i]->rc->size, roff);

      lowersolve_nt_hmatrix_amatrix(aunit, a->son[i + i * sons], xp1);

      roff2 = roff + a->son[i]->rc->size;
      for (j = i + 1; j < sons; j++) {
	xp2 = init_sub_amatrix(&tmp2, xp, xp->rows, 0, a->son[j]->rc->size,
			       roff2);

	addmul_hmatrix_amatrix_amatrix(-1.0, false, a->son[j + i * sons],
				       true, xp1, true, xp2);

	uninit_amatrix(xp2);

	roff2 += a->son[j]->rc->size;
      }
      assert(roff2 == a->rc->size);

      uninit_amatrix(xp1);

      roff += a->son[i]->rc->size;
    }
    assert(roff == a->rc->size);
  }
}

static void
lowersolve_tn_hmatrix_amatrix(bool aunit, pchmatrix a, pamatrix xp)
{
  amatrix   tmp1, tmp2;
  pamatrix  xp1, xp2;
  uint      sons;
  uint      roff, roff2;
  uint      i, j;

  assert(a->rc == a->cc);
  assert(a->cc->size == xp->rows);

  if (a->f)
    triangularsolve_amatrix(true, aunit, true, a->f, false, xp);
  else {
    assert(a->son != 0);
    assert(a->rsons == a->csons);

    sons = a->rsons;

    roff = a->rc->size;
    for (i = sons; i-- > 0;) {
      roff -= a->son[i]->rc->size;

      xp1 =
	init_sub_amatrix(&tmp1, xp, a->son[i]->rc->size, roff, xp->cols, 0);

      lowersolve_tn_hmatrix_amatrix(aunit, a->son[i + i * sons], xp1);

      roff2 = roff;
      for (j = i; j-- > 0;) {
	roff2 -= a->son[j]->rc->size;

	xp2 =
	  init_sub_amatrix(&tmp2, xp, a->son[j]->rc->size, roff2, xp->cols,
			   0);

	addmul_hmatrix_amatrix_amatrix(-1.0, true, a->son[i + j * sons],
				       false, xp1, false, xp2);

	uninit_amatrix(xp2);
      }
      assert(roff2 == 0);

      uninit_amatrix(xp1);
    }
    assert(roff == 0);
  }
}

static void
lowersolve_tt_hmatrix_amatrix(bool aunit, pchmatrix a, pamatrix xp)
{
  amatrix   tmp1, tmp2;
  pamatrix  xp1, xp2;
  uint      sons;
  uint      roff, roff2;
  uint      i, j;

  assert(a->rc == a->cc);
  assert(a->cc->size == xp->cols);

  if (a->f)
    triangularsolve_amatrix(true, aunit, true, a->f, true, xp);
  else {
    assert(a->son != 0);
    assert(a->rsons == a->csons);

    sons = a->rsons;

    roff = a->rc->size;
    for (i = sons; i-- > 0;) {
      roff -= a->son[i]->rc->size;

      xp1 =
	init_sub_amatrix(&tmp1, xp, xp->rows, 0, a->son[i]->rc->size, roff);

      lowersolve_tt_hmatrix_amatrix(aunit, a->son[i + i * sons], xp1);

      roff2 = roff;
      for (j = i; j-- > 0;) {
	roff2 -= a->son[j]->rc->size;

	xp2 = init_sub_amatrix(&tmp2, xp, xp->rows, 0, a->son[j]->rc->size,
			       roff2);

	addmul_hmatrix_amatrix_amatrix(-1.0, true, a->son[i + j * sons], true,
				       xp1, true, xp2);

	uninit_amatrix(xp2);
      }
      assert(roff2 == 0);

      uninit_amatrix(xp1);
    }
    assert(roff == 0);
  }
}

static void
lowersolve_hmatrix_amatrix(bool aunit, bool atrans, pchmatrix a,
			   bool xtrans, pamatrix xp)
{
  if (atrans) {
    if (xtrans)
      lowersolve_tt_hmatrix_amatrix(aunit, a, xp);
    else
      lowersolve_tn_hmatrix_amatrix(aunit, a, xp);
  }
  else {
    if (xtrans)
      lowersolve_nt_hmatrix_amatrix(aunit, a, xp);
    else
      lowersolve_nn_hmatrix_amatrix(aunit, a, xp);
  }
}

static void
uppersolve_nn_hmatrix_amatrix(bool aunit, pchmatrix a, pamatrix xp)
{
  amatrix   tmp1, tmp2;
  pamatrix  xp1, xp2;
  uint      sons;
  uint      roff, roff2;
  uint      i, j;

  assert(a->rc == a->cc);
  assert(a->cc->size == xp->rows);

  if (a->f)
    triangularsolve_amatrix(false, aunit, false, a->f, false, xp);
  else {
    assert(a->son != 0);
    assert(a->rsons == a->csons);

    sons = a->rsons;

    roff = a->rc->size;
    for (i = sons; i-- > 0;) {
      roff -= a->son[i]->rc->size;

      xp1 =
	init_sub_amatrix(&tmp1, xp, a->son[i]->rc->size, roff, xp->cols, 0);

      uppersolve_nn_hmatrix_amatrix(aunit, a->son[i + i * sons], xp1);

      roff2 = roff;
      for (j = i; j-- > 0;) {
	roff2 -= a->son[j]->rc->size;

	xp2 =
	  init_sub_amatrix(&tmp2, xp, a->son[j]->rc->size, roff2, xp->cols,
			   0);

	addmul_hmatrix_amatrix_amatrix(-1.0, false, a->son[j + i * sons],
				       false, xp1, false, xp2);

	uninit_amatrix(xp2);
      }
      assert(roff2 == 0);

      uninit_amatrix(xp1);
    }
    assert(roff == 0);
  }
}

static void
uppersolve_nt_hmatrix_amatrix(bool aunit, pchmatrix a, pamatrix xp)
{
  amatrix   tmp1, tmp2;
  pamatrix  xp1, xp2;
  uint      sons;
  uint      roff, roff2;
  uint      i, j;

  assert(a->rc == a->cc);
  assert(a->cc->size == xp->cols);

  if (a->f)
    triangularsolve_amatrix(false, aunit, false, a->f, true, xp);
  else {
    assert(a->son != 0);
    assert(a->rsons == a->csons);

    sons = a->rsons;

    roff = a->rc->size;
    for (i = sons; i-- > 0;) {
      roff -= a->son[i]->rc->size;

      xp1 =
	init_sub_amatrix(&tmp1, xp, xp->rows, 0, a->son[i]->rc->size, roff);

      uppersolve_nt_hmatrix_amatrix(aunit, a->son[i + i * sons], xp1);

      roff2 = roff;
      for (j = i; j-- > 0;) {
	roff2 -= a->son[j]->rc->size;

	xp2 = init_sub_amatrix(&tmp2, xp, xp->rows, 0, a->son[j]->rc->size,
			       roff2);

	addmul_hmatrix_amatrix_amatrix(-1.0, false, a->son[j + i * sons],
				       true, xp1, true, xp2);

	uninit_amatrix(xp2);
      }
      assert(roff2 == 0);

      uninit_amatrix(xp1);
    }
    assert(roff == 0);
  }
}

static void
uppersolve_tn_hmatrix_amatrix(bool aunit, pchmatrix a, pamatrix xp)
{
  amatrix   tmp1, tmp2;
  pamatrix  xp1, xp2;
  uint      sons;
  uint      roff, roff2;
  uint      i, j;

  assert(a->rc == a->cc);
  assert(a->cc->size == xp->rows);

  if (a->f)
    triangularsolve_amatrix(false, aunit, true, a->f, false, xp);
  else {
    assert(a->son != 0);
    assert(a->rsons == a->csons);

    sons = a->rsons;

    roff = 0;
    for (i = 0; i < sons; i++) {
      xp1 =
	init_sub_amatrix(&tmp1, xp, a->son[i]->rc->size, roff, xp->cols, 0);

      uppersolve_tn_hmatrix_amatrix(aunit, a->son[i + i * sons], xp1);

      roff2 = roff + a->son[i]->rc->size;
      for (j = i + 1; j < sons; j++) {
	xp2 =
	  init_sub_amatrix(&tmp2, xp, a->son[j]->rc->size, roff2, xp->cols,
			   0);

	addmul_hmatrix_amatrix_amatrix(-1.0, true, a->son[i + j * sons],
				       false, xp1, false, xp2);

	uninit_amatrix(xp2);

	roff2 += a->son[j]->rc->size;
      }
      assert(roff2 == a->rc->size);

      uninit_amatrix(xp1);

      roff += a->son[i]->rc->size;
    }
    assert(roff == a->rc->size);
  }
}

static void
uppersolve_tt_hmatrix_amatrix(bool aunit, pchmatrix a, pamatrix xp)
{
  amatrix   tmp1, tmp2;
  pamatrix  xp1, xp2;
  uint      sons;
  uint      roff, roff2;
  uint      i, j;

  assert(a->rc == a->cc);
  assert(a->cc->size == xp->cols);

  if (a->f)
    triangularsolve_amatrix(false, aunit, true, a->f, true, xp);
  else {
    assert(a->son != 0);
    assert(a->rsons == a->csons);

    sons = a->rsons;

    roff = 0;
    for (i = 0; i < sons; i++) {
      xp1 =
	init_sub_amatrix(&tmp1, xp, xp->rows, 0, a->son[i]->rc->size, roff);

      uppersolve_tt_hmatrix_amatrix(aunit, a->son[i + i * sons], xp1);

      roff2 = roff + a->son[i]->rc->size;
      for (j = i + 1; j < sons; j++) {
	xp2 = init_sub_amatrix(&tmp2, xp, xp->rows, 0, a->son[j]->rc->size,
			       roff2);

	addmul_hmatrix_amatrix_amatrix(-1.0, true, a->son[i + j * sons], true,
				       xp1, true, xp2);

	uninit_amatrix(xp2);

	roff2 += a->son[j]->rc->size;
      }
      assert(roff2 == a->rc->size);

      uninit_amatrix(xp1);

      roff += a->son[i]->rc->size;
    }
    assert(roff == a->rc->size);
  }
}

static void
uppersolve_hmatrix_amatrix(bool aunit, bool atrans, pchmatrix a,
			   bool xtrans, pamatrix xp)
{
  if (atrans) {
    if (xtrans)
      uppersolve_tt_hmatrix_amatrix(aunit, a, xp);
    else
      uppersolve_tn_hmatrix_amatrix(aunit, a, xp);
  }
  else {
    if (xtrans)
      uppersolve_nt_hmatrix_amatrix(aunit, a, xp);
    else
      uppersolve_nn_hmatrix_amatrix(aunit, a, xp);
  }
}

void
triangularinvmul_hmatrix_amatrix(bool alower, bool aunit, bool atrans,
				 pchmatrix a, bool xtrans, pamatrix xp)
{
  if (alower)
    lowersolve_hmatrix_amatrix(aunit, atrans, a, xtrans, xp);
  else
    uppersolve_hmatrix_amatrix(aunit, atrans, a, xtrans, xp);
}

static void
lowersolve_nn_hmatrix(bool aunit, pchmatrix a, pctruncmode tm,
		      real eps, phmatrix xp)
{
  phmatrix  atmp;
  uint      sons;
  uint      i, j, k;

  assert(a->rc == a->cc);
  assert(a->cc == xp->rc);

  if (xp->f)
    lowersolve_nn_hmatrix_amatrix(aunit, a, xp->f);
  else if (xp->r)
    lowersolve_nn_hmatrix_amatrix(aunit, a, &xp->r->A);
  else {
    if (a->f) {
      atmp = split_sub_amatrix(a->f, a->rc, a->cc, (xp->son[0]->rc != xp->rc),
			       (xp->son[0]->rc != xp->rc));

      assert(xp->rsons == atmp->rsons);

      sons = atmp->rsons;

      for (k = 0; k < xp->csons; k++)
	for (i = 0; i < sons; i++) {
	  lowersolve_nn_hmatrix(aunit, atmp->son[i + i * sons], tm, eps,
				xp->son[i + k * xp->rsons]);

	  for (j = i + 1; j < sons; j++)
	    addmul_hmatrix(-1.0, false, atmp->son[j + i * sons], false,
			   xp->son[i + k * xp->rsons], tm, eps,
			   xp->son[j + k * xp->rsons]);
	}

      del_hmatrix(atmp);
    }
    else {
      assert(a->rsons == a->csons);
      assert(xp->rsons == a->rsons);

      sons = a->rsons;

      for (k = 0; k < xp->csons; k++)
	for (i = 0; i < sons; i++) {
	  lowersolve_nn_hmatrix(aunit, a->son[i + i * sons], tm, eps,
				xp->son[i + k * xp->rsons]);

	  for (j = i + 1; j < sons; j++)
	    addmul_hmatrix(-1.0, false, a->son[j + i * sons], false,
			   xp->son[i + k * xp->rsons], tm, eps,
			   xp->son[j + k * xp->rsons]);
	}
    }
  }
}

static void
lowersolve_nt_hmatrix(bool aunit, pchmatrix a, pctruncmode tm,
		      real eps, phmatrix xp)
{
  phmatrix  atmp;
  uint      sons;
  uint      i, j, k;

  assert(a->rc == a->cc);
  assert(a->cc == xp->cc);

  if (xp->f)
    lowersolve_nt_hmatrix_amatrix(aunit, a, xp->f);
  else if (xp->r)
    lowersolve_nn_hmatrix_amatrix(aunit, a, &xp->r->B);
  else {
    if (a->f) {
      atmp = split_sub_amatrix(a->f, a->rc, a->cc, (xp->son[0]->cc != xp->cc),
			       (xp->son[0]->cc != xp->cc));

      assert(xp->csons == atmp->rsons);

      sons = atmp->rsons;

      for (k = 0; k < xp->rsons; k++)
	for (i = 0; i < sons; i++) {
	  lowersolve_nt_hmatrix(aunit, atmp->son[i + i * sons], tm, eps,
				xp->son[k + i * xp->rsons]);

	  for (j = i + 1; j < sons; j++)
	    addmul_hmatrix(-1.0, false, xp->son[k + i * xp->rsons], true,
			   atmp->son[j + i * sons], tm, eps,
			   xp->son[k + j * xp->rsons]);
	}

      del_hmatrix(atmp);
    }
    else {
      assert(a->rsons == a->csons);
      assert(xp->csons == a->rsons);

      sons = a->rsons;

      for (k = 0; k < xp->rsons; k++)
	for (i = 0; i < sons; i++) {
	  lowersolve_nt_hmatrix(aunit, a->son[i + i * sons], tm, eps,
				xp->son[k + i * xp->rsons]);

	  for (j = i + 1; j < sons; j++)
	    addmul_hmatrix(-1.0, false, xp->son[k + i * xp->rsons], true,
			   a->son[j + i * sons], tm, eps,
			   xp->son[k + j * xp->rsons]);
	}
    }
  }
}

static void
lowersolve_tn_hmatrix(bool aunit, pchmatrix a, pctruncmode tm,
		      real eps, phmatrix xp)
{
  phmatrix  atmp;
  uint      sons;
  uint      i, j, k;

  assert(a->rc == a->cc);
  assert(a->cc == xp->rc);

  if (xp->f)
    lowersolve_tn_hmatrix_amatrix(aunit, a, xp->f);
  else if (xp->r)
    lowersolve_tn_hmatrix_amatrix(aunit, a, &xp->r->A);
  else {
    if (a->f) {
      atmp = split_sub_amatrix(a->f, a->rc, a->cc, (xp->son[0]->rc != xp->rc),
			       (xp->son[0]->rc != xp->rc));

      assert(xp->rsons == atmp->rsons);

      sons = atmp->rsons;

      for (k = 0; k < xp->csons; k++)
	for (i = sons; i-- > 0;) {
	  lowersolve_tn_hmatrix(aunit, atmp->son[i + i * sons], tm, eps,
				xp->son[i + k * xp->rsons]);

	  for (j = i; j-- > 0;)
	    addmul_hmatrix(-1.0, true, atmp->son[i + j * sons], false,
			   xp->son[i + k * xp->rsons], tm, eps,
			   xp->son[j + k * xp->rsons]);
	}

      del_hmatrix(atmp);
    }
    else {
      assert(a->rsons == a->csons);
      assert(xp->rsons == a->rsons);

      sons = a->rsons;

      for (k = 0; k < xp->csons; k++)
	for (i = sons; i-- > 0;) {
	  lowersolve_tn_hmatrix(aunit, a->son[i + i * sons], tm, eps,
				xp->son[i + k * xp->rsons]);

	  for (j = i; j-- > 0;)
	    addmul_hmatrix(-1.0, true, a->son[i + j * sons], false,
			   xp->son[i + k * xp->rsons], tm, eps,
			   xp->son[j + k * xp->rsons]);
	}
    }
  }
}

static void
lowersolve_tt_hmatrix(bool aunit, pchmatrix a, pctruncmode tm,
		      real eps, phmatrix xp)
{
  phmatrix  atmp;
  uint      sons;
  uint      i, j, k;

  assert(a->rc == a->cc);
  assert(a->cc == xp->cc);

  if (xp->f)
    lowersolve_tt_hmatrix_amatrix(aunit, a, xp->f);
  else if (xp->r)
    lowersolve_tn_hmatrix_amatrix(aunit, a, &xp->r->B);
  else {
    if (a->f) {
      atmp = split_sub_amatrix(a->f, a->rc, a->cc, (xp->son[0]->cc != xp->cc),
			       (xp->son[0]->cc != xp->cc));

      assert(xp->csons == atmp->rsons);

      sons = atmp->rsons;

      for (k = 0; k < xp->rsons; k++)
	for (i = sons; i-- > 0;) {
	  lowersolve_tt_hmatrix(aunit, atmp->son[i + i * sons], tm, eps,
				xp->son[k + i * xp->rsons]);

	  for (j = i; j-- > 0;)
	    addmul_hmatrix(-1.0, false, xp->son[k + i * xp->rsons], false,
			   atmp->son[i + j * sons], tm, eps,
			   xp->son[k + j * xp->rsons]);
	}

      del_hmatrix(atmp);
    }
    else {
      assert(a->rsons == a->csons);
      assert(xp->csons == a->rsons);

      sons = a->rsons;

      for (k = 0; k < xp->rsons; k++)
	for (i = sons; i-- > 0;) {
	  lowersolve_tt_hmatrix(aunit, a->son[i + i * sons], tm, eps,
				xp->son[k + i * xp->rsons]);

	  for (j = i; j-- > 0;)
	    addmul_hmatrix(-1.0, false, xp->son[k + i * xp->rsons], false,
			   a->son[i + j * sons], tm, eps,
			   xp->son[k + j * xp->rsons]);
	}
    }
  }
}

static void
lowersolve_hmatrix(bool aunit, bool atrans, pchmatrix a,
		   pctruncmode tm, real eps, bool xtrans, phmatrix xp)
{
  if (atrans) {
    if (xtrans)
      lowersolve_tt_hmatrix(aunit, a, tm, eps, xp);
    else
      lowersolve_tn_hmatrix(aunit, a, tm, eps, xp);
  }
  else {
    if (xtrans)
      lowersolve_nt_hmatrix(aunit, a, tm, eps, xp);
    else
      lowersolve_nn_hmatrix(aunit, a, tm, eps, xp);
  }
}

static void
uppersolve_nn_hmatrix(bool aunit, pchmatrix a, pctruncmode tm,
		      real eps, phmatrix xp)
{
  phmatrix  atmp;
  uint      sons;
  uint      i, j, k;

  assert(a->rc == a->cc);
  assert(a->cc == xp->rc);

  if (xp->f)
    uppersolve_nn_hmatrix_amatrix(aunit, a, xp->f);
  else if (xp->r)
    uppersolve_nn_hmatrix_amatrix(aunit, a, &xp->r->A);
  else {
    if (a->f) {
      atmp = split_sub_amatrix(a->f, a->rc, a->cc, (xp->son[0]->rc != xp->rc),
			       (xp->son[0]->rc != xp->rc));

      assert(xp->rsons == atmp->rsons);

      sons = atmp->rsons;

      for (k = 0; k < xp->csons; k++)
	for (i = sons; i-- > 0;) {
	  uppersolve_nn_hmatrix(aunit, atmp->son[i + i * sons], tm, eps,
				xp->son[i + k * xp->rsons]);

	  for (j = i; j-- > 0;)
	    addmul_hmatrix(-1.0, false, atmp->son[j + i * sons], false,
			   xp->son[i + k * xp->rsons], tm, eps,
			   xp->son[j + k * xp->rsons]);
	}

      del_hmatrix(atmp);
    }
    else {
      assert(a->rsons == a->csons);
      assert(xp->rsons == a->rsons);

      sons = a->rsons;

      for (k = 0; k < xp->csons; k++)
	for (i = sons; i-- > 0;) {
	  uppersolve_nn_hmatrix(aunit, a->son[i + i * sons], tm, eps,
				xp->son[i + k * xp->rsons]);

	  for (j = i; j-- > 0;)
	    addmul_hmatrix(-1.0, false, a->son[j + i * sons], false,
			   xp->son[i + k * xp->rsons], tm, eps,
			   xp->son[j + k * xp->rsons]);
	}
    }
  }
}

static void
uppersolve_nt_hmatrix(bool aunit, pchmatrix a, pctruncmode tm,
		      real eps, phmatrix xp)
{
  phmatrix  atmp;
  uint      sons;
  uint      i, j, k;

  assert(a->rc == a->cc);
  assert(a->cc == xp->cc);

  if (xp->f)
    uppersolve_nt_hmatrix_amatrix(aunit, a, xp->f);
  else if (xp->r)
    uppersolve_nn_hmatrix_amatrix(aunit, a, &xp->r->B);
  else {
    if (a->f) {
      atmp = split_sub_amatrix(a->f, a->rc, a->cc, (xp->son[0]->cc != xp->cc),
			       (xp->son[0]->cc != xp->cc));

      assert(xp->csons == atmp->rsons);

      sons = atmp->rsons;

      for (k = 0; k < xp->rsons; k++)
	for (i = sons; i-- > 0;) {
	  uppersolve_nt_hmatrix(aunit, atmp->son[i + i * sons], tm, eps,
				xp->son[k + i * xp->rsons]);

	  for (j = i; j-- > 0;)
	    addmul_hmatrix(-1.0, false, xp->son[k + i * xp->rsons], true,
			   atmp->son[j + i * sons], tm, eps,
			   xp->son[k + j * xp->rsons]);
	}

      del_hmatrix(atmp);
    }
    else {
      assert(a->rsons == a->csons);
      assert(xp->csons == a->rsons);

      sons = a->rsons;

      for (k = 0; k < xp->rsons; k++)
	for (i = sons; i-- > 0;) {
	  uppersolve_nt_hmatrix(aunit, a->son[i + i * sons], tm, eps,
				xp->son[k + i * xp->rsons]);

	  for (j = i; j-- > 0;)
	    addmul_hmatrix(-1.0, false, xp->son[k + i * xp->rsons], true,
			   a->son[j + i * sons], tm, eps,
			   xp->son[k + j * xp->rsons]);
	}
    }
  }
}

static void
uppersolve_tn_hmatrix(bool aunit, pchmatrix a, pctruncmode tm,
		      real eps, phmatrix xp)
{
  phmatrix  atmp;
  uint      sons;
  uint      i, j, k;

  assert(a->rc == a->cc);
  assert(a->cc == xp->rc);

  if (xp->f)
    uppersolve_tn_hmatrix_amatrix(aunit, a, xp->f);
  else if (xp->r)
    uppersolve_tn_hmatrix_amatrix(aunit, a, &xp->r->A);
  else {
    if (a->f) {
      atmp = split_sub_amatrix(a->f, a->rc, a->cc, (xp->son[0]->rc != xp->rc),
			       (xp->son[0]->rc != xp->rc));

      assert(xp->rsons == atmp->rsons);

      sons = atmp->rsons;

      for (k = 0; k < xp->csons; k++)
	for (i = 0; i < sons; i++) {
	  uppersolve_tn_hmatrix(aunit, atmp->son[i + i * sons], tm, eps,
				xp->son[i + k * xp->rsons]);

	  for (j = i + 1; j < sons; j++)
	    addmul_hmatrix(-1.0, true, atmp->son[i + j * sons], false,
			   xp->son[i + k * xp->rsons], tm, eps,
			   xp->son[j + k * xp->rsons]);
	}

      del_hmatrix(atmp);
    }
    else {
      assert(a->rsons == a->csons);
      assert(xp->rsons == a->rsons);

      sons = a->rsons;

      for (k = 0; k < xp->csons; k++)
	for (i = 0; i < sons; i++) {
	  uppersolve_tn_hmatrix(aunit, a->son[i + i * sons], tm, eps,
				xp->son[i + k * xp->rsons]);

	  for (j = i + 1; j < sons; j++)
	    addmul_hmatrix(-1.0, true, a->son[i + j * sons], false,
			   xp->son[i + k * xp->rsons], tm, eps,
			   xp->son[j + k * xp->rsons]);
	}
    }
  }
}

static void
uppersolve_tt_hmatrix(bool aunit, pchmatrix a, pctruncmode tm,
		      real eps, phmatrix xp)
{
  phmatrix  atmp;
  uint      sons;
  uint      i, j, k;

  assert(a->rc == a->cc);
  assert(a->cc == xp->cc);

  if (xp->f)
    uppersolve_tt_hmatrix_amatrix(aunit, a, xp->f);
  else if (xp->r)
    uppersolve_tn_hmatrix_amatrix(aunit, a, &xp->r->B);
  else {
    if (a->f) {
      atmp = split_sub_amatrix(a->f, a->rc, a->cc, (xp->son[0]->cc != xp->cc),
			       (xp->son[0]->cc != xp->cc));

      assert(xp->csons == atmp->rsons);

      sons = atmp->rsons;

      for (k = 0; k < xp->rsons; k++)
	for (i = 0; i < sons; i++) {
	  uppersolve_tt_hmatrix(aunit, atmp->son[i + i * sons], tm, eps,
				xp->son[k + i * xp->rsons]);

	  for (j = i + 1; j < sons; j++)
	    addmul_hmatrix(-1.0, false, xp->son[k + i * xp->rsons], false,
			   atmp->son[i + j * sons], tm, eps,
			   xp->son[k + j * xp->rsons]);
	}

      del_hmatrix(atmp);
    }
    else {
      assert(a->rsons == a->csons);
      assert(xp->csons == a->rsons);

      sons = a->rsons;

      for (k = 0; k < xp->rsons; k++)
	for (i = 0; i < sons; i++) {
	  uppersolve_tt_hmatrix(aunit, a->son[i + i * sons], tm, eps,
				xp->son[k + i * xp->rsons]);

	  for (j = i + 1; j < sons; j++)
	    addmul_hmatrix(-1.0, false, xp->son[k + i * xp->rsons], false,
			   a->son[i + j * sons], tm, eps,
			   xp->son[k + j * xp->rsons]);
	}
    }
  }
}

static void
uppersolve_hmatrix(bool aunit, bool atrans, pchmatrix a,
		   pctruncmode tm, real eps, bool xtrans, phmatrix xp)
{
  if (atrans) {
    if (xtrans)
      uppersolve_tt_hmatrix(aunit, a, tm, eps, xp);
    else
      uppersolve_tn_hmatrix(aunit, a, tm, eps, xp);
  }
  else {
    if (xtrans)
      uppersolve_nt_hmatrix(aunit, a, tm, eps, xp);
    else
      uppersolve_nn_hmatrix(aunit, a, tm, eps, xp);
  }
}

void
triangularinvmul_hmatrix(bool alower, bool aunit, bool atrans, pchmatrix a,
			 pctruncmode tm, real eps, bool xtrans, phmatrix xp)
{
  if (alower)
    lowersolve_hmatrix(aunit, atrans, a, tm, eps, xtrans, xp);
  else
    uppersolve_hmatrix(aunit, atrans, a, tm, eps, xtrans, xp);
}

static void
lowersolve_nn_amatrix_hmatrix(bool aunit, pcamatrix a,
			      pctruncmode tm, real eps, phmatrix xp)
{
  hmatrix   tmp;

  init_hmatrix(&tmp, xp->rc, xp->rc);
  tmp.f = (pamatrix) a;

  lowersolve_nn_hmatrix(aunit, &tmp, tm, eps, xp);
}

static void
lowersolve_nt_amatrix_hmatrix(bool aunit, pcamatrix a,
			      pctruncmode tm, real eps, phmatrix xp)
{
  hmatrix   tmp;

  init_hmatrix(&tmp, xp->cc, xp->cc);
  tmp.f = (pamatrix) a;

  lowersolve_nt_hmatrix(aunit, &tmp, tm, eps, xp);
}

static void
lowersolve_tn_amatrix_hmatrix(bool aunit, pcamatrix a,
			      pctruncmode tm, real eps, phmatrix xp)
{
  hmatrix   tmp;

  init_hmatrix(&tmp, xp->rc, xp->rc);
  tmp.f = (pamatrix) a;

  lowersolve_tn_hmatrix(aunit, &tmp, tm, eps, xp);
}

static void
lowersolve_tt_amatrix_hmatrix(bool aunit, pcamatrix a,
			      pctruncmode tm, real eps, phmatrix xp)
{
  hmatrix   tmp;

  init_hmatrix(&tmp, xp->cc, xp->cc);
  tmp.f = (pamatrix) a;

  lowersolve_tt_hmatrix(aunit, &tmp, tm, eps, xp);
}

static void
lowersolve_amatrix_hmatrix(bool aunit, bool atrans, pcamatrix a,
			   pctruncmode tm, real eps, bool xtrans, phmatrix xp)
{
  if (atrans) {
    if (xtrans)
      lowersolve_tt_amatrix_hmatrix(aunit, a, tm, eps, xp);
    else
      lowersolve_tn_amatrix_hmatrix(aunit, a, tm, eps, xp);
  }
  else {
    if (xtrans)
      lowersolve_nt_amatrix_hmatrix(aunit, a, tm, eps, xp);
    else
      lowersolve_nn_amatrix_hmatrix(aunit, a, tm, eps, xp);
  }
}

static void
uppersolve_nn_amatrix_hmatrix(bool aunit, pcamatrix a,
			      pctruncmode tm, real eps, phmatrix xp)
{
  hmatrix   tmp;

  init_hmatrix(&tmp, xp->rc, xp->rc);
  tmp.f = (pamatrix) a;

  uppersolve_nn_hmatrix(aunit, &tmp, tm, eps, xp);
}

static void
uppersolve_nt_amatrix_hmatrix(bool aunit, pcamatrix a,
			      pctruncmode tm, real eps, phmatrix xp)
{
  hmatrix   tmp;

  init_hmatrix(&tmp, xp->cc, xp->cc);
  tmp.f = (pamatrix) a;

  uppersolve_nt_hmatrix(aunit, &tmp, tm, eps, xp);
}

static void
uppersolve_tn_amatrix_hmatrix(bool aunit, pcamatrix a,
			      pctruncmode tm, real eps, phmatrix xp)
{
  hmatrix   tmp;

  init_hmatrix(&tmp, xp->rc, xp->rc);
  tmp.f = (pamatrix) a;

  uppersolve_tn_hmatrix(aunit, &tmp, tm, eps, xp);
}

static void
uppersolve_tt_amatrix_hmatrix(bool aunit, pcamatrix a,
			      pctruncmode tm, real eps, phmatrix xp)
{
  hmatrix   tmp;

  init_hmatrix(&tmp, xp->cc, xp->cc);
  tmp.f = (pamatrix) a;

  uppersolve_tt_hmatrix(aunit, &tmp, tm, eps, xp);
}

static void
uppersolve_amatrix_hmatrix(bool aunit, bool atrans, pcamatrix a,
			   pctruncmode tm, real eps, bool xtrans, phmatrix xp)
{
  if (atrans) {
    if (xtrans)
      uppersolve_tt_amatrix_hmatrix(aunit, a, tm, eps, xp);
    else
      uppersolve_tn_amatrix_hmatrix(aunit, a, tm, eps, xp);
  }
  else {
    if (xtrans)
      uppersolve_nt_amatrix_hmatrix(aunit, a, tm, eps, xp);
    else
      uppersolve_nn_amatrix_hmatrix(aunit, a, tm, eps, xp);
  }
}

void
triangularinvmul_amatrix_hmatrix(bool alower, bool aunit, bool atrans,
				 pcamatrix a, pctruncmode tm, real eps,
				 bool xtrans, phmatrix xp)
{
  if (alower)
    lowersolve_amatrix_hmatrix(aunit, atrans, a, tm, eps, xtrans, xp);
  else
    uppersolve_amatrix_hmatrix(aunit, atrans, a, tm, eps, xtrans, xp);
}

/* ------------------------------------------------------------
 * Evaluate triangular matrices
 * ------------------------------------------------------------ */

static void
lowereval_n_hmatrix_avector(bool aunit, pchmatrix a, pavector xp)
{
  avector   tmp1, tmp2;
  pavector  xp1, xp2;
  uint      sons;
  uint      roff, roff2;
  uint      i, j;

  assert(a->rc == a->cc);
  assert(a->cc->size == xp->dim);

  if (a->f)
    triangulareval_amatrix_avector(true, aunit, false, a->f, xp);
  else {
    assert(a->son != 0);
    assert(a->rsons == a->csons);

    sons = a->rsons;

    roff = a->rc->size;
    for (i = sons; i-- > 0;) {
      roff -= a->son[i]->rc->size;

      xp1 = init_sub_avector(&tmp1, xp, a->son[i]->rc->size, roff);

      roff2 = roff + a->son[i]->rc->size;
      for (j = i + 1; j < sons; j++) {
	xp2 = init_sub_avector(&tmp2, xp, a->son[j]->rc->size, roff2);

	fastaddeval_hmatrix_avector(1.0, a->son[j + i * sons], xp1, xp2);

	uninit_avector(xp2);

	roff2 += a->son[j]->rc->size;
      }
      assert(roff2 == a->rc->size);

      lowereval_n_hmatrix_avector(aunit, a->son[i + i * sons], xp1);

      uninit_avector(xp1);
    }
    assert(roff == 0);
  }
}

static void
lowereval_t_hmatrix_avector(bool aunit, pchmatrix a, pavector xp)
{
  avector   tmp1, tmp2;
  pavector  xp1, xp2;
  uint      sons;
  uint      roff, roff2;
  uint      i, j;

  assert(a->rc == a->cc);
  assert(a->cc->size == xp->dim);

  if (a->f)
    triangulareval_amatrix_avector(true, aunit, true, a->f, xp);
  else {
    assert(a->son != 0);
    assert(a->rsons == a->csons);

    sons = a->rsons;

    roff = 0;
    for (i = 0; i < sons; i++) {
      xp1 = init_sub_avector(&tmp1, xp, a->son[i]->rc->size, roff);

      roff2 = roff;
      for (j = i; j-- > 0;) {
	roff2 -= a->son[j]->rc->size;

	xp2 = init_sub_avector(&tmp2, xp, a->son[j]->rc->size, roff2);

	fastaddevaltrans_hmatrix_avector(1.0, a->son[i + j * sons], xp1, xp2);

	uninit_avector(xp2);
      }
      assert(roff2 == 0);

      lowereval_t_hmatrix_avector(aunit, a->son[i + i * sons], xp1);

      uninit_avector(xp1);

      roff += a->son[i]->rc->size;
    }
    assert(roff == a->rc->size);
  }
}

static void
lowereval_hmatrix_avector(bool aunit, bool atrans, pchmatrix a, pavector xp)
{
  if (atrans)
    lowereval_t_hmatrix_avector(aunit, a, xp);
  else
    lowereval_n_hmatrix_avector(aunit, a, xp);
}

static void
uppereval_n_hmatrix_avector(bool aunit, pchmatrix a, pavector xp)
{
  avector   tmp1, tmp2;
  pavector  xp1, xp2;
  uint      sons;
  uint      roff, roff2;
  uint      i, j;

  assert(a->rc == a->cc);
  assert(a->cc->size == xp->dim);

  if (a->f)
    triangulareval_amatrix_avector(false, aunit, false, a->f, xp);
  else {
    assert(a->son != 0);
    assert(a->rsons == a->csons);

    sons = a->rsons;

    roff = 0;
    for (i = 0; i < sons; i++) {
      xp1 = init_sub_avector(&tmp1, xp, a->son[i]->rc->size, roff);

      roff2 = roff;
      for (j = i; j-- > 0;) {
	roff2 -= a->son[j]->rc->size;

	xp2 = init_sub_avector(&tmp2, xp, a->son[j]->rc->size, roff2);

	fastaddeval_hmatrix_avector(1.0, a->son[j + i * sons], xp1, xp2);

	uninit_avector(xp2);
      }
      assert(roff2 == 0);

      uppereval_n_hmatrix_avector(aunit, a->son[i + i * sons], xp1);

      uninit_avector(xp1);

      roff += a->son[i]->rc->size;
    }
    assert(roff == a->rc->size);
  }
}

static void
uppereval_t_hmatrix_avector(bool aunit, pchmatrix a, pavector xp)
{
  avector   tmp1, tmp2;
  pavector  xp1, xp2;
  uint      sons;
  uint      roff, roff2;
  uint      i, j;

  assert(a->rc == a->cc);
  assert(a->cc->size == xp->dim);

  if (a->f)
    triangulareval_amatrix_avector(false, aunit, true, a->f, xp);
  else {
    assert(a->son != 0);
    assert(a->rsons == a->csons);

    sons = a->rsons;

    roff = a->rc->size;
    for (i = sons; i-- > 0;) {
      roff -= a->son[i]->rc->size;

      xp1 = init_sub_avector(&tmp1, xp, a->son[i]->rc->size, roff);

      roff2 = roff + a->son[i]->rc->size;
      for (j = i + 1; j < sons; j++) {
	xp2 = init_sub_avector(&tmp2, xp, a->son[j]->rc->size, roff2);

	fastaddevaltrans_hmatrix_avector(1.0, a->son[i + j * sons], xp1, xp2);

	uninit_avector(xp2);

	roff2 += a->son[j]->rc->size;
      }
      assert(roff2 == a->rc->size);

      uppereval_t_hmatrix_avector(aunit, a->son[i + i * sons], xp1);

      uninit_avector(xp1);
    }
    assert(roff == 0);
  }
}

static void
uppereval_hmatrix_avector(bool aunit, bool atrans, pchmatrix a, pavector xp)
{
  if (atrans)
    uppereval_t_hmatrix_avector(aunit, a, xp);
  else
    uppereval_n_hmatrix_avector(aunit, a, xp);
}

void
triangularmul_hmatrix_avector(bool alower, bool aunit, bool atrans,
			      pchmatrix a, pavector xp)
{
  if (alower)
    lowereval_hmatrix_avector(aunit, atrans, a, xp);
  else
    uppereval_hmatrix_avector(aunit, atrans, a, xp);
}

void
triangulareval_hmatrix_avector(bool alower, bool aunit, bool atrans,
			       pchmatrix a, pavector x)
{
  avector   tmp;
  pavector  xp;
  const uint *idx;
  uint      i, n;

  assert(a->rc == a->cc);
  assert(x->dim == a->rc->size);

  n = a->rc->size;
  idx = a->rc->idx;

  xp = init_avector(&tmp, n);
  for (i = 0; i < n; i++)
    xp->v[i] = x->v[idx[i]];

  triangularmul_hmatrix_avector(alower, aunit, atrans, a, xp);

  for (i = 0; i < n; i++)
    x->v[idx[i]] = xp->v[i];
  uninit_avector(xp);
}

static void
lowereval_n_hmatrix_amatrix(bool aunit, pchmatrix a, bool xtrans, pamatrix xp)
{
  amatrix   tmp1, tmp2;
  pamatrix  xp1, xp2;
  uint      sons;
  uint      roff, roff2;
  uint      i, j;

  assert(a->rc == a->cc);
  if (xtrans)
    assert(a->cc->size == xp->cols);
  else
    assert(a->cc->size == xp->rows);

  if (a->f)
    triangulareval_amatrix(true, aunit, false, a->f, xtrans, xp);
  else {
    assert(a->son != 0);
    assert(a->rsons == a->csons);

    sons = a->rsons;

    roff = a->rc->size;
    for (i = sons; i-- > 0;) {
      roff -= a->son[i]->rc->size;

      xp1 = (xtrans ?
	     init_sub_amatrix(&tmp1, xp, xp->rows, 0, a->son[i]->rc->size,
			      roff) :
	     init_sub_amatrix(&tmp1, xp, a->son[i]->rc->size, roff, xp->cols,
			      0));

      roff2 = roff + a->son[i]->rc->size;
      for (j = i + 1; j < sons; j++) {
	xp2 = (xtrans ?
	       init_sub_amatrix(&tmp2, xp, xp->rows, 0, a->son[j]->rc->size,
				roff2) :
	       init_sub_amatrix(&tmp2, xp, a->son[j]->rc->size, roff2,
				xp->cols, 0));

	addmul_hmatrix_amatrix_amatrix(1.0, false, a->son[j + i * sons],
				       xtrans, xp1, xtrans, xp2);

	uninit_amatrix(xp2);

	roff2 += a->son[j]->rc->size;
      }
      assert(roff2 == a->rc->size);

      lowereval_n_hmatrix_amatrix(aunit, a->son[i + i * sons], xtrans, xp1);

      uninit_amatrix(xp1);
    }
    assert(roff == 0);
  }
}

static void
lowereval_t_hmatrix_amatrix(bool aunit, pchmatrix a, bool xtrans, pamatrix xp)
{
  amatrix   tmp1, tmp2;
  pamatrix  xp1, xp2;
  uint      sons;
  uint      roff, roff2;
  uint      i, j;

  assert(a->rc == a->cc);
  if (xtrans)
    assert(a->cc->size == xp->cols);
  else
    assert(a->cc->size == xp->rows);

  if (a->f)
    triangulareval_amatrix(true, aunit, true, a->f, xtrans, xp);
  else {
    assert(a->son != 0);
    assert(a->rsons == a->csons);

    sons = a->rsons;

    roff = 0;
    for (i = 0; i < sons; i++) {
      xp1 = (xtrans ?
	     init_sub_amatrix(&tmp1, xp, xp->rows, 0, a->son[i]->rc->size,
			      roff) :
	     init_sub_amatrix(&tmp1, xp, a->son[i]->rc->size, roff, xp->cols,
			      0));

      roff2 = roff;
      for (j = i; j-- > 0;) {
	roff2 -= a->son[j]->rc->size;

	xp2 = (xtrans ?
	       init_sub_amatrix(&tmp2, xp, xp->rows, 0, a->son[j]->rc->size,
				roff2) :
	       init_sub_amatrix(&tmp2, xp, a->son[j]->rc->size, roff2,
				xp->cols, 0));

	addmul_hmatrix_amatrix_amatrix(1.0, true, a->son[i + j * sons],
				       xtrans, xp1, xtrans, xp2);

	uninit_amatrix(xp2);
      }
      assert(roff2 == 0);

      lowereval_t_hmatrix_amatrix(aunit, a->son[i + i * sons], xtrans, xp1);

      uninit_amatrix(xp1);

      roff += a->son[i]->rc->size;
    }
    assert(roff == a->rc->size);
  }
}

static void
lowereval_hmatrix_amatrix(bool aunit, bool atrans, pchmatrix a,
			  bool xtrans, pamatrix xp)
{
  if (atrans)
    lowereval_t_hmatrix_amatrix(aunit, a, xtrans, xp);
  else
    lowereval_n_hmatrix_amatrix(aunit, a, xtrans, xp);
}

static void
uppereval_n_hmatrix_amatrix(bool aunit, pchmatrix a, bool xtrans, pamatrix xp)
{
  amatrix   tmp1, tmp2;
  pamatrix  xp1, xp2;
  uint      sons;
  uint      roff, roff2;
  uint      i, j;

  assert(a->rc == a->cc);
  if (xtrans)
    assert(a->cc->size == xp->cols);
  else
    assert(a->cc->size == xp->rows);

  if (a->f)
    triangulareval_amatrix(false, aunit, false, a->f, xtrans, xp);
  else {
    assert(a->son != 0);
    assert(a->rsons == a->csons);

    sons = a->rsons;

    roff = 0;
    for (i = 0; i < sons; i++) {
      xp1 = (xtrans ?
	     init_sub_amatrix(&tmp1, xp, xp->rows, 0, a->son[i]->rc->size,
			      roff) :
	     init_sub_amatrix(&tmp1, xp, a->son[i]->rc->size, roff, xp->cols,
			      0));

      roff2 = roff;
      for (j = i; j-- > 0;) {
	roff2 -= a->son[j]->rc->size;

	xp2 = (xtrans ?
	       init_sub_amatrix(&tmp2, xp, xp->rows, 0, a->son[j]->rc->size,
				roff2) :
	       init_sub_amatrix(&tmp2, xp, a->son[j]->rc->size, roff2,
				xp->cols, 0));

	addmul_hmatrix_amatrix_amatrix(1.0, false, a->son[j + i * sons],
				       xtrans, xp1, xtrans, xp2);

	uninit_amatrix(xp2);
      }
      assert(roff2 == 0);

      uppereval_n_hmatrix_amatrix(aunit, a->son[i + i * sons], xtrans, xp1);

      uninit_amatrix(xp1);

      roff += a->son[i]->rc->size;
    }
    assert(roff == a->rc->size);
  }
}

static void
uppereval_t_hmatrix_amatrix(bool aunit, pchmatrix a, bool xtrans, pamatrix xp)
{
  amatrix   tmp1, tmp2;
  pamatrix  xp1, xp2;
  uint      sons;
  uint      roff, roff2;
  uint      i, j;

  assert(a->rc == a->cc);
  if (xtrans)
    assert(a->cc->size == xp->cols);
  else
    assert(a->cc->size == xp->rows);

  if (a->f)
    triangulareval_amatrix(false, aunit, true, a->f, xtrans, xp);
  else {
    assert(a->son != 0);
    assert(a->rsons == a->csons);

    sons = a->rsons;

    roff = a->rc->size;
    for (i = sons; i-- > 0;) {
      roff -= a->son[i]->rc->size;

      xp1 = (xtrans ?
	     init_sub_amatrix(&tmp1, xp, xp->rows, 0, a->son[i]->rc->size,
			      roff) :
	     init_sub_amatrix(&tmp1, xp, a->son[i]->rc->size, roff, xp->cols,
			      0));

      roff2 = roff + a->son[i]->rc->size;
      for (j = i + 1; j < sons; j++) {
	xp2 = (xtrans ?
	       init_sub_amatrix(&tmp2, xp, xp->rows, 0, a->son[j]->rc->size,
				roff2) :
	       init_sub_amatrix(&tmp2, xp, a->son[j]->rc->size, roff2,
				xp->cols, 0));

	addmul_hmatrix_amatrix_amatrix(1.0, true, a->son[i + j * sons],
				       xtrans, xp1, xtrans, xp2);

	uninit_amatrix(xp2);

	roff2 += a->son[j]->rc->size;
      }
      assert(roff2 == a->rc->size);

      uppereval_t_hmatrix_amatrix(aunit, a->son[i + i * sons], xtrans, xp1);

      uninit_amatrix(xp1);
    }
    assert(roff == 0);
  }
}

static void
uppereval_hmatrix_amatrix(bool aunit, bool atrans, pchmatrix a,
			  bool xtrans, pamatrix xp)
{
  if (atrans)
    uppereval_t_hmatrix_amatrix(aunit, a, xtrans, xp);
  else
    uppereval_n_hmatrix_amatrix(aunit, a, xtrans, xp);
}

void
triangularmul_hmatrix_amatrix(bool alower, bool aunit, bool atrans,
			      pchmatrix a, bool xtrans, pamatrix xp)
{
  if (alower)
    lowereval_hmatrix_amatrix(aunit, atrans, a, xtrans, xp);
  else
    uppereval_hmatrix_amatrix(aunit, atrans, a, xtrans, xp);
}

static void
lowereval_nn_hmatrix(bool aunit, pchmatrix a, pctruncmode tm,
		     real eps, phmatrix xp)
{
  phmatrix  atmp;
  uint      sons;
  uint      i, j, k;

  assert(a->rc == a->cc);
  assert(a->cc == xp->rc);

  if (xp->f)
    lowereval_n_hmatrix_amatrix(aunit, a, false, xp->f);
  else if (xp->r)
    lowereval_n_hmatrix_amatrix(aunit, a, false, &xp->r->A);
  else {
    if (a->f) {
      atmp = split_sub_amatrix(a->f, a->rc, a->cc, (xp->son[0]->rc != xp->rc),
			       (xp->son[0]->rc != xp->rc));

      assert(xp->rsons == atmp->rsons);

      sons = atmp->rsons;

      for (k = 0; k < xp->csons; k++)
	for (i = sons; i-- > 0;) {
	  for (j = i + 1; j < sons; j++)
	    addmul_hmatrix(1.0, false, atmp->son[j + i * sons], false,
			   xp->son[i + k * xp->rsons], tm, eps,
			   xp->son[j + k * xp->rsons]);

	  lowereval_nn_hmatrix(aunit, atmp->son[i + i * sons], tm, eps,
			       xp->son[i + k * xp->rsons]);
	}

      del_hmatrix(atmp);
    }
    else {
      assert(a->rsons == a->csons);
      assert(xp->rsons == a->rsons);

      sons = a->rsons;

      for (k = 0; k < xp->csons; k++)
	for (i = sons; i-- > 0;) {
	  for (j = i + 1; j < sons; j++)
	    addmul_hmatrix(1.0, false, a->son[j + i * sons], false,
			   xp->son[i + k * xp->rsons], tm, eps,
			   xp->son[j + k * xp->rsons]);

	  lowereval_nn_hmatrix(aunit, a->son[i + i * sons], tm, eps,
			       xp->son[i + k * xp->rsons]);
	}
    }
  }
}

static void
lowereval_nt_hmatrix(bool aunit, pchmatrix a, pctruncmode tm,
		     real eps, phmatrix xp)
{
  phmatrix  atmp;
  uint      sons;
  uint      i, j, k;

  assert(a->rc == a->cc);
  assert(a->cc == xp->cc);

  if (xp->f)
    lowereval_n_hmatrix_amatrix(aunit, a, true, xp->f);
  else if (xp->r)
    lowereval_n_hmatrix_amatrix(aunit, a, false, &xp->r->B);
  else {
    if (a->f) {
      atmp = split_sub_amatrix(a->f, a->rc, a->cc, (xp->son[0]->cc != xp->cc),
			       (xp->son[0]->cc != xp->cc));

      assert(xp->csons == atmp->rsons);

      sons = atmp->rsons;

      for (k = 0; k < xp->rsons; k++)
	for (i = sons; i-- > 0;) {
	  for (j = i + 1; j < sons; j++)
	    addmul_hmatrix(1.0, false, xp->son[k + i * xp->rsons], true,
			   atmp->son[j + i * sons], tm, eps,
			   xp->son[k + j * xp->rsons]);

	  lowereval_nt_hmatrix(aunit, atmp->son[i + i * sons], tm, eps,
			       xp->son[k + i * xp->rsons]);
	}

      del_hmatrix(atmp);
    }
    else {
      assert(a->rsons == a->csons);
      assert(xp->rsons == a->rsons);

      sons = a->rsons;

      for (k = 0; k < xp->rsons; k++)
	for (i = sons; i-- > 0;) {
	  for (j = i + 1; j < sons; j++)
	    addmul_hmatrix(1.0, false, xp->son[k + i * xp->rsons], true,
			   a->son[j + i * sons], tm, eps,
			   xp->son[k + j * xp->rsons]);

	  lowereval_nt_hmatrix(aunit, a->son[i + i * sons], tm, eps,
			       xp->son[k + i * xp->rsons]);
	}
    }
  }
}

static void
lowereval_tn_hmatrix(bool aunit, pchmatrix a, pctruncmode tm,
		     real eps, phmatrix xp)
{
  phmatrix  atmp;
  uint      sons;
  uint      i, j, k;

  assert(a->rc == a->cc);
  assert(a->cc == xp->rc);

  if (xp->f)
    lowereval_t_hmatrix_amatrix(aunit, a, false, xp->f);
  else if (xp->r)
    lowereval_t_hmatrix_amatrix(aunit, a, false, &xp->r->A);
  else {
    if (a->f) {
      atmp = split_sub_amatrix(a->f, a->rc, a->cc, (xp->son[0]->rc != xp->rc),
			       (xp->son[0]->rc != xp->rc));

      assert(xp->rsons == atmp->rsons);

      sons = atmp->rsons;

      for (k = 0; k < xp->csons; k++)
	for (i = 0; i < sons; i++) {
	  for (j = i; j-- > 0;)
	    addmul_hmatrix(1.0, true, atmp->son[i + j * sons], false,
			   xp->son[i + k * xp->rsons], tm, eps,
			   xp->son[j + k * xp->rsons]);

	  lowereval_tn_hmatrix(aunit, atmp->son[i + i * sons], tm, eps,
			       xp->son[i + k * xp->rsons]);
	}

      del_hmatrix(atmp);
    }
    else {
      assert(a->rsons == a->csons);
      assert(xp->rsons == a->rsons);

      sons = a->rsons;

      for (k = 0; k < xp->csons; k++)
	for (i = 0; i < sons; i++) {
	  for (j = i; j-- > 0;)
	    addmul_hmatrix(1.0, true, a->son[i + j * sons], false,
			   xp->son[i + k * xp->rsons], tm, eps,
			   xp->son[j + k * xp->rsons]);

	  lowereval_tn_hmatrix(aunit, a->son[i + i * sons], tm, eps,
			       xp->son[i + k * xp->rsons]);
	}
    }
  }
}

static void
lowereval_tt_hmatrix(bool aunit, pchmatrix a, pctruncmode tm,
		     real eps, phmatrix xp)
{
  phmatrix  atmp;
  uint      sons;
  uint      i, j, k;

  assert(a->rc == a->cc);
  assert(a->cc == xp->cc);

  if (xp->f)
    lowereval_t_hmatrix_amatrix(aunit, a, true, xp->f);
  else if (xp->r)
    lowereval_t_hmatrix_amatrix(aunit, a, false, &xp->r->B);
  else {
    if (a->f) {
      atmp = split_sub_amatrix(a->f, a->rc, a->cc, (xp->son[0]->cc != xp->cc),
			       (xp->son[0]->cc != xp->cc));

      assert(xp->csons == atmp->rsons);

      sons = atmp->rsons;

      for (k = 0; k < xp->rsons; k++)
	for (i = 0; i < sons; i++) {
	  for (j = i; j-- > 0;)
	    addmul_hmatrix(1.0, false, xp->son[k + i * xp->rsons], false,
			   atmp->son[i + j * sons], tm, eps,
			   xp->son[k + j * xp->rsons]);

	  lowereval_tt_hmatrix(aunit, atmp->son[i + i * sons], tm, eps,
			       xp->son[k + i * xp->rsons]);
	}

      del_hmatrix(atmp);
    }
    else {
      assert(a->rsons == a->csons);
      assert(xp->csons == a->rsons);

      sons = a->rsons;

      for (k = 0; k < xp->rsons; k++)
	for (i = 0; i < sons; i++) {
	  for (j = i; j-- > 0;)
	    addmul_hmatrix(1.0, false, xp->son[k + i * xp->rsons], false,
			   a->son[i + j * sons], tm, eps,
			   xp->son[k + j * xp->rsons]);

	  lowereval_tt_hmatrix(aunit, a->son[i + i * sons], tm, eps,
			       xp->son[k + i * xp->rsons]);
	}
    }
  }
}

static void
lowereval_hmatrix(bool aunit, bool atrans, pchmatrix a,
		  pctruncmode tm, real eps, bool xtrans, phmatrix xp)
{
  if (atrans) {
    if (xtrans)
      lowereval_tt_hmatrix(aunit, a, tm, eps, xp);
    else
      lowereval_tn_hmatrix(aunit, a, tm, eps, xp);
  }
  else {
    if (xtrans)
      lowereval_nt_hmatrix(aunit, a, tm, eps, xp);
    else
      lowereval_nn_hmatrix(aunit, a, tm, eps, xp);
  }
}

static void
uppereval_nn_hmatrix(bool aunit, pchmatrix a, pctruncmode tm,
		     real eps, phmatrix xp)
{
  phmatrix  atmp;
  uint      sons;
  uint      i, j, k;

  assert(a->rc == a->cc);
  assert(a->cc == xp->rc);

  if (xp->f)
    uppereval_n_hmatrix_amatrix(aunit, a, false, xp->f);
  else if (xp->r)
    uppereval_n_hmatrix_amatrix(aunit, a, false, &xp->r->A);
  else {
    if (a->f) {
      atmp = split_sub_amatrix(a->f, a->rc, a->cc, (xp->son[0]->rc != xp->rc),
			       (xp->son[0]->rc != xp->rc));

      assert(xp->rsons == atmp->rsons);

      sons = atmp->rsons;

      for (k = 0; k < xp->csons; k++)
	for (i = 0; i < sons; i++) {
	  for (j = i; j-- > 0;)
	    addmul_hmatrix(1.0, false, atmp->son[j + i * sons], false,
			   xp->son[i + k * xp->rsons], tm, eps,
			   xp->son[j + k * xp->rsons]);

	  uppereval_nn_hmatrix(aunit, atmp->son[i + i * sons], tm, eps,
			       xp->son[i + k * xp->rsons]);
	}

      del_hmatrix(atmp);
    }
    else {
      assert(a->rsons == a->csons);
      assert(xp->rsons == a->rsons);

      sons = a->rsons;

      for (k = 0; k < xp->csons; k++)
	for (i = 0; i < sons; i++) {
	  for (j = i; j-- > 0;)
	    addmul_hmatrix(1.0, false, a->son[j + i * sons], false,
			   xp->son[i + k * xp->rsons], tm, eps,
			   xp->son[j + k * xp->rsons]);

	  uppereval_nn_hmatrix(aunit, a->son[i + i * sons], tm, eps,
			       xp->son[i + k * xp->rsons]);
	}
    }
  }
}

static void
uppereval_nt_hmatrix(bool aunit, pchmatrix a, pctruncmode tm,
		     real eps, phmatrix xp)
{
  phmatrix  atmp;
  uint      sons;
  uint      i, j, k;

  assert(a->rc == a->cc);
  assert(a->cc == xp->cc);

  if (xp->f)
    uppereval_n_hmatrix_amatrix(aunit, a, true, xp->f);
  else if (xp->r)
    uppereval_n_hmatrix_amatrix(aunit, a, false, &xp->r->B);
  else {
    if (a->f) {
      atmp = split_sub_amatrix(a->f, a->rc, a->cc, (xp->son[0]->cc != xp->cc),
			       (xp->son[0]->cc != xp->cc));

      assert(xp->csons == atmp->rsons);

      sons = atmp->rsons;

      for (k = 0; k < xp->rsons; k++)
	for (i = 0; i < sons; i++) {
	  for (j = i; j-- > 0;)
	    addmul_hmatrix(1.0, false, xp->son[k + i * xp->rsons], true,
			   atmp->son[j + i * sons], tm, eps,
			   xp->son[k + j * xp->rsons]);

	  uppereval_nt_hmatrix(aunit, atmp->son[i + i * sons], tm, eps,
			       xp->son[k + i * xp->rsons]);
	}

      del_hmatrix(atmp);
    }
    else {
      assert(a->rsons == a->csons);
      assert(xp->csons == a->rsons);

      sons = a->rsons;

      for (k = 0; k < xp->rsons; k++)
	for (i = 0; i < sons; i++) {
	  for (j = i; j-- > 0;)
	    addmul_hmatrix(1.0, false, xp->son[k + i * xp->rsons], true,
			   a->son[j + i * sons], tm, eps,
			   xp->son[k + j * xp->rsons]);

	  uppereval_nt_hmatrix(aunit, a->son[i + i * sons], tm, eps,
			       xp->son[k + i * xp->rsons]);
	}
    }
  }
}

static void
uppereval_tn_hmatrix(bool aunit, pchmatrix a, pctruncmode tm,
		     real eps, phmatrix xp)
{
  phmatrix  atmp;
  uint      sons;
  uint      i, j, k;

  assert(a->rc == a->cc);
  assert(a->cc == xp->rc);

  if (xp->f)
    uppereval_t_hmatrix_amatrix(aunit, a, false, xp->f);
  else if (xp->r)
    uppereval_t_hmatrix_amatrix(aunit, a, false, &xp->r->A);
  else {
    if (a->f) {
      atmp = split_sub_amatrix(a->f, a->rc, a->cc, (xp->son[0]->rc != xp->rc),
			       (xp->son[0]->rc != xp->rc));

      assert(xp->rsons == atmp->rsons);

      sons = atmp->rsons;

      for (k = 0; k < xp->csons; k++)
	for (i = sons; i-- > 0;) {
	  for (j = i + 1; j < sons; j++)
	    addmul_hmatrix(1.0, true, atmp->son[i + j * sons], false,
			   xp->son[i + k * xp->rsons], tm, eps,
			   xp->son[j + k * xp->rsons]);

	  uppereval_tn_hmatrix(aunit, atmp->son[i + i * sons], tm, eps,
			       xp->son[i + k * xp->rsons]);
	}

      del_hmatrix(atmp);
    }
    else {
      assert(a->rsons == a->csons);
      assert(xp->rsons == a->rsons);

      sons = a->rsons;

      for (k = 0; k < xp->csons; k++)
	for (i = sons; i-- > 0;) {
	  for (j = i + 1; j < sons; j++)
	    addmul_hmatrix(1.0, true, a->son[i + j * sons], false,
			   xp->son[i + k * xp->rsons], tm, eps,
			   xp->son[j + k * xp->rsons]);

	  uppereval_tn_hmatrix(aunit, a->son[i + i * sons], tm, eps,
			       xp->son[i + k * xp->rsons]);
	}
    }
  }
}

static void
uppereval_tt_hmatrix(bool aunit, pchmatrix a, pctruncmode tm,
		     real eps, phmatrix xp)
{
  phmatrix  atmp;
  uint      sons;
  uint      i, j, k;

  assert(a->rc == a->cc);
  assert(a->cc == xp->cc);

  if (xp->f)
    uppereval_t_hmatrix_amatrix(aunit, a, true, xp->f);
  else if (xp->r)
    uppereval_t_hmatrix_amatrix(aunit, a, false, &xp->r->B);
  else {
    if (a->f) {
      atmp = split_sub_amatrix(a->f, a->rc, a->cc, (xp->son[0]->cc != xp->cc),
			       (xp->son[0]->cc != xp->cc));

      assert(xp->csons == atmp->rsons);

      sons = atmp->rsons;

      for (k = 0; k < xp->rsons; k++)
	for (i = sons; i-- > 0;) {
	  for (j = i + 1; j < sons; j++)
	    addmul_hmatrix(1.0, false, xp->son[k + i * xp->rsons], false,
			   atmp->son[i + j * sons], tm, eps,
			   xp->son[k + j * xp->rsons]);

	  uppereval_tt_hmatrix(aunit, atmp->son[i + i * sons], tm, eps,
			       xp->son[k + i * xp->rsons]);
	}

      del_hmatrix(atmp);
    }
    else {
      assert(a->rsons == a->csons);
      assert(xp->csons == a->rsons);

      sons = a->rsons;

      for (k = 0; k < xp->rsons; k++)
	for (i = sons; i-- > 0;) {
	  for (j = i + 1; j < sons; j++)
	    addmul_hmatrix(1.0, false, xp->son[k + i * xp->rsons], false,
			   a->son[i + j * sons], tm, eps,
			   xp->son[k + j * xp->rsons]);

	  uppereval_tt_hmatrix(aunit, a->son[i + i * sons], tm, eps,
			       xp->son[k + i * xp->rsons]);
	}
    }
  }
}

static void
uppereval_hmatrix(bool aunit, bool atrans, pchmatrix a,
		  pctruncmode tm, real eps, bool xtrans, phmatrix xp)
{
  if (atrans) {
    if (xtrans)
      uppereval_tt_hmatrix(aunit, a, tm, eps, xp);
    else
      uppereval_tn_hmatrix(aunit, a, tm, eps, xp);
  }
  else {
    if (xtrans)
      uppereval_nt_hmatrix(aunit, a, tm, eps, xp);
    else
      uppereval_nn_hmatrix(aunit, a, tm, eps, xp);
  }
}

void
triangularmul_hmatrix(bool alower, bool aunit, bool atrans, pchmatrix a,
		      pctruncmode tm, real eps, bool xtrans, phmatrix xp)
{
  if (alower)
    lowereval_hmatrix(aunit, atrans, a, tm, eps, xtrans, xp);
  else
    uppereval_hmatrix(aunit, atrans, a, tm, eps, xtrans, xp);
}

/* ------------------------------------------------------------
 * Triangular factorizations
 * ------------------------------------------------------------ */

void
lrdecomp_hmatrix(phmatrix a, pctruncmode tm, real eps)
{
  uint      sons;
  uint      i, j, k;
  uint      res;

  assert(a->rc == a->cc);

  if (a->f) {
    res = lrdecomp_amatrix(a->f);
    assert(res == 0);
  }
  else {
    assert(a->son != 0);
    assert(a->rsons == a->csons);

    sons = a->rsons;

    for (k = 0; k < sons; k++) {
      lrdecomp_hmatrix(a->son[k + k * sons], tm, eps);

      for (j = k + 1; j < sons; j++)
	lowersolve_hmatrix(true, false, a->son[k + k * sons], tm, eps, false,
			   a->son[k + j * sons]);

      for (i = k + 1; i < sons; i++)
	uppersolve_hmatrix(false, true, a->son[k + k * sons], tm, eps, true,
			   a->son[i + k * sons]);

      for (j = k + 1; j < sons; j++)
	for (i = k + 1; i < sons; i++)
	  addmul_hmatrix(-1.0, false, a->son[i + k * sons], false,
			 a->son[k + j * sons], tm, eps, a->son[i + j * sons]);
    }
  }
}

void
lrsolve_n_hmatrix_avector(pchmatrix a, pavector x)
{
  avector   tmp;
  pavector  xp;
  const uint *idx;
  uint      i, n;

  assert(x->dim == a->rc->size);

  n = a->rc->size;
  idx = a->rc->idx;

  xp = init_avector(&tmp, n);
  for (i = 0; i < n; i++)
    xp->v[i] = x->v[idx[i]];

  lowersolve_hmatrix_avector(true, false, a, xp);
  uppersolve_hmatrix_avector(false, false, a, xp);

  for (i = 0; i < n; i++)
    x->v[idx[i]] = xp->v[i];
  uninit_avector(xp);
}

void
lrsolve_t_hmatrix_avector(pchmatrix a, pavector x)
{
  avector   tmp;
  pavector  xp;
  const uint *idx;
  uint      i, n;

  assert(x->dim == a->rc->size);

  n = a->rc->size;
  idx = a->rc->idx;

  xp = init_avector(&tmp, n);
  for (i = 0; i < n; i++)
    xp->v[i] = x->v[idx[i]];

  uppersolve_hmatrix_avector(false, true, a, xp);
  lowersolve_hmatrix_avector(true, true, a, xp);

  for (i = 0; i < n; i++)
    x->v[idx[i]] = xp->v[i];
  uninit_avector(xp);
}

void
lrsolve_hmatrix_avector(bool atrans, pchmatrix a, pavector x)
{
  if (atrans) {
    lrsolve_t_hmatrix_avector(a, x);
  }
  else {
    lrsolve_n_hmatrix_avector(a, x);
  }
}

void
lreval_n_hmatrix_avector(pchmatrix a, pavector x)
{
  avector   tmp;
  pavector  xp;
  const uint *idx;
  uint      i, n;

  assert(x->dim == a->rc->size);

  n = a->rc->size;
  idx = a->rc->idx;

  xp = init_avector(&tmp, n);
  for (i = 0; i < n; i++)
    xp->v[i] = x->v[idx[i]];

  uppereval_hmatrix_avector(false, false, a, xp);
  lowereval_hmatrix_avector(true, false, a, xp);

  for (i = 0; i < n; i++)
    x->v[idx[i]] = xp->v[i];
  uninit_avector(xp);
}

void
lreval_t_hmatrix_avector(pchmatrix a, pavector x)
{
  avector   tmp;
  pavector  xp;
  const uint *idx;
  uint      i, n;

  assert(x->dim == a->rc->size);

  n = a->rc->size;
  idx = a->rc->idx;

  xp = init_avector(&tmp, n);
  for (i = 0; i < n; i++)
    xp->v[i] = x->v[idx[i]];

  lowereval_hmatrix_avector(true, true, a, xp);
  uppereval_hmatrix_avector(false, true, a, xp);

  for (i = 0; i < n; i++)
    x->v[idx[i]] = xp->v[i];
  uninit_avector(xp);
}

void
lreval_hmatrix_avector(bool atrans, pchmatrix a, pavector x)
{
  if (atrans) {
    lreval_t_hmatrix_avector(a, x);
  }
  else {
    lreval_n_hmatrix_avector(a, x);
  }
}

void
choldecomp_hmatrix(phmatrix a, pctruncmode tm, real eps)
{
  uint      sons;
  uint      i, j, k;
  uint      res;

  assert(a->rc == a->cc);

  if (a->f) {
    res = choldecomp_amatrix(a->f);
    assert(res == 0);
  }
  else {
    assert(a->son != 0);
    assert(a->rsons == a->csons);

    sons = a->rsons;

    for (k = 0; k < sons; k++) {
      choldecomp_hmatrix(a->son[k + k * sons], tm, eps);

      for (i = k + 1; i < sons; i++)
	lowersolve_hmatrix(false, false, a->son[k + k * sons], tm, eps, true,
			   a->son[i + k * sons]);

      for (j = k + 1; j < sons; j++) {
	addmul_lower_hmatrix(-1.0, false, a->son[j + k * sons], true,
			     a->son[j + k * sons], tm, eps,
			     a->son[j + j * sons]);
	for (i = j + 1; i < sons; i++)
	  addmul_hmatrix(-1.0, false, a->son[i + k * sons], true,
			 a->son[j + k * sons], tm, eps, a->son[i + j * sons]);
      }
    }
  }
}

void
cholsolve_hmatrix_avector(pchmatrix a, pavector x)
{
  avector   tmp;
  pavector  xp;
  const uint *idx;
  uint      i, n;

  assert(x->dim == a->rc->size);

  n = a->rc->size;
  idx = a->rc->idx;

  xp = init_avector(&tmp, n);
  for (i = 0; i < n; i++)
    xp->v[i] = x->v[idx[i]];

  lowersolve_hmatrix_avector(false, false, a, xp);
  lowersolve_hmatrix_avector(false, true, a, xp);

  for (i = 0; i < n; i++)
    x->v[idx[i]] = xp->v[i];
  uninit_avector(xp);
}

void
choleval_hmatrix_avector(pchmatrix a, pavector x)
{
  avector   tmp;
  pavector  xp;
  const uint *idx;
  uint      i, n;

  assert(x->dim == a->rc->size);

  n = a->rc->size;
  idx = a->rc->idx;

  xp = init_avector(&tmp, n);
  for (i = 0; i < n; i++)
    xp->v[i] = x->v[idx[i]];

  lowereval_hmatrix_avector(false, true, a, xp);
  lowereval_hmatrix_avector(false, false, a, xp);

  for (i = 0; i < n; i++)
    x->v[idx[i]] = xp->v[i];
  uninit_avector(xp);
}
