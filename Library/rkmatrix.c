/* ------------------------------------------------------------
 This is the file "rkmatrix.c" of the H2Lib package.
 All rights reserved, Steffen Boerm 2010
 ------------------------------------------------------------ */

#include "rkmatrix.h"

#include "basic.h"

/* ------------------------------------------------------------
 Constructors and destructors
 ------------------------------------------------------------ */

prkmatrix
init_rkmatrix(prkmatrix r, uint rows, uint cols, uint k)
{
  init_amatrix(&r->A, rows, k);
  init_amatrix(&r->B, cols, k);
  r->k = k;

  return r;
}

pcrkmatrix
init_sub_rkmatrix(prkmatrix r, pcrkmatrix src, uint rows, uint roff,
		  uint cols, uint coff)
{
  prkmatrix wsrc = (prkmatrix) src;
  uint      k = src->k;

  init_sub_amatrix(&r->A, &wsrc->A, rows, roff, k, 0);
  init_sub_amatrix(&r->B, &wsrc->B, cols, coff, k, 0);
  r->k = k;

  return r;
}

void
uninit_rkmatrix(prkmatrix r)
{
  uninit_amatrix(&r->B);
  uninit_amatrix(&r->A);
}

prkmatrix
new_rkmatrix(uint rows, uint cols, uint k)
{
  prkmatrix r;

  r = (prkmatrix) allocmem(sizeof(rkmatrix));

  return init_rkmatrix(r, rows, cols, k);
}

pcrkmatrix
new_sub_rkmatrix(pcrkmatrix src, uint rows, uint roff, uint cols, uint coff)
{
  prkmatrix r;

  r = (prkmatrix) allocmem(sizeof(rkmatrix));

  return init_sub_rkmatrix(r, src, rows, roff, cols, coff);
}

void
del_rkmatrix(prkmatrix r)
{
  uninit_rkmatrix(r);

  freemem(r);
}

/* ------------------------------------------------------------
 Change rank
 ------------------------------------------------------------ */

void
setrank_rkmatrix(prkmatrix r, uint k)
{
  resize_amatrix(&r->A, r->A.rows, k);
  resize_amatrix(&r->B, r->B.rows, k);
  r->k = k;
}

void
resize_rkmatrix(prkmatrix r, uint rows, uint cols, uint k)
{
  resize_amatrix(&r->A, rows, k);
  resize_amatrix(&r->B, cols, k);
  r->k = k;
}

/* ------------------------------------------------------------
 Statistics
 ------------------------------------------------------------ */

size_t
getsize_rkmatrix(pcrkmatrix r)
{
  size_t    sz;

  sz = sizeof(rkmatrix);
  sz += getsize_heap_amatrix(&r->A);
  sz += getsize_heap_amatrix(&r->B);

  return sz;
}

size_t
getsize_heap_rkmatrix(pcrkmatrix r)
{
  size_t    sz;

  sz = getsize_heap_amatrix(&r->A);
  sz += getsize_heap_amatrix(&r->B);

  return sz;
}

/* ------------------------------------------------------------
 Simple utility functions
 ------------------------------------------------------------ */

prkmatrix
clone_rkmatrix(pcrkmatrix r)
{
  prkmatrix rnew;

  rnew = new_rkmatrix(r->A.rows, r->B.rows, r->k);

  copy_rkmatrix(false, r, rnew);

  return rnew;
}

void
copy_rkmatrix(bool atrans, pcrkmatrix a, prkmatrix b)
{
  if (atrans) {
    assert(a->A.rows == b->B.rows);
    assert(a->B.rows == b->A.rows);

    if (b->k != a->k)
      setrank_rkmatrix(b, a->k);

    copy_amatrix(false, &a->A, &b->B);
    copy_amatrix(false, &a->B, &b->A);
  }
  else {
    assert(a->A.rows == b->A.rows);
    assert(a->B.rows == b->B.rows);

    if (b->k != a->k)
      setrank_rkmatrix(b, a->k);

    copy_amatrix(false, &a->A, &b->A);
    copy_amatrix(false, &a->B, &b->B);
  }
}

void
scale_rkmatrix(field alpha, prkmatrix r)
{
  scale_amatrix(alpha, &r->A);
}

void
random_rkmatrix(prkmatrix r, uint kmax)
{
  setrank_rkmatrix(r, kmax);
  random_amatrix(&r->A);
  random_amatrix(&r->B);
}

/* ------------------------------------------------------------
 Matrix-vector multiplication
 ------------------------------------------------------------ */

void
addeval_rkmatrix_avector(field alpha, pcrkmatrix r, pcavector x, pavector y)
{
  pavector  ac, bc;
  avector   atmp, btmp;
  field     beta;
  uint      nu;

  assert(y->dim == r->A.rows);
  assert(x->dim == r->B.rows);
  assert(r->k <= r->A.cols);
  assert(r->k <= r->B.cols);

  for (nu = 0; nu < r->k; nu++) {
    ac = init_column_avector(&atmp, (pamatrix) &r->A, nu);
    bc = init_column_avector(&btmp, (pamatrix) &r->B, nu);

    beta = dotprod_avector(bc, x);
    add_avector(alpha * beta, ac, y);

    uninit_avector(bc);
    uninit_avector(ac);
  }
}

void
addevaltrans_rkmatrix_avector(field alpha, pcrkmatrix r, pcavector x,
			      pavector y)
{
  pavector  ac, bc;
  avector   atmp, btmp;
  field     beta;
  uint      nu;

  assert(x->dim == r->A.rows);
  assert(y->dim == r->B.rows);
  assert(r->k <= r->A.cols);
  assert(r->k <= r->B.cols);

  for (nu = 0; nu < r->k; nu++) {
    ac = init_column_avector(&atmp, (pamatrix) &r->A, nu);
    bc = init_column_avector(&btmp, (pamatrix) &r->B, nu);

    beta = dotprod_avector(ac, x);
    add_avector(alpha * beta, bc, y);

    uninit_avector(bc);
    uninit_avector(ac);
  }
}

void
mvm_rkmatrix_avector(field alpha, bool rtrans, pcrkmatrix r, pcavector x,
		     pavector y)
{
  if (rtrans)
    addevaltrans_rkmatrix_avector(alpha, r, x, y);
  else
    addeval_rkmatrix_avector(alpha, r, x, y);
}

real
norm2_rkmatrix(pcrkmatrix R)
{
  return norm2_matrix((mvm_t) mvm_rkmatrix_avector, (void *) R, R->A.rows,
		      R->B.rows);
}

real
norm2diff_rkmatrix(pcrkmatrix ra, pcrkmatrix rb)
{
  return norm2diff_matrix((mvm_t) mvm_rkmatrix_avector, (void *) ra,
			  (mvm_t) mvm_rkmatrix_avector, (void *) rb,
			  ra->A.rows, ra->B.rows);
}
