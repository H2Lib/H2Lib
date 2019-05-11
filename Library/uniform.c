
/* ------------------------------------------------------------
   This is the file "uniform.c" of the H2Lib package.
   All rights reserved, Steffen Boerm 2009
   ------------------------------------------------------------ */

#include <math.h>

#include "uniform.h"
#include "basic.h"

/* ------------------------------------------------------------
   Constructors and destructors
   ------------------------------------------------------------ */

puniform
new_uniform(pclusterbasis rb, pclusterbasis cb)
{
  puniform  u;

  u = allocmem(sizeof(uniform));

  u->rb = u->cb = 0;
  u->rnext = u->rprev = u->cnext = u->cprev = 0;

  ref_row_uniform(u, rb);
  ref_col_uniform(u, cb);

  init_amatrix(&u->S, rb->k, cb->k);

  return u;
}

void
del_uniform(puniform u)
{
  assert(u != 0);

  uninit_amatrix(&u->S);

  unref_row_uniform(u);
  unref_col_uniform(u);

  freemem(u);
}

/* ------------------------------------------------------------
   Reference management
   ------------------------------------------------------------ */

void
ref_row_uniform(puniform u, pclusterbasis rb)
{
  if (u->rb == rb)
    return;

  if (u->rb)
    unref_row_uniform(u);

  if (rb == 0)
    return;

  ref_clusterbasis(&u->rb, rb);

  u->rnext = rb->rlist;
  u->rprev = 0;
  if (rb->rlist)
    rb->rlist->rprev = u;
  rb->rlist = u;
}

void
ref_col_uniform(puniform u, pclusterbasis cb)
{
  if (u->cb == cb)
    return;

  if (u->cb)
    unref_col_uniform(u);

  if (cb == 0)
    return;

  ref_clusterbasis(&u->cb, cb);

  u->cnext = cb->clist;
  u->cprev = 0;
  if (cb->clist)
    cb->clist->cprev = u;
  cb->clist = u;
}

void
unref_row_uniform(puniform u)
{
  if (u->rb == 0) {
    assert(u->rnext == 0);
    assert(u->rprev == 0);
    return;
  }

  if (u->rnext)
    u->rnext->rprev = u->rprev;

  if (u->rprev)
    u->rprev->rnext = u->rnext;

  if (u->rb->rlist == u)
    u->rb->rlist = u->rnext;

  unref_clusterbasis(u->rb);

  u->rb = 0;
}

void
unref_col_uniform(puniform u)
{
  if (u->cb == 0) {
    assert(u->cnext == 0);
    assert(u->cprev == 0);
    return;
  }

  if (u->cnext)
    u->cnext->cprev = u->cprev;

  if (u->cprev)
    u->cprev->cnext = u->cnext;

  if (u->cb->clist == u)
    u->cb->clist = u->cnext;

  unref_clusterbasis(u->cb);

  u->cb = 0;
}

/* ------------------------------------------------------------
   Statistics
   ------------------------------------------------------------ */

size_t
getsize_uniform(pcuniform u)
{
  size_t    sz;

  sz = sizeof(uniform);
  sz += getsize_heap_amatrix(&u->S);

  return sz;
}

/* ------------------------------------------------------------
   Arithmetic operations
   ------------------------------------------------------------ */

void
clear_uniform(puniform u)
{
  clear_amatrix(&u->S);
}

void
copy_uniform(bool trans, pcuniform src, puniform trg)
{
  if (trans) {
    assert(trg->rb == src->cb);
    assert(trg->cb == src->rb);

    copy_amatrix(true, &src->S, &trg->S);
  }
  else {
    assert(trg->rb == src->rb);
    assert(trg->cb == src->cb);

    copy_amatrix(false, &src->S, &trg->S);
  }
}

puniform
clone_uniform(pcuniform src)
{
  puniform  u;

  u = new_uniform(src->rb, src->cb);
  copy_uniform(false, src, u);

  return u;
}

void
scale_uniform(field alpha, puniform u)
{
  scale_amatrix(alpha, &u->S);
}

void
random_uniform(puniform u)
{
  random_amatrix(&u->S);
}

/* ------------------------------------------------------------
   Matrix-vector multiplication
   ------------------------------------------------------------ */

void
mvm_uniform_avector(field alpha, bool trans, pcuniform u,
		    pcavector x, pavector y)
{
  avector   tmp1, tmp2;
  pavector  xt, yt;

  if (trans) {
    xt = init_avector(&tmp1, u->rb->kbranch);
    yt = init_avector(&tmp2, u->cb->kbranch);

    compress_clusterbasis_avector(u->rb, x, xt);

    clear_avector(yt);

    mvm_amatrix_avector(alpha, true, &u->S, xt, yt);

    expand_clusterbasis_avector(u->cb, yt, y);

    uninit_avector(yt);
    uninit_avector(xt);
  }
  else {
    xt = init_avector(&tmp1, u->cb->kbranch);
    yt = init_avector(&tmp2, u->rb->kbranch);

    compress_clusterbasis_avector(u->cb, x, xt);

    clear_avector(yt);

    mvm_amatrix_avector(alpha, false, &u->S, xt, yt);

    expand_clusterbasis_avector(u->rb, yt, y);

    uninit_avector(yt);
    uninit_avector(xt);
  }
}

/* ------------------------------------------------------------
   Conversion operations
   ------------------------------------------------------------ */

void
add_projected_uniform(pcuniform u,
		      pcclusteroperator ro, pcclusteroperator co,
		      puniform unew)
{
  amatrix   tmp;
  pamatrix  XS;

  if (u->rb == unew->rb) {
    if (u->cb == unew->cb)
      add_amatrix(1.0, false, &u->S, &unew->S);
    else {
      assert(co->krow == unew->cb->k);
      assert(co->kcol == u->cb->k);

      addmul_amatrix(1.0, false, &u->S, true, &co->C, &unew->S);
    }
  }
  else {
    if (u->cb == unew->cb) {
      assert(ro->krow == unew->rb->k);
      assert(ro->kcol == u->rb->k);

      addmul_amatrix(1.0, false, &ro->C, false, &u->S, &unew->S);
    }
    else {
      assert(ro->krow == unew->rb->k);
      assert(ro->kcol == u->rb->k);
      assert(co->krow == unew->cb->k);
      assert(co->kcol == u->cb->k);

      XS = init_amatrix(&tmp, unew->rb->k, u->cb->k);
      clear_amatrix(XS);
      addmul_amatrix(1.0, false, &ro->C, false, &u->S, XS);
      addmul_amatrix(1.0, false, XS, true, &co->C, &unew->S);
      uninit_amatrix(XS);
    }
  }
}

void
project_inplace_uniform(puniform u,
			pclusterbasis rb, pcclusteroperator ro,
			pclusterbasis cb, pcclusteroperator co)
{
  amatrix   tmp;
  pamatrix  X;

  if (u->rb == rb) {
    if (u->cb == cb) {
      /* Nothing to do */
    }
    else {
      assert(co);
      assert(co->krow == cb->k);
      assert(co->kcol == u->cb->k);

      /* Transform coupling matrix */
      X = init_amatrix(&tmp, u->rb->k, cb->k);
      clear_amatrix(X);
      addmul_amatrix(1.0, false, &u->S, true, &co->C, X);

      resize_amatrix(&u->S, X->rows, X->cols);
      copy_amatrix(false, X, &u->S);

      uninit_amatrix(X);

      /* Update pointer */
      ref_col_uniform(u, cb);
    }
  }
  else {
    assert(ro);
    assert(ro->krow == rb->k);
    assert(ro->kcol == u->rb->k);

    if (u->cb == cb) {
      /* Transform coupling matrix */
      X = init_amatrix(&tmp, rb->k, u->cb->k);
      clear_amatrix(X);
      addmul_amatrix(1.0, false, &ro->C, false, &u->S, X);

      resize_amatrix(&u->S, X->rows, X->cols);
      copy_amatrix(false, X, &u->S);

      uninit_amatrix(X);

      /* Update pointer */
      ref_row_uniform(u, rb);
    }
    else {
      assert(co);
      assert(co->krow == cb->k);
      assert(co->kcol == u->cb->k);

      /* Transform coupling matrix for row... */
      X = init_amatrix(&tmp, rb->k, u->cb->k);
      clear_amatrix(X);
      addmul_amatrix(1.0, false, &ro->C, false, &u->S, X);

      /* ... and column basis */
      resize_amatrix(&u->S, rb->k, cb->k);
      clear_amatrix(&u->S);
      addmul_amatrix(1.0, false, X, true, &co->C, &u->S);

      uninit_amatrix(X);

      /* Update pointers */
      ref_row_uniform(u, rb);
      ref_col_uniform(u, cb);
    }
  }
}

void
add_rkmatrix_uniform(pcrkmatrix r, puniform unew)
{
  amatrix   tmp1, tmp2, tmp3, tmp4;
  pamatrix  At, Bt, Ac, Bc;
  uint      k;

  assert(r->A.cols == r->B.cols);

  k = r->A.cols;

  At = init_amatrix(&tmp1, unew->rb->kbranch, k);
  compress_clusterbasis_amatrix(unew->rb, &r->A, At);

  Bt = init_amatrix(&tmp2, unew->cb->kbranch, k);
  compress_clusterbasis_amatrix(unew->cb, &r->B, Bt);

  Ac = init_sub_amatrix(&tmp3, At, unew->rb->k, 0, k, 0);
  Bc = init_sub_amatrix(&tmp4, Bt, unew->cb->k, 0, k, 0);

  addmul_amatrix(1.0, false, Ac, true, Bc, &unew->S);

  uninit_amatrix(Bc);
  uninit_amatrix(Ac);
  uninit_amatrix(Bt);
  uninit_amatrix(At);
}
