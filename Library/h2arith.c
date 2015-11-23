
/* ------------------------------------------------------------
 * This is the file "h2arith.c" of the H2Lib package.
 * All rights reserved, Knut Reimer 2012
 * ------------------------------------------------------------ */

#include <stdio.h>

#include "h2arith.h"

#include "basic.h"
#include "amatrix.h"
#include "factorizations.h"
#include "hmatrix.h"
#include "truncation.h"
#include "harith.h"
#include "clusteroperator.h"
#include "h2matrix.h"
#include "h2compression.h"
#include "h2update.h"

ph2matrix
build_from_block_lower_h2matrix(pcblock b, pclusterbasis rb, pclusterbasis cb)
{
  ph2matrix h, h1;
  pcblock   b1;
  pclusterbasis rb1, cb1;
  uint      rsons, csons;
  uint      i, j;

  h = NULL;

  if (b->son) {
    rsons = b->rsons;
    csons = b->csons;

    h = new_super_h2matrix(rb, cb, rsons, csons);

    for (j = 0; j < csons; j++) {
      for (i = j; i < rsons; i++) {
	b1 = b->son[i + j * rsons];

	rb1 = rb;
	if (b1->rc != b->rc) {
	  assert(rb->sons == rsons);
	  rb1 = rb->son[i];
	}

	cb1 = cb;
	if (b1->cc != b->cc) {
	  assert(cb->sons == csons);
	  cb1 = cb->son[j];
	}

	h1 = build_from_block_lower_h2matrix(b1, rb1, cb1);

	ref_h2matrix(h->son + i + j * rsons, h1);
      }
    }

    if (rb->t == cb->t) {
      for (j = 0; j < csons; j++) {
	for (i = 0; i < UINT_MIN(j, rsons); i++) {
	  b1 = b->son[i + j * rsons];

	  rb1 = rb;
	  if (b1->rc != b->rc) {
	    assert(rb->sons == rsons);
	    rb1 = rb->son[i];
	  }

	  cb1 = cb;
	  if (b1->cc != b->cc) {
	    assert(cb->sons == csons);
	    cb1 = cb->son[j];
	  }

	  h1 = new_zero_h2matrix(rb1, cb1);

	  ref_h2matrix(h->son + i + j * rsons, h1);
	}
      }
    }
    else {
      for (j = 0; j < csons; j++) {
	for (i = 0; i < UINT_MIN(j, rsons); i++) {
	  b1 = b->son[i + j * rsons];

	  rb1 = rb;
	  if (b1->rc != b->rc) {
	    assert(rb->sons == rsons);
	    rb1 = rb->son[i];
	  }

	  cb1 = cb;
	  if (b1->cc != b->cc) {
	    assert(cb->sons == csons);
	    cb1 = cb->son[j];
	  }

	  h1 = build_from_block_lower_h2matrix(b1, rb1, cb1);

	  ref_h2matrix(h->son + i + j * rsons, h1);
	}
      }
    }
  }
  else if (b->a > 0)
    h = new_zero_h2matrix(rb, cb);
  else {
    assert(b->a == 0);
    h = new_full_h2matrix(rb, cb);
    clear_amatrix(h->f);
  }

  update_h2matrix(h);

  return h;
}

ph2matrix
build_from_block_upper_h2matrix(pcblock b, pclusterbasis rb, pclusterbasis cb)
{
  ph2matrix h, h1;
  pcblock   b1;
  pclusterbasis rb1, cb1;
  uint      rsons, csons;
  uint      i, j;

  h = NULL;

  if (b->son) {
    rsons = b->rsons;
    csons = b->csons;

    h = new_super_h2matrix(rb, cb, rsons, csons);

    for (j = 0; j < csons; j++) {
      for (i = 0; i < UINT_MIN(j + 1, rsons); i++) {
	b1 = b->son[i + j * rsons];

	rb1 = rb;
	if (b1->rc != b->rc) {
	  assert(rb->sons == rsons);
	  rb1 = rb->son[i];
	}

	cb1 = cb;
	if (b1->cc != b->cc) {
	  assert(cb->sons == csons);
	  cb1 = cb->son[j];
	}

	h1 = build_from_block_upper_h2matrix(b1, rb1, cb1);

	ref_h2matrix(h->son + i + j * rsons, h1);
      }
    }

    if (rb->t == cb->t) {
      for (j = 0; j < csons; j++) {
	for (i = j + 1; i < rsons; i++) {
	  b1 = b->son[i + j * rsons];

	  rb1 = rb;
	  if (b1->rc != b->rc) {
	    assert(rb->sons == rsons);
	    rb1 = rb->son[i];
	  }

	  cb1 = cb;
	  if (b1->cc != b->cc) {
	    assert(cb->sons == csons);
	    cb1 = cb->son[j];
	  }

	  h1 = new_zero_h2matrix(rb1, cb1);

	  ref_h2matrix(h->son + i + j * rsons, h1);
	}
      }
    }
    else {
      for (j = 0; j < csons; j++) {
	for (i = j + 1; i < rsons; i++) {
	  b1 = b->son[i + j * rsons];

	  rb1 = rb;
	  if (b1->rc != b->rc) {
	    assert(rb->sons == rsons);
	    rb1 = rb->son[i];
	  }

	  cb1 = cb;
	  if (b1->cc != b->cc) {
	    assert(cb->sons == csons);
	    cb1 = cb->son[j];
	  }

	  h1 = build_from_block_upper_h2matrix(b1, rb1, cb1);

	  ref_h2matrix(h->son + i + j * rsons, h1);
	}
      }
    }

  }
  else if (b->a > 0)
    h = new_zero_h2matrix(rb, cb);
  else {
    assert(b->a == 0);
    h = new_full_h2matrix(rb, cb);
    clear_amatrix(h->f);
  }

  update_h2matrix(h);

  return h;
}

void
tolower_h2matrix(ph2matrix A)
{
  uint      sons = A->rb->sons;

  ph2matrix h2;
  uint      i, j;

  assert(A->rb->t == A->cb->t);

  for (i = 0; i < sons; i++) {
    tolower_h2matrix(A->son[i + i * sons]);
    for (j = i + 1; j < sons; j++) {
      h2 = new_zero_h2matrix(A->rb->son[i], A->cb->son[j]);
      ref_h2matrix(&A->son[i + j * sons], h2);
    }
  }
}

/* ------------------------------------------------------------
 * Sum
 * ------------------------------------------------------------ */

void
add_hmatrix_h2matrix(pchmatrix Gh, ph2matrix Gh2, pclusteroperator rwf,
		     pclusteroperator cwf, ptruncmode tm, real eps)
{
  uint      rsons = Gh->rsons;
  uint      csons = Gh->csons;

  pclusteroperator rw, cw;
  uint      i, j;

  assert(Gh2->rb->t == Gh->rc);
  assert(Gh2->cb->t == Gh->cc);

  if (Gh->son) {
    rw = identify_son_clusterweight_clusteroperator(rwf, Gh->rc);
    cw = identify_son_clusterweight_clusteroperator(cwf, Gh->cc);

    assert(Gh2->son);
    for (i = 0; i < rsons; i++) {
      for (j = 0; j < csons; j++) {
	add_hmatrix_h2matrix(Gh->son[i + j * rsons], Gh2->son[i + j * rsons],
			     rw, cw, tm, eps);
      }
    }
  }
  else if (Gh->f) {
    assert(Gh2->f);
    add_amatrix(1.0, false, Gh->f, Gh2->f);
  }
  else if (Gh->r) {
    rkupdate_h2matrix(Gh->r, Gh2, rwf, cwf, tm, eps);
  }
}

/* ------------------------------------------------------------
 * Convert functions
 * ------------------------------------------------------------ */

pamatrix
convert_h2matrix_amatrix(bool h2trans, pch2matrix h2)
{
  pcclusterbasis rb = (h2trans ? h2->cb : h2->rb);
  pcclusterbasis cb = (h2trans ? h2->rb : h2->cb);
  uint      rows = rb->t->size;
  uint      cols = cb->t->size;

  uint      i;
  pavector  e_i, xt, yt;
  avector   v;
  pamatrix  A;

  A = new_zero_amatrix(rows, cols);
  e_i = new_avector(cols);
  clear_avector(e_i);
  xt = new_coeffs_clusterbasis_avector(cb);
  yt = new_coeffs_clusterbasis_avector(rb);

  for (i = 0; i < cols; i++) {
    setentry_avector(e_i, i, 1.0);
    init_column_avector(&v, A, i);

    clear_avector(yt);
    forward_nopermutation_clusterbasis_avector(cb, e_i, xt);
    if (h2trans)
      fastaddevaltrans_h2matrix_avector(1.0, h2, xt, yt);
    else
      fastaddeval_h2matrix_avector(1.0, h2, xt, yt);
    backward_nopermutation_clusterbasis_avector(rb, yt, &v);

    uninit_avector(&v);
    setentry_avector(e_i, i, 0.0);
  }

  del_avector(yt);
  del_avector(xt);
  del_avector(e_i);

  return A;
}

void
fastaddmul_clusterbasis_amatrix(pcclusterbasis cb, pamatrix Yt, pamatrix Y)
{
  amatrix   loc1, loc2, loc3;
  pamatrix  Y1, Yt1, Yc, Yp;
  uint      i, yoff;

  assert(Y->rows == cb->t->size);
  assert(Yt->rows == cb->kbranch);
  assert(Y->cols == Yt->cols);

  Yc = init_sub_amatrix(&loc1, Yt, cb->k, 0, Yt->cols, 0);

  if (cb->sons > 0) {
    yoff = 0;
    for (i = 0; i < cb->sons; i++) {
      Yt1 = init_sub_amatrix(&loc2, Yt, cb->son[i]->k, cb->k, Yt->cols, 0);
      addmul_amatrix(1.0, false, &cb->son[i]->E, false, Yc, Yt1);
      uninit_amatrix(Yt1);

      Y1 = init_sub_amatrix(&loc2, Y, cb->t->son[i]->size, yoff, Y->cols, 0);
      Yt1 = init_sub_amatrix(&loc3, Yt, cb->son[i]->kbranch, cb->k, Yt->cols,
			     0);
      fastaddmul_clusterbasis_amatrix(cb->son[i], Yt1, Y1);
      uninit_amatrix(Yt1);
      uninit_amatrix(Y1);

      yoff += cb->t->son[i]->size;
    }
    assert(yoff == cb->t->size);
  }
  else {
    Yp = init_sub_amatrix(&loc2, Yt, cb->t->size, cb->k, Yt->cols, 0);
    addmul_amatrix(1.0, false, &cb->V, false, Yc, Yp);
    add_amatrix(1.0, false, Yp, Y);
    uninit_amatrix(Yp);
  }
  uninit_amatrix(Yc);
  clear_amatrix(Yt);
}

prkmatrix
convert_uniform_rkmatrix(bool utrans, pcuniform u)
{
  pclusterbasis rb = u->rb;
  pclusterbasis cb = u->cb;
  uint      rows = u->rb->t->size;
  uint      cols = u->cb->t->size;
  uint      krow = u->rb->k;
  uint      kcol = u->cb->k;

  prkmatrix r;

  pamatrix  Yt;
  amatrix   tmp;
  uint      k;

  if (krow <= kcol) {
    k = krow;
    r = new_rkmatrix(rows, cols, k);
    /* r->A = rb->V */
    Yt = init_identity_amatrix(&tmp, rb->kbranch, k);
    clear_amatrix(&r->A);
    fastaddmul_clusterbasis_amatrix(rb, Yt, &r->A);
    uninit_amatrix(Yt);
    /* r->B = cb->W * u->S^T */
    Yt = init_zero_amatrix(&tmp, cb->kbranch, k);
    copy_sub_amatrix(true, &u->S, Yt);
    clear_amatrix(&r->B);
    fastaddmul_clusterbasis_amatrix(cb, Yt, &r->B);
    uninit_amatrix(Yt);
  }
  else {
    k = kcol;
    r = new_rkmatrix(rows, cols, k);
    /* r->A = rb->V * u->S */
    Yt = init_zero_amatrix(&tmp, rb->kbranch, k);
    copy_sub_amatrix(false, &u->S, Yt);
    clear_amatrix(&r->A);
    fastaddmul_clusterbasis_amatrix(rb, Yt, &r->A);
    uninit_amatrix(Yt);
    /* r->B = cb->W */
    Yt = init_identity_amatrix(&tmp, cb->kbranch, k);
    clear_amatrix(&r->B);
    fastaddmul_clusterbasis_amatrix(cb, Yt, &r->B);
    uninit_amatrix(Yt);
  }

  if (utrans) {
    tmp = r->A;
    r->A = r->B;
    r->B = tmp;
  }

  return r;
}

phmatrix
convert_h2matrix_hmatrix(pch2matrix h2)
{
  phmatrix  hm, hm1;

  uint      i;

  if (h2->u) {
    hm = new_hmatrix(h2->rb->t, h2->cb->t);
    hm->r = convert_uniform_rkmatrix(false, h2->u);
  }
  else if (h2->f) {
    hm = new_full_hmatrix(h2->rb->t, h2->cb->t);
    copy_amatrix(false, h2->f, hm->f);
  }
  else {
    assert(h2->son);
    hm = new_super_hmatrix(h2->rb->t, h2->cb->t, h2->rsons, h2->csons);
    for (i = 0; i < h2->rsons * h2->csons; i++) {
      hm1 = convert_h2matrix_hmatrix(h2->son[i]);
      ref_hmatrix(hm->son + i, hm1);
    }
  }

  update_hmatrix(hm);

  return hm;
}

/* ------------------------------------------------------------
 * Multiplication
 * ------------------------------------------------------------ */

prkmatrix
mul_h2matrix_1_rkmatrix(pch2matrix A, pch2matrix B, real tol)
{
  uint      rows = A->rb->t->size;
  uint      s = A->cb->t->size;
  uint      cols = B->cb->t->size;
  uint      rsons = A->rsons;
  uint      ssons = A->csons;
  uint      csons = B->csons;

  prkmatrix R;

  pamatrix  tmp;
  phmatrix  htmp;
  prkmatrix rtmp;
  prkmatrix p;
  pccluster rc1, cc1;

  uint      i, j, l;

  assert(A->cb->t == B->rb->t);

  /*  first case: A is a uniform matrix  */
  if (A->u) {
    p = convert_uniform_rkmatrix(false, A->u);
    R = new_rkmatrix(rows, cols, p->k);
    copy_amatrix(false, &p->A, &R->A);
    clear_amatrix(&R->B);
    addmul_h2matrix_amatrix_amatrix(1.0, true, B, false, &p->B, &R->B);
    del_rkmatrix(p);
  }
  /*  second case: B is a uniform matrix  */
  else if (B->u) {
    p = convert_uniform_rkmatrix(false, B->u);
    R = new_rkmatrix(rows, cols, p->k);
    copy_amatrix(false, &p->B, &R->B);
    clear_amatrix(&R->A);
    addmul_h2matrix_amatrix_amatrix(1.0, false, A, false, &p->A, &R->A);
    del_rkmatrix(p);
  }
  /*  third case: A is an amatrix */
  else if (A->f) {
    if (rows <= s) {
      R = new_rkmatrix(rows, cols, rows);
      identity_amatrix(&R->A);
      clear_amatrix(&R->B);
      addmul_h2matrix_amatrix_amatrix(1.0, true, B, true, A->f, &R->B);
    }
    else {
      /* R = (prkmatrix) allocmem(sizeof(rkmatrix));
         R->k = s; */
      R = new_rkmatrix(rows, cols, s);
      copy_amatrix(false, A->f, &R->A);
      /*      init_sub_amatrix(&R->A, A->f, rows, 0, s, 0); */
      tmp = convert_h2matrix_amatrix(true, B);
      copy_amatrix(false, tmp, &R->B);
      /*      R->B = *tmp; */
      del_amatrix(tmp);
    }
  }
  /*  fourth case: B is an amatrix */
  else if (B->f) {
    if (cols <= s) {
      R = new_rkmatrix(rows, cols, cols);
      identity_amatrix(&R->B);
      clear_amatrix(&R->A);
      addmul_h2matrix_amatrix_amatrix(1.0, false, A, false, B->f, &R->A);
    }
    else {
      /* R = (prkmatrix) allocmem(sizeof(rkmatrix));
         R->k = s; */
      R = new_rkmatrix(rows, cols, s);
      tmp = convert_h2matrix_amatrix(false, A);
      /* R->A = *tmp; */
      copy_amatrix(false, tmp, &R->A);
      del_amatrix(tmp);

      /*      init_amatrix(&R->B, cols, s); */
      copy_amatrix(true, B->f, &R->B);
    }
  }
  /*  fifth case: A and B have sons */
  else if (A->son && B->son) {
    htmp = new_super_hmatrix(A->rb->t, B->cb->t, rsons, csons);

    for (i = 0; i < rsons; i++) {
      rc1 = (A->son[i]->rb == A->rb ? A->rb->t : A->rb->t->son[i]);

      for (j = 0; j < csons; j++) {
	cc1 = (B->son[j * ssons]->cb == B->cb ? B->cb->t : B->cb->t->son[j]);

	ref_hmatrix(htmp->son + i + j * rsons, new_hmatrix(rc1, cc1));
	/*      htmp->son[i+j*rsons] = new_hmatrix(rc1, cc1); */
	rtmp = mul_h2matrix_1_rkmatrix(A->son[i], B->son[j * ssons], tol);

	for (l = 1; l < ssons; l++) {
	  p = mul_h2matrix_1_rkmatrix(A->son[i + l * rsons],
				      B->son[l + j * ssons], tol);

	  add_rkmatrix(1.0, p, 0, tol, rtmp);

	  del_rkmatrix(p);
	}
	htmp->son[i + j * rsons]->r = rtmp;
      }
    }

    R = merge_hmatrix_rkmatrix(htmp, 0, tol);

    del_hmatrix(htmp);
  }
  else {
    R = new_rkmatrix(rows, cols, 0);
    clear_amatrix(&R->A);
    clear_amatrix(&R->B);
  }
  assert(R->A.owner == NULL);
  return R;
}

prkmatrix
mul_h2matrix_2_rkmatrix(pch2matrix A, pch2matrix B, real tol)
{
  uint      rows = A->rb->t->size;
  uint      s = A->cb->t->size;
  uint      cols = B->rb->t->size;
  uint      rsons = A->rsons;
  uint      ssons = A->csons;
  uint      csons = B->rsons;

  prkmatrix R;

  pamatrix  tmp;
  phmatrix  htmp;
  prkmatrix rtmp;
  prkmatrix p;
  pccluster rc1, cc1;

  uint      i, j, l;

  assert(A->cb->t == B->cb->t);

  /*  first case: A is a uniform matrix  */
  if (A->u) {
    p = convert_uniform_rkmatrix(false, A->u);
    R = new_rkmatrix(rows, cols, p->k);
    copy_amatrix(false, &p->A, &R->A);
    clear_amatrix(&R->B);
    addmul_h2matrix_amatrix_amatrix(1.0, false, B, false, &p->B, &R->B);
    del_rkmatrix(p);
  }
  /*  second case: B is a uniform matrix  */
  else if (B->u) {
    p = convert_uniform_rkmatrix(false, B->u);
    R = new_rkmatrix(rows, cols, p->k);
    copy_amatrix(false, &p->A, &R->B);
    clear_amatrix(&R->A);
    addmul_h2matrix_amatrix_amatrix(1.0, false, A, false, &p->B, &R->A);
    del_rkmatrix(p);
  }
  /*  third case: A is an amatrix */
  else if (A->f) {
    if (rows <= s) {
      R = new_rkmatrix(rows, cols, rows);
      identity_amatrix(&R->A);
      clear_amatrix(&R->B);
      addmul_h2matrix_amatrix_amatrix(1.0, false, B, true, A->f, &R->B);
    }
    else {
      /*      R = (prkmatrix) allocmem(sizeof(rkmatrix));
         R->k = s; */
      R = new_rkmatrix(rows, cols, s);
      copy_amatrix(false, A->f, &R->A);

      /*      init_sub_amatrix(&R->A, A->f, rows, 0, s, 0); */
      tmp = convert_h2matrix_amatrix(false, B);
      copy_amatrix(false, tmp, &R->B);
      /*      R->B = *tmp; */
      del_amatrix(tmp);
    }
  }
  /*  fourth case: B is an amatrix */
  else if (B->f) {
    if (cols <= s) {
      R = new_rkmatrix(rows, cols, cols);
      identity_amatrix(&R->B);
      clear_amatrix(&R->A);
      addmul_h2matrix_amatrix_amatrix(1.0, false, A, true, B->f, &R->A);
    }
    else {
      /*      R = (prkmatrix) allocmem(sizeof(rkmatrix));
         R->k = s; */
      R = new_rkmatrix(rows, cols, s);
      tmp = convert_h2matrix_amatrix(false, A);
      /*      R->A = *tmp; */
      copy_amatrix(false, tmp, &R->A);
      del_amatrix(tmp);
      /*      init_amatrix(&R->B, cols, s) */ ;
      copy_amatrix(false, B->f, &R->B);
    }
  }
  /*  fifth case: A and B have sons */
  else if (A->son && B->son) {
    htmp = new_super_hmatrix(A->rb->t, B->rb->t, rsons, csons);

    for (i = 0; i < rsons; i++) {
      rc1 = (A->son[i]->rb == A->rb ? A->rb->t : A->rb->t->son[i]);

      for (j = 0; j < csons; j++) {
	cc1 = (B->son[j]->rb == B->rb ? B->rb->t : B->rb->t->son[j]);

	/*      htmp->son[i+j*rsons] = new_hmatrix(rc1, cc1); */
	ref_hmatrix(htmp->son + i + j * rsons, new_hmatrix(rc1, cc1));
	rtmp = mul_h2matrix_2_rkmatrix(A->son[i], B->son[j], tol);

	for (l = 1; l < ssons; l++) {
	  p = mul_h2matrix_2_rkmatrix(A->son[i + l * rsons],
				      B->son[j + l * csons], tol);

	  add_rkmatrix(1.0, p, 0, tol, rtmp);

	  del_rkmatrix(p);
	}
	htmp->son[i + j * rsons]->r = rtmp;
      }
    }

    R = merge_hmatrix_rkmatrix(htmp, 0, tol);

    del_hmatrix(htmp);
  }
  else {
    R = new_rkmatrix(rows, cols, 0);
    clear_amatrix(&R->A);
    clear_amatrix(&R->B);
  }

  return R;
}

prkmatrix
mul_h2matrix_rkmatrix(pch2matrix A, bool btrans, pch2matrix B, real tol)
{
  prkmatrix R;
  if (!btrans)
    R = mul_h2matrix_1_rkmatrix(A, B, tol);
  else
    R = mul_h2matrix_2_rkmatrix(A, B, tol);

  return R;
}

void
addmul_1_h2matrix(field alpha, pch2matrix A, pch2matrix B, ph2matrix C,
		  pclusteroperator rwf, pclusteroperator cwf, ptruncmode tm,
		  real tol)
{
  uint      rows = A->rb->t->size;
  uint      s = A->cb->t->size;
  uint      cols = B->cb->t->size;
  uint      rsons = A->rsons;
  uint      ssons = A->csons;
  uint      csons = B->csons;

  prkmatrix R, p;
  pamatrix  X;
  pclusteroperator rw, cw;
  uint      i, j, l;

  assert(C->rb->t == A->rb->t);
  assert(C->cb->t == B->cb->t);
  assert(A->cb->t == B->rb->t);
  assert(rwf->son);
  assert(rwf->sons > 0);
  assert(cwf->son);
  assert(cwf->sons > 0);

  /*  first case: a is a uniform matrix  */
  if (A->u) {
    p = convert_uniform_rkmatrix(false, A->u);
    R = new_rkmatrix(rows, cols, p->k);
    copy_amatrix(false, &p->A, &R->A);
    clear_amatrix(&R->B);
    addmul_h2matrix_amatrix_amatrix(alpha, true, B, false, &p->B, &R->B);
    del_rkmatrix(p);

    trunc_rkmatrix(0, tol, R);
    rkupdate_h2matrix(R, C, rwf, cwf, tm, tol);

    del_rkmatrix(R);
  }
  /*  second case: b is a uniform matrix  */
  else if (B->u) {
    p = convert_uniform_rkmatrix(false, B->u);
    R = new_rkmatrix(rows, cols, p->k);
    copy_amatrix(false, &p->B, &R->B);
    clear_amatrix(&R->A);
    addmul_h2matrix_amatrix_amatrix(alpha, false, A, false, &p->A, &R->A);
    del_rkmatrix(p);

    trunc_rkmatrix(0, tol, R);
    rkupdate_h2matrix(R, C, rwf, cwf, tm, tol);

    del_rkmatrix(R);
  }
  /*  third case: a is an amatrix */
  else if (A->f) {
    if (rows <= s) {
      R = new_rkmatrix(rows, cols, rows);
      identity_amatrix(&R->A);
      clear_amatrix(&R->B);
      addmul_h2matrix_amatrix_amatrix(alpha, true, B, true, A->f, &R->B);
    }
    else {
      /*      R = (prkmatrix) allocmem(sizeof(rkmatrix));
         R->k = s; */
      R = new_rkmatrix(rows, cols, s);
      copy_amatrix(false, A->f, &R->A);
      /*      init_sub_amatrix(&R->A, A->f, rows, 0, s, 0); */
      X = convert_h2matrix_amatrix(true, B);
      copy_amatrix(false, X, &R->B);
      del_amatrix(X);
      /*      R->B = *X;
         freemem(X); */
      scale_amatrix(alpha, &R->B);
    }

    rkupdate_h2matrix(R, C, rwf, cwf, tm, tol);

    del_rkmatrix(R);
  }
  /*  fourth case: b is an amatrix */
  else if (B->f) {
    if (cols <= s && false) {
      R = new_rkmatrix(rows, cols, cols);
      identity_amatrix(&R->B);
      clear_amatrix(&R->A);
      addmul_h2matrix_amatrix_amatrix(alpha, false, A, false, B->f, &R->A);
    }
    else {
      /*      R = (prkmatrix) allocmem(sizeof(rkmatrix));
         R->k = s; */
      R = new_rkmatrix(rows, cols, s);
      X = convert_h2matrix_amatrix(false, A);
      /*      R->A = *X; */
      copy_amatrix(false, X, &R->A);
      /*      freemem(X); */
      del_amatrix(X);
      /*      init_amatrix(&R->B, cols, s); */
      copy_amatrix(true, B->f, &R->B);
      scale_amatrix(alpha, &R->B);
    }

    rkupdate_h2matrix(R, C, rwf, cwf, tm, tol);

    del_rkmatrix(R);
  }
  /*  fifth case: C is a uniform matrix  */
  else if (C->u) {
    R = mul_h2matrix_rkmatrix(A, false, B, tol);
    scale_amatrix(alpha, &R->A);
    rkupdate_h2matrix(R, C, rwf, cwf, tm, tol);

    del_rkmatrix(R);
  }
  /*  sixth case: C is an amatrix */
  else if (C->f) {
    X = convert_h2matrix_amatrix(false, B);
    addmul_h2matrix_amatrix_amatrix(alpha, false, A, false, X, C->f);
    del_amatrix(X);
  }
  /* seventh case: all three clusters have sons and c is not a rkmatrix  */
  else if (C->son && A->son && B->son) {
    rw = rwf;
    i = 0;
    while (rw->t != A->rb->t && i < rwf->sons) {
      rw = rwf->son[i];
      i++;
    }
    assert(rw->t == A->rb->t);
    cw = cwf;
    i = 0;
    while (cw->t != B->cb->t && i < cwf->sons) {
      cw = cwf->son[i];
      i++;
    }
    assert(cw->t == B->cb->t);

    if (rw->sons == 0)
      rw = rwf;
    if (cw->sons == 0)
      cw = cwf;

    for (i = 0; i < rsons; i++) {
      for (j = 0; j < csons; j++) {
	for (l = 0; l < ssons; l++) {
	  addmul_1_h2matrix(alpha, A->son[i + l * rsons],
			    B->son[l + j * ssons], C->son[i + j * rsons], rw,
			    cw, tm, tol);
	}
      }
    }
    /*orthogonal_block_h2matrix(C, rw, cw); */
    update_clusterbasis(C->rb);
    update_clusterbasis(C->cb);
  }
  /*eighth case: C is admissible zero block and A and B are not */
  else if (A->son && B->son) {
    C->u = new_uniform(C->rb, C->cb);
    clear_amatrix(&C->u->S);

    R = mul_h2matrix_rkmatrix(A, false, B, tol);
    scale_amatrix(alpha, &R->A);
    rkupdate_h2matrix(R, C, rwf, cwf, tm, tol);

    del_rkmatrix(R);
  }
}

void
addmul_2_h2matrix(field alpha, pch2matrix A, pch2matrix B, ph2matrix C,
		  pclusteroperator rwf, pclusteroperator cwf, ptruncmode tm,
		  real tol)
{
  uint      rows = A->rb->t->size;
  uint      s = A->cb->t->size;
  uint      cols = B->rb->t->size;
  uint      rsons = A->rsons;
  uint      ssons = A->csons;
  uint      csons = B->rsons;

  prkmatrix R, p;
  pamatrix  X;
  pclusteroperator rw, cw;
  uint      i, j, l;

  assert(C->rb->t == A->rb->t);
  assert(C->cb->t == B->rb->t);
  assert(A->cb->t == B->cb->t);
  assert(rwf->son);
  assert(rwf->sons > 0);
  assert(cwf->son);
  assert(cwf->sons > 0);

  /*  first case: a is a uniform matrix  */
  if (A->u) {
    p = convert_uniform_rkmatrix(false, A->u);
    R = new_rkmatrix(rows, cols, p->k);
    copy_amatrix(false, &p->A, &R->A);
    clear_amatrix(&R->B);
    addmul_h2matrix_amatrix_amatrix(alpha, false, B, false, &p->B, &R->B);
    del_rkmatrix(p);

    trunc_rkmatrix(0, tol, R);
    rkupdate_h2matrix(R, C, rwf, cwf, tm, tol);

    del_rkmatrix(R);
  }
  /*  second case: b is a uniform matrix  */
  else if (B->u) {
    p = convert_uniform_rkmatrix(false, B->u);
    R = new_rkmatrix(rows, cols, p->k);
    copy_amatrix(false, &p->A, &R->B);
    clear_amatrix(&R->A);
    addmul_h2matrix_amatrix_amatrix(alpha, false, A, false, &p->B, &R->A);
    del_rkmatrix(p);

    trunc_rkmatrix(0, tol, R);
    rkupdate_h2matrix(R, C, rwf, cwf, tm, tol);

    del_rkmatrix(R);
  }
  /*  third case: a is an amatrix */
  else if (A->f) {
    if (rows <= s) {
      R = new_rkmatrix(rows, cols, rows);
      identity_amatrix(&R->A);
      clear_amatrix(&R->B);
      addmul_h2matrix_amatrix_amatrix(alpha, false, B, true, A->f, &R->B);
    }
    else {
      /*       R = (prkmatrix) allocmem(sizeof(rkmatrix));
         R->k = s; */
      R = new_rkmatrix(rows, cols, s);
      copy_amatrix(false, A->f, &R->A);
      /*     init_sub_amatrix(&R->A, A->f, rows, 0, s, 0); */
      X = convert_h2matrix_amatrix(false, B);
      copy_amatrix(false, X, &R->B);
      /*      R->B = *X;
         freemem(X); */
      del_amatrix(X);
      scale_amatrix(alpha, &R->B);
    }

    rkupdate_h2matrix(R, C, rwf, cwf, tm, tol);

    del_rkmatrix(R);
  }
  /*  fourth case: b is an amatrix */
  else if (B->f) {
    if (cols <= s && false) {
      R = new_rkmatrix(rows, cols, cols);
      identity_amatrix(&R->B);
      clear_amatrix(&R->A);
      addmul_h2matrix_amatrix_amatrix(alpha, false, A, true, B->f, &R->A);
    }
    else {
      /*       R = (prkmatrix) allocmem(sizeof(rkmatrix));
         R->k = s; */
      R = new_rkmatrix(rows, cols, s);

      X = convert_h2matrix_amatrix(false, A);
      /*      R->A = *X; */
      copy_amatrix(false, X, &R->A);
      /*      freemem(X); */
      del_amatrix(X);
      /*   init_amatrix(&R->B, cols, s); */
      copy_amatrix(false, B->f, &R->B);
      scale_amatrix(alpha, &R->B);
    }

    rkupdate_h2matrix(R, C, rwf, cwf, tm, tol);

    del_rkmatrix(R);
  }
  /*  fifth case: C is a uniform matrix  */
  else if (C->u) {
    R = mul_h2matrix_rkmatrix(A, true, B, tol);
    scale_amatrix(alpha, &R->A);
    rkupdate_h2matrix(R, C, rwf, cwf, tm, tol);

    del_rkmatrix(R);
  }
  /*  sixth case: C is an amatrix */
  else if (C->f) {
    X = convert_h2matrix_amatrix(true, B);
    addmul_h2matrix_amatrix_amatrix(alpha, false, A, false, X, C->f);
    del_amatrix(X);
  }
  /* seventh case: all three clusters have sons and c is not a rkmatrix  */
  else if (C->son && A->son && B->son) {
    rw = identify_son_clusterweight_clusteroperator(rwf, C->rb->t);
    cw = identify_son_clusterweight_clusteroperator(cwf, C->cb->t);

    for (i = 0; i < rsons; i++) {
      for (j = 0; j < csons; j++) {
	for (l = 0; l < ssons; l++) {
	  addmul_2_h2matrix(alpha, A->son[i + l * rsons],
			    B->son[j + l * csons], C->son[i + j * rsons], rw,
			    cw, tm, tol);
	}
      }
    }
    /*orthogonal_block_h2matrix(C, rw, cw); */
    update_clusterbasis(C->rb);
    update_clusterbasis(C->cb);
  }
  /*eighth case: C is admissible zero block and A and B are not */
  else if (A->son && B->son) {
    C->u = new_uniform(C->rb, C->cb);
    clear_amatrix(&C->u->S);

    R = mul_h2matrix_rkmatrix(A, true, B, tol);
    scale_amatrix(alpha, &R->A);
    rkupdate_h2matrix(R, C, rwf, cwf, tm, tol);

    del_rkmatrix(R);
  }
}

void
addmul_h2matrix(field alpha, pch2matrix A, bool btrans, pch2matrix B,
		ph2matrix C, pclusteroperator rwf, pclusteroperator cwf,
		ptruncmode tm, real tol)
{
  if (!btrans)
    addmul_1_h2matrix(alpha, A, B, C, rwf, cwf, tm, tol);
  else
    addmul_2_h2matrix(alpha, A, B, C, rwf, cwf, tm, tol);
}

ph2matrix
mul_h2matrix(field alpha, ph2matrix A, ph2matrix B, ph2matrix h2,
	     ptruncmode tm, real tol)
{
  pccluster rc = A->rb->t;
  pccluster cc = B->cb->t;

  pclusterbasis rb, cb;
  ph2matrix C;
  pclusteroperator rwf, cwf;

  rb = build_from_cluster_clusterbasis(rc);
  cb = build_from_cluster_clusterbasis(cc);
  C = clonestructure_h2matrix(h2, rb, cb);

  rwf = prepare_row_clusteroperator(C->rb, C->cb, tm);
  cwf = prepare_col_clusteroperator(C->rb, C->cb, tm);

  addmul_h2matrix(alpha, A, false, B, C, rwf, cwf, tm, tol);

  del_clusteroperator(rwf);
  del_clusteroperator(cwf);
  update_tree_clusterbasis(C->rb);
  update_tree_clusterbasis(C->cb);
  return C;
}

/* ------------------------------------------------------------
 Invert an H2-matrix
 ------------------------------------------------------------ */

void
invert_h2matrix(ph2matrix h2, pclusteroperator rwf, pclusteroperator cwf,
		ph2matrix inv, pclusteroperator rwfinv,
		pclusteroperator cwfinv, ptruncmode tm, real tol)
{

  uint      i, j, l;

  pccluster t = h2->rb->t;
  uint      sons = h2->rsons;
  pclusteroperator rw, cw, rwinv, cwinv;
  ph2matrix *work;

  assert(t == h2->cb->t);

  if (h2->son) {
    assert(h2->csons == sons);
    /* select weights for current block */
    rw = identify_son_clusterweight_clusteroperator(rwf, t);
    rwinv = identify_son_clusterweight_clusteroperator(rwfinv, t);
    cw = identify_son_clusterweight_clusteroperator(cwf, t);
    cwinv = identify_son_clusterweight_clusteroperator(cwfinv, t);

    work = (ph2matrix *) allocmem((sons * sons) * sizeof(ph2matrix));
    for (l = 0; l < sons - 1; l++) {
      invert_h2matrix(h2->son[l + l * sons], rw, cw, inv->son[l + l * sons],
		      rwinv, cwinv, tm, tol);

      for (i = l + 1; i < sons; i++) {
	work[i + l * sons] = mul_h2matrix(1.0, h2->son[i + l * sons],
					  inv->son[l + l * sons],
					  h2->son[i + l * sons], tm, tol);
	work[l + i * sons] =
	  mul_h2matrix(1.0, inv->son[l + l * sons], h2->son[l + i * sons],
		       h2->son[l + i * sons], tm, tol);
	for (j = l + 1; j < sons; j++)
	  addmul_h2matrix(-1.0, work[i + l * sons], false,
			  h2->son[l + j * sons], h2->son[i + j * sons], rw,
			  cw, tm, tol);
      }
    }

    invert_h2matrix(h2->son[sons * sons - 1], rw, cw,
		    inv->son[sons * sons - 1], rwinv, cwinv, tm, tol);

    for (l = sons - 1; l > 0; l--) {
      for (i = l; i < sons; i++) {
	for (j = l; j < sons; j++) {
	  addmul_h2matrix(-1.0, inv->son[i + j * sons], false,
			  work[j + (l - 1) * sons],
			  inv->son[i + (l - 1) * sons], rwinv, cwinv, tm,
			  tol);
	  addmul_h2matrix(-1.0, work[(l - 1) + j * sons], false,
			  inv->son[j + i * sons],
			  inv->son[(l - 1) + i * sons], rwinv, cwinv, tm,
			  tol);
	}
	addmul_h2matrix(-1.0, work[(l - 1) + i * sons], false,
			inv->son[i + (l - 1) * sons],
			inv->son[(l - 1) + (l - 1) * sons], rwinv, cwinv, tm,
			tol);
      }
      for (i = l; i < sons; i++) {
	del_h2matrix(work[i + (l - 1) * sons]);
	del_h2matrix(work[(l - 1) + i * sons]);
      }
    }
    freemem(work);

    update_clusterbasis(inv->rb);
    update_clusterbasis(inv->cb);
    update_clusterbasis(h2->rb);
    update_clusterbasis(h2->cb);
  }
  else {
    assert(h2->f);
    copy_amatrix(false, h2->f, inv->f);
    qrinvert_amatrix(inv->f);
  }
}

/* ------------------------------------------------------------
 * Forward / backward substitution
 * ------------------------------------------------------------ */

/* ------------------------------------------------------------
 * lower-/uppersolve_h2matrix_avector
 * ------------------------------------------------------------ */

/* x = L^{-1}*x
 * L is lower triangular part of a */
void
lowersolve_h2matrix_1_avector(bool unit, pch2matrix a, pavector x)
{
  uint      sons = a->rsons;
  pccluster t = a->rb->t;

  pavector  x1, x2, x1t, x2t;
  avector   sub1, sub2, tmp1, tmp2;
  uint      roff1, roff2;
  uint      i, j;

  assert(sons == a->csons);
  assert(t == a->cb->t);
  assert(t->size == x->dim);

  if (a->son != 0) {
    roff1 = 0;
    for (i = 0; i < sons; i++) {
      /* compute the i-th row-block */
      x1 = init_sub_avector(&sub1, x, t->son[i]->size, roff1);
      lowersolve_h2matrix_1_avector(unit, a->son[i + i * sons], x1);

      /* update the lower blocks */
      roff2 = roff1;
      for (j = i + 1; j < sons; j++) {
	roff2 += t->son[j - 1]->size;
	x2 = init_sub_avector(&sub2, x, t->son[j]->size, roff2);

	x1t = init_avector(&tmp1, a->son[j + i * sons]->cb->ktree);
	x2t = init_avector(&tmp2, a->son[j + i * sons]->rb->ktree);
	clear_avector(x2t);
	forward_nopermutation_clusterbasis_avector(a->son[j + i * sons]->cb,
						   x1, x1t);
	fastaddeval_h2matrix_avector(-1.0, a->son[j + i * sons], x1t, x2t);
	backward_nopermutation_clusterbasis_avector(a->son[j + i * sons]->rb,
						    x2t, x2);
	uninit_avector(x1t);
	uninit_avector(x2t);

	uninit_avector(x2);
      }

      uninit_avector(x1);
      roff1 += t->son[i]->size;
    }
  }
  else {
    assert(a->f != 0);
    triangularsolve_amatrix_avector(true, unit, false, a->f, x);
  }
}

/* x = R^{-T}*x 
 * R is upper triangular part of a */
static void
lowersolve_h2matrix_2_avector(bool unit, pch2matrix a, pavector x)
{
  uint      sons = a->rsons;
  pccluster t = a->rb->t;

  pavector  x1, x2, x1t, x2t;
  avector   sub1, sub2, tmp1, tmp2;
  uint      coff1, coff2;
  uint      i, j;

  assert(t == a->cb->t);
  assert(t->size == x->dim);

  if (a->son != 0) {
    coff1 = 0;
    for (i = 0; i < sons; i++) {
      /* compute the i-th row-block */
      x1 = init_sub_avector(&sub1, x, t->son[i]->size, coff1);
      lowersolve_h2matrix_2_avector(unit, a->son[i + i * sons], x1);

      /* update the lower blocks */
      coff2 = coff1;
      for (j = i + 1; j < sons; j++) {
	coff2 += t->son[j - 1]->size;
	x2 = init_sub_avector(&sub2, x, t->son[j]->size, coff2);

	x1t = init_avector(&tmp1, a->son[i + j * sons]->rb->ktree);
	x2t = init_avector(&tmp2, a->son[i + j * sons]->cb->ktree);
	clear_avector(x2t);
	forward_nopermutation_clusterbasis_avector(a->son[i + j * sons]->rb,
						   x1, x1t);
	fastaddevaltrans_h2matrix_avector(-1.0, a->son[i + j * sons], x1t,
					  x2t);
	backward_nopermutation_clusterbasis_avector(a->son[i + j * sons]->cb,
						    x2t, x2);
	uninit_avector(x1t);
	uninit_avector(x2t);

	uninit_avector(x2);
      }

      uninit_avector(x1);
      coff1 += t->son[i]->size;
    }
  }
  else {
    assert(a->f != 0);
    triangularsolve_amatrix_avector(false, unit, true, a->f, x);
  }
}

/* x = T^{-1}*x 
 * T is lower triangular matrix
 * lower triangular or transposed upper triangular part of a */
void
lowersolve_h2matrix_avector(bool unit, bool atrans, pch2matrix a, pavector x)
{
  if (!atrans)
    lowersolve_h2matrix_1_avector(unit, a, x);
  else
    lowersolve_h2matrix_2_avector(unit, a, x);
}

/* x = R^{-1}*x 
 * R is upper triangular part of a */
static void
uppersolve_h2matrix_1_avector(bool unit, pch2matrix a, pavector x)
{
  uint      sons = a->rsons;
  pccluster t = a->rb->t;

  pavector  x1, x2, x1t, x2t;
  avector   sub1, sub2, tmp1, tmp2;
  uint      roff1, roff2;
  uint      i, j;

  assert(t == a->cb->t);
  assert(t->size == x->dim);

  if (a->son != 0) {
    roff1 = t->size;
    for (i = sons; i-- > 0;) {
      roff1 -= t->son[i]->size;
      /* compute the i-th row-block */
      x1 = init_sub_avector(&sub1, x, t->son[i]->size, roff1);
      uppersolve_h2matrix_1_avector(unit, a->son[i + i * sons], x1);

      /* update the upper blocks */
      roff2 = roff1;
      for (j = i; j-- > 0;) {
	roff2 -= t->son[j]->size;
	x2 = init_sub_avector(&sub2, x, t->son[j]->size, roff2);

	x1t = init_avector(&tmp1, a->son[j + i * sons]->cb->ktree);
	x2t = init_avector(&tmp2, a->son[j + i * sons]->rb->ktree);
	clear_avector(x2t);
	forward_nopermutation_clusterbasis_avector(a->son[j + i * sons]->cb,
						   x1, x1t);
	fastaddeval_h2matrix_avector(-1.0, a->son[j + i * sons], x1t, x2t);
	backward_nopermutation_clusterbasis_avector(a->son[j + i * sons]->rb,
						    x2t, x2);
	uninit_avector(x1t);
	uninit_avector(x2t);

	uninit_avector(x2);
      }

      uninit_avector(x1);
    }
  }
  else {
    assert(a->f != 0);
    triangularsolve_amatrix_avector(false, unit, false, a->f, x);
  }
}

/* x = L^{-T}*x 
 * L is lower triangular part of a */
static void
uppersolve_h2matrix_2_avector(bool unit, pch2matrix a, pavector x)
{
  uint      sons = a->rsons;
  pccluster t = a->rb->t;

  pavector  x1, x2, x1t, x2t;
  avector   sub1, sub2, tmp1, tmp2;
  uint      coff1, coff2;
  uint      i, j;

  assert(t == a->cb->t);
  assert(t->size == x->dim);

  if (a->son != 0) {
    coff1 = t->size;
    for (i = sons; i-- > 0;) {
      coff1 -= t->son[i]->size;
      /* compute the i-th row-block */
      x1 = init_sub_avector(&sub1, x, t->son[i]->size, coff1);
      uppersolve_h2matrix_2_avector(unit, a->son[i + i * sons], x1);

      /* update the upper blocks */
      coff2 = coff1;
      for (j = i; j-- > 0;) {
	coff2 -= t->son[j]->size;
	x2 = init_sub_avector(&sub2, x, t->son[j]->size, coff2);

	x1t = init_avector(&tmp1, a->son[i + j * sons]->rb->ktree);
	x2t = init_avector(&tmp2, a->son[i + j * sons]->cb->ktree);
	clear_avector(x2t);
	forward_nopermutation_clusterbasis_avector(a->son[i + j * sons]->rb,
						   x1, x1t);
	fastaddevaltrans_h2matrix_avector(-1.0, a->son[i + j * sons], x1t,
					  x2t);
	backward_nopermutation_clusterbasis_avector(a->son[i + j * sons]->cb,
						    x2t, x2);
	uninit_avector(x1t);
	uninit_avector(x2t);

	uninit_avector(x2);
      }

      uninit_avector(x1);
    }
  }
  else {
    assert(a->f != 0);
    triangularsolve_amatrix_avector(true, unit, true, a->f, x);
  }
}

/* x = T^{-1}*x 
 * T is upper triangular matrix
 * upper triangular or transposed lower triangular part of a */
void
uppersolve_h2matrix_avector(bool unit, bool atrans, pch2matrix a, pavector x)
{
  if (!atrans)
    uppersolve_h2matrix_1_avector(unit, a, x);
  else
    uppersolve_h2matrix_2_avector(unit, a, x);
}

/* ------------------------------------------------------------
 * lower-/uppersolve_h2matrix_amatrix
 * ------------------------------------------------------------ */

/* X = T^{-1}*X 
 * T is lower triangular part of A */
static void
lowersolve_h2matrix_1_amatrix(bool unit, pch2matrix a, bool xtrans,
			      pamatrix x)
{
  uint      sons = a->rsons;
  pccluster t = a->rb->t;

  pamatrix  x1, x2;
  amatrix   sub1, sub2;
  uint      roff1, roff2;
  uint      i, j;

  assert(sons == a->csons);
  assert(t == a->cb->t);
  assert((xtrans && (t->size == x->cols))
	 || (!xtrans && (t->size == x->rows)));

  if (a->son != 0) {
    roff1 = 0;
    for (i = 0; i < sons; i++) {
      /* compute the i-th row-block */
      if (xtrans == true)
	x1 = init_sub_amatrix(&sub1, x, x->rows, 0, t->son[i]->size, roff1);
      else
	x1 = init_sub_amatrix(&sub1, x, t->son[i]->size, roff1, x->cols, 0);

      lowersolve_h2matrix_1_amatrix(unit, a->son[i + i * sons], xtrans, x1);

      /* update the lower blocks */
      roff2 = roff1;
      for (j = i + 1; j < sons; j++) {
	roff2 += t->son[j - 1]->size;
	if (xtrans == true) {
	  x2 = init_sub_amatrix(&sub2, x, x->rows, 0, t->son[j]->size, roff2);
	  addmul_amatrix_h2matrix_amatrix(-1.0, false, x1, true,
					  a->son[j + i * sons], x2);
	}
	else {
	  x2 = init_sub_amatrix(&sub2, x, t->son[j]->size, roff2, x->cols, 0);
	  addmul_h2matrix_amatrix_amatrix(-1.0, false, a->son[j + i * sons],
					  false, x1, x2);
	}
	uninit_amatrix(x2);
      }

      uninit_amatrix(x1);
      roff1 += t->son[i]->size;
    }
  }
  else {
    assert(a->f != 0);
    triangularsolve_amatrix(true, unit, false, a->f, xtrans, x);
  }
}

/* X = T^{-T}*X
 * T is upper triangular part of A */
static void
lowersolve_h2matrix_2_amatrix(bool unit, pch2matrix a, bool xtrans,
			      pamatrix x)
{
  uint      sons = a->rsons;
  pccluster t = a->rb->t;

  pamatrix  x1, x2;
  amatrix   sub1, sub2;
  uint      coff1, coff2;
  uint      i, j;

  assert(sons == a->csons);
  assert(t == a->cb->t);
  assert((xtrans && (t->size == x->cols))
	 || (!xtrans && (t->size == x->rows)));

  if (a->son != 0) {
    coff1 = 0;
    for (i = 0; i < sons; i++) {
      /* compute the i-th row-block */
      if (!xtrans)
	x1 = init_sub_amatrix(&sub1, x, t->son[i]->size, coff1, x->cols, 0);
      else
	x1 = init_sub_amatrix(&sub1, x, x->rows, 0, t->son[i]->size, coff1);

      lowersolve_h2matrix_2_amatrix(unit, a->son[i + i * sons], xtrans, x1);

      /* update the lower block */
      coff2 = coff1;
      for (j = i + 1; j < sons; j++) {
	coff2 += t->son[j - 1]->size;
	if (!xtrans) {
	  x2 = init_sub_amatrix(&sub2, x, t->son[j]->size, coff2, x->cols, 0);
	  addmul_h2matrix_amatrix_amatrix(-1.0, true, a->son[i + j * sons],
					  false, x1, x2);
	}
	else {
	  x2 = init_sub_amatrix(&sub2, x, x->rows, 0, t->son[j]->size, coff2);
	  addmul_amatrix_h2matrix_amatrix(-1.0, false, x1, false,
					  a->son[i + j * sons], x2);
	}
	uninit_amatrix(x2);
      }
      uninit_amatrix(x1);
      coff1 += t->son[i]->size;
    }
  }
  else {
    assert(a->f != 0);
    triangularsolve_amatrix(false, unit, true, a->f, xtrans, x);
  }
}

void
lowersolve_h2matrix_amatrix(bool unit, bool atrans, pch2matrix A,
			    bool xtrans, pamatrix X)
{
  if (!atrans)
    lowersolve_h2matrix_1_amatrix(unit, A, xtrans, X);
  else
    lowersolve_h2matrix_2_amatrix(unit, A, xtrans, X);
}

/* X = T^{-1}*X 
 * T is upper triangular part of A */
static void
uppersolve_h2matrix_1_amatrix(bool unit, pch2matrix a, bool xtrans,
			      pamatrix x)
{
  uint      sons = a->rsons;
  pccluster t = a->rb->t;

  pamatrix  x1, x2;
  amatrix   sub1, sub2;
  uint      roff1, roff2;
  uint      i, j;

  assert(t == a->cb->t);
  if (xtrans == true) {
    assert(t->size == x->cols);
  }
  else {
    assert(xtrans == false);
    assert(t->size == x->rows);
  }

  if (a->son != 0) {
    roff1 = t->size;
    for (i = sons; i-- > 0;) {
      roff1 -= t->son[i]->size;
      /* compute the i-th row-block */
      if (xtrans == true) {
	x1 = init_sub_amatrix(&sub1, x, x->rows, 0, t->son[i]->size, roff1);
      }
      else {
	x1 = init_sub_amatrix(&sub1, x, t->son[i]->size, roff1, x->cols, 0);
      }
      uppersolve_h2matrix_1_amatrix(unit, a->son[i + i * sons], xtrans, x1);

      /* update the lower blocks */
      roff2 = roff1;
      for (j = i; j-- > 0;) {
	roff2 -= t->son[j]->size;
	if (xtrans == true) {
	  x2 = init_sub_amatrix(&sub2, x, x->rows, 0, t->son[j]->size, roff2);
	  addmul_amatrix_h2matrix_amatrix(-1.0, false, x1, true,
					  a->son[j + i * sons], x2);
	}
	else {
	  x2 = init_sub_amatrix(&sub2, x, t->son[j]->size, roff2, x->cols, 0);
	  addmul_h2matrix_amatrix_amatrix(-1.0, false, a->son[j + i * sons],
					  false, x1, x2);
	}
	uninit_amatrix(x2);
      }

      uninit_amatrix(x1);
    }
  }
  else {
    assert(a->f != 0);
    triangularsolve_amatrix(false, unit, false, a->f, xtrans, x);
  }
}

/* X = T^{-T}*X
 * T is upper triangular part of A */
static void
uppersolve_h2matrix_2_amatrix(bool unit, pch2matrix a, bool xtrans,
			      pamatrix x)
{
  uint      sons = a->rsons;
  pccluster t = a->rb->t;

  pamatrix  x1, x2;
  amatrix   sub1, sub2;
  uint      coff1, coff2;
  uint      i, j;

  assert(t == a->cb->t);
  if (xtrans == true) {
    assert(t->size == x->rows);
  }
  else {
    assert(xtrans == false);
    assert(t->size == x->cols);
  }

  if (a->son != 0) {
    coff1 = t->size;
    for (i = sons; i-- > 0;) {
      coff1 -= t->son[i]->size;
      /* compute the i-th row-block */
      if (xtrans == true) {
	x1 = init_sub_amatrix(&sub1, x, t->son[i]->size, coff1, x->cols, 0);
      }
      else {
	x1 = init_sub_amatrix(&sub1, x, x->rows, 0, t->son[i]->size, coff1);
      }
      uppersolve_h2matrix_2_amatrix(unit, a->son[i + i * sons], xtrans, x1);

      /* update the lower block */
      coff2 = coff1;
      for (j = i; j-- > 0;) {
	coff2 -= t->son[j]->size;
	if (xtrans == true) {
	  x2 = init_sub_amatrix(&sub2, x, t->son[j]->size, coff2, x->cols, 0);
	  addmul_h2matrix_amatrix_amatrix(-1.0, true, a->son[i + j * sons],
					  false, x1, x2);
	}
	else {
	  x2 = init_sub_amatrix(&sub2, x, x->rows, 0, t->son[j]->size, coff2);
	  addmul_amatrix_h2matrix_amatrix(-1.0, false, x1, false,
					  a->son[i + j * sons], x2);
	}
	uninit_amatrix(x2);
      }
      uninit_amatrix(x1);
    }
  }
  else {
    assert(a->f != 0);
    triangularsolve_amatrix(true, unit, true, a->f, xtrans, x);
  }
}

void
uppersolve_h2matrix_amatrix(bool unit, bool atrans, pch2matrix A,
			    bool xtrans, pamatrix X)
{
  if (!atrans)
    uppersolve_h2matrix_1_amatrix(unit, A, xtrans, X);
  else
    uppersolve_h2matrix_2_amatrix(unit, A, xtrans, X);
}

/* ------------------------------------------------------------
 * lower-/uppersolve_amatrix_h2matrix
 * ------------------------------------------------------------ */

/* R = T^{-1}*X 
 * T is lower triangular part of L */
static void
lowersolve_amatrix_1_h2matrix(bool unit, pcamatrix L, pch2matrix X,
			      ph2matrix R, pclusteroperator rwfup,
			      pclusteroperator cwfup, ptruncmode tm, real tol)
{
  uint      rsons = X->rsons;
  uint      csons = X->csons;

  prkmatrix r;
  pclusteroperator rwup, cwup;
  uint      j;

  assert(L->rows == L->cols);
  assert(L->cols == X->rb->t->size);

  if (X->son != 0) {
    assert(rsons == 1);
    assert(X->rb->t->son == 0);

    rwup = identify_son_clusterweight_clusteroperator(rwfup, X->rb->t);
    cwup = identify_son_clusterweight_clusteroperator(cwfup, X->cb->t);

    for (j = 0; j < csons; j++)
      lowersolve_amatrix_1_h2matrix(unit, L, X->son[j], R->son[j], rwup, cwup,
				    tm, tol);

    update_clusterbasis(R->rb);
    update_clusterbasis(R->cb);
  }
  else if (X->u != 0) {
    r = convert_uniform_rkmatrix(false, X->u);
    triangularsolve_amatrix(true, unit, false, L, false, &r->A);
    rkupdate_h2matrix(r, R, rwfup, cwfup, tm, tol);
    del_rkmatrix(r);
  }
  else if (X->f != 0) {
    copy_amatrix(false, X->f, R->f);
    triangularsolve_amatrix(true, unit, false, L, false, R->f);
  }
}

/* R^T = T^{-1}*X^T 
 * T is lower triangular part of L */
static void
lowersolve_amatrix_2_h2matrix(bool unit, pcamatrix L, pch2matrix X,
			      ph2matrix R, pclusteroperator rwfup,
			      pclusteroperator cwfup, ptruncmode tm, real tol)
{
  uint      rsons = X->csons;
  uint      csons = X->rsons;

  prkmatrix r;
  pclusteroperator rwup, cwup;
  uint      j;

  assert(L->rows == L->cols);
  assert(L->cols == X->cb->t->size);

  if (X->son != 0) {
    assert(rsons == 1);
    assert(X->cb->t->son == 0);

    rwup = identify_son_clusterweight_clusteroperator(rwfup, X->rb->t);
    cwup = identify_son_clusterweight_clusteroperator(cwfup, X->cb->t);

    for (j = 0; j < csons; j++)
      lowersolve_amatrix_2_h2matrix(unit, L, X->son[j], R->son[j], rwup, cwup,
				    tm, tol);

    update_clusterbasis(R->rb);
    update_clusterbasis(R->cb);
  }
  else if (X->u != 0) {
    r = convert_uniform_rkmatrix(false, X->u);
    clear_amatrix(&X->u->S);
    triangularsolve_amatrix(true, unit, false, L, false, &r->B);
    rkupdate_h2matrix(r, R, rwfup, cwfup, tm, tol);
    del_rkmatrix(r);
  }
  else if (X->f != 0) {
    copy_amatrix(false, X->f, R->f);
    triangularsolve_amatrix(true, unit, false, L, true, R->f);
  }
}

/* L = T^{-T}*X 
 * T is upper triangular part of R */
static void
lowersolve_amatrix_3_h2matrix(bool unit, pcamatrix L, pch2matrix X,
			      ph2matrix R, pclusteroperator rwfup,
			      pclusteroperator cwfup, ptruncmode tm, real tol)
{
  (void) unit;
  (void) L;
  (void) X;
  (void) R;
  (void) rwfup;
  (void) cwfup;
  (void) tm;
  (void) tol;

  printf("The function lowersolve_amatrix_h2matrix_3 has not yet "
	 "been implemented\n");
  abort();
}

/* L^T = T^{-T}*X^T 
 * T is upper triangular part of R */
static void
lowersolve_amatrix_4_h2matrix(bool unit, pcamatrix R, pch2matrix X,
			      ph2matrix L, pclusteroperator rwflow,
			      pclusteroperator cwflow, ptruncmode tm,
			      real tol)
{
  uint      rsons = X->rsons;
  uint      csons = X->csons;

  prkmatrix r;
  pclusteroperator rwlow, cwlow;
  uint      j;

  assert(R->rows == R->cols);
  assert(R->rows == X->cb->t->size);

  if (X->son != 0) {
    assert(csons == 1);
    assert(X->cb->t->son == 0);

    cwlow = cwflow;
    rwlow = identify_son_clusterweight_clusteroperator(rwflow, X->rb->t);

    for (j = 0; j < rsons; j++)
      lowersolve_amatrix_h2matrix(unit, true, R, true, X->son[j], L->son[j],
				  rwlow, cwlow, tm, tol);

    update_clusterbasis(L->rb);
    update_clusterbasis(L->cb);
  }
  else if (X->u != 0) {
    r = convert_uniform_rkmatrix(false, X->u);
    triangularsolve_amatrix(false, unit, true, R, false, &r->B);
    rkupdate_h2matrix(r, L, rwflow, cwflow, tm, tol);
    del_rkmatrix(r);
  }
  else if (X->f != 0) {
    copy_amatrix(false, X->f, L->f);
    triangularsolve_amatrix(false, unit, true, R, true, L->f);
  }
}

void
lowersolve_amatrix_h2matrix(bool unit, bool atrans, pcamatrix A,
			    bool xytrans, pch2matrix X, ph2matrix Y,
			    pclusteroperator rwf, pclusteroperator cwf,
			    ptruncmode tm, real tol)
{
  if (!atrans) {
    if (!xytrans)
      lowersolve_amatrix_1_h2matrix(unit, A, X, Y, rwf, cwf, tm, tol);
    else
      lowersolve_amatrix_2_h2matrix(unit, A, X, Y, rwf, cwf, tm, tol);
  }
  else {
    if (!xytrans)
      lowersolve_amatrix_3_h2matrix(unit, A, X, Y, rwf, cwf, tm, tol);
    else
      lowersolve_amatrix_4_h2matrix(unit, A, X, Y, rwf, cwf, tm, tol);
  }
}

/* R = T^{-1}*X 
 * T is lower triangular part of L */
static void
uppersolve_amatrix_1_h2matrix(bool unit, pcamatrix a, ph2matrix x,
			      pclusteroperator rwf, pclusteroperator cwf,
			      ptruncmode tm, real tol)
{
  uint      rsons = x->rsons;
  uint      csons = x->csons;

  prkmatrix r, p;
  pclusteroperator rw, cw;
  uint      i, j;

  assert(a->rows == a->cols);
  assert(a->cols == x->rb->t->size);

  if (x->son != 0) {
    assert(rsons == 1);
    assert(x->rb->t->son == 0);
    rw = rwf;

    cw = cwf;
    i = 0;
    while (cw->t != x->cb->t && i < cwf->sons) {
      cw = cwf->son[i];
      i++;
    }
    assert(cw->t == x->cb->t);
    if (cw->sons == 0)
      cw = cwf;

    for (j = 0; j < csons; j++)
      uppersolve_amatrix_1_h2matrix(unit, a, x->son[j], rw, cw, tm, tol);

    update_clusterbasis(x->rb);
    update_clusterbasis(x->cb);
  }
  else if (x->u != 0) {
    r = convert_uniform_rkmatrix(false, x->u);
    p = new_rkmatrix(r->A.rows, r->B.rows, r->k);
    copy_amatrix(false, &r->A, &p->A);
    copy_amatrix(false, &r->B, &p->B);
    triangularsolve_amatrix(false, unit, false, a, false, &r->A);
    add_rkmatrix(-1.0, p, 0, tol, r);
    rkupdate_h2matrix(r, x, rwf, cwf, tm, tol);
    del_rkmatrix(r);
  }
  else if (x->f != 0) {
    triangularsolve_amatrix(false, unit, false, a, false, x->f);
  }
}

/* R^T = T^{-1}*X^T 
 * T is upper triangular part of L */
static void
uppersolve_amatrix_2_h2matrix(bool unit, pcamatrix L, ph2matrix R,
			      pclusteroperator rwfup, pclusteroperator cwfup,
			      ptruncmode tm, real tol)
{
  (void) unit;
  (void) L;
  (void) R;
  (void) rwfup;
  (void) cwfup;
  (void) tm;
  (void) tol;

  printf("The function uppersolve_amatrix_h2matrix_2 has not yet "
	 "been implemented\n");
  abort();
}

/* L = T^{-T}*X 
 * T is lower triangular part of R */
static void
uppersolve_amatrix_3_h2matrix(bool unit, pcamatrix L, ph2matrix R,
			      pclusteroperator rwfup, pclusteroperator cwfup,
			      ptruncmode tm, real tol)
{
  (void) unit;
  (void) L;
  (void) R;
  (void) rwfup;
  (void) cwfup;
  (void) tm;
  (void) tol;

  printf("The function uppersolve_amatrix_h2matrix_3 has not yet "
	 "been implemented\n");
  abort();
}

/* L^T = T^{-T}*X^T 
 * T is lower triangular part of R */
static void
uppersolve_amatrix_4_h2matrix(bool unit, pcamatrix a, ph2matrix x,
			      pclusteroperator rwf, pclusteroperator cwf,
			      ptruncmode tm, real tol)
{
  uint      rsons = x->rsons;
  uint      csons = x->csons;

  prkmatrix r, p;
  pclusteroperator rw, cw;
  uint      i, j;

  assert(a->rows == a->cols);
  assert(a->rows == x->cb->t->size);

  if (x->son != 0) {
    assert(csons == 1);
    assert(x->cb->t->son == 0);
    cw = cwf;

    rw = rwf;
    i = 0;
    while (rw->t != x->rb->t && i < rwf->sons) {
      rw = rwf->son[i];
      i++;
    }
    assert(rw->t == x->rb->t);
    if (rw->sons == 0)
      rw = rwf;

    for (j = 0; j < rsons; j++)
      uppersolve_amatrix_4_h2matrix(unit, a, x->son[j], rw, cw, tm, tol);

    update_clusterbasis(x->rb);
    update_clusterbasis(x->cb);
  }
  else if (x->u != 0) {
    r = convert_uniform_rkmatrix(false, x->u);
    p = new_rkmatrix(r->A.rows, r->B.rows, r->k);
    copy_amatrix(false, &r->A, &p->A);
    copy_amatrix(false, &r->B, &p->B);
    triangularsolve_amatrix(true, unit, true, a, false, &r->B);
    add_rkmatrix(-1.0, p, 0, tol, r);
    rkupdate_h2matrix(r, x, rwf, cwf, tm, tol);
    del_rkmatrix(r);
  }
  else if (x->f != 0) {
    triangularsolve_amatrix(true, unit, true, a, true, x->f);
  }
}

void
uppersolve_amatrix_h2matrix(bool unit, bool atrans, pcamatrix A,
			    bool ytrans, ph2matrix Y, pclusteroperator rwf,
			    pclusteroperator cwf, ptruncmode tm, real tol)
{
  if (!atrans) {
    if (!ytrans)
      uppersolve_amatrix_1_h2matrix(unit, A, Y, rwf, cwf, tm, tol);
    else
      uppersolve_amatrix_2_h2matrix(unit, A, Y, rwf, cwf, tm, tol);
  }
  else {
    if (!ytrans)
      uppersolve_amatrix_3_h2matrix(unit, A, Y, rwf, cwf, tm, tol);
    else
      uppersolve_amatrix_4_h2matrix(unit, A, Y, rwf, cwf, tm, tol);
  }
}

/* ------------------------------------------------------------
 lower-/uppersolve_h2matrix
 ------------------------------------------------------------ */

/* R = T^{-1}*X 
 * T is lower triangular part of L */
static void
lowersolve_h2matrix_1_h2matrix(bool unit, pch2matrix L, ph2matrix X,
			       pclusteroperator rwf, pclusteroperator cwf,
			       ph2matrix R, pclusteroperator rwfup,
			       pclusteroperator cwfup, ptruncmode tm,
			       real tol)
{
  uint      rsons = L->rsons;
  uint      csons = X->csons;
  pccluster t = L->rb->t;

  prkmatrix r;
  pclusteroperator rw, cw, rwup, cwup;
  uint      i, j, l;

  assert(X->rb->t == R->rb->t);
  assert(X->cb->t == R->cb->t);
  assert(X->rsons == R->rsons);
  assert(X->csons == R->csons);
  assert(rsons == L->csons);
  assert(t == L->cb->t);
  assert(t == X->rb->t);

  if (L->son != 0) {
    if (X->son != 0) {
      assert(rsons == X->rsons);
      rw = identify_son_clusterweight_clusteroperator(rwf, X->rb->t);
      cw = identify_son_clusterweight_clusteroperator(cwf, X->cb->t);
      rwup = identify_son_clusterweight_clusteroperator(rwfup, R->rb->t);
      cwup = identify_son_clusterweight_clusteroperator(cwfup, R->cb->t);

      for (i = 0; i < rsons; i++) {
	for (j = 0; j < csons; j++) {
	  /* compute the i-th row */
	  lowersolve_h2matrix_1_h2matrix(unit, L->son[i + i * rsons],
					 X->son[i + j * rsons], rw, cw,
					 R->son[i + j * rsons], rwup, cwup,
					 tm, tol);

	  /* update the lower rows */
	  for (l = i + 1; l < rsons; l++) {
	    addmul_h2matrix(-1.0, L->son[l + i * rsons], false,
			    R->son[i + j * rsons], X->son[l + j * rsons], rw,
			    cw, tm, tol);
	  }
	}
      }
      update_clusterbasis(X->rb);
      update_clusterbasis(X->cb);
      update_clusterbasis(R->rb);
      update_clusterbasis(R->cb);
    }
    else if (X->u != 0) {
      r = convert_uniform_rkmatrix(false, X->u);
      lowersolve_h2matrix_amatrix(unit, false, L, false, &r->A);
      rkupdate_h2matrix(r, R, rwfup, cwfup, tm, tol);
      del_rkmatrix(r);
    }
    else if (X->f != 0) {
      copy_amatrix(false, X->f, R->f);
      lowersolve_h2matrix_amatrix(unit, false, L, false, R->f);
    }
  }
  else {
    assert(L->f != 0);
    lowersolve_amatrix_h2matrix(unit, false, L->f, false, X, R, rwfup, cwfup,
				tm, tol);
  }
}

/* R^* = T^{-1}*X^* 
 T is lower triangular part of L */
static void
lowersolve_h2matrix_2_h2matrix(bool unit, pch2matrix L, ph2matrix X,
			       pclusteroperator rwf, pclusteroperator cwf,
			       ph2matrix R, pclusteroperator rwfup,
			       pclusteroperator cwfup, ptruncmode tm,
			       real tol)
{
  uint      rsons = L->rsons;
  uint      csons = X->rsons;
  pccluster t = L->rb->t;

  prkmatrix r;
  pclusteroperator rw, cw, rwup, cwup;
  uint      i, j, l;

  assert(X->rb->t == R->rb->t);
  assert(X->cb->t == R->cb->t);
  assert(X->rsons == R->rsons);
  assert(X->csons == R->csons);
  assert(rsons == L->csons);
  assert(t == L->cb->t);
  assert(t == X->cb->t);

  if (L->son != 0) {
    if (X->son != 0) {
      assert(rsons == X->csons);
      rw = identify_son_clusterweight_clusteroperator(rwf, X->rb->t);
      cw = identify_son_clusterweight_clusteroperator(cwf, X->cb->t);
      rwup = identify_son_clusterweight_clusteroperator(rwfup, R->rb->t);
      cwup = identify_son_clusterweight_clusteroperator(cwfup, R->cb->t);

      for (i = 0; i < rsons; i++) {
	for (j = 0; j < csons; j++) {
	  /* compute the i-th row */
	  lowersolve_h2matrix_2_h2matrix(unit, L->son[i + i * rsons],
					 X->son[j + i * csons], rw, cw,
					 R->son[j + i * csons], rwup, cwup,
					 tm, tol);

	  /* update the lower rows */
	  for (l = i + 1; l < rsons; l++) {
	    addmul_h2matrix(-1.0, R->son[j + i * csons], true,
			    L->son[l + i * rsons], X->son[j + l * csons], rw,
			    cw, tm, tol);
	  }
	}
      }
      update_clusterbasis(X->rb);
      update_clusterbasis(X->cb);
      update_clusterbasis(R->rb);
      update_clusterbasis(R->cb);
    }
    else if (X->u != 0) {
      r = convert_uniform_rkmatrix(false, X->u);
      lowersolve_h2matrix_amatrix(unit, false, L, false, &r->B);
      rkupdate_h2matrix(r, R, rwfup, cwfup, tm, tol);
      del_rkmatrix(r);
    }
    else if (X->f != 0) {
      copy_amatrix(false, X->f, R->f);
      lowersolve_h2matrix_amatrix(unit, false, L, true, R->f);
    }
  }
  else {
    assert(L->f != 0);
    lowersolve_amatrix_h2matrix(unit, false, L->f, true, X, R, rwfup, cwfup,
				tm, tol);
  }
}

/* L = T^{-T}X
 * T is upper triangular part of R */
static void
lowersolve_h2matrix_3_h2matrix(bool unit, pch2matrix R, ph2matrix X,
			       pclusteroperator rwf, pclusteroperator cwf,
			       ph2matrix L, pclusteroperator rwflow,
			       pclusteroperator cwflow, ptruncmode tm,
			       real tol)
{
  (void) unit;
  (void) L;
  (void) X;
  (void) R;
  (void) rwf;
  (void) cwf;
  (void) rwflow;
  (void) cwflow;
  (void) tm;
  (void) tol;

  printf("The function lowersolve_h2matrix_3_h2matrix has not yet "
	 "been implemented\n");
  abort();
}

/* L^* = T^{-T}X^*
 * T is upper triangular part of R */
static void
lowersolve_h2matrix_4_h2matrix(bool unit, pch2matrix R, ph2matrix X,
			       pclusteroperator rwf, pclusteroperator cwf,
			       ph2matrix L, pclusteroperator rwflow,
			       pclusteroperator cwflow, ptruncmode tm,
			       real tol)
{
  uint      rsons = X->rsons;
  uint      csons = R->csons;
  pccluster t = R->rb->t;

  prkmatrix r;
  pclusteroperator rw, cw, rwlow, cwlow;
  uint      i, j, l;

  assert(X->rb->t == L->rb->t);
  assert(X->cb->t == L->cb->t);
  assert(X->rsons == L->rsons);
  assert(X->csons == L->csons);
  assert(csons == R->rsons);
  assert(t == R->cb->t);
  assert(t == X->cb->t);

  if (R->son != 0) {
    if (X->son != 0) {
      assert(csons == X->csons);
      rw = identify_son_clusterweight_clusteroperator(rwf, X->rb->t);
      cw = identify_son_clusterweight_clusteroperator(cwf, X->cb->t);
      rwlow = identify_son_clusterweight_clusteroperator(rwflow, L->rb->t);
      cwlow = identify_son_clusterweight_clusteroperator(cwflow, L->cb->t);

      for (i = 0; i < csons; i++) {
	for (j = 0; j < rsons; j++) {
	  /* computes the i-th row */
	  lowersolve_h2matrix_4_h2matrix(unit, R->son[i + i * csons],
					 X->son[j + i * rsons], rw, cw,
					 L->son[j + i * rsons], rwlow, cwlow,
					 tm, tol);

	  /* update the lower block */
	  for (l = i + 1; l < csons; l++) {
	    addmul_h2matrix(-1.0, L->son[j + i * rsons], false,
			    R->son[i + l * csons], X->son[j + l * rsons], rw,
			    cw, tm, tol);
	  }
	}
      }
      update_clusterbasis(X->rb);
      update_clusterbasis(X->cb);
      update_clusterbasis(L->rb);
      update_clusterbasis(L->cb);
    }
    else if (X->u != 0) {
      r = convert_uniform_rkmatrix(false, X->u);
      lowersolve_h2matrix_amatrix(unit, true, R, false, &r->B);
      rkupdate_h2matrix(r, L, rwflow, cwflow, tm, tol);
      del_rkmatrix(r);
    }
    else if (X->f != 0) {
      copy_amatrix(false, X->f, L->f);
      lowersolve_h2matrix_amatrix(unit, true, R, true, L->f);
    }
  }
  else {
    assert(R->f != 0);
    lowersolve_amatrix_h2matrix(unit, true, R->f, true, X, L, rwflow, cwflow,
				tm, tol);
  }
}

void
lowersolve_h2matrix(bool aunit, bool atrans, pch2matrix A,
		    bool xytrans, ph2matrix X, pclusteroperator xrwf,
		    pclusteroperator xcwf, ph2matrix Y, pclusteroperator yrwf,
		    pclusteroperator ycwf, ptruncmode tm, real tol)
{
  if (!atrans) {
    if (!xytrans)
      lowersolve_h2matrix_1_h2matrix(aunit, A, X, xrwf, xcwf, Y, yrwf, ycwf,
				     tm, tol);
    else
      lowersolve_h2matrix_2_h2matrix(aunit, A, X, xrwf, xcwf, Y, yrwf, ycwf,
				     tm, tol);
  }
  else {
    if (!xytrans)
      lowersolve_h2matrix_3_h2matrix(aunit, A, X, xrwf, xcwf, Y, yrwf, ycwf,
				     tm, tol);
    else
      lowersolve_h2matrix_4_h2matrix(aunit, A, X, xrwf, xcwf, Y, yrwf, ycwf,
				     tm, tol);
  }
}

/* x = T^{-1}*x 
 * T is upper triangular part of a */
static void
uppersolve_h2matrix_1_h2matrix(bool unit, pch2matrix a, ph2matrix x,
			       pclusteroperator rwf, pclusteroperator cwf,
			       ptruncmode tm, real tol)
{
  uint      rsons = a->rsons;
  uint      csons = x->csons;
  pccluster t = a->rb->t;

  prkmatrix r, p;
  pclusteroperator rw, cw;
  uint      i, j, l;

  assert(t == a->cb->t);
  assert(t == x->rb->t);

  if (a->son != 0) {
    if (x->son != 0) {
      rw = rwf;
      i = 0;
      while (rw->t != x->rb->t && i < rwf->sons) {
	rw = rwf->son[i];
	i++;
      }
      assert(rw->t == x->rb->t);
      cw = cwf;
      i = 0;
      while (cw->t != x->cb->t && i < cwf->sons) {
	cw = cwf->son[i];
	i++;
      }
      assert(cw->t == x->cb->t);

      if (rw->sons == 0)
	rw = rwf;
      if (cw->sons == 0)
	cw = cwf;

      for (i = rsons; i-- > 0;) {
	for (j = csons; j-- > 0;) {
	  /* compute the i-th row */
	  uppersolve_h2matrix_1_h2matrix(unit, a->son[i + i * rsons],
					 x->son[i + j * rsons], rw, cw, tm,
					 tol);

	  /* update the upper rows */
	  for (l = i; l-- > 0;) {
	    addmul_h2matrix(-1.0, a->son[l + i * rsons], false,
			    x->son[i + j * rsons], x->son[l + j * rsons], rw,
			    cw, tm, tol);
	  }
	}
      }
      update_clusterbasis(x->rb);
      update_clusterbasis(x->cb);
    }
    else if (x->u != 0) {
      r = convert_uniform_rkmatrix(false, x->u);
      p = new_rkmatrix(r->A.rows, r->B.rows, r->k);
      copy_amatrix(false, &r->A, &p->A);
      copy_amatrix(false, &r->B, &p->B);
      uppersolve_h2matrix_amatrix(unit, false, a, false, &r->A);
      add_rkmatrix(-1.0, p, 0, tol, r);
      rkupdate_h2matrix(r, x, rwf, cwf, tm, tol);
      del_rkmatrix(r);
    }
    else if (x->f != 0) {
      uppersolve_h2matrix_amatrix(unit, false, a, false, x->f);
    }
  }
  else {
    assert(a->f != 0);
    uppersolve_amatrix_h2matrix(unit, false, a->f, false, x, rwf, cwf, tm,
				tol);
  }
}

/* x^* = T^{-1}*x^*
 * T is upper triangular part of a */
static void
uppersolve_h2matrix_2_h2matrix(bool unit, pch2matrix a, ph2matrix x,
			       pclusteroperator rwf, pclusteroperator cwf,
			       ptruncmode tm, real tol)
{
  (void) unit;
  (void) a;
  (void) x;
  (void) rwf;
  (void) cwf;
  (void) tm;
  (void) tol;

  printf("The function uppersolve_h2matrix_2_h2matrix has not yet "
	 "been implemented\n");
  abort();
}

/* x = T^{-T}*x
 * T is lower triangular part of a */
static void
uppersolve_h2matrix_3_h2matrix(bool unit, pch2matrix a, ph2matrix x,
			       pclusteroperator rwf, pclusteroperator cwf,
			       ptruncmode tm, real tol)
{
  (void) unit;
  (void) a;
  (void) x;
  (void) rwf;
  (void) cwf;
  (void) tm;
  (void) tol;

  printf("The function uppersolve_h2matrix_3_h2matrix has not yet "
	 "been implemented\n");
  abort();
}

/* x^* = T^{-T}*x^*
 * T is lower triangular part of a */
static void
uppersolve_h2matrix_4_h2matrix(bool unit, pch2matrix a, ph2matrix x,
			       pclusteroperator rwf, pclusteroperator cwf,
			       ptruncmode tm, real tol)
{
  uint      rsons = x->rsons;
  uint      csons = a->csons;
  pccluster t = a->rb->t;

  prkmatrix r, p;
  pclusteroperator rw, cw;
  uint      i, j, l;

  assert(t == a->cb->t);
  assert(t == x->cb->t);

  if (a->son != 0) {
    if (x->son != 0) {
      rw = rwf;
      i = 0;
      while (rw->t != x->rb->t && i < rwf->sons) {
	rw = rwf->son[i];
	i++;
      }
      assert(rw->t == x->rb->t);
      cw = cwf;
      i = 0;
      while (cw->t != x->cb->t && i < cwf->sons) {
	cw = cwf->son[i];
	i++;
      }
      assert(cw->t == x->cb->t);

      if (rw->sons == 0)
	rw = rwf;
      if (cw->sons == 0)
	cw = cwf;

      for (i = csons; i-- > 0;) {
	for (j = rsons; j-- > 0;) {
	  /* computes the i-th row */
	  uppersolve_h2matrix_4_h2matrix(unit, a->son[i + i * csons],
					 x->son[j + i * rsons], rw, cw, tm,
					 tol);

	  /* update the upper block */
	  for (l = i + 1; l < csons; l++) {
	    addmul_h2matrix(-1.0, x->son[j + i * rsons], false,
			    a->son[i + l * csons], x->son[j + l * rsons], rw,
			    cw, tm, tol);
	  }
	}
      }
      update_clusterbasis(x->rb);
      update_clusterbasis(x->cb);
    }
    else if (x->u != 0) {
      r = convert_uniform_rkmatrix(false, x->u);
      p = new_rkmatrix(r->A.rows, r->B.rows, r->k);
      copy_amatrix(false, &r->A, &p->A);
      copy_amatrix(false, &r->B, &p->B);
      uppersolve_h2matrix_amatrix(unit, true, a, true, &r->B);
      add_rkmatrix(-1.0, p, 0, tol, r);
      rkupdate_h2matrix(r, x, rwf, cwf, tm, tol);
      del_rkmatrix(r);
    }
    else if (x->f != 0) {
      uppersolve_h2matrix_amatrix(unit, true, a, false, x->f);
    }
  }
  else {
    assert(a->f != 0);
    uppersolve_amatrix_h2matrix(unit, true, a->f, true, x, rwf, cwf, tm, tol);
  }
}

void
uppersolve_h2matrix(bool aunit, bool atrans, pch2matrix A,
		    bool ytrans, ph2matrix Y, pclusteroperator yrwf,
		    pclusteroperator ycwf, ptruncmode tm, real tol)
{
  if (!atrans) {
    if (!ytrans)
      uppersolve_h2matrix_1_h2matrix(aunit, A, Y, yrwf, ycwf, tm, tol);
    else
      uppersolve_h2matrix_2_h2matrix(aunit, A, Y, yrwf, ycwf, tm, tol);
  }
  else {
    if (!ytrans)
      uppersolve_h2matrix_3_h2matrix(aunit, A, Y, yrwf, ycwf, tm, tol);
    else
      uppersolve_h2matrix_4_h2matrix(aunit, A, Y, yrwf, ycwf, tm, tol);
  }
}

/* ------------------------------------------------------------
 LR decomposition
 ------------------------------------------------------------ */

void
lrdecomp_h2matrix(ph2matrix X, pclusteroperator rwf, pclusteroperator cwf,
		  ph2matrix L, pclusteroperator rwflow,
		  pclusteroperator cwflow, ph2matrix R,
		  pclusteroperator rwfup, pclusteroperator cwfup,
		  ptruncmode tm, real tol)
{
  uint      sons = X->rsons;
  pccluster t = X->rb->t;

  pclusteroperator rw, cw, rwlow, cwlow, rwup, cwup;
  uint      i, j, k /*, size */ ;

  assert(t == X->cb->t);
  assert(t == L->rb->t);
  assert(t == L->cb->t);
  assert(t == R->rb->t);
  assert(t == R->cb->t);

  if (X->son != 0) {
    rw = identify_son_clusterweight_clusteroperator(rwf, t);
    cw = identify_son_clusterweight_clusteroperator(cwf, t);
    rwlow = identify_son_clusterweight_clusteroperator(rwflow, t);
    cwlow = identify_son_clusterweight_clusteroperator(cwflow, t);
    rwup = identify_son_clusterweight_clusteroperator(rwfup, t);
    cwup = identify_son_clusterweight_clusteroperator(cwfup, t);

    for (i = 0; i < sons; i++) {

      /* compute the i-th diagonal block */
      lrdecomp_h2matrix(X->son[i + i * sons], rw, cw, L->son[i + i * sons],
			rwlow, cwlow, R->son[i + i * sons], rwup, cwup, tm,
			tol);

      for (j = i + 1; j < sons; j++) {

	/* compute the i-th row */
	lowersolve_h2matrix(true, false, L->son[i + i * sons], false,
			    X->son[i + j * sons], rw, cw,
			    R->son[i + j * sons], rwup, cwup, tm, tol);

	/* compute the i-th column */
	lowersolve_h2matrix(false, true, R->son[i + i * sons], true,
			    X->son[j + i * sons], rw, cw,
			    L->son[j + i * sons], rwlow, cwlow, tm, tol);
      }

      /* update the lower right block */
      for (j = i + 1; j < sons; j++) {
	for (k = i + 1; k < sons; k++) {
	  addmul_h2matrix(-1.0, L->son[j + i * sons], false,
			  R->son[i + k * sons], X->son[j + k * sons], rw, cw,
			  tm, tol);
	}
      }
    }
    update_clusterbasis(X->rb);
    update_clusterbasis(X->cb);
    update_clusterbasis(L->rb);
    update_clusterbasis(L->cb);
    update_clusterbasis(R->rb);
    update_clusterbasis(R->cb);
  }
  else {
    assert(X->f != 0);
    lrdecomp_amatrix(X->f);
    copy_lower_amatrix(X->f, true, L->f);
    copy_upper_amatrix(X->f, false, R->f);
  }
}

void
lrsolve_h2matrix_avector(pch2matrix L, pch2matrix R, pavector x)
{
  pavector  xp;
  avector   xtmp;
  uint      i, ip;

  /*assert(x->dim == h2->cb->t->size); */

  /* Permutation of x */
  xp = init_avector(&xtmp, x->dim);
  for (i = 0; i < xp->dim; i++) {
    ip = L->cb->t->idx[i];
    assert(ip < x->dim);
    xp->v[i] = x->v[ip];
  }

  lowersolve_h2matrix_avector(true, false, L, xp);
  uppersolve_h2matrix_avector(false, false, R, xp);

  /* Reverse permutation of x */
  for (i = 0; i < xp->dim; i++) {
    ip = R->rb->t->idx[i];
    assert(ip < x->dim);
    x->v[ip] = xp->v[i];
  }
  uninit_avector(xp);
}

/* ------------------------------------------------------------
 Cholesky decomposition
 ------------------------------------------------------------ */

void
init_cholesky_h2matrix(ph2matrix A, pclusteroperator * prwf,
		       pclusteroperator * pcwf, ptruncmode tm)
{
  tolower_h2matrix(A);
  if (*prwf)
    del_clusteroperator(*prwf);
  if (*pcwf)
    del_clusteroperator(*pcwf);
  *prwf = prepare_row_clusteroperator(A->rb, A->cb, tm);
  *pcwf = prepare_col_clusteroperator(A->rb, A->cb, tm);
}

void
choldecomp_h2matrix(ph2matrix A, pclusteroperator rwf,
		    pclusteroperator cwf, ph2matrix L,
		    pclusteroperator rwflow, pclusteroperator cwflow,
		    ptruncmode tm, real tol)
{
  uint      sons = A->rb->sons;
  pccluster t = A->rb->t;

  pclusteroperator rw, cw, rwlow, cwlow;
  uint      i, j, l;

  assert(t == A->cb->t);

  if (A->son != 0) {

    rw = identify_son_clusterweight_clusteroperator(rwf, t);
    cw = identify_son_clusterweight_clusteroperator(cwf, t);
    rwlow = identify_son_clusterweight_clusteroperator(rwflow, t);
    cwlow = identify_son_clusterweight_clusteroperator(cwflow, t);

    for (i = 0; i < sons; i++) {

      /* compute the i-th diagonal block */
      choldecomp_h2matrix(A->son[i + i * sons], rw, cw, L->son[i + i * sons],
			  rwlow, cwlow, tm, tol);

      for (j = i + 1; j < sons; j++) {

	/* compute the i-th column */
	lowersolve_h2matrix(false, false, L->son[i + i * sons], true,
			    A->son[j + i * sons], rw, cw,
			    L->son[j + i * sons], rwlow, cwlow, tm, tol);
      }

      /* update the lower right block */
      for (j = i + 1; j < sons; j++) {
	for (l = i + 1; l <= j; l++) {
	  addmul_h2matrix(-1.0, L->son[j + i * sons], true,
			  L->son[l + i * sons], A->son[j + l * sons], rw, cw,
			  tm, tol);
	}
      }
    }
    update_clusterbasis(A->rb);
    update_clusterbasis(A->cb);
    update_clusterbasis(L->rb);
    update_clusterbasis(L->cb);
  }
  else {
    assert(A->f != 0);
    choldecomp_amatrix(A->f);
    copy_lower_amatrix(A->f, false, L->f);
  }
}

void
cholsolve_h2matrix_avector(pch2matrix h2, pavector x)
{
  pavector  xp;
  avector   xtmp;
  uint      i, ip;

  assert(x->dim == h2->cb->t->size);

  /* Permutation of x */
  xp = init_avector(&xtmp, x->dim);
  for (i = 0; i < xp->dim; i++) {
    ip = h2->cb->t->idx[i];
    assert(ip < x->dim);
    xp->v[i] = x->v[ip];
  }

  lowersolve_h2matrix_avector(false, false, h2, xp);
  uppersolve_h2matrix_avector(false, true, h2, xp);

  /* Reverse permutation of x */
  for (i = 0; i < xp->dim; i++) {
    ip = h2->rb->t->idx[i];
    assert(ip < x->dim);
    x->v[ip] = xp->v[i];
  }
  uninit_avector(xp);
}
