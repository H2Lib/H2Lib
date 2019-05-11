#include "h2update.h"

#include <stdio.h>

#include "factorizations.h"
#include "h2compression.h"
#include "basic.h"

#include "laplacebem2d.h"
#include "laplacebem3d.h"

real
norm2_rkupdate_uniform(puniform u, uint k)
{
  pclusterbasis rb = u->rb;
  pclusterbasis cb = u->cb;

  amatrix   tmp1, tmp2, tmp3, tmp4;
  pamatrix  uz, zuz, Z1, Z2;
  real      norm = 0.0;

  /* zuz = rb->Z * u->S * cb->Z */
  if (cb->Z) {
    uz = init_zero_amatrix(&tmp1, u->S.rows, cb->Z->rows);
    Z1 = init_sub_amatrix(&tmp3, cb->Z, cb->Z->rows, 0, u->S.cols, 0);
    addmul_amatrix(1.0, false, &u->S, true, Z1, uz);
    uninit_amatrix(Z1);
  }
  else
    uz = &u->S;

  if (rb->Z) {
    zuz = init_zero_amatrix(&tmp2, rb->Z->rows, uz->cols);
    Z2 = init_sub_amatrix(&tmp4, rb->Z, rb->Z->rows, 0, u->S.rows, 0);
    addmul_amatrix(1.0, false, Z2, false, uz, zuz);
    uninit_amatrix(Z2);
  }
  else
    zuz = uz;

  if (rb->Z && cb->Z) {
    Z1 = init_sub_amatrix(&tmp3, cb->Z, cb->Z->rows, 0, k, u->S.cols);
    Z2 = init_sub_amatrix(&tmp4, rb->Z, rb->Z->rows, 0, k, u->S.rows);
    addmul_amatrix(1.0, false, Z2, true, Z1, zuz);
    uninit_amatrix(Z1);
    uninit_amatrix(Z2);
  }

  /* norm of zuz via vector iteration */
  norm = norm2_amatrix(zuz);

  if (rb->Z)
    uninit_amatrix(zuz);
  if (cb->Z)
    uninit_amatrix(uz);

  return REAL_SQRT(norm);
}

real
normfrob_rkupdate_uniform(puniform u, uint k)
{
  pclusterbasis rb = u->rb;
  pclusterbasis cb = u->cb;

  amatrix   tmp1, tmp2, tmp3, tmp4;
  pamatrix  uz, zuz, Z1, Z2;
  real      norm = 0.0;

  /* zuz = rb->Z * u->S * cb->Z */
  if (cb->Z) {
    uz = init_zero_amatrix(&tmp1, u->S.rows, cb->Z->rows);
    Z1 = init_sub_amatrix(&tmp3, cb->Z, cb->Z->rows, 0, u->S.cols, 0);
    addmul_amatrix(1.0, false, &u->S, true, Z1, uz);
    uninit_amatrix(Z1);
  }
  else
    uz = &u->S;

  if (rb->Z) {
    zuz = init_zero_amatrix(&tmp2, rb->Z->rows, uz->cols);
    Z2 = init_sub_amatrix(&tmp4, rb->Z, rb->Z->rows, 0, u->S.rows, 0);
    addmul_amatrix(1.0, false, Z2, false, uz, zuz);
    uninit_amatrix(Z2);
  }
  else
    zuz = uz;

  if (rb->Z && cb->Z) {
    Z1 = init_sub_amatrix(&tmp3, cb->Z, cb->Z->rows, 0, k, u->S.cols);
    Z2 = init_sub_amatrix(&tmp4, rb->Z, rb->Z->rows, 0, k, u->S.rows);
    addmul_amatrix(1.0, false, Z2, true, Z1, zuz);
    uninit_amatrix(Z1);
    uninit_amatrix(Z2);
  }

  /* norm of zuz */
  norm = normfrob_amatrix(zuz);

  if (rb->Z)
    uninit_amatrix(zuz);
  if (cb->Z)
    uninit_amatrix(uz);

  return norm;
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
/* compute the R of QR decomposition of the extended clusterbasis (V A) */
static void
orthoweight_rkupdate_clusterbasis(pclusterbasis cb, pamatrix A)
{
  uint      sons = cb->sons;
  pclusterbasis *son = cb->son;
  uint      k = A->cols;

  amatrix   tmp1, tmp2, tmp4;
  avector   tmp3;
  pamatrix  Vhat, Vhat1, A1, R;
  uint      m, off, roff;
  pavector  tau;
  uint      i, refl;

  assert(cb->t->size == A->rows);

  if (sons > 0) {
    /* compute weights for sons */
    m = 0;
    roff = 0;
    for (i = 0; i < sons; i++) {
      A1 = init_sub_amatrix(&tmp2, A, son[i]->t->size, roff, k, 0);
      orthoweight_rkupdate_clusterbasis(son[i], A1);
      m += son[i]->Z->rows;
      roff += son[i]->t->size;

      uninit_amatrix(A1);
    }
    assert(m <= roff);
    assert(roff == cb->t->size);

    /* compute Vhat for current cluster */
    Vhat = init_amatrix(&tmp1, m, cb->k + k);
    off = 0;
    for (i = 0; i < sons; i++) {
      Vhat1 = init_sub_amatrix(&tmp2, Vhat, son[i]->Z->rows, off, cb->k, 0);
      R =
	init_sub_amatrix(&tmp4, son[i]->Z, son[i]->Z->rows, 0, cb->son[i]->k,
			 0);
      clear_amatrix(Vhat1);
      addmul_amatrix(1.0, false, R, false, &cb->son[i]->E, Vhat1);
      uninit_amatrix(Vhat1);
      uninit_amatrix(R);

      assert(son[i]->Z->cols == cb->son[i]->k + k);
      Vhat1 = init_sub_amatrix(&tmp2, Vhat, son[i]->Z->rows, off, k, cb->k);
      R = init_sub_amatrix(&tmp4, son[i]->Z, son[i]->Z->rows, 0, k,
			   cb->son[i]->k);
      copy_amatrix(false, R, Vhat1);
      uninit_amatrix(Vhat1);
      uninit_amatrix(R);

      off += son[i]->Z->rows;
    }
    assert(off == m);
  }
  /* compute Vhat for current leaf cluster */
  else {
    assert(sons == 0);
    m = cb->t->size;

    Vhat = init_amatrix(&tmp1, m, cb->k + k);

    Vhat1 = init_sub_amatrix(&tmp2, Vhat, m, 0, cb->k, 0);
    copy_amatrix(false, &cb->V, Vhat1);
    uninit_amatrix(Vhat1);

    Vhat1 = init_sub_amatrix(&tmp2, Vhat, m, 0, k, cb->k);
    copy_amatrix(false, A, Vhat1);
    uninit_amatrix(Vhat1);
  }

  /* compute weight of current cluster */
  refl = UINT_MIN(m, cb->k + k);
  tau = init_avector(&tmp3, refl);
  qrdecomp_amatrix(Vhat, tau);
  assert(cb->Z == NULL);
  cb->Z = new_amatrix(refl, cb->k + k);
  copy_upper_amatrix(Vhat, false, cb->Z);
  uninit_avector(tau);
  uninit_amatrix(Vhat);
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
/* computes the totalweights for all sons of the extended row clusterbasis */
static void
rowweight_rkupdate_clusteroperator(pclusterbasis rb,
				   pclusteroperator rw, int k, ptruncmode tm)
{
  uint      sons = rw->sons;

  pamatrix  Yhat, Yhat1;
  pamatrix  Z, Z1;		/* Z = Z or Z = R_{si} */
  amatrix   tmp1, tmp2, tmp3;
  pavector  tau;
  avector   tmp4;
  pclusterbasis son;
  puniform  u;
  real      zeta_age, norm, alpha;
  uint      rows, cols, refl;	/* size of Yhat */
  uint      off;
  uint      i;

  assert(rb->t == rw->t);

  zeta_age = (tm ? tm->zeta_age : 1.0);

  if (sons > 0) {
    for (i = 0; i < sons; i++) {
      son = rb->son[i];

      /* rows of Yhat */
      rows = rw->krow;
      u = son->rlist;
      while (u != NULL) {
	Z = u->cb->Z;
	/* u is a subblock of AB* */
	if (Z != NULL) {
	  rows += Z->rows;
	}
	/* u is no subblock of AB* */
	else {
	  rows += u->S.cols;
	}
	u = u->rnext;
      }
      /* cols of Yhat */
      cols = son->k + k;
      Yhat = init_amatrix(&tmp1, rows, cols);
      /* compute Yhat */
      assert(rw->kcol == son->E.cols + k);
      /* columns associated with V */
      Yhat1 = init_sub_amatrix(&tmp2, Yhat, rw->krow, 0, son->k, 0);
      Z1 = init_sub_amatrix(&tmp3, &rw->C, rw->krow, 0, son->E.cols, 0);
      clear_amatrix(Yhat1);
      addmul_amatrix(zeta_age, false, Z1, true, &son->E, Yhat1);
      uninit_amatrix(Yhat1);
      uninit_amatrix(Z1);

      /* columns associated with A */
      Yhat1 = init_sub_amatrix(&tmp2, Yhat, rw->krow, 0, k, son->k);
      Z1 = init_sub_amatrix(&tmp3, &rw->C, rw->krow, 0, k, son->E.cols);
      copy_amatrix(false, Z1, Yhat1);
      scale_amatrix(zeta_age, Yhat1);
      uninit_amatrix(Yhat1);
      uninit_amatrix(Z1);

      off = rw->krow;
      u = son->rlist;
      while (u) {
	/* Compute block weight if required */
	alpha = 1.0;
	if (tm && tm->blocks) {
	  if (tm->frobenius)
	    norm = normfrob_rkupdate_uniform(u, k);
	  else
	    norm = norm2_rkupdate_uniform(u, k);

	  alpha = (norm > 0.0 ? 1.0 / norm : 1.0);
	}
	Z = u->cb->Z;
	/* u is a subblock of AB* */
	if (Z != NULL) {
	  assert(Z->cols == u->S.cols + k);
	  /* columns associated with V */
	  Yhat1 = init_sub_amatrix(&tmp2, Yhat, Z->rows, off, son->k, 0);
	  Z1 = init_sub_amatrix(&tmp3, Z, Z->rows, 0, u->S.cols, 0);
	  clear_amatrix(Yhat1);
	  addmul_amatrix(alpha, false, Z1, true, &u->S, Yhat1);
	  uninit_amatrix(Yhat1);
	  uninit_amatrix(Z1);

	  /* columns associated with A */
	  Yhat1 = init_sub_amatrix(&tmp2, Yhat, Z->rows, off, k, son->k);
	  Z1 = init_sub_amatrix(&tmp3, Z, Z->rows, 0, k, u->S.cols);
	  copy_amatrix(false, Z1, Yhat1);
	  scale_amatrix(alpha, Yhat1);
	  uninit_amatrix(Yhat1);
	  uninit_amatrix(Z1);

	  off += Z->rows;
	}
	/* u is no subblock of AB* */
	else {
	  /* columns associated with V */
	  Yhat1 = init_sub_amatrix(&tmp2, Yhat, u->S.cols, off, son->k, 0);
	  copy_amatrix(true, &u->S, Yhat1);
	  scale_amatrix(alpha, Yhat1);
	  uninit_amatrix(Yhat1);

	  /* columns associated with A */
	  Yhat1 = init_sub_amatrix(&tmp2, Yhat, u->S.cols, off, k, son->k);
	  clear_amatrix(Yhat1);
	  uninit_amatrix(Yhat1);

	  off += u->S.cols;
	}
	u = u->rnext;
      }
      assert(off == rows);

      /* compute the weight */
      refl = UINT_MIN(rows, cols);
      tau = init_avector(&tmp4, refl);
      qrdecomp_amatrix(Yhat, tau);
      resize_clusteroperator(rw->son[i], refl, cols);
      copy_upper_amatrix(Yhat, false, &rw->son[i]->C);
      uninit_avector(tau);
      uninit_amatrix(Yhat);

      /* compute the weight for the son */
      rowweight_rkupdate_clusteroperator(son, rw->son[i], k, tm);
    }
  }
  else {
    assert(sons == 0);
  }
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
/* computes the totalweights for all sons of the extended col clusterbasis */
static void
colweight_rkupdate_clusteroperator(pclusterbasis cb,
				   pclusteroperator cw, int k, ptruncmode tm)
{
  uint      sons = cw->sons;

  pamatrix  Yhat, Yhat1;
  pamatrix  Z, Z1;		/* Z = Z or Z = R_{si} */
  amatrix   tmp1, tmp2, tmp3;
  pavector  tau;
  avector   tmp4;
  pclusterbasis son;
  puniform  u;
  real      zeta_age, norm, alpha;
  uint      rows, cols, refl;	/* size of Yhat */
  uint      off;
  uint      i;

  assert(cb->t == cw->t);

  zeta_age = (tm ? tm->zeta_age : 1.0);

  if (sons > 0) {
    for (i = 0; i < sons; i++) {
      son = cb->son[i];

      /* rows of Yhat */
      rows = cw->krow;
      u = son->clist;
      while (u != NULL) {
	Z = u->rb->Z;
	/* u is a subblock of AB* */
	if (Z != NULL) {
	  rows += Z->rows;
	}
	/* u is no subblock of AB* */
	else {
	  rows += u->S.rows;
	}
	u = u->cnext;
      }
      /* cols of Yhat */
      cols = son->k + k;
      Yhat = init_amatrix(&tmp1, rows, cols);

      /* compute Yhat */
      assert(cw->kcol == son->E.cols + k);
      /* columns associated with W */
      Yhat1 = init_sub_amatrix(&tmp2, Yhat, cw->krow, 0, son->k, 0);
      Z1 = init_sub_amatrix(&tmp3, &cw->C, cw->krow, 0, son->E.cols, 0);
      clear_amatrix(Yhat1);
      addmul_amatrix(zeta_age, false, Z1, true, &son->E, Yhat1);
      uninit_amatrix(Yhat1);
      uninit_amatrix(Z1);

      /* columns associated with B */
      Yhat1 = init_sub_amatrix(&tmp2, Yhat, cw->krow, 0, k, son->k);
      Z1 = init_sub_amatrix(&tmp3, &cw->C, cw->krow, 0, k, son->E.cols);
      copy_amatrix(false, Z1, Yhat1);
      scale_amatrix(zeta_age, Yhat1);
      uninit_amatrix(Yhat1);
      uninit_amatrix(Z1);

      off = cw->krow;
      u = son->clist;
      while (u != NULL) {
	/* Compute block weight if required */
	alpha = 1.0;
	if (tm && tm->blocks) {
	  if (tm->frobenius)
	    norm = normfrob_rkupdate_uniform(u, k);
	  else
	    norm = norm2_rkupdate_uniform(u, k);

	  alpha = (norm > 0.0 ? 1.0 / norm : 1.0);
	}
	Z = u->rb->Z;
	/* u is a subblock of AB* */
	if (Z != NULL) {
	  assert(Z->cols == u->S.rows + k);
	  /* columns associated with V */
	  Yhat1 = init_sub_amatrix(&tmp2, Yhat, Z->rows, off, son->k, 0);
	  Z1 = init_sub_amatrix(&tmp3, Z, Z->rows, 0, u->S.rows, 0);
	  clear_amatrix(Yhat1);
	  addmul_amatrix(alpha, false, Z1, false, &u->S, Yhat1);
	  uninit_amatrix(Yhat1);
	  uninit_amatrix(Z1);

	  /* columns associated with A */
	  Yhat1 = init_sub_amatrix(&tmp2, Yhat, Z->rows, off, k, son->k);
	  Z1 = init_sub_amatrix(&tmp3, Z, Z->rows, 0, k, u->S.rows);
	  copy_amatrix(false, Z1, Yhat1);
	  scale_amatrix(alpha, Yhat1);
	  uninit_amatrix(Yhat1);
	  uninit_amatrix(Z1);

	  off += Z->rows;
	}
	/* u is no subblock of AB* */
	else {
	  /* columns associated with V */
	  Yhat1 = init_sub_amatrix(&tmp2, Yhat, u->S.rows, off, son->k, 0);
	  copy_amatrix(false, &u->S, Yhat1);
	  scale_amatrix(alpha, Yhat1);
	  uninit_amatrix(Yhat1);

	  /* columns associated with A */
	  Yhat1 = init_sub_amatrix(&tmp2, Yhat, u->S.rows, off, k, son->k);
	  clear_amatrix(Yhat1);
	  uninit_amatrix(Yhat1);

	  off += u->S.rows;
	}
	u = u->cnext;
      }
      assert(off == rows);

      /* compute the weight */
      refl = UINT_MIN(rows, cols);
      tau = init_avector(&tmp4, refl);
      qrdecomp_amatrix(Yhat, tau);
      resize_clusteroperator(cw->son[i], refl, cols);
      copy_upper_amatrix(Yhat, false, &cw->son[i]->C);
      uninit_avector(tau);
      uninit_amatrix(Yhat);

      /* compute the weight for the son */
      colweight_rkupdate_clusteroperator(son, cw->son[i], k, tm);
    }
  }
  else {
    assert(sons == 0);
  }
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
/* compute adaptive clusterbasis for extended clusterbasis (V A) */
static void
truncate_rkupdate_clusterbasis(pclusterbasis cb, pamatrix A,
			       pclusteroperator cw, pctruncmode tm, real eps)
{
  amatrix   tmp1, tmp2, tmp3;
  realavector tmp4;
  pamatrix  Vhat, Vhat1, VhatZ, Q, Q1, A1, C1, E1;
  prealavector sigma;
  pclusteroperator cw1;
  real      zeta_level;
  uint      i, off, m, k;

  zeta_level = (tm ? tm->zeta_level : 1.0);

  Vhat = 0;
  if (cb->sons == 0) {
    /* In son clusters, we have Vhat = (V A) */
    m = cb->t->size;
    Vhat = init_amatrix(&tmp1, cb->t->size, cb->k + A->cols);

    Vhat1 = init_sub_amatrix(&tmp2, Vhat, cb->t->size, 0, cb->k, 0);
    copy_amatrix(false, &cb->V, Vhat1);
    uninit_amatrix(Vhat1);
    Vhat1 = init_sub_amatrix(&tmp2, Vhat, cb->t->size, 0, A->cols, cb->k);
    copy_amatrix(false, A, Vhat1);
    uninit_amatrix(Vhat1);
  }
  else {
    /* Save original transfer matrices */
    E1 = allocmem(sizeof(amatrix) * cb->sons);

    /* Compute cluster bases for son clusters recursively */
    m = 0;
    off = 0;
    assert(cb->sons == cw->sons);
    for (i = 0; i < cb->sons; i++) {
      cw1 = cw->son[i];

      init_amatrix(E1+i, cb->son[i]->E.rows, cb->son[i]->E.cols);
      copy_amatrix(false, &cb->son[i]->E, E1+i);

      A1 = init_sub_amatrix(&tmp2, A, cb->son[i]->t->size, off, A->cols, 0);
      truncate_rkupdate_clusterbasis(cb->son[i], A1, cw1, tm,
				     eps * zeta_level);
      uninit_amatrix(A1);

      off += cb->son[i]->t->size;
      m += cb->son[i]->k;
    }
    assert(off == A->rows);

    /* compute Vhat */
    Vhat = init_amatrix(&tmp1, m, cb->k + A->cols);
    /* Blocks of Vhat are (R_{t'|...}E_{t'}   (R_{t'|...}) for sons t' */
    off = 0;
    for (i = 0; i < cb->sons; i++) {
      Vhat1 = init_sub_amatrix(&tmp2, Vhat, cb->son[i]->k, off, cb->k, 0);
      C1 = init_sub_amatrix(&tmp3, &cw->son[i]->C, cb->son[i]->k, 0,
			    E1[i].rows, 0);
      clear_amatrix(Vhat1);
      addmul_amatrix(1.0, false, C1, false, E1+i, Vhat1);
      uninit_amatrix(Vhat1);
      uninit_amatrix(C1);

      Vhat1 =
	init_sub_amatrix(&tmp2, Vhat, cb->son[i]->k, off, A->cols, cb->k);
      C1 =
	init_sub_amatrix(&tmp3, &cw->son[i]->C, cb->son[i]->k, 0,
			 A->cols, E1[i].rows);
      copy_amatrix(false, C1, Vhat1);
      uninit_amatrix(Vhat1);
      uninit_amatrix(C1);

      uninit_amatrix(E1+i);

      off += cb->son[i]->k;
    }
    assert(off == m);
    
    freemem(E1);
  }

  /* Multiply by weight matrix */
  VhatZ = init_amatrix(&tmp2, Vhat->rows, cw->krow);
  clear_amatrix(VhatZ);
  addmul_amatrix(1.0, false, Vhat, true, &cw->C, VhatZ);

  /* Compute singular value decomposition of VhatZ */
  k = UINT_MIN(VhatZ->rows, VhatZ->cols);
  Q = init_amatrix(&tmp3, VhatZ->rows, VhatZ->cols);
  sigma = init_realavector(&tmp4, k);
  svd_amatrix(VhatZ, sigma, Q, 0);
  uninit_amatrix(VhatZ);

  /* Find appropriate rank */
  k = findrank_truncmode(tm, eps, sigma);
  uninit_realavector(sigma);

  /* Set rank of new cluster basis */
  resize_clusterbasis(cb, k);

  if (cb->sons == 0) {
    assert(cb->V.rows == Q->rows);
    assert(cb->V.cols == k);

    /* Copy cluster basis matrix to new cluster basis */
    Q1 = init_sub_amatrix(&tmp2, Q, Q->rows, 0, k, 0);
    copy_amatrix(false, Q1, &cb->V);
    uninit_amatrix(Q1);
  }
  else {
    off = 0;
    for (i = 0; i < cb->sons; i++) {
      assert(cb->son[i]->E.rows == cb->son[i]->k);
      assert(cb->son[i]->E.cols == k);

      /* Copy transfer matrix to new cluster basis */
      Q1 = init_sub_amatrix(&tmp2, Q, cb->son[i]->k, off, k, 0);
      copy_amatrix(false, Q1, &cb->son[i]->E);
      uninit_amatrix(Q1);

      off += cb->son[i]->k;
    }
    assert(off == m);
  }

  /* Compute transformation from old to new basis */
  Q1 = init_sub_amatrix(&tmp2, Q, Q->rows, 0, k, 0);
  resize_clusteroperator(cw, k, Vhat->cols);
  clear_amatrix(&cw->C);
  addmul_amatrix(1.0, true, Q1, false, Vhat, &cw->C);
  uninit_amatrix(Q1);

  /* Clean up */
  uninit_amatrix(Q);
  uninit_amatrix(Vhat);

  update_clusterbasis(cb);
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
/* adapts the coupling matrices of subblocks */
static void
rkupdate_inside_h2matrix(ph2matrix Gh2, pamatrix A, pamatrix B,
			 pclusteroperator rw, pclusteroperator cw)
{
  uint      rsons = Gh2->rsons;
  uint      csons = Gh2->csons;

  pclusteroperator rw1, cw1;
  amatrix   tmp1, tmp2, tmp3, tmp4;
  pamatrix  A1, B1, S, S1, S2;
  uint      roff, coff;
  uint      k, rsize, csize;
  uint      i, j;

  assert(rw->t == Gh2->rb->t);
  assert(cw->t == Gh2->cb->t);
  assert(rw->t->size == A->rows);
  assert(cw->t->size == B->rows);

  k = A->cols;
  assert(k == B->cols);

  /* update the sons recursively */
  if (Gh2->son) {
    roff = 0;
    for (i = 0; i < rsons; i++) {
      if (rw->sons == 0)
	rw1 = rw;
      else
	rw1 = rw->son[i];
      rsize = rw1->t->size;
      A1 = init_sub_amatrix(&tmp1, A, rsize, roff, k, 0);

      coff = 0;
      for (j = 0; j < csons; j++) {
	if (cw->sons == 0)
	  cw1 = cw;
	else
	  cw1 = cw->son[j];
	csize = cw1->t->size;
	B1 = init_sub_amatrix(&tmp2, B, csize, coff, k, 0);
	rkupdate_inside_h2matrix(Gh2->son[i + j * rsons], A1, B1, rw1, cw1);

	coff += csize;
	uninit_amatrix(B1);
      }
      assert(coff == cw->t->size);

      roff += rsize;
      uninit_amatrix(A1);
    }
    assert(roff == rw->t->size);
  }
  /* in inadmissible leafs just add AB^T */
  else if (Gh2->f) {
    addmul_amatrix(1.0, false, A, true, B, Gh2->f);
  }
  else if (Gh2->u) {
    S = &Gh2->u->S;
    /* left projection of old coupling matrix */
    S2 = init_amatrix(&tmp4, Gh2->rb->k, S->cols);
    A1 = init_sub_amatrix(&tmp1, &rw->C, rw->krow, 0, S->rows, 0);
    clear_amatrix(S2);
    addmul_amatrix(1.0, false, A1, false, S, S2);
    uninit_amatrix(A1);
    /* right projection of old coupling matrix */
    S1 = init_amatrix(&tmp3, Gh2->rb->k, Gh2->cb->k);
    B1 = init_sub_amatrix(&tmp2, &cw->C, cw->krow, 0, S->cols, 0);
    clear_amatrix(S1);
    addmul_amatrix(1.0, false, S2, true, B1, S1);
    uninit_amatrix(B1);
    uninit_amatrix(S2);
    /* add the coupling part corresponding to AB^T */
    assert(rw->kcol == S->rows + k);
    A1 = init_sub_amatrix(&tmp1, &rw->C, rw->krow, 0, k, S->rows);
    assert(cw->kcol == S->cols + k);
    B1 = init_sub_amatrix(&tmp2, &cw->C, cw->krow, 0, k, S->cols);
    addmul_amatrix(1.0, false, A1, true, B1, S1);
    uninit_amatrix(A1);
    uninit_amatrix(B1);

    uninit_amatrix(S);
    Gh2->u->S = *S1;
    assert(Gh2->u->rb->k == Gh2->u->S.rows);
    assert(Gh2->u->cb->k == Gh2->u->S.cols);
  }
  else {
    /* Gh2 is admissible zero block */
    Gh2->u = new_uniform(Gh2->rb, Gh2->cb);
    S = &Gh2->u->S;
    clear_amatrix(S);
    /* add the coupling part corresponding to AB^T */
    A1 = init_sub_amatrix(&tmp1, &rw->C, rw->krow, 0, k, rw->kcol - k);
    B1 = init_sub_amatrix(&tmp2, &cw->C, cw->krow, 0, k, cw->kcol - k);
    addmul_amatrix(1.0, false, A1, true, B1, S);
    uninit_amatrix(A1);
    uninit_amatrix(B1);
  }
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
/* adapts the coupling matrices outside the block of AB* to new row clusterbasis */
static void
rkupdate_rowout_h2matrix(pclusterbasis rb, pclusteroperator rw)
{
  puniform  u;
  amatrix   tmp1, tmp2;
  pamatrix  S, S1, R1;
  uint      i;

  assert(rb->t == rw->t);
  assert(rb->k == rw->krow);

  u = rb->rlist;
  while (u) {
    assert(u->rb == rb);
    /* u is no subblock of AB* */
    if (u->cb->Z == 0) {
      /* left projection of coupling matrix */
      S = &u->S;
      S1 = init_amatrix(&tmp1, rb->k, S->cols);
      R1 = init_sub_amatrix(&tmp2, &rw->C, rw->krow, 0, S->rows, 0);
      clear_amatrix(S1);
      addmul_amatrix(1.0, false, R1, false, S, S1);
      uninit_amatrix(R1);

      uninit_amatrix(S);
      u->S = *S1;
    }
    assert(u->rb->k == u->S.rows);
    assert(u->cb->k == u->S.cols);
    u = u->rnext;
  }

  /* update the sons recursively */
  for (i = 0; i < rb->sons; i++) {
    rkupdate_rowout_h2matrix(rb->son[i], rw->son[i]);
  }
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
/* adapts the coupling matrices outside the block of AB* to new col clusterbasis */
static void
rkupdate_colout_h2matrix(pclusterbasis cb, pclusteroperator cw)
{
  puniform  u;
  amatrix   tmp1, tmp2;
  pamatrix  S, S1, R1;
  uint      i;

  assert(cb->t == cw->t);
  assert(cb->k == cw->krow);

  u = cb->clist;
  while (u) {
    assert(u->cb == cb);
    /* u is no subblock of AB* */
    if (u->rb->Z == 0) {
      /* right projection of coupling matrix */
      S = &u->S;
      S1 = init_amatrix(&tmp1, S->rows, cb->k);
      R1 = init_sub_amatrix(&tmp2, &cw->C, cw->krow, 0, S->cols, 0);
      clear_amatrix(S1);
      addmul_amatrix(1.0, false, S, true, R1, S1);
      uninit_amatrix(R1);

      uninit_amatrix(S);
      u->S = *S1;
    }
    assert(u->rb->k == u->S.rows);
    assert(u->cb->k == u->S.cols);
    u = u->cnext;
  }

  /* update the sons recursively */
  for (i = 0; i < cb->sons; i++) {
    rkupdate_colout_h2matrix(cb->son[i], cw->son[i]);
  }
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
/* compute the new total weights for the row cluster */
static void
totalweights_row_clusteroperator(pclusterbasis rb,
				 pclusteroperator rw, ptruncmode tm)
{
  uint      sons = rw->sons;

  pamatrix  Yhat, Yhat1;
  amatrix   tmp1, tmp2;
  pavector  tau;
  avector   tmp3;
  pclusterbasis son;
  puniform  u;
  real      zeta_age, norm, alpha;
  uint      rows, cols, refl;	/* size of Yhat */
  uint      off;
  uint      i;

  zeta_age = (tm ? tm->zeta_age : 1.0);

  if (rw->son != NULL) {
    for (i = 0; i < sons; i++) {
      son = rb->son[i];

      /* rows of Yhat */
      rows = rw->krow;
      u = son->rlist;
      while (u != NULL) {
	rows += u->S.cols;
	u = u->rnext;
      }
      /* cols of Yhat */
      cols = son->k;
      Yhat = init_amatrix(&tmp1, rows, cols);

      /* transfer the weight of the father */
      Yhat1 = init_sub_amatrix(&tmp2, Yhat, rw->krow, 0, son->k, 0);
      clear_amatrix(Yhat1);
      addmul_amatrix(zeta_age, false, &rw->C, true, &son->E, Yhat1);
      uninit_amatrix(Yhat1);

      off = rw->krow;
      u = son->rlist;
      while (u) {
	/* Compute block weight if required */
	alpha = 1.0;
	if (tm && tm->blocks) {
	  if (tm->frobenius)
	    norm = normfrob_rkupdate_uniform(u, 0);
	  else
	    norm = norm2_rkupdate_uniform(u, 0);

	  alpha = (norm > 0.0 ? 1.0 / norm : 1.0);
	}
	/* coupling matrices for admissible blocks */
	Yhat1 = init_sub_amatrix(&tmp2, Yhat, u->S.cols, off, son->k, 0);
	copy_amatrix(true, &u->S, Yhat1);
	scale_amatrix(alpha, Yhat1);
	uninit_amatrix(Yhat1);

	off += u->S.cols;
	u = u->rnext;
      }
      assert(off == rows);

      /* compute the weight for current cluster */
      refl = UINT_MIN(rows, cols);
      tau = init_avector(&tmp3, refl);
      qrdecomp_amatrix(Yhat, tau);
      resize_clusteroperator(rw->son[i], refl, cols);
      copy_upper_amatrix(Yhat, false, &rw->son[i]->C);
      uninit_avector(tau);
      uninit_amatrix(Yhat);

      /* compute the weights for the son recursively */
      totalweights_row_clusteroperator(son, rw->son[i], tm);
    }
  }
}

/* compute the new total weights for the column cluster */
static void
totalweights_col_clusteroperator(pclusterbasis cb,
				 pclusteroperator cw, ptruncmode tm)
{
  uint      sons = cw->sons;

  pamatrix  Yhat, Yhat1;
  amatrix   tmp1, tmp2;
  pavector  tau;
  avector   tmp3;
  pclusterbasis son;
  puniform  u;
  real      zeta_age, norm, alpha;
  uint      rows, cols, refl;	/* size of Yhat */
  uint      off;
  uint      i;

  zeta_age = (tm ? tm->zeta_age : 1.0);

  if (cw->son != NULL) {
    for (i = 0; i < sons; i++) {
      son = cb->son[i];

      /* rows of Yhat */
      rows = cw->krow;
      u = son->clist;
      while (u != NULL) {
	rows += u->S.rows;
	u = u->cnext;
      }
      /* cols of Yhat */
      cols = son->k;
      Yhat = init_amatrix(&tmp1, rows, cols);

      /* transfer the weight of the father */
      Yhat1 = init_sub_amatrix(&tmp2, Yhat, cw->krow, 0, son->k, 0);
      clear_amatrix(Yhat1);
      addmul_amatrix(zeta_age, false, &cw->C, true, &son->E, Yhat1);
      uninit_amatrix(Yhat1);

      off = cw->krow;
      u = son->clist;
      while (u) {
	/* Compute block weight if required */
	alpha = 1.0;
	if (tm && tm->blocks) {
	  if (tm->frobenius)
	    norm = normfrob_rkupdate_uniform(u, 0);
	  else
	    norm = norm2_rkupdate_uniform(u, 0);

	  alpha = (norm > 0.0 ? 1.0 / norm : 1.0);
	}
	/* coupling matrices for admissible blocks */
	Yhat1 = init_sub_amatrix(&tmp2, Yhat, u->S.rows, off, son->k, 0);
	scale_amatrix(alpha, Yhat1);
	copy_amatrix(false, &u->S, Yhat1);
	uninit_amatrix(Yhat1);

	off += u->S.rows;
	u = u->cnext;
      }
      assert(off == rows);

      /* compute the weight current cluster */
      refl = UINT_MIN(rows, cols);
      tau = init_avector(&tmp3, refl);
      qrdecomp_amatrix(Yhat, tau);
      resize_clusteroperator(cw->son[i], refl, cols);
      copy_upper_amatrix(Yhat, false, &cw->son[i]->C);
      uninit_avector(tau);
      uninit_amatrix(Yhat);

      /* compute the weights for the son recursively */
      totalweights_col_clusteroperator(son, cw->son[i], tm);
    }
  }
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
/* adds a uniform block to the h2matrix structure */
static void
rkupdate_adduniform_h2matrix(ph2matrix h2)
{
  uint      i;

  if (h2->son)
    for (i = 0; i < h2->rsons * h2->csons; i++)
      rkupdate_adduniform_h2matrix(h2->son[i]);

  if (!h2->u && !h2->f && !h2->son) {
    h2->u = new_uniform(h2->rb, h2->cb);
    clear_amatrix(&h2->u->S);
  }
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
/* Gh2 = Gh2 + R */
void
rkupdate_h2matrix(prkmatrix R, ph2matrix Gh2, pclusteroperator rwf,
		  pclusteroperator cwf, ptruncmode tm, real eps)
{
  pclusterbasis rb = Gh2->rb;
  pclusterbasis cb = Gh2->cb;

  pclusteroperator rw, cw;
  amatrix   tmp1, tmp2, tmp3;
  pamatrix  Yhat, Yhat1;
  pamatrix  Z, Z1, E1;
  avector   tmp4;
  pavector  tau;
  puniform  u;
  real      zeta_age, norm, alpha;
  uint      i, k;
  uint      rows, cols, refl;
  uint      off;

  assert(rwf->son);
  assert(cwf->son);

  zeta_age = (tm ? tm->zeta_age : 1.0);

  rkupdate_adduniform_h2matrix(Gh2);

  /* calculate orthogonal weights for (V_t A) */
  orthoweight_rkupdate_clusterbasis(rb, &R->A);

  /* calculate orthogonal weights for (W_s B) */
  orthoweight_rkupdate_clusterbasis(cb, &R->B);

  /* search for row weight of the active block */
  rw = NULL;
  for (i = 0; i < rwf->sons; i++) {
    if (rwf->son[i]->t == rb->t) {
      rw = rwf->son[i];
      break;
    }
  }
  assert(rw != NULL);

  /* search for col weight of the active block */
  cw = NULL;
  for (i = 0; i < cwf->sons; i++) {
    if (cwf->son[i]->t == cb->t) {
      cw = cwf->son[i];
      break;
    }
  }
  assert(cw != NULL);

  /* -------------------------------------------------- */
  /* compute total weight for current row cluster (V A) */

  /* rows of Yhat */
  rows = rwf->krow;
  u = rb->rlist;
  while (u != NULL) {
    Z = u->cb->Z;
    /* u is a subblock of AB* */
    if (Z != NULL) {
      rows += Z->rows;
    }
    /* u is no subblock of AB* */
    else {
      rows += u->S.cols;
    }
    u = u->rnext;
  }
  /* cols of Yhat */
  cols = rb->k + R->k;
  Yhat = init_amatrix(&tmp1, rows, cols);

  /* columns associated with V */
  Yhat1 = init_sub_amatrix(&tmp2, Yhat, rwf->krow, 0, rb->k, 0);
  clear_amatrix(Yhat1);
  addmul_amatrix(zeta_age, false, &rwf->C, true, &rb->E, Yhat1);
  uninit_amatrix(Yhat1);

  /* columns associated with A */
  Yhat1 = init_sub_amatrix(&tmp2, Yhat, rwf->krow, 0, R->k, rb->k);
  clear_amatrix(Yhat1);
  uninit_amatrix(Yhat1);

  off = rwf->krow;
  u = rb->rlist;
  while (u != NULL) {
    /* Compute block weight if required */
    alpha = 1.0;
    if (tm && tm->blocks) {
      if (tm->frobenius)
	norm = normfrob_rkupdate_uniform(u, R->k);
      else
	norm = norm2_rkupdate_uniform(u, R->k);

      alpha = (norm > 0.0 ? 1.0 / norm : 1.0);
    }
    Z = u->cb->Z;
    /* u is a subblock of AB* */
    if (Z != NULL) {
      /* columns associated with V */
      Yhat1 = init_sub_amatrix(&tmp2, Yhat, Z->rows, off, rb->k, 0);
      Z1 = init_sub_amatrix(&tmp3, Z, Z->rows, 0, u->S.cols, 0);
      clear_amatrix(Yhat1);
      addmul_amatrix(alpha, false, Z1, true, &u->S, Yhat1);
      uninit_amatrix(Yhat1);
      uninit_amatrix(Z1);

      /* columns associated with A */
      Yhat1 = init_sub_amatrix(&tmp2, Yhat, Z->rows, off, R->k, rb->k);
      assert(R->k == Z->cols - u->S.cols);
      Z1 = init_sub_amatrix(&tmp3, Z, Z->rows, 0, R->k, u->S.cols);
      copy_amatrix(false, Z1, Yhat1);
      scale_amatrix(alpha, Yhat1);
      uninit_amatrix(Yhat1);
      uninit_amatrix(Z1);

      off += Z->rows;
    }
    /* u is no subblock of AB* */
    else {
      /* columns associated with V */
      Yhat1 = init_sub_amatrix(&tmp2, Yhat, u->S.cols, off, rb->k, 0);
      copy_amatrix(true, &u->S, Yhat1);
      scale_amatrix(alpha, Yhat1);
      uninit_amatrix(Yhat1);

      /* columns associated with A */
      Yhat1 = init_sub_amatrix(&tmp2, Yhat, u->S.cols, off, R->k, rb->k);
      clear_amatrix(Yhat1);
      uninit_amatrix(Yhat1);

      off += u->S.cols;
    }
    u = u->rnext;
  }
  assert(off == rows);

  /* compute total weight */
  k = UINT_MIN(rows, cols);
  tau = init_avector(&tmp4, k);
  qrdecomp_amatrix(Yhat, tau);
  resize_clusteroperator(rw, k, cols);
  copy_upper_amatrix(Yhat, false, &rw->C);

  uninit_avector(tau);
  uninit_amatrix(Yhat);

  /* compute total weight for son cluster of (V A) */
  rowweight_rkupdate_clusteroperator(rb, rw, R->k, tm);

  /* -------------------------------------------------- */
  /* compute total weight for current column cluster (W B) */

  /* rows of Yhat */
  rows = cwf->krow;
  u = cb->clist;
  while (u != NULL) {
    Z = u->rb->Z;
    /* u is a subblock of AB* */
    if (Z != NULL) {
      rows += Z->rows;
    }
    /* u is no subblock of AB* */
    else {
      rows += u->S.rows;
    }
    u = u->cnext;
  }
  /* cols of Yhat */
  cols = cb->k + R->k;
  Yhat = init_amatrix(&tmp1, rows, cols);

  /* columns associated with V */
  Yhat1 = init_sub_amatrix(&tmp2, Yhat, cwf->krow, 0, cb->k, 0);
  clear_amatrix(Yhat1);
  addmul_amatrix(zeta_age, false, &cwf->C, true, &cb->E, Yhat1);
  uninit_amatrix(Yhat1);

  /* columns associated with A */
  Yhat1 = init_sub_amatrix(&tmp2, Yhat, cwf->krow, 0, R->k, cb->k);
  clear_amatrix(Yhat1);
  uninit_amatrix(Yhat1);

  off = cwf->krow;
  u = cb->clist;
  while (u) {
    /* Compute block weight if required */
    alpha = 1.0;
    if (tm && tm->blocks) {
      if (tm->frobenius)
	norm = normfrob_rkupdate_uniform(u, R->k);
      else
	norm = norm2_rkupdate_uniform(u, R->k);

      alpha = (norm > 0.0 ? 1.0 / norm : 1.0);
    }
    Z = u->rb->Z;
    /* u is a subblock of AB* */
    if (Z != NULL) {
      /* columns associated with V */
      Yhat1 = init_sub_amatrix(&tmp2, Yhat, Z->rows, off, cb->k, 0);
      Z1 = init_sub_amatrix(&tmp3, Z, Z->rows, 0, u->S.rows, 0);
      clear_amatrix(Yhat1);
      addmul_amatrix(alpha, false, Z1, false, &u->S, Yhat1);
      uninit_amatrix(Yhat1);
      uninit_amatrix(Z1);

      /* columns associated with A */
      Yhat1 = init_sub_amatrix(&tmp2, Yhat, Z->rows, off, R->k, cb->k);
      assert(R->k == Z->cols - u->S.rows);
      Z1 = init_sub_amatrix(&tmp3, Z, Z->rows, 0, R->k, u->S.rows);
      copy_amatrix(false, Z1, Yhat1);
      scale_amatrix(alpha, Yhat1);
      uninit_amatrix(Yhat1);
      uninit_amatrix(Z1);

      off += Z->rows;
    }
    /* u is no subblock of AB* */
    else {
      /* columns associated with V */
      Yhat1 = init_sub_amatrix(&tmp2, Yhat, u->S.rows, off, cb->k, 0);
      copy_amatrix(false, &u->S, Yhat1);
      scale_amatrix(alpha, Yhat1);
      uninit_amatrix(Yhat1);

      /* columns associated with A */
      Yhat1 = init_sub_amatrix(&tmp2, Yhat, u->S.rows, off, R->k, cb->k);
      clear_amatrix(Yhat1);
      uninit_amatrix(Yhat1);

      off += u->S.rows;
    }
    u = u->cnext;
  }
  assert(off == rows);

  /* compute total weight */
  k = UINT_MIN(rows, cols);
  tau = init_avector(&tmp4, k);
  qrdecomp_amatrix(Yhat, tau);
  resize_clusteroperator(cw, k, cols);
  copy_upper_amatrix(Yhat, false, &cw->C);

  uninit_avector(tau);
  uninit_amatrix(Yhat);

  /* compute total weight for son cluster of (W B) */
  colweight_rkupdate_clusteroperator(cb, cw, R->k, tm);

  /* -------------------------------------------------- */
  /* truncate the row clusterbasis (V A) */
  k = rb->k;
  E1 = init_amatrix(&tmp2, k, rb->E.cols);
  copy_amatrix(false, &rb->E, E1);
  truncate_rkupdate_clusterbasis(rb, &R->A, rw, tm, eps);
  Z1 = init_sub_amatrix(&tmp3, &rw->C, rb->k, 0, k, 0);
  clear_amatrix(&rb->E);
  addmul_amatrix(1.0, false, Z1, false, E1, &rb->E);
  uninit_amatrix(Z1);
  uninit_amatrix(E1);

  /* -------------------------------------------------- */
  /* truncate the col clusterbasis (W B) */
  k = cb->k;
  E1 = init_amatrix(&tmp2, k, cb->E.cols);
  copy_amatrix(false, &cb->E, E1);
  truncate_rkupdate_clusterbasis(cb, &R->B, cw, tm, eps);
  Z1 = init_sub_amatrix(&tmp3, &cw->C, cb->k, 0, k, 0);
  clear_amatrix(&cb->E);
  addmul_amatrix(1.0, false, Z1, false, E1, &cb->E);
  uninit_amatrix(Z1);
  uninit_amatrix(E1);

  /* -------------------------------------------------- */
  /* update the subblocks of Gh2 */
  rkupdate_inside_h2matrix(Gh2, &R->A, &R->B, rw, cw);

  /* update the blocks outside of Gh2 */
  rkupdate_rowout_h2matrix(rb, rw);
  rkupdate_colout_h2matrix(cb, cw);

  /* -------------------------------------------------- */
  /* update the totalweights of current row clusterbasis */

  /* rows of Yhat */
  rows = rwf->krow;
  u = rb->rlist;
  while (u != NULL) {
    rows += u->S.cols;
    u = u->rnext;
  }
  /* cols of Yhat */
  cols = rb->k;
  Yhat = init_amatrix(&tmp1, rows, cols);

  /* transfer the weight of the father */
  Yhat1 = init_sub_amatrix(&tmp2, Yhat, rwf->krow, 0, rb->k, 0);
  clear_amatrix(Yhat1);
  addmul_amatrix(zeta_age, false, &rwf->C, true, &rb->E, Yhat1);
  uninit_amatrix(Yhat1);

  off = rwf->krow;
  u = rb->rlist;
  while (u) {
    /* Compute block weight if required */
    alpha = 1.0;
    if (tm && tm->blocks) {
      if (tm->frobenius)
	norm = normfrob_rkupdate_uniform(u, 0);
      else
	norm = norm2_rkupdate_uniform(u, 0);

      alpha = (norm > 0.0 ? 1.0 / norm : 1.0);
    }
    /* coupling matrices for admissible blocks */
    Yhat1 = init_sub_amatrix(&tmp2, Yhat, u->S.cols, off, rb->k, 0);
    copy_amatrix(true, &u->S, Yhat1);
    scale_amatrix(alpha, Yhat1);
    uninit_amatrix(Yhat1);

    off += u->S.cols;
    u = u->rnext;
  }
  assert(off == rows);

  /* compute the weight */
  refl = UINT_MIN(rows, cols);
  tau = init_avector(&tmp4, refl);
  qrdecomp_amatrix(Yhat, tau);
  resize_clusteroperator(rw, refl, cols);
  copy_upper_amatrix(Yhat, false, &rw->C);
  uninit_avector(tau);
  uninit_amatrix(Yhat);

  /* update the totalweights of sons of current row clusterbasis */
  totalweights_row_clusteroperator(rb, rw, tm);

  /* -------------------------------------------------- */
  /* update the totalweights of current col clusterbasis */

  /* rows of Yhat */
  rows = cwf->krow;
  u = cb->clist;
  while (u != NULL) {
    rows += u->S.rows;
    u = u->cnext;
  }
  /* cols of Yhat */
  cols = cb->k;
  Yhat = init_amatrix(&tmp1, rows, cols);

  /* transfer the weight of the father */
  Yhat1 = init_sub_amatrix(&tmp2, Yhat, cwf->krow, 0, cb->k, 0);
  clear_amatrix(Yhat1);
  addmul_amatrix(zeta_age, false, &cwf->C, true, &cb->E, Yhat1);
  uninit_amatrix(Yhat1);

  off = cwf->krow;
  u = cb->clist;
  while (u) {
    /* Compute block weight if required */
    alpha = 1.0;
    if (tm && tm->blocks) {
      if (tm->frobenius)
	norm = normfrob_rkupdate_uniform(u, 0);
      else
	norm = norm2_rkupdate_uniform(u, 0);

      alpha = (norm > 0.0 ? 1.0 / norm : 1.0);
    }
    /* coupling matrices for admissible blocks */
    Yhat1 = init_sub_amatrix(&tmp2, Yhat, u->S.rows, off, cb->k, 0);
    copy_amatrix(false, &u->S, Yhat1);
    scale_amatrix(alpha, Yhat1);
    uninit_amatrix(Yhat1);

    off += u->S.rows;
    u = u->cnext;
  }
  assert(off == rows);

  /* compute the weight */
  refl = UINT_MIN(rows, cols);
  tau = init_avector(&tmp4, refl);
  qrdecomp_amatrix(Yhat, tau);
  resize_clusteroperator(cw, refl, cols);
  copy_upper_amatrix(Yhat, false, &cw->C);
  uninit_avector(tau);
  uninit_amatrix(Yhat);

  /* update the totalweights of sons of current row clusterbasis */
  totalweights_col_clusteroperator(cb, cw, tm);

  /* clean up weights in clusterbasis */
  clear_weight_clusterbasis(rb);
  clear_weight_clusterbasis(cb);
}

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* builds a clusteroperator and computes the total weights for the row cluster basis */
pclusteroperator
prepare_row_clusteroperator(pclusterbasis rb, pclusterbasis cb, ptruncmode tm)
{
  pclusteroperator rw, rwf;

  rwf = (pclusteroperator) allocmem(sizeof(clusteroperator));
  init_leaf_clusteroperator(rwf, NULL);
  rwf->sons = 1;
  rwf->son = (pclusteroperator *) allocmem(sizeof(pclusteroperator));
  *rwf->son = NULL;

  rw = build_from_clusterbasis_clusteroperator(rb);
  orthoweight_clusterbasis(cb);
  totalweight_row_clusteroperator(rb, rw, tm);
  clear_weight_clusterbasis(cb);

  ref_clusteroperator(rwf->son, rw);

  return rwf;
}

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* builds a clusteroperator and computes the total weights for the column cluster basis */
pclusteroperator
prepare_col_clusteroperator(pclusterbasis rb, pclusterbasis cb, ptruncmode tm)
{
  pclusteroperator cw, cwf;

  cwf = (pclusteroperator) allocmem(sizeof(clusteroperator));
  init_leaf_clusteroperator(cwf, NULL);
  cwf->sons = 1;
  cwf->son = (pclusteroperator *) allocmem(sizeof(pclusteroperator));
  *cwf->son = NULL;

  cw = build_from_clusterbasis_clusteroperator(cb);
  orthoweight_clusterbasis(rb);
  totalweight_col_clusteroperator(cb, cw, tm);
  clear_weight_clusterbasis(rb);

  ref_clusteroperator(cwf->son, cw);

  return cwf;
}
