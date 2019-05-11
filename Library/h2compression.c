
/* ------------------------------------------------------------
 * This is the file "h2compression.c" of the H2Lib package.
 * All rights reserved, Steffen Boerm 2011
 * ------------------------------------------------------------ */

#include "h2compression.h"

#include "h2update.h"
#include "eigensolvers.h"
#include "factorizations.h"
#include "basic.h"

/* ------------------------------------------------------------
 * High-level compression functions
 * ------------------------------------------------------------ */

ph2matrix
compress_amatrix_h2matrix(pcamatrix G, pcblock b, pctruncmode tm, real eps)
{
  pclusterbasis rb, cb;
  ph2matrix G2;

  rb = buildrowbasis_amatrix(G, b, tm, eps);
  cb = buildcolbasis_amatrix(G, b, tm, eps);

  G2 = build_projected_amatrix_h2matrix(G, b, rb, cb);

  return G2;
}

ph2matrix
compress_hmatrix_h2matrix(pchmatrix G, pctruncmode tm, real eps)
{
  pclusterbasis rb, cb;
  ph2matrix G2;

  rb = buildrowbasis_hmatrix(G, tm, eps);
  cb = buildcolbasis_hmatrix(G, tm, eps);

  G2 = build_projected_hmatrix_h2matrix(G, rb, cb);

  return G2;
}

ph2matrix
compress_h2matrix_h2matrix(pch2matrix G, bool rbortho, bool cbortho,
			   pctruncmode tm, real eps)
{
  pclusteroperator rbw, cbw;
  pclusteroperator rlw, clw;
  pclusterbasis rb, cb;
  pclusteroperator ro, co;
  ph2matrix Gp;

  rbw = 0;
  if (!rbortho) {
    rbw = build_from_clusterbasis_clusteroperator(G->rb);
    weight_clusterbasis_clusteroperator(G->rb, rbw);
  }

  cbw = 0;
  if (!cbortho) {
    cbw = build_from_clusterbasis_clusteroperator(G->cb);
    weight_clusterbasis_clusteroperator(G->cb, cbw);
  }

  rlw = build_from_clusterbasis_clusteroperator(G->rb);
  clw = build_from_clusterbasis_clusteroperator(G->cb);
  localweights_h2matrix(G, rbw, cbw, tm, rlw, clw);

  rb = clonestructure_clusterbasis(G->rb);
  ro = build_from_clusterbasis_clusteroperator(G->rb);
  truncate_clusterbasis(G->rb, 0, clw, tm, eps, rb, ro);

  cb = clonestructure_clusterbasis(G->cb);
  co = build_from_clusterbasis_clusteroperator(G->cb);
  truncate_clusterbasis(G->cb, 0, rlw, tm, eps, cb, co);

  Gp = build_projected_h2matrix(G, rb, ro, cb, co);

  del_clusteroperator(co);
  del_clusteroperator(ro);
  del_clusteroperator(clw);
  del_clusteroperator(rlw);
  if (cbw)
    del_clusteroperator(cbw);
  if (rbw)
    del_clusteroperator(rbw);

  return Gp;
}

ph2matrix
compress_symmetric_h2matrix_h2matrix(pch2matrix G, bool rbortho,
				     pctruncmode tm, real eps)
{
  pclusteroperator rbw;
  pclusteroperator rlw;
  pclusterbasis rb;
  pclusteroperator ro;
  ph2matrix Gp;

  rbw = 0;
  if (!rbortho) {
    rbw = build_from_clusterbasis_clusteroperator(G->rb);
    weight_clusterbasis_clusteroperator(G->rb, rbw);
  }

  rlw = build_from_clusterbasis_clusteroperator(G->rb);

  rowweights_h2matrix(G, rbw, rbw, tm, rlw);

  rb = clonestructure_clusterbasis(G->rb);
  ro = build_from_clusterbasis_clusteroperator(G->rb);
  truncate_clusterbasis(G->rb, 0, rlw, tm, eps, rb, ro);

  Gp = build_projected_h2matrix(G, rb, ro, rb, ro);

  del_clusteroperator(ro);
  del_clusteroperator(rlw);
  if (rbw)
    del_clusteroperator(rbw);

  return Gp;
}

/* ------------------------------------------------------------
 * Compute local and total weights for H^2-matrices
 * ------------------------------------------------------------ */

typedef struct _weightsdata weightsdata;
typedef weightsdata *pweightsdata;
struct _weightsdata {
  pclusterbasis *rbn;
  pclusteroperator *rbwn;
  pclusteroperator *cbwn;
  pclusteroperator *rwn;
  pclusteroperator *cwn;
  pctruncmode tm;
  bool      Gtrans;
};

static    real
norm2_fast_uniform(pcuniform u, pcclusteroperator rw, pcclusteroperator cw)
{
  amatrix   tmp1, tmp2;
  pamatrix  ur, urc;
  real      norm;

  if (rw) {
    if (cw) {
      ur = init_amatrix(&tmp1, rw->krow, u->cb->k);
      urc = init_amatrix(&tmp2, rw->krow, cw->krow);

      clear_amatrix(ur);
      addmul_amatrix(1.0, false, &rw->C, false, &u->S, ur);
      clear_amatrix(urc);
      addmul_amatrix(1.0, false, ur, true, &cw->C, urc);

      uninit_amatrix(ur);
    }
    else {
      urc = init_amatrix(&tmp2, rw->krow, u->cb->k);

      clear_amatrix(urc);
      addmul_amatrix(1.0, false, &rw->C, false, &u->S, urc);
    }
  }
  else {
    if (cw) {
      urc = init_amatrix(&tmp2, u->rb->k, cw->krow);

      clear_amatrix(urc);
      addmul_amatrix(1.0, false, &u->S, true, &cw->C, urc);
    }
    else
      urc = (pamatrix) &u->S;
  }

  norm = norm2_amatrix(urc);

  if (urc != &u->S)
    uninit_amatrix(urc);

  return norm;
}

static    real
normfrob_fast_uniform(pcuniform u, pcclusteroperator rw, pcclusteroperator cw)
{
  amatrix   tmp1, tmp2;
  pamatrix  ur, urc;
  real      norm;

  if (rw) {
    if (cw) {
      ur = init_amatrix(&tmp1, rw->krow, u->cb->k);
      urc = init_amatrix(&tmp2, rw->krow, cw->krow);

      clear_amatrix(ur);
      addmul_amatrix(1.0, false, &rw->C, false, &u->S, ur);
      clear_amatrix(urc);
      addmul_amatrix(1.0, false, ur, true, &cw->C, urc);

      uninit_amatrix(ur);
    }
    else {
      urc = init_amatrix(&tmp2, rw->krow, u->cb->k);

      clear_amatrix(urc);
      addmul_amatrix(1.0, false, &rw->C, false, &u->S, urc);
    }
  }
  else {
    if (cw) {
      urc = init_amatrix(&tmp2, u->rb->k, cw->krow);

      clear_amatrix(urc);
      addmul_amatrix(1.0, false, &u->S, true, &cw->C, urc);
    }
    else
      urc = (pamatrix) &u->S;
  }

  norm = normfrob_amatrix(urc);

  if (urc != &u->S)
    uninit_amatrix(urc);

  return norm;
}

static void
localweights1(pccluster t, uint tname, uint pardepth,
	      pch2matrixlist hl, void *data)
{
  amatrix   tmp1, tmp2;
  avector   tmp3;
  pweightsdata wd = (pweightsdata) data;
  pcclusterbasis rb = wd->rbn[tname];
  pclusteroperator *rbwn = wd->rbwn;
  pclusteroperator *cbwn = wd->cbwn;
  pcclusteroperator rbw, cbw;
  pclusteroperator rw = wd->rwn[tname];
  pctruncmode tm = wd->tm;
  bool      Gtrans = wd->Gtrans;
  pch2matrix G;
  pamatrix  X, X1;
  pavector  tau;
  pch2matrixlist hl0;
  real      norm, alpha;
  uint      n, k, kmax;
  uint      off;

  /* Stop gcc from complaining about unused parameters */
  (void) t;
  (void) pardepth;

  k = rb->k;

  /*
   * Set up matrix X containing contributions from all admissible blocks
   */

  n = 0;
  X = 0;
  if (Gtrans) {
    /* Determine number of rows */
    for (hl0 = hl; hl0; hl0 = hl0->next) {
      G = hl0->G;
      rbw = rbwn[hl0->rname];

      assert(G->cb == rb);

      if (G->u) {
	if (rbw) {
	  assert(rbw->kcol == G->u->rb->k);

	  n += rbw->krow;
	}
	else
	  n += G->rb->k;
      }
    }

    /* Allocate matrix */
    X = init_amatrix(&tmp1, n, k);

    off = 0;
    for (hl0 = hl; hl0; hl0 = hl0->next) {
      G = hl0->G;
      rbw = rbwn[hl0->rname];
      cbw = cbwn[hl0->cname];

      if (G->u) {
	/* Compute weight factor */
	alpha = 1.0;
	if (tm && tm->blocks) {
	  if (tm->frobenius)
	    norm = normfrob_fast_uniform(G->u, rbw, cbw);
	  else
	    norm = norm2_fast_uniform(G->u, rbw, cbw);

	  alpha = (norm > 0.0 ? 1.0 / norm : 1.0);
	}

	if (rbw) {
	  /* Submatrix alpha R_t S_{ts} if a row weight R_t is given */
	  X1 = init_sub_amatrix(&tmp2, X, rbw->krow, off, k, 0);

	  clear_amatrix(X1);
	  addmul_amatrix(alpha, false, &rbw->C, false, &G->u->S, X1);

	  uninit_amatrix(X1);

	  off += rbw->krow;
	}
	else {
	  /* Submatrix alpha S_{ts} if no weight is given */
	  X1 = init_sub_amatrix(&tmp2, X, G->rb->k, off, k, 0);

	  clear_amatrix(X1);
	  add_amatrix(alpha, false, &G->u->S, X1);

	  uninit_amatrix(X1);

	  off += G->rb->k;
	}
      }
    }
    assert(off == n);
  }
  else {
    /* Determine number of rows */
    for (hl0 = hl; hl0; hl0 = hl0->next) {
      G = hl0->G;
      cbw = cbwn[hl0->cname];

      assert(G->rb == rb);

      if (G->u) {
	if (cbw) {
	  assert(cbwn[hl0->cname]->kcol == G->u->cb->k);

	  n += cbwn[hl0->cname]->krow;
	}
	else
	  n += G->cb->k;
      }
    }

    /* Allocate matrix */
    X = init_amatrix(&tmp1, n, k);

    off = 0;
    for (hl0 = hl; hl0; hl0 = hl0->next) {
      G = hl0->G;
      rbw = rbwn[hl0->rname];
      cbw = cbwn[hl0->cname];

      if (G->u) {
	/* Compute weight factor */
	alpha = 1.0;
	if (tm && tm->blocks) {
	  if (tm->frobenius)
	    norm = normfrob_fast_uniform(G->u, rbw, cbw);
	  else
	    norm = norm2_fast_uniform(G->u, rbw, cbw);

	  alpha = (norm > 0.0 ? 1.0 / norm : 1.0);
	}

	if (cbw) {
	  /* Submatrix alpha R_s S_{ts}^* if a column weight R_s is given */
	  X1 = init_sub_amatrix(&tmp2, X, cbw->krow, off, k, 0);

	  clear_amatrix(X1);
	  addmul_amatrix(alpha, false, &cbw->C, true, &G->u->S, X1);

	  uninit_amatrix(X1);

	  off += cbw->krow;
	}
	else {
	  /* Submatrix alpha S_{ts}^* if no weight is given */
	  X1 = init_sub_amatrix(&tmp2, X, G->cb->k, off, k, 0);

	  clear_amatrix(X1);
	  add_amatrix(alpha, true, &G->u->S, X1);

	  uninit_amatrix(X1);

	  off += G->cb->k;
	}
      }
    }
    assert(off == n);
  }

  /* Compute QR factorization */
  kmax = UINT_MIN(k, n);
  tau = init_avector(&tmp3, kmax);
  qrdecomp_amatrix(X, tau);

  /* Copy R into weight matrix */
  resize_clusteroperator(rw, kmax, k);
  copy_upper_amatrix(X, false, &rw->C);

  /* Clean up */
  uninit_avector(tau);
  uninit_amatrix(X);
}

void
rowweights_h2matrix(pch2matrix G, pcclusteroperator rbw,
		    pcclusteroperator cbw, pctruncmode tm,
		    pclusteroperator rlw)
{
  weightsdata wd;

  wd.rbn = enumerate_clusterbasis(G->rb->t, (pclusterbasis) G->rb);
  wd.rbwn = enumerate_clusteroperator(G->rb->t, (pclusteroperator) rbw);
  wd.cbwn = enumerate_clusteroperator(G->cb->t, (pclusteroperator) cbw);
  wd.rwn = enumerate_clusteroperator(G->rb->t, rlw);
  wd.cwn = 0;
  wd.tm = tm;
  wd.Gtrans = false;

  iterate_rowlist_h2matrix((ph2matrix) G, 0, 0, 0, max_pardepth,
			   localweights1, 0, &wd);

  /* Clean up */
  freemem(wd.rwn);
  freemem(wd.cbwn);
  freemem(wd.rbwn);
  freemem(wd.rbn);
}

void
colweights_h2matrix(pch2matrix G, pcclusteroperator rbw,
		    pcclusteroperator cbw, pctruncmode tm,
		    pclusteroperator clw)
{
  weightsdata wd;

  wd.rbn = enumerate_clusterbasis(G->cb->t, (pclusterbasis) G->cb);
  wd.rbwn = enumerate_clusteroperator(G->rb->t, (pclusteroperator) rbw);
  wd.cbwn = enumerate_clusteroperator(G->cb->t, (pclusteroperator) cbw);
  wd.rwn = enumerate_clusteroperator(G->cb->t, clw);
  wd.cwn = 0;
  wd.tm = tm;
  wd.Gtrans = true;

  iterate_collist_h2matrix((ph2matrix) G, 0, 0, 0, max_pardepth,
			   localweights1, 0, &wd);

  /* Clean up */
  freemem(wd.rwn);
  freemem(wd.cbwn);
  freemem(wd.rbwn);
  freemem(wd.rbn);
}

static void
localweights2(ph2matrix G, uint bname, uint rname, uint cname,
	      uint pardepth, void *data)
{
  pweightsdata wd = (pweightsdata) data;
  pcclusteroperator rbw = wd->rbwn[rname];
  pcclusteroperator cbw = wd->cbwn[cname];
  pclusteroperator rlw = wd->rwn[rname];
  pclusteroperator clw = wd->cwn[cname];
  pctruncmode tm = wd->tm;
  pamatrix  W, W1;
  pavector  tau;
  amatrix   tmp1, tmp2;
  avector   tmp3;
  pcclusterbasis rb, cb;
  real      norm, alpha;
  uint      n, k;

  /* Stop gcc from complaining about unused parameters */
  (void) bname;
  (void) pardepth;

  /* Quick exit if inadmissible submatrix */
  if (!G->u)
    return;

  rb = G->rb;
  cb = G->cb;

  assert(rb == G->u->rb);
  assert(cb == G->u->cb);
  assert(rlw->t == rb->t);
  assert(clw->t == cb->t);

  /*
   * Compute block weight if required
   */

  alpha = 1.0;
  if (tm && tm->blocks) {
    if (tm->frobenius)
      norm = normfrob_fast_uniform(G->u, rbw, cbw);
    else
      norm = norm2_fast_uniform(G->u, rbw, cbw);

    if (norm > 0.0)
      alpha = 1.0 / norm;
  }

  /*
   * Update the row weight
   */

  /* W combines the old row weight and the current block */
  n = (cbw ? rlw->krow + cbw->krow : rlw->krow + cb->k);
  W = init_amatrix(&tmp1, n, rb->k);

  /* Copy old row weight to upper block of W */
  if (rlw->krow > 0) {
    assert(rlw->kcol == rb->k);

    W1 = init_sub_amatrix(&tmp2, W, rlw->krow, 0, rb->k, 0);
    copy_amatrix(false, &rlw->C, W1);
    uninit_amatrix(W1);
  }

  /* Add W_s S_b^* or S_b^* to lower block of W */
  if (cbw) {
    W1 = init_sub_amatrix(&tmp2, W, cbw->krow, rlw->krow, rb->k, 0);
    clear_amatrix(W1);
    addmul_amatrix(alpha, false, &cbw->C, true, &G->u->S, W1);
    uninit_amatrix(W1);
  }
  else {
    W1 = init_sub_amatrix(&tmp2, W, cb->k, rlw->krow, rb->k, 0);
    clear_amatrix(W1);
    add_amatrix(alpha, true, &G->u->S, W1);
    uninit_amatrix(W1);
  }

  /* Compute QR factorization */
  k = UINT_MIN(rb->k, n);
  tau = init_avector(&tmp3, k);
  qrdecomp_amatrix(W, tau);
  uninit_avector(tau);

  /* Use R as new row weight */
  resize_clusteroperator(rlw, k, rb->k);
  copy_upper_amatrix(W, false, &rlw->C);
  uninit_amatrix(W);

  /*
   * Update the column weight
   */

  /* W combines the old row weight and the current block */
  n = (rbw ? clw->krow + rbw->krow : clw->krow + rb->k);
  W = init_amatrix(&tmp1, n, cb->k);

  /* Copy old column weight to upper block of W */
  if (clw->krow > 0) {
    assert(clw->kcol == cb->k);

    W1 = init_sub_amatrix(&tmp2, W, clw->krow, 0, cb->k, 0);
    copy_amatrix(false, &clw->C, W1);
    uninit_amatrix(W1);
  }

  /* Add V_t S_b or S_b to lower block of W */
  if (rbw) {
    W1 = init_sub_amatrix(&tmp2, W, rbw->krow, clw->krow, cb->k, 0);
    clear_amatrix(W1);
    addmul_amatrix(alpha, false, &rbw->C, false, &G->u->S, W1);
    uninit_amatrix(W1);
  }
  else {
    W1 = init_sub_amatrix(&tmp2, W, rb->k, clw->krow, cb->k, 0);
    clear_amatrix(W1);
    add_amatrix(alpha, false, &G->u->S, W1);
    uninit_amatrix(W1);
  }

  /* Compute QR factorization */
  k = UINT_MIN(cb->k, n);
  tau = init_avector(&tmp3, k);
  qrdecomp_amatrix(W, tau);
  uninit_avector(tau);

  /* Use R as new row weight */
  resize_clusteroperator(clw, k, cb->k);
  copy_upper_amatrix(W, false, &clw->C);
  uninit_amatrix(W);
}

void
localweights_h2matrix(pch2matrix G, pcclusteroperator rbw,
		      pcclusteroperator cbw, pctruncmode tm,
		      pclusteroperator rw, pclusteroperator cw)
{
  weightsdata wd;

  wd.rbn = 0;
  wd.rbwn = enumerate_clusteroperator(G->rb->t, (pclusteroperator) rbw);
  wd.cbwn = enumerate_clusteroperator(G->cb->t, (pclusteroperator) cbw);
  wd.rwn = enumerate_clusteroperator(G->rb->t, rw);
  wd.cwn = enumerate_clusteroperator(G->cb->t, cw);
  wd.tm = tm;

  iterate_h2matrix((ph2matrix) G, 0, 0, 0, max_pardepth, localweights2, 0,
		   &wd);

  freemem(wd.cwn);
  freemem(wd.rwn);
  freemem(wd.cbwn);
  freemem(wd.rbwn);
}

static void
accumulate(pccluster t, uint tname, void *data)
{
  pweightsdata wd = (pweightsdata) data;
  pcclusterbasis cb = wd->rbn[tname];
  pclusteroperator lw = wd->rwn[tname];
  pctruncmode tm = wd->tm;
  amatrix   tmp1, tmp2;
  avector   tmp3;
  pamatrix  w, w1;
  pavector  tau;
  real      zeta_age;
  uint      i, k;

  (void) t;

  assert(lw->kcol == cb->k);

  zeta_age = (tm ? tm->zeta_age : 1.0);

  if (lw->sons == 0)
    return;
  else {
    assert(cb->sons == lw->sons);

    for (i = 0; i < lw->sons; i++) {
      assert(lw->son[i]->kcol == cb->son[i]->k);

      /* w combines the son's local weight and the father's contribution */
      w = init_amatrix(&tmp1, lw->son[i]->krow + lw->krow, cb->son[i]->k);

      /* Copy son's local weight to upper block of w */
      w1 = init_sub_amatrix(&tmp2, w, lw->son[i]->krow, 0, cb->son[i]->k, 0);
      copy_amatrix(false, &lw->son[i]->C, w1);
      uninit_amatrix(w1);

      /* Add father's contribution using transfer matrix to lower block */
      w1 =
	init_sub_amatrix(&tmp2, w, lw->krow, lw->son[i]->krow, cb->son[i]->k,
			 0);
      clear_amatrix(w1);
      addmul_amatrix(zeta_age, false, &lw->C, true, &cb->son[i]->E, w1);
      uninit_amatrix(w1);

      /* Compute QR factorization */
      k = UINT_MIN(w->rows, w->cols);
      tau = init_avector(&tmp3, k);
      qrdecomp_amatrix(w, tau);
      uninit_avector(tau);

      /* Use R as new total weight */
      resize_clusteroperator(lw->son[i], k, cb->son[i]->k);
      copy_upper_amatrix(w, false, &lw->son[i]->C);

      uninit_amatrix(w);
    }
  }
}

void
accumulate_clusteroperator(pcclusterbasis cb, pctruncmode tm,
			   pclusteroperator lw)
{
  weightsdata wd;

  wd.rbn = enumerate_clusterbasis(cb->t, (pclusterbasis) cb);
  wd.rbwn = 0;
  wd.cbwn = 0;
  wd.rwn = enumerate_clusteroperator(cb->t, lw);
  wd.cwn = 0;
  wd.tm = tm;

  iterate_parallel_cluster(cb->t, 0, max_pardepth, accumulate, 0, &wd);

  freemem(wd.rwn);
  freemem(wd.rbn);
}

void
totalweights_h2matrix(pch2matrix G, bool rbortho, bool cbortho,
		      pctruncmode tm, pclusteroperator rw,
		      pclusteroperator cw)
{
  pclusteroperator rbw, cbw;

  /* Compute weights for row cluster basis */
  rbw = 0;
  if (!rbortho) {
    rbw = build_from_clusterbasis_clusteroperator(G->rb);
    weight_clusterbasis_clusteroperator(G->rb, rbw);
  }

  /* Compute weights for column cluster basis */
  cbw = 0;
  if (!cbortho) {
    cbw = rbw;
    if (G->rb != G->cb) {
      cbw = build_from_clusterbasis_clusteroperator(G->cb);
      weight_clusterbasis_clusteroperator(G->cb, cbw);
    }
  }

  /* Assemble local weights */
  localweights_h2matrix(G, rbw, cbw, tm, rw, cw);

  /* Assemble total weights */
  if (rw)
    accumulate_clusteroperator(G->rb, tm, rw);
  if (cw)
    accumulate_clusteroperator(G->cb, tm, cw);

  /* Clean up */
  if (cbw && cbw != rbw)
    del_clusteroperator(cbw);
  if (rbw)
    del_clusteroperator(rbw);
}

/* ------------------------------------------------------------
 * Compute truncated cluster basis
 * ------------------------------------------------------------ */

struct _truncate_data {
  pclusteroperator *cw;		/* Total cluster weights */
  pclusteroperator *clw;	/* Local cluster weights */
  pamatrix  Wn;			/* Accumulated local weights */

  pclusterbasis *cbold;		/* Old cluster basis */
  pclusterbasis *cbnew;		/* New cluster basis */
  pclusteroperator *old2new;	/* Basis change */

  pctruncmode tm;		/* Truncation mode */
  preal     eps;		/* Truncation accuracies */
};

static void
truncate_pre(pccluster t, uint tname, void *data)
{
  amatrix   tmp1, tmp2;
  avector   tmp3;
  struct _truncate_data *td = (struct _truncate_data *) data;
  pcclusteroperator clw = (td->clw ? td->clw[tname] : 0);
  pamatrix  Wn = td->Wn;
  pcclusterbasis cbold = td->cbold[tname];
  pctruncmode tm = td->tm;
  preal     eps = td->eps;
  pamatrix  W, W1;
  pavector  tau;
  uint      tname1;
  real      zeta_age, zeta_level;
  uint      kold, k, n, i;

  assert(cbold->t == t);

  /* Stop gcc from complaining about unused parameters */
  (void) t;

  /* Inherited weight is multiplied by this factor */
  zeta_age = (tm ? tm->zeta_age : 1.0);

  /* eps is multiplied by this factor if the level increases */
  zeta_level = (tm ? tm->zeta_level : 1.0);

  if (Wn) {
    kold = cbold->k;
    k = Wn[tname].rows;

    /* Inherited weight */
    assert(Wn[tname].cols == kold);

    /* Add local weight if available */
    if (clw) {
      assert(clw->t == t);
      assert(clw->kcol == kold);

      /* Combine inherited and local weight */
      n = k + clw->krow;

      W = init_amatrix(&tmp1, n, kold);

      W1 = init_sub_amatrix(&tmp2, W, k, 0, kold, 0);
      copy_amatrix(false, Wn + tname, W1);
      uninit_amatrix(W1);

      W1 = init_sub_amatrix(&tmp2, W, clw->krow, k, kold, 0);
      copy_amatrix(false, &clw->C, W1);
      uninit_amatrix(W1);

      /* Find QR decomposition */
      k = UINT_MIN(kold, n);
      tau = init_avector(&tmp3, k);
      qrdecomp_amatrix(W, tau);

      /* Copy R factor to Wn[tname] */
      resize_amatrix(Wn + tname, k, kold);
      copy_upper_amatrix(W, false, Wn + tname);

      /* Clean up */
      uninit_avector(tau);
      uninit_amatrix(W);
    }

    /* Propagate to sons if necessary */
    tname1 = tname + 1;
    for (i = 0; i < cbold->sons; i++) {
      init_amatrix(Wn + tname1, k, cbold->son[i]->k);

      clear_amatrix(Wn + tname1);
      addmul_amatrix(zeta_age, false, Wn + tname, true, &cbold->son[i]->E,
		     Wn + tname1);

      tname1 += t->son[i]->desc;
    }
    assert(tname1 == tname + t->desc);
  }

  /* Compute accuracies for sons */
  tname1 = tname + 1;
  for (i = 0; i < cbold->sons; i++) {
    eps[tname1] = eps[tname] * zeta_level;

    tname1 += t->son[i]->desc;
  }
  assert(tname1 == tname + t->desc);
}

static void
truncate_post(pccluster t, uint tname, void *data)
{
  amatrix   tmp1, tmp2, tmp3;
  realavector tmp4;
  avector   tmp5;
  struct _truncate_data *td = (struct _truncate_data *) data;
  pcclusterbasis cbold = td->cbold[tname];
  pclusterbasis cbnew = td->cbnew[tname];
  pclusteroperator old2new = td->old2new[tname];
  pcclusteroperator cw = (td->cw ? td->cw[tname] : 0);
  pamatrix  Wn = td->Wn;
  pctruncmode tm = td->tm;
  preal     eps = td->eps;
  pamatrix  Vhat, Vhat1, VhatZ, X, X1, Z, Q, Q1;
  pavector  tau;
  prealavector sigma;
  uint      m, kold, k, kmax;
  uint      tname1;
  uint      i, off;

  assert(cbold->t == t);
  assert(cbnew->t == t);
  assert(old2new->t == t);

  /* Stop gcc from complaining about unused parameters */
  (void) t;

  kold = cbold->k;

  /*
   * Construct Vhat
   */

  m = 0;
  Vhat = 0;
  if (cbold->sons > 0) {
    assert(cbnew->sons == cbold->sons);

    /* Determine number of rows */
    for (i = 0; i < cbnew->sons; i++)
      m += cbnew->son[i]->k;

    /* Allocate Vhat */
    Vhat = init_amatrix(&tmp1, m, kold);

    /* Fill submatrices */
    off = 0;
    for (i = 0; i < cbnew->sons; i++) {
      Vhat1 = init_sub_amatrix(&tmp2, Vhat, cbnew->son[i]->k, off, kold, 0);

      clear_amatrix(Vhat1);
      addmul_amatrix(1.0, false, &old2new->son[i]->C, false,
		     &cbold->son[i]->E, Vhat1);

      uninit_amatrix(Vhat1);

      off += cbnew->son[i]->k;
    }
    assert(off == m);
  }
  else {
    /* For leaves, the number of rows equals the size of t */
    m = cbold->t->size;

    /* Allocate Vhat */
    Vhat = init_amatrix(&tmp1, m, kold);

    /* Copy directly from original basis */
    copy_amatrix(false, &cbold->V, Vhat);
  }

  /*
   * Construct weighted matrix VhatZ
   */

  VhatZ = 0;
  if (cw) {
    assert(cw->kcol == kold);

    if (Wn) {
      /* Total and local weights present, merge them into Z */
      assert(Wn != 0);
      assert(Wn[tname].cols == kold);

      /* Copy cw->C and Wn[tname] into a block matrix */
      X = init_amatrix(&tmp2, cw->krow + Wn[tname].rows, kold);

      X1 = init_sub_amatrix(&tmp3, X, cw->krow, 0, kold, 0);
      copy_amatrix(false, &cw->C, X1);
      uninit_amatrix(X1);

      X1 = init_sub_amatrix(&tmp3, X, Wn[tname].rows, cw->krow, kold, 0);
      copy_amatrix(false, Wn + tname, X1);
      uninit_amatrix(X1);

      /* Compute QR factorization */
      k = UINT_MIN(X->rows, X->cols);
      tau = init_avector(&tmp5, k);
      qrdecomp_amatrix(X, tau);

      /* R factor is the new weight */
      Z = init_amatrix(&tmp3, k, kold);
      copy_upper_amatrix(X, false, Z);

      /* Clean up intermediate matrices */
      uninit_avector(tau);
      uninit_amatrix(X);

      /* Multiply Vhat by Z */
      VhatZ = init_amatrix(&tmp2, m, k);
      clear_amatrix(VhatZ);
      addmul_amatrix(1.0, false, Vhat, true, Z, VhatZ);

      /* Clean up */
      uninit_amatrix(Z);
    }
    else {
      /* Multiply Vhat by cw->C */
      VhatZ = init_amatrix(&tmp2, m, cw->krow);
      clear_amatrix(VhatZ);
      addmul_amatrix(1.0, false, Vhat, true, &cw->C, VhatZ);
    }
  }
  else {
    if (Wn) {
      assert(Wn[tname].cols == kold);

      /* Multiply Vhat by Wn[tname] */
      VhatZ = init_amatrix(&tmp2, m, Wn[tname].rows);
      clear_amatrix(VhatZ);
      addmul_amatrix(1.0, false, Vhat, true, Wn + tname, VhatZ);
    }
    else {
      /* No weighting, VhatZ equals Vhat */
      VhatZ = init_amatrix(&tmp2, m, kold);
      copy_amatrix(false, Vhat, VhatZ);
    }
  }
  assert(VhatZ->rows == m);

  /*
   * Find new basis by SVD
   */

  kmax = UINT_MIN(m, VhatZ->cols);
  Q = init_amatrix(&tmp3, m, kmax);
  sigma = init_realavector(&tmp4, kmax);
  svd_amatrix(VhatZ, sigma, Q, 0);

  /* Determine new rank */
  k = findrank_truncmode(tm, eps[tname], sigma);

  /* DEBUGGING
     printf("%u:", k);
     for(i=0; i<kmax; i++)
     printf(" %.4e", sigma->v[i]);
     printf("\n");
   */

  /* Clean up intermediate matrices */
  uninit_realavector(sigma);
  uninit_amatrix(VhatZ);

  /* Set rank of new cluster basis */
  resize_clusterbasis(cbnew, k);

  /* Compute basis change matrix */
  resize_clusteroperator(old2new, k, kold);
  Q1 = init_sub_amatrix(&tmp2, Q, m, 0, k, 0);
  clear_amatrix(&old2new->C);
  addmul_amatrix(1.0, true, Q1, false, Vhat, &old2new->C);
  uninit_amatrix(Q1);

  /* Copy result to cluster basis */
  if (cbold->sons > 0) {
    off = 0;
    for (i = 0; i < cbold->sons; i++) {
      Q1 = init_sub_amatrix(&tmp2, Q, cbnew->son[i]->k, off, k, 0);
      copy_amatrix(false, Q1, &cbnew->son[i]->E);
      uninit_amatrix(Q1);

      off += cbnew->son[i]->k;
    }
    assert(off == m);
  }
  else {
    Q1 = init_sub_amatrix(&tmp2, Q, t->size, 0, k, 0);
    copy_amatrix(false, Q1, &cbnew->V);
    uninit_amatrix(Q1);
  }

  /* Clean up */
  uninit_amatrix(Q);
  uninit_amatrix(Vhat);
  tname1 = tname + 1;
  if (Wn) {
    for (i = 0; i < cbold->sons; i++) {
      uninit_amatrix(Wn + tname1);
      tname1 += cbold->son[i]->t->desc;
    }
    assert(tname1 == tname + cbold->t->desc);
  }
}

void
truncate_clusterbasis(pcclusterbasis cb, pcclusteroperator cw,
		      pcclusteroperator clw, pctruncmode tm, real eps,
		      pclusterbasis cbnew, pclusteroperator old2new)
{
  struct _truncate_data td;

  td.cw = (cw ? enumerate_clusteroperator(cw->t, (pclusteroperator) cw) : 0);
  td.clw =
    (clw ? enumerate_clusteroperator(clw->t, (pclusteroperator) clw) : 0);
  td.Wn = (clw ? (pamatrix) allocmem(sizeof(amatrix) * cb->t->desc) : 0);
  td.cbold = enumerate_clusterbasis(cb->t, (pclusterbasis) cb);
  td.cbnew = enumerate_clusterbasis(cbnew->t, cbnew);
  td.old2new = enumerate_clusteroperator(cb->t, old2new);
  td.tm = tm;
  td.eps = allocreal(cb->t->desc);

  td.eps[0] = eps;
  if (td.Wn)
    init_amatrix(td.Wn, 0, cb->k);

  iterate_cluster(cb->t, 0, truncate_pre, truncate_post, &td);

  if (td.Wn)
    uninit_amatrix(td.Wn);
  freemem(td.eps);
  freemem(td.old2new);
  freemem(td.cbnew);
  freemem(td.cbold);
  if (clw) {
    freemem(td.Wn);
    freemem(td.clw);
  }
  if (cw)
    freemem(td.cw);
}

pclusterbasis
buildrowbasis_h2matrix(pch2matrix G, bool rbortho, bool cbortho,
		       pctruncmode tm, real eps, pclusteroperator old2new)
{
  pclusterbasis rb;
  pclusteroperator rbw, cbw;
  pclusteroperator rlw;

  /* Build weights for the row cluster basis */
  rbw = 0;
  if (!rbortho) {
    rbw = build_from_clusterbasis_clusteroperator(G->rb);
    weight_clusterbasis_clusteroperator(G->rb, rbw);
  }

  /* Build weights for the column cluster basis */
  cbw = 0;
  if (G->rb == G->cb)
    cbw = rbw;
  else if (!cbortho) {
    cbw = build_from_clusterbasis_clusteroperator(G->cb);
    weight_clusterbasis_clusteroperator(G->cb, cbw);
  }

  /* Build local weights */
  rlw = build_from_clusterbasis_clusteroperator(G->rb);
  localweights_h2matrix(G, rbw, cbw, tm, rlw, 0);

  /* Build new cluster basis by weighted truncation */
  rb = clonestructure_clusterbasis(G->rb);
  truncate_clusterbasis(G->rb, 0, rlw, tm, eps, rb, old2new);

  /* Clean up */
  del_clusteroperator(rlw);
  if (cbw && cbw != rbw)
    del_clusteroperator(cbw);
  if (rbw)
    del_clusteroperator(rbw);

  return rb;
}

pclusterbasis
buildcolbasis_h2matrix(pch2matrix G, bool rbortho, bool cbortho,
		       pctruncmode tm, real eps, pclusteroperator old2new)
{
  pclusterbasis cb;
  pclusteroperator rbw, cbw;
  pclusteroperator clw;

  /* Build weights for the row cluster basis */
  rbw = 0;
  if (!rbortho) {
    rbw = build_from_clusterbasis_clusteroperator(G->rb);
    weight_clusterbasis_clusteroperator(G->rb, rbw);
  }

  /* Build weights for the column cluster basis */
  cbw = 0;
  if (G->rb == G->cb)
    cbw = rbw;
  else if (!cbortho) {
    cbw = build_from_clusterbasis_clusteroperator(G->cb);
    weight_clusterbasis_clusteroperator(G->cb, cbw);
  }

  /* Build local weights */
  clw = build_from_clusterbasis_clusteroperator(G->cb);
  localweights_h2matrix(G, rbw, cbw, tm, 0, clw);

  /* Build new cluster basis by weighted truncation */
  cb = clonestructure_clusterbasis(G->cb);
  truncate_clusterbasis(G->cb, 0, clw, tm, eps, cb, old2new);

  /* Clean up */
  del_clusteroperator(clw);
  if (cbw && cbw != rbw)
    del_clusteroperator(cbw);
  if (rbw)
    del_clusteroperator(rbw);

  return cb;
}

struct _computebasis_data {
  pamatrix  rbw;		/* Weights for row basis */
  pamatrix  cbw;		/* Weights for column basis */

  pamatrix  Z;			/* Total basis weights */

  pclusterbasis *cbold;		/* Old cluster basis */
  pclusterbasis *cbnew;		/* New cluster basis */
  pclusteroperator *old2new;	/* Basis change */

  pctruncmode tm;		/* Truncation mode */
  preal     eps;		/* Truncation accuracies */
  bool      Gtrans;		/* Use S^* instead of S? */
};

/* Compute spectral norm of a product S_b R_s^* or S_b^* R_t^* */
static    real
norm2_product(bool Strans, pcamatrix S, pcamatrix R)
{
  amatrix   tmp;
  pamatrix  Z;
  real      norm;

  norm = 0.0;
  if (R) {
    if (Strans) {
      assert(S->rows == R->cols);

      Z = init_amatrix(&tmp, S->cols, R->rows);
      clear_amatrix(Z);
      addmul_amatrix(1.0, true, S, true, R, Z);
    }
    else {
      assert(S->cols == R->cols);

      Z = init_amatrix(&tmp, S->rows, R->rows);
      clear_amatrix(Z);
      addmul_amatrix(1.0, false, S, true, R, Z);
    }
    norm = norm2_amatrix(Z);
    uninit_amatrix(Z);
  }
  else
    norm = norm2_amatrix(S);

  return norm;
}

/* Compute Frobenius norm of a product S_b R_s^* or S_b^* R_t^* */
static    real
normfrob_product(bool Strans, pcamatrix S, pcamatrix R)
{
  amatrix   tmp;
  pamatrix  Z;
  real      norm;

  norm = 0.0;
  if (R) {
    if (Strans) {
      assert(S->rows == R->cols);

      Z = init_amatrix(&tmp, S->cols, R->rows);
      clear_amatrix(Z);
      addmul_amatrix(1.0, true, S, true, R, Z);
    }
    else {
      assert(S->cols == R->cols);

      Z = init_amatrix(&tmp, S->rows, R->rows);
      clear_amatrix(Z);
      addmul_amatrix(1.0, false, S, true, R, Z);
    }
    norm = normfrob_amatrix(Z);
    uninit_amatrix(Z);
  }
  else
    norm = normfrob_amatrix(S);

  return norm;
}

/* Top-down phase of basis construction: prepare cluster weights */
static void
computebasis_pre(pccluster t, uint tname, uint pardepth,
		 pch2matrixlist hl, void *data)
{
  amatrix   tmp1, tmp2;
  avector   tmp3;
  struct _computebasis_data *cd = (struct _computebasis_data *) data;
  pamatrix  rbw = cd->rbw;
  pamatrix  cbw = cd->cbw;
  pamatrix  Z = cd->Z;
  pclusterbasis cbold = cd->cbold[tname];
  bool      Gtrans = cd->Gtrans;
  pctruncmode tm = cd->tm;
  preal     eps = cd->eps;
  real      zeta_level, zeta_age;
  pamatrix  Zhat, Zhat1;
  pavector  tau;
  pch2matrixlist hl0;
  real      norm;
  uint      n, k;
  uint      i, off, tname1;

  /* Stop gcc from complaining about unused parameters */
  (void) pardepth;

  /* Quick exit */
  if (cbold == 0)
    return;

  /* eps is multiplied by this factor if the level increases */
  zeta_level = (tm ? tm->zeta_level : 1.0);

  /* Inherited weight is multiplied by this factor */
  zeta_age = (tm ? tm->zeta_age : 1.0);

  /* Inherited weight */
  assert(Z[tname].cols == cbold->k);

  if (Gtrans) {
    /* Determine size of Zhat */
    n = Z[tname].rows;
    if (rbw) {
      for (hl0 = hl; hl0; hl0 = hl0->next) {
	assert(hl0->G->cb->t == t);

	if (hl0->G->u) {
	  assert(rbw[hl0->rname].cols == hl0->G->rb->k);
	  n += rbw[hl0->rname].rows;
	}
      }
    }
    else {
      for (hl0 = hl; hl0; hl0 = hl0->next) {
	assert(hl0->G->cb->t == t);

	if (hl0->G->u)
	  n += hl0->G->rb->k;
      }
    }

    /* Create Zhat */
    Zhat = init_amatrix(&tmp1, n, cbold->k);

    /* Copy inherited part */
    Zhat1 = init_sub_amatrix(&tmp2, Zhat, Z[tname].rows, 0, Z[tname].cols, 0);
    copy_amatrix(false, Z + tname, Zhat1);
    uninit_amatrix(Zhat1);

    /* Add admissible submatrices */
    off = Z[tname].rows;
    if (rbw) {
      for (hl0 = hl; hl0; hl0 = hl0->next)
	if (hl0->G->u) {
	  Zhat1 = init_sub_amatrix(&tmp2, Zhat, rbw[hl0->rname].rows, off,
				   cbold->k, 0);
	  clear_amatrix(Zhat1);
	  addmul_amatrix(1.0, false, rbw + hl0->rname, false, &hl0->G->u->S,
			 Zhat1);

	  /* Scale block if required */
	  if (tm && tm->blocks) {
	    if (tm->frobenius)
	      norm = normfrob_product(false, Zhat1,
				      (cbw ? cbw + hl0->cname : 0));
	    else
	      norm =
		norm2_product(false, Zhat1, (cbw ? cbw + hl0->cname : 0));

	    if (norm > 0.0)
	      scale_amatrix(1.0 / norm, Zhat1);
	  }
	  uninit_amatrix(Zhat1);

	  off += rbw[hl0->rname].rows;
	}
    }
    else {
      for (hl0 = hl; hl0; hl0 = hl0->next)
	if (hl0->G->u) {
	  Zhat1 = init_sub_amatrix(&tmp2, Zhat, hl0->G->rb->k, off, cbold->k,
				   0);
	  copy_amatrix(false, &hl0->G->u->S, Zhat1);

	  /* Scale block if required */
	  if (tm && tm->blocks) {
	    if (tm->frobenius)
	      norm = normfrob_product(false, Zhat1,
				      (cbw ? cbw + hl0->cname : 0));
	    else
	      norm =
		norm2_product(false, Zhat1, (cbw ? cbw + hl0->cname : 0));

	    if (norm > 0.0)
	      scale_amatrix(1.0 / norm, Zhat1);
	  }
	  uninit_amatrix(Zhat1);

	  off += hl0->G->rb->k;
	}
    }
    assert(off == n);
  }
  else {
    /* Determine size of Zhat */
    n = Z[tname].rows;
    if (cbw) {
      for (hl0 = hl; hl0; hl0 = hl0->next) {
	assert(hl0->G->rb->t == t);

	if (hl0->G->u) {
	  assert(cbw[hl0->cname].cols == hl0->G->cb->k);
	  n += cbw[hl0->cname].rows;
	}
      }
    }
    else {
      for (hl0 = hl; hl0; hl0 = hl0->next) {
	assert(hl0->G->rb->t == t);

	if (hl0->G->u)
	  n += hl0->G->cb->k;
      }
    }

    /* Create Zhat */
    Zhat = init_amatrix(&tmp1, n, cbold->k);

    /* Copy inherited part */
    Zhat1 = init_sub_amatrix(&tmp2, Zhat, Z[tname].rows, 0, Z[tname].cols, 0);
    copy_amatrix(false, Z + tname, Zhat1);
    uninit_amatrix(Zhat1);

    /* Add admissible submatrices */
    off = Z[tname].rows;
    if (cbw) {
      for (hl0 = hl; hl0; hl0 = hl0->next)
	if (hl0->G->u) {
	  Zhat1 = init_sub_amatrix(&tmp2, Zhat, cbw[hl0->cname].rows, off,
				   cbold->k, 0);
	  clear_amatrix(Zhat1);
	  addmul_amatrix(1.0, false, cbw + hl0->cname, true, &hl0->G->u->S,
			 Zhat1);

	  /* Scale block if required */
	  if (tm && tm->blocks) {
	    if (tm->frobenius)
	      norm = normfrob_product(false, Zhat1,
				      (rbw ? rbw + hl0->rname : 0));
	    else
	      norm =
		norm2_product(false, Zhat1, (rbw ? rbw + hl0->rname : 0));

	    if (norm > 0.0)
	      scale_amatrix(1.0 / norm, Zhat1);
	  }
	  uninit_amatrix(Zhat1);

	  off += cbw[hl0->cname].rows;
	}
    }
    else {
      for (hl0 = hl; hl0; hl0 = hl0->next)
	if (hl0->G->u) {
	  Zhat1 = init_sub_amatrix(&tmp2, Zhat, hl0->G->cb->k, off, cbold->k,
				   0);
	  copy_amatrix(true, &hl0->G->u->S, Zhat1);

	  /* Scale block if required */
	  if (tm && tm->blocks) {
	    if (tm->frobenius)
	      norm = normfrob_product(false, Zhat1,
				      (rbw ? rbw + hl0->rname : 0));
	    else
	      norm =
		norm2_product(false, Zhat1, (rbw ? rbw + hl0->rname : 0));

	    if (norm > 0.0)
	      scale_amatrix(1.0 / norm, Zhat1);
	  }
	  uninit_amatrix(Zhat1);

	  off += hl0->G->cb->k;
	}
    }
    assert(off == n);
  }

  /* Compute QR decomposition of Zhat */
  k = UINT_MIN(n, cbold->k);
  tau = init_avector(&tmp3, k);
  qrdecomp_amatrix(Zhat, tau);
  uninit_avector(tau);

  /* Copy into Z */
  resize_amatrix(Z + tname, k, cbold->k);
  copy_upper_amatrix(Zhat, false, Z + tname);

  /* Initialize sons */
  if (cbold->sons > 0) {
    assert(cbold->sons == t->sons);

    tname1 = tname + 1;
    for (i = 0; i < t->sons; i++) {
      eps[tname1] = eps[tname] * zeta_level;

      init_amatrix(Z + tname1, k, cbold->son[i]->k);
      clear_amatrix(Z + tname1);
      addmul_amatrix(zeta_age, false, Z + tname, true, &cbold->son[i]->E,
		     Z + tname1);

      tname1 += t->son[i]->desc;
    }
    assert(tname1 == tname + t->desc);
  }

  /* Clean up */
  uninit_amatrix(Zhat);
}

/* Bottom-up phase of basis construction: find nested basis */
static void
computebasis_post(pccluster t, uint tname, uint pardepth,
		  pch2matrixlist hl, void *data)
{
  amatrix   tmp1, tmp2, tmp3, tmp4;
  realavector tmp5;
  struct _computebasis_data *cd = (struct _computebasis_data *) data;
  pamatrix  Z = cd->Z;
  pclusterbasis cbold = cd->cbold[tname];
  pclusterbasis cbnew = cd->cbnew[tname];
  pclusteroperator old2new = cd->old2new[tname];
  pctruncmode tm = cd->tm;
  preal     eps = cd->eps;
  pamatrix  Vhat, Vhat1, VhatZ, Q, Q1;
  prealavector sigma;
  uint      tname1;
  uint      m, k, kmax;
  uint      i, off;

  /* Stop gcc from complaining about unused parameters */
  (void) hl;
  (void) pardepth;

  /* Quick exit */
  if (cbold->k == 0 || Z[tname].rows == 0) {
    resize_clusterbasis(cbnew, 0);
    resize_clusteroperator(old2new, 0, 0);

    if (cbold->sons > 0) {
      assert(cbold->sons == t->sons);

      tname1 = tname + 1;
      for (i = 0; i < cbold->sons; i++) {
	uninit_amatrix(Z + tname1);

	tname1 += cbold->son[i]->t->desc;
      }
      assert(tname1 == tname + cbold->t->desc);
    }

    return;
  }

  /* Construct Vhat */
  if (cbold->sons > 0) {
    assert(cbold->sons == cbnew->sons);
    assert(cbold->sons == old2new->sons);

    m = 0;
    for (i = 0; i < cbold->sons; i++)
      m += cbnew->son[i]->k;

    Vhat = init_amatrix(&tmp1, m, cbold->k);

    off = 0;
    for (i = 0; i < cbold->sons; i++) {
      Vhat1 =
	init_sub_amatrix(&tmp2, Vhat, cbnew->son[i]->k, off, cbold->k, 0);

      clear_amatrix(Vhat1);
      addmul_amatrix(1.0, false, &old2new->son[i]->C, false,
		     &cbold->son[i]->E, Vhat1);

      uninit_amatrix(Vhat1);

      off += cbnew->son[i]->k;
    }
    assert(off == m);
  }
  else {
    m = t->size;

    Vhat = init_amatrix(&tmp1, m, cbold->k);

    copy_amatrix(false, &cbold->V, Vhat);
  }

  /* Multiply by weight matrix */
  assert(Z[tname].cols == cbold->k);
  VhatZ = init_amatrix(&tmp2, m, Z[tname].rows);

  clear_amatrix(VhatZ);
  addmul_amatrix(1.0, false, Vhat, true, Z + tname, VhatZ);

  /* Compute singular value decomposition */
  kmax = UINT_MIN(m, Z[tname].rows);
  Q = init_amatrix(&tmp3, m, kmax);
  sigma = init_realavector(&tmp5, kmax);
  svd_amatrix(VhatZ, sigma, Q, 0);

  /* Determine new rank */
  k = findrank_truncmode(tm, eps[tname], sigma);

  /* DEBUGGING
     printf("%u:", k);
     for(i=0; i<kmax; i++)
     printf(" %.4e", sigma->v[i]);
     printf("\n");
   */

  /* Set rank of new cluster basis */
  resize_clusterbasis(cbnew, k);

  /* Compute basis change matrix */
  resize_clusteroperator(old2new, k, cbold->k);
  Q1 = init_sub_amatrix(&tmp4, Q, m, 0, k, 0);
  clear_amatrix(&old2new->C);
  addmul_amatrix(1.0, true, Q1, false, Vhat, &old2new->C);
  uninit_amatrix(Q1);

  /* Copy result to cluster basis */
  if (cbold->sons > 0) {
    off = 0;
    for (i = 0; i < cbold->sons; i++) {
      Q1 = init_sub_amatrix(&tmp4, Q, cbnew->son[i]->k, off, k, 0);
      copy_amatrix(false, Q1, &cbnew->son[i]->E);
      uninit_amatrix(Q1);

      off += cbnew->son[i]->k;
    }
    assert(off == m);
  }
  else {
    Q1 = init_sub_amatrix(&tmp4, Q, t->size, 0, k, 0);
    copy_amatrix(false, Q1, &cbnew->V);
    uninit_amatrix(Q1);
  }

  /* Clean up */
  uninit_realavector(sigma);
  uninit_amatrix(Q);
  uninit_amatrix(VhatZ);
  uninit_amatrix(Vhat);
  if (cbold->sons > 0) {
    assert(cbold->sons == t->sons);

    tname1 = tname + 1;
    for (i = 0; i < cbold->sons; i++) {
      uninit_amatrix(Z + tname1);

      tname1 += cbold->son[i]->t->desc;
    }
    assert(tname1 == tname + cbold->t->desc);
  }
}

void
truncrowbasis_h2matrix(ph2matrix G, bool rbortho, bool cbortho,
		       pctruncmode tm, real eps, pclusterbasis rbnew,
		       pclusteroperator old2new)
{
  struct _computebasis_data cd;
  uint      i;

  cd.rbw = (rbortho ? 0 : weight_enum_clusterbasis_clusteroperator(G->rb));
  cd.cbw = (cbortho ? 0 : weight_enum_clusterbasis_clusteroperator(G->cb));
  cd.Z = (pamatrix) allocmem(sizeof(amatrix) * G->rb->t->desc);
  cd.cbold = enumerate_clusterbasis(G->rb->t, G->rb);
  cd.cbnew = enumerate_clusterbasis(G->rb->t, rbnew);
  cd.old2new = enumerate_clusteroperator(G->rb->t, old2new);
  cd.tm = tm;
  cd.eps = allocreal(G->rb->t->desc);
  cd.Gtrans = false;

  cd.eps[0] = eps;
  init_amatrix(cd.Z, 0, G->rb->k);

  iterate_rowlist_h2matrix(G, 0, 0, 0, max_pardepth, computebasis_pre,
			   computebasis_post, &cd);

  uninit_amatrix(cd.Z);

  freemem(cd.eps);
  freemem(cd.old2new);
  freemem(cd.cbnew);
  freemem(cd.cbold);
  freemem(cd.Z);
  if (cd.cbw) {
    for (i = 0; i < G->cb->t->desc; i++)
      uninit_amatrix(cd.cbw + i);
    freemem(cd.cbw);
  }
  if (cd.rbw) {
    for (i = 0; i < G->rb->t->desc; i++)
      uninit_amatrix(cd.rbw + i);
    freemem(cd.rbw);
  }
}

void
trunccolbasis_h2matrix(ph2matrix G, bool rbortho, bool cbortho,
		       pctruncmode tm, real eps, pclusterbasis cbnew,
		       pclusteroperator old2new)
{
  struct _computebasis_data cd;
  uint      i;

  cd.rbw = (rbortho ? 0 : weight_enum_clusterbasis_clusteroperator(G->rb));
  cd.cbw = (cbortho ? 0 : weight_enum_clusterbasis_clusteroperator(G->cb));
  cd.Z = (pamatrix) allocmem(sizeof(amatrix) * G->cb->t->desc);
  cd.cbold = enumerate_clusterbasis(G->cb->t, G->cb);
  cd.cbnew = enumerate_clusterbasis(G->cb->t, cbnew);
  cd.old2new = enumerate_clusteroperator(G->cb->t, old2new);
  cd.tm = tm;
  cd.eps = allocreal(G->cb->t->desc);
  cd.Gtrans = true;

  cd.eps[0] = eps;
  init_amatrix(cd.Z, 0, G->cb->k);

  iterate_collist_h2matrix(G, 0, 0, 0, max_pardepth, computebasis_pre,
			   computebasis_post, &cd);

  uninit_amatrix(cd.Z);

  freemem(cd.eps);
  freemem(cd.old2new);
  freemem(cd.cbnew);
  freemem(cd.cbold);
  freemem(cd.Z);
  if (cd.cbw) {
    for (i = 0; i < G->cb->t->desc; i++)
      uninit_amatrix(cd.cbw + i);
    freemem(cd.cbw);
  }
  if (cd.rbw) {
    for (i = 0; i < G->rb->t->desc; i++)
      uninit_amatrix(cd.rbw + i);
    freemem(cd.rbw);
  }
}

ph2matrix
build_projected_h2matrix(pch2matrix h2, pclusterbasis rb,
			 pcclusteroperator ro, pclusterbasis cb,
			 pcclusteroperator co)
{
  ph2matrix h2new, h2new1;
  pclusterbasis rb1, cb1;
  pcclusteroperator ro1, co1;
  uint      rsons, csons;
  uint      i, j;

  h2new = 0;
  if (h2->son) {
    rsons = h2->rsons;
    csons = h2->csons;

    h2new = new_super_h2matrix(rb, cb, rsons, csons);

    for (j = 0; j < csons; j++) {
      cb1 = cb;
      co1 = co;
      if (h2->son[j * rsons]->cb != h2->cb) {
	assert(j < cb->sons);
	cb1 = cb->son[j];

	co1 = 0;
	if (co) {
	  assert(j < co->sons);
	  co1 = co->son[j];
	}
      }
      for (i = 0; i < rsons; i++) {
	rb1 = rb;
	ro1 = ro;
	if (h2->son[i]->rb != h2->rb) {
	  assert(i < rb->sons);
	  rb1 = rb->son[i];

	  ro1 = 0;
	  if (ro) {
	    assert(i < ro->sons);
	    ro1 = ro->son[i];
	  }
	}

	h2new1 =
	  build_projected_h2matrix(h2->son[i + j * rsons], rb1, ro1, cb1,
				   co1);
	ref_h2matrix(h2new->son + i + j * rsons, h2new1);
      }
    }
  }
  else if (h2->u) {
    h2new = new_uniform_h2matrix(rb, cb);

    clear_amatrix(&h2new->u->S);
    add_projected_uniform(h2->u, ro, co, h2new->u);
  }
  else if (h2->f) {
    h2new = new_full_h2matrix(rb, cb);

    copy_amatrix(false, h2->f, h2new->f);
  }
  else
    h2new = new_zero_h2matrix(rb, cb);

  update_h2matrix(h2new);

  return h2new;
}

typedef struct _projectiondata projectiondata;
typedef projectiondata *pprojectiondata;
struct _projectiondata {
  pclusterbasis *rbn;
  pclusterbasis *cbn;
  pclusteroperator *ron;
  pclusteroperator *con;
};

static void
project_inplace(ph2matrix G, uint mname, uint rname, uint cname,
		uint pardepth, void *data)
{
  pprojectiondata pd = (pprojectiondata) data;
  pclusterbasis rb = pd->rbn[rname];
  pclusterbasis cb = pd->cbn[cname];
  pcclusteroperator ro = pd->ron[rname];
  pcclusteroperator co = pd->con[cname];

  /* To stop gcc from complaining about unused parameters */
  (void) mname;
  (void) pardepth;

  if (G->u) {
    assert(G->rb == G->u->rb);
    assert(G->cb == G->u->cb);

    project_inplace_uniform(G->u, rb, ro, cb, co);
  }

  ref_clusterbasis(&G->rb, rb);
  ref_clusterbasis(&G->cb, cb);
}

static void
project_parallel_inplace_h2matrix(ph2matrix G, uint pardepth,
				  pclusterbasis rb, pcclusteroperator ro,
				  pclusterbasis cb, pcclusteroperator co)
{
  projectiondata pd;

  pd.rbn = enumerate_clusterbasis(G->rb->t, rb);
  pd.cbn = enumerate_clusterbasis(G->cb->t, cb);
  pd.ron = enumerate_clusteroperator(G->rb->t, (pclusteroperator) ro);
  pd.con = enumerate_clusteroperator(G->cb->t, (pclusteroperator) co);

  iterate_h2matrix(G, 0, 0, 0, pardepth, 0, project_inplace, &pd);

  freemem(pd.con);
  freemem(pd.ron);
  freemem(pd.cbn);
  freemem(pd.rbn);
}

void
project_inplace_h2matrix(ph2matrix G, pclusterbasis rb,
			 pcclusteroperator ro, pclusterbasis cb,
			 pcclusteroperator co)
{
  project_parallel_inplace_h2matrix(G, max_pardepth, rb, ro, cb, co);
}

/* ------------------------------------------------------------
 * Specialized recompression routines for H^2-matrix arithmetic algorithms
 * ------------------------------------------------------------ */

/* compute the R of QR decomposition of the clusterbasis */
void
orthoweight_clusterbasis(pclusterbasis cb)
{
  uint      sons = cb->sons;
  pclusterbasis *son = cb->son;

  amatrix   tmp1, tmp2;
  avector   tmp3;
  pamatrix  Vhat, Vhat1;
  uint      m, off, roff;
  pavector  tau;
  uint      i, refl;

  if (sons > 0) {
    m = 0;
    roff = 0;
    for (i = 0; i < sons; i++) {
      orthoweight_clusterbasis(son[i]);
      m += son[i]->Z->rows;
      roff += son[i]->t->size;
    }
    assert(m <= roff);
    assert(roff == cb->t->size);

    Vhat = init_amatrix(&tmp1, m, cb->k);

    off = 0;
    for (i = 0; i < sons; i++) {
      Vhat1 = init_sub_amatrix(&tmp2, Vhat, son[i]->Z->rows, off, cb->k, 0);

      clear_amatrix(Vhat1);
      addmul_amatrix(1.0, false, son[i]->Z, false, &cb->son[i]->E, Vhat1);

      uninit_amatrix(Vhat1);

      off += son[i]->Z->rows;
    }
    assert(off == m);
  }
  else {
    assert(sons == 0);
    m = cb->t->size;

    Vhat = init_amatrix(&tmp1, m, cb->k);

    copy_amatrix(false, &cb->V, Vhat);
  }

  refl = UINT_MIN(m, cb->k);

  tau = init_avector(&tmp3, refl);
  qrdecomp_amatrix(Vhat, tau);

  assert(cb->Z == NULL);
  cb->Z = new_amatrix(refl, cb->k);
  copy_upper_amatrix(Vhat, false, cb->Z);

  uninit_avector(tau);
  uninit_amatrix(Vhat);
}

/* compute the totalweights of the row clusterbasis */
void
totalweight_row_clusteroperator(pclusterbasis rb, pclusteroperator rw,
				pctruncmode tm)
{
  uint      sons = rw->sons;

  pamatrix  Yhat, Yhat1;
  pamatrix  Z;			/* Z = Z or Z = R_{si} */
  amatrix   tmp1, tmp2;
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
	rows += Z->rows;
	u = u->rnext;
      }

      /* cols of Yhat */
      cols = son->k;

      Yhat = init_amatrix(&tmp1, rows, cols);

      assert(rw->kcol == son->E.cols);
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
	Z = u->cb->Z;

	assert(Z->cols == u->S.cols);
	Yhat1 = init_sub_amatrix(&tmp2, Yhat, Z->rows, off, son->k, 0);
	clear_amatrix(Yhat1);
	addmul_amatrix(alpha, false, Z, true, &u->S, Yhat1);
	uninit_amatrix(Yhat1);

	off += Z->rows;
	u = u->rnext;
      }
      assert(off == rows);

      refl = UINT_MIN(rows, cols);

      tau = init_avector(&tmp4, refl);
      qrdecomp_amatrix(Yhat, tau);

      resize_clusteroperator(rw->son[i], refl, cols);
      copy_upper_amatrix(Yhat, false, &rw->son[i]->C);

      uninit_avector(tau);
      uninit_amatrix(Yhat);

      totalweight_row_clusteroperator(son, rw->son[i], tm);
    }
  }
  else {
    assert(sons == 0);
  }
}

/* compute the totalweights of the extended col clusterbasis */
void
totalweight_col_clusteroperator(pclusterbasis cb, pclusteroperator cw,
				pctruncmode tm)
{
  uint      sons = cw->sons;

  pamatrix  Yhat, Yhat1;
  pamatrix  Z;			/* Z = Z or Z = R_{si} */
  amatrix   tmp1, tmp2;
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
	rows += Z->rows;
	u = u->cnext;
      }

      /* cols of Yhat */
      cols = son->k;

      Yhat = init_amatrix(&tmp1, rows, cols);

      assert(cw->kcol == son->E.cols);
      Yhat1 = init_sub_amatrix(&tmp2, Yhat, cw->krow, 0, son->k, 0);
      clear_amatrix(Yhat1);
      addmul_amatrix(zeta_age, false, &cw->C, true, &son->E, Yhat1);
      uninit_amatrix(Yhat1);

      off = cw->krow;
      u = son->clist;
      while (u != NULL) {
	/* Compute block weight if required */
	alpha = 1.0;
	if (tm && tm->blocks) {
	  if (tm->frobenius)
	    norm = normfrob_rkupdate_uniform(u, 0);
	  else
	    norm = norm2_rkupdate_uniform(u, 0);

	  alpha = (norm > 0.0 ? 1.0 / norm : 1.0);
	}
	Z = u->rb->Z;

	assert(Z->cols == u->S.rows);
	Yhat1 = init_sub_amatrix(&tmp2, Yhat, Z->rows, off, son->k, 0);
	clear_amatrix(Yhat1);
	addmul_amatrix(alpha, false, Z, false, &u->S, Yhat1);
	uninit_amatrix(Yhat1);

	off += Z->rows;
	u = u->cnext;
      }
      assert(off == rows);

      refl = UINT_MIN(rows, cols);

      tau = init_avector(&tmp4, refl);
      qrdecomp_amatrix(Yhat, tau);

      resize_clusteroperator(cw->son[i], refl, cols);
      copy_upper_amatrix(Yhat, false, &cw->son[i]->C);

      uninit_avector(tau);
      uninit_amatrix(Yhat);

      totalweight_col_clusteroperator(son, cw->son[i], tm);
    }
  }
  else {
    assert(sons == 0);
  }
}

/* compute adaptive clusterbasis for extended clusterbasis */
void
truncate_inplace_clusterbasis(pclusterbasis cb, pclusteroperator cw,
			      pctruncmode tm, real eps)
{
  amatrix   tmp1, tmp2, tmp3;
  realavector tmp4;
  pamatrix  Vhat, Vhat1, VhatZ, Q, Q1;
  prealavector sigma;
  pclusteroperator cw1;
  real      zeta_level;
  uint      i, off, m, k;

  zeta_level = (tm ? tm->zeta_level : 1.0);

  Vhat = 0;
  if (cb->sons == 0) {
    /* In son clusters, we have Vhat = (V A) */
    m = cb->t->size;
    Vhat = init_amatrix(&tmp1, cb->t->size, cb->k);
    copy_amatrix(false, &cb->V, Vhat);
  }
  else {
    /* Compute cluster bases for son clusters recursively */
    m = 0;
    off = 0;
    for (i = 0; i < cb->sons; i++) {
      assert(i < cw->sons);
      cw1 = cw->son[i];

      truncate_inplace_clusterbasis(cb->son[i], cw1, tm, eps * zeta_level);

      m += cb->son[i]->k;
    }

    Vhat = init_amatrix(&tmp1, m, cb->k);

    /* Blocks of Vhat are (R_{t'|...}E_{t'}   (R_{t'|...}) for sons t' */
    off = 0;
    for (i = 0; i < cb->sons; i++) {
      Vhat1 = init_sub_amatrix(&tmp2, Vhat, cb->son[i]->k, off, cb->k, 0);
      clear_amatrix(Vhat1);
      addmul_amatrix(1.0, false, &cw->son[i]->C, false, &cb->son[i]->E,
		     Vhat1);
      uninit_amatrix(Vhat1);

      off += cb->son[i]->k;
    }
    assert(off == m);
  }

  VhatZ = 0;

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
}

void
recompress_inplace_h2matrix(ph2matrix G, pctruncmode tm, real eps)
{
  pclusterbasis rb = G->rb;
  pclusterbasis cb = G->cb;

  pclusterbasis rbnew, cbnew;
  pclusteroperator rw, cw;

  if (rb == cb) {
    rbnew = clone_clusterbasis(rb);
    cbnew = rbnew;

    rw = build_from_clusterbasis_clusteroperator(rb);
    cw = rw;

    orthoweight_clusterbasis(rb);

    totalweight_row_clusteroperator(rb, rw, tm);

    clear_weight_clusterbasis(rb);

    truncate_inplace_clusterbasis(rbnew, rw, tm, eps);
  }
  else {
    rw = build_from_clusterbasis_clusteroperator(rb);
    cw = build_from_clusterbasis_clusteroperator(cb);

    rbnew = clone_clusterbasis(rb);
    cbnew = clone_clusterbasis(cb);

    orthoweight_clusterbasis(rb);
    orthoweight_clusterbasis(cb);

    totalweight_row_clusteroperator(rb, rw, tm);
    totalweight_col_clusteroperator(cb, cw, tm);

    clear_weight_clusterbasis(rb);
    clear_weight_clusterbasis(cb);

    truncate_inplace_clusterbasis(rbnew, rw, tm, eps);
    truncate_inplace_clusterbasis(cbnew, cw, tm, eps);
  }

  project_inplace_h2matrix(G, rbnew, rw, cbnew, cw);

  del_clusteroperator(rw);
  if (rw != cw) {
    del_clusteroperator(cw);
  }
}

/* ------------------------------------------------------------
 * Unify multiple cluster bases
 * ------------------------------------------------------------ */

ptruncblock
new_truncblock(pcclusterbasis cb, pcclusteroperator cw, ptruncblock next)
{
  ptruncblock tb;

  tb = (ptruncblock) allocmem((size_t) sizeof(truncblock));
  tb->cb = cb;
  tb->cw = cw;
  tb->old2new = 0;
  tb->cw_factor = 1.0;
  tb->next = next;

  return tb;
}

static    ptruncblock
new_sub_truncblock(pcclusterbasis cb, pcclusteroperator cw,
		   field cw_factor, pccluster t, uint off, ptruncblock next)
{
  ptruncblock tb;

  tb = (ptruncblock) allocmem((size_t) sizeof(truncblock));
  tb->cb = init_sub_clusterbasis(&tb->tmp_cb, (pclusterbasis) cb, t, off);
  tb->cw = cw;
  tb->old2new = 0;
  tb->cw_factor = cw_factor;
  tb->next = next;

  return tb;
}

void
del_truncblock(ptruncblock tb)
{
  ptruncblock next;

  while (tb) {
    next = tb->next;

    if (tb->cb == &tb->tmp_cb)
      uninit_clusterbasis((pclusterbasis) tb->cb);

    del_clusteroperator(tb->old2new);

    freemem(tb);

    tb = next;
  }
}

static void
del_partial_truncblock(ptruncblock tb)
{
  ptruncblock next;

  while (tb) {
    next = tb->next;

    if (tb->cb == &tb->tmp_cb)
      uninit_clusterbasis((pclusterbasis) tb->cb);

    freemem(tb);

    tb = next;
  }
}

ptruncblock
reverse_truncblock(ptruncblock tb)
{
  ptruncblock tb1, next;

  tb1 = 0;
  while (tb) {
    next = tb->next;
    tb->next = tb1;
    tb1 = tb;
    tb = next;
  }

  return tb1;
}

static    pclusterbasis
unify_parallel_clusterbasis(pccluster t, ptruncblock tb,
			    pctruncmode tm, real eps, uint pardepth,
			    pclusteroperator * cw)
{
  pclusterbasis cb, *cb1;
  pclusteroperator *cw1;
  amatrix   tmp1, tmp2, tmp3, tmp4;
  avector   tmp5;
  realavector tmp6;
  ptruncblock *tb1, tb2, tb3;
  pamatrix  Vhat, Vhat1, VhatZ, VhatZ1, W, W1, Q, Q1;
  pavector  tau;
  prealavector sigma;
  real      zeta_age, zeta_level;
  uint      sons, n, nw, m;
#ifdef USE_OPENMP
  uint      nthreads;		/* HACK: Solaris workaround */
#endif
  uint      pardepth1;
  uint      i, roff, coff, coffw, k, kw;

  zeta_age = (tm ? tm->zeta_age : 1.0);
  zeta_level = (tm ? tm->zeta_level : 1.0);

  pardepth1 = (pardepth > 0 ? pardepth - 1 : 0);

  /* Do we have to treat sons? How many columns in Vhat and VhatZ? */
  sons = 0;
  n = nw = 0;
  for (tb2 = tb; tb2; tb2 = tb2->next) {
    assert(tb2->cb->t == t);

    if (tb2->cb->sons > 0) {
      assert(tb2->cb->sons == t->sons);
      sons = tb2->cb->sons;
    }

    n += tb2->cb->k;

    kw = (tb2->cw ? tb2->cw->krow : tb2->cb->k);
    nw += kw;
  }

  /* Set up clusterbasis */
  cb = (sons > 0 ? new_clusterbasis(t) : new_leaf_clusterbasis(t));
  assert(cb->sons == sons);

  /* Set up clusteroperator for total weights */
  if (cw) {
    *cw = (sons > 0 ? new_clusteroperator(t) : new_leaf_clusteroperator(t));
    assert((*cw)->sons == sons);
  }

  /* Set up clusteroperators for basis change */
  if (sons > 0)
    for (tb2 = tb; tb2; tb2 = tb2->next)
      tb2->old2new = new_clusteroperator(t);
  else
    for (tb2 = tb; tb2; tb2 = tb2->next)
      tb2->old2new = new_leaf_clusteroperator(t);

  Vhat = 0;
  if (sons > 0) {
    /* Build truncblock structures for all sons */
    tb1 = (ptruncblock *) allocmem((size_t) sizeof(ptruncblock) * sons);
    for (i = 0; i < sons; i++)
      tb1[i] = 0;

    for (tb2 = tb; tb2; tb2 = tb2->next) {
      if (tb2->cb->son) {
	assert(tb2->cb->sons == sons);
	assert(tb2->cw == 0 || tb2->cw->sons == sons);

	for (i = 0; i < sons; i++)
	  tb1[i] = new_truncblock(tb2->cb->son[i],
				  (tb2->cw == 0 ? 0 : tb2->cw->son[i]),
				  tb1[i]);
      }
      else {
	roff = 0;
	for (i = 0; i < sons; i++) {
	  tb1[i] = new_sub_truncblock(tb2->cb, tb2->cw,
				      tb2->cw_factor * zeta_age, t->son[i],
				      roff, tb1[i]);

	  roff += t->son[i]->size;
	}
	assert(roff == t->size);
      }
    }

    /* Reverse block order */
    for (i = 0; i < sons; i++)
      tb1[i] = reverse_truncblock(tb1[i]);

    /* Compute cluster bases for son clusters recursively */
    cb1 = (pclusterbasis *) allocmem((size_t) sizeof(pclusterbasis) * sons);
    cw1 = 0;
    if (cw)
      cw1 = (pclusteroperator *) allocmem((size_t) sizeof(pclusteroperator) *
					  sons);

#ifdef USE_OPENMP
    nthreads = sons;
    (void) nthreads;
#pragma omp parallel for if(pardepth > 0), num_threads(nthreads)
#endif
    for (i = 0; i < sons; i++)
      cb1[i] = unify_parallel_clusterbasis(t->son[i], tb1[i], tm,
					   eps * zeta_level, pardepth1,
					   (cw ? cw1 + i : 0));

    /* Set up sons of clusterbasis, determine number of rows of Vhat */
    m = 0;
    for (i = 0; i < sons; i++) {
      ref_clusterbasis(cb->son + i, cb1[i]);

      m += cb1[i]->k;
    }

    /* Set up sons of weight clusteroperator */
    if (cw) {
      assert(cw1);

      for (i = 0; i < sons; i++)
	ref_clusteroperator((*cw)->son + i, cw1[i]);

      freemem(cw1);
    }

    /* Set up sons of basis change clusteroperator */
    for (i = 0; i < sons; i++) {
      for (tb2 = tb, tb3 = tb1[i]; tb2; tb2 = tb2->next, tb3 = tb3->next) {
	assert(tb2->cb->son == 0 || tb2->cb->son[i] == tb3->cb);
	assert(i < tb2->old2new->sons);
	ref_clusteroperator(tb2->old2new->son + i, tb3->old2new);
      }
      assert(tb3 == 0);
    }

    /* Clean up auxiliary data */
    freemem(cb1);
    for (i = 0; i < sons; i++)
      del_partial_truncblock(tb1[i]);
    freemem(tb1);

    /* Compute Vhat */
    Vhat = init_amatrix(&tmp1, m, n);
    coff = 0;
    for (tb2 = tb; tb2; tb2 = tb2->next) {
      if (tb2->cb->sons == 0) {
	/* tb2->cb is a leaf, so we have to remove the virtual sons
	   created by new_sub_truncblock */
	roff = 0;
	for (i = 0; i < sons; i++) {
	  Vhat1 =
	    init_sub_amatrix(&tmp2, Vhat, cb->son[i]->k, roff, tb2->cb->k,
			     coff);
	  copy_amatrix(false, &tb2->old2new->son[i]->C, Vhat1);
	  uninit_amatrix(Vhat1);

	  roff += cb->son[i]->k;
	}
	assert(roff == m);

	removesons_clusteroperator(tb2->old2new);
      }
      else {
	/* tb2->cb is not a leaf, so we have to take the transfer
	   matrices into account when assembling Vhat */
	roff = 0;
	for (i = 0; i < sons; i++) {
	  Vhat1 =
	    init_sub_amatrix(&tmp2, Vhat, cb->son[i]->k, roff, tb2->cb->k,
			     coff);
	  clear_amatrix(Vhat1);

	  assert(cb->son[i]->k == tb2->old2new->son[i]->krow);
	  assert(tb2->cb->son[i]->k == tb2->old2new->son[i]->kcol);
	  addmul_amatrix(1.0, false, &tb2->old2new->son[i]->C, false,
			 &tb2->cb->son[i]->E, Vhat1);
	  uninit_amatrix(Vhat1);

	  roff += cb->son[i]->k;
	}
	assert(roff == m);
      }

      coff += tb2->cb->k;
    }
    assert(coff == n);
  }
  else {
    /* In leaf clusters, we have Vhat = V */
    m = t->size;
    Vhat = init_amatrix(&tmp1, m, n);
    coff = 0;
    for (tb2 = tb; tb2; tb2 = tb2->next) {
      Vhat1 = init_sub_amatrix(&tmp2, Vhat, m, 0, tb2->cb->k, coff);
      copy_amatrix(false, &tb2->cb->V, Vhat1);
      uninit_amatrix(Vhat1);

      coff += tb2->cb->k;
    }
    assert(coff == n);
  }

  /* Multiply by weight matrices */
  VhatZ = init_amatrix(&tmp2, m, nw);
  coff = coffw = 0;
  for (tb2 = tb; tb2; tb2 = tb2->next) {
    Vhat1 = init_sub_amatrix(&tmp3, Vhat, m, 0, tb2->cb->k, coff);
    if (tb2->cw) {
      assert(tb2->cw->kcol == tb2->cb->k);

      VhatZ1 = init_sub_amatrix(&tmp4, VhatZ, m, 0, tb2->cw->krow, coffw);
      clear_amatrix(VhatZ1);

      addmul_amatrix(tb2->cw_factor, false, Vhat1, true, &tb2->cw->C, VhatZ1);

      uninit_amatrix(VhatZ1);

      coffw += tb2->cw->krow;
    }
    else {
      VhatZ1 = init_sub_amatrix(&tmp4, VhatZ, m, 0, tb2->cb->k, coffw);

      copy_amatrix(false, Vhat1, VhatZ1);

      uninit_amatrix(VhatZ1);

      coffw += tb2->cb->k;
    }
    uninit_amatrix(Vhat1);

    coff += tb2->cb->k;
  }
  assert(coff == n);
  assert(coffw == nw);

  /* Compute singular value decomposition of VhatZ */
  k = UINT_MIN(m, nw);
  Q = init_amatrix(&tmp3, m, k);
  sigma = init_realavector(&tmp6, k);
  svd_amatrix(VhatZ, sigma, Q, 0);
  uninit_amatrix(VhatZ);

  /* Find appropriate rank */
  k = findrank_truncmode(tm, eps, sigma);
  uninit_realavector(sigma);

  /* Set rank of new cluster basis */
  resize_clusterbasis(cb, k);
  assert(cb->sons == sons);

  if (sons == 0) {
    assert(cb->V.rows == Q->rows);
    assert(cb->V.cols == k);

    /* Copy cluster basis matrix to new cluster basis */
    Q1 = init_sub_amatrix(&tmp2, Q, Q->rows, 0, k, 0);
    copy_amatrix(false, Q1, &cb->V);
    uninit_amatrix(Q1);
  }
  else {
    roff = 0;
    for (i = 0; i < sons; i++) {
      assert(cb->son[i]->E.rows == cb->son[i]->k);
      assert(cb->son[i]->E.cols == k);

      /* Copy transfer matrix to new cluster basis */
      Q1 = init_sub_amatrix(&tmp2, Q, cb->son[i]->k, roff, k, 0);
      copy_amatrix(false, Q1, &cb->son[i]->E);
      uninit_amatrix(Q1);

      roff += cb->son[i]->k;
    }
    assert(roff == m);
  }

  /* Compute transformation from old to new basis */
  Q1 = init_sub_amatrix(&tmp2, Q, m, 0, k, 0);
  coff = 0;
  for (tb2 = tb; tb2; tb2 = tb2->next) {
    Vhat1 = init_sub_amatrix(&tmp4, Vhat, m, 0, tb2->cb->k, coff);
    resize_clusteroperator(tb2->old2new, k, tb2->cb->k);
    clear_amatrix(&tb2->old2new->C);
    addmul_amatrix(1.0, true, Q1, false, Vhat1, &tb2->old2new->C);
    uninit_amatrix(Vhat1);

    coff += tb2->cb->k;
  }
  assert(coff == n);

  uninit_amatrix(Q1);
  uninit_amatrix(Q);
  uninit_amatrix(Vhat);

  /* Compute total weight matrix */
  if (cw) {
    W = init_amatrix(&tmp1, nw, k);
    roff = 0;
    for (tb2 = tb; tb2; tb2 = tb2->next)
      if (tb2->cw) {
	W1 = init_sub_amatrix(&tmp2, W, tb2->cw->krow, roff, k, 0);
	clear_amatrix(W1);
	addmul_amatrix(tb2->cw_factor, false, &tb2->cw->C, true,
		       &tb2->old2new->C, W1);
	uninit_amatrix(W1);
	roff += tb2->cw->krow;
      }
      else {
	W1 = init_sub_amatrix(&tmp2, W, tb2->cb->k, roff, k, 0);
	copy_amatrix(true, &tb2->old2new->C, W1);
	uninit_amatrix(W1);
	roff += tb2->cb->k;
      }
    assert(roff == nw);

    tau = init_avector(&tmp5, k);
    qrdecomp_amatrix(W, tau);

    resize_clusteroperator(*cw, UINT_MIN(k, nw), k);
    copy_upper_amatrix(W, false, &(*cw)->C);

    uninit_avector(tau);
    uninit_amatrix(W);
  }

  return cb;
}

pclusterbasis
unify_clusterbasis(pccluster t, ptruncblock tb, pctruncmode tm,
		   real eps, pclusteroperator * cw)
{
  return unify_parallel_clusterbasis(t, tb, tm, eps, max_pardepth, cw);
}

/* ------------------------------------------------------------
 * Combine H^2-submatrices
 * ------------------------------------------------------------ */

void
unify_parallel_h2matrix(ph2matrix G, uint pardepth, pclusteroperator * rw1,
			pclusteroperator * cw1, pctruncmode tm, real eps,
			pclusteroperator * rw, pclusteroperator * cw)
{
  ptruncblock tb, tb1;
  pclusterbasis rb, cb;
  pclusterbasis *rb1;		/* Unified row bases */
  pclusterbasis *cb1;		/* Unified column bases */
  pclusteroperator *rw2;	/* Total weights for row bases */
  pclusteroperator *cw2;	/* Total weights for column bases */
  pclusteroperator *ro;		/* Basis change for row bases */
  pclusteroperator *co;		/* Basis change for column bases */
  uint      rsons, csons;
  uint      i, j;

  rsons = G->rsons;
  csons = G->csons;

  /* Construct unified row bases */
  rb1 = (pclusterbasis *) allocmem((size_t) sizeof(pclusterbasis) * rsons);
  rw2 = (pclusteroperator *) allocmem((size_t) sizeof(pclusteroperator) *
				      rsons);
  ro =
    (pclusteroperator *) allocmem((size_t) sizeof(pclusteroperator) * rsons *
				  csons);
  for (i = 0; i < rsons; i++) {
    tb = 0;
    for (j = 0; j < csons; j++)
      tb = new_truncblock(G->son[i + j * rsons]->rb, rw1[i + j * rsons], tb);
    tb = reverse_truncblock(tb);

    rb1[i] = unify_parallel_clusterbasis(G->son[i]->rb->t, tb, tm, eps,
					 pardepth, rw2 + i);
    assert(rb1[i]->k == rw2[i]->kcol);

    for (j = 0, tb1 = tb; j < csons; j++, tb1 = tb1->next)
      ro[i + j * rsons] = tb1->old2new;

    del_partial_truncblock(tb);
  }

  /* Construct unified column bases */
  cb1 = (pclusterbasis *) allocmem((size_t) sizeof(pclusterbasis) * csons);
  cw2 = (pclusteroperator *) allocmem((size_t) sizeof(pclusteroperator) *
				      csons);
  co =
    (pclusteroperator *) allocmem((size_t) sizeof(pclusteroperator) * rsons *
				  csons);
  for (j = 0; j < csons; j++) {
    tb = 0;
    for (i = 0; i < rsons; i++)
      tb = new_truncblock(G->son[i + j * rsons]->cb, cw1[i + j * rsons], tb);
    tb = reverse_truncblock(tb);

    cb1[j] =
      unify_parallel_clusterbasis(G->son[j * rsons]->cb->t, tb, tm, eps,
				  pardepth, cw2 + j);
    assert(cb1[j]->k == cw2[j]->kcol);

    for (i = 0, tb1 = tb; i < rsons; i++, tb1 = tb1->next)
      co[i + j * rsons] = tb1->old2new;

    del_partial_truncblock(tb);
  }

  /* Change bases */
  if (pardepth > 0) {
    for (j = 0; j < csons; j++)
      for (i = 0; i < rsons; i++)
	project_parallel_inplace_h2matrix(G->son[i + j * rsons],
					  (pardepth > 0 ? pardepth - 1 : 0),
					  rb1[i], ro[i + j * rsons], cb1[j],
					  co[i + j * rsons]);
  }
  else {
    for (j = 0; j < csons; j++)
      for (i = 0; i < rsons; i++)
	project_inplace_h2matrix(G->son[i + j * rsons], rb1[i],
				 ro[i + j * rsons], cb1[j],
				 co[i + j * rsons]);
  }

  /* Create row basis for root cluster */
  rb = 0;
  if (G->son[0]->rb->t == G->rb->t) {
    assert(rsons == 1);
    rb = rb1[0];

    if (rw)
      *rw = rw2[0];
  }
  else {
    rb = new_clusterbasis(G->rb->t);
    assert(rb->sons == rsons);
    for (i = 0; i < rsons; i++)
      ref_clusterbasis(rb->son + i, rb1[i]);
    resize_clusterbasis(rb, 0);

    if (rw) {
      *rw = new_clusteroperator(rb->t);
      assert((*rw)->sons == rsons);
      for (i = 0; i < rsons; i++)
	ref_clusteroperator((*rw)->son + i, rw2[i]);
      resize_clusteroperator(*rw, 0, rb->k);
    }
    else
      for (i = 0; i < rsons; i++)
	del_clusteroperator(rw2[i]);
  }
  ref_clusterbasis(&G->rb, rb);

  /* Create column basis for root cluster */
  cb = 0;
  if (G->son[0]->cb->t == G->cb->t) {
    assert(csons == 1);
    cb = cb1[0];

    if (cw)
      *cw = cw2[0];
  }
  else {
    cb = new_clusterbasis(G->cb->t);
    assert(cb->sons == csons);
    for (j = 0; j < csons; j++)
      ref_clusterbasis(cb->son + j, cb1[j]);
    resize_clusterbasis(cb, 0);

    if (cw) {
      *cw = new_clusteroperator(cb->t);
      assert((*cw)->sons == csons);
      for (j = 0; j < csons; j++)
	ref_clusteroperator((*cw)->son + j, cw2[j]);
      resize_clusteroperator(*cw, 0, cb->k);
    }
    else
      for (j = 0; j < csons; j++)
	del_clusteroperator(cw2[j]);
  }
  ref_clusterbasis(&G->cb, cb);

  /* Clean up */
  for (j = 0; j < csons; j++)
    for (i = 0; i < rsons; i++) {
      del_clusteroperator(ro[i + j * rsons]);
      del_clusteroperator(co[i + j * rsons]);
    }
  freemem(co);
  freemem(cw2);
  freemem(cb1);
  freemem(ro);
  freemem(rw2);
  freemem(rb1);
}

void
unify_h2matrix(ph2matrix G, pclusteroperator * rw1, pclusteroperator * cw1,
	       pctruncmode tm, real eps, pclusteroperator * rw,
	       pclusteroperator * cw)
{
  ptruncblock tb, tb1;
  pclusterbasis rb, cb;
  pclusterbasis *rb1;		/* Unified row bases */
  pclusterbasis *cb1;		/* Unified column bases */
  pclusteroperator *rw2;	/* Total weights for row bases */
  pclusteroperator *cw2;	/* Total weights for column bases */
  pclusteroperator *ro;		/* Basis change for row bases */
  pclusteroperator *co;		/* Basis change for column bases */
  uint      rsons, csons;
  uint      i, j;

  rsons = G->rsons;
  csons = G->csons;

  /* Construct unified row bases */
  rb1 = (pclusterbasis *) allocmem((size_t) sizeof(pclusterbasis) * rsons);
  rw2 = (pclusteroperator *) allocmem((size_t) sizeof(pclusteroperator) *
				      rsons);
  ro =
    (pclusteroperator *) allocmem((size_t) sizeof(pclusteroperator) * rsons *
				  csons);
  for (i = 0; i < rsons; i++) {
    tb = 0;
    for (j = 0; j < csons; j++)
      tb = new_truncblock(G->son[i + j * rsons]->rb, rw1[i + j * rsons], tb);
    tb = reverse_truncblock(tb);

    rb1[i] = unify_clusterbasis(G->son[i]->rb->t, tb, tm, eps, rw2 + i);
    assert(rb1[i]->k == rw2[i]->kcol);

    for (j = 0, tb1 = tb; j < csons; j++, tb1 = tb1->next)
      ro[i + j * rsons] = tb1->old2new;

    del_partial_truncblock(tb);
  }

  /* Construct unified column bases */
  cb1 = (pclusterbasis *) allocmem((size_t) sizeof(pclusterbasis) * csons);
  cw2 = (pclusteroperator *) allocmem((size_t) sizeof(pclusteroperator) *
				      csons);
  co =
    (pclusteroperator *) allocmem((size_t) sizeof(pclusteroperator) * rsons *
				  csons);
  for (j = 0; j < csons; j++) {
    tb = 0;
    for (i = 0; i < rsons; i++)
      tb = new_truncblock(G->son[i + j * rsons]->cb, cw1[i + j * rsons], tb);
    tb = reverse_truncblock(tb);

    cb1[j] =
      unify_clusterbasis(G->son[j * rsons]->cb->t, tb, tm, eps, cw2 + j);
    assert(cb1[j]->k == cw2[j]->kcol);

    for (i = 0, tb1 = tb; i < rsons; i++, tb1 = tb1->next)
      co[i + j * rsons] = tb1->old2new;

    del_partial_truncblock(tb);
  }

  /* Change bases */
  for (j = 0; j < csons; j++)
    for (i = 0; i < rsons; i++)
      project_inplace_h2matrix(G->son[i + j * rsons], rb1[i],
			       ro[i + j * rsons], cb1[j], co[i + j * rsons]);

  /* Create row basis for root cluster */
  rb = 0;
  if (G->son[0]->rb->t == G->rb->t) {
    assert(rsons == 1);
    rb = rb1[0];

    if (rw)
      *rw = rw2[0];
  }
  else {
    rb = new_clusterbasis(G->rb->t);
    assert(rb->sons == rsons);
    for (i = 0; i < rsons; i++)
      ref_clusterbasis(rb->son + i, rb1[i]);
    resize_clusterbasis(rb, 0);

    if (rw) {
      *rw = new_clusteroperator(rb->t);
      assert((*rw)->sons == rsons);
      for (i = 0; i < rsons; i++)
	ref_clusteroperator((*rw)->son + i, rw2[i]);
      resize_clusteroperator(*rw, 0, rb->k);
    }
    else
      for (i = 0; i < rsons; i++)
	del_clusteroperator(rw2[i]);
  }
  ref_clusterbasis(&G->rb, rb);

  /* Create column basis for root cluster */
  cb = 0;
  if (G->son[0]->cb->t == G->cb->t) {
    assert(csons == 1);
    cb = cb1[0];

    if (cw)
      *cw = cw2[0];
  }
  else {
    cb = new_clusterbasis(G->cb->t);
    assert(cb->sons == csons);
    for (j = 0; j < csons; j++)
      ref_clusterbasis(cb->son + j, cb1[j]);
    resize_clusterbasis(cb, 0);

    if (cw) {
      *cw = new_clusteroperator(cb->t);
      assert((*cw)->sons == csons);
      for (j = 0; j < csons; j++)
	ref_clusteroperator((*cw)->son + j, cw2[j]);
      resize_clusteroperator(*cw, 0, cb->k);
    }
    else
      for (j = 0; j < csons; j++)
	del_clusteroperator(cw2[j]);
  }
  ref_clusterbasis(&G->cb, cb);

  /* Clean up */
  for (j = 0; j < csons; j++)
    for (i = 0; i < rsons; i++) {
      del_clusteroperator(ro[i + j * rsons]);
      del_clusteroperator(co[i + j * rsons]);
    }
  freemem(co);
  freemem(cw2);
  freemem(cb1);
  freemem(ro);
  freemem(rw2);
  freemem(rb1);
}

/* ------------------------------------------------------------
 * Lossless conversion of an rkmatrix into a uniform matrix
 * ------------------------------------------------------------ */

void
convert_rkmatrix_uniform(pcrkmatrix r, puniform u, pctruncmode tm,
			 pclusteroperator * rw, pclusteroperator * cw)
{
  pccluster rc = u->rb->t;
  pccluster cc = u->cb->t;
  pclusterbasis rb, cb;
  amatrix   tmp1, tmp2, tmp3;
  avector   tmp4;
  pamatrix  Ac, Bc, R, C, W;
  pavector  tau;
  uint      rows, cols;
  uint      k, kr, kc;
  real      norm;

  rows = rc->size;
  cols = cc->size;

  assert(rows == r->A.rows);
  assert(cols == r->B.rows);
  assert(r->k == r->A.cols);
  assert(r->k == r->B.cols);

  kr = UINT_MIN(rows, r->k);
  kc = UINT_MIN(cols, r->k);
  k = UINT_MIN(kr, kc);

  rb = new_leaf_clusterbasis(rc);
  cb = new_leaf_clusterbasis(cc);
  resize_clusterbasis(rb, k);
  resize_clusterbasis(cb, k);
  ref_row_uniform(u, rb);
  ref_col_uniform(u, cb);
  resize_amatrix(&u->S, k, k);

  tau = init_avector(&tmp4, r->k);

  if (kr <= kc) {
    /* Copy A */
    Ac = init_amatrix(&tmp1, r->A.rows, r->A.cols);
    copy_amatrix(false, &r->A, Ac);

    /* Compute A = Q_A R_A */
    qrdecomp_amatrix(Ac, tau);
    R = init_sub_amatrix(&tmp3, Ac, kr, 0, r->A.cols, 0);

    /* Q_A is the row basis */
    qrexpand_amatrix(Ac, tau, &rb->V);

    /* Copy B */
    Bc = init_amatrix(&tmp2, r->B.rows, r->B.cols);
    copy_amatrix(false, &r->B, Bc);

    /* Compute C^* = R_A B^*, i.e., C = B R_A^* */
    triangulareval_amatrix(false, false, false, R, true, Bc);
    uninit_amatrix(R);
    uninit_amatrix(Ac);
    C = init_sub_amatrix(&tmp1, Bc, r->B.rows, 0, kr, 0);

    /* Compute C = Q_C R_C */
    qrdecomp_amatrix(C, tau);

    /* Q_C is the column basis */
    qrexpand_amatrix(C, tau, &cb->V);

    /* R_C is the row weight matrix... */
    W = init_amatrix(&tmp3, k, k);
    copy_upper_amatrix(C, false, W);
    if (rw) {
      *rw = new_leaf_clusteroperator(rc);
      resize_clusteroperator(*rw, k, k);
      copy_amatrix(false, W, &(*rw)->C);
    }

    /* ... and its adjoint is the column weight and the coupling matrix */
    copy_amatrix(true, W, &u->S);
    if (cw) {
      *cw = new_leaf_clusteroperator(cc);
      resize_clusteroperator(*cw, k, k);
      copy_amatrix(true, W, &(*cw)->C);
    }

    /* Clean up */
    uninit_amatrix(W);
    uninit_amatrix(C);
    uninit_amatrix(Bc);
  }
  else {
    /* Copy B */
    Bc = init_amatrix(&tmp1, r->B.rows, r->B.cols);
    copy_amatrix(false, &r->B, Bc);

    /* Compute B = Q_B R_B */
    qrdecomp_amatrix(Bc, tau);
    R = init_sub_amatrix(&tmp3, Bc, kc, 0, r->B.cols, 0);

    /* Q_B is the column basis */
    qrexpand_amatrix(Bc, tau, &cb->V);

    /* Copy A */
    Ac = init_amatrix(&tmp2, r->A.rows, r->A.cols);
    copy_amatrix(false, &r->A, Ac);

    /* Compute C^* = R_B A^*, i.e., C = A R_B^* */
    triangulareval_amatrix(false, false, false, R, true, Ac);
    uninit_amatrix(R);
    uninit_amatrix(Bc);
    C = init_sub_amatrix(&tmp1, Ac, r->A.rows, 0, kc, 0);

    /* Compute C = Q_C R_C */
    qrdecomp_amatrix(C, tau);

    /* Q_C is the row basis */
    qrexpand_amatrix(C, tau, &rb->V);

    /* R_C is the column weight matrix and the coupling matrix... */
    copy_upper_amatrix(C, false, &u->S);
    if (cw) {
      *cw = new_leaf_clusteroperator(cc);
      resize_clusteroperator(*cw, k, k);
      copy_amatrix(false, &u->S, &(*cw)->C);
    }

    /* ... and its adjoint is the column weight matrix */
    if (rw) {
      *rw = new_leaf_clusteroperator(rc);
      resize_clusteroperator(*rw, k, k);
      copy_amatrix(true, &u->S, &(*rw)->C);
    }

    /* Clean up */
    uninit_amatrix(C);
    uninit_amatrix(Ac);
  }
  uninit_avector(tau);

  /* Scale weight matrices */
  if (tm && tm->blocks) {
    norm = (tm->frobenius ? normfrob_amatrix(&u->S) : norm2_amatrix(&u->S));
    if (norm > 0.0) {
      if (*rw)
	scale_amatrix(1.0 / norm, &(*rw)->C);
      if (*cw)
	scale_amatrix(1.0 / norm, &(*cw)->C);
    }
  }
}

/* ------------------------------------------------------------
 * H-matrix blocks
 * ------------------------------------------------------------ */

typedef struct _hcompactive hcompactive;
typedef hcompactive *phcompactive;

struct _hcompactive {
  pchmatrix hm;
  amatrix   A;
  real      weight;
  phcompactive next;
};

typedef struct _hcomppassive hcomppassive;
typedef hcomppassive *phcomppassive;

struct _hcomppassive {
  pchmatrix hm;
  phcomppassive next;
};

static void
del_hcompactive(phcompactive ha)
{
  phcompactive next;

  while (ha) {
    next = ha->next;

    uninit_amatrix(&ha->A);
    freemem(ha);

    ha = next;
  }
}

static void
del_hcomppassive(phcomppassive hp)
{
  phcomppassive next;

  while (hp) {
    next = hp->next;

    freemem(hp);

    hp = next;
  }
}

static void
addrow_hcomp(pccluster rc, pchmatrix hm, pctruncmode tm,
	     phcompactive * active, phcomppassive * passive)
{
  phcompactive ha;
  phcomppassive hp;
  amatrix   tmp1, tmp2, tmp3;
  avector   tmp4;
  pamatrix  R, R1, Ahat;
  pavector  tau;
  real      norm, weight;
  uint      rsons, csons;
  uint      i, j, k;

  if (hm->son) {
    rsons = hm->rsons;
    csons = hm->csons;

    /* Check whether there is a son matching the cluster rc */
    i = 0;
    while (i < rsons && hm->son[i]->rc != rc)
      i++;

    /* If there is, check the sons recursively */
    if (i < rsons) {
      for (j = 0; j < csons; j++)
	addrow_hcomp(rc, hm->son[i + j * rsons], tm, active, passive);
    }
    else {			/* Otherwise, this matrix block is passive */
      assert(hm->rc == rc);

      hp = (phcomppassive) allocmem(sizeof(hcomppassive));
      hp->hm = hm;
      hp->next = *passive;
      *passive = hp;
    }
  }
  else if (hm->r) {
    assert(hm->r->A.rows == hm->rc->size);
    assert(hm->r->A.cols == hm->r->k);
    assert(hm->r->B.rows == hm->cc->size);
    assert(hm->r->B.cols == hm->r->k);

    /* Copy B part of the rkmatrix */
    R = init_amatrix(&tmp1, hm->cc->size, hm->r->k);
    copy_amatrix(false, &hm->r->B, R);

    /* Compute QR factorization B=QR */
    k = UINT_MIN(hm->cc->size, hm->r->k);
    tau = init_avector(&tmp4, k);
    qrdecomp_amatrix(R, tau);
    uninit_avector(tau);

    /* Multiply A by R^* */
    R1 = init_sub_amatrix(&tmp2, R, k, 0, hm->r->k, 0);
    Ahat = init_amatrix(&tmp3, hm->rc->size, hm->r->k);
    copy_amatrix(false, &hm->r->A, Ahat);
    triangulareval_amatrix(false, false, false, R1, true, Ahat);
    uninit_amatrix(R1);
    uninit_amatrix(R);

    /* Compute weight factor if necessary */
    weight = 1.0;
    if (tm && tm->blocks) {
      norm = (tm->frobenius ? normfrob_amatrix(Ahat) : norm2_amatrix(Ahat));
      if (norm > 0.0)
	weight = 1.0 / norm;
    }

    /* Copy the result to ha->A */
    ha = (phcompactive) allocmem(sizeof(hcompactive));
    ha->hm = hm;
    init_amatrix(&ha->A, hm->rc->size, k);
    copy_sub_amatrix(false, Ahat, &ha->A);
    ha->weight = weight;
    ha->next = *active;
    *active = ha;
    uninit_amatrix(Ahat);
  }
}

static void
addcol_hcomp(pccluster cc, pchmatrix hm, pctruncmode tm,
	     phcompactive * active, phcomppassive * passive)
{
  phcompactive ha;
  phcomppassive hp;
  amatrix   tmp1, tmp2, tmp3;
  avector   tmp4;
  pamatrix  R, R1, Bhat;
  pavector  tau;
  real      norm, weight;
  uint      rsons, csons;
  uint      i, j, k;

  if (hm->son) {
    rsons = hm->rsons;
    csons = hm->csons;

    /* Check whether there is a son matching the cluster cc */
    j = 0;
    while (j < csons && hm->son[j * rsons]->cc != cc)
      j++;

    /* If there is, check the sons recursively */
    if (j < csons) {
      for (i = 0; i < rsons; i++)
	addcol_hcomp(cc, hm->son[i + j * rsons], tm, active, passive);
    }
    else {			/* Otherwise, this matrix block is passive */
      assert(hm->cc == cc);

      hp = (phcomppassive) allocmem(sizeof(hcomppassive));
      hp->hm = hm;
      hp->next = *passive;
      *passive = hp;
    }
  }
  else if (hm->r) {
    assert(hm->r->A.rows == hm->rc->size);
    assert(hm->r->A.cols == hm->r->k);
    assert(hm->r->B.rows == hm->cc->size);
    assert(hm->r->B.cols == hm->r->k);

    /* Copy A part of the rkmatrix */
    R = init_amatrix(&tmp1, hm->rc->size, hm->r->k);
    copy_amatrix(false, &hm->r->A, R);

    /* Compute QR factorization A=QR */
    k = UINT_MIN(hm->rc->size, hm->r->k);
    tau = init_avector(&tmp4, k);
    qrdecomp_amatrix(R, tau);
    uninit_avector(tau);

    /* Multiply B by R^* */
    R1 = init_sub_amatrix(&tmp2, R, k, 0, hm->r->k, 0);
    Bhat = init_amatrix(&tmp3, hm->cc->size, hm->r->k);
    copy_amatrix(false, &hm->r->B, Bhat);
    triangulareval_amatrix(false, false, false, R1, true, Bhat);
    uninit_amatrix(R1);
    uninit_amatrix(R);

    /* Compute weight factor if necessary */
    weight = 1.0;
    if (tm && tm->blocks) {
      norm = (tm->frobenius ? normfrob_amatrix(Bhat) : norm2_amatrix(Bhat));
      if (norm > 0.0)
	weight = 1.0 / norm;
    }

    /* Copy the result to ha->A */
    ha = (phcompactive) allocmem(sizeof(hcompactive));
    ha->hm = hm;
    init_amatrix(&ha->A, hm->cc->size, k);
    copy_sub_amatrix(false, Bhat, &ha->A);
    ha->weight = weight;
    ha->next = *active;
    *active = ha;
    uninit_amatrix(Bhat);
  }
}

static    pclusterbasis
buildbasis_hcomp(pccluster t, bool colbasis,
		 phcompactive active, phcomppassive passive, pctruncmode tm,
		 real eps)
{
  pclusterbasis cb, cb1;
  amatrix   tmp1, tmp2, tmp3, tmp4;
  realavector tmp5;
  pamatrix  Ahat, Ahat0, Ahat1;
  pamatrix  Q, Q1;
  prealavector sigma;
  phcompactive active1, ha, ha1;
  phcomppassive passive1, hp;
  real      zeta_age, zeta_level;
  uint      i, off, m, n, k;

  zeta_age = (tm ? tm->zeta_age : 1.0);
  zeta_level = (tm ? tm->zeta_level : 1.0);

  cb = new_clusterbasis(t);

  if (cb->sons > 0) {
    assert(cb->sons == t->sons);

    off = 0;
    m = 0;
    for (i = 0; i < t->sons; i++) {
      active1 = 0;
      passive1 = 0;

      /* Check for passive blocks that become active in son */
      if (colbasis) {
	for (hp = passive; hp; hp = hp->next)
	  addcol_hcomp(t->son[i], hp->hm, tm, &active1, &passive1);
      }
      else {
	for (hp = passive; hp; hp = hp->next)
	  addrow_hcomp(t->son[i], hp->hm, tm, &active1, &passive1);
      }

      /* Add submatrices for already active blocks to list */
      for (ha = active; ha; ha = ha->next) {
	ha1 = (phcompactive) allocmem(sizeof(hcompactive));
	ha1->hm = ha->hm;
	init_sub_amatrix(&ha1->A, &ha->A, t->son[i]->size, off, ha->A.cols,
			 0);
	ha1->weight = ha->weight * zeta_age;
	ha1->next = active1;
	active1 = ha1;
      }

      /* Create cluster basis for son */
      cb1 = buildbasis_hcomp(t->son[i], colbasis, active1, passive1, tm,
			     eps * zeta_level);
      ref_clusterbasis(cb->son + i, cb1);

      /* Clean up block lists */
      del_hcompactive(active1);
      del_hcomppassive(passive1);

      off += t->son[i]->size;
      m += cb1->k;
    }
    assert(off == t->size);

    for (ha = active; ha; ha = ha->next) {
      Ahat = init_amatrix(&tmp1, m, ha->A.cols);

      /* Combine projections of son's submatrices */
      off = 0;
      m = 0;
      for (i = 0; i < t->sons; i++) {
	/* Matrix prepared by recursive call */
	Ahat1 =
	  init_sub_amatrix(&tmp2, &ha->A, cb->son[i]->k, off, Ahat->cols, 0);

	/* Appropriate submatrix in Ahat */
	Ahat0 =
	  init_sub_amatrix(&tmp3, Ahat, cb->son[i]->k, m, Ahat->cols, 0);

	/* Copy submatrix to its place */
	copy_amatrix(false, Ahat1, Ahat0);
	uninit_amatrix(Ahat1);
	uninit_amatrix(Ahat0);

	off += t->son[i]->size;
	m += cb->son[i]->k;
      }
      assert(off == t->size);
      assert(m == Ahat->rows);

      copy_sub_amatrix(false, Ahat, &ha->A);
      uninit_amatrix(Ahat);
    }
  }
  else
    m = t->size;

  /* Determine number of columns of SVD matrix */
  n = 0;
  for (ha = active; ha; ha = ha->next)
    n += ha->A.cols;

  /* Quick exit if zero columns */
  if (n == 0) {
    resize_clusterbasis(cb, 0);

    return cb;
  }

  /* Combine submatrices */
  Ahat = init_amatrix(&tmp1, m, n);
  n = 0;
  for (ha = active; ha; ha = ha->next) {
    Ahat0 = init_sub_amatrix(&tmp2, Ahat, m, 0, ha->A.cols, n);

    copy_sub_amatrix(false, &ha->A, Ahat0);
    if (ha->weight != 1.0)
      scale_amatrix(ha->weight, Ahat0);
    uninit_amatrix(Ahat0);

    n += ha->A.cols;
  }
  assert(n == Ahat->cols);

  /* Compute singular value decomposition */
  k = UINT_MIN(m, n);
  Q = init_amatrix(&tmp2, m, n);
  sigma = init_realavector(&tmp5, k);
  svd_amatrix(Ahat, sigma, Q, 0);
  uninit_amatrix(Ahat);

  /* Find appropriate rank */
  k = findrank_truncmode(tm, eps, sigma);
  uninit_realavector(sigma);

  /* Set rank of new cluster basis */
  resize_clusterbasis(cb, k);

  if (cb->sons == 0) {
    assert(cb->V.rows == Q->rows);
    assert(cb->V.cols == k);

    /* Copy cluster basis matrix to new cluster basis */
    Q1 = init_sub_amatrix(&tmp1, Q, Q->rows, 0, k, 0);
    copy_amatrix(false, Q1, &cb->V);
    uninit_amatrix(Q1);
  }
  else {
    off = 0;
    for (i = 0; i < cb->sons; i++) {
      assert(cb->son[i]->E.rows == cb->son[i]->k);
      assert(cb->son[i]->E.cols == k);

      /* Copy transfer matrix to new cluster basis */
      Q1 = init_sub_amatrix(&tmp1, Q, cb->son[i]->k, off, k, 0);
      copy_amatrix(false, Q1, &cb->son[i]->E);
      uninit_amatrix(Q1);

      off += cb->son[i]->k;
    }
    assert(off == m);
  }

  /* Prepare projection matrices for father */
  Q1 = init_sub_amatrix(&tmp1, Q, m, 0, k, 0);
  for (ha = active; ha; ha = ha->next) {
    /* Copy original matrix */
    Ahat = init_amatrix(&tmp3, m, ha->A.cols);
    copy_sub_amatrix(false, &ha->A, Ahat);

    /* Pick submatrix for result */
    Ahat0 = init_sub_amatrix(&tmp4, &ha->A, k, 0, ha->A.cols, 0);
    clear_amatrix(Ahat0);

    /* Compute projection */
    addmul_amatrix(1.0, true, Q1, false, Ahat, Ahat0);
    uninit_amatrix(Ahat0);
    uninit_amatrix(Ahat);
  }
  uninit_amatrix(Q1);
  uninit_amatrix(Q);

  /* Now we're done with this subtree */
  update_clusterbasis(cb);

  return cb;
}

pclusterbasis
buildrowbasis_hmatrix(pchmatrix G, pctruncmode tm, real eps)
{
  pclusterbasis rb;
  phcompactive active;
  phcomppassive passive;

  active = 0;
  passive = 0;
  addrow_hcomp(G->rc, G, tm, &active, &passive);

  rb = buildbasis_hcomp(G->rc, false, active, passive, tm, eps);

  del_hcompactive(active);
  del_hcomppassive(passive);

  return rb;
}

pclusterbasis
buildcolbasis_hmatrix(pchmatrix G, pctruncmode tm, real eps)
{
  pclusterbasis cb;
  phcompactive active;
  phcomppassive passive;

  active = 0;
  passive = 0;
  addcol_hcomp(G->cc, G, tm, &active, &passive);

  cb = buildbasis_hcomp(G->cc, true, active, passive, tm, eps);

  del_hcompactive(active);
  del_hcomppassive(passive);

  return cb;
}

/* ------------------------------------------------------------
 * Approximate H-matrix in new cluster bases
 * ------------------------------------------------------------ */

ph2matrix
build_projected_hmatrix_h2matrix(pchmatrix G, pclusterbasis rb,
				 pclusterbasis cb)
{
  ph2matrix G2, G21;
  pclusterbasis rb1, cb1;
  uint      rsons, csons;
  uint      i, j;

  assert(G->rc == rb->t);
  assert(G->cc == cb->t);

  if (G->son) {
    rsons = G->rsons;
    csons = G->csons;

    G2 = new_super_h2matrix(rb, cb, rsons, csons);

    for (j = 0; j < csons; j++) {
      cb1 = cb;
      if (G->son[j * rsons]->cc != cb->t) {
	assert(j < cb->sons);
	cb1 = cb->son[j];
      }

      for (i = 0; i < rsons; i++) {
	rb1 = rb;
	if (G->son[i]->rc != rb->t) {
	  assert(i < rb->sons);
	  rb1 = rb->son[i];
	}

	G21 =
	  build_projected_hmatrix_h2matrix(G->son[i + j * rsons], rb1, cb1);
	ref_h2matrix(G2->son + i + j * rsons, G21);
      }
    }
  }
  else if (G->f) {
    G2 = new_full_h2matrix(rb, cb);

    copy_amatrix(false, G->f, G2->f);
  }
  else if (G->r && G->r->A.cols > 0) {
    G2 = new_uniform_h2matrix(rb, cb);

    clear_uniform(G2->u);
    add_rkmatrix_uniform(G->r, G2->u);
  }
  else
    G2 = new_zero_h2matrix(rb, cb);

  update_h2matrix(G2);

  return G2;
}

/* ------------------------------------------------------------
 * Dense matrix blocks
 * ------------------------------------------------------------ */

typedef struct _compactive compactive;
typedef compactive *pcompactive;

struct _compactive {
  pcamatrix G;
  pcblock   b;
  amatrix   A;
  real      weight;
  pcompactive next;
};

typedef struct _comppassive comppassive;
typedef comppassive *pcomppassive;

struct _comppassive {
  pcamatrix G;
  pcblock   b;
  pcomppassive next;
};

static void
del_compactive(pcompactive ca)
{
  pcompactive next;

  while (ca) {
    next = ca->next;

    uninit_amatrix(&ca->A);
    freemem(ca);

    ca = next;
  }
}

static void
del_comppassive(pcomppassive cp)
{
  pcomppassive next;

  while (cp) {
    next = cp->next;

    freemem(cp);

    cp = next;
  }
}

static void
addrow_comp(pccluster rc, pcamatrix G, pcblock b, pctruncmode tm,
	    pcompactive * active, pcomppassive * passive)
{
  pcompactive ca;
  pcomppassive cp;
  pamatrix  Ahat;
  const uint *ridx, *cidx;
  longindex ldAhat, ldG;
  uint      rsize, csize;
  real      norm, weight;
  uint      rsons, csons;
  uint      i, j;

  if (b->son) {
    rsons = b->rsons;
    csons = b->csons;

    /* Check whether there is a son matching the cluster rc */
    i = 0;
    while (i < rsons && b->son[i]->rc != rc)
      i++;

    /* If there is, check the sons recursively */
    if (i < rsons) {
      for (j = 0; j < csons; j++)
	addrow_comp(rc, G, b->son[i + j * rsons], tm, active, passive);
    }
    else {			/* Otherwise, this matrix block is passive */
      assert(b->rc == rc);

      cp = (pcomppassive) allocmem(sizeof(comppassive));
      cp->G = G;
      cp->b = b;
      cp->next = *passive;
      *passive = cp;
    }
  }
  else if (b->a) {
    assert(b->rc == rc);

    /* Create an active block */
    rsize = b->rc->size;
    csize = b->cc->size;
    ca = (pcompactive) allocmem(sizeof(compactive));
    ca->G = G;
    ca->b = b;
    Ahat = init_amatrix(&ca->A, rsize, csize);
    ca->next = *active;
    *active = ca;

    /* Copy entries from original matrix G */
    ldAhat = Ahat->ld;
    ldG = G->ld;
    ridx = b->rc->idx;
    cidx = b->cc->idx;
    for (j = 0; j < csize; j++)
      for (i = 0; i < rsize; i++)
	Ahat->a[i + j * ldAhat] = G->a[ridx[i] + cidx[j] * ldG];

    /* Compute weight factor if necessary */
    weight = 1.0;
    if (tm && tm->blocks) {
      norm = (tm->frobenius ? normfrob_amatrix(Ahat) : norm2_amatrix(Ahat));
      if (norm > 0.0)
	weight = 1.0 / norm;
    }
    ca->weight = weight;
  }
}

static void
addcol_comp(pccluster cc, pcamatrix G, pcblock b, pctruncmode tm,
	    pcompactive * active, pcomppassive * passive)
{
  pcompactive ca;
  pcomppassive cp;
  pamatrix  Bhat;
  const uint *ridx, *cidx;
  longindex ldBhat, ldG;
  uint      rsize, csize;
  real      norm, weight;
  uint      rsons, csons;
  uint      i, j;

  if (b->son) {
    rsons = b->rsons;
    csons = b->csons;

    /* Check whether there is a son matching the cluster cc */
    j = 0;
    while (j < csons && b->son[j * rsons]->cc != cc)
      j++;

    /* If there is, check the sons recursively */
    if (j < csons) {
      for (i = 0; i < rsons; i++)
	addcol_comp(cc, G, b->son[i + j * rsons], tm, active, passive);
    }
    else {			/* Otherwise, this matrix block is passive */
      assert(b->cc == cc);

      cp = (pcomppassive) allocmem(sizeof(comppassive));
      cp->G = G;
      cp->b = b;
      cp->next = *passive;
      *passive = cp;
    }
  }
  else if (b->a) {
    assert(b->cc == cc);

    /* Create an active block */
    rsize = b->rc->size;
    csize = b->cc->size;
    ca = (pcompactive) allocmem(sizeof(compactive));
    ca->G = G;
    ca->b = b;
    Bhat = init_amatrix(&ca->A, csize, rsize);
    ca->next = *active;
    *active = ca;

    /* Copy entries from original matrix G */
    ldBhat = Bhat->ld;
    ldG = G->ld;
    ridx = b->rc->idx;
    cidx = b->cc->idx;
    for (j = 0; j < rsize; j++)
      for (i = 0; i < csize; i++)
	Bhat->a[i + j * ldBhat] = G->a[ridx[j] + cidx[i] * ldG];

    /* Compute weight factor if necessary */
    weight = 1.0;
    if (tm && tm->blocks) {
      norm = (tm->frobenius ? normfrob_amatrix(Bhat) : norm2_amatrix(Bhat));
      if (norm > 0.0)
	weight = 1.0 / norm;
    }
    ca->weight = weight;
  }
}

static    pclusterbasis
buildbasis_comp(pccluster t, bool colbasis,
		pcompactive active, pcomppassive passive, pctruncmode tm,
		real eps)
{
  pclusterbasis cb, cb1;
  amatrix   tmp1, tmp2, tmp3, tmp4;
  realavector tmp5;
  pamatrix  Ahat, Ahat0, Ahat1;
  pamatrix  Q, Q1;
  prealavector sigma;
  pcompactive active1, ca, ca1;
  pcomppassive passive1, cp;
  real      zeta_age, zeta_level;
  uint      i, off, m, n, k;

  zeta_age = (tm ? tm->zeta_age : 1.0);
  zeta_level = (tm ? tm->zeta_level : 1.0);

  cb = new_clusterbasis(t);

  if (cb->sons > 0) {
    assert(cb->sons == t->sons);

    off = 0;
    m = 0;
    for (i = 0; i < t->sons; i++) {
      active1 = 0;
      passive1 = 0;

      /* Check for passive blocks that become active in son */
      if (colbasis) {
	for (cp = passive; cp; cp = cp->next)
	  addcol_comp(t->son[i], cp->G, cp->b, tm, &active1, &passive1);
      }
      else {
	for (cp = passive; cp; cp = cp->next)
	  addrow_comp(t->son[i], cp->G, cp->b, tm, &active1, &passive1);
      }

      /* Add submatrices for already active blocks to list */
      for (ca = active; ca; ca = ca->next) {
	ca1 = (pcompactive) allocmem(sizeof(compactive));
	ca1->G = ca->G;
	ca1->b = ca->b;
	init_sub_amatrix(&ca1->A, &ca->A, t->son[i]->size, off, ca->A.cols,
			 0);
	ca1->weight = ca->weight * zeta_age;
	ca1->next = active1;
	active1 = ca1;
      }

      /* Create cluster basis for son */
      cb1 = buildbasis_comp(t->son[i], colbasis, active1, passive1, tm,
			    eps * zeta_level);
      ref_clusterbasis(cb->son + i, cb1);

      /* Clean up block lists */
      del_compactive(active1);
      del_comppassive(passive1);

      off += t->son[i]->size;
      m += cb1->k;
    }
    assert(off == t->size);

    for (ca = active; ca; ca = ca->next) {
      Ahat = init_amatrix(&tmp1, m, ca->A.cols);

      /* Combine projections of son's submatrices */
      off = 0;
      m = 0;
      for (i = 0; i < t->sons; i++) {
	/* Matrix prepared by recursive call */
	Ahat1 =
	  init_sub_amatrix(&tmp2, &ca->A, cb->son[i]->k, off, Ahat->cols, 0);

	/* Appropriate submatrix in Ahat */
	Ahat0 =
	  init_sub_amatrix(&tmp3, Ahat, cb->son[i]->k, m, Ahat->cols, 0);

	/* Copy submatrix to its place */
	copy_amatrix(false, Ahat1, Ahat0);
	uninit_amatrix(Ahat1);
	uninit_amatrix(Ahat0);

	off += t->son[i]->size;
	m += cb->son[i]->k;
      }
      assert(off == t->size);
      assert(m == Ahat->rows);

      copy_sub_amatrix(false, Ahat, &ca->A);
      uninit_amatrix(Ahat);
    }
  }
  else
    m = t->size;

  /* Determine number of columns of SVD matrix */
  n = 0;
  for (ca = active; ca; ca = ca->next)
    n += ca->A.cols;

  /* Quick exit if zero columns */
  if (n == 0) {
    resize_clusterbasis(cb, 0);

    return cb;
  }

  /* Combine submatrices */
  Ahat = init_amatrix(&tmp1, m, n);
  n = 0;
  for (ca = active; ca; ca = ca->next) {
    Ahat0 = init_sub_amatrix(&tmp2, Ahat, m, 0, ca->A.cols, n);

    copy_sub_amatrix(false, &ca->A, Ahat0);
    if (ca->weight != 1.0)
      scale_amatrix(ca->weight, Ahat0);
    uninit_amatrix(Ahat0);

    n += ca->A.cols;
  }
  assert(n == Ahat->cols);

  /* Compute singular value decomposition */
  k = UINT_MIN(m, n);
  Q = init_amatrix(&tmp2, m, n);
  sigma = init_realavector(&tmp5, k);
  svd_amatrix(Ahat, sigma, Q, 0);
  uninit_amatrix(Ahat);

  /* Find appropriate rank */
  k = findrank_truncmode(tm, eps, sigma);
  uninit_realavector(sigma);

  /* Set rank of new cluster basis */
  resize_clusterbasis(cb, k);

  if (cb->sons == 0) {
    assert(cb->V.rows == Q->rows);
    assert(cb->V.cols == k);

    /* Copy cluster basis matrix to new cluster basis */
    Q1 = init_sub_amatrix(&tmp1, Q, Q->rows, 0, k, 0);
    copy_amatrix(false, Q1, &cb->V);
    uninit_amatrix(Q1);
  }
  else {
    off = 0;
    for (i = 0; i < cb->sons; i++) {
      assert(cb->son[i]->E.rows == cb->son[i]->k);
      assert(cb->son[i]->E.cols == k);

      /* Copy transfer matrix to new cluster basis */
      Q1 = init_sub_amatrix(&tmp1, Q, cb->son[i]->k, off, k, 0);
      copy_amatrix(false, Q1, &cb->son[i]->E);
      uninit_amatrix(Q1);

      off += cb->son[i]->k;
    }
    assert(off == m);
  }

  /* Prepare projection matrices for father */
  Q1 = init_sub_amatrix(&tmp3, Q, m, 0, k, 0);
  for (ca = active; ca; ca = ca->next) {
    /* Copy original matrix */
    Ahat = init_amatrix(&tmp1, m, ca->A.cols);
    copy_sub_amatrix(false, &ca->A, Ahat);

    /* Pick submatrix for result */
    Ahat0 = init_sub_amatrix(&tmp4, &ca->A, k, 0, ca->A.cols, 0);
    clear_amatrix(Ahat0);

    /* Compute projection */
    addmul_amatrix(1.0, true, Q1, false, Ahat, Ahat0);
    uninit_amatrix(Ahat0);
    uninit_amatrix(Ahat);
  }
  uninit_amatrix(Q1);
  uninit_amatrix(Q);

  /* Now we're done with this subtree */
  update_clusterbasis(cb);

  return cb;
}

pclusterbasis
buildrowbasis_amatrix(pcamatrix G, pcblock b, pctruncmode tm, real eps)
{
  pclusterbasis rb;
  pcompactive active;
  pcomppassive passive;

  active = 0;
  passive = 0;
  addrow_comp(b->rc, G, b, tm, &active, &passive);

  rb = buildbasis_comp(b->rc, false, active, passive, tm, eps);

  del_compactive(active);
  del_comppassive(passive);

  return rb;
}

pclusterbasis
buildcolbasis_amatrix(pcamatrix G, pcblock b, pctruncmode tm, real eps)
{
  pclusterbasis cb;
  pcompactive active;
  pcomppassive passive;

  active = 0;
  passive = 0;
  addcol_comp(b->cc, G, b, tm, &active, &passive);

  cb = buildbasis_comp(b->cc, true, active, passive, tm, eps);

  del_compactive(active);
  del_comppassive(passive);

  return cb;
}

/* ------------------------------------------------------------
 * Orthogonal projection
 * ------------------------------------------------------------ */

ph2matrix
build_projected_amatrix_h2matrix(pcamatrix G, pcblock b,
				 pclusterbasis rb, pclusterbasis cb)
{
  ph2matrix h2, h21;
  pamatrix  f;
  pclusterbasis rb1, cb1;
  pccluster rc, cc;
  const uint *ridx, *cidx;
  longindex ldf, ldG;
  uint      rsize, csize;
  uint      rsons, csons;
  uint      i, j;

  if (b->son) {
    rsons = b->rsons;
    csons = b->csons;
    h2 = new_super_h2matrix(rb, cb, rsons, csons);

    for (j = 0; j < csons; j++) {
      cb1 = cb;
      if (b->son[0]->cc != b->cc) {
	assert(j < cb->sons);
	cb1 = cb->son[j];
      }

      for (i = 0; i < rsons; i++) {
	rb1 = rb;
	if (b->son[0]->rc != b->rc) {
	  assert(i < rb->sons);
	  rb1 = rb->son[i];
	}

	h21 = build_projected_amatrix_h2matrix(G, b->son[i + j * rsons], rb1,
					       cb1);
	ref_h2matrix(h2->son + i + j * rsons, h21);
      }
    }
  }
  else if (b->a) {
    h2 = new_uniform_h2matrix(rb, cb);

    collectdense_h2matrix(G, rb, cb, &h2->u->S);
  }
  else {
    h2 = new_full_h2matrix(rb, cb);
    f = h2->f;

    rc = rb->t;
    cc = cb->t;
    ridx = rc->idx;
    cidx = cc->idx;
    rsize = rc->size;
    csize = cc->size;

    assert(h2->f->rows == rc->size);
    assert(h2->f->cols == cc->size);

    ldf = f->ld;
    ldG = G->ld;
    for (j = 0; j < csize; j++)
      for (i = 0; i < rsize; i++)
	f->a[i + j * ldf] = G->a[ridx[i] + cidx[j] * ldG];
  }

  update_h2matrix(h2);

  return h2;
}
