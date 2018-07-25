
/* ------------------------------------------------------------
 * This is the file "clusterbasis.c" of the H2Lib package.
 * All rights reserved, Steffen Boerm 2009
 * ------------------------------------------------------------ */

#include "clusterbasis.h"

#include "basic.h"
#include "factorizations.h"

#ifdef USE_NETCDF
#include <stdio.h>
#include <string.h>
#include <netcdf.h>
#endif

static uint active_clusterbasis = 0;

/* ------------------------------------------------------------
 * Constructors and destructors
 * ------------------------------------------------------------ */

static void
init_raw_clusterbasis(pclusterbasis cb, pccluster t)
{
  assert(cb != NULL);

  cb->t = t;
  cb->k = 0;
  cb->ktree = 0;
  cb->kbranch = 0;

  cb->sons = 0;
  cb->son = 0;

  cb->refs = 0;
  cb->rlist = 0;
  cb->clist = 0;

  cb->Z = NULL;

#ifdef USE_OPENMP
#pragma omp atomic
#endif
  active_clusterbasis++;
}

pclusterbasis
init_clusterbasis(pclusterbasis cb, pccluster t)
{
  uint      i, sons;

  (void) init_leaf_clusterbasis(cb, t);

  sons = t->sons;
  if (sons > 0) {
    cb->sons = sons;
    cb->son =
      (pclusterbasis *) allocmem((size_t) sizeof(pclusterbasis) * sons);
    for (i = 0; i < sons; i++)
      cb->son[i] = NULL;
  }

  return cb;
}

pclusterbasis
init_leaf_clusterbasis(pclusterbasis cb, pccluster t)
{
  assert(cb != NULL);

  init_raw_clusterbasis(cb, t);

  init_amatrix(&cb->V, 0, 0);
  init_amatrix(&cb->E, 0, 0);

  return cb;
}

pclusterbasis
init_sub_clusterbasis(pclusterbasis cb, pclusterbasis src,
		      pccluster t, uint off)
{
  assert(cb != 0);

  init_raw_clusterbasis(cb, t);

  init_sub_amatrix(&cb->V, &src->V, t->size, off, src->k, 0);
  init_amatrix(&cb->E, 0, 0);

  cb->k = cb->ktree = cb->kbranch = src->k;

  return cb;
}

void
uninit_clusterbasis(pclusterbasis cb)
{
  uint      i;

  assert(cb->refs == 0);
  assert(cb->rlist == 0);
  assert(cb->clist == 0);

  if (cb->sons > 0) {
    for (i = 0; i < cb->sons; i++)
      unref_clusterbasis(cb->son[i]);
    freemem(cb->son);
  }

  uninit_amatrix(&cb->V);
  uninit_amatrix(&cb->E);

  assert(active_clusterbasis > 0);

#ifdef USE_OPENMP
#pragma omp atomic
#endif
  active_clusterbasis--;
}

pclusterbasis
new_clusterbasis(pccluster t)
{
  pclusterbasis cb;

  cb = allocmem(sizeof(clusterbasis));

  init_clusterbasis(cb, t);

  return cb;
}

pclusterbasis
new_leaf_clusterbasis(pccluster t)
{
  pclusterbasis cb;

  cb = allocmem(sizeof(clusterbasis));

  init_leaf_clusterbasis(cb, t);

  return cb;
}

void
del_clusterbasis(pclusterbasis cb)
{
  uninit_clusterbasis(cb);

  freemem(cb);
}

/* ------------------------------------------------------------
 * Reference counting
 * ------------------------------------------------------------ */

void
ref_clusterbasis(pclusterbasis * ptr, pclusterbasis cb)
{
  if (*ptr)
    unref_clusterbasis(*ptr);

  *ptr = cb;

  if (cb)
    cb->refs++;
}

void
unref_clusterbasis(pclusterbasis cb)
{
  assert(cb->refs > 0);

  cb->refs--;

  if (cb->refs == 0)
    del_clusterbasis(cb);
}

/* ------------------------------------------------------------
 * Low-level management
 * ------------------------------------------------------------ */

void
update_clusterbasis(pclusterbasis cb)
{
  uint      stree, sbranch;
  uint      i;

  stree = sbranch = 0;

  if (cb->sons == 0) {
    stree += cb->t->size;
    sbranch = cb->t->size;
  }
  else {
    for (i = 0; i < cb->sons; i++) {
      stree += cb->son[i]->ktree;
      if (cb->son[i]->kbranch > sbranch)
	sbranch = cb->son[i]->kbranch;
    }
  }

  cb->ktree = stree + cb->k;
  cb->kbranch = sbranch + cb->k;
}

void
update_tree_clusterbasis(pclusterbasis cb)
{
  uint      i;

  for (i = 0; i < cb->sons; i++) {
    update_tree_clusterbasis(cb->son[i]);
  }
  update_clusterbasis(cb);
}

void
resize_clusterbasis(pclusterbasis cb, uint k)
{
  setrank_clusterbasis(cb, k);
}

void
setrank_clusterbasis(pclusterbasis cb, uint k)
{
  uint      i;

  if (cb->sons > 0) {
    for (i = 0; i < cb->sons; i++)
      resize_amatrix(&cb->son[i]->E, cb->son[i]->k, k);
  }
  else
    resize_amatrix(&cb->V, cb->t->size, k);

  resize_amatrix(&cb->E, k, cb->E.cols);

  cb->k = k;

  update_clusterbasis(cb);
}

/* ------------------------------------------------------------
 * Build clusterbasis based on cluster
 * ------------------------------------------------------------ */

pclusterbasis
build_from_cluster_clusterbasis(pccluster t)
{
  pclusterbasis cb, cb1;
  uint      i;

  cb = new_clusterbasis(t);

  if (cb->sons > 0) {
    for (i = 0; i < cb->sons; i++) {
      cb1 = build_from_cluster_clusterbasis(t->son[i]);
      ref_clusterbasis(cb->son + i, cb1);
    }
  }

  resize_clusterbasis(cb, 0);

  return cb;
}

/* ------------------------------------------------------------
 * Clone a cluster basis
 * ------------------------------------------------------------ */

pclusterbasis
clone_clusterbasis(pcclusterbasis cb)
{
  pclusterbasis cbnew, cbnew1;
  uint      i;

  cbnew = 0;

  if (cb->sons > 0) {
    cbnew = new_clusterbasis(cb->t);
    assert(cbnew->sons == cb->sons);

    for (i = 0; i < cb->sons; i++) {
      cbnew1 = clone_clusterbasis(cb->son[i]);

      ref_clusterbasis(cbnew->son + i, cbnew1);
    }
  }
  else
    cbnew = new_leaf_clusterbasis(cb->t);

  resize_clusterbasis(cbnew, cb->k);

  if (cb->sons > 0) {
    for (i = 0; i < cb->sons; i++)
      copy_amatrix(false, &cb->son[i]->E, &cbnew->son[i]->E);
    if (cb->Z != NULL) {
      copy_amatrix(false, cb->son[i]->Z, cbnew->son[i]->Z);
    }
  }
  else
    copy_amatrix(false, &cb->V, &cbnew->V);

  return cbnew;
}

pclusterbasis
clonestructure_clusterbasis(pcclusterbasis cb)
{
  pclusterbasis cbnew, cbnew1;
  uint      i;

  cbnew = 0;

  if (cb->sons > 0) {
    cbnew = new_clusterbasis(cb->t);
    assert(cbnew->sons == cb->sons);

    for (i = 0; i < cb->sons; i++) {
      cbnew1 = clonestructure_clusterbasis(cb->son[i]);
      ref_clusterbasis(cbnew->son + i, cbnew1);
    }

    update_clusterbasis(cbnew);
  }
  else
    cbnew = new_leaf_clusterbasis(cb->t);

  return cbnew;
}

/* ------------------------------------------------------------
 * Statistics
 * ------------------------------------------------------------ */

uint
getactives_clusterbasis()
{
  return active_clusterbasis;
}

size_t
getsize_clusterbasis(pcclusterbasis cb)
{
  size_t    sz;
  uint      i;

  sz = (size_t) sizeof(clusterbasis);
  sz += getsize_heap_amatrix(&cb->V);
  sz += getsize_heap_amatrix(&cb->E);

  if (cb->sons > 0) {
    sz += (size_t) sizeof(pclusterbasis) * cb->sons;

    for (i = 0; i < cb->sons; i++)
      sz += getsize_clusterbasis(cb->son[i]);
  }

  return sz;
}

/* ------------------------------------------------------------
 * Simple utility functions
 * ------------------------------------------------------------ */

void
clear_weight_clusterbasis(pclusterbasis cb)
{
  uint      i;

  if (cb->son != NULL) {
    for (i = 0; i < cb->sons; i++) {
      clear_weight_clusterbasis(cb->son[i]);
    }
  }
  if (cb->Z != NULL) {
    del_amatrix(cb->Z);
  }
  cb->Z = NULL;
}

/* ------------------------------------------------------------
 * Hierarchical iterator
 * ------------------------------------------------------------ */

void
iterate_clusterbasis(pcclusterbasis cb, uint cbname,
		     void (*pre) (pcclusterbasis cb, uint cbname, void *data),
		     void (*post) (pcclusterbasis cb, uint cbname,
				   void *data), void *data)
{
  uint      cbname1;
  uint      i;

  if (pre)
    pre(cb, cbname, data);

  if (cb->sons > 0) {
    assert(cb->sons == cb->t->sons);

    cbname1 = cbname + 1;
    for (i = 0; i < cb->t->sons; i++) {
      iterate_clusterbasis(cb->son[i], cbname1, pre, post, data);

      cbname1 += cb->t->son[i]->desc;
    }
    assert(cbname1 == cbname + cb->t->desc);
  }

  if (post)
    post(cb, cbname, data);
}

void
iterate_parallel_clusterbasis(pcclusterbasis cb, uint cbname,
			      uint pardepth, void (*pre) (pcclusterbasis cb,
							  uint cbname,
							  void *data),
			      void (*post) (pcclusterbasis cb, uint cbname,
					    void *data), void *data)
{
  uint     *cbname1, cbname2;
#ifdef USE_OPENMP
  uint      nthreads;		/* HACK: Solaris workaround */
#endif
  uint      i;

  if (pre)
    pre(cb, cbname, data);

  if (cb->sons > 0) {
    assert(cb->sons == cb->t->sons);

    cbname1 = (uint *) allocmem((size_t) sizeof(uint) * cb->sons);

    cbname2 = cbname + 1;
    for (i = 0; i < cb->t->sons; i++) {
      cbname1[i] = cbname2;
      cbname2 += cb->son[i]->t->desc;
    }
    assert(cbname2 == cbname + cb->t->desc);

#ifdef USE_OPENMP
    nthreads = cb->sons;
    (void) nthreads;
#pragma omp parallel for if(pardepth > 0), num_threads(nthreads)
#endif
    for (i = 0; i < cb->sons; i++)
      iterate_parallel_clusterbasis(cb->son[i], cbname1[i],
				    (pardepth > 0 ? pardepth - 1 : 0), pre,
				    post, data);

    freemem(cbname1);
  }

  if (post)
    post(cb, cbname, data);
}

/* ------------------------------------------------------------
 * Enumeration
 * ------------------------------------------------------------ */

static void
enumerate(pccluster t, uint tname, pclusterbasis cb, pclusterbasis * cbn)
{
  uint      tname1;
  uint      i;
  cbn[tname] = cb;

  tname1 = tname + 1;

  if (cb == 0 || cb->son == 0) {
    assert(cb == 0 || cb->t == t);
    for (i = 0; i < t->sons; i++) {
      enumerate(t->son[i], tname1, 0, cbn);

      tname1 += t->son[i]->desc;
    }
  }
  else {
    assert(cb->t == t);
    assert(t->sons == cb->sons);

    for (i = 0; i < t->sons; i++) {
      enumerate(t->son[i], tname1, cb->son[i], cbn);

      tname1 += t->son[i]->desc;
    }
  }
  assert(tname1 == tname + t->desc);
}

pclusterbasis *
enumerate_clusterbasis(pccluster t, pclusterbasis cb)
{
  pclusterbasis *cbn;

  cbn = (pclusterbasis *) allocmem((size_t) sizeof(pclusterbasis) * t->desc);

  enumerate(t, 0, cb, cbn);

  return cbn;
}

/* ------------------------------------------------------------
 * Forward and backward transformation
 * ------------------------------------------------------------ */

pavector
new_coeffs_clusterbasis_avector(pcclusterbasis cb)
{
  pavector  xt;

  xt = new_avector(cb->ktree);

  return xt;
}

void
forward_clusterbasis_avector(pcclusterbasis cb, pcavector x, pavector xt)
{
#ifdef USE_OPENMP
  forward_parallel_clusterbasis_avector(cb, x, xt, max_pardepth);
#else
  avector   loc1, loc2;
  pavector  xt1, xc, xp;
  uint      i, xtoff;

  assert(xt->dim == cb->ktree);

  /* This part of xt contains the coefficients for the current cluster */
  xc = init_sub_avector(&loc1, xt, cb->k, 0);
  clear_avector(xc);

  if (cb->sons > 0) {
    xtoff = cb->k;
    for (i = 0; i < cb->sons; i++) {
      /* This part corresponds to the subtree rooted in the i-th son */
      xt1 = init_sub_avector(&loc2, xt, cb->son[i]->ktree, xtoff);

      /* Compute coefficients in the subtree */
      forward_clusterbasis_avector(cb->son[i], x, xt1);

      /* Multiply by transfer matrix */
      mvm_amatrix_avector(1.0, true, &cb->son[i]->E, xt1, xc);

      uninit_avector(xt1);

      xtoff += cb->son[i]->ktree;
    }
    assert(xtoff == cb->ktree);
  }
  else {
    /* Permuted entries of x */
    xp = init_sub_avector(&loc2, xt, cb->t->size, cb->k);

    /* Find and copy entries */
    for (i = 0; i < cb->t->size; i++)
      xp->v[i] = x->v[cb->t->idx[i]];

    /* Multiply by leaf matrix */
    mvm_amatrix_avector(1.0, true, &cb->V, xp, xc);

    uninit_avector(xp);
  }

  uninit_avector(xc);
#endif
}

void
forward_parallel_clusterbasis_avector(pcclusterbasis cb, pcavector x,
				      pavector xt, uint pardepth)
{
  avector   loc1, loc2;
  pavector *xt1, xc, xp;
#ifdef USE_OPENMP
  uint      nthreads;		/* HACK: Solaris workaround */
#endif
  uint      i, xtoff;

  assert(xt->dim == cb->ktree);

  xc = init_sub_avector(&loc1, xt, cb->k, 0);
  clear_avector(xc);

  if (cb->sons > 0) {
    xt1 = (pavector *) allocmem((size_t) sizeof(pavector) * cb->sons);
    xtoff = cb->k;
    for (i = 0; i < cb->sons; i++) {
      xt1[i] = new_sub_avector(xt, cb->son[i]->ktree, xtoff);

      xtoff += cb->son[i]->ktree;
    }
    assert(xtoff == cb->ktree);

#ifdef USE_OPENMP
    nthreads = cb->sons;
    (void) nthreads;
#pragma omp parallel for if(pardepth>0), num_threads(nthreads)
#endif
    for (i = 0; i < cb->sons; i++)
      forward_parallel_clusterbasis_avector(cb->son[i], x, xt1[i],
					    (pardepth >
					     0 ? pardepth - 1 : 0));

    for (i = 0; i < cb->sons; i++) {
      mvm_amatrix_avector(1.0, true, &cb->son[i]->E, xt1[i], xc);

      del_avector(xt1[i]);
    }
    freemem(xt1);
  }
  else {
    xp = init_sub_avector(&loc2, xt, cb->t->size, cb->k);

    for (i = 0; i < cb->t->size; i++)
      setentry_avector(xp, i, getentry_avector(x, cb->t->idx[i]));

    mvm_amatrix_avector(1.0, true, &cb->V, xp, xc);

    uninit_avector(xp);
  }

  uninit_avector(xc);
}

void
forward_nopermutation_clusterbasis_avector(pcclusterbasis cb, pcavector xp,
					   pavector xt)
{
  avector   loc1, loc2, loc3;
  pavector  xp1, xt1, xc;
  uint      i, xpoff, xtoff;

  assert(xp->dim == cb->t->size);
  assert(xt->dim == cb->ktree);

  xc = init_sub_avector(&loc1, xt, cb->k, 0);
  clear_avector(xc);

  if (cb->sons > 0) {
    xpoff = 0;
    xtoff = cb->k;
    for (i = 0; i < cb->sons; i++) {
      xp1 =
	init_sub_avector(&loc2, (pavector) xp, cb->t->son[i]->size, xpoff);
      xt1 = init_sub_avector(&loc3, xt, cb->son[i]->ktree, xtoff);

      forward_nopermutation_clusterbasis_avector(cb->son[i], xp1, xt1);

      uninit_avector(xt1);
      uninit_avector(xp1);

      xt1 = init_sub_avector(&loc3, xt, cb->son[i]->k, xtoff);
      addevaltrans_amatrix_avector(1.0, &cb->son[i]->E, xt1, xc);
      uninit_avector(xt1);

      xpoff += cb->t->son[i]->size;
      xtoff += cb->son[i]->ktree;
    }
    assert(xpoff == cb->t->size);
    assert(xtoff == cb->ktree);
  }
  else {
    xt1 = init_sub_avector(&loc2, xt, cb->t->size, cb->k);

    copy_avector(xp, xt1);

    addevaltrans_amatrix_avector(1.0, &cb->V, xt1, xc);

    uninit_avector(xt1);
  }

  uninit_avector(xc);
}

static void
forward_notransfer(pcclusterbasis cb, pcavector x, pavector xt, pavector xp)
{
  avector   loc1, loc2, loc3;
  pavector  xt1, xp1, xc;
  uint      i, xtoff, xpoff;

  assert(xp->dim == cb->t->size);
  assert(xt->dim == cb->ktree);

  /* This part of xt contains the coefficients for the current cluster */
  xc = init_sub_avector(&loc1, xt, cb->k, 0);
  clear_avector(xc);

  if (cb->sons > 0) {
    xtoff = cb->k;
    xpoff = 0;
    for (i = 0; i < cb->sons; i++) {
      /* This part corresponds to the subtree rooted in the i-th son */
      xt1 = init_sub_avector(&loc2, xt, cb->son[i]->ktree, xtoff);
      xp1 = init_sub_avector(&loc3, xp, cb->son[i]->t->size, xpoff);

      /* Compute coefficients in the subtree */
      forward_notransfer(cb->son[i], x, xt1, xp1);

      /* Release subvectors */
      uninit_avector(xp1);
      uninit_avector(xt1);

      /* Switch to next subvector */
      xtoff += cb->son[i]->ktree;
      xpoff += cb->son[i]->t->size;
    }
    assert(xtoff == cb->ktree);
    assert(xpoff == cb->t->size);
  }
  else {
    xt1 = init_sub_avector(&loc2, xt, cb->t->size, cb->k);
    /* Find and copy entries */
    for (i = 0; i < cb->t->size; i++) {
      assert(cb->t->idx[i] < x->dim);
      xp->v[i] = x->v[cb->t->idx[i]];
      xt1->v[i] = xp->v[i];
    }
  }

  /* Multiply by leaf matrix */
  mvm_amatrix_avector(1.0, true, &cb->V, xp, xc);

  uninit_avector(xc);
}

void
forward_notransfer_clusterbasis_avector(pcclusterbasis cb, pcavector x,
					pavector xt)
{
  pavector  xp;
  avector   tmp;

  xp = init_avector(&tmp, cb->t->size);

  forward_notransfer(cb, x, xt, xp);

  uninit_avector(xp);
}

void
backward_clusterbasis_avector(pcclusterbasis cb, pavector yt, pavector y)
{
#ifdef USE_OPENMP
  backward_parallel_clusterbasis_avector(cb, yt, y, max_pardepth);
#else
  avector   loc1, loc2;
  pavector  yt1, yc, yp;
  uint      i, ytoff;

  assert(yt->dim == cb->ktree);

  /* This part of yt contains the coefficients for the current cluster */
  yc = init_sub_avector(&loc1, yt, cb->k, 0);

  if (cb->sons > 0) {
    ytoff = cb->k;
    for (i = 0; i < cb->sons; i++) {
      /* This part corresponds to the subtree rooted in the i-th son */
      yt1 = init_sub_avector(&loc2, yt, cb->son[i]->ktree, ytoff);

      /* Multiply by transfer matrix */
      mvm_amatrix_avector(1.0, false, &cb->son[i]->E, yc, yt1);

      /* Treat coefficients in the subtree */
      backward_clusterbasis_avector(cb->son[i], yt1, y);

      uninit_avector(yt1);

      ytoff += cb->son[i]->ktree;
    }
    assert(ytoff == cb->ktree);
  }
  else {
    /* Permuted entries of x */
    yp = init_sub_avector(&loc2, yt, cb->t->size, cb->k);

    /* Multiply by leaf matrix */
    mvm_amatrix_avector(1.0, false, &cb->V, yc, yp);

    /* Find and copy entries */
    for (i = 0; i < cb->t->size; i++)
      y->v[cb->t->idx[i]] += yp->v[i];

    uninit_avector(yp);
  }

  uninit_avector(yc);
#endif
}

void
backward_parallel_clusterbasis_avector(pcclusterbasis cb, pavector yt,
				       pavector y, uint pardepth)
{
  avector   loc1, loc2;
  pavector *yt1, yc, yp;
#ifdef USE_OPENMP
  uint      nthreads;		/* HACK: Solaris workaround */
#endif
  uint      i, ytoff;

  assert(yt->dim == cb->ktree);

  yc = init_sub_avector(&loc1, yt, cb->k, 0);

  if (cb->sons > 0) {
    yt1 = (pavector *) allocmem((size_t) sizeof(pavector) * cb->sons);
    ytoff = cb->k;
    for (i = 0; i < cb->sons; i++) {
      yt1[i] = new_sub_avector(yt, cb->son[i]->ktree, ytoff);

      ytoff += cb->son[i]->ktree;
    }
    assert(ytoff == cb->ktree);

#ifdef USE_OPENMP
    nthreads = cb->sons;
    (void) nthreads;
#pragma omp parallel for if(pardepth>0), num_threads(nthreads)
#endif
    for (i = 0; i < cb->sons; i++) {
      mvm_amatrix_avector(1.0, false, &cb->son[i]->E, yc, yt1[i]);

      backward_parallel_clusterbasis_avector(cb->son[i], yt1[i], y,
					     (pardepth >
					      0 ? pardepth - 1 : 0));

      del_avector(yt1[i]);
    }
    freemem(yt1);
  }
  else {
    yp = init_sub_avector(&loc2, yt, cb->t->size, cb->k);

    mvm_amatrix_avector(1.0, false, &cb->V, yc, yp);

    for (i = 0; i < cb->t->size; i++)
      addentry_avector(y, cb->t->idx[i], getentry_avector(yp, i));

    uninit_avector(yp);
  }

  uninit_avector(yc);
}

void
backward_nopermutation_clusterbasis_avector(pcclusterbasis cb, pavector yt,
					    pavector yp)
{
  avector   loc1, loc2, loc3;
  pavector  yt1, yp1, yc;
  uint      i, ytoff, ypoff;

  assert(yp->dim == cb->t->size);
  assert(yt->dim == cb->ktree);

  yc = init_sub_avector(&loc1, yt, cb->k, 0);

  if (cb->sons > 0) {
    ypoff = 0;
    ytoff = cb->k;
    for (i = 0; i < cb->sons; i++) {
      yt1 = init_sub_avector(&loc3, yt, cb->son[i]->k, ytoff);
      addeval_amatrix_avector(1.0, &cb->son[i]->E, yc, yt1);
      uninit_avector(yt1);

      yp1 = init_sub_avector(&loc2, yp, cb->t->son[i]->size, ypoff);
      yt1 = init_sub_avector(&loc3, yt, cb->son[i]->ktree, ytoff);

      backward_nopermutation_clusterbasis_avector(cb->son[i], yt1, yp1);

      uninit_avector(yt1);
      uninit_avector(yp1);

      ypoff += cb->t->son[i]->size;
      ytoff += cb->son[i]->ktree;
    }
    assert(ypoff == cb->t->size);
    assert(ytoff == cb->ktree);
  }
  else {
    yt1 = init_sub_avector(&loc2, yt, cb->t->size, cb->k);

    addeval_amatrix_avector(1.0, &cb->V, yc, yt1);

    add_avector(1.0, yt1, yp);

    uninit_avector(yt1);
  }

  uninit_avector(yc);
}

static void
backward_notransfer(pcclusterbasis cb, pavector yt, pavector y, pavector yp)
{
  avector   loc1, loc2, loc3;
  pavector  yt1, yp1, yc;
  uint      i, ytoff, ypoff;

  assert(yp->dim == cb->t->size);
  assert(yt->dim == cb->ktree);

  /* This part of xt contains the coefficients for the current cluster */
  yc = init_sub_avector(&loc1, yt, cb->k, 0);

  /* Multiply by leaf matrix */
  mvm_amatrix_avector(1.0, false, &cb->V, yc, yp);

  if (cb->sons > 0) {
    ytoff = cb->k;
    ypoff = 0;
    for (i = 0; i < cb->sons; i++) {
      /* This part corresponds to the subtree rooted in the i-th son */
      yt1 = init_sub_avector(&loc2, yt, cb->son[i]->ktree, ytoff);
      yp1 = init_sub_avector(&loc3, yp, cb->son[i]->t->size, ypoff);

      /* Compute coefficients in the subtree */
      backward_notransfer(cb->son[i], yt1, y, yp1);

      /* Release subvectors */
      uninit_avector(yp1);
      uninit_avector(yt1);

      /* Switch to next subvector */
      ytoff += cb->son[i]->ktree;
      ypoff += cb->son[i]->t->size;
    }
    assert(ytoff == cb->ktree);
    assert(ypoff == cb->t->size);
  }
  else {
    yt1 = init_sub_avector(&loc2, yt, cb->t->size, cb->k);
    /* Add entries to result */
    for (i = 0; i < cb->t->size; i++) {
      assert(cb->t->idx[i] < y->dim);
      y->v[cb->t->idx[i]] += yp->v[i] + yt1->v[i];
    }
  }

  uninit_avector(yc);
}

void
backward_notransfer_clusterbasis_avector(pcclusterbasis cb, pavector yt,
					 pavector y)
{
  pavector  yp;
  avector   tmp;

  yp = init_avector(&tmp, cb->t->size);
  clear_avector(yp);

  backward_notransfer(cb, yt, y, yp);

  uninit_avector(yp);
}

/* ------------------------------------------------------------
 * Forward and backward transformation for the root only
 * ------------------------------------------------------------ */

void
compress_clusterbasis_avector(pcclusterbasis cb, pcavector xp, pavector xt)
{
  avector   tmp1, tmp2, tmp3;
  pavector  xt1, xp1, xc;
  uint      i, off;

  assert(xt->dim == cb->kbranch);

  /* This part of xt contains the coefficients for the current cluster */
  xc = init_sub_avector(&tmp1, xt, cb->k, 0);
  clear_avector(xc);

  if (cb->sons > 0) {
    off = 0;
    for (i = 0; i < cb->sons; i++) {
      /* These parts correspond to the subtree rooted in the i-th son */
      xt1 = init_sub_avector(&tmp2, xt, cb->son[i]->kbranch, cb->k);
      xp1 = init_sub_avector(&tmp3, (pavector) xp, cb->son[i]->t->size, off);

      /* Compute coefficients in the subtree */
      compress_clusterbasis_avector(cb->son[i], xp1, xt1);

      /* Multiply by transfer matrix */
      mvm_amatrix_avector(1.0, true, &cb->son[i]->E, xt1, xc);

      uninit_avector(xp1);
      uninit_avector(xt1);

      off += cb->son[i]->t->size;
    }
    assert(off == cb->t->size);
  }
  else {
    /* Multiply by leaf matrix */
    mvm_amatrix_avector(1.0, true, &cb->V, xp, xc);
  }
  uninit_avector(xc);
}

void
expand_clusterbasis_avector(pcclusterbasis cb, pavector yt, pavector yp)
{
  avector   tmp1, tmp2, tmp3;
  pavector  yt1, yp1, yc;
  uint      i, off;

  assert(yt->dim == cb->kbranch);

  /* This part of yt contains the coefficients for the current cluster */
  yc = init_sub_avector(&tmp1, yt, cb->k, 0);

  if (cb->sons > 0) {
    off = 0;
    for (i = 0; i < cb->sons; i++) {
      /* This part corresponds to the i-th son */
      yt1 = init_sub_avector(&tmp2, yt, cb->son[i]->k, cb->k);
      clear_avector(yt1);

      /* Multiply by transfer matrix */
      mvm_amatrix_avector(1.0, false, &cb->son[i]->E, yc, yt1);
      uninit_avector(yt1);

      /* These parts correspond to the subtree rooted in the i-th son */
      yt1 = init_sub_avector(&tmp2, yt, cb->son[i]->kbranch, cb->k);
      yp1 = init_sub_avector(&tmp3, yp, cb->son[i]->t->size, off);

      /* Treat coefficients in the subtree */
      expand_clusterbasis_avector(cb->son[i], yt1, yp1);

      uninit_avector(yp1);
      uninit_avector(yt1);

      off += cb->son[i]->t->size;
    }
    assert(off == cb->t->size);
  }
  else {
    /* Multiply by leaf matrix */
    mvm_amatrix_avector(1.0, false, &cb->V, yc, yp);
  }
  uninit_avector(yc);
}

void
compress_clusterbasis_amatrix(pcclusterbasis cb, pcamatrix Xp, pamatrix Xt)
{
  amatrix   loc1, loc2, loc3;
  pamatrix  Xp1, Xt1, Xc;
  uint      i, xoff;

  assert(Xp->rows >= cb->t->size);
  assert(Xt->rows >= cb->kbranch);

  Xc = init_sub_amatrix(&loc1, Xt, cb->k, 0, Xt->cols, 0);
  clear_amatrix(Xc);

  if (cb->sons > 0) {
    xoff = 0;
    for (i = 0; i < cb->sons; i++) {
      Xp1 = init_sub_amatrix(&loc2, (pamatrix) Xp, cb->t->son[i]->size, xoff,
			     Xp->cols, 0);
      Xt1 = init_sub_amatrix(&loc3, Xt, cb->son[i]->kbranch, cb->k, Xt->cols,
			     0);

      compress_clusterbasis_amatrix(cb->son[i], Xp1, Xt1);
      uninit_amatrix(Xt1);
      uninit_amatrix(Xp1);

      Xt1 = init_sub_amatrix(&loc2, Xt, cb->son[i]->k, cb->k, Xt->cols, 0);

      addmul_amatrix(1.0, true, &cb->son[i]->E, false, Xt1, Xc);
      uninit_amatrix(Xt1);

      xoff += cb->t->son[i]->size;
    }
    assert(xoff == cb->t->size);
  }
  else
    addmul_amatrix(1.0, true, &cb->V, false, Xp, Xc);

  uninit_amatrix(Xc);
}

void
compress_parallel_clusterbasis_amatrix(pcclusterbasis cb, pcamatrix Xp,
				       pamatrix Xt, uint pardepth)
{
  amatrix   loc1;
  pamatrix *Xp1, *Xt1, *Xr1, Xc;
#ifdef USE_OPENMP
  uint      nthreads;		/* HACK: Solaris workaround */
#endif
  uint      i, xoff;

  assert(Xp->rows >= cb->t->size);
  assert(Xt->rows >= cb->kbranch);

  Xc = init_sub_amatrix(&loc1, Xt, cb->k, 0, Xt->cols, 0);
  clear_amatrix(Xc);

  if (cb->sons > 0) {
    Xp1 = (pamatrix *) allocmem((size_t) sizeof(pamatrix) * cb->sons);
    Xt1 = (pamatrix *) allocmem((size_t) sizeof(pamatrix) * cb->sons);
    Xr1 = (pamatrix *) allocmem((size_t) sizeof(pamatrix) * cb->sons);
    xoff = 0;
    for (i = 0; i < cb->sons; i++) {
      Xp1[i] = new_sub_amatrix((pamatrix) Xp, cb->t->son[i]->size, xoff,
			       Xp->cols, 0);
      Xt1[i] = new_sub_amatrix(Xt, cb->son[i]->kbranch, cb->k, Xt->cols, 0);
      Xr1[i] = new_sub_amatrix(Xt, cb->son[i]->k, cb->k, Xt->cols, 0);

      xoff += cb->t->son[i]->size;
    }
    assert(xoff == cb->t->size);

#ifdef USE_OPENMP
    nthreads = cb->sons;
    (void) nthreads;
#pragma omp parallel for if(pardepth>0), num_threads(nthreads)
#endif
    for (i = 0; i < cb->sons; i++)
      compress_parallel_clusterbasis_amatrix(cb->son[i], Xp1[i], Xt1[i],
					     (pardepth >
					      0 ? pardepth - 1 : 0));

    for (i = 0; i < cb->sons; i++) {
      addmul_amatrix(1.0, true, &cb->son[i]->E, false, Xr1[i], Xc);

      del_amatrix(Xr1[i]);
      del_amatrix(Xt1[i]);
      del_amatrix(Xp1[i]);
    }
    freemem(Xr1);
    freemem(Xt1);
    freemem(Xp1);
  }
  else
    addmul_amatrix(1.0, true, &cb->V, false, Xp, Xc);

  uninit_amatrix(Xc);
}

void
expand_clusterbasis_amatrix(pcclusterbasis cb, pamatrix Yt, pamatrix Yp)
{
  amatrix   tmp1, tmp2, tmp3;
  pamatrix  Yt1, Yp1, Yc;
  uint      i, off;

  assert(Yt->rows == cb->kbranch);

  /* This part of Yt contains the coefficients for the current cluster */
  Yc = init_sub_amatrix(&tmp1, Yt, cb->k, 0, Yt->cols, 0);

  if (cb->sons > 0) {
    off = 0;
    for (i = 0; i < cb->sons; i++) {
      /* This part corresponds to the i-th son */
      Yt1 = init_sub_amatrix(&tmp2, Yt, cb->son[i]->k, cb->k, Yt->cols, 0);
      clear_amatrix(Yt1);

      addmul_amatrix(1.0, false, &cb->son[i]->E, false, Yc, Yt1);
      uninit_amatrix(Yt1);

      /* These parts correspond to the subtree rooted in the i-th son */
      Yt1 = init_sub_amatrix(&tmp2, Yt, cb->son[i]->kbranch, cb->k, Yt->cols, 0);
      Yp1 = init_sub_amatrix(&tmp3, Yp, cb->son[i]->t->size, off, Yt->cols, 0);

      /* Treat coefficients in the subtree */
      expand_clusterbasis_amatrix(cb->son[i], Yt1, Yp1);

      uninit_amatrix(Yp1);
      uninit_amatrix(Yt1);

      off += cb->son[i]->t->size;
    }
    assert(off == cb->t->size);
  }
  else {
    /* Multiply by leaf matrix */
    addmul_amatrix(1.0, false, &cb->V, false, Yc, Yp);
  }
  uninit_amatrix(Yc);
}

/* ------------------------------------------------------------
 * Forward and backward transformation for matrices
 * ------------------------------------------------------------ */

void
forward_clusterbasis_amatrix(pcclusterbasis cb, pcamatrix Xp, pamatrix Xt)
{
  amatrix   loc1, loc2, loc3;
  pamatrix  Xp1, Xt1, Xc;
  uint      i, xpoff, xtoff;

  assert(Xp->rows == cb->t->size);
  assert(Xt->rows == cb->ktree);
  assert(Xp->cols == Xt->cols);

  Xc = init_sub_amatrix(&loc1, Xt, cb->k, 0, Xt->cols, 0);
  clear_amatrix(Xc);

  if (cb->sons > 0) {
    xpoff = 0;
    xtoff = cb->k;
    for (i = 0; i < cb->sons; i++) {
      Xp1 = init_sub_amatrix(&loc2, (pamatrix) Xp, cb->t->son[i]->size, xpoff,
			     Xp->cols, 0);
      Xt1 =
	init_sub_amatrix(&loc3, Xt, cb->son[i]->ktree, xtoff, Xt->cols, 0);

      forward_clusterbasis_amatrix(cb->son[i], Xp1, Xt1);

      uninit_amatrix(Xt1);
      uninit_amatrix(Xp1);

      Xt1 = init_sub_amatrix(&loc3, Xt, cb->son[i]->k, xtoff, Xt->cols, 0);

      addmul_amatrix(1.0, true, &cb->son[i]->E, false, Xt1, Xc);
      uninit_amatrix(Xt1);

      xpoff += cb->t->son[i]->size;
      xtoff += cb->son[i]->ktree;
    }
    assert(xpoff == cb->t->size);
    assert(xtoff == cb->ktree);
  }
  else {
    Xt1 = init_sub_amatrix(&loc2, Xt, cb->t->size, cb->k, Xt->cols, 0);

    copy_amatrix(false, Xp, Xt1);

    addmul_amatrix(1.0, true, &cb->V, false, Xt1, Xc);

    uninit_amatrix(Xt1);
  }

  uninit_amatrix(Xc);
}

void
forward_clusterbasis_trans_amatrix(pcclusterbasis cb, pcamatrix Xp,
				   pamatrix Xt)
{
  amatrix   loc1, loc2, loc3;
  pamatrix  Xp1, Xt1, Xc;
  uint      i, xpoff, xtoff;

  assert(Xp->cols == cb->t->size);
  assert(Xt->rows == cb->ktree);
  assert(Xp->rows == Xt->cols);

  Xc = init_sub_amatrix(&loc1, Xt, cb->k, 0, Xt->cols, 0);
  clear_amatrix(Xc);

  if (cb->sons > 0) {
    xpoff = 0;
    xtoff = cb->k;
    for (i = 0; i < cb->sons; i++) {
      Xp1 = init_sub_amatrix(&loc2, (pamatrix) Xp, Xp->rows, 0,
			     cb->t->son[i]->size, xpoff);
      Xt1 =
	init_sub_amatrix(&loc3, Xt, cb->son[i]->ktree, xtoff, Xt->cols, 0);

      forward_clusterbasis_trans_amatrix(cb->son[i], Xp1, Xt1);

      uninit_amatrix(Xt1);
      uninit_amatrix(Xp1);

      Xt1 = init_sub_amatrix(&loc2, Xt, cb->son[i]->k, xtoff, Xt->cols, 0);

      addmul_amatrix(1.0, true, &cb->son[i]->E, false, Xt1, Xc);
      uninit_amatrix(Xt1);

      xpoff += cb->t->son[i]->size;
      xtoff += cb->son[i]->ktree;
    }
    assert(xpoff == cb->t->size);
    assert(xtoff == cb->ktree);
  }
  else {
    Xt1 = init_sub_amatrix(&loc2, Xt, cb->t->size, cb->k, Xt->cols, 0);

    copy_amatrix(true, Xp, Xt1);

    addmul_amatrix(1.0, true, &cb->V, false, Xt1, Xc);

    uninit_amatrix(Xt1);
  }

  uninit_amatrix(Xc);
}

void
backward_clusterbasis_amatrix(pcclusterbasis cb, pamatrix Yt, pamatrix Yp)
{
  amatrix   loc1, loc2, loc3;
  pamatrix  Yt1, Yp1, Yc;
  uint      i, ypoff, ytoff;

  assert(Yp->rows == cb->t->size);
  assert(Yt->rows == cb->ktree);
  assert(Yp->cols == Yt->cols);

  Yc = init_sub_amatrix(&loc1, Yt, cb->k, 0, Yt->cols, 0);

  if (cb->sons > 0) {
    ypoff = 0;
    ytoff = cb->k;

    for (i = 0; i < cb->sons; i++) {
      Yt1 = init_sub_amatrix(&loc2, Yt, cb->son[i]->k, ytoff, Yt->cols, 0);
      addmul_amatrix(1.0, false, &cb->son[i]->E, false, Yc, Yt1);
      uninit_amatrix(Yt1);

      Yp1 =
	init_sub_amatrix(&loc2, Yp, cb->t->son[i]->size, ypoff, Yp->cols, 0);
      Yt1 =
	init_sub_amatrix(&loc3, Yt, cb->son[i]->ktree, ytoff, Yt->cols, 0);

      backward_clusterbasis_amatrix(cb->son[i], Yt1, Yp1);

      uninit_amatrix(Yt1);
      uninit_amatrix(Yp1);

      ypoff += cb->t->son[i]->size;
      ytoff += cb->son[i]->ktree;
    }
    assert(ypoff == cb->t->size);
    assert(ytoff == cb->ktree);
  }
  else {
    Yt1 = init_sub_amatrix(&loc2, Yt, cb->t->size, cb->k, Yt->cols, 0);

    addmul_amatrix(1.0, false, &cb->V, false, Yc, Yt1);

    add_amatrix(1.0, false, Yt1, Yp);

    uninit_amatrix(Yt1);
  }

  uninit_amatrix(Yc);
}

void
backward_clusterbasis_trans_amatrix(pcclusterbasis cb, pamatrix Yt,
				    pamatrix Yp)
{
  amatrix   loc1, loc2, loc3;
  pamatrix  Yt1, Yp1, Yc;
  uint      i, ypoff, ytoff;

  assert(Yp->cols == cb->t->size);
  assert(Yt->rows == cb->ktree);
  assert(Yp->rows == Yt->cols);

  Yc = init_sub_amatrix(&loc1, Yt, cb->k, 0, Yt->cols, 0);

  if (cb->sons > 0) {
    ypoff = 0;
    ytoff = cb->k;

    for (i = 0; i < cb->sons; i++) {
      Yt1 = init_sub_amatrix(&loc2, Yt, cb->son[i]->k, ytoff, Yt->cols, 0);
      addmul_amatrix(1.0, false, &cb->son[i]->E, false, Yc, Yt1);
      uninit_amatrix(Yt1);

      Yp1 =
	init_sub_amatrix(&loc2, Yp, Yp->rows, 0, cb->t->son[i]->size, ypoff);
      Yt1 =
	init_sub_amatrix(&loc3, Yt, cb->son[i]->ktree, ytoff, Yt->cols, 0);

      backward_clusterbasis_amatrix(cb->son[i], Yt1, Yp1);

      uninit_amatrix(Yt1);
      uninit_amatrix(Yp1);

      ypoff += cb->t->son[i]->size;
      ytoff += cb->son[i]->ktree;
    }
    assert(ypoff == cb->t->size);
    assert(ytoff == cb->ktree);
  }
  else {
    Yt1 = init_sub_amatrix(&loc2, Yt, cb->t->size, cb->k, Yt->cols, 0);

    addmul_amatrix(1.0, false, &cb->V, false, Yc, Yt1);

    add_amatrix(1.0, true, Yt1, Yp);

    uninit_amatrix(Yt1);
  }

  uninit_amatrix(Yc);
}

/* ------------------------------------------------------------
 * Simple computations
 * ------------------------------------------------------------ */

void
addeval_clusterbasis_avector(field alpha, pcclusterbasis cb,
			     pcavector yc, pavector yp)
{
  pavector  yt, yt1;
  avector   tmp1, tmp2;

  assert(yc->dim == cb->k);
  assert(yp->dim == cb->t->size);

  yt = init_avector(&tmp1, cb->kbranch);

  yt1 = init_sub_avector(&tmp2, yt, cb->k, 0);
  copy_avector(yc, yt1);
  scale_avector(alpha, yt1);
  uninit_avector(yt1);

  expand_clusterbasis_avector(cb, yt, yp);

  uninit_avector(yt);

  /*
     pavector x1, y1;
     avector tmp1, tmp2;
     uint i, off;

     assert(x->dim == cb->k);
     assert(y->dim == cb->t->size);

     if (cb->son == 0)
     addeval_amatrix_avector(1.0, &cb->V, x, y);
     else {
     assert(cb->sons > 0);
     off = 0;
     for (i = 0; i < cb->sons; i++) {
     x1 = init_avector(&tmp1, cb->son[i]->k);
     clear_avector(x1);
     addeval_amatrix_avector(1.0, &cb->son[i]->E, x, x1);
     y1 = init_sub_avector(&tmp2, y, cb->son[i]->t->size, off);
     addeval_clusterbasis_avector(cb->son[i], x1, y1);
     uninit_avector(x1);
     uninit_avector(y1);
     off += cb->son[i]->t->size;
     }
     assert(off == cb->t->size);
     }
   */
}

void
addevaltrans_clusterbasis_avector(field alpha, pcclusterbasis cb,
				  pcavector xp, pavector xc)
{
  pavector  xt, xt1;
  avector   tmp1, tmp2;

  assert(xp->dim == cb->t->size);
  assert(xc->dim == cb->k);

  xt = init_avector(&tmp1, cb->kbranch);

  compress_clusterbasis_avector(cb, xp, xt);

  xt1 = init_sub_avector(&tmp2, xt, cb->k, 0);
  add_avector(alpha, xt1, xc);
  uninit_avector(xt1);

  uninit_avector(xt);

  /*
     pavector x1, y1;
     avector tmp1, tmp2;
     uint i, off;

     assert(x->dim == cb->t->size);
     assert(y->dim == cb->k);

     if (cb->son == 0)
     addevaltrans_amatrix_avector(1.0, &cb->V, x, y);
     else {
     assert(cb->sons > 0);
     off = 0;
     for (i = 0; i < cb->sons; i++) {
     x1 = init_sub_avector(&tmp1, x, cb->son[i]->t->size, off);
     y1 = init_avector(&tmp2, cb->son[i]->k);
     clear_avector(y1);
     addevaltrans_clusterbasis_avector(cb->son[i], x1, y1);
     addevaltrans_amatrix_avector(1.0, &cb->son[i]->E, y1, y);
     uninit_avector(x1);
     uninit_avector(y1);
     off += cb->son[i]->t->size;
     }
     assert(off == cb->t->size);
     }
   */
}

/* ------------------------------------------------------------
 * Orthogonalization
 * ------------------------------------------------------------ */

void
ortho_clusterbasis(pclusterbasis cb, pclusteroperator old2new)
{
  amatrix   tmp1, tmp2, tmp3;
  avector   tmp4;
  pamatrix  Vhat, Eold, Qhat, Vhat1, Eold1, Qhat1;
  pavector  tau;
  uint      sons = cb->sons;
  uint      i, k, kold, m, off;

  assert(old2new == 0 || sons == old2new->sons);
  assert(cb->rlist == 0);
  assert(cb->clist == 0);

  kold = cb->k;

  if (sons > 0) {
    /* Orthogonalize sons' bases */
    m = 0;
    for (i = 0; i < cb->sons; i++) {
      ortho_clusterbasis(cb->son[i], (old2new ? old2new->son[i] : 0));

      m += cb->son[i]->k;
    }

    /* Prepare Vhat */
    Vhat = init_amatrix(&tmp1, m, kold);

    /* Copy transfer matrices into Vhat */
    off = 0;
    for (i = 0; i < cb->sons; i++) {
      Vhat1 = init_sub_amatrix(&tmp2, Vhat, cb->son[i]->k, off, kold, 0);

      copy_amatrix(false, &cb->son[i]->E, Vhat1);

      uninit_amatrix(Vhat1);

      off += cb->son[i]->k;
    }
    assert(off == m);

    /* Determine new rank */
    k = UINT_MIN(m, kold);

    /* Compute Householder factorization Vhat = Qhat R */
    tau = init_avector(&tmp4, k);
    qrdecomp_amatrix(Vhat, tau);

    /* Clone original transfer matrix */
    Eold = init_amatrix(&tmp2, kold, cb->E.cols);
    copy_amatrix(false, &cb->E, Eold);

    /* Switch to new rank */
    setrank_clusterbasis(cb, k);

    /* Replace transfer matrix E by R E, restoring consistency */
    Vhat1 = init_sub_amatrix(&tmp3, Vhat, k, 0, kold, 0);
    triangulareval_amatrix(false, false, false, Vhat1, false, Eold);
    uninit_amatrix(Vhat1);
    Eold1 = init_sub_amatrix(&tmp3, Eold, k, 0, Eold->cols, 0);
    copy_amatrix(false, Eold1, &cb->E);
    uninit_amatrix(Eold1);
    uninit_amatrix(Eold);

    /* Copy R into "old2new" cluster operator if required */
    if (old2new) {
      resize_clusteroperator(old2new, k, kold);

      copy_upper_amatrix(Vhat, false, &old2new->C);
    }

    /* Expand matrix Qhat of Householder factorization */
    Qhat = init_amatrix(&tmp2, m, k);
    qrexpand_amatrix(Vhat, tau, Qhat);

    /* Blocks of Qhat are the new transfer matrices */
    off = 0;
    for (i = 0; i < cb->sons; i++) {
      Qhat1 = init_sub_amatrix(&tmp3, Qhat, cb->son[i]->k, off, k, 0);

      assert(cb->son[i]->E.rows == cb->son[i]->k);
      assert(cb->son[i]->E.cols == cb->k);
      copy_amatrix(false, Qhat1, &cb->son[i]->E);

      uninit_amatrix(Qhat1);

      off += cb->son[i]->k;
    }
    assert(off == m);

    /* Clean up */
    uninit_amatrix(Qhat);
    uninit_avector(tau);
    uninit_amatrix(Vhat);
  }
  else {
    m = cb->t->size;

    /* Copy leaf matrix into Vhat */
    Vhat = init_amatrix(&tmp1, m, kold);

    copy_amatrix(false, &cb->V, Vhat);

    /* Determine new rank */
    k = UINT_MIN(m, kold);

    /* Compute Householder factorization Vhat = Qhat R */
    tau = init_avector(&tmp4, k);
    qrdecomp_amatrix(Vhat, tau);

    /* Clone original transfer matrix */
    Eold = init_amatrix(&tmp2, kold, cb->E.cols);
    copy_amatrix(false, &cb->E, Eold);

    /* Switch to new rank */
    setrank_clusterbasis(cb, k);

    /* Replace transfer matrix E by R E, restoring consistency */
    Vhat1 = init_sub_amatrix(&tmp3, Vhat, k, 0, kold, 0);
    triangulareval_amatrix(false, false, false, Vhat1, false, Eold);
    uninit_amatrix(Vhat1);
    Eold1 = init_sub_amatrix(&tmp3, Eold, k, 0, Eold->cols, 0);
    copy_amatrix(false, Eold1, &cb->E);
    uninit_amatrix(Eold1);
    uninit_amatrix(Eold);

    /* Copy R into "old2new" cluster operator if required */
    if (old2new) {
      resize_clusteroperator(old2new, k, kold);

      copy_upper_amatrix(Vhat, false, &old2new->C);
    }

    /* Expand matrix Qhat of Householder factorization */
    qrexpand_amatrix(Vhat, tau, &cb->V);

    /* Clean up */
    uninit_avector(tau);
    uninit_amatrix(Vhat);
  }
}

real
check_ortho_clusterbasis(pcclusterbasis cb)
{
  amatrix   tmp;
  pamatrix  G;
  real      norm, norm1;
  uint      i;

  G = init_amatrix(&tmp, cb->k, cb->k);
  clear_amatrix(G);

  norm = 0.0;

  if (cb->sons > 0) {
    for (i = 0; i < cb->sons; i++) {
      norm1 = check_ortho_clusterbasis(cb->son[i]);
      if (norm1 > norm)
	norm = norm1;

      addmul_amatrix(1.0, true, &cb->son[i]->E, false, &cb->son[i]->E, G);
    }
  }
  else
    addmul_amatrix(1.0, true, &cb->V, false, &cb->V, G);

  for (i = 0; i < cb->k; i++)
    addentry_amatrix(G, i, i, -1.0);

  norm1 = normfrob_amatrix(G);

  if (norm1 > norm)
    norm = norm1;

  uninit_amatrix(G);

  return norm;
}

pclusteroperator
weight_clusterbasis_clusteroperator(pcclusterbasis cb, pclusteroperator co)
{
  amatrix   tmp1, tmp2;
  avector   tmp3;
  pamatrix  Vhat, Vhat1;
  pavector  tau;
  uint      i, k, m, off;

  assert(cb->sons == co->sons);

  if (cb->sons > 0) {
    m = 0;
    for (i = 0; i < cb->sons; i++) {
      weight_clusterbasis_clusteroperator(cb->son[i], co->son[i]);
      m += co->son[i]->krow;
    }

    Vhat = init_amatrix(&tmp1, m, cb->k);
    clear_amatrix(Vhat);

    off = 0;
    for (i = 0; i < cb->sons; i++) {
      Vhat1 = init_sub_amatrix(&tmp2, Vhat, co->son[i]->krow, off, cb->k, 0);

      assert(co->son[i]->kcol == cb->son[i]->E.rows);
      addmul_amatrix(1.0, false, &co->son[i]->C, false, &cb->son[i]->E,
		     Vhat1);

      uninit_amatrix(Vhat1);

      off += co->son[i]->krow;
    }
    assert(off == m);

    k = UINT_MIN(m, cb->k);

    tau = init_avector(&tmp3, k);

    qrdecomp_amatrix(Vhat, tau);

    resize_clusteroperator(co, k, cb->k);

    copy_upper_amatrix(Vhat, false, &co->C);

    uninit_avector(tau);
    uninit_amatrix(Vhat);
  }
  else {
    m = cb->t->size;

    Vhat = init_amatrix(&tmp1, m, cb->k);

    copy_amatrix(false, &cb->V, Vhat);

    k = UINT_MIN(m, cb->k);

    tau = init_avector(&tmp3, k);

    qrdecomp_amatrix(Vhat, tau);

    resize_clusteroperator(co, k, cb->k);

    copy_upper_amatrix(Vhat, false, &co->C);

    uninit_avector(tau);
    uninit_amatrix(Vhat);
  }

  return co;
}

struct _computeweight_data {
  pclusterbasis *V;
  pamatrix  R;
};

static void
computeweight(pccluster t, uint tname, void *data)
{
  amatrix   tmp1, tmp2;
  avector   tmp3;
  struct _computeweight_data *wd = (struct _computeweight_data *) data;
  pcclusterbasis cb = wd->V[tname];
  pamatrix  R = wd->R;
  pamatrix  Vhat, Vhat1;
  pavector  tau;
  uint      sons;
  uint      m, k;
  uint      i, tname1;

  /* Quick exit if there is no cluster basis for this cluster */
  if (cb == 0)
    return;

  /* Quick exit if the rank is zero */
  if (cb->k == 0) {
    init_amatrix(R + tname, 0, 0);
    return;
  }

  sons = cb->sons;

  if (sons > 0) {
    assert(t->sons == sons);

    /* Determine ranks of sons */
    m = 0;
    tname1 = tname + 1;
    for (i = 0; i < sons; i++) {
      m += R[tname1].rows;
      tname1 += t->son[i]->desc;
    }
    assert(tname1 == tname + t->desc);

    /* Assemble half-compressed matrix Vhat */
    Vhat = init_amatrix(&tmp1, m, cb->k);
    m = 0;
    tname1 = tname + 1;
    for (i = 0; i < sons; i++) {
      Vhat1 = init_sub_amatrix(&tmp2, Vhat, R[tname1].rows, m, cb->k, 0);

      clear_amatrix(Vhat1);
      addmul_amatrix(1.0, false, R + tname1, false, &cb->son[i]->E, Vhat1);

      uninit_amatrix(Vhat1);

      m += R[tname1].rows;
      tname1 += t->son[i]->desc;
    }
    assert(tname1 == tname + t->desc);
    assert(m == Vhat->rows);
  }
  else {
    m = cb->t->size;

    assert(m == cb->V.rows);

    /* Copy V to Vhat */
    Vhat = init_amatrix(&tmp1, m, cb->k);
    copy_amatrix(false, &cb->V, Vhat);
  }

  /* Find QR decomposition */
  k = UINT_MIN(cb->k, m);
  tau = init_avector(&tmp3, k);
  qrdecomp_amatrix(Vhat, tau);

  /* Copy result */
  init_amatrix(R + tname, k, cb->k);
  copy_upper_amatrix(Vhat, false, R + tname);

  /* Clean up */
  uninit_avector(tau);
  uninit_amatrix(Vhat);
}

pamatrix
weight_enum_clusterbasis_clusteroperator(pcclusterbasis cb)
{
  pamatrix  R;
  struct _computeweight_data wd;

  R = (pamatrix) allocmem(sizeof(amatrix) * cb->t->desc);

  wd.V = enumerate_clusterbasis(cb->t, (pclusterbasis) cb);
  wd.R = R;

  iterate_cluster(cb->t, 0, 0, computeweight, &wd);

  freemem(wd.V);

  return wd.R;
}

/* ------------------------------------------------------------
 * File I/O
 * ------------------------------------------------------------ */

#ifdef USE_NETCDF
static void
write_count(pcclusterbasis cb, size_t * clusters, size_t * coeffs)
{
  uint      sons;
  uint      i;

  /* Increase cluster counter */
  (*clusters)++;

  /* Add number of coefficients of transfer matrix */
  (*coeffs) += cb->E.rows * cb->E.cols;

  if (cb->sons > 0) {
    sons = cb->sons;

    /* Handle sons (and check that their transfer matrices are of
     * the expected size */
    for (i = 0; i < sons; i++) {
      assert(cb->son[i]->E.rows == cb->son[i]->k);
      assert(cb->son[i]->E.cols == cb->k);

      write_count(cb->son[i], clusters, coeffs);
    }
  }
  else
    (*coeffs) += cb->V.rows * cb->V.cols;
}

static void
write_cdf(pcclusterbasis cb,
	  size_t clusters, size_t coeffs,
	  size_t * clusteridx, size_t * coeffidx,
	  int nc_file, int nc_sons, int nc_size, int nc_rank, int nc_coeff)
{
  size_t    start, count;
  ptrdiff_t stride;
  int       val, result;
  uint      i, j;

  assert(*clusteridx <= clusters);

  /* Write number of sons to nc_sons[*clusteridx] */
  start = *clusteridx;
  count = 1;
  stride = 1;
  val = cb->sons;
  result = nc_put_vars(nc_file, nc_sons, &start, &count, &stride, &val);
  assert(result == NC_NOERR);

  /* Write size of cluster to nc_size[*clusteridx] */
  val = cb->t->size;
  result = nc_put_vars(nc_file, nc_size, &start, &count, &stride, &val);
  assert(result == NC_NOERR);

  /* Write rank to nc_rank[*clusteridx] */
  val = cb->k;
  result = nc_put_vars(nc_file, nc_rank, &start, &count, &stride, &val);
  assert(result == NC_NOERR);

  /* Increment cluster counter */
  (*clusteridx)++;

  if (cb->sons > 0) {
    /* Take care of sons */
    for (i = 0; i < cb->sons; i++)
      write_cdf(cb->son[i],
		clusters, coeffs,
		clusteridx, coeffidx,
		nc_file, nc_sons, nc_size, nc_rank, nc_coeff);

    /* Write transfer matrices */
    if(coeffs > 0) {
      start = *coeffidx;
      stride = 1;
      for (i = 0; i < cb->sons; i++) {
	count = cb->son[i]->k;
	assert(start + cb->son[i]->k * cb->k <= coeffs);
	
	for (j = 0; j < cb->k; j++) {
	  result =
	    nc_put_vars(nc_file, nc_coeff, &start, &count, &stride,
			cb->son[i]->E.a + j * cb->son[i]->E.ld);
	  assert(result == NC_NOERR);
	  
	  start += cb->son[i]->k;
	}
      }
      (*coeffidx) = start;
    }
  }
  else {
    /* Write leaf matrix to nc_coeff */
    if(coeffs > 0) {
      start = *coeffidx;
      count = cb->V.rows;
      assert(start + cb->V.rows * cb->V.cols <= coeffs);
      for (j = 0; j < cb->V.cols; j++) {
	result =
	  nc_put_vars(nc_file, nc_coeff, &start, &count, &stride,
		      cb->V.a + j * cb->V.ld);
	assert(result == NC_NOERR);
	
	start += cb->V.rows;
      }
      *coeffidx = start;
    }
  }
}

void
write_cdf_clusterbasis(pcclusterbasis cb, const char *name)
{
  size_t    clusters, clusteridx;
  size_t    coeffs, coeffidx;
  int       nc_file, nc_sons, nc_size, nc_rank, nc_coeff;
  int       nc_clusters, nc_coeffs;
  int       result;

  /* Count number of clusters and coefficients */
  clusters = 0;
  coeffs = 0;
  write_count(cb, &clusters, &coeffs);

  /* Create NetCDF file */
  result = nc_create(name, NC_64BIT_OFFSET, &nc_file);
  assert(result == NC_NOERR);

  /* Define "clusters" dimension */
  result = nc_def_dim(nc_file, "clusters", clusters, &nc_clusters);
  assert(result == NC_NOERR);

  /* Define "coeffs" dimension */
  nc_coeffs = 0;
  if(coeffs > 0) {
    result = nc_def_dim(nc_file, "coeffs", coeffs, &nc_coeffs);
    assert(result == NC_NOERR);
  }

  /* Define "sons" variable */
  result = nc_def_var(nc_file, "sons", NC_INT, 1, &nc_clusters, &nc_sons);
  assert(result == NC_NOERR);

  /* Define "size" variable */
  result = nc_def_var(nc_file, "size", NC_INT, 1, &nc_clusters, &nc_size);
  assert(result == NC_NOERR);

  /* Define "rank" variable */
  result = nc_def_var(nc_file, "rank", NC_INT, 1, &nc_clusters, &nc_rank);
  assert(result == NC_NOERR);

  /* Define "coeff" variable */
  if(coeffs > 0) {
    result = nc_def_var(nc_file, "coeff", NC_DOUBLE, 1, &nc_coeffs, &nc_coeff);
    assert(result == NC_NOERR);
  }

  /* Finish NetCDF define mode */
  result = nc_enddef(nc_file);
  assert(result == NC_NOERR);

  /* Write matrices to NetCDF variables */
  clusteridx = 0;
  coeffidx = 0;
  write_cdf(cb, clusters, coeffs, &clusteridx, &coeffidx,
	    nc_file, nc_sons, nc_size, nc_rank, nc_coeff);
  assert(clusteridx == clusters);
  assert(coeffidx == coeffs);

  /* Close file */
  result = nc_close(nc_file);
  assert(result == NC_NOERR);
}

static void
prefix_name(char *buf, int bufsize, const char *prefix, const char *name)
{
  if (prefix)
    snprintf(buf, bufsize, "%s_%s", prefix, name);
  else
    strncpy(buf, name, bufsize);
}

void
write_cdfpart_clusterbasis(pcclusterbasis cb, int nc_file, const char *prefix)
{
  size_t    clusters, clusteridx;
  size_t    coeffs, coeffidx;
  char     *buf;
  int       bufsize;
  int       nc_sons, nc_size, nc_rank, nc_coeff;
  int       nc_clusters, nc_coeffs;
  int       result;

  /* Prepare buffer for prefixed names */
  bufsize = strlen(prefix) + 16;
  buf = (char *) allocmem(sizeof(char) * bufsize);

  /* Count number of clusters and coefficients */
  clusters = 0;
  coeffs = 0;
  write_count(cb, &clusters, &coeffs);

  /* Switch NetCDF file to define mode */
  result = nc_redef(nc_file);
  assert(result == NC_NOERR || result == NC_EINDEFINE);

  /* Define "clusters" dimension */
  prefix_name(buf, bufsize, prefix, "clusters");
  result = nc_def_dim(nc_file, buf, clusters, &nc_clusters);
  assert(result == NC_NOERR);

  /* Define "coeffs" dimension */
  nc_coeffs = 0;
  if(coeffs > 0) {
    prefix_name(buf, bufsize, prefix, "coeffs");
    result = nc_def_dim(nc_file, buf, coeffs, &nc_coeffs);
    assert(result == NC_NOERR);
  }

  /* Define "sons" variable */
  prefix_name(buf, bufsize, prefix, "sons");
  result = nc_def_var(nc_file, buf, NC_INT, 1, &nc_clusters, &nc_sons);
  assert(result == NC_NOERR);

  /* Define "size" variable */
  prefix_name(buf, bufsize, prefix, "size");
  result = nc_def_var(nc_file, buf, NC_INT, 1, &nc_clusters, &nc_size);
  assert(result == NC_NOERR);

  /* Define "rank" variable */
  prefix_name(buf, bufsize, prefix, "rank");
  result = nc_def_var(nc_file, buf, NC_INT, 1, &nc_clusters, &nc_rank);
  assert(result == NC_NOERR);

  /* Define "coeff" variable */
  if(coeffs > 0) {
    prefix_name(buf, bufsize, prefix, "coeff");
    result = nc_def_var(nc_file, buf, NC_DOUBLE, 1, &nc_coeffs, &nc_coeff);
    assert(result == NC_NOERR);
  }

  /* Finish NetCDF define mode */
  result = nc_enddef(nc_file);
  assert(result == NC_NOERR);

  /* Write matrices to NetCDF variables */
  clusteridx = 0;
  coeffidx = 0;
  write_cdf(cb, clusters, coeffs, &clusteridx, &coeffidx,
	    nc_file, nc_sons, nc_size, nc_rank, nc_coeff);
  assert(clusteridx == clusters);
  assert(coeffidx == coeffs);

  /* Clean up */
  nc_sync(nc_file);
  freemem(buf);
}

static    pclusterbasis
read_cdf_part(int nc_file, size_t clusters, size_t coeffs,
	      int nc_sons, int nc_size, int nc_rank,
	      int nc_coeff,
	      pcluster *t, uint * idx, size_t * clusteridx, size_t * coeffidx)
{
  pclusterbasis cb, cb1;
  uint     *idx1;
  uint      size;
  uint      sons;
  uint      k;
  uint      i, j;
  size_t    start, count;
  ptrdiff_t stride;
  int       val, result;

  /* Get number of sons */
  start = *clusteridx;
  count = 1;
  stride = 1;
  result = nc_get_vars(nc_file, nc_sons, &start, &count, &stride, &val);
  assert(result == NC_NOERR);
  sons = val;

  /* Get size of cluster */
  result = nc_get_vars(nc_file, nc_size, &start, &count, &stride, &val);
  assert(result == NC_NOERR);
  size = val;

  /* Get rank */
  result = nc_get_vars(nc_file, nc_rank, &start, &count, &stride, &val);
  assert(result == NC_NOERR);
  k = val;

  if (sons > 0) {
    /* If we have not inherited a cluster, create one */
    if (*t == 0) {
      /* If we have not inherited an index, create one */
      if (idx == 0) {
	idx = (uint *) allocmem(sizeof(uint) * size);
	for (i = 0; i < size; i++)
	  idx[i] = i;
      }

      *t = new_cluster(size, idx, sons, 1);
    }
    else {
      assert((*t)->size == size);
      assert((*t)->sons == sons);

      idx = (*t)->idx;
    }

    /* Create the new cluster basis */
    cb = new_clusterbasis(*t);

    /* Increment cluster counter */
    (*clusteridx)++;

    /* Handle sons */
    idx1 = idx;
    for (i = 0; i < sons; i++) {
      cb1 = read_cdf_part(nc_file, clusters, coeffs,
			  nc_sons, nc_size, nc_rank,
			  nc_coeff,
			  (*t)->son + i, idx1, clusteridx, coeffidx);
      ref_clusterbasis(cb->son + i, cb1);

      idx1 += cb1->t->size;
    }
    assert(idx1 == idx + cb->t->size);

    /* Set rank */
    resize_clusterbasis(cb, k);

    /* Read transfer matrices */
    if(coeffs > 0) {
      start = (*coeffidx);
      stride = 1;
      for (i = 0; i < sons; i++) {
	assert(cb->son[i]->E.rows == cb->son[i]->k);
	assert(cb->son[i]->E.cols == k);
	
	count = cb->son[i]->k;
	for (j = 0; j < k; j++) {
	  result =
	    nc_get_vars(nc_file, nc_coeff, &start, &count, &stride,
			cb->son[i]->E.a + j * cb->son[i]->E.ld);
	  assert(result == NC_NOERR);
	  
	  start += cb->son[i]->k;
	}
      }
      (*coeffidx) = start;
    }
  }
  else {
    /* If we have not inherited a cluster, create one */
    if (*t == 0) {
      /* If we have not inherited an index, create one */
      if (idx == 0) {
	idx = (uint *) allocmem(sizeof(uint) * size);
	for (i = 0; i < size; i++)
	  idx[i] = i;
      }

      *t = new_cluster(size, idx, 0, 1);
    }
    else {
      assert((*t)->size == size);

      idx = (*t)->idx;
    }

    /* Create the new cluster basis */
    cb = new_leaf_clusterbasis(*t);

    /* Set rank */
    resize_clusterbasis(cb, k);

    /* Increment cluster counter */
    (*clusteridx)++;

    /* Read leaf matrix */
    if(coeffs > 0) {
      start = (*coeffidx);
      count = size;
      stride = 1;
      for (j = 0; j < k; j++) {
	result =
	  nc_get_vars(nc_file, nc_coeff, &start, &count, &stride,
		      cb->V.a + j * cb->V.ld);
	assert(result == NC_NOERR);
	
	start += size;
      }
      (*coeffidx) = start;
    }
  }

  update_clusterbasis(cb);

  return cb;
}

pclusterbasis
read_cdf_clusterbasis(const char *name, pccluster t)
{
  pclusterbasis cb;
  pcluster  t0;
  size_t    clusters, clusteridx;
  size_t    coeffs, coeffidx;
  char      dimname[NC_MAX_NAME + 1];
  int       nc_file, nc_sons, nc_size, nc_rank, nc_coeff;
  int       nc_clusters, nc_coeffs;
  int       result;

  /* Open NetCDF file */
  result = nc_open(name, NC_NOWRITE, &nc_file);
  assert(result == NC_NOERR);

  /* Get "clusters" dimension */
  result = nc_inq_dimid(nc_file, "clusters", &nc_clusters);
  assert(result == NC_NOERR);
  result = nc_inq_dim(nc_file, nc_clusters, dimname, &clusters);
  assert(result == NC_NOERR);

  /* Get "coeffs" dimension */
  result = nc_inq_dimid(nc_file, "coeffs", &nc_coeffs);
  coeffs = 0;
  /* Error code NC_EBADDIM means that "coeffs" is not present in
   * the file, and by our convention this means that it is zero. */
  if(result != NC_EBADDIM) {
    assert(result == NC_NOERR);
    result = nc_inq_dim(nc_file, nc_coeffs, dimname, &coeffs);
    assert(result == NC_NOERR);
  }

  /* Get "sons" variable */
  result = nc_inq_varid(nc_file, "sons", &nc_sons);
  assert(result == NC_NOERR);

  /* Get "size" variable */
  result = nc_inq_varid(nc_file, "size", &nc_size);
  assert(result == NC_NOERR);

  /* Get "rank" variable */
  result = nc_inq_varid(nc_file, "rank", &nc_rank);
  assert(result == NC_NOERR);

  /* Get "coeff" variable */
  result = nc_inq_varid(nc_file, "coeff", &nc_coeff);
  assert(result == NC_NOERR);

  /* Read cluster basis from NetCDF variables */
  t0 = (pcluster) t;
  clusteridx = 0;
  coeffidx = 0;

  cb = read_cdf_part(nc_file, clusters, coeffs,
		     nc_sons, nc_size, nc_rank, nc_coeff,
		     &t0, 0, &clusteridx, &coeffidx);
  assert(clusteridx == clusters);
  assert(coeffidx == coeffs);

  /* Closer NetCDF file */
  nc_close(nc_file);

  return cb;
}

pclusterbasis
read_cdfpart_clusterbasis(int nc_file, const char *prefix, pccluster t)
{
  pclusterbasis cb;
  pcluster  t0;
  size_t    clusters, clusteridx;
  size_t    coeffs, coeffidx;
  char      dimname[NC_MAX_NAME + 1];
  char     *buf;
  int       bufsize;
  int       nc_sons, nc_size, nc_rank, nc_coeff;
  int       nc_clusters, nc_coeffs;
  int       result;

  /* Prepare buffer for prefixed names */
  bufsize = strlen(prefix) + 16;
  buf = (char *) allocmem(sizeof(char) * bufsize);

  /* Get "clusters" dimension */
  prefix_name(buf, bufsize, prefix, "clusters");
  result = nc_inq_dimid(nc_file, buf, &nc_clusters);
  assert(result == NC_NOERR);
  result = nc_inq_dim(nc_file, nc_clusters, dimname, &clusters);
  assert(result == NC_NOERR);

  /* Get "coeffs" dimension */
  prefix_name(buf, bufsize, prefix, "coeffs");
  result = nc_inq_dimid(nc_file, buf, &nc_coeffs);
  coeffs = 0;
  /* Error code NC_EBADDIM means that "coeffs" is not present in
   * the file, and by our convention this means that it is zero. */
  if(result != NC_EBADDIM) {
    assert(result == NC_NOERR);
    result = nc_inq_dim(nc_file, nc_coeffs, dimname, &coeffs);
    assert(result == NC_NOERR);
  }

  /* Get "sons" variable */
  prefix_name(buf, bufsize, prefix, "sons");
  result = nc_inq_varid(nc_file, buf, &nc_sons);
  assert(result == NC_NOERR);

  /* Get "size" variable */
  prefix_name(buf, bufsize, prefix, "size");
  result = nc_inq_varid(nc_file, buf, &nc_size);
  assert(result == NC_NOERR);

  /* Get "rank" variable */
  prefix_name(buf, bufsize, prefix, "rank");
  result = nc_inq_varid(nc_file, buf, &nc_rank);
  assert(result == NC_NOERR);

  /* Get "coeff" variable */
  if(coeffs > 0) {
    prefix_name(buf, bufsize, prefix, "coeff");
    result = nc_inq_varid(nc_file, buf, &nc_coeff);
    assert(result == NC_NOERR);
  }

  /* Read cluster basis from NetCDF variables */
  t0 = (pcluster) t;
  clusteridx = 0;
  coeffidx = 0;

  cb = read_cdf_part(nc_file, clusters, coeffs,
		     nc_sons, nc_size, nc_rank, nc_coeff,
		     &t0, 0, &clusteridx, &coeffidx);
  assert(clusteridx == clusters);
  assert(coeffidx == coeffs);

  /* Clean up */
  freemem(buf);

  return cb;
}
#endif
