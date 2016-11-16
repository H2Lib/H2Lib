
/* ------------------------------------------------------------
 * This is the file "dclusterbasis.c" of the H2Lib package.
 * All rights reserved, Steffen Boerm 2015
 * ------------------------------------------------------------ */

#include <stdio.h>

#include "factorizations.h"
#include "dclusterbasis.h"

static uint active_dclusterbasis = 0;

/* ------------------------------------------------------------
 * Constructors and destructors
 * ------------------------------------------------------------ */

pdclusterbasis
init_dclusterbasis(pdclusterbasis cb, pcdcluster t)
{
  uint      directions;
  uint      iota, j;

  assert(cb != NULL);

  directions = (t->directions == 0 ? 1 : t->directions);

  cb->t = t;
  cb->directions = directions;
  cb->ktree = 0;
  cb->kbranch = 0;

  cb->k = (uint *) allocmem(sizeof(uint) * directions);
  cb->koff = (uint *) allocmem(sizeof(uint) * (directions + 1));
  for (iota = 0; iota < directions; iota++) {
    cb->k[iota] = 0;
    cb->koff[iota] = 0;
  }
  cb->koff[iota] = 0;

  if (t->sons == 0) {
    cb->V = (pamatrix) allocmem(sizeof(amatrix) * directions);
    for (iota = 0; iota < directions; iota++)
      init_amatrix(cb->V + iota, t->size, 0);

    cb->E = 0;

    cb->sons = 0;
    cb->son = 0;
    cb->dirson = 0;
  }
  else {
    cb->E = (pamatrix *) allocmem(sizeof(pamatrix) * t->sons);
    cb->son = (pdclusterbasis *) allocmem(sizeof(pdclusterbasis) * t->sons);
    cb->dirson = (uint **) allocmem(sizeof(uint *) * t->sons);
    for (j = 0; j < t->sons; j++) {
      cb->son[j] = 0;
      cb->dirson[j] = (uint *) allocmem(sizeof(uint) * directions);
      cb->E[j] = (pamatrix) allocmem(sizeof(amatrix) * directions);
      for (iota = 0; iota < directions; iota++) {
	cb->dirson[j][iota] = (t->dirson ? t->dirson[j][iota] : 0);
	init_amatrix(cb->E[j] + iota, 0, 0);
      }
    }

    cb->V = 0;

    cb->sons = t->sons;
  }

#ifdef USE_OPENMP
#pragma omp atomic
#endif
  active_dclusterbasis++;

  return cb;
}

void
uninit_dclusterbasis(pdclusterbasis cb)
{
  uint      directions;
  uint      i, j;

  directions = cb->directions;

  freemem(cb->k);
  freemem(cb->koff);

  if (cb->sons > 0) {
    for (i = 0; i < cb->sons; i++) {
      del_dclusterbasis(cb->son[i]);
      freemem(cb->dirson[i]);
    }
    freemem(cb->son);
    freemem(cb->dirson);
  }

  if (cb->V) {
    for (i = 0; i < directions; i++)
      uninit_amatrix(cb->V + i);
    freemem(cb->V);
    cb->V = 0;
  }

  if (cb->E) {
    for (j = 0; j < cb->sons; j++) {
      for (i = 0; i < directions; i++)
	uninit_amatrix(cb->E[j] + i);
      freemem(cb->E[j]);
    }
    freemem(cb->E);
    cb->E = 0;
  }

  assert(active_dclusterbasis > 0);

#ifdef USE_OPENMP
#pragma omp atomic
#endif
  active_dclusterbasis--;
}

pdclusterbasis
new_dclusterbasis(pcdcluster t)
{
  pdclusterbasis cb;

  cb = allocmem(sizeof(dclusterbasis));

  init_dclusterbasis(cb, t);

  return cb;
}

void
del_dclusterbasis(pdclusterbasis cb)
{
  uninit_dclusterbasis(cb);

  freemem(cb);
}

/* ------------------------------------------------------------
 * Low-level management
 * ------------------------------------------------------------ */

void
update_dclusterbasis(pdclusterbasis cb)
{
  uint      ksum, stree, sbranch;
  uint      directions;
  uint      i;

  directions = cb->directions;

  ksum = 0;
  for (i = 0; i < directions; i++) {
    cb->koff[i] = ksum;
    ksum += cb->k[i];
  }
  cb->koff[i] = ksum;

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

  cb->ktree = stree + ksum;
  cb->kbranch = sbranch + ksum;
}

void
setrank_dclusterbasis(pdclusterbasis cb, uint iota, uint k)
{
  uint      j, iota1;

  assert(iota < cb->directions);

  if (cb->sons > 0) {
    for (j = 0; j < cb->sons; j++) {
      /* Direction iota in the current cluster corresponds
       * to direction iota1 in the j-th son */
      iota1 = cb->dirson[j][iota];

      /* Set dimensions of transfer matrix */
      resize_amatrix(cb->E[j] + iota, cb->son[j]->k[iota1], k);
    }
  }
  else {
    /* Set dimensions of leaf matrix */
    resize_amatrix(cb->V + iota, cb->t->size, k);
  }

  cb->k[iota] = k;
}

void
initmatrices_dclusterbasis(pdclusterbasis cb)
{
  uint      j, iota, iota1;

  if (cb->sons > 0) {
    for (j = 0; j < cb->sons; j++) {
      initmatrices_dclusterbasis(cb->son[j]);

      for (iota = 0; iota < cb->directions; iota++) {
	iota1 = cb->dirson[j][iota];

	resize_amatrix(cb->E[j] + iota, cb->son[j]->k[iota1], cb->k[iota]);
      }
    }
  }
  else {
    for (iota = 0; iota < cb->directions; iota++)
      resize_amatrix(cb->V + iota, cb->t->size, cb->k[iota]);
  }
}

static void
scan_blocks(uint rank, pcdblock b, pdclusterbasis rb, pdclusterbasis cb)
{

  uint      rsons, csons;
  uint      i, j;

  assert(b->rc == rb->t);
  assert(b->cc == cb->t);

  /* Mark directions in row and column cluster basis */
  if (b->adm) {
    rb->k[b->rd] = rank;
    cb->k[b->cd] = rank;
  }

  if (b->son) {
    rsons = b->rsons;
    csons = b->csons;

    if (b->son[0]->cc == b->cc) {
      assert(csons == 1);

      if (b->son[0]->rc == b->rc) {
	assert(rsons == 1);

	/* Row and column have no sons, but the block does.
	 * This should never happen, but we are prepared for it, anyway. */
	scan_blocks(rank, b->son[0], rb, cb);
      }
      else {
	assert(rsons == rb->sons);

	/* Row has sons, column does not */
	for (i = 0; i < rsons; i++)
	  scan_blocks(rank, b->son[i], rb->son[i], cb);
      }
    }
    else {
      assert(csons == cb->sons);

      if (b->son[0]->rc == b->rc) {
	assert(rsons == 1);

	/* Column has sons, row does not */
	for (j = 0; j < csons; j++)
	  scan_blocks(rank, b->son[j], rb, cb->son[j]);
      }
      else {
	assert(rsons == rb->sons);

	/* Row and column have sons */
	for (j = 0; j < csons; j++)
	  for (i = 0; i < rsons; i++)
	    scan_blocks(rank, b->son[i + j * rsons], rb->son[i], cb->son[j]);
      }
    }
  }
}

static void
scan_sons(pdclusterbasis cb)
{
  uint      directions = cb->directions;
  uint      iota, iota1, j;

  if (cb->sons > 0) {
    for (j = 0; j < cb->sons; j++) {
      for (iota = 0; iota < directions; iota++)
	if (cb->k[iota] > 0) {
	  /* This direction in the son corresponds to iota in the
	   * current cluster */
	  iota1 =
	    (cb->son[j]->t->directions > 0 ? cb->t->dirson[j][iota] : 0);
	  assert(iota1 < cb->son[j]->directions);

	  /* If a direction is used in the current cluster, the
	   * corresponding direction has to be present in the son */
	  cb->son[j]->k[iota1] = cb->k[iota];
	}

      scan_sons(cb->son[j]);
    }
  }

  update_dclusterbasis(cb);
}


void
findranks_dclusterbasis(uint rank, pcdblock b, pdclusterbasis rb,
			pdclusterbasis cb)
{
  scan_blocks(rank, b, rb, cb);

  scan_sons(rb);
  scan_sons(cb);
}

/* ------------------------------------------------------------
 * Statistics
 * ------------------------------------------------------------ */

uint
getactives_dclusterbasis()
{
  return active_dclusterbasis;
}

size_t
getsize_dclusterbasis(pcdclusterbasis cb)
{
  size_t    sz;
  uint      i;

  sz = getsize_nonrecursive_dclusterbasis(cb);

  for (i = 0; i < cb->sons; i++)
    sz += getsize_dclusterbasis(cb->son[i]);

  return sz;
}

size_t
getsize_nonrecursive_dclusterbasis(pcdclusterbasis cb)
{
  uint      directions;
  size_t    sz;
  uint      iota, j;

  directions = cb->directions;

  sz = (size_t) sizeof(dclusterbasis);

  if (cb->V) {
    for (iota = 0; iota < directions; iota++)
      sz += getsize_heap_amatrix(cb->V + iota);
    sz += sizeof(amatrix) * directions;
  }

  if (cb->E) {
    for (j = 0; j < cb->sons; j++) {
      for (iota = 0; iota < directions; iota++)
	sz += getsize_heap_amatrix(cb->E[j] + iota);
      sz += sizeof(amatrix) * directions;
    }
    sz += sizeof(pamatrix) * cb->sons;
  }

  sz += sizeof(uint) * directions;
  sz += sizeof(uint) * (directions + 1);

  if (cb->sons > 0)
    sz += (size_t) sizeof(pdclusterbasis) * cb->sons;

  return sz;
}

uint
getmaxrank_dclusterbasis(pcdclusterbasis cb)
{
  uint      kmax, iota, i;

  kmax = 0;
  for (iota = 0; iota < cb->directions; iota++)
    kmax = UINT_MAX(kmax, cb->k[iota]);

  for (i = 0; i < cb->sons; i++)
    kmax = UINT_MAX(kmax, getmaxrank_dclusterbasis(cb->son[i]));

  return kmax;
}

uint
getactivedirections_dclusterbasis(pcdclusterbasis cb)
{
  uint      d;
  uint      j, iota;

  d = 0;
  for (iota = 0; iota < cb->directions; iota++)
    if (cb->k[iota] > 0)
      d++;

  for (j = 0; j < cb->sons; j++)
    d += getactivedirections_dclusterbasis(cb->son[j]);

  return d;
}

static void
print_tree(pcdclusterbasis cb, uint level)
{
  uint      directions;
  uint      i, j;

  directions = cb->directions;

  for (i = 0; i < level; i++)
    printf("  ");
  printf("k=(%u", cb->k[0]);
  for (i = 1; i < directions; i++)
    printf(",%u", cb->k[i]);
  printf("), ktree=%u, size=%u\n", cb->ktree, cb->t->size);

  for (j = 0; j < cb->sons; j++)
    print_tree(cb->son[j], level + 1);
}

void
print_tree_dclusterbasis(pcdclusterbasis cb)
{
  print_tree(cb, 0);
}

/* ------------------------------------------------------------
 * Build cluster basis based on cluster tree
 * ------------------------------------------------------------ */

pdclusterbasis
buildfromdcluster_dclusterbasis(pcdcluster t)
{
  pdclusterbasis cb, cb1;
  uint      i;

  cb = new_dclusterbasis(t);

  for (i = 0; i < t->sons; i++) {
    cb1 = buildfromdcluster_dclusterbasis(t->son[i]);
    cb->son[i] = cb1;
  }

  update_dclusterbasis(cb);

  return cb;
}

/* ------------------------------------------------------------
 * Forward and backward transformation
 * ------------------------------------------------------------ */

pavector
newcoeffs_dclusterbasis(pcdclusterbasis cb)
{
  pavector  xt;

  xt = new_avector(cb->ktree);

  return xt;
}

void
forward_dclusterbasis(pcdclusterbasis cb, pcavector x, pavector xt)
{
  avector   tmp1, tmp2;
  pavector  xt1, xc, xp;
  pcdclusterbasis cb1;
  uint      directions;
  uint      iota, iota1, j, xtoff;

  assert(xt->dim == cb->ktree);

  directions = cb->directions;

  if (cb->sons > 0) {
    /* This part of xt contains the coefficients for all directions
     * in the current cluster */
    xc = init_sub_avector(&tmp1, xt, cb->koff[directions], 0);
    clear_avector(xc);
    uninit_avector(xc);

    xtoff = cb->koff[directions];
    for (j = 0; j < cb->sons; j++) {
      cb1 = cb->son[j];

      /* This part of xt corresponds to the subtree rooted in the
       * j-th son */
      xt1 = init_sub_avector(&tmp2, xt, cb1->ktree, xtoff);

      /* Compute coefficients in the subtree */
      forward_dclusterbasis(cb1, x, xt1);

      /* Clean up subvector */
      uninit_avector(xt1);

      for (iota = 0; iota < directions; iota++) {
	/* This part of xt corresponds to iota in the current cluster */
	xc = init_sub_avector(&tmp1, xt, cb->k[iota], cb->koff[iota]);

	/* This direction corresponds to iota in the j-th son */
	iota1 = cb->dirson[j][iota];

	/* This part of xt corresponds to the direction iota1 in
	 * the j-th son */
	xt1 = init_sub_avector(&tmp2, xt,
			       cb1->k[iota1], xtoff + cb1->koff[iota1]);

	/* Multiply by adjoint transfer matrix */
	mvm_amatrix_avector(1.0, true, cb->E[j] + iota, xt1, xc);

	/* Clean up source vector */
	uninit_avector(xt1);

	/* Clean up target vector */
	uninit_avector(xc);
      }

      /* Get offset for next son */
      xtoff += cb1->ktree;
    }
    assert(xtoff == cb->ktree);
  }
  else {
    /* Permuted entries of x */
    xp = init_sub_avector(&tmp2, xt, cb->t->size, cb->koff[directions]);

    /* Find and copy entries */
    for (j = 0; j < cb->t->size; j++)
      xp->v[j] = x->v[cb->t->idx[j]];

    for (iota = 0; iota < directions; iota++) {
      /* This part of xt contains the coefficients for the direction
       * iota in the current cluster */
      xc = init_sub_avector(&tmp1, xt, cb->k[iota], cb->koff[iota]);
      clear_avector(xc);

      /* Multiply by leaf matrix for direction iota */
      mvm_amatrix_avector(1.0, true, cb->V + iota, xp, xc);

      /* Clean up target subvector */
      uninit_avector(xc);
    }

    /* Clean up source subvector */
    uninit_avector(xp);
  }
}

void
backward_dclusterbasis(pcdclusterbasis cb, pavector yt, pavector y)
{
  avector   tmp1, tmp2;
  pavector  yt1, yc, yp;
  pcdclusterbasis cb1;
  uint      directions;
  uint      iota, iota1, j, ytoff;

  assert(yt->dim == cb->ktree);

  directions = cb->directions;

  if (cb->sons > 0) {
    ytoff = cb->koff[directions];
    for (j = 0; j < cb->sons; j++) {
      cb1 = cb->son[j];

      for (iota = 0; iota < directions; iota++) {
	/* This part of yt corresponds to iota in the current cluster */
	yc = init_sub_avector(&tmp1, yt, cb->k[iota], cb->koff[iota]);

	/* This direction corresponds to iota in the j-th son */
	iota1 = cb->dirson[j][iota];

	/* This part of yt corresponds to the direction iota1 in
	 * the j-th son */
	yt1 = init_sub_avector(&tmp2, yt,
			       cb1->k[iota1], ytoff + cb1->koff[iota1]);

	/* Multiply by transfer matrix */
	mvm_amatrix_avector(1.0, false, cb->E[j] + iota, yc, yt1);

	/* Clean up target subvector */
	uninit_avector(yt1);

	/* Clean up source subvector */
	uninit_avector(yc);
      }

      /* This part of yt corresponds to the subtree rooted in the
       * j-th son */
      yt1 = init_sub_avector(&tmp2, yt, cb1->ktree, ytoff);

      /* Propagate coefficients in the subtree */
      backward_dclusterbasis(cb1, yt1, y);

      /* Clean up source vector */
      uninit_avector(yt1);

      /* Get offset for next son */
      ytoff += cb1->ktree;
    }
    assert(ytoff == cb->ktree);
  }
  else {
    /* Permuted entries of y */
    yp = init_sub_avector(&tmp1, yt, cb->t->size, cb->koff[directions]);

    for (iota = 0; iota < directions; iota++) {
      /* This part of yt contains the coefficients for the direction
       * iota in the current cluster */
      yc = init_sub_avector(&tmp2, yt, cb->k[iota], cb->koff[iota]);

      /* Multiply by leaf matrix for direction iota */
      mvm_amatrix_avector(1.0, false, cb->V + iota, yc, yp);

      /* Clean up source subvector */
      uninit_avector(yc);
    }

    /* Add entries to target vector */
    for (j = 0; j < cb->t->size; j++)
      y->v[cb->t->idx[j]] += yp->v[j];

    /* Clean up source vector */
    uninit_avector(yp);
  }
}

void
slowforward_dclusterbasis(pcdclusterbasis cb, pcavector x, pavector xt)
{
  amatrix   tmp1, tmp2;
  avector   tmp3, tmp4;
  pamatrix  Vt, Vc, Vfull;
  pavector  xp, xc, xt1;
  uint      directions = cb->directions;
  uint      iota, j, xtoff;

  for (iota = 0; iota < directions; iota++) {
    Vt = init_amatrix(&tmp1, cb->kbranch, cb->k[iota]);
    clear_amatrix(Vt);

    Vc = init_sub_amatrix(&tmp2, Vt, cb->k[iota], cb->koff[iota],
			  cb->k[iota], 0);
    identity_amatrix(Vc);
    uninit_amatrix(Vc);

    Vfull = init_amatrix(&tmp2, cb->t->size, cb->k[iota]);
    clear_amatrix(Vfull);

    blockexpand_dclusterbasis(1.0, cb, Vt, Vfull);

    xp = init_avector(&tmp3, cb->t->size);
    for (j = 0; j < cb->t->size; j++)
      xp->v[j] = x->v[cb->t->idx[j]];

    xc = init_sub_avector(&tmp4, xt, cb->k[iota], cb->koff[iota]);
    clear_avector(xc);
    mvm_amatrix_avector(1.0, true, Vfull, xp, xc);
    uninit_avector(xc);

    if (cb->sons == 0) {
      xc = init_sub_avector(&tmp4, xt, cb->t->size, cb->koff[directions]);
      copy_avector(xp, xc);
      uninit_avector(xc);
    }

    uninit_avector(xp);
    uninit_amatrix(Vfull);
    uninit_amatrix(Vt);
  }

  if (cb->sons > 0) {
    xtoff = cb->koff[directions];
    for (j = 0; j < cb->sons; j++) {
      xt1 = init_sub_avector(&tmp3, xt, cb->son[j]->ktree, xtoff);

      slowforward_dclusterbasis(cb->son[j], x, xt1);

      uninit_avector(xt1);

      xtoff += cb->son[j]->ktree;
    }
    assert(xtoff == cb->ktree);
  }
}

void
slowbackward_dclusterbasis(pcdclusterbasis cb, pcavector yt, pavector y)
{
  amatrix   tmp1, tmp2;
  avector   tmp3, tmp4;
  pamatrix  Vt, Vc, Vfull;
  pavector  yp, yc, yt1;
  uint      directions = cb->directions;
  uint      iota, j, ytoff;

  for (iota = 0; iota < directions; iota++) {
    Vt = init_amatrix(&tmp1, cb->kbranch, cb->k[iota]);
    clear_amatrix(Vt);

    Vc = init_sub_amatrix(&tmp2, Vt, cb->k[iota], cb->koff[iota],
			  cb->k[iota], 0);
    identity_amatrix(Vc);
    uninit_amatrix(Vc);

    Vfull = init_amatrix(&tmp2, cb->t->size, cb->k[iota]);
    clear_amatrix(Vfull);

    blockexpand_dclusterbasis(1.0, cb, Vt, Vfull);

    yp = init_avector(&tmp3, cb->t->size);

    if (cb->sons == 0) {
      yc =
	init_sub_avector(&tmp4, (pavector) yt, cb->t->size,
			 cb->koff[directions]);
      copy_avector(yc, yp);
      uninit_avector(yc);
    }
    else
      clear_avector(yp);

    yc = init_sub_avector(&tmp4, (pavector) yt, cb->k[iota], cb->koff[iota]);
    mvm_amatrix_avector(1.0, false, Vfull, yc, yp);
    uninit_avector(yc);

    for (j = 0; j < cb->t->size; j++)
      y->v[cb->t->idx[j]] += yp->v[j];

    uninit_avector(yp);
    uninit_amatrix(Vfull);
    uninit_amatrix(Vt);
  }

  if (cb->sons > 0) {
    ytoff = cb->koff[directions];
    for (j = 0; j < cb->sons; j++) {
      yt1 = init_sub_avector(&tmp3, (pavector) yt, cb->son[j]->ktree, ytoff);

      slowbackward_dclusterbasis(cb->son[j], yt1, y);

      uninit_avector(yt1);

      ytoff += cb->son[j]->ktree;
    }
    assert(ytoff == cb->ktree);
  }
}

/* ------------------------------------------------------------
 * Forward transformation for the root only
 * ------------------------------------------------------------ */

void
compress_dclusterbasis(pcdclusterbasis cb, pcavector xp, pavector xt)
{
  avector   tmp1, tmp2;
  pavector  xt1, xp1, xc, xc1;
  pcdclusterbasis cb1;
  uint      directions;
  uint      iota, iota1, j, xpoff;

  assert(xt->dim == cb->kbranch);
  assert(xp->dim == cb->t->size);

  directions = cb->directions;

  /* Clear target vector */
  xc = init_sub_avector(&tmp1, xt, cb->koff[directions], 0);
  clear_avector(xc);
  uninit_avector(xc);

  if (cb->sons > 0) {
    xpoff = 0;
    for (j = 0; j < cb->sons; j++) {
      cb1 = cb->son[j];

      /* This part of xp corresponds to the j-th son */
      xp1 = init_sub_avector(&tmp1, (pavector) xp, cb1->t->size, xpoff);

      /* This part of xt corresponds to the subtree rooted in the
       * j-th son */
      xt1 = init_sub_avector(&tmp2, xt, cb1->kbranch, cb->koff[directions]);

      /* Compute coefficients in the subtree */
      compress_dclusterbasis(cb1, xp1, xt1);

      /* Clean up target subvector */
      uninit_avector(xt1);

      /* Clean up source subvector */
      uninit_avector(xp1);

      for (iota = 0; iota < directions; iota++) {
	/* This direction corresponds to iota in the j-th son */
	iota1 = cb->dirson[j][iota];

	/* This part of xt corresponds to iota in the current cluster */
	xc = init_sub_avector(&tmp1, xt, cb->k[iota], cb->koff[iota]);

	/* This part of xt corresponds to iota1 in the j-th son */
	xc1 = init_sub_avector(&tmp2, xt,
			       cb1->k[iota1],
			       cb->koff[directions] + cb1->koff[iota1]);

	/* Multiply by transfer matrix */
	mvm_amatrix_avector(1.0, true, cb->E[j] + iota, xc1, xc);

	/* Clean up source subvector */
	uninit_avector(xc1);

	/* Clean up target subvector */
	uninit_avector(xc);
      }

      /* Get offset for next son */
      xpoff += cb1->t->size;
    }
    assert(xpoff == cb->t->size);
  }
  else {
    for (iota = 0; iota < directions; iota++) {
      /* This part of xt contains the coefficients for the direction
       * iota in the current cluster */
      xc = init_sub_avector(&tmp1, xt, cb->k[iota], cb->koff[iota]);

      /* Multiply by leaf matrix for direction iota */
      mvm_amatrix_avector(1.0, true, cb->V + iota, xp, xc);

      /* Clean up subvector */
      uninit_avector(xc);
    }
  }
}

void
expand_dclusterbasis(field alpha, pcdclusterbasis cb,
		     pavector yt, pavector yp)
{
  avector   tmp1, tmp2;
  pavector  yt1, yp1, yc, yc1;
  pcdclusterbasis cb1;
  uint      directions;
  uint      iota, iota1, j, ypoff;

  assert(yt->dim == cb->kbranch);
  assert(yp->dim == cb->t->size);

  directions = cb->directions;

  if (cb->sons > 0) {
    ypoff = 0;
    for (j = 0; j < cb->sons; j++) {
      cb1 = cb->son[j];

      /* Clear target vector */
      yt1 = init_sub_avector(&tmp1, yt,
			     cb1->koff[cb1->directions],
			     cb->koff[directions]);
      clear_avector(yt1);
      uninit_avector(yt1);

      for (iota = 0; iota < directions; iota++) {
	/* This directions corresponds to iota in the j-th son */
	iota1 = cb->dirson[j][iota];

	/* This part of yt corresponds to iota in the current cluster */
	yc = init_sub_avector(&tmp1, yt, cb->k[iota], cb->koff[iota]);

	/* This part of yt corresponds to iota1 in the j-th son */
	yc1 = init_sub_avector(&tmp2, yt,
			       cb1->k[iota1],
			       cb->koff[directions] + cb1->koff[iota1]);

	/* Multiply by transfer matrix */
	mvm_amatrix_avector(1.0, false, cb->E[j] + iota, yc, yc1);

	/* Clean up target subvector */
	uninit_avector(yc1);

	/* Clean up source subvector */
	uninit_avector(yc);
      }

      /* This part of yp corresponds to the j-th son */
      yp1 = init_sub_avector(&tmp1, yp, cb1->t->size, ypoff);

      /* This part of yt corresponds to the subtree rooted in the
       * j-th son */
      yt1 = init_sub_avector(&tmp2, yt, cb1->kbranch, cb->koff[directions]);

      /* Expand subtree */
      expand_dclusterbasis(alpha, cb1, yt1, yp1);

      /* Clean up source subvector */
      uninit_avector(yt1);

      /* Clean up target subvector */
      uninit_avector(yp1);

      /* Get offset for next son */
      ypoff += cb1->t->size;
    }
    assert(ypoff == cb->t->size);
  }
  else {
    for (iota = 0; iota < directions; iota++) {
      /* This part of yt contains the coefficients for the direction
       * iota in the current cluster */
      yc = init_sub_avector(&tmp1, yt, cb->k[iota], cb->koff[iota]);

      /* Multiply by leaf matrix for direction iota */
      mvm_amatrix_avector(alpha, false, cb->V + iota, yc, yp);

      /* Clean up subvector */
      uninit_avector(yc);
    }
  }
}

void
blockcompress_dclusterbasis(pcdclusterbasis cb, pcamatrix Xp, pamatrix Xt)
{
  amatrix   tmp1, tmp2;
  pamatrix  Xp1, Xt1, Xc, Xc1;
  uint      directions;
  uint      iota, iota1, j, xpoff;

  assert(Xp->rows >= cb->t->size);
  assert(Xt->rows >= cb->kbranch);
  assert(Xp->cols == Xt->cols);

  directions = cb->directions;

  /* Clear target matrix */
  Xc = init_sub_amatrix(&tmp1, Xt, cb->koff[directions], 0, Xt->cols, 0);
  clear_amatrix(Xc);
  uninit_amatrix(Xc);

  if (cb->sons > 0) {
    xpoff = 0;
    for (j = 0; j < cb->sons; j++) {
      /* This part of Xp corresponds to the j-th son */
      Xp1 = init_sub_amatrix(&tmp1,
			     (pamatrix) Xp, cb->t->son[j]->size, xpoff,
			     Xp->cols, 0);

      /* This part of Xt corresponds to the subtree rooted in the
       * j-th son */
      Xt1 = init_sub_amatrix(&tmp2,
			     Xt, cb->son[j]->kbranch, cb->koff[directions],
			     Xt->cols, 0);

      /* Compute coefficients in the subtree */
      blockcompress_dclusterbasis(cb->son[j], Xp1, Xt1);

      /* Clean up target submatrix */
      uninit_amatrix(Xt1);

      /* Clean up source submatrix */
      uninit_amatrix(Xp1);

      for (iota = 0; iota < directions; iota++) {
	/* This direction corresponds to iota in the j-th son */
	iota1 = cb->dirson[j][iota];

	/* This part of Xt corresponds to iota in the current cluster */
	Xc = init_sub_amatrix(&tmp1, Xt,
			      cb->k[iota], cb->koff[iota], Xt->cols, 0);

	/* This part of Xt corresponds to iota1 in the j-th son */
	Xc1 = init_sub_amatrix(&tmp2, Xt,
			       cb->son[j]->k[iota1],
			       cb->koff[directions] + cb->son[j]->koff[iota1],
			       Xt->cols, 0);

	/* Multiply by transfer matrix */
	addmul_amatrix(1.0, true, cb->E[j] + iota, false, Xc1, Xc);

	/* Clean up source subvector */
	uninit_amatrix(Xc1);

	/* Clean up target subvector */
	uninit_amatrix(Xc);
      }

      xpoff += cb->son[j]->t->size;
    }
    assert(xpoff == cb->t->size);
  }
  else {
    for (iota = 0; iota < directions; iota++) {
      /* This part of Xt contains the coefficients for the direction
       * iota in the current cluster */
      Xc = init_sub_amatrix(&tmp1, Xt, cb->k[iota], cb->koff[iota],
			    Xt->cols, 0);

      /* Multiply by leaf matrix for direction iota */
      addmul_amatrix(1.0, true, cb->V + iota, false, Xp, Xc);

      /* Clean up submatrix */
      uninit_amatrix(Xc);
    }
  }
}

void
blockexpand_dclusterbasis(field alpha, pcdclusterbasis cb,
			  pamatrix Xt, pamatrix Xp)
{
  amatrix   tmp1, tmp2;
  pamatrix  Xp1, Xt1, Xc, Xc1;
  pcdclusterbasis cb1;
  uint      directions;
  uint      iota, iota1, j, xpoff;

  assert(Xp->rows == cb->t->size);
  assert(Xt->rows == cb->kbranch);
  assert(Xp->cols == Xt->cols);

  directions = cb->directions;

  if (cb->sons > 0) {
    xpoff = 0;
    for (j = 0; j < cb->sons; j++) {
      cb1 = cb->son[j];

      /* Clear target matrix */
      Xc = init_sub_amatrix(&tmp1, Xt,
			    cb1->koff[cb1->directions], cb->koff[directions],
			    Xt->cols, 0);
      clear_amatrix(Xc);
      uninit_amatrix(Xc);

      for (iota = 0; iota < directions; iota++) {
	/* This direction corresponds to iota in the j-th son */
	iota1 = cb->dirson[j][iota];

	/* This part of Xt corresponds to iota in the current cluster */
	Xc = init_sub_amatrix(&tmp1, Xt,
			      cb->k[iota], cb->koff[iota], Xt->cols, 0);

	/* This part of Xt corresponds to iota1 in the j-th son */
	Xc1 = init_sub_amatrix(&tmp2, Xt,
			       cb1->k[iota1],
			       cb->koff[directions] + cb1->koff[iota1],
			       Xt->cols, 0);

	/* Multiply by transfer matrix */
	addmul_amatrix(1.0, false, cb->E[j] + iota, false, Xc, Xc1);

	/* Clean up target subvector */
	uninit_amatrix(Xc1);

	/* Clean up source subvector */
	uninit_amatrix(Xc);
      }

      /* This part of Xp corresponds to the j-th son */
      Xp1 = init_sub_amatrix(&tmp1,
			     Xp, cb->t->son[j]->size, xpoff, Xp->cols, 0);

      /* This part of Xt corresponds to the subtree rooted in the
       * j-th son */
      Xt1 = init_sub_amatrix(&tmp2,
			     Xt, cb1->kbranch, cb->koff[directions],
			     Xt->cols, 0);

      /* Expand subtree */
      blockexpand_dclusterbasis(alpha, cb->son[j], Xt1, Xp1);

      /* Clean up source submatrix */
      uninit_amatrix(Xt1);

      /* Clean up target submatrix */
      uninit_amatrix(Xp1);

      /* Update offset */
      xpoff += cb->son[j]->t->size;
    }
    assert(xpoff == cb->t->size);
  }
  else {
    for (iota = 0; iota < directions; iota++) {
      /* This part of Xt contains the coefficients for the direction
       * iota in the current cluster */
      Xc = init_sub_amatrix(&tmp1, Xt, cb->k[iota], cb->koff[iota],
			    Xt->cols, 0);

      /* Multiply by leaf matrix for direction iota */
      addmul_amatrix(alpha, false, cb->V + iota, false, Xc, Xp);

      /* Clean up submatrix */
      uninit_amatrix(Xc);
    }
  }
}

/* ------------------------------------------------------------
 * Enumeration
 * ------------------------------------------------------------ */

static void
enumerate(pcdcluster t, uint tname, pdclusterbasis cb, pdclusterbasis * cbn)
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

pdclusterbasis *
enumerate_dclusterbasis(pcdcluster t, pdclusterbasis cb)
{
  pdclusterbasis *cbn;

  cbn =
    (pdclusterbasis *) allocmem((size_t) sizeof(pdclusterbasis) * t->desc);

  enumerate(t, 0, cb, cbn);

  return cbn;
}

/* ------------------------------------------------------------
 * Hierarchical iterator
 * ------------------------------------------------------------ */

void
iterate_dclusterbasis(pdclusterbasis cb, uint tname,
		      uint pardepth,
		      void (*pre) (pdclusterbasis cb, uint tname,
				   uint pardepth, void *data),
		      void (*post) (pdclusterbasis cb, uint tname,
				    uint pardepth, void *data), void *data)
{
  uint     *tname1, tname2;
  pcdcluster t = cb->t;
  uint      i;
  int       p, nthreads;

  if (pre)
    pre(cb, tname, pardepth, data);

  if (cb->sons > 0) {
    tname1 = (uint *) allocmem((size_t) sizeof(uint) * cb->sons);

    tname2 = tname + 1;
    for (i = 0; i < t->sons; i++) {
      tname1[i] = tname2;
      tname2 += t->son[i]->desc;
    }
    assert(tname2 == tname + t->desc);

    nthreads = cb->sons;
#ifdef USE_OPENMP
#pragma omp parallel for if(pardepth > 0), num_threads(nthreads)
#endif
    for (p = 0; p < nthreads; p++)
      iterate_dclusterbasis(cb->son[p], tname1[p],
			    (pardepth > 0 ? pardepth - 1 : 0),
			    pre, post, data);

    freemem(tname1);
  }

  if (post)
    post(cb, tname, pardepth, data);
}


/* ------------------------------------------------------------
 * Orthogonalization
 * ------------------------------------------------------------ */

typedef struct _orthodata orthodata;
typedef orthodata *porthodata;
struct _orthodata {

  pdclusteroperator *nco;

};

static void
ortho_dclusterbasis_operator(pdclusterbasis cb, uint tname,
			     uint pardepth, void *data)
{

  porthodata odata = (porthodata) data;
  pdclusteroperator co = odata->nco[tname];
  amatrix   A1, A2;
  avector   b1;
  pamatrix  Vhat, Qhat, Vhat1, Qhat1;
  pavector  tau;
  uint      i, j, roff, rows, dim;
  uint      temp;
  uint      direction = cb->directions;

  assert(cb->t == co->t);

  if (cb->sons > 0) {
    /* For every single direction */
    for (j = 0; j < direction; j++) {
      rows = 0;
      for (i = 0; i < cb->sons; i++) {
	rows += co->son[i]->krow[cb->dirson[i][j]];
      }

      if (rows > 0 && cb->k[j] > 0) {
	Vhat = init_amatrix(&A1, rows, cb->k[j]);	/* Initialize matrix with required size */
	clear_amatrix(Vhat);

	roff = 0;
	/* Compute new Vhat_t = R_tE_t for every single son and save this matrix as part of the Vhat */
	for (i = 0; i < cb->sons; i++) {
	  temp = cb->dirson[i][j];	/* For one term save the needed direction */

	  Vhat1 =
	    init_sub_amatrix(&A2, Vhat, cb->son[i]->k[temp], roff, cb->k[j],
			     0);
	  clear_amatrix(Vhat1);
	  assert(co->son[i]->krow[temp] == cb->son[i]->k[temp]);	/* It should be the same size */
	  assert(co->son[i]->kcol[temp] == cb->E[i][j].rows);
	  addmul_amatrix(1.0, false, &co->son[i]->C[temp], false,
			 &cb->E[i][j], Vhat1);
	  uninit_amatrix(Vhat1);
	  roff += cb->son[i]->k[temp];
	}
	assert(roff == rows);
	dim = UINT_MIN(rows, cb->k[j]);

	tau = init_avector(&b1, dim);
	qrdecomp_amatrix(Vhat, tau);	/* Solve qr decomposition */

	resize_dclusteroperator(co, dim, Vhat->cols, j);	/* Prepare for restructure  Vhat->cols should be cb->k[j] */
	copy_upper_amatrix(Vhat, false, &co->C[j]);	/* Rewrite C[j] with our R_t */

	setrank_dclusterbasis(cb, j, dim);	/* Next step will be Qhat so, that we get the E_t */
	Qhat = init_amatrix(&A2, rows, dim);
	qrexpand_amatrix(Vhat, tau, Qhat);
	uninit_amatrix(Vhat);
	uninit_avector(tau);

	roff = 0;		/* Write every single E_t for the direction j of this clusterbasis */
	for (i = 0; i < cb->sons; i++) {
	  temp = cb->dirson[i][j];
	  Qhat1 =
	    init_sub_amatrix(&A1, Qhat, cb->son[i]->k[temp], roff, dim, 0);
	  assert(cb->E[i][j].rows == cb->son[i]->k[temp]);
	  assert(cb->E[i][j].cols == cb->k[j]);
	  copy_amatrix(false, Qhat1, &cb->E[i][j]);
	  uninit_amatrix(Qhat1);
	  roff += cb->son[i]->k[temp];
	}
	assert(roff == rows);

	uninit_amatrix(Qhat);
	update_dclusterbasis(cb);

      }
      else {
	resize_dclusteroperator(co, 0, cb->k[j], j);
	setrank_dclusterbasis(cb, j, 0);
	update_dclusterbasis(cb);
      }
    }
  }
  else {			/* Leaf */

    rows = cb->t->size;		/* Get size for all V[i] */

    for (i = 0; i < direction; i++) {	/* Run through all directions and evaluate QR decomposition */
      if (cb->k[i] > 0) {
	Vhat = init_amatrix(&A1, rows, cb->k[i]);	/* Set up Vhat rows = size of cluster */
	/* Column = rank of underlying direction */
	copy_amatrix(false, &cb->V[i], Vhat);

	dim = UINT_MIN(rows, cb->k[i]);	/* Set up tau */
	tau = init_avector(&b1, dim);

	qrdecomp_amatrix(Vhat, tau);

	resize_dclusteroperator(co, dim, cb->k[i], i);	/* Change size of clusteroperator (Attention co->V will be destroyed!) */
	setrank_dclusterbasis(cb, i, dim);	/* Change size of clusterbasis */

	copy_upper_amatrix(Vhat, false, &co->C[i]);	/* Reearn R */
	qrexpand_amatrix(Vhat, tau, &cb->V[i]);	/* Reearn Q */

	uninit_avector(tau);
	uninit_amatrix(Vhat);
      }
      else {
	resize_dclusteroperator(co, 0, cb->k[i], i);	/* Change size of clusteroperator (Attention co->V will be destroyed!) */
	setrank_dclusterbasis(cb, i, 0);
      }
    }
    update_dclusterbasis(cb);
  }
}


void
ortho_dclusterbasis(pdclusterbasis cb, pdclusteroperator co)
{

  orthodata odata;

  odata.nco = enumerate_dclusteroperator(cb->t, co);
  iterate_dclusterbasis(cb, 0, max_pardepth, NULL,
			ortho_dclusterbasis_operator, (void *) &odata);

  freemem(odata.nco);

}


static void
weight(pdclusterbasis cb, uint tname, uint pardepth, void *data)
{

  porthodata odata = (porthodata) data;
  pdclusteroperator co = odata->nco[tname];
  amatrix   A1, A2;
  avector   b1;
  pamatrix  Vhat, Vhat1;
  pavector  tau;
  uint      i, j, roff, rows, dim;
  uint      temp;
  uint      direction = cb->directions;

  assert(cb->t == co->t);


  if (cb->sons > 0) {
    /* Every single direction */
    for (j = 0; j < direction; j++) {
      rows = 0;
      for (i = 0; i < cb->sons; i++) {
	rows += co->son[i]->krow[cb->dirson[i][j]];
      }

      Vhat = init_amatrix(&A1, rows, cb->k[j]);	/* Initialize matrix with required size */
      clear_amatrix(Vhat);

      roff = 0;
      /* Compute new Vhat_t = R_tE_t for every single son and save this matrix as part of the Vhat */
      for (i = 0; i < cb->sons; i++) {
	temp = cb->dirson[i][j];	/* For one term save the needed direction */

	Vhat1 =
	  init_sub_amatrix(&A2, Vhat, co->son[i]->krow[temp], roff, cb->k[j],
			   0);
	clear_amatrix(Vhat1);
	addmul_amatrix(1.0, false, &co->son[i]->C[temp], false, &cb->E[i][j],
		       Vhat1);
	uninit_amatrix(Vhat1);
	roff += co->son[i]->krow[temp];
      }
      assert(roff == rows);
      dim = UINT_MIN(rows, cb->k[j]);

      tau = init_avector(&b1, dim);
      qrdecomp_amatrix(Vhat, tau);	/* Solve qr decomposition */

      resize_dclusteroperator(co, dim, Vhat->cols, j);	/* Prepare for restructure  Vhat->cols should be cb->k[j] */
      copy_upper_amatrix(Vhat, false, &co->C[j]);	/* Rewrite C[j] with our R_t */

      uninit_amatrix(Vhat);
      uninit_avector(tau);
    }
  }
  else {			/* Leaf */

    rows = cb->t->size;		/* Get size for all V[i] */

    for (i = 0; i < direction; i++) {	/* Run through all directions and evaluate QR decomposition */
      Vhat = init_amatrix(&A1, rows, cb->k[i]);	/* Set up Vhat rows = size of cluster */
      /* Column = rank of underlying direction */
      copy_amatrix(false, &cb->V[i], Vhat);

      dim = UINT_MIN(rows, cb->k[i]);	/* Set up tau */
      tau = init_avector(&b1, dim);

      qrdecomp_amatrix(Vhat, tau);
      resize_dclusteroperator(co, dim, cb->k[i], i);	/* Change size of clusteroperator (Attention co->V will be destroyed!) */
      copy_upper_amatrix(Vhat, false, &co->C[i]);	/* Reearn R */

      uninit_avector(tau);
      uninit_amatrix(Vhat);
    }
  }
}


void
weight_dclusterbasis_dclusteroperator(pdclusterbasis cb, pdclusteroperator co)
{

  orthodata odata;

  odata.nco = enumerate_dclusteroperator(cb->t, co);
  iterate_dclusterbasis(cb, 0, max_pardepth, NULL, weight, (void *) &odata);

  freemem(odata.nco);

}


real
check_ortho_dclusterbasis(pcdclusterbasis cb)
{
  /* The idea is to evaluate the Frobenius Norm of E and V */

  amatrix   tmp;
  pamatrix  A;
  real      norm = 0.0;
  real      tmp_norm = 0.0;
  uint      i, j, l;
  uint      k = 0;

  /* Find the largest rank */
  for (j = 0; j < cb->directions; j++) {
    k = (k > cb->k[j] ? k : cb->k[j]);
  }

  A = init_amatrix(&tmp, k, k);

  /* If the cluster has sons I should walk through until I reach the root and on my way back I evaluate the norm of the E matrices  */

  if (cb->sons > 0) {
    for (i = 0; i < cb->sons; i++) {
      tmp_norm = check_ortho_dclusterbasis(cb->son[i]);
      norm = (norm > tmp_norm ? norm : tmp_norm);
    }
    /* For a given orthogonal matrix D is D*D = Id true, so reconstruct the whole Q and evaluate Q*Q */
    for (j = 0; j < cb->directions; j++) {
      clear_amatrix(A);
      for (i = 0; i < cb->sons; i++) {
	addmul_amatrix(1.0, true, &cb->E[i][j], false, &cb->E[i][j], A);

      }
      for (l = 0; l < cb->k[j]; l++) {
	addentry_amatrix(A, l, l, -1.0);
      }
      tmp_norm = normfrob_amatrix(A);
      norm = (norm > tmp_norm ? norm : tmp_norm);
    }

  }
  else {
    for (j = 0; j < cb->directions; j++) {
      clear_amatrix(A);
      addmul_amatrix(1.0, true, &cb->V[j], false, &cb->V[j], A);

      for (i = 0; i < cb->V[j].cols; i++) {
	addentry_amatrix(A, i, i, -1.0);
      }
      tmp_norm = normfrob_amatrix(A);
      norm = (norm > tmp_norm ? norm : tmp_norm);
    }
  }
  uninit_amatrix(A);
  return norm;
}

/* --------------------------------------- */
/*  Copy                                   */
/* --------------------------------------- */



static void
copy_dclusterbasis(pcdclusterbasis cb, pdclusterbasis copy)
{

  uint      i, j;
  uint      directions = cb->directions;

  assert(cb->sons == copy->sons);
  assert(cb->t == copy->t);
  assert(cb->directions == copy->directions);

  if (cb->sons > 0) {
    for (i = 0; i < cb->sons; i++) {
      copy_dclusterbasis(cb->son[i], copy->son[i]);

      for (j = 0; j < directions; j++) {
	//printf("non-leaf\n");        
	setrank_dclusterbasis(copy, j, cb->k[j]);
	//printf(" rows %u %u, cols %u %u  %u \n", cb->E[i][j].rows, copy->E[i][j].rows, cb->E[i][j].cols, copy->E[i][j].cols, copy->son[i]->k[copy->dirson[i][j]]);          
	copy_amatrix(false, &cb->E[i][j], &copy->E[i][j]);
      }
    }
  }
  else {
    for (j = 0; j < directions; j++) {
      //printf("leaf\n");        
      setrank_dclusterbasis(copy, j, cb->k[j]);
      copy_amatrix(false, &cb->V[j], &copy->V[j]);
    }
  }

  update_dclusterbasis(copy);
}

pdclusterbasis
clone_structure_dclusterbasis(pcdclusterbasis cb)
{

  uint      i;
  pdclusterbasis clone, clone1;

  clone = new_dclusterbasis(cb->t);

  if (cb->sons > 0) {

    assert(clone->sons == cb->sons);
    assert(clone->directions == cb->directions);

    for (i = 0; i < cb->sons; i++) {
      clone1 = clone_structure_dclusterbasis(cb->son[i]);
      clone->son[i] = clone1;
    }
  }
  update_dclusterbasis(clone);

  return clone;
}

pdclusterbasis
duplicate_dclusterbasis(pcdclusterbasis cb)
{

  pdclusterbasis cb2;

  cb2 = clone_structure_dclusterbasis(cb);
  initmatrices_dclusterbasis(cb2);
  copy_dclusterbasis(cb, cb2);

  return cb2;
}
