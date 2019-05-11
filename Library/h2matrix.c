/* ------------------------------------------------------------
 * This is the file "h2matrix.c" of the H2Lib package.
 * All rights reserved, Steffen Boerm 2009
 * ------------------------------------------------------------ */

#include <stdio.h>

#include "h2matrix.h"

#include "basic.h"

#ifdef USE_NETCDF
#include <stdio.h>
#include <string.h>
#include <netcdf.h>
#endif

/* ------------------------------------------------------------
 * Constructors and destructors
 * ------------------------------------------------------------ */

ph2matrix
new_h2matrix(pclusterbasis rb, pclusterbasis cb)
{
  ph2matrix h2;

  h2 = allocmem(sizeof(h2matrix));

  h2->rb = h2->cb = NULL;

  ref_clusterbasis(&h2->rb, rb);
  ref_clusterbasis(&h2->cb, cb);
  h2->u = NULL;
  h2->f = NULL;

  h2->son = NULL;
  h2->rsons = 0;
  h2->csons = 0;

  h2->refs = 0;
  h2->desc = 0;

  return h2;
}

ph2matrix
new_uniform_h2matrix(pclusterbasis rb, pclusterbasis cb)
{
  ph2matrix h2;

  h2 = new_h2matrix(rb, cb);

  h2->u = new_uniform(rb, cb);

  h2->desc = 1;

  return h2;
}

ph2matrix
new_full_h2matrix(pclusterbasis rb, pclusterbasis cb)
{
  ph2matrix h2;

  h2 = new_h2matrix(rb, cb);

  h2->f = new_amatrix(rb->t->size, cb->t->size);

  h2->desc = 1;

  return h2;
}

ph2matrix
new_super_h2matrix(pclusterbasis rb, pclusterbasis cb, uint rsons, uint csons)
{
  ph2matrix h2;
  uint      i, j;

  h2 = new_h2matrix(rb, cb);

  h2->rsons = rsons;
  h2->csons = csons;

  h2->son =
    (ph2matrix *) allocmem((size_t) sizeof(ph2matrix) * rsons * csons);
  for (j = 0; j < csons; j++)
    for (i = 0; i < rsons; i++)
      h2->son[i + j * rsons] = NULL;

  return h2;
}

ph2matrix
new_zero_h2matrix(pclusterbasis rb, pclusterbasis cb)
{
  ph2matrix h2;

  h2 = new_h2matrix(rb, cb);

  h2->u = 0;
  h2->f = 0;
  h2->son = 0;
  h2->desc = 1;

  return h2;
}

ph2matrix
clonestructure_h2matrix(pch2matrix h2, pclusterbasis rb, pclusterbasis cb)
{
  ph2matrix h, h1;
  ph2matrix h21;
  pclusterbasis rb1, cb1;
  uint      rsons, csons;
  uint      i, j;

  h = NULL;

  if (h2->son) {
    rsons = h2->rsons;
    csons = h2->csons;

    h = new_super_h2matrix(rb, cb, rsons, csons);

    for (j = 0; j < csons; j++)
      for (i = 0; i < rsons; i++) {
	h21 = h2->son[i + j * rsons];

	rb1 = rb;
	if (h21->rb->t != h2->rb->t) {
	  assert(rb->sons == rsons);
	  rb1 = rb->son[i];
	}

	cb1 = cb;
	if (h21->cb->t != h2->cb->t) {
	  assert(cb->sons == csons);
	  cb1 = cb->son[j];
	}

	h1 = clonestructure_h2matrix(h21, rb1, cb1);

	ref_h2matrix(h->son + i + j * rsons, h1);
      }

  }
  else if (h2->u) {
    h = new_uniform_h2matrix(rb, cb);
    clear_amatrix(&h->u->S);
  }
  else if (h2->f) {
    h = new_full_h2matrix(rb, cb);
    clear_amatrix(h->f);
  }
  else
    h = new_zero_h2matrix(rb, cb);

  update_h2matrix(h);

  return h;
}

ph2matrix
clone_h2matrix(pch2matrix h2, pclusterbasis rb, pclusterbasis cb)
{
  uint      rsons, csons;

  ph2matrix h2clone;

  ph2matrix tmp;
  pclusterbasis rb1, cb1;
  uint      i, j;

  assert(h2->rb->t == rb->t);
  assert(h2->cb->t == cb->t);

  if (h2->son != NULL) {
    rsons = h2->rsons;
    csons = h2->csons;
    h2clone = new_super_h2matrix(rb, cb, rsons, csons);

    for (j = 0; j < csons; j++) {
      if (h2->cb == h2->son[j * rsons]->cb)
	cb1 = cb;
      else
	cb1 = cb->son[j];

      for (i = 0; i < rsons; i++) {
	if (h2->rb == h2->son[i]->rb)
	  rb1 = rb;
	else
	  rb1 = rb->son[i];

	tmp = clone_h2matrix(h2->son[i + j * rsons], rb1, cb1);
	ref_h2matrix(h2clone->son + i + j * rsons, tmp);
      }
    }
  }
  else if (h2->u != NULL) {
    h2clone = new_uniform_h2matrix(rb, cb);
    copy_amatrix(false, &h2->u->S, &h2clone->u->S);
  }
  else if (h2->f != NULL) {
    h2clone = new_full_h2matrix(rb, cb);
    copy_amatrix(false, h2->f, h2clone->f);
  }
  else {
    h2clone = new_zero_h2matrix(rb, cb);
  }

  update_h2matrix(h2clone);

  return h2clone;
}

void
update_h2matrix(ph2matrix h2)
{
  uint      desc;
  uint      rsons, csons;
  uint      i, j;

  desc = 1;

  if (h2->son) {
    rsons = h2->rsons;
    csons = h2->csons;

    for (j = 0; j < csons; j++)
      for (i = 0; i < rsons; i++)
	desc += h2->son[i + j * rsons]->desc;
  }

  h2->desc = desc;
}

void
del_h2matrix(ph2matrix h2)
{
  uint      rsons = h2->rsons;
  uint      csons = h2->csons;
  uint      i, j;

  assert(h2->refs == 0);

  if (h2->son) {
    for (j = 0; j < csons; j++)
      for (i = 0; i < rsons; i++)
	unref_h2matrix(h2->son[i + j * rsons]);
    freemem(h2->son);
  }

  if (h2->f)
    del_amatrix(h2->f);

  if (h2->u)
    del_uniform(h2->u);

  unref_clusterbasis(h2->cb);
  unref_clusterbasis(h2->rb);

  freemem(h2);
}

/* ------------------------------------------------------------
 * Reference counting
 * ------------------------------------------------------------ */

void
ref_h2matrix(ph2matrix *ptr, ph2matrix h2)
{
  if (*ptr)
    unref_h2matrix(*ptr);

  *ptr = h2;

  if (h2)
    h2->refs++;
}

void
unref_h2matrix(ph2matrix h2)
{
  assert(h2->refs > 0);

  h2->refs--;

  if (h2->refs == 0)
    del_h2matrix(h2);
}

/* ------------------------------------------------------------
 * Statistics
 * ------------------------------------------------------------ */

size_t
getsize_h2matrix(pch2matrix h2)
{
  size_t    sz;
  uint      rsons = h2->rsons;
  uint      csons = h2->csons;
  uint      i, j;

  sz = (size_t) sizeof(h2matrix);

  if (h2->u)
    sz += getsize_uniform(h2->u);

  if (h2->f)
    sz += getsize_amatrix(h2->f);

  for (j = 0; j < csons; j++)
    for (i = 0; i < rsons; i++)
      sz += getsize_h2matrix(h2->son[i + j * rsons]);

  return sz;
}

size_t
gettotalsize_h2matrix(pch2matrix h2)
{
  size_t    sz;

  sz = getsize_h2matrix(h2);
  sz += getsize_clusterbasis(h2->rb);
  if (h2->rb != h2->cb)
    sz += getsize_clusterbasis(h2->cb);

  return sz;
}

size_t
getnearsize_h2matrix(pch2matrix h2)
{
  size_t    sz;
  uint      rsons = h2->rsons;
  uint      csons = h2->csons;
  uint      i, j;

  sz = 0;

  if (h2->f)
    sz += getsize_amatrix(h2->f);

  for (j = 0; j < csons; j++)
    for (i = 0; i < rsons; i++)
      sz += getnearsize_h2matrix(h2->son[i + j * rsons]);

  return sz;
}

size_t
getfarsize_h2matrix(pch2matrix h2)
{
  size_t    sz;
  uint      rsons = h2->rsons;
  uint      csons = h2->csons;
  uint      i, j;

  sz = 0;

  if (h2->u)
    sz += getsize_uniform(h2->u);

  for (j = 0; j < csons; j++)
    for (i = 0; i < rsons; i++)
      sz += getfarsize_h2matrix(h2->son[i + j * rsons]);

  return sz;
}

/* ------------------------------------------------------------
 * Simple utility functions
 * ------------------------------------------------------------ */

void
clear_h2matrix(ph2matrix h2)
{
  uint      rsons, csons;
  uint      i, j;

  if (h2->son) {
    rsons = h2->rsons;
    csons = h2->csons;

    for (j = 0; j < csons; j++)
      for (i = 0; i < rsons; i++)
	clear_h2matrix(h2->son[i + j * rsons]);
  }
  else if (h2->u)
    clear_amatrix(&h2->u->S);
  else if (h2->f)
    clear_amatrix(h2->f);
}

void
scale_h2matrix(field alpha, ph2matrix h2)
{
  uint      rsons = h2->rsons;
  uint      csons = h2->csons;
  uint      i, j;

  if (h2->son) {
    for (j = 0; j < csons; j++)
      for (i = 0; i < rsons; i++)
	scale_h2matrix(alpha, h2->son[i + j * rsons]);
  }
  else if (h2->u)
    scale_uniform(alpha, h2->u);
  else if (h2->f)
    scale_amatrix(alpha, h2->f);
}

void
random_h2matrix(ph2matrix h2)
{
  uint      rsons = h2->rsons;
  uint      csons = h2->csons;
  uint      i, j;

  if (h2->son) {
    for (j = 0; j < csons; j++)
      for (i = 0; i < rsons; i++)
	random_h2matrix(h2->son[i + j * rsons]);
  }
  else if (h2->u)
    random_uniform(h2->u);
  else if (h2->f)
    random_amatrix(h2->f);
}

/* ------------------------------------------------------------
 * Build H^2-matrix based on block tree
 * ------------------------------------------------------------ */

ph2matrix
build_from_block_h2matrix(pcblock b, pclusterbasis rb, pclusterbasis cb)
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

    for (j = 0; j < csons; j++)
      for (i = 0; i < rsons; i++) {
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

	h1 = build_from_block_h2matrix(b1, rb1, cb1);

	ref_h2matrix(h->son + i + j * rsons, h1);
      }
  }
  else if (b->a > 0)
    h = new_uniform_h2matrix(rb, cb);
  else
    h = new_full_h2matrix(rb, cb);

  update_h2matrix(h);

  return h;
}

/* ------------------------------------------------------------
 * Build block tree from H^2-matrix
 * ------------------------------------------------------------ */

pblock
build_from_h2matrix_block(pch2matrix G)
{
  uint      i, j, rsons, csons;
  pblock    b;

  b = 0;
  if (G->son) {
    rsons = G->rsons;
    csons = G->csons;

    b = new_block((pcluster) G->rb->t, (pcluster) G->cb->t, false, rsons,
		  csons);

    for (j = 0; j < csons; j++) {
      for (i = 0; i < rsons; i++) {
	b->son[i + j * rsons] =
	  build_from_h2matrix_block(G->son[i + j * rsons]);
      }
    }
  }
  else if (G->u) {
    b = new_block((pcluster) G->rb->t, (pcluster) G->cb->t, true, 0, 0);
  }
  else {
    assert(G->f);
    b = new_block((pcluster) G->rb->t, (pcluster) G->cb->t, false, 0, 0);
  }

  update_block(b);

  return b;
}

/* ------------------------------------------------------------
 * Enumeration by block number
 * ------------------------------------------------------------ */

static void
enumerate(uint bname, ph2matrix h2, ph2matrix *h2n)
{
  uint      bname1;
  uint      i, j;

  h2n[bname] = h2;

  bname1 = bname + 1;

  for (j = 0; j < h2->csons; j++) {
    for (i = 0; i < h2->rsons; i++) {
      enumerate(bname1, h2->son[i + j * h2->rsons], h2n);

      bname1 += h2->son[i + j * h2->rsons]->desc;
    }
  }

  assert(bname1 == bname + h2->desc);
}

ph2matrix *
enumerate_h2matrix(ph2matrix h2)
{
  ph2matrix *h2n;

  h2n = (ph2matrix *) allocmem((size_t) sizeof(ph2matrix) * h2->desc);

  enumerate(0, h2, h2n);

  return h2n;
}

/* ------------------------------------------------------------
 * Hierarchical iterators
 * ------------------------------------------------------------ */

void
iterate_h2matrix(ph2matrix G, uint mname, uint rname, uint cname,
		 uint pardepth,
		 void (*pre) (ph2matrix G, uint mname, uint rname, uint cname,
			      uint pardepth, void *data),
		 void (*post) (ph2matrix G, uint mname, uint rname,
			       uint cname, uint pardepth, void *data),
		 void *data)
{
  pccluster rc, cc;
  uint      rsons, csons;
  ph2matrix *Gn;
  uint     *mnamen, *rnamen, *cnamen;
  uint      mname1, rname1, cname1;
#ifdef USE_OPENMP
  uint      nthreads;		/* HACK: Solaris cc requires an lvalue */
#endif
  uint      i, j, k, l;

  /* Call a priori callback function */
  if (pre)
    pre(G, mname, rname, cname, pardepth, data);

  /* Handle sons */
  if (G->son) {
    rc = G->rb->t;
    cc = G->cb->t;

    rsons = G->rsons;
    csons = G->csons;

    if (G->rb == G->son[0]->rb) {
      if (G->cb == G->son[0]->cb) {
	/* Just one son */
	iterate_h2matrix(G->son[0], mname + 1, rname, cname,
			 (pardepth > 0 ? pardepth - 1 : 0), pre, post, data);
      }
      else {
	/* No row son, but column sons */
	mname1 = mname + 1;
	cname1 = cname + 1;
	for (j = 0; j < csons; j++) {
	  iterate_h2matrix(G->son[j * rsons], mname1, rname, cname1,
			   (pardepth > 0 ? pardepth - 1 : 0), pre, post,
			   data);
	  mname1 += G->son[j * rsons]->desc;
	  cname1 += cc->son[j]->desc;
	}
	assert(mname1 == mname + G->desc);
	assert(cname1 == cname + cc->desc);
      }
    }
    else {
      if (G->cb == G->son[0]->cb) {
	/* No column son, but row sons */
	mname1 = mname + 1;
	rname1 = rname + 1;
	for (i = 0; i < rsons; i++) {
	  iterate_h2matrix(G->son[i], mname1, rname1, cname,
			   (pardepth > 0 ? pardepth - 1 : 0), pre, post,
			   data);
	  mname1 += G->son[i]->desc;
	  rname1 += rc->son[i]->desc;
	}
	assert(mname1 == mname + G->desc);
	assert(rname1 == rname + rc->desc);
      }
      else {
	/* Row and column sons, so we can parallelize something */
	if (rsons >= csons) {
	  /* Allocate arrays for submatrices and names */
	  Gn = (ph2matrix *) allocmem(sizeof(ph2matrix) * csons);
	  mnamen = (uint *) allocmem(sizeof(uint) * csons);
	  rnamen = (uint *) allocmem(sizeof(uint) * csons);
	  cnamen = (uint *) allocmem(sizeof(uint) * csons);

	  /* Consider all parallelizable subsets, in this case we use
	     row-wise shifted cyclic diagonals */
	  for (k = 0; k < rsons; k++) {
	    cname1 = cname + 1;

	    for (j = 0; j < csons; j++) {
	      i = (k + j) % rsons;

	      Gn[j] = G->son[i + j * rsons];

	      /* Compute submatrix number */
	      mname1 = mname + 1;
	      for (l = 0; l < i + j * rsons; l++)
		mname1 += G->son[l]->desc;

	      /* Compute row number */
	      rname1 = rname;
	      if (G->rb != G->son[0]->rb) {
		rname1 = rname + 1;
		for (l = 0; l < i; l++)
		  rname1 += rc->son[l]->desc;
	      }

	      mnamen[j] = mname1;
	      rnamen[j] = rname1;
	      cnamen[j] = cname1;

	      cname1 += cc->son[j]->desc;
	    }
	    assert(cname1 == cname + cc->desc);

	    /* Recursive call */
#ifdef USE_OPENMP
	    nthreads = csons;
	    (void) nthreads;
#pragma omp parallel for if(pardepth>0),num_threads(nthreads)
#endif
	    for (j = 0; j < csons; j++)
	      iterate_h2matrix(Gn[j], mnamen[j], rnamen[j], cnamen[j],
			       (pardepth > 0 ? pardepth - 1 : 0), pre, post,
			       data);
	  }

	  /* Clean up */
	  freemem(cnamen);
	  freemem(rnamen);
	  freemem(mnamen);
	  freemem(Gn);
	}
	else {
	  /* Allocate arrays for submatrices and names */
	  Gn = (ph2matrix *) allocmem(sizeof(ph2matrix) * rsons);
	  mnamen = (uint *) allocmem(sizeof(uint) * rsons);
	  rnamen = (uint *) allocmem(sizeof(uint) * rsons);
	  cnamen = (uint *) allocmem(sizeof(uint) * rsons);

	  /* Consider all parallelizable subsets, in this case we use
	     column-wise shifted cyclic diagonals */
	  for (k = 0; k < csons; k++) {
	    rname1 = rname + 1;

	    for (i = 0; i < rsons; i++) {
	      j = (k + i) % csons;

	      Gn[i] = G->son[i + j * rsons];

	      /* Compute submatrix number */
	      mname1 = mname + 1;
	      for (l = 0; l < i + j * rsons; l++)
		mname1 += G->son[l]->desc;

	      /* Compute colum number */
	      cname1 = cname;
	      if (G->cb != G->son[0]->cb) {
		cname1 = cname + 1;
		for (l = 0; l < j; l++)
		  cname1 += cc->son[l]->desc;
	      }

	      mnamen[i] = mname1;
	      rnamen[i] = rname1;
	      cnamen[i] = cname1;

	      rname1 += rc->son[i]->desc;
	    }
	    assert(rname1 == rname + rc->desc);

	    /* Recursive call */
#ifdef USE_OPENMP
	    nthreads = rsons;
	    (void) nthreads;
#pragma omp parallel for if(pardepth>0),num_threads(nthreads)
#endif
	    for (j = 0; j < rsons; j++)
	      iterate_h2matrix(Gn[j], mnamen[j], rnamen[j], cnamen[j],
			       (pardepth > 0 ? pardepth - 1 : 0), pre, post,
			       data);
	  }

	  /* Clean up */
	  freemem(cnamen);
	  freemem(rnamen);
	  freemem(mnamen);
	  freemem(Gn);
	}
      }
    }
  }

  /* Call a posteriori callback function */
  if (post)
    post(G, mname, rname, cname, pardepth, data);
}

static void
del_h2matrixlist(ph2matrixlist hl)
{
  ph2matrixlist next;

  while (hl) {
    next = hl->next;
    freemem(hl);
    hl = next;
  }
}

static    ph2matrixlist
addrow_iterate(ph2matrix G, uint mname, uint rname,
	       uint cname, pch2matrixlist father, ph2matrixlist next)
{
  ph2matrixlist hl;
  pch2matrixlist hlf;
  pccluster cc;
  uint      mname1, cname1;
  uint      j;

  hl = (ph2matrixlist) allocmem(sizeof(h2matrixlist));
  hl->G = G;
  hl->mname = mname;
  hl->rname = rname;
  hl->cname = cname;
  hl->father = father;
  hl->next = next;

  hlf = hl;

  if (G->son && G->son[0]->rb == G->rb) {
    cc = G->cb->t;

    assert(G->csons == cc->sons);

    mname1 = mname + 1;
    cname1 = cname + 1;

    for (j = 0; j < G->csons; j++) {
      hl =
	addrow_iterate(G->son[j * G->rsons], mname1, rname, cname1, hlf, hl);

      mname1 += G->son[j * G->rsons]->desc;
      cname1 += cc->son[j]->desc;
    }
    assert(mname1 == mname + G->desc);
    assert(cname1 == cname + cc->desc);
  }

  return hl;
}

static void
iterate_rowlist(pccluster rc, uint rname, ph2matrixlist hl,
		uint pardepth,
		void (*pre) (pccluster t, uint tname, uint pardepth,
			     pch2matrixlist hl, void *data),
		void (*post) (pccluster t, uint tname, uint pardepth,
			      pch2matrixlist hl, void *data), void *data)
{
  ph2matrixlist hl0, *hl1;
  ph2matrix G;
  uint      rsons, csons;
  uint      mname1, rname1, cname1;
  uint     *rnames;
#ifdef USE_OPENMP
  uint      nthreads;		/* HACK: Solaris workaround */
#endif
  uint      i, j;

  /* Quick exit */
  if (hl == 0)
    return;

  /* Use "pre" callback, if one is given */
  if (pre)
    pre(rc, rname, pardepth, hl, data);

  /* Consider sons */
  if (rc->son) {
    /* Initialize matrix lists for all sons */
    hl1 = (ph2matrixlist *) allocmem((size_t) sizeof(ph2matrixlist) *
				     rc->sons);
    for (i = 0; i < rc->sons; i++)
      hl1[i] = 0;

    /* Initialize names for all sons */
    rnames = (uint *) allocmem((size_t) sizeof(uint) * rc->sons);
    rname1 = rname + 1;
    for (i = 0; i < rc->sons; i++) {
      rnames[i] = rname1;
      rname1 += rc->son[i]->desc;
    }
    assert(rname1 == rname + rc->desc);

    /* Fill matrix lists for all sons */
    for (hl0 = hl; hl0; hl0 = hl0->next) {
      assert(hl0->rname == rname);

      G = hl0->G;

      assert(G->rb->t == rc);

      /* Consider submatrices */
      if (G->son && G->son[0]->rb->t != rc) {
	assert(G->rsons == rc->sons);

	rsons = G->rsons;
	csons = G->csons;

	mname1 = hl0->mname + 1;
	cname1 = (G->son[0]->cb == G->cb ? hl0->cname : hl0->cname + 1);
	for (j = 0; j < csons; j++) {
	  rname1 = rname + 1;

	  for (i = 0; i < rsons; i++) {
	    assert(G->son[i + j * rsons]->rb->t == rc->son[i]);

	    hl1[i] = addrow_iterate(G->son[i + j * rsons], mname1, rname1,
				    cname1, hl0, hl1[i]);

	    mname1 += G->son[i + j * rsons]->desc;
	    rname1 += rc->son[i]->desc;
	  }
	  assert(rname1 == rname + rc->desc);

	  cname1 += G->son[j * rsons]->cb->t->desc;
	}
	assert(cname1 == hl0->cname + G->cb->t->desc);
      }
    }

    /* Recursively handle sons */
#ifdef USE_OPENMP
    nthreads = rc->sons;	/* HACK: Solaris workaround */
    (void) nthreads;
#pragma omp parallel for if(pardepth > 0), num_threads(nthreads)
#endif
    for (i = 0; i < rc->sons; i++)
      iterate_rowlist(rc->son[i], rnames[i], hl1[i],
		      (pardepth > 0 ? pardepth - 1 : 0), pre, post, data);

    /* Clean up */
    for (i = 0; i < rc->sons; i++)
      del_h2matrixlist(hl1[i]);
    freemem(hl1);
    freemem(rnames);
  }

  /* Use "post" callback, if one is given */
  if (post)
    post(rc, rname, pardepth, hl, data);
}

void
iterate_rowlist_h2matrix(ph2matrix G, uint mname, uint rname, uint cname,
			 uint pardepth,
			 void (*pre) (pccluster t, uint tname, uint pardepth,
				      pch2matrixlist hl, void *data),
			 void (*post) (pccluster t, uint tname, uint pardepth,
				       pch2matrixlist hl, void *data),
			 void *data)
{
  ph2matrixlist hl;

  hl = addrow_iterate(G, mname, rname, cname, 0, 0);

  iterate_rowlist(G->rb->t, rname, hl, pardepth, pre, post, data);

  del_h2matrixlist(hl);
}

static    ph2matrixlist
addcol_iterate(ph2matrix G, uint mname, uint rname,
	       uint cname, pch2matrixlist father, ph2matrixlist next)
{
  ph2matrixlist hl;
  pch2matrixlist hlf;
  pccluster rc;
  uint      mname1, rname1;
  uint      i;

  hl = (ph2matrixlist) allocmem(sizeof(h2matrixlist));
  hl->G = G;
  hl->mname = mname;
  hl->rname = rname;
  hl->cname = cname;
  hl->father = father;
  hl->next = next;

  hlf = hl;

  if (G->son && G->son[0]->cb == G->cb) {
    rc = G->rb->t;

    assert(G->rsons == rc->sons);

    mname1 = mname + 1;
    rname1 = rname + 1;

    for (i = 0; i < G->rsons; i++) {
      hl = addcol_iterate(G->son[i], mname1, rname1, cname, hlf, hl);

      mname1 += G->son[i]->desc;
      rname1 += rc->son[i]->desc;
    }
    assert(mname1 == mname + G->desc);
    assert(rname1 == rname + rc->desc);
  }

  return hl;
}

static void
iterate_collist(pccluster cc, uint cname, ph2matrixlist hl,
		uint pardepth,
		void (*pre) (pccluster t, uint tname, uint pardepth,
			     pch2matrixlist hl, void *data),
		void (*post) (pccluster t, uint tname, uint pardepth,
			      pch2matrixlist hl, void *data), void *data)
{
  ph2matrixlist hl0, *hl1;
  pch2matrix G;
  uint      rsons, csons;
  uint      mname1, rname1, cname1;
  uint     *cnames;
#ifdef USE_OPENMP
  uint      nthreads;		/* HACK: Solaris workaround */
#endif
  uint      i, j;

  /* Quick exit */
  if (hl == 0)
    return;

  /* Use "pre" callback, if one is given */
  if (pre)
    pre(cc, cname, pardepth, hl, data);

  /* Consider sons */
  if (cc->son) {
    /* Initialize matrix lists for all sons */
    hl1 = (ph2matrixlist *) allocmem((size_t) sizeof(ph2matrixlist) *
				     cc->sons);
    for (j = 0; j < cc->sons; j++)
      hl1[j] = 0;

    /* Initialize names for all sons */
    cnames = (uint *) allocmem((size_t) sizeof(uint) * cc->sons);
    cname1 = cname + 1;
    for (j = 0; j < cc->sons; j++) {
      cnames[j] = cname1;
      cname1 += cc->son[j]->desc;
    }
    assert(cname1 == cname + cc->desc);

    /* Fill matrix lists for all sons */
    for (hl0 = hl; hl0; hl0 = hl0->next) {
      assert(hl0->cname == cname);

      G = hl0->G;

      assert(G->cb->t == cc);

      /* Consider submatrices */
      if (G->son && G->son[0]->cb->t != cc) {
	assert(G->csons == cc->sons);

	rsons = G->rsons;
	csons = G->csons;

	mname1 = hl0->mname + 1;
	cname1 = cname + 1;
	for (j = 0; j < csons; j++) {
	  rname1 = (G->son[0]->rb == G->rb ? hl0->rname : hl0->rname + 1);

	  for (i = 0; i < rsons; i++) {
	    assert(G->son[i + j * rsons]->cb->t == cc->son[j]);

	    hl1[j] = addcol_iterate(G->son[i + j * rsons], mname1, rname1,
				    cname1, hl0, hl1[j]);

	    mname1 += G->son[i + j * rsons]->desc;
	    rname1 += G->son[i + j * rsons]->rb->t->desc;
	  }
	  assert(rname1 == hl0->rname + G->rb->t->desc);

	  cname1 += G->son[j * rsons]->cb->t->desc;
	}
	assert(mname1 == hl0->mname + G->desc);
	assert(cname1 == cname + G->cb->t->desc);
      }
    }

    /* Recursively handle sons */
#ifdef USE_OPENMP
    nthreads = cc->sons;
    (void) nthreads;
#pragma omp parallel for if(pardepth > 0), num_threads(nthreads)
#endif
    for (j = 0; j < cc->sons; j++)
      iterate_collist(cc->son[j], cnames[j], hl1[j],
		      (pardepth > 0 ? pardepth - 1 : 0), pre, post, data);

    /* Clean up */
    for (j = 0; j < cc->sons; j++)
      del_h2matrixlist(hl1[j]);
    freemem(hl1);
    freemem(cnames);
  }

  /* Use "post" callback, if one is given */
  if (post)
    post(cc, cname, pardepth, hl, data);
}

void
iterate_collist_h2matrix(ph2matrix G, uint mname, uint rname, uint cname,
			 uint pardepth,
			 void (*pre) (pccluster t, uint tname, uint pardepth,
				      pch2matrixlist hl, void *data),
			 void (*post) (pccluster t, uint tname, uint pardepth,
				       pch2matrixlist hl, void *data),
			 void *data)
{
  ph2matrixlist hl;

  hl = addcol_iterate(G, mname, rname, cname, 0, 0);

  iterate_collist(G->cb->t, cname, hl, pardepth, pre, post, data);

  del_h2matrixlist(hl);
}

struct _listdata {
  void      (*pre) (ph2matrix G, uint mname, uint rname, uint cname,
		    uint pardepth, void *data);
  void      (*post) (ph2matrix G, uint mname, uint rname, uint cname,
		     uint pardepth, void *data);
  void     *data;
};

static void
call_reversed_h2matrixlist(pch2matrixlist hl, uint pardepth,
			   void (*pre) (ph2matrix G, uint mname, uint rname,
					uint cname, uint pardepth,
					void *data), void *data)
{
  if (hl->next)
    call_reversed_h2matrixlist(hl->next, pardepth, pre, data);

  pre(hl->G, hl->mname, hl->rname, hl->cname, pardepth, data);
}

static void
pre_h2matrixlist(pccluster t, uint tname, uint pardepth,
		 pch2matrixlist hl, void *data)
{
  struct _listdata *ld = (struct _listdata *) data;

  (void) t;
  (void) tname;

  if (ld->pre)
    call_reversed_h2matrixlist(hl, pardepth, ld->pre, ld->data);
}

static void
call_inorder_h2matrixlist(pch2matrixlist hl, uint pardepth,
			  void (*post) (ph2matrix G, uint mname, uint rname,
					uint cname, uint pardepth,
					void *data), void *data)
{
  post(hl->G, hl->mname, hl->rname, hl->cname, pardepth, data);

  if (hl->next)
    call_inorder_h2matrixlist(hl->next, pardepth, post, data);
}

static void
post_h2matrixlist(pccluster t, uint tname, uint pardepth,
		  pch2matrixlist hl, void *data)
{
  struct _listdata *ld = (struct _listdata *) data;

  (void) t;
  (void) tname;

  if (ld->post)
    call_inorder_h2matrixlist(hl, pardepth, ld->post, ld->data);
}

void
iterate_byrow_h2matrix(ph2matrix G, uint mname, uint rname, uint cname,
		       uint pardepth,
		       void (*pre) (ph2matrix G, uint mname, uint rname,
				    uint cname, uint pardepth, void *data),
		       void (*post) (ph2matrix G, uint mname, uint rname,
				     uint cname, uint pardepth, void *data),
		       void *data)
{
  struct _listdata ld;
  ph2matrixlist hl;

  ld.pre = pre;
  ld.post = post;
  ld.data = data;

  hl = addrow_iterate(G, mname, rname, cname, 0, 0);

  iterate_rowlist(G->rb->t, rname, hl, pardepth, pre_h2matrixlist,
		  post_h2matrixlist, &ld);

  del_h2matrixlist(hl);
}

void
iterate_bycol_h2matrix(ph2matrix G, uint mname, uint rname, uint cname,
		       uint pardepth,
		       void (*pre) (ph2matrix G, uint mname, uint rname,
				    uint cname, uint pardepth, void *data),
		       void (*post) (ph2matrix G, uint mname, uint rname,
				     uint cname, uint pardepth, void *data),
		       void *data)
{
  struct _listdata ld;
  ph2matrixlist hl;

  ld.pre = pre;
  ld.post = post;
  ld.data = data;

  hl = addcol_iterate(G, mname, rname, cname, 0, 0);

  iterate_collist(G->cb->t, cname, hl, pardepth, pre_h2matrixlist,
		  post_h2matrixlist, &ld);

  del_h2matrixlist(hl);
}

/* ------------------------------------------------------------
 * Matrix-vector multiplication
 * ------------------------------------------------------------ */

void
mvm_h2matrix_avector(field alpha, bool h2trans, pch2matrix h2, pcavector x,
		     pavector y)
{
  if (h2trans)
    addevaltrans_h2matrix_avector(alpha, h2, x, y);
  else
    addeval_h2matrix_avector(alpha, h2, x, y);
}

void
fastaddeval_h2matrix_avector(field alpha, pch2matrix h2, pavector xt,
			     pavector yt)
{
  avector   loc1, loc2;
  pavector  xp, yp, xt1, yt1;
  pcclusterbasis rb = h2->rb;
  pcclusterbasis cb = h2->cb;
  uint      rsons = h2->rsons;
  uint      csons = h2->csons;
  uint      xtoff, ytoff;
  uint      i, j;

  if (h2->u) {
    addeval_amatrix_avector(alpha, &h2->u->S, xt, yt);
  }
  else if (h2->f) {
    xp = init_sub_avector(&loc1, xt, cb->t->size, cb->k);
    yp = init_sub_avector(&loc2, yt, rb->t->size, rb->k);

    addeval_amatrix_avector(alpha, h2->f, xp, yp);

    uninit_avector(yp);
    uninit_avector(xp);
  }
  else if (h2->son) {
    xtoff = cb->k;
    for (j = 0; j < csons; j++) {
      assert(csons == 1 || cb->sons > 0);
      xt1 =
	(cb->sons >
	 0 ? init_sub_avector(&loc1, xt, cb->son[j]->ktree,
			      xtoff) : init_sub_avector(&loc1, xt, cb->ktree,
							0));

      ytoff = rb->k;
      for (i = 0; i < rsons; i++) {
	assert(rsons == 1 || rb->sons > 0);
	yt1 = (rb->sons > 0 ?
	       init_sub_avector(&loc2, yt, rb->son[i]->ktree, ytoff) :
	       init_sub_avector(&loc2, yt, rb->ktree, 0));

	fastaddeval_h2matrix_avector(alpha, h2->son[i + j * rsons], xt1, yt1);

	uninit_avector(yt1);

	ytoff += (rb->sons > 0 ? rb->son[i]->ktree : rb->t->size);
      }
      assert(ytoff == rb->ktree);

      uninit_avector(xt1);

      xtoff += (cb->sons > 0 ? cb->son[j]->ktree : cb->t->size);
    }
    assert(xtoff == cb->ktree);
  }
}

void
addeval_h2matrix_avector(field alpha, pch2matrix h2, pcavector x, pavector y)
{
  pavector  xt, yt;

  xt = new_coeffs_clusterbasis_avector(h2->cb);
  yt = new_coeffs_clusterbasis_avector(h2->rb);

  clear_avector(yt);

  forward_clusterbasis_avector(h2->cb, x, xt);

  fastaddeval_h2matrix_avector(alpha, h2, xt, yt);

  backward_clusterbasis_avector(h2->rb, yt, y);

  del_avector(yt);
  del_avector(xt);
}

void
fastaddevaltrans_h2matrix_avector(field alpha, pch2matrix h2, pavector xt,
				  pavector yt)
{
  avector   loc1, loc2;
  pavector  xp, yp, xt1, yt1;
  pcclusterbasis rb = h2->rb;
  pcclusterbasis cb = h2->cb;
  uint      rsons = h2->rsons;
  uint      csons = h2->csons;
  uint      xtoff, ytoff;
  uint      i, j;

  if (h2->u) {
    addevaltrans_amatrix_avector(alpha, &h2->u->S, xt, yt);
  }
  else if (h2->f) {
    xp = init_sub_avector(&loc1, xt, rb->t->size, rb->k);
    yp = init_sub_avector(&loc2, yt, cb->t->size, cb->k);

    addevaltrans_amatrix_avector(alpha, h2->f, xp, yp);

    uninit_avector(yp);
    uninit_avector(xp);
  }
  else if (h2->son) {
    ytoff = cb->k;
    for (j = 0; j < csons; j++) {
      yt1 =
	(cb->sons >
	 0 ? init_sub_avector(&loc2, yt, cb->son[j]->ktree,
			      ytoff) : init_sub_avector(&loc2, yt, cb->ktree,
							0));

      xtoff = rb->k;
      for (i = 0; i < rsons; i++) {
	xt1 = (rb->sons > 0 ?
	       init_sub_avector(&loc1, xt, rb->son[i]->ktree, xtoff) :
	       init_sub_avector(&loc1, xt, rb->ktree, 0));

	fastaddevaltrans_h2matrix_avector(alpha, h2->son[i + j * rsons], xt1,
					  yt1);

	uninit_avector(xt1);

	xtoff += (rb->sons > 0 ? rb->son[i]->ktree : rb->t->size);
      }
      assert(xtoff == rb->ktree);

      uninit_avector(yt1);

      ytoff += (cb->sons > 0 ? cb->son[j]->ktree : cb->t->size);
    }
    assert(ytoff == cb->ktree);
  }
}

void
addevaltrans_h2matrix_avector(field alpha, pch2matrix h2, pcavector x,
			      pavector y)
{
  pavector  xt, yt;

  xt = new_coeffs_clusterbasis_avector(h2->rb);
  yt = new_coeffs_clusterbasis_avector(h2->cb);

  clear_avector(yt);

  forward_clusterbasis_avector(h2->rb, x, xt);

  fastaddevaltrans_h2matrix_avector(alpha, h2, xt, yt);

  backward_clusterbasis_avector(h2->cb, yt, y);

  del_avector(yt);
  del_avector(xt);
}

static void
addevalsymm_offdiag(field alpha, pch2matrix h2, pavector xt,
		    pavector xta, pavector yt, pavector yta)
{
  avector   tmp1, tmp2, tmp3, tmp4;
  pavector  xp, yp;
  pavector  xt1, xta1, yt1, yta1;
  pcclusterbasis rb = h2->rb;
  pcclusterbasis cb = h2->cb;
  uint      rsons, csons;
  uint      xtoff, ytoff;
  uint      i, j;

  assert(xt->dim == h2->cb->ktree);
  assert(xta->dim == h2->rb->ktree);
  assert(yt->dim == h2->rb->ktree);
  assert(yta->dim == h2->cb->ktree);

  if (h2->f) {
    xp = init_sub_avector(&tmp1, xt, cb->t->size, cb->k);
    yp = init_sub_avector(&tmp2, yt, rb->t->size, rb->k);

    addeval_amatrix_avector(alpha, h2->f, xp, yp);

    uninit_avector(yp);
    uninit_avector(xp);

    xp = init_sub_avector(&tmp1, xta, rb->t->size, rb->k);
    yp = init_sub_avector(&tmp2, yta, cb->t->size, cb->k);

    addevaltrans_amatrix_avector(alpha, h2->f, xp, yp);

    uninit_avector(yp);
    uninit_avector(xp);
  }
  else if (h2->u) {
    addeval_amatrix_avector(alpha, &h2->u->S, xt, yt);
    addevaltrans_amatrix_avector(alpha, &h2->u->S, xta, yta);
  }
  else {
    assert(h2->son != 0);

    rsons = h2->rsons;
    csons = h2->csons;

    xtoff = cb->k;
    for (j = 0; j < csons; j++) {
      xt1 = 0;
      yta1 = 0;
      if (h2->son[j * rsons]->cb == cb) {
	xt1 = init_sub_avector(&tmp1, xt, cb->ktree, 0);
	yta1 = init_sub_avector(&tmp2, yta, cb->ktree, 0);
	xtoff += cb->t->size;
      }
      else {
	assert(j < cb->sons);

	xt1 = init_sub_avector(&tmp1, xt, cb->son[j]->ktree, xtoff);
	yta1 = init_sub_avector(&tmp2, yta, cb->son[j]->ktree, xtoff);
	xtoff += cb->son[j]->ktree;
      }

      ytoff = rb->k;
      for (i = 0; i < rsons; i++) {
	yt1 = 0;
	xta1 = 0;
	if (h2->son[i]->rb == rb) {
	  yt1 = init_sub_avector(&tmp3, yt, rb->ktree, 0);
	  xta1 = init_sub_avector(&tmp4, xta, rb->ktree, 0);
	  ytoff += rb->t->size;
	}
	else {
	  assert(i < rb->sons);

	  yt1 = init_sub_avector(&tmp3, yt, rb->son[i]->ktree, ytoff);
	  xta1 = init_sub_avector(&tmp4, xta, rb->son[i]->ktree, ytoff);
	  ytoff += rb->son[i]->ktree;
	}

	addevalsymm_offdiag(alpha, h2->son[i + j * rsons], xt1, xta1, yt1,
			    yta1);

	uninit_avector(xta1);
	uninit_avector(yt1);
      }
      assert(ytoff == rb->ktree);

      uninit_avector(yta1);
      uninit_avector(xt1);
    }
    assert(xtoff == cb->ktree);
  }
}

static void
addevalsymm_diag(field alpha, pch2matrix h2, pavector xt,
		 pavector xta, pavector yt, pavector yta)
{
  avector   tmp1, tmp2, tmp3, tmp4;
  pavector  xt1, xta1, yt1, yta1;
  pavector  xp, yp;
  pcclusterbasis rb = h2->rb;
  pcclusterbasis cb = h2->cb;
  pfield    aa;
  uint      lda, sons;
  uint      xtoff, ytoff;
  uint      n;
  uint      i, j;

  assert(h2->rb->t == h2->cb->t);
  assert(xt->dim == h2->cb->ktree);
  assert(xta->dim == h2->rb->ktree);
  assert(yt->dim == h2->rb->ktree);
  assert(yta->dim == h2->cb->ktree);

  if (h2->f) {
    aa = h2->f->a;
    lda = h2->f->ld;

    n = rb->t->size;
    xp = init_sub_avector(&tmp1, xt, n, cb->k);
    yp = init_sub_avector(&tmp2, yt, n, rb->k);

    for (j = 0; j < n; j++) {
      yp->v[j] += alpha * aa[j + j * lda] * xp->v[j];
      for (i = j + 1; i < n; i++) {
	yp->v[i] += alpha * aa[i + j * lda] * xp->v[j];
	yp->v[j] += alpha * CONJ(aa[i + j * lda]) * xp->v[i];
      }
    }

    uninit_avector(yp);
    uninit_avector(xp);
  }
  else {
    assert(h2->son != 0);
    assert(h2->rsons == h2->csons);

    sons = h2->rsons;

    xtoff = cb->k;
    for (j = 0; j < sons; j++) {
      xt1 = 0;
      yta1 = 0;
      if (h2->son[j * sons]->cb == cb) {
	xt1 = init_sub_avector(&tmp1, xt, cb->ktree, 0);
	yta1 = init_sub_avector(&tmp2, yta, cb->ktree, 0);
	xtoff += cb->t->size;
      }
      else {
	assert(j < cb->sons);

	xt1 = init_sub_avector(&tmp1, xt, cb->son[j]->ktree, xtoff);
	yta1 = init_sub_avector(&tmp2, yta, cb->son[j]->ktree, xtoff);
	xtoff += cb->son[j]->ktree;
      }

      ytoff = rb->k;

      yt1 = 0;
      xta1 = 0;
      if (h2->son[0]->rb == rb) {
	yt1 = init_sub_avector(&tmp3, yt, rb->ktree, 0);
	xta1 = init_sub_avector(&tmp4, xta, rb->ktree, 0);
	ytoff += rb->t->size;
      }
      else {
	/* Skip superdiagonal blocks */
	for (i = 0; i < j; i++)
	  ytoff += rb->son[i]->ktree;

	/* Subvectors for diagonal block */
	yt1 = init_sub_avector(&tmp3, yt, rb->son[i]->ktree, ytoff);
	xta1 = init_sub_avector(&tmp4, xta, rb->son[i]->ktree, ytoff);
	ytoff += rb->son[i]->ktree;
      }

      addevalsymm_diag(alpha, h2->son[j + j * sons], xt1, xta1, yt1, yta1);

      uninit_avector(xta1);
      uninit_avector(yt1);

      for (i = j + 1; i < sons; i++) {
	assert(i < rb->sons);

	yt1 = init_sub_avector(&tmp3, yt, rb->son[i]->ktree, ytoff);
	xta1 = init_sub_avector(&tmp4, xta, rb->son[i]->ktree, ytoff);
	ytoff += rb->son[i]->ktree;

	addevalsymm_offdiag(alpha, h2->son[i + j * sons], xt1, xta1, yt1,
			    yta1);

	uninit_avector(xta1);
	uninit_avector(yt1);
      }
      assert(ytoff == rb->ktree);

      uninit_avector(yta1);
      uninit_avector(xt1);
    }
    assert(xtoff == cb->ktree);
  }
}

void
addevalsymm_h2matrix_avector(field alpha, pch2matrix h2, pcavector x,
			     pavector y)
{
  pavector  xt, yt, xta, yta;

  assert(h2->rb->t == h2->cb->t);

  /* Transformed coefficients */
  xt = new_coeffs_clusterbasis_avector(h2->cb);
  xta = new_coeffs_clusterbasis_avector(h2->rb);
  yt = new_coeffs_clusterbasis_avector(h2->rb);
  yta = new_coeffs_clusterbasis_avector(h2->cb);

  /* Clear row coefficients */
  clear_avector(yt);
  clear_avector(yta);

  /* Column coefficients filled by forward transformation */
  forward_clusterbasis_avector(h2->cb, x, xt);
  forward_clusterbasis_avector(h2->rb, x, xta);

  /* Multiplication step */
  addevalsymm_diag(alpha, h2, xt, xta, yt, yta);

  /* Row coefficients added to result by backward transformation */
  backward_clusterbasis_avector(h2->rb, yt, y);
  backward_clusterbasis_avector(h2->cb, yta, y);

  /* Clean up */
  uninit_avector(yta);
  uninit_avector(yt);
  uninit_avector(xta);
  uninit_avector(xt);
}

/* ------------------------------------------------------------
 * Addmul H2-Matrices and Amatrix
 * ------------------------------------------------------------ */

void
fastaddmul_h2matrix_amatrix_amatrix(field alpha, bool h2trans,
				    pch2matrix h2, pcamatrix Xt, pamatrix Yt)
{
  amatrix   loc1, loc2;
  pamatrix  Xp, Yp, Xt1, Yt1;
  pcclusterbasis rb = (h2trans ? h2->cb : h2->rb);
  pcclusterbasis cb = (h2trans ? h2->rb : h2->cb);
  uint      rsons = (h2trans ? h2->csons : h2->rsons);
  uint      csons = (h2trans ? h2->rsons : h2->csons);
  uint      cols = Xt->cols;
  uint      xtoff, ytoff;
  uint      i, j, k;

  if (h2->u) {
    Xt1 = init_sub_amatrix(&loc1, (pamatrix) Xt, cb->k, 0, Xt->cols, 0);
    Yt1 = init_sub_amatrix(&loc2, Yt, rb->k, 0, Yt->cols, 0);
    addmul_amatrix(alpha, h2trans, &h2->u->S, false, Xt1, Yt1);
    uninit_amatrix(Yt1);
    uninit_amatrix(Xt1);
  }
  else if (h2->f) {
    Xp = init_sub_amatrix(&loc1, (pamatrix) Xt, cb->t->size, cb->k, Xt->cols,
			  0);
    Yp = init_sub_amatrix(&loc2, Yt, rb->t->size, rb->k, Yt->cols, 0);
    addmul_amatrix(alpha, h2trans, h2->f, false, Xp, Yp);
    uninit_amatrix(Yp);
    uninit_amatrix(Xp);
  }
  else if (h2->son) {
    xtoff = cb->k;
    for (j = 0; j < csons; j++) {
      assert(csons == 1 || cb->sons > 0);
      Xt1 = (cb->sons > 0 ?
	     init_sub_amatrix(&loc1, (pamatrix) Xt, cb->son[j]->ktree, xtoff,
			      cols, 0) :
	     init_sub_amatrix(&loc1, (pamatrix) Xt, cb->ktree, 0, cols, 0));

      ytoff = rb->k;
      for (i = 0; i < rsons; i++) {
	assert(rsons == 1 || rb->sons > 0);
	Yt1 = (rb->sons > 0 ?
	       init_sub_amatrix(&loc2, Yt, rb->son[i]->ktree, ytoff, cols,
				0) : init_sub_amatrix(&loc2, Yt, rb->ktree, 0,
						      cols, 0));
	k = (h2trans ? i * csons + j : i + j * rsons);
	fastaddmul_h2matrix_amatrix_amatrix(alpha, h2trans, h2->son[k], Xt1,
					    Yt1);
	uninit_amatrix(Yt1);

	ytoff += (rb->sons > 0 ? rb->son[i]->ktree : rb->t->size);
      }
      assert(ytoff == rb->ktree);
      uninit_amatrix(Xt1);

      xtoff += (cb->sons > 0 ? cb->son[j]->ktree : cb->t->size);
    }
    assert(xtoff == cb->ktree);
  }
}

void
addmul_h2matrix_amatrix_amatrix(field alpha, bool h2trans, pch2matrix h2,
				bool xtrans, pcamatrix X, pamatrix Y)
{
  pamatrix  Xt, Yt;

  if (!h2trans) {
    if (!xtrans) {
      assert(X->cols == Y->cols);
      assert(X->rows == h2->cb->t->size);
      assert(Y->rows == h2->rb->t->size);

      Xt = new_amatrix(h2->cb->ktree, X->cols);
      Yt = new_amatrix(h2->rb->ktree, Y->cols);
      clear_amatrix(Yt);
      forward_clusterbasis_amatrix(h2->cb, X, Xt);
      fastaddmul_h2matrix_amatrix_amatrix(alpha, false, h2, Xt, Yt);
      backward_clusterbasis_amatrix(h2->rb, Yt, Y);
      del_amatrix(Yt);
      del_amatrix(Xt);
    }
    else {
      assert(xtrans);
      assert(X->rows == Y->cols);
      assert(X->cols == h2->cb->t->size);
      assert(Y->rows == h2->rb->t->size);

      Xt = new_amatrix(h2->cb->ktree, X->rows);
      Yt = new_amatrix(h2->rb->ktree, Y->cols);
      clear_amatrix(Yt);
      forward_clusterbasis_trans_amatrix(h2->cb, X, Xt);
      fastaddmul_h2matrix_amatrix_amatrix(alpha, false, h2, Xt, Yt);
      backward_clusterbasis_amatrix(h2->rb, Yt, Y);
      del_amatrix(Yt);
      del_amatrix(Xt);
    }
  }
  else {
    assert(h2trans);
    if (!xtrans) {
      assert(X->cols == Y->cols);
      assert(X->rows == h2->rb->t->size);
      assert(Y->rows == h2->cb->t->size);

      Xt = new_amatrix(h2->rb->ktree, X->cols);
      Yt = new_amatrix(h2->cb->ktree, Y->cols);
      clear_amatrix(Yt);
      forward_clusterbasis_amatrix(h2->rb, X, Xt);
      fastaddmul_h2matrix_amatrix_amatrix(alpha, true, h2, Xt, Yt);
      backward_clusterbasis_amatrix(h2->cb, Yt, Y);
      del_amatrix(Yt);
      del_amatrix(Xt);
    }
    else {
      assert(xtrans);
      assert(X->rows == Y->cols);
      assert(X->cols == h2->rb->t->size);
      assert(Y->rows == h2->cb->t->size);

      Xt = new_amatrix(h2->rb->ktree, X->rows);
      Yt = new_amatrix(h2->cb->ktree, Y->cols);
      clear_amatrix(Yt);
      forward_clusterbasis_trans_amatrix(h2->rb, X, Xt);
      fastaddmul_h2matrix_amatrix_amatrix(alpha, true, h2, Xt, Yt);
      backward_clusterbasis_amatrix(h2->cb, Yt, Y);
      del_amatrix(Yt);
      del_amatrix(Xt);
    }
  }
}

void
addmul_amatrix_h2matrix_amatrix(field alpha, bool xtrans, pcamatrix X,
				bool h2trans, pch2matrix h2, pamatrix Y)
{
  pamatrix  Xt, Yt;

  if (h2trans) {
    if (xtrans) {
      assert(X->cols == Y->rows);
      assert(X->rows == h2->cb->t->size);
      assert(Y->cols == h2->rb->t->size);

      Xt = new_amatrix(h2->cb->ktree, X->cols);
      Yt = new_amatrix(h2->rb->ktree, Y->rows);
      clear_amatrix(Yt);
      forward_clusterbasis_amatrix(h2->cb, X, Xt);
      fastaddmul_h2matrix_amatrix_amatrix(alpha, false, h2, Xt, Yt);
      backward_clusterbasis_amatrix(h2->rb, Yt, Y);
      del_amatrix(Yt);
      del_amatrix(Xt);
    }
    else {
      assert(!xtrans);
      assert(X->rows == Y->rows);
      assert(X->cols == h2->cb->t->size);
      assert(Y->cols == h2->rb->t->size);

      Xt = new_amatrix(h2->cb->ktree, X->rows);
      Yt = new_amatrix(h2->rb->ktree, Y->rows);
      clear_amatrix(Yt);
      forward_clusterbasis_trans_amatrix(h2->cb, X, Xt);
      fastaddmul_h2matrix_amatrix_amatrix(alpha, false, h2, Xt, Yt);
      backward_clusterbasis_trans_amatrix(h2->rb, Yt, Y);
      del_amatrix(Yt);
      del_amatrix(Xt);
    }
  }
  else {
    assert(!h2trans);
    if (xtrans) {
      assert(X->cols == Y->rows);
      assert(X->rows == h2->rb->t->size);
      assert(Y->cols == h2->cb->t->size);

      Xt = new_amatrix(h2->rb->ktree, X->cols);
      Yt = new_amatrix(h2->cb->ktree, Y->rows);
      clear_amatrix(Yt);
      forward_clusterbasis_amatrix(h2->rb, X, Xt);
      fastaddmul_h2matrix_amatrix_amatrix(alpha, true, h2, Xt, Yt);
      backward_clusterbasis_trans_amatrix(h2->cb, Yt, Y);
      del_amatrix(Yt);
      del_amatrix(Xt);
    }
    else {
      assert(!xtrans);
      assert(X->rows == Y->rows);
      assert(X->cols == h2->rb->t->size);
      assert(Y->cols == h2->cb->t->size);

      Xt = new_amatrix(h2->rb->ktree, X->rows);
      Yt = new_amatrix(h2->cb->ktree, Y->rows);
      clear_amatrix(Yt);
      forward_clusterbasis_trans_amatrix(h2->rb, X, Xt);
      fastaddmul_h2matrix_amatrix_amatrix(alpha, true, h2, Xt, Yt);
      backward_clusterbasis_trans_amatrix(h2->cb, Yt, Y);
      del_amatrix(Yt);
      del_amatrix(Xt);
    }
  }
}

/* ------------------------------------------------------------
 * Orthogonal projection
 * ------------------------------------------------------------ */

void
collectdense_h2matrix(pcamatrix a, pcclusterbasis rb, pcclusterbasis cb,
		      pamatrix s)
{
  amatrix   tmp1, tmp2;
  pamatrix  s1, s2;
  pccluster row, col;
  uint      i, j;

  assert(s->rows == rb->k);
  assert(s->cols == cb->k);

  clear_amatrix(s);

  if (rb->sons > 0) {
    if (cb->sons > 0) {		/* rb has sons, cb has sons */
      for (j = 0; j < cb->sons; j++) {
	s1 = init_amatrix(&tmp1, rb->k, cb->son[j]->k);
	clear_amatrix(s1);

	for (i = 0; i < rb->sons; i++) {
	  s2 = init_amatrix(&tmp2, rb->son[i]->k, cb->son[j]->k);

	  collectdense_h2matrix(a, rb->son[i], cb->son[j], s2);

	  addmul_amatrix(1.0, true, &rb->son[i]->E, false, s2, s1);

	  uninit_amatrix(s2);
	}
	addmul_amatrix(1.0, false, s1, false, &cb->son[j]->E, s);

	uninit_amatrix(s1);
      }
    }
    else {			/* rb has sons, cb is a leaf */
      for (i = 0; i < rb->sons; i++) {
	s2 = init_amatrix(&tmp2, rb->son[i]->k, cb->k);

	collectdense_h2matrix(a, rb->son[i], cb, s2);

	addmul_amatrix(1.0, true, &rb->son[i]->E, false, s2, s);

	uninit_amatrix(s2);
      }
    }
  }
  else {
    if (cb->sons > 0) {		/* rb is a leaf, cb has sons */
      for (j = 0; j < cb->sons; j++) {
	s1 = init_amatrix(&tmp1, rb->k, cb->son[j]->k);

	collectdense_h2matrix(a, rb, cb->son[j], s1);

	addmul_amatrix(1.0, false, s1, false, &cb->son[j]->E, s);

	uninit_amatrix(s1);
      }
    }
    else {			/* rb is a leaf, cb is a leaf */
      row = rb->t;
      col = cb->t;

      s1 = init_amatrix(&tmp1, rb->k, col->size);
      clear_amatrix(s1);

      s2 = init_amatrix(&tmp2, row->size, col->size);

      for (j = 0; j < col->size; j++)
	for (i = 0; i < row->size; i++)
	  s2->a[i + j * s2->ld] = a->a[row->idx[i] + col->idx[j] * a->ld];

      addmul_amatrix(1.0, true, &rb->V, false, s2, s1);

      addmul_amatrix(1.0, false, s1, false, &cb->V, s);

      uninit_amatrix(s2);
      uninit_amatrix(s1);
    }
  }
}

void
project_amatrix_h2matrix(ph2matrix h2, pcamatrix a)
{
  pccluster row, col;
  uint      i, j;

  if (h2->son) {
    for (j = 0; j < h2->csons; j++)
      for (i = 0; i < h2->rsons; i++)
	project_amatrix_h2matrix(h2->son[i + j * h2->rsons], a);
  }
  else if (h2->f) {
    row = h2->rb->t;
    col = h2->cb->t;

    assert(h2->f->rows == row->size);
    assert(h2->f->cols == col->size);

    for (j = 0; j < col->size; j++)
      for (i = 0; i < row->size; i++)
	setentry_amatrix(h2->f, i, j,
			 getentry_amatrix(a, row->idx[i], col->idx[j]));
  }
  else
    collectdense_h2matrix(a, h2->rb, h2->cb, &h2->u->S);
}

void
project_hmatrix_h2matrix(ph2matrix h2, phmatrix h)
{
  uint      i;
  amatrix   loca, locb;
  pamatrix  A, B;

  if (h2->son) {
    assert(h->son);
    for (i = 0; i < h2->rsons * h2->csons; i++) {
      project_hmatrix_h2matrix(h2->son[i], h->son[i]);
    }
  }
  else if (h2->f) {
    assert(h->f);
    copy_amatrix(false, h->f, h2->f);
  }
  else {
    assert(h2->u);
    assert(h->r);

    A = init_amatrix(&loca, h2->rb->kbranch, h->r->k);
    compress_clusterbasis_amatrix(h2->rb, &h->r->A, A);
    A->rows = h2->rb->k;
    B = init_amatrix(&locb, h2->cb->kbranch, h->r->k);
    compress_clusterbasis_amatrix(h2->cb, &h->r->B, B);
    B->rows = h2->cb->k;

    clear_amatrix(&h2->u->S);
    addmul_amatrix(1.0, false, A, true, B, &h2->u->S);

    uninit_amatrix(B);
    uninit_amatrix(A);
  }
}

/* ------------------------------------------------------------
 * Spectral norm
 * ------------------------------------------------------------ */

real
norm2_h2matrix(pch2matrix H2)
{
  return norm2_matrix((mvm_t) mvm_h2matrix_avector, (void *) H2,
		      H2->rb->t->size, H2->cb->t->size);
}

real
norm2diff_h2matrix(pch2matrix a, pch2matrix b)
{
  return norm2diff_matrix((mvm_t) mvm_h2matrix_avector, (void *) a,
			  (mvm_t) mvm_h2matrix_avector, (void *) b,
			  a->rb->t->size, a->cb->t->size);
}

/* ------------------------------------------------------------
 File I/O
 ------------------------------------------------------------ */

#ifdef USE_NETCDF
static void
write_count(pch2matrix G, size_t * blocks, size_t * coeffs)
{
  uint      rsons, csons;
  uint      i, j;

  /* Increase submatrix counter */
  (*blocks)++;

  if (G->f) {
    /* If it's an amatrix: add number of coefficients to coeffs */
    (*coeffs) += G->f->rows * G->f->cols;
  }
  else if (G->u) {
    /* If it's a uniform matrix: also add number of coefficients to coeffs */
    (*coeffs) += G->u->rb->k * G->u->cb->k;
  }
  else if (G->son) {
    rsons = G->rsons;
    csons = G->csons;

    /* If it's subdivided: check sons */
    for (j = 0; j < csons; j++)
      for (i = 0; i < rsons; i++)
	write_count(G->son[i + j * rsons], blocks, coeffs);
  }
}

static void
write_cdf(pch2matrix G,
	  size_t blocks, size_t coeffs,
	  size_t * blockidx, size_t * coeffidx,
	  int nc_file, int nc_type, int nc_rsons, int nc_csons, int nc_coeff)
{
  size_t    start, count;
  ptrdiff_t stride;
  pcamatrix f;
  pcuniform u;
  uint      rsons, csons;
  uint      rows, cols;
  uint      kr, kc;
  int       val;
  uint      i, j;

  assert(*blockidx <= blocks);

  if (G->f) {
    f = G->f;
    rows = f->rows;
    cols = f->cols;

    /* Write type 1 to nc_type[blockidx] */
    start = *blockidx;
    count = 1;
    stride = 1;
    val = 1;
    nc_put_vars(nc_file, nc_type, &start, &count, &stride, &val);

    /* An amatrix has no rsons or csons */
    val = 0;
    nc_put_vars(nc_file, nc_rsons, &start, &count, &stride, &val);
    nc_put_vars(nc_file, nc_csons, &start, &count, &stride, &val);

    /* Write coefficients */
    start = *coeffidx;
    count = rows;
    assert(start + rows * cols <= coeffs);
    for (j = 0; j < cols; j++) {
      nc_put_vars(nc_file, nc_coeff, &start, &count, &stride,
		  f->a + j * f->ld);
      start += rows;
    }
    *coeffidx = start;

    /* Increase submatrices counter */
    (*blockidx)++;
  }
  else if (G->u) {
    u = G->u;
    kr = u->rb->k;
    kc = u->cb->k;

    /* Write type 2 to nc_type[blockidx] */
    start = *blockidx;
    count = 1;
    stride = 1;
    val = 2;
    nc_put_vars(nc_file, nc_type, &start, &count, &stride, &val);

    /* An rkmatrix has no rsons or csons */
    val = 0;
    nc_put_vars(nc_file, nc_rsons, &start, &count, &stride, &val);
    nc_put_vars(nc_file, nc_csons, &start, &count, &stride, &val);

    /* Write coefficients */
    start = *coeffidx;
    count = kr;
    assert(start + kr * kc <= coeffs);
    for (j = 0; j < kc; j++) {
      nc_put_vars(nc_file, nc_coeff, &start, &count, &stride,
		  u->S.a + j * u->S.ld);
      start += kr;
    }
    (*coeffidx) = start;

    /* Increase submatrix counter */
    (*blockidx)++;
  }
  else if (G->son) {
    rsons = G->rsons;
    csons = G->csons;

    /* Write type 3 to nc_type */
    start = *blockidx;
    count = 1;
    stride = 1;
    val = 3;
    nc_put_vars(nc_file, nc_type, &start, &count, &stride, &val);

    /* Write rsons and csons.
     * rsons == 0 means that the sons use the same row cluster basis
     * as the father. The same convention is used for csons. */
    val = (G->rb == G->son[0]->rb ? 0 : rsons);
    nc_put_vars(nc_file, nc_rsons, &start, &count, &stride, &val);
    val = (G->cb == G->son[0]->cb ? 0 : csons);
    nc_put_vars(nc_file, nc_csons, &start, &count, &stride, &val);

    /* Increase submatrix counter */
    (*blockidx)++;

    /* Take care of sons */
    for (j = 0; j < csons; j++)
      for (i = 0; i < rsons; i++)
	write_cdf(G->son[i + j * rsons], blocks, coeffs,
		  blockidx, coeffidx,
		  nc_file, nc_type, nc_rsons, nc_csons, nc_coeff);
  }
}

void
write_cdf_h2matrix(pch2matrix G, const char *name)
{
  size_t    blocks, blockidx;
  size_t    coeffs, coeffidx;
  int       nc_file, nc_type, nc_rsons, nc_csons, nc_coeff;
  int       nc_blocks, nc_coeffs;
  int       result;

  /* Count number of submatrices and coefficients */
  blocks = 0;
  coeffs = 0;
  write_count(G, &blocks, &coeffs);

  /* Create NetCDF file */
  result = nc_create(name, NC_64BIT_OFFSET, &nc_file);
  assert(result == NC_NOERR);

  /* Define "blocks" dimension */
  result = nc_def_dim(nc_file, "blocks", blocks, &nc_blocks);
  assert(result == NC_NOERR);

  /* Define "coeffs" dimension */
  result = nc_def_dim(nc_file, "coeffs", coeffs, &nc_coeffs);
  assert(result == NC_NOERR);

  /* Define "type" variable */
  result = nc_def_var(nc_file, "type", NC_INT, 1, &nc_blocks, &nc_type);
  assert(result == NC_NOERR);

  /* Define "rsons" variable */
  result = nc_def_var(nc_file, "rsons", NC_INT, 1, &nc_blocks, &nc_rsons);
  assert(result == NC_NOERR);

  /* Define "csons" variable */
  result = nc_def_var(nc_file, "csons", NC_INT, 1, &nc_blocks, &nc_csons);
  assert(result == NC_NOERR);

  /* Define "coeff" variable */
  result = nc_def_var(nc_file, "coeff", NC_DOUBLE, 1, &nc_coeffs, &nc_coeff);
  assert(result == NC_NOERR);

  /* Finish NetCDF define mode */
  result = nc_enddef(nc_file);
  assert(result == NC_NOERR);

  /* Write matrices to NetCDF variables */
  blockidx = 0;
  coeffidx = 0;
  write_cdf(G, blocks, coeffs, &blockidx, &coeffidx,
	    nc_file, nc_type, nc_rsons, nc_csons, nc_coeff);
  assert(blockidx == blocks);
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
write_cdfpart_h2matrix(pch2matrix G, int nc_file, const char *prefix)
{
  size_t    blocks, blockidx;
  size_t    coeffs, coeffidx;
  char     *buf;
  int       bufsize;
  int       nc_type, nc_rsons, nc_csons, nc_coeff;
  int       nc_blocks, nc_coeffs;
  int       result;

  /* Prepare buffer for prefixed names */
  bufsize = strlen(prefix) + 16;
  buf = (char *) allocmem(sizeof(char) * bufsize);

  /* Count number of submatrices and coefficients */
  blocks = 0;
  coeffs = 0;
  write_count(G, &blocks, &coeffs);

  /* Switch NetCDF file to define mode */
  result = nc_redef(nc_file);
  assert(result == NC_NOERR || result == NC_EINDEFINE);

  /* Define "blocks" dimension */
  prefix_name(buf, bufsize, prefix, "blocks");
  result = nc_def_dim(nc_file, buf, blocks, &nc_blocks);
  assert(result == NC_NOERR);

  /* Define "coeffs" dimension */
  prefix_name(buf, bufsize, prefix, "coeffs");
  result = nc_def_dim(nc_file, buf, coeffs, &nc_coeffs);
  assert(result == NC_NOERR);

  /* Define "type" variable */
  prefix_name(buf, bufsize, prefix, "type");
  result = nc_def_var(nc_file, buf, NC_INT, 1, &nc_blocks, &nc_type);
  assert(result == NC_NOERR);

  /* Define "rsons" variable */
  prefix_name(buf, bufsize, prefix, "rsons");
  result = nc_def_var(nc_file, buf, NC_INT, 1, &nc_blocks, &nc_rsons);
  assert(result == NC_NOERR);

  /* Define "csons" variable */
  prefix_name(buf, bufsize, prefix, "csons");
  result = nc_def_var(nc_file, buf, NC_INT, 1, &nc_blocks, &nc_csons);
  assert(result == NC_NOERR);

  /* Define "coeff" variable */
  prefix_name(buf, bufsize, prefix, "coeff");
  result = nc_def_var(nc_file, buf, NC_DOUBLE, 1, &nc_coeffs, &nc_coeff);
  assert(result == NC_NOERR);

  /* Finish NetCDF define mode */
  result = nc_enddef(nc_file);
  assert(result == NC_NOERR);

  /* Write matrices to NetCDF variables */
  blockidx = 0;
  coeffidx = 0;
  write_cdf(G, blocks, coeffs, &blockidx, &coeffidx,
	    nc_file, nc_type, nc_rsons, nc_csons, nc_coeff);
  assert(blockidx == blocks);
  assert(coeffidx == coeffs);

  /* Clean up */
  nc_sync(nc_file);
  freemem(buf);
}

static ph2matrix
read_cdf_part(int nc_file, pclusterbasis rb, pclusterbasis cb,
	      int nc_type, int nc_rsons, int nc_csons, int nc_coeffs,
	      size_t * blockidx, size_t * coeffidx)
{
  ph2matrix G, G1;
  pamatrix  f;
  puniform  u;
  uint      rsons, csons;
  uint      rows, cols;
  uint      kr, kc;
  uint      i, j;
  size_t    start, count;
  ptrdiff_t stride;
  int       val, result;

  /* Get type */
  start = *blockidx;
  count = 1;
  stride = 1;
  result = nc_get_vars(nc_file, nc_type, &start, &count, &stride, &val);
  assert(result == NC_NOERR);

  switch (val) {
  case 1:
    /* Create full matrix */
    G = new_full_h2matrix(rb, cb);
    f = G->f;

    rows = f->rows;
    cols = f->cols;

    /* Read coefficients */
    start = (*coeffidx);
    count = rows;
    stride = 1;
    for (j = 0; j < cols; j++) {
      result =
	nc_get_vars(nc_file, nc_coeffs, &start, &count, &stride,
		    f->a + j * f->ld);
      assert(result == NC_NOERR);

      start += rows;
    }
    (*coeffidx) = start;

    /* Increment block index */
    (*blockidx)++;
    break;

  case 2:
    /* Create uniform matrix */
    G = new_uniform_h2matrix(rb, cb);
    u = G->u;

    kr = u->rb->k;
    kc = u->cb->k;

    /* Read coefficients */
    start = (*coeffidx);
    count = kr;
    stride = 1;
    for (j = 0; j < kc; j++) {
      result =
	nc_get_vars(nc_file, nc_coeffs, &start, &count, &stride,
		    u->S.a + j * u->S.ld);
      assert(result == NC_NOERR);

      start += kr;
    }
    (*coeffidx) = start;

    /* Increment block index */
    (*blockidx)++;
    break;

  case 3:
    /* Get rsons */
    result = nc_get_vars(nc_file, nc_rsons, &start, &count, &stride, &val);
    assert(result == NC_NOERR);
    rsons = val;

    /* Get csons */
    result = nc_get_vars(nc_file, nc_csons, &start, &count, &stride, &val);
    assert(result == NC_NOERR);
    csons = val;

    /* Increment block index */
    (*blockidx)++;

    /* Handle submatrices */
    if (csons == 0) {
      /* Re-using cb is only allowed if cb is a leaf */
      assert(cb->sons == 0);

      if (rsons == 0) {
	/* Re-using rb is only allowed if rb is a leaf */
	assert(rb->sons == 0);

	/* This should never happen */
	fprintf(stderr, "Warning: son of h2matrix with leaf cluster bases\n");

	/* Create subdivided matrix */
	G = new_super_h2matrix(rb, cb, 1, 1);

	/* Read son */
	G1 = read_cdf_part(nc_file, rb, cb,
			   nc_type, nc_rsons, nc_csons, nc_coeffs,
			   blockidx, coeffidx);
	ref_h2matrix(G->son, G1);
      }
      else {
	assert(rb->sons == rsons);

	/* Create subdivided matrix */
	G = new_super_h2matrix(rb, cb, rsons, 1);

	/* Read sons */
	for (i = 0; i < rsons; i++) {
	  G1 = read_cdf_part(nc_file, rb->son[i], cb,
			     nc_type, nc_rsons, nc_csons, nc_coeffs,
			     blockidx, coeffidx);
	  ref_h2matrix(G->son + i, G1);
	}
      }
    }
    else {
      assert(cb->sons == csons);

      if (rsons == 0) {
	/* Re-using rb is only allowed if rb is a leaf */
	assert(rb->sons == 0);

	/* Create subdivided matrix */
	G = new_super_h2matrix(rb, cb, 1, csons);

	/* Read sons */
	for (j = 0; j < csons; j++) {
	  G1 = read_cdf_part(nc_file, rb, cb->son[j],
			     nc_type, nc_rsons, nc_csons, nc_coeffs,
			     blockidx, coeffidx);
	  ref_h2matrix(G->son + j, G1);
	}
      }
      else {
	assert(rb->sons == rsons);

	/* Create subdivided matrix */
	G = new_super_h2matrix(rb, cb, rsons, csons);

	/* Read sons */
	for (j = 0; j < csons; j++)
	  for (i = 0; i < rsons; i++) {
	    G1 = read_cdf_part(nc_file, rb->son[i], cb->son[j],
			       nc_type, nc_rsons, nc_csons, nc_coeffs,
			       blockidx, coeffidx);
	    ref_h2matrix(G->son + i + j * rsons, G1);
	  }
      }
    }
    break;

  default:
    fprintf(stderr, "Unknown block type %d\n", val);
    abort();
  }

  update_h2matrix(G);

  return G;
}

ph2matrix
read_cdf_h2matrix(const char *name, pclusterbasis rb, pclusterbasis cb)
{
  ph2matrix G;
  size_t    blocks, blockidx;
  size_t    coeffs, coeffidx;
  char      dimname[NC_MAX_NAME + 1];
  int       nc_file, nc_type, nc_rsons, nc_csons, nc_coeff;
  int       nc_blocks, nc_coeffs;
  int       result;

  /* Open NetCDF file */
  result = nc_open(name, NC_NOWRITE, &nc_file);
  assert(result == NC_NOERR);

  /* Get "blocks" dimension */
  result = nc_inq_dimid(nc_file, "blocks", &nc_blocks);
  assert(result == NC_NOERR);
  result = nc_inq_dim(nc_file, nc_blocks, dimname, &blocks);
  assert(result == NC_NOERR);

  /* Get "coeffs" dimension */
  result = nc_inq_dimid(nc_file, "coeffs", &nc_coeffs);
  assert(result == NC_NOERR);
  result = nc_inq_dim(nc_file, nc_coeffs, dimname, &coeffs);
  assert(result == NC_NOERR);

  /* Get "type" variable */
  result = nc_inq_varid(nc_file, "type", &nc_type);
  assert(result == NC_NOERR);

  /* Get "rsons" variable */
  result = nc_inq_varid(nc_file, "rsons", &nc_rsons);
  assert(result == NC_NOERR);

  /* Get "csons" variable */
  result = nc_inq_varid(nc_file, "csons", &nc_csons);
  assert(result == NC_NOERR);

  /* Get "coeff" variable */
  result = nc_inq_varid(nc_file, "coeff", &nc_coeff);
  assert(result == NC_NOERR);

  /* Read matrices from NetCDF variables */
  blockidx = 0;
  coeffidx = 0;

  G = read_cdf_part(nc_file, rb, cb,
		    nc_type, nc_rsons, nc_csons, nc_coeff,
		    &blockidx, &coeffidx);
  assert(blockidx == blocks);
  assert(coeffidx == coeffs);

  return G;
}

ph2matrix
read_cdfpart_h2matrix(int nc_file, const char *prefix,
		      pclusterbasis rb, pclusterbasis cb)
{
  ph2matrix G;
  size_t    blocks, blockidx;
  size_t    coeffs, coeffidx;
  char      dimname[NC_MAX_NAME + 1];
  char     *buf;
  int       bufsize;
  int       nc_type, nc_rsons, nc_csons, nc_coeff;
  int       nc_blocks, nc_coeffs;
  int       result;

  /* Prepare buffer for prefixed names */
  bufsize = strlen(prefix) + 16;
  buf = (char *) allocmem(sizeof(char) * bufsize);

  /* Get "blocks" dimension */
  prefix_name(buf, bufsize, prefix, "blocks");
  result = nc_inq_dimid(nc_file, buf, &nc_blocks);
  assert(result == NC_NOERR);
  result = nc_inq_dim(nc_file, nc_blocks, dimname, &blocks);
  assert(result == NC_NOERR);

  /* Get "coeffs" dimension */
  prefix_name(buf, bufsize, prefix, "coeffs");
  result = nc_inq_dimid(nc_file, buf, &nc_coeffs);
  assert(result == NC_NOERR);
  result = nc_inq_dim(nc_file, nc_coeffs, dimname, &coeffs);
  assert(result == NC_NOERR);

  /* Get "type" variable */
  prefix_name(buf, bufsize, prefix, "type");
  result = nc_inq_varid(nc_file, buf, &nc_type);
  assert(result == NC_NOERR);

  /* Get "rsons" variable */
  prefix_name(buf, bufsize, prefix, "rsons");
  result = nc_inq_varid(nc_file, buf, &nc_rsons);
  assert(result == NC_NOERR);

  /* Get "csons" variable */
  prefix_name(buf, bufsize, prefix, "csons");
  result = nc_inq_varid(nc_file, buf, &nc_csons);
  assert(result == NC_NOERR);

  /* Get "coeff" variable */
  prefix_name(buf, bufsize, prefix, "coeff");
  result = nc_inq_varid(nc_file, buf, &nc_coeff);
  assert(result == NC_NOERR);

  /* Read matrices from NetCDF variables */
  blockidx = 0;
  coeffidx = 0;

  G = read_cdf_part(nc_file, rb, cb,
		    nc_type, nc_rsons, nc_csons, nc_coeff,
		    &blockidx, &coeffidx);
  assert(blockidx == blocks);
  assert(coeffidx == coeffs);

  /* Clean up */
  freemem(buf);

  return G;
}

void
write_cdfcomplete_h2matrix(pch2matrix G, const char *name)
{
  int       nc_file;
  int       result;

  /* Create NetCDF file */
  result = nc_create(name, NC_64BIT_OFFSET, &nc_file);
  assert(result == NC_NOERR);

  /* Write row cluster tree */
  write_cdfpart_cluster(G->rb->t, nc_file, "rc");

  /* Write row cluster basis */
  write_cdfpart_clusterbasis(G->rb, nc_file, "rb");

  if (G->cb->t != G->rb->t) {
    /* Write column cluster tree */
    write_cdfpart_cluster(G->cb->t, nc_file, "cc");
  }

  if (G->cb != G->rb) {
    /* Write column cluster basis */
    write_cdfpart_clusterbasis(G->cb, nc_file, "cb");
  }

  /* Write h2matrix */
  write_cdfpart_h2matrix(G, nc_file, "ma");

  /* Close NetCDF file */
  nc_close(nc_file);
}

ph2matrix
read_cdfcomplete_h2matrix(const char *name)
{
  pcluster  rc, cc;
  pclusterbasis rb, cb;
  ph2matrix G;
  int       nc_file, nc_sons;
  int       result;

  /* Open NetCDF file */
  result = nc_open(name, NC_NOWRITE, &nc_file);
  assert(result == NC_NOERR);

  /* Read row cluster tree */
  rc = read_cdfpart_cluster(nc_file, "rc");

  /* Read row cluster basis */
  rb = read_cdfpart_clusterbasis(nc_file, "rb", rc);

  /* Check whether there is a separate column cluster tree */
  result = nc_inq_varid(nc_file, "cc_sons", &nc_sons);
  if (result == NC_NOERR) {
    /* If there is, read it */
    cc = read_cdfpart_cluster(nc_file, "cc");
  }
  else {
    /* Otherwise assume it's the row cluster tree */
    cc = rc;
  }

  /* Check whether there is a separate column cluster basis */
  result = nc_inq_varid(nc_file, "cb_sons", &nc_sons);
  if (result == NC_NOERR) {
    /* If these is, read it */
    cb = read_cdfpart_clusterbasis(nc_file, "cb", cc);
  }
  else {
    /* Otherwise assume it's the row cluster basis */
    cb = rb;
  }

  /* Read h2matrix */
  G = read_cdfpart_h2matrix(nc_file, "ma", rb, cb);

  /* Close NetCDF file */
  nc_close(nc_file);

  return G;
}
#endif

/* ------------------------------------------------------------
 * Drawing
 * ------------------------------------------------------------ */

#ifdef USE_CAIRO
static void
cairodraw(cairo_t * cr, pch2matrix G, bool storage, uint levels)
{
  cairo_text_extents_t extents;
  char      buf[11];
  uint      rsons, csons;
  uint      rsize, csize;
  uint      roff, coff;
  uint      i, j;

  if (G->son && levels != 1) {
    rsons = G->rsons;
    csons = G->csons;

    coff = 0;
    for (j = 0; j < csons; j++) {
      roff = 0;
      for (i = 0; i < rsons; i++) {
	cairo_save(cr);
	cairo_translate(cr, coff, roff);
	cairodraw(cr, G->son[i + j * rsons], storage, levels - 1);
	cairo_restore(cr);

	roff += G->son[i + j * rsons]->rb->t->size;
      }
      assert(roff == G->rb->t->size);

      coff += G->son[j * rsons]->cb->t->size;
    }
    assert(coff == G->cb->t->size);
  }
  else {
    rsize = G->rb->t->size;
    csize = G->cb->t->size;

    if (G->son) {
      cairo_rectangle(cr, 0.0, 0.0, csize, rsize);
      cairo_save(cr);
      cairo_set_source_rgb(cr, 0.9, 0.9, 1.0);
      cairo_fill_preserve(cr);
      cairo_restore(cr);
      cairo_stroke(cr);
    }
    else if (G->u) {
      if (storage) {
	cairo_rectangle(cr, 0.0, 0.0, (G->cb->k > csize ? csize : G->cb->k),
			(G->rb->k > rsize ? rsize : G->rb->k));
	cairo_save(cr);
	cairo_set_source_rgb(cr, 1.0, 0.0, 1.0);
	cairo_fill(cr);
	cairo_restore(cr);

	if (G->rb->k > 0 && G->cb->k > 0) {
	  sprintf(buf, "%dx%d", G->rb->k, G->cb->k);
	  cairo_set_font_size(cr, UINT_MIN(rsize, csize) / 4.75);
	  cairo_text_extents(cr, buf, &extents);

	  cairo_save(cr);
	  cairo_translate(cr, (csize - extents.width) / 2,
			  (rsize + extents.height) / 2);
	  cairo_show_text(cr, buf);
	  cairo_restore(cr);
	}

	cairo_rectangle(cr, 0.0, 0.0, csize, rsize);
	cairo_stroke(cr);
      }
      else {
	cairo_rectangle(cr, 0.0, 0.0, csize, rsize);
	cairo_save(cr);
	cairo_set_source_rgb(cr, 0.2, 0.2, 1.0);
	cairo_fill_preserve(cr);
	cairo_restore(cr);
	cairo_stroke(cr);
      }
    }
    else if (G->f) {
      cairo_rectangle(cr, 0.0, 0.0, csize, rsize);
      cairo_save(cr);
      cairo_set_source_rgb(cr, 1.0, 0.0, 0.0);
      cairo_fill_preserve(cr);
      cairo_restore(cr);
      cairo_stroke(cr);
    }
    else {
      cairo_rectangle(cr, 0.0, 0.0, csize, rsize);
      cairo_stroke(cr);
    }
  }
}

void
draw_cairo_h2matrix(cairo_t * cr, pch2matrix G, bool storage, uint levels)
{
  double    sx, sy, ex, ey;
  uint      rsize, csize;
  double    scalex, scaley, scale;

  /* Save Cairo state */
  cairo_save(cr);

  /* Obtain size of block */
  rsize = G->rb->t->size;
  csize = G->cb->t->size;

  /* Obtain size of current Cairo bounding box */
  cairo_clip_extents(cr, &sx, &sy, &ex, &ey);

  /* Compute scaling factor */
  scalex = (ex - sx) / rsize;
  scaley = (ey - sy) / csize;
  scale = (scalex < scaley ? scalex : scaley);

  /* Center block in bounding box */
  cairo_translate(cr, 0.5 * (ex - sx - scale * rsize),
		  0.5 * (ey - sy - scale * csize));

  /* Scale coordinates */
  cairo_scale(cr, scale, scale);
  cairo_set_line_width(cr, cairo_get_line_width(cr) / scale);

  /* Start drawing */
  cairodraw(cr, G, storage, levels);

  /* Restore Cairo state */
  cairo_restore(cr);
}
#endif
