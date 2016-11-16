/* ------------------------------------------------------------
 * This is the file "dh2matrix.c" of the H2Lib package.
 * All rights reserved, Steffen Boerm 2015
 * ------------------------------------------------------------ */

#include "dh2matrix.h"

#include "basic.h"

/* DEBUGGING */
#include <stdio.h>

/* ------------------------------------------------------------
 * Constructors and destructors
 * ------------------------------------------------------------ */

pdh2matrix
new_dh2matrix(pdclusterbasis rb, pdclusterbasis cb)
{
  pdh2matrix h2;

  h2 = allocmem(sizeof(dh2matrix));

  h2->rb = rb;
  h2->cb = cb;

  h2->u = NULL;
  h2->f = NULL;

  h2->son = NULL;
  h2->rsons = 0;
  h2->csons = 0;

  h2->desc = 0;

  return h2;
}

pdh2matrix
new_uniform_dh2matrix(pdclusterbasis rb, uint rd, pdclusterbasis cb, uint cd)
{
  pdh2matrix h2;

  h2 = new_dh2matrix(rb, cb);

  h2->u = new_duniform(rb, rd, cb, cd);

  h2->desc = 1;

  return h2;
}

pdh2matrix
new_full_dh2matrix(pdclusterbasis rb, pdclusterbasis cb)
{
  pdh2matrix h2;

  h2 = new_dh2matrix(rb, cb);

  h2->f = new_amatrix(rb->t->size, cb->t->size);

  h2->desc = 1;

  return h2;
}

pdh2matrix
new_super_dh2matrix(pdclusterbasis rb, pdclusterbasis cb, uint rsons,
		    uint csons)
{
  pdh2matrix h2;
  uint      i, j;

  h2 = new_dh2matrix(rb, cb);

  h2->rsons = rsons;
  h2->csons = csons;

  h2->son = (pdh2matrix *) allocmem((size_t) sizeof(pdh2matrix) * rsons *
				    csons);
  for (j = 0; j < csons; j++)
    for (i = 0; i < rsons; i++)
      h2->son[i + j * rsons] = 0;

  return h2;
}

void
update_dh2matrix(pdh2matrix h2)
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
del_dh2matrix(pdh2matrix h2)
{
  uint      rsons = h2->rsons;
  uint      csons = h2->csons;
  uint      i, j;

  if (h2->son) {
    for (j = 0; j < csons; j++)
      for (i = 0; i < rsons; i++)
	del_dh2matrix(h2->son[i + j * rsons]);
    freemem(h2->son);
  }

  if (h2->f)
    del_amatrix(h2->f);

  if (h2->u)
    del_duniform(h2->u);

  freemem(h2);
}

pdh2matrix
clone_dh2matrix(pcdh2matrix dh2, pdclusterbasis rb, pdclusterbasis cb)
{

  uint      i, j;
  uint      rsons, csons;
  pdh2matrix Newdh2;

  assert(dh2->rb->t == rb->t);
  assert(dh2->cb->t == cb->t);

  rsons = dh2->rsons;
  csons = dh2->csons;

  if (rsons + csons > 0) {

    Newdh2 = new_super_dh2matrix(rb, cb, rsons, csons);

    for (j = 0; j < csons; j++) {
      for (i = 0; i < rsons; i++) {
	Newdh2->son[i + j * rsons] =
	  clone_dh2matrix(dh2->son[i + j * rsons], rb->son[i], cb->son[j]);
      }
    }
  }
  else {
    if (dh2->u) {
      Newdh2 = new_uniform_dh2matrix(rb, dh2->u->rd, cb, dh2->u->cd);
      copy_amatrix(false, &dh2->u->S, &Newdh2->u->S);
    }
    else {
      assert(dh2->f);
      Newdh2 = new_full_dh2matrix(rb, cb);
      copy_amatrix(false, dh2->f, Newdh2->f);
    }
  }

  update_dh2matrix(Newdh2);

  return Newdh2;
}

/* ------------------------------------------------------------
 Statistics
 ------------------------------------------------------------ */

size_t
getsize_dh2matrix(pcdh2matrix h2)
{
  size_t    sz;
  uint      rsons = h2->rsons;
  uint      csons = h2->csons;
  uint      i, j;

  sz = (size_t) sizeof(dh2matrix);

  if (h2->u)
    sz += getsize_duniform(h2->u);

  if (h2->f)
    sz += getsize_amatrix(h2->f);

  for (j = 0; j < csons; j++)
    for (i = 0; i < rsons; i++)
      sz += getsize_dh2matrix(h2->son[i + j * rsons]);

  return sz;
}

size_t
getnearsize_dh2matrix(pcdh2matrix h2)
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
      sz += getnearsize_dh2matrix(h2->son[i + j * rsons]);

  return sz;
}

size_t
getfarsize_dh2matrix(pcdh2matrix h2)
{
  size_t    sz;
  uint      rsons = h2->rsons;
  uint      csons = h2->csons;
  uint      i, j;

  sz = 0;

  if (h2->u)
    sz += getsize_duniform(h2->u);

  for (j = 0; j < csons; j++)
    for (i = 0; i < rsons; i++)
      sz += getfarsize_dh2matrix(h2->son[i + j * rsons]);

  return sz;
}

size_t
gettotalsize_dh2matrix(pcdh2matrix h2)
{
  size_t    sz;

  sz = getsize_dh2matrix(h2);
  sz += getsize_dclusterbasis(h2->rb);
  if (h2->rb != h2->cb)
    sz += getsize_dclusterbasis(h2->cb);

  return sz;
}

/* ------------------------------------------------------------
 * Build H^2-matrix based on block tree
 * ------------------------------------------------------------ */

pdh2matrix
buildfromblock_dh2matrix(pcdblock b, pdclusterbasis rb, pdclusterbasis cb)
{
  pdh2matrix G, G1;
  pcdblock  b1;
  pdclusterbasis rb1, cb1;
  uint      rsons, csons;
  uint      i, j;

  G = 0;

  if (b->son) {
    rsons = b->rsons;
    csons = b->csons;

    G = new_super_dh2matrix(rb, cb, rsons, csons);

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

	G1 = buildfromblock_dh2matrix(b1, rb1, cb1);

	G->son[i + j * rsons] = G1;
      }
  }
  else if (b->adm)
    G = new_uniform_dh2matrix(rb, b->rd, cb, b->cd);
  else
    G = new_full_dh2matrix(rb, cb);

  update_dh2matrix(G);

  return G;
}

/* ------------------------------------------------------------
 * Copy nearfield matrices from a given dense matrix
 * ------------------------------------------------------------ */

typedef struct {
  pcamatrix G;
} copynear_data;

static void
copynear_pre(pdh2matrix Gh, uint mname, uint rname, uint cname,
	     uint pardepth, void *data)
{
  const uint *ridx = Gh->rb->t->idx;
  const uint *cidx = Gh->cb->t->idx;
  copynear_data *cd = (copynear_data *) data;
  pcamatrix G = cd->G;
  longindex ldG = G->ld;
  pamatrix  f;
  longindex ldf;
  uint      i, j, ii, jj;

  (void) mname;
  (void) rname;
  (void) cname;
  (void) pardepth;

  if (Gh->f) {
    f = Gh->f;
    ldf = f->ld;

    for (j = 0; j < f->cols; j++) {
      jj = cidx[j];
      for (i = 0; i < f->rows; i++) {
	ii = ridx[i];

	f->a[i + j * ldf] = G->a[ii + jj * ldG];
      }
    }
  }
}

void
copynear_dh2matrix(pcamatrix G, pdh2matrix Gh)
{
  copynear_data cd;

  cd.G = G;

  iterate_dh2matrix(Gh, 0, 0, 0, max_pardepth, copynear_pre, 0, &cd);
}

/* ------------------------------------------------------------
 * Matrix-vector multiplication
 * ------------------------------------------------------------ */

/* Virtual subdivision of non-strict inadmissible leaves */
static void
addeval_subamatrix(field alpha, pcdclusterbasis rb,
		   pcdclusterbasis cb, pcamatrix f, pcavector xt, pavector yt)
{
  amatrix   tmp1;
  avector   tmp2, tmp3;
  pamatrix  f1;
  pavector  xt1, yt1;
  uint      xoff, yoff, xtoff, ytoff;
  uint      i, j;

  if (cb->sons > 0) {
    if (rb->sons > 0) {
      xoff = 0;
      xtoff = cb->koff[cb->directions];
      for (j = 0; j < cb->sons; j++) {
	xt1 =
	  init_sub_avector(&tmp2, (pavector) xt, cb->son[j]->ktree, xtoff);

	yoff = 0;
	ytoff = rb->koff[rb->directions];
	for (i = 0; i < rb->sons; i++) {
	  yt1 = init_sub_avector(&tmp3, yt, rb->son[i]->ktree, ytoff);

	  f1 =
	    init_sub_amatrix(&tmp1, (pamatrix) f, rb->son[i]->t->size, yoff,
			     cb->son[j]->t->size, xoff);

	  addeval_subamatrix(alpha, rb->son[i], cb->son[j], f1, xt1, yt1);

	  uninit_amatrix(f1);
	  uninit_avector(yt1);

	  yoff += rb->son[i]->t->size;
	  ytoff += rb->son[i]->ktree;
	}
	assert(yoff == rb->t->size);
	assert(ytoff == rb->ktree);

	uninit_avector(xt1);

	xoff += cb->son[j]->t->size;
	xtoff += cb->son[j]->ktree;
      }
      assert(xoff == cb->t->size);
      assert(xtoff == cb->ktree);
    }
    else {
      xoff = 0;
      xtoff = cb->koff[cb->directions];
      for (j = 0; j < cb->sons; j++) {
	xt1 =
	  init_sub_avector(&tmp2, (pavector) xt, cb->son[j]->ktree, xtoff);

	f1 = init_sub_amatrix(&tmp1, (pamatrix) f, rb->t->size, 0,
			      cb->son[j]->t->size, xoff);

	addeval_subamatrix(alpha, rb, cb->son[j], f1, xt1, yt);

	uninit_amatrix(f1);
	uninit_avector(xt1);

	xoff += cb->son[j]->t->size;
	xtoff += cb->son[j]->ktree;
      }
      assert(xoff == cb->t->size);
      assert(xtoff == cb->ktree);
    }
  }
  else {
    if (rb->sons > 0) {
      yoff = 0;
      ytoff = rb->koff[rb->directions];
      for (i = 0; i < rb->sons; i++) {
	yt1 = init_sub_avector(&tmp3, yt, rb->son[i]->ktree, ytoff);

	f1 = init_sub_amatrix(&tmp1, (pamatrix) f, rb->son[i]->t->size, yoff,
			      cb->t->size, 0);

	addeval_subamatrix(alpha, rb->son[i], cb, f1, xt, yt1);

	uninit_amatrix(f1);
	uninit_avector(yt1);

	yoff += rb->son[i]->t->size;
	ytoff += rb->son[i]->ktree;
      }
      assert(yoff == rb->t->size);
      assert(ytoff == rb->ktree);
    }
    else {
      xt1 = init_sub_avector(&tmp2, (pavector) xt, cb->t->size,
			     cb->koff[cb->directions]);
      yt1 =
	init_sub_avector(&tmp3, yt, rb->t->size, rb->koff[rb->directions]);

      mvm_amatrix_avector(alpha, false, f, xt1, yt1);

      uninit_avector(yt1);
      uninit_avector(xt1);
    }
  }
}

/* Virtual subdivision of non-strict inadmissible leaves */
static void
addevaltrans_subamatrix(field alpha, pcdclusterbasis rb,
			pcdclusterbasis cb, pcamatrix f, pcavector xt,
			pavector yt)
{
  amatrix   tmp1;
  avector   tmp2, tmp3;
  pamatrix  f1;
  pavector  xt1, yt1;
  uint      xoff, yoff, xtoff, ytoff;
  uint      i, j;

  if (cb->sons > 0) {
    if (rb->sons > 0) {
      yoff = 0;
      ytoff = cb->koff[cb->directions];
      for (j = 0; j < cb->sons; j++) {
	yt1 = init_sub_avector(&tmp2, yt, cb->son[j]->ktree, ytoff);

	xoff = 0;
	xtoff = rb->koff[rb->directions];
	for (i = 0; i < rb->sons; i++) {
	  xt1 = init_sub_avector(&tmp3, (pavector) xt, rb->son[i]->ktree,
				 xtoff);

	  f1 =
	    init_sub_amatrix(&tmp1, (pamatrix) f, rb->son[i]->t->size, xoff,
			     cb->son[j]->t->size, yoff);

	  addevaltrans_subamatrix(alpha, rb->son[i], cb->son[j], f1, xt1,
				  yt1);

	  uninit_amatrix(f1);
	  uninit_avector(xt1);

	  xoff += rb->son[i]->t->size;
	  xtoff += rb->son[i]->ktree;
	}
	assert(xoff == rb->t->size);
	assert(xtoff == rb->ktree);

	uninit_avector(yt1);

	yoff += cb->son[j]->t->size;
	ytoff += cb->son[j]->ktree;
      }
      assert(yoff == cb->t->size);
      assert(ytoff == cb->ktree);
    }
    else {
      yoff = 0;
      ytoff = cb->koff[cb->directions];
      for (j = 0; j < cb->sons; j++) {
	yt1 = init_sub_avector(&tmp2, yt, cb->son[j]->ktree, ytoff);

	f1 = init_sub_amatrix(&tmp1, (pamatrix) f, rb->t->size, 0,
			      cb->son[j]->t->size, yoff);

	addevaltrans_subamatrix(alpha, rb, cb->son[j], f1, xt, yt1);

	uninit_amatrix(f1);
	uninit_avector(yt1);

	yoff += cb->son[j]->t->size;
	ytoff += cb->son[j]->ktree;
      }
      assert(yoff == cb->t->size);
      assert(ytoff == cb->ktree);
    }
  }
  else {
    if (rb->sons > 0) {
      xoff = 0;
      xtoff = rb->koff[rb->directions];
      for (i = 0; i < rb->sons; i++) {
	xt1 =
	  init_sub_avector(&tmp3, (pavector) xt, rb->son[i]->ktree, xtoff);

	f1 = init_sub_amatrix(&tmp1, (pamatrix) f, rb->son[i]->t->size, xoff,
			      cb->t->size, 0);

	addevaltrans_subamatrix(alpha, rb->son[i], cb, f1, xt1, yt);

	uninit_amatrix(f1);
	uninit_avector(xt1);

	xoff += rb->son[i]->t->size;
	xtoff += rb->son[i]->ktree;
      }
      assert(xoff == rb->t->size);
      assert(xtoff == rb->ktree);
    }
    else {
      xt1 = init_sub_avector(&tmp2, (pavector) xt, rb->t->size,
			     rb->koff[rb->directions]);
      yt1 =
	init_sub_avector(&tmp3, yt, cb->t->size, cb->koff[cb->directions]);

      mvm_amatrix_avector(alpha, true, f, xt1, yt1);

      uninit_avector(yt1);
      uninit_avector(xt1);
    }
  }
}

void
fastaddeval_dh2matrix_avector(field alpha, pcdh2matrix h2, pcavector xt,
			      pavector yt)
{
  avector   tmp1, tmp2;
  pavector  xt1, yt1;
  pcdclusterbasis rb = h2->rb;
  pcdclusterbasis cb = h2->cb;
  uint      rsons = h2->rsons;
  uint      csons = h2->csons;
  uint      xtoff, ytoff;
  uint      i, j;

  assert(xt->dim == cb->ktree);
  assert(yt->dim == rb->ktree);

  if (h2->u) {
    fastaddeval_duniform_avector(alpha, h2->u, xt, yt);
  }
  else if (h2->f) {
    addeval_subamatrix(alpha, rb, cb, h2->f, xt, yt);
  }
  else if (h2->son) {
    if (h2->son[0]->cb == h2->cb) {
      assert(h2->csons == 1);

      if (h2->son[0]->rb == h2->rb) {
	assert(h2->rsons == 1);

	fastaddeval_dh2matrix_avector(alpha, h2->son[0], xt, yt);
      }
      else {
	assert(h2->rsons == rb->sons);

	ytoff = rb->koff[rb->directions];
	for (i = 0; i < rsons; i++) {
	  yt1 = init_sub_avector(&tmp1, yt, rb->son[i]->ktree, ytoff);

	  fastaddeval_dh2matrix_avector(alpha, h2->son[i], xt, yt1);

	  uninit_avector(yt1);

	  ytoff += rb->son[i]->ktree;
	}
	assert(ytoff == rb->ktree);
      }
    }
    else {
      assert(h2->csons == cb->sons);

      if (h2->son[0]->rb == h2->rb) {
	assert(h2->rsons == 1);

	xtoff = cb->koff[cb->directions];
	for (j = 0; j < csons; j++) {
	  xt1 = init_sub_avector(&tmp2, (pavector) xt, cb->son[j]->ktree,
				 xtoff);

	  fastaddeval_dh2matrix_avector(alpha, h2->son[j], xt1, yt);

	  uninit_avector(xt1);

	  xtoff += cb->son[j]->ktree;
	}
	assert(xtoff == cb->ktree);
      }
      else {
	assert(h2->rsons == rb->sons);

	xtoff = cb->koff[cb->directions];
	for (j = 0; j < csons; j++) {
	  xt1 = init_sub_avector(&tmp2, (pavector) xt, cb->son[j]->ktree,
				 xtoff);

	  ytoff = rb->koff[rb->directions];
	  for (i = 0; i < rsons; i++) {
	    yt1 = init_sub_avector(&tmp1, yt, rb->son[i]->ktree, ytoff);

	    fastaddeval_dh2matrix_avector(alpha, h2->son[i + j * rsons], xt1,
					  yt1);

	    uninit_avector(yt1);

	    ytoff += rb->son[i]->ktree;
	  }
	  assert(ytoff == rb->ktree);

	  uninit_avector(xt1);

	  xtoff += cb->son[j]->ktree;
	}
	assert(xtoff == cb->ktree);
      }
    }
  }
}

void
fastaddevaltrans_dh2matrix_avector(field alpha, pcdh2matrix h2,
				   pcavector xt, pavector yt)
{
  avector   tmp1, tmp2;
  pavector  xt1, yt1;
  pcdclusterbasis rb = h2->rb;
  pcdclusterbasis cb = h2->cb;
  uint      rsons = h2->rsons;
  uint      csons = h2->csons;
  uint      xtoff, ytoff;
  uint      i, j;

  assert(xt->dim == rb->ktree);
  assert(yt->dim == cb->ktree);

  if (h2->u) {
    fastaddevaltrans_duniform_avector(alpha, h2->u, xt, yt);
  }
  else if (h2->f) {
    addevaltrans_subamatrix(alpha, rb, cb, h2->f, xt, yt);
  }
  else if (h2->son) {
    if (h2->son[0]->cb == h2->cb) {
      assert(h2->csons == 1);

      if (h2->son[0]->rb == h2->rb) {
	assert(h2->rsons == 1);

	fastaddevaltrans_dh2matrix_avector(alpha, h2->son[0], xt, yt);
      }
      else {
	assert(h2->rsons == rb->sons);

	xtoff = rb->koff[rb->directions];
	for (i = 0; i < rsons; i++) {
	  xt1 = init_sub_avector(&tmp1, (pavector) xt, rb->son[i]->ktree,
				 xtoff);

	  fastaddevaltrans_dh2matrix_avector(alpha, h2->son[i], xt1, yt);

	  uninit_avector(xt1);

	  xtoff += rb->son[i]->ktree;
	}
	assert(xtoff == rb->ktree);
      }
    }
    else {
      assert(h2->csons == cb->sons);

      if (h2->son[0]->rb == h2->rb) {
	assert(h2->rsons == 1);

	ytoff = cb->koff[cb->directions];
	for (j = 0; j < csons; j++) {
	  yt1 = init_sub_avector(&tmp2, yt, cb->son[j]->ktree, ytoff);

	  fastaddevaltrans_dh2matrix_avector(alpha, h2->son[j], xt, yt1);

	  uninit_avector(yt1);

	  ytoff += cb->son[j]->ktree;
	}
	assert(ytoff == cb->ktree);
      }
      else {
	assert(h2->rsons == rb->sons);

	ytoff = cb->koff[cb->directions];
	for (j = 0; j < csons; j++) {
	  yt1 = init_sub_avector(&tmp2, yt, cb->son[j]->ktree, ytoff);

	  xtoff = rb->koff[rb->directions];
	  for (i = 0; i < rsons; i++) {
	    xt1 = init_sub_avector(&tmp1, (pavector) xt, rb->son[i]->ktree,
				   xtoff);

	    fastaddevaltrans_dh2matrix_avector(alpha, h2->son[i + j * rsons],
					       xt1, yt1);

	    uninit_avector(xt1);

	    xtoff += rb->son[i]->ktree;
	  }
	  assert(xtoff == rb->ktree);

	  uninit_avector(yt1);

	  ytoff += cb->son[j]->ktree;
	}
	assert(ytoff == cb->ktree);
      }
    }
  }
}

void
addeval_dh2matrix_avector(field alpha, pcdh2matrix h2, pcavector x,
			  pavector y)
{
  pavector  xt, yt;

  xt = newcoeffs_dclusterbasis(h2->cb);
  yt = newcoeffs_dclusterbasis(h2->rb);

  forward_dclusterbasis(h2->cb, x, xt);

  clear_avector(yt);
  fastaddeval_dh2matrix_avector(alpha, h2, xt, yt);

  backward_dclusterbasis(h2->rb, yt, y);

  del_avector(yt);
  del_avector(xt);
}

void
addevaltrans_dh2matrix_avector(field alpha, pcdh2matrix h2, pcavector x,
			       pavector y)
{
  pavector  xt, yt;

  xt = newcoeffs_dclusterbasis(h2->rb);
  yt = newcoeffs_dclusterbasis(h2->cb);

  forward_dclusterbasis(h2->rb, x, xt);

  clear_avector(yt);
  fastaddevaltrans_dh2matrix_avector(alpha, h2, xt, yt);

  backward_dclusterbasis(h2->cb, yt, y);

  del_avector(yt);
  del_avector(xt);
}

void
mvm_dh2matrix_avector(field alpha, bool h2trans, pcdh2matrix h2,
		      pcavector x, pavector y)
{
  if (h2trans)
    addevaltrans_dh2matrix_avector(alpha, h2, x, y);
  else
    addeval_dh2matrix_avector(alpha, h2, x, y);
}

/* ------------------------------------------------------------
 * Parallel matrix-vector multiplication
 * ------------------------------------------------------------ */

typedef struct _addevalblock addevalblock;
struct _addevalblock {
  pcdh2matrix h2;

  uint      xoff;

  addevalblock *next;
};

typedef struct {
  field     alpha;

  addevalblock **bn;

  pcavector xt;
  pavector  yt;

  uint     *yoff;

} addeval_data;

static addevalblock *
splitrow_addeval(pcdh2matrix h2, uint xoff, addevalblock * next)
{
  addevalblock *b;
  uint      j, off;

  b = 0;

  if (h2->son && h2->son[0]->rb == h2->rb) {
    assert(h2->rsons == 1);
    assert(h2->son[0]->cb != h2->cb);

    b = next;
    off = xoff + h2->cb->koff[h2->cb->directions];
    for (j = 0; j < h2->csons; j++) {
      b = splitrow_addeval(h2->son[j], off, b);

      off += h2->son[j]->cb->ktree;
    }
    assert(off == xoff + h2->cb->ktree);
  }
  else {
    b = (addevalblock *) allocmem(sizeof(addevalblock));
    b->h2 = h2;
    b->xoff = xoff;
    b->next = next;
  }

  return b;
}

static void
addeval_pre(pdclusterbasis rb, uint tname, uint pardepth, void *data)
{
  avector   tmp1, tmp2;
  pavector  xt1, yt1;
  addeval_data *ad = (addeval_data *) data;
  pcavector xt = ad->xt;
  pavector  yt = ad->yt;
  addevalblock **bn = ad->bn;
  addevalblock *b = ad->bn[tname];
  uint      yoff = ad->yoff[tname];
  field     alpha = ad->alpha;
  pcdh2matrix h2;
  pcdclusterbasis cb;
  uint      rsons, csons;
  uint      xoff;
  uint      i, j, off, tname1;

  (void) pardepth;

  if (rb->sons > 0) {
    tname1 = tname + 1;
    off = yoff + rb->koff[rb->directions];
    for (i = 0; i < rb->sons; i++) {
      ad->bn[tname1] = 0;
      ad->yoff[tname1] = off;

      tname1 += rb->son[i]->t->desc;
      off += rb->son[i]->ktree;
    }
    assert(tname1 == tname + rb->t->desc);
    assert(off == yoff + rb->ktree);
  }

  while (b) {
    h2 = b->h2;
    xoff = b->xoff;
    cb = h2->cb;

    if (h2->son) {
      rsons = h2->rsons;
      csons = h2->csons;

      tname1 = tname + 1;
      for (i = 0; i < rsons; i++) {
	assert(h2->son[i]->rb == rb->son[i]);

	if (h2->son[0]->cb == cb) {
	  assert(csons == 1);

	  bn[tname1] = splitrow_addeval(h2->son[i], xoff, bn[tname1]);
	}
	else {
	  off = xoff + cb->koff[cb->directions];
	  for (j = 0; j < csons; j++) {
	    bn[tname1] = splitrow_addeval(h2->son[i + j * rsons], off,
					  bn[tname1]);

	    off += cb->son[j]->ktree;
	  }
	  assert(off == xoff + cb->ktree);
	}

	tname1 += rb->son[i]->t->desc;
      }
      assert(tname1 == tname + rb->t->desc);
    }
    else {
      xt1 = init_sub_avector(&tmp1, (pavector) xt, cb->ktree, xoff);
      yt1 = init_sub_avector(&tmp2, yt, rb->ktree, yoff);
      fastaddeval_dh2matrix_avector(alpha, h2, xt1, yt1);
      uninit_avector(yt1);
      uninit_avector(xt1);
    }

    b = b->next;
  }
}

void
addeval_parallel_dh2matrix_avector(field alpha, pcdh2matrix h2,
				   pcavector x, pavector y)
{
  addeval_data ad;
  addevalblock *b, *bnext;
  pavector  xt, yt;
  uint      i;

  xt = newcoeffs_dclusterbasis(h2->cb);
  yt = newcoeffs_dclusterbasis(h2->rb);

  forward_dclusterbasis(h2->cb, x, xt);

  clear_avector(yt);

  ad.bn =
    (addevalblock **) allocmem(sizeof(addevalblock *) * h2->rb->t->desc);
  ad.yoff = (uint *) allocmem(sizeof(uint) * h2->rb->t->desc);
  ad.xt = xt;
  ad.yt = yt;
  ad.alpha = alpha;
  ad.bn[0] = splitrow_addeval(h2, 0, 0);
  ad.yoff[0] = 0;
  iterate_dclusterbasis(h2->rb, 0, max_pardepth, addeval_pre, 0, &ad);

  for (i = 0; i < h2->rb->t->desc; i++) {
    b = ad.bn[i];
    while (b) {
      bnext = b->next;
      freemem(b);
      b = bnext;
    }
  }

  freemem(ad.yoff);
  freemem(ad.bn);

  backward_dclusterbasis(h2->rb, yt, y);

  del_avector(yt);
  del_avector(xt);
}

static addevalblock *
splitcol_addeval(pcdh2matrix h2, uint xoff, addevalblock * next)
{
  addevalblock *b;
  uint      i, off;

  b = 0;

  if (h2->son && h2->son[0]->cb == h2->cb) {
    assert(h2->csons == 1);
    assert(h2->son[0]->rb != h2->rb);

    b = next;
    off = xoff + h2->rb->koff[h2->rb->directions];
    for (i = 0; i < h2->rsons; i++) {
      b = splitcol_addeval(h2->son[i], off, b);

      off += h2->son[i]->rb->ktree;
    }
    assert(off == xoff + h2->rb->ktree);
  }
  else {
    b = (addevalblock *) allocmem(sizeof(addevalblock));
    b->h2 = h2;
    b->xoff = xoff;
    b->next = next;
  }

  return b;
}

static void
addevaltrans_pre(pdclusterbasis cb, uint tname, uint pardepth, void *data)
{
  avector   tmp1, tmp2;
  pavector  xt1, yt1;
  addeval_data *ad = (addeval_data *) data;
  pcavector xt = ad->xt;
  pavector  yt = ad->yt;
  addevalblock **bn = ad->bn;
  addevalblock *b = ad->bn[tname];
  uint      yoff = ad->yoff[tname];
  field     alpha = ad->alpha;
  pcdh2matrix h2;
  pcdclusterbasis rb;
  uint      rsons, csons;
  uint      xoff;
  uint      i, j, off, tname1;

  (void) pardepth;

  if (cb->sons > 0) {
    tname1 = tname + 1;
    off = yoff + cb->koff[cb->directions];
    for (j = 0; j < cb->sons; j++) {
      ad->bn[tname1] = 0;
      ad->yoff[tname1] = off;

      tname1 += cb->son[j]->t->desc;
      off += cb->son[j]->ktree;
    }
    assert(tname1 == tname + cb->t->desc);
    assert(off == yoff + cb->ktree);
  }

  while (b) {
    h2 = b->h2;
    xoff = b->xoff;
    rb = h2->rb;

    if (h2->son) {
      rsons = h2->rsons;
      csons = h2->csons;

      tname1 = tname + 1;
      for (j = 0; j < csons; j++) {
	assert(h2->son[j * rsons]->cb == cb->son[j]);

	if (h2->son[0]->rb == rb) {
	  assert(rsons == 1);

	  bn[tname1] = splitcol_addeval(h2->son[j], xoff, bn[tname1]);
	}
	else {
	  off = xoff + rb->koff[rb->directions];
	  for (i = 0; i < rsons; i++) {
	    bn[tname1] = splitcol_addeval(h2->son[i + j * rsons], off,
					  bn[tname1]);

	    off += rb->son[i]->ktree;
	  }
	  assert(off == xoff + rb->ktree);
	}

	tname1 += cb->son[j]->t->desc;
      }
      assert(tname1 == tname + cb->t->desc);
    }
    else {
      xt1 = init_sub_avector(&tmp1, (pavector) xt, rb->ktree, xoff);
      yt1 = init_sub_avector(&tmp2, yt, cb->ktree, yoff);
      fastaddevaltrans_dh2matrix_avector(alpha, h2, xt1, yt1);
      uninit_avector(yt1);
      uninit_avector(xt1);
    }

    b = b->next;
  }
}

void
addevaltrans_parallel_dh2matrix_avector(field alpha, pcdh2matrix h2,
					pcavector x, pavector y)
{
  addeval_data ad;
  addevalblock *b, *bnext;
  pavector  xt, yt;
  uint      i;

  xt = newcoeffs_dclusterbasis(h2->rb);
  yt = newcoeffs_dclusterbasis(h2->cb);

  forward_dclusterbasis(h2->rb, x, xt);

  clear_avector(yt);

  ad.bn =
    (addevalblock **) allocmem(sizeof(addevalblock *) * h2->cb->t->desc);
  ad.yoff = (uint *) allocmem(sizeof(uint) * h2->cb->t->desc);
  ad.xt = xt;
  ad.yt = yt;
  ad.alpha = alpha;
  ad.bn[0] = splitcol_addeval(h2, 0, 0);
  ad.yoff[0] = 0;
  iterate_dclusterbasis(h2->cb, 0, max_pardepth, addevaltrans_pre, 0, &ad);

  for (i = 0; i < h2->cb->t->desc; i++) {
    b = ad.bn[i];
    while (b) {
      bnext = b->next;
      freemem(b);
      b = bnext;
    }
  }

  freemem(ad.yoff);
  freemem(ad.bn);

  backward_dclusterbasis(h2->cb, yt, y);

  del_avector(yt);
  del_avector(xt);
}

/* ------------------------------------------------------------
 * Slow direct matrix-vector multiplication,
 * for debugging purposes
 * ------------------------------------------------------------ */

void
slowaddeval_dh2matrix_avector(field alpha, pcdh2matrix h2, pcavector x,
			      pavector y)
{
  avector   tmp1, tmp2;
  pavector  xp, yp;
  pcdclusterbasis rb = h2->rb;
  pcdclusterbasis cb = h2->cb;
  uint      rsons = h2->rsons;
  uint      csons = h2->csons;
  uint      i, j;

  if (h2->u) {
    slowaddeval_duniform_avector(alpha, h2->u, x, y);
  }
  else if (h2->f) {
    xp = init_avector(&tmp1, cb->t->size);
    yp = init_avector(&tmp2, rb->t->size);

    for (j = 0; j < cb->t->size; j++)
      xp->v[j] = x->v[cb->t->idx[j]];

    clear_avector(yp);
    mvm_amatrix_avector(alpha, false, h2->f, xp, yp);

    for (i = 0; i < rb->t->size; i++)
      y->v[rb->t->idx[i]] += yp->v[i];

    uninit_avector(yp);
    uninit_avector(xp);
  }
  else if (h2->son) {
    for (j = 0; j < csons; j++)
      for (i = 0; i < rsons; i++)
	slowaddeval_dh2matrix_avector(alpha, h2->son[i + j * rsons], x, y);
  }
}

void
slowaddevaltrans_dh2matrix_avector(field alpha, pcdh2matrix h2,
				   pcavector x, pavector y)
{
  avector   tmp1, tmp2;
  pavector  xp, yp;
  pcdclusterbasis rb = h2->rb;
  pcdclusterbasis cb = h2->cb;
  uint      rsons = h2->rsons;
  uint      csons = h2->csons;
  uint      i, j;

  if (h2->u) {
    slowaddevaltrans_duniform_avector(alpha, h2->u, x, y);
  }
  else if (h2->f) {
    xp = init_avector(&tmp1, rb->t->size);
    yp = init_avector(&tmp2, cb->t->size);

    for (i = 0; i < rb->t->size; i++)
      xp->v[i] = x->v[rb->t->idx[i]];

    clear_avector(yp);
    mvm_amatrix_avector(alpha, true, h2->f, xp, yp);

    for (j = 0; j < cb->t->size; j++)
      y->v[cb->t->idx[j]] += yp->v[j];

    uninit_avector(yp);
    uninit_avector(xp);
  }
  else if (h2->son) {
    for (j = 0; j < csons; j++)
      for (i = 0; i < rsons; i++)
	slowaddevaltrans_dh2matrix_avector(alpha, h2->son[i + j * rsons], x,
					   y);
  }
}

/* ------------------------------------------------------------
 * Conversion to a full matrix
 * ------------------------------------------------------------ */

typedef struct {
  field     alpha;
  pamatrix  G;
} expand_data;

static void
expand_pre(pdh2matrix h2, uint mname, uint rname, uint cname,
	   uint pardepth, void *data)
{
  expand_data *ed = (expand_data *) data;
  pamatrix  G = ed->G;
  field     alpha = ed->alpha;
  longindex ldG = G->ld;
  longindex ldf;
  pcamatrix f;
  const uint *ridx;
  const uint *cidx;
  uint      i, j, ii, jj;

  (void) mname;
  (void) rname;
  (void) cname;
  (void) pardepth;

  if (h2->u)
    expand_duniform(alpha, h2->u, G);
  else if (h2->f) {
    f = h2->f;
    ldf = f->ld;
    ridx = h2->rb->t->idx;
    cidx = h2->cb->t->idx;

    for (j = 0; j < f->cols; j++) {
      jj = cidx[j];
      for (i = 0; i < f->rows; i++) {
	ii = ridx[i];
	G->a[ii + jj * ldG] += alpha * f->a[i + j * ldf];
      }
    }
  }
}

void
expand_dh2matrix(field alpha, pcdh2matrix h2, pamatrix G)
{
  expand_data ed;

  ed.alpha = alpha;
  ed.G = G;

  iterate_dh2matrix((pdh2matrix) h2, 0, 0, 0, max_pardepth, expand_pre, 0,
		    &ed);
}

/* ------------------------------------------------------------
 Hierarchical iterators
 ------------------------------------------------------------ */

static void
enumerate(uint bname, pdh2matrix h2, pdh2matrix *h2n)
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

pdh2matrix *
enumerate_dh2matrix(pdh2matrix h2)
{
  pdh2matrix *h2n;

  h2n = (pdh2matrix *) allocmem((size_t) sizeof(pdh2matrix) * h2->desc);

  enumerate(0, h2, h2n);

  return h2n;
}

void
iterate_dh2matrix(pdh2matrix G, uint mname, uint rname, uint cname,
		  uint pardepth,
		  void (*pre) (pdh2matrix G, uint mname, uint rname,
			       uint cname, uint pardepth, void *data),
		  void (*post) (pdh2matrix G, uint mname, uint rname,
				uint cname, uint pardepth, void *data),
		  void *data)
{
  pcdcluster rc, cc;
  uint      rsons, csons;
  pdh2matrix *Gn;
  uint     *mnamen, *rnamen, *cnamen;
  uint      mname1, rname1, cname1;
#ifdef USE_OPENMP
  uint      nthreads;		/* HACK: Solaris cc requires an lvalue */
#endif
  uint      i, j, k, l;
  int       p;

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
	assert(rsons == 1);
	assert(csons == 1);

	iterate_dh2matrix(G->son[0], mname + 1, rname, cname,
			  (pardepth > 0 ? pardepth - 1 : 0), pre, post, data);
      }
      else {
	/* No row son, but column sons */
	assert(rsons == 1);

	mname1 = mname + 1;
	cname1 = cname + 1;
	for (j = 0; j < csons; j++) {
	  iterate_dh2matrix(G->son[j * rsons], mname1, rname, cname1,
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
	assert(csons == 1);

	mname1 = mname + 1;
	rname1 = rname + 1;
	for (i = 0; i < rsons; i++) {
	  iterate_dh2matrix(G->son[i], mname1, rname1, cname,
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
	  Gn = (pdh2matrix *) allocmem(sizeof(pdh2matrix) * csons);
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
	      assert(G->rb != G->son[0]->rb);
	      rname1 = rname + 1;
	      for (l = 0; l < i; l++)
		rname1 += rc->son[l]->desc;

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
	    for (p = 0; p < (int) csons; p++)
	      iterate_dh2matrix(Gn[p], mnamen[p], rnamen[p], cnamen[p],
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
	  Gn = (pdh2matrix *) allocmem(sizeof(pdh2matrix) * rsons);
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
	      assert(G->cb != G->son[0]->cb);
	      cname1 = cname + 1;
	      for (l = 0; l < j; l++)
		cname1 += cc->son[l]->desc;

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
	    for (p = 0; p < (int) rsons; p++)
	      iterate_dh2matrix(Gn[p], mnamen[p], rnamen[p], cnamen[p],
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

/* ------------------------------------------------------------
 Matrix norms
 ------------------------------------------------------------ */

real
norm2_dh2matrix(pcdh2matrix DH2)
{
  return norm2_matrix((mvm_t) mvm_dh2matrix_avector, (void *) DH2,
		      DH2->rb->t->size, DH2->cb->t->size);
}


real
norm2diff_dh2matrix(pcdh2matrix a, pcdh2matrix b)
{
  pavector  x, y;
  avector   tmp1, tmp2;
  uint      rows, cols;
  uint      i;
  real      norm;

  rows = getrows_dh2matrix(a);
  cols = getcols_dh2matrix(a);

  assert(getrows_dh2matrix(b) == rows);
  assert(getcols_dh2matrix(b) == cols);

  x = init_avector(&tmp1, cols);
  y = init_avector(&tmp2, rows);

  random_avector(x);
  norm = norm2_avector(x);

  for (i = 0; i < NORM_STEPS && norm > 0.0; i++) {
    scale_avector(1.0 / norm, x);

    clear_avector(y);
    addeval_dh2matrix_avector(1.0, a, x, y);
    addeval_dh2matrix_avector(-1.0, b, x, y);

    clear_avector(x);
    addevaltrans_dh2matrix_avector(1.0, a, y, x);
    addevaltrans_dh2matrix_avector(-1.0, b, y, x);

    norm = norm2_avector(x);
  }

  uninit_avector(y);
  uninit_avector(x);

  return REAL_SQRT(norm);
}

/* ------------------------------------------------------------
 * Drawing
 * ------------------------------------------------------------ */

#ifdef USE_CAIRO
static void
cairodraw(cairo_t * cr, pcdh2matrix G, bool storage, bool ranks, uint levels)
{
  cairo_text_extents_t extents;
  char      buf[11];
  uint      rsons, csons;
  uint      rsize, csize;
  uint      roff, coff;
  uint      kr, kc;
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
	cairodraw(cr, G->son[i + j * rsons], storage, ranks, levels - 1);
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
	kr = G->rb->k[G->u->rd];
	if (kr > rsize)
	  kr = rsize;

	kc = G->cb->k[G->u->cd];
	if (kc > csize)
	  kc = csize;

	cairo_rectangle(cr, 0.0, 0.0, kc, kr);
	cairo_save(cr);
	cairo_set_source_rgb(cr, 1.0, 0.0, 1.0);
	cairo_fill(cr);
	cairo_restore(cr);

	if (ranks && kr > 0 && kc > 0) {
	  sprintf(buf, "%dx%d", kr, kc);
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
draw_cairo_dh2matrix(cairo_t * cr, pcdh2matrix G, bool storage, bool ranks,
		     uint levels)
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
  cairodraw(cr, G, storage, ranks, levels);

  /* Restore Cairo state */
  cairo_restore(cr);
}

#endif

typedef struct _admisblock admisblock;

/* List for admissible blocks */
struct _admisblock {

  uint      name;		/* number of the corresponding block */
  uint      rname;		/* number of row cluster */
  uint      cname;		/* number of col cluster */
  uint      father;		/* number of father */
  uint      son;		/* son number */
  uint      length;

  struct _admisblock *next;	/* Next admissible block */
};


typedef struct _compdata compdata;
typedef compdata *pcompdata;
struct _compdata {

  pdclusteroperator *nco;
  pdclusteroperator *nro;
  pdclusteroperator *nbco;

  pdclusteroperator *noro;
  pdclusteroperator *noco;
  pdh2matrix *nG;

  pdclusterbasis *ncb;
  pdclusterbasis *nrb;

  pctruncmode tm;
  preal     eps;
  bool      rows;

  admisblock **cblock;
  admisblock **rblock;

};

static admisblock *
create_newadmisblock(uint bname, uint rname, uint cname, uint father)
{

  admisblock *ab = (admisblock *) allocmem(sizeof(admisblock));

  ab->name = bname;
  ab->rname = rname;
  ab->cname = cname;
  ab->father = father;
  ab->next = NULL;
  return ab;
}

static void
resize_coupling(pdh2matrix A, uint mname, uint rname, uint cname,
		uint pardepth, void *data)
{

  pcompdata cdata = (pcompdata) data;
  pdclusteroperator co = cdata->nco[cname];
  pdclusteroperator ro = cdata->nro[rname];
  amatrix   tmp;
  pamatrix  T;
  uint      cdir, rdir;


  if (A->u) {
    assert(co->dir == ro->dir);
    /* Take directions for S_b from duniform */
    cdir = A->u->cd;
    rdir = A->u->rd;

    assert(ro->C[rdir].cols == co->C[cdir].cols);

    T = init_amatrix(&tmp, A->u->S.rows, co->C[cdir].rows);
    clear_amatrix(T);
    /* First compute T = S_{b} * R^*_{s}, than resize S_{b} and save R_{t} * T in S_{b} */
    addmul_amatrix(1.0, false, &A->u->S, true, &co->C[cdir], T);
    resize_amatrix(&A->u->S, ro->C[rdir].rows, co->C[cdir].rows);
    clear_amatrix(&A->u->S);
    addmul_amatrix(1.0, false, &ro->C[rdir], false, T, &A->u->S);
    uninit_amatrix(T);

  }
}

void
resize_coupling_dh2matrix(pdh2matrix A, pdclusteroperator ro,
			  pdclusteroperator co)
{

  compdata  cdata;

  cdata.nco = enumerate_dclusteroperator(A->cb->t, co);
  cdata.nro = enumerate_dclusteroperator(A->rb->t, ro);

  iterate_dh2matrix(A, 0, 0, 0, max_pardepth, resize_coupling, 0, &cdata);

  freemem(cdata.nco);
  freemem(cdata.nro);
}

static void
build_dh2matrix_pre(pdh2matrix Gh, uint mname, uint rname, uint cname,
		    uint pardepth, void *data)
{

  pcompdata cdata = (pcompdata) data;
  pdclusteroperator co = cdata->nco[cname];
  pdclusteroperator ro = cdata->nro[rname];
  pdclusterbasis cb = cdata->ncb[cname];
  pdclusterbasis rb = cdata->nrb[rname];
  uint      rsons, csons;

  assert(Gh->cb->t == cb->t);
  assert(Gh->rb->t == rb->t);
  assert(cb->t == co->t);
  assert(rb->t == ro->t);

  /* Matrix is still subdivided */
  if (Gh->son) {
    rsons = Gh->rsons;
    csons = Gh->csons;

    cdata->nG[mname] = new_super_dh2matrix(rb, cb, rsons, csons);
  }
  else {			/* Two possible cases */
    if (Gh->u) {		/* Farfield */
      cdata->nG[mname] = new_uniform_dh2matrix(rb, Gh->u->rd, cb, Gh->u->cd);
      clear_amatrix(&cdata->nG[mname]->u->S);
      add_projected_duniform(Gh->u, ro, co, cdata->nG[mname]->u);
    }
    else {			/* Nearfield */
      cdata->nG[mname] = new_full_dh2matrix(rb, cb);
      copy_amatrix(false, Gh->f, cdata->nG[mname]->f);
    }
  }
}


static void
build_dh2matrix_post(pdh2matrix Gh, uint mname, uint rname, uint cname,
		     uint pardepth, void *data)
{

  pcompdata cdata = (pcompdata) data;
  uint      mname1;
  uint      rsons, csons;
  uint      i, j;

  /* Now we are searching for the same way as in the pre function to find the right sons ...  */
  if (Gh->son) {
    rsons = Gh->rsons;
    csons = Gh->csons;

    if (Gh->rb == Gh->son[0]->rb) {
      if (Gh->cb == Gh->son[0]->cb) {
	/* Just one son */
	assert(rsons == 1);
	assert(csons == 1);

	mname1 = mname + 1;
	cdata->nG[mname]->son[0] = cdata->nG[mname1];
      }
      else {
	/* No row son, but column sons */
	assert(rsons == 1);

	mname1 = mname + 1;
	for (j = 0; j < csons; j++) {
	  cdata->nG[mname]->son[j] = cdata->nG[mname1];
	  mname1 += Gh->son[j * rsons]->desc;
	}
	assert(mname1 == mname + Gh->desc);
      }
    }
    else {
      if (Gh->cb == Gh->son[0]->cb) {
	/* No column son, but row sons */
	assert(csons == 1);

	mname1 = mname + 1;
	for (i = 0; i < rsons; i++) {
	  cdata->nG[mname]->son[i] = cdata->nG[mname1];
	  mname1 += Gh->son[i]->desc;
	}
	assert(mname1 == mname + Gh->desc);
      }
      else {
	mname1 = mname + 1;
	for (j = 0; j < csons; j++) {
	  for (i = 0; i < rsons; i++) {
	    cdata->nG[mname]->son[i + j * rsons] = cdata->nG[mname1];
	    mname1 += Gh->son[i + j * rsons]->desc;
	  }
	}
      }
    }
  }
  update_dh2matrix(cdata->nG[mname]);
}

pdh2matrix
build_projected_dh2matrix(pcdh2matrix dh2, pdclusterbasis rb,
			  pdclusterbasis cb, pdclusteroperator ro,
			  pdclusteroperator co)
{

  pdh2matrix *newdh2;
  compdata  cdata;

  cdata.nrb = enumerate_dclusterbasis(rb->t, rb);
  cdata.ncb = enumerate_dclusterbasis(cb->t, cb);
  cdata.nro = enumerate_dclusteroperator(ro->t, ro);
  cdata.nco = enumerate_dclusteroperator(co->t, co);

  newdh2 = (pdh2matrix *) allocmem(sizeof(pdh2matrix) * dh2->desc);
  cdata.nG = newdh2;		/*This will be the new matrix */

  iterate_dh2matrix((pdh2matrix) dh2, 0, 0, 0, max_pardepth,
		    build_dh2matrix_pre, build_dh2matrix_post, &cdata);
  /* cleaning up */
  freemem(cdata.nrb);
  freemem(cdata.ncb);
  freemem(cdata.nro);
  freemem(cdata.nco);

  return newdh2[0];
}

static void
create_lists(uint bname, uint rname, uint cname, uint rfather, uint cfather,
	     uint rson, uint cson, void *data, admisblock ** rlist,
	     admisblock ** clist)
{

  pcompdata cdata = (pcompdata) data;
  uint      i, j;
  uint      bname1, rname1, cname1;
  pdh2matrix G = cdata->nG[bname];
  admisblock *abr = *rlist;
  admisblock *abc = *clist;

  if (G->son != NULL) {

    abr->father = rfather;
    abc->father = cfather;
    abr->son = rson;
    abc->son = cson;
    bname1 = bname + 1;
    cname1 = (G->csons > 0 ? cname + 1 : cname);

    for (j = 0; j < G->csons; j++) {
      rname1 = (G->rsons > 0 ? rname + 1 : rname);

      for (i = 0; i < G->rsons; i++) {
	create_lists(bname1, rname1, cname1, rname, cname, i, j, data,
		     &(cdata->rblock[rname1]), &(cdata->cblock[cname1]));
	rname1 += cdata->nG[bname1]->rb->t->desc;
	bname1 += cdata->nG[bname1]->desc;
      }
      assert(rname1 == rname + G->rb->t->desc);

      cname1 += G->son[j * G->rsons]->cb->t->desc;
    }
    assert(cname1 == cname + G->cb->t->desc);
    assert(bname1 == bname + G->desc);
  }
  else {
    abr->father = rfather;
    abc->father = cfather;
    abr->son = rson;
    abc->son = cson;

    if (G->u) {
      admisblock *newr = create_newadmisblock(bname, rname, cname, rfather);
      newr->son = rson;
      /* save number bname in list for rows */
      if (abr->length == 0) {
	abr->name = bname;
	abr->rname = rname;
	abr->cname = cname;
	abr->father = rfather;
	abr->son = rson;
	abr->length = 1;
      }
      else {
	while (abr->next != NULL) {
	  abr = abr->next;
	}
	abr->next = newr;
	cdata->rblock[rname]->length += 1;
      }
      /* save number bname in list for cols */
      admisblock *newc = create_newadmisblock(bname, rname, cname, cfather);
      newr->son = cson;
      if (abc->length == 0) {
	abc->name = bname;
	abc->rname = rname;
	abc->cname = cname;
	abc->father = cfather;
	abc->son = cson;
	abc->length = 1;
      }
      else {
	while (abc->next != NULL) {
	  abc = abc->next;
	}
	abc->next = newc;
	cdata->cblock[cname]->length += 1;
      }
    }
  }
}


static void
dweights(pdclusterbasis cb, uint tname, uint pardepth, void *data)
{

  pcompdata cdata = (pcompdata) data;
  pdclusteroperator yco = cdata->nco[tname];
  pdclusteroperator yro = cdata->nro[tname];
  pdclusteroperator oco;
  pdclusteroperator oro;
  pctruncmode tm = cdata->tm;
  bool      rows = cdata->rows;
  pdh2matrix G;
  uint    **information;

  uint      i, j;
  uint      n, off;
  uint      all, size, length;
  uint      rd, cd;
  amatrix   tmp1, tmp2;
  pamatrix  What, What1;
  uint      dim, m;
  avector   b1;
  pavector  tau;
  real      alpha, norm;
  admisblock *ab;

  pamatrix  Yhat, Yhat1;
  real      zeta_age;

  /* Find out if is the row or col cluster basis and take corresponing list */
  if (rows == true) {
    ab = cdata->rblock[tname];
  }
  else {
    ab = cdata->cblock[tname];
  }
  /* If the list isn't empty we have to do something */
  if (ab->length > 0) {

    /* Set up basis change matrices from orthogonalization */
    if (cdata->noco != NULL) {
      oco = cdata->noco[ab->cname];
    }
    else {
      oco = NULL;
    }
    if (cdata->noro != NULL) {
      oro = cdata->noro[ab->rname];
    }
    else {
      oro = NULL;
    }
    alpha = 1.0;

    /* Two different versions (row or col) */
    if (rows == true) {		/* Row version */
      length = ab->length;
      information = (uint **) allocmem(sizeof(uint *) * (length));
      for (j = 0; j < length; j++) {
	information[j] = (uint *) allocmem(sizeof(uint) * 2);
	information[j][0] = 0;	/* direction */
	information[j][1] = 0;	/* all clusteroperators or basis together ....for counting */
      }

      all = 1;			/* number of different cluster for computing weights */

      /* First entry */
      information[0][0] = cdata->nG[ab->name]->u->rd;
      information[0][1] =
	(oco ? oco->krow[cdata->nG[ab->name]->u->cd] : cdata->
	 nG[ab->name]->u->cb->k[cdata->nG[ab->name]->u->cd]);

      while (ab->next != NULL) {
	ab = ab->next;
	/* Find direction if it is already used or next empty array */
	j = 0;
	while ((information[j][1] > 0)
	       && (information[j][0] != cdata->nG[ab->name]->u->cd)) {
	  j += 1;
	  assert(j < length);
	}
	oco = (oco ? cdata->noco[ab->cname] : NULL);
	information[j][0] = cdata->nG[ab->name]->u->rd;
	information[j][1] +=
	  (oco ? oco->krow[cdata->nG[ab->name]->u->cd] : cdata->
	   nG[ab->name]->u->cb->k[cdata->nG[ab->name]->u->cd]);

	all += 1;
      }
      /* Set up matrices */
      j = 0;
      if (oco) {		/* Case 1 the original cluster basis was not(!) orthogonalized What <- R_s S*_{t,s} */
	while (all > 0) {
	  /* size of matrix */
	  rd = information[j][0];
	  size = information[j][1];
	  off = 0;
	  ab = cdata->rblock[tname];

	  /* Compute */
	  while ((off < size) && (ab != NULL)) {
	    if (cdata->nG[ab->name]->u->rd == rd) {	/* element of list belongs to this direction */
	      cd = cdata->nG[ab->name]->u->cd;
	      G = cdata->nG[ab->name];
	      oco = cdata->noco[ab->cname];
	      assert(G->u->cd == cd);
	      assert(G->u->rd == rd);
	      assert(tname == ab->rname);
	      n = yro->krow[rd] + oco->krow[cd];
	      What = init_amatrix(&tmp1, n, G->rb->k[rd]);
	      clear_amatrix(What);

	      /* Factor for Truncation, if necessary */
	      if (tm && tm->blocks) {
		if (tm->frobenius) {	// case Frobenius norm 
		  norm = normfrob_fast_duniform(G->u, oro, oco);
		}
		else {		// case Euclidean norm 
		  norm = norm2_fast_duniform(G->u, oro, oco);
		}
		alpha = (norm > 0.0 ? 1.0 / norm : 1.0);
	      }
	      /* If we already have something */
	      if (yro->krow[rd] > 0) {
		What1 =
		  init_sub_amatrix(&tmp2, What, yro->krow[rd], 0,
				   yro->C[rd].cols, 0);
		clear_amatrix(What1);
		copy_amatrix(false, &yro->C[rd], What1);
		uninit_amatrix(What1);
	      }
	      /* Compute local weight */
	      What1 =
		init_sub_amatrix(&tmp2, What, oco->krow[cd], yro->krow[rd],
				 G->rb->k[rd], 0);
	      clear_amatrix(What1);
	      addmul_amatrix(alpha, false, &oco->C[cd], true, &G->u->S,
			     What1);
	      uninit_amatrix(What1);

	      off += oco->krow[cd];

	      /* QR for smaller matrix */
	      dim = UINT_MIN(G->rb->k[rd], n);
	      tau = init_avector(&b1, dim);
	      qrdecomp_amatrix(What, tau);
	      uninit_avector(tau);

	      resize_dclusteroperator(yro, dim, G->rb->k[rd], rd);
	      copy_upper_amatrix(What, false, &yro->C[rd]);
	      uninit_amatrix(What);
	    }
	    ab = ab->next;
	  }
	  /* prepare for next round */
	  all -= 1;
	  j += 1;		/* next element of array */
	}
      }
      else {			/* Case 2 the original cluster basis was already orthogonalized What <- S*_{t,s} */
	while (all > 0) {
	  rd = information[j][0];
	  size = information[j][1];
	  off = 0;
	  ab = cdata->rblock[tname];
	  /* Compute */
	  while ((off < size) && (ab != NULL)) {
	    if (cdata->nG[ab->name]->u->rd == rd) {	/* element of list belongs to this direction */
	      cd = cdata->nG[ab->name]->u->cd;
	      G = cdata->nG[ab->name];
	      assert(G->u->cd == cd);
	      assert(G->u->rd == rd);
	      n = yro->krow[rd] + G->cb->k[cd];
	      What = init_amatrix(&tmp1, n, G->rb->k[rd]);
	      clear_amatrix(What);

	      /* Factor for Truncation, if necessary */
	      if (tm && tm->blocks) {
		if (tm->frobenius) {	// case Frobenius norm 
		  norm = normfrob_fast_duniform(G->u, oro, oco);
		}
		else {		// case Euclidean norm 
		  norm = norm2_fast_duniform(G->u, oro, oco);
		}
		alpha = (norm > 0.0 ? 1.0 / norm : 1.0);
	      }
	      if (yro->krow[rd] > 0) {
		What1 =
		  init_sub_amatrix(&tmp2, What, yro->krow[rd], 0,
				   yro->C[rd].cols, 0);
		clear_amatrix(What1);
		copy_amatrix(false, &yro->C[rd], What1);
		uninit_amatrix(What1);
	      }

	      What1 =
		init_sub_amatrix(&tmp2, What, G->cb->k[cd], yro->krow[rd],
				 G->rb->k[rd], 0);
	      clear_amatrix(What1);
	      add_amatrix(alpha, true, &G->u->S, What1);
	      uninit_amatrix(What1);

	      off += G->u->cb->k[cd];

	      /* QR */
	      dim = UINT_MIN(G->rb->k[rd], n);
	      tau = init_avector(&b1, dim);
	      qrdecomp_amatrix(What, tau);
	      uninit_avector(tau);

	      resize_dclusteroperator(yro, dim, G->rb->k[rd], rd);
	      copy_upper_amatrix(What, false, &yro->C[rd]);
	      uninit_amatrix(What);
	    }
	    ab = ab->next;
	  }
	  /* prepare for next round */
	  all -= 1;
	  j += 1;		/* next element of array */
	}
      }
    }
    else {
      length = ab->length;
      information = (uint **) allocmem(sizeof(uint *) * (length));
      for (j = 0; j < length; j++) {
	information[j] = (uint *) allocmem(sizeof(uint) * 2);
	information[j][0] = 0;	/* direction */
	information[j][1] = 0;	/* size clusteroperator or clusterbasis */
      }
      all = 1;			/* number of different cluster for computing weights */
      /* First entry */
      information[0][0] = cdata->nG[ab->name]->u->cd;
      information[0][1] =
	(oro ? oro->krow[cdata->nG[ab->name]->u->rd] : cdata->
	 nG[ab->name]->u->rb->k[cdata->nG[ab->name]->u->rd]);

      while (ab->next != NULL) {
	ab = ab->next;
	/* Find direction if already used or next empty array */
	j = 0;
	while ((information[j][1] > 0)
	       && (information[j][0] != cdata->nG[ab->name]->u->rd)) {
	  j += 1;
	  assert(j < length);
	}
	oro = (oro ? cdata->noro[ab->rname] : NULL);
	information[j][0] = cdata->nG[ab->name]->u->cd;
	information[j][1] +=
	  (oro ? oro->krow[cdata->nG[ab->name]->u->rd] : cdata->
	   nG[ab->name]->u->rb->k[cdata->nG[ab->name]->u->rd]);
	all += 1;
      }
      /* Set up matrices */
      j = 0;
      if (oro) {		/* Case 1 the original cluster basis was not(!) orthogonalized What <- R_s S*_{t,s} */
	while (all > 0) {
	  /* size of matrix */
	  cd = information[j][0];
	  size = information[j][1];
	  off = 0;
	  ab = cdata->cblock[tname];
	  /* Compute */
	  while ((off < size) && (ab != NULL)) {
	    if (cdata->nG[ab->name]->u->cd == cd) {	/* element of list belongs to this direction */
	      G = cdata->nG[ab->name];
	      rd = G->u->rd;
	      oro = cdata->noro[ab->rname];
	      assert(G->u->cd == cd);
	      assert(G->u->rd == rd);
	      assert(tname == ab->cname);
	      n = oro->krow[rd] + yco->krow[cd];
	      What = init_amatrix(&tmp1, n, G->cb->k[cd]);
	      clear_amatrix(What);
	      /* Factor for Truncation, if necessary */
	      if (tm && tm->blocks) {
		if (tm->frobenius) {	// case Frobenius norm 
		  norm = normfrob_fast_duniform(G->u, oro, oco);
		}
		else {		// case Euclidean norm 
		  norm = norm2_fast_duniform(G->u, oro, oco);
		}
		alpha = (norm > 0.0 ? 1.0 / norm : 1.0);
	      }

	      if (yco->krow[cd] > 0) {
		What1 =
		  init_sub_amatrix(&tmp2, What, yco->krow[cd], 0,
				   G->cb->k[cd], 0);
		copy_amatrix(false, &yco->C[cd], What1);
		uninit_amatrix(What1);
	      }

	      What1 =
		init_sub_amatrix(&tmp2, What, oro->krow[rd], yco->krow[cd],
				 G->cb->k[cd], 0);
	      clear_amatrix(What1);
	      addmul_amatrix(alpha, false, &oro->C[rd], false, &G->u->S,
			     What1);
	      uninit_amatrix(What1);

	      off += oro->krow[rd];

	      /* QR */
	      dim = UINT_MIN(G->cb->k[cd], n);
	      tau = init_avector(&b1, dim);
	      qrdecomp_amatrix(What, tau);
	      uninit_avector(tau);

	      resize_dclusteroperator(yco, dim, G->cb->k[cd], cd);
	      copy_upper_amatrix(What, false, &yco->C[cd]);
	      uninit_amatrix(What);
	    }
	    ab = ab->next;
	  }
	  /* prepare for next round */
	  all -= 1;
	  j += 1;		/* next element of array */
	}
      }
      else {			/* Case 2 the original cluster basis was already orthogonalized What <- S*_{t,s} */
	while (all > 0) {
	  /* size of matrix */
	  cd = information[j][0];
	  size = information[j][1];
	  off = 0;
	  ab = cdata->cblock[tname];
	  /* Compute */
	  while ((off < size) && (ab != NULL)) {
	    if (cdata->nG[ab->name]->u->cd == cd) {	/* element of list belongs to this direction */
	      G = cdata->nG[ab->name];
	      rd = G->u->rd;
	      assert(G->u->cd == cd);
	      assert(G->u->rd == rd);
	      n = G->rb->k[rd] + yco->krow[cd];
	      What = init_amatrix(&tmp1, n, G->cb->k[cd]);
	      clear_amatrix(What);

	      /* Factor for Truncation, if necessary */
	      if (tm && tm->blocks) {
		if (tm->frobenius) {	// case Frobenius norm 
		  norm = normfrob_fast_duniform(G->u, oro, oco);
		}
		else {		// case Euclidean norm 
		  norm = norm2_fast_duniform(G->u, oro, oco);
		}
		alpha = (norm > 0.0 ? 1.0 / norm : 1.0);
	      }

	      if (yco->krow[cd] > 0) {
		What1 =
		  init_sub_amatrix(&tmp2, What, yco->krow[cd], 0,
				   G->cb->k[cd], 0);
		copy_amatrix(false, &yco->C[cd], What1);
		uninit_amatrix(What1);
	      }

	      What1 =
		init_sub_amatrix(&tmp2, What, G->rb->k[rd], yco->krow[cd],
				 G->cb->k[cd], 0);
	      clear_amatrix(What1);
	      add_amatrix(alpha, false, &G->u->S, What1);
	      uninit_amatrix(What1);

	      off += G->u->rb->k[rd];

	      /* QR */
	      dim = UINT_MIN(G->cb->k[cd], n);
	      tau = init_avector(&b1, dim);
	      qrdecomp_amatrix(What, tau);
	      uninit_avector(tau);

	      resize_dclusteroperator(yco, dim, G->cb->k[cd], cd);
	      copy_upper_amatrix(What, false, &yco->C[cd]);
	      uninit_amatrix(What);

	    }
	    ab = ab->next;
	  }
	  /* prepare for next round */
	  all -= 1;
	  j += 1;		/* next element of array */
	}
      }
    }
    /* first cleaning up */
    for (j = 0; j < length; j++) {
      freemem(information[j]);
    }
    freemem(information);
  }

  /* Now combine to total weights */
  if ((rows == true) && (tname != 0)) {
    ab = cdata->rblock[tname];
    zeta_age = (tm ? tm->zeta_age : 1.0);
    for (i = 0; i < cb->directions; i++) {
      /* Finding all corresponding father directions and save size */
      m = 0;
      for (j = 0; j < cdata->nrb[ab->father]->directions; j++) {
	m =
	  (cdata->nrb[ab->father]->dirson[ab->son][j] ==
	   i ? m + cdata->nro[ab->father]->C[j].rows : m);
      }
      /* Set up matrix */
      Yhat = init_amatrix(&tmp1, m + yro->C[i].rows, yro->C[i].cols);
      /* Evaluate ZE* from father and collect in upper part of Yhat */
      off = 0;
      if (m > 0) {
	for (j = 0; j < cdata->nrb[ab->father]->directions; j++) {
	  if (cdata->nrb[ab->father]->dirson[ab->son][j] == i) {
	    Yhat1 =
	      init_sub_amatrix(&tmp2, Yhat, cdata->nro[ab->father]->C[j].rows,
			       off,
			       cdata->nrb[ab->father]->E[ab->son][j].rows, 0);
	    clear_amatrix(Yhat1);
	    assert(cdata->nro[ab->father]->C[j].cols ==
		   cdata->nrb[ab->father]->E[ab->son][j].cols);
	    addmul_amatrix(zeta_age, false, &cdata->nro[ab->father]->C[j],
			   true, &cdata->nrb[ab->father]->E[ab->son][j],
			   Yhat1);
	    off += cdata->nro[ab->father]->C[j].rows;
	    uninit_amatrix(Yhat1);
	  }
	}
	assert(off == m);
      }
      Yhat1 =
	init_sub_amatrix(&tmp2, Yhat, yro->C[i].rows, m, yro->C[i].cols, 0);
      copy_amatrix(false, &yro->C[i], Yhat1);
      uninit_amatrix(Yhat1);

      dim = UINT_MIN(m + yro->C[i].rows, yro->C[i].cols);
      tau = init_avector(&b1, dim);
      qrdecomp_amatrix(Yhat, tau);
      uninit_avector(tau);

      resize_dclusteroperator(yro, dim, yro->C[i].cols, i);
      copy_upper_amatrix(Yhat, false, &yro->C[i]);
      uninit_amatrix(Yhat);
    }
  }
  else {
    if (tname != 0) {
      zeta_age = (tm ? tm->zeta_age : 1.0);
      ab = cdata->cblock[tname];
      /* For all his direction take a look at the cluster operator */
      for (i = 0; i < cb->directions; i++) {
	/* Finding all father directions and save size */
	m = 0;
	for (j = 0; j < cdata->ncb[ab->father]->directions; j++) {
	  m =
	    (cdata->ncb[ab->father]->dirson[ab->son][j] ==
	     i ? m + cdata->nco[ab->father]->C[j].rows : m);
	}
	/* Set up matrix */
	Yhat = init_amatrix(&tmp1, m + yco->C[i].rows, yco->C[i].cols);
	/* Evaluate ZE* from father and collect in upper part of Yhat */
	off = 0;
	if (m > 0) {
	  for (j = 0; j < cdata->ncb[ab->father]->directions; j++) {
	    if (cdata->ncb[ab->father]->dirson[ab->son][j] == i) {
	      Yhat1 =
		init_sub_amatrix(&tmp2, Yhat,
				 cdata->nco[ab->father]->C[j].rows, off,
				 cdata->ncb[ab->father]->E[ab->son][j].rows,
				 0);
	      clear_amatrix(Yhat1);
	      assert(cdata->nco[ab->father]->C[j].cols ==
		     cdata->ncb[ab->father]->E[ab->son][j].cols);
	      addmul_amatrix(zeta_age, false, &cdata->nco[ab->father]->C[j],
			     true, &cdata->ncb[ab->father]->E[ab->son][j],
			     Yhat1);
	      off += cdata->nco[ab->father]->C[j].rows;
	      uninit_amatrix(Yhat1);
	    }
	  }
	  assert(off == m);
	}
	Yhat1 =
	  init_sub_amatrix(&tmp2, Yhat, yco->C[i].rows, m, yco->C[i].cols, 0);
	copy_amatrix(false, &yco->C[i], Yhat1);
	uninit_amatrix(Yhat1);

	dim = UINT_MIN(m + yco->C[i].rows, yco->C[i].cols);
	tau = init_avector(&b1, dim);
	qrdecomp_amatrix(Yhat, tau);
	uninit_avector(tau);

	resize_dclusteroperator(yco, dim, yco->C[i].cols, i);
	copy_upper_amatrix(Yhat, false, &yco->C[i]);
	uninit_amatrix(Yhat);
      }
    }
  }
}


static void
compress_dweights(pcdh2matrix G, pdclusteroperator oro, pdclusteroperator oco,
		  pdclusteroperator yro, pdclusteroperator yco,
		  pctruncmode tm)
{

  compdata  cdata;
  uint      i;
  uint      maxcname = G->cb->t->desc;
  uint      maxrname = G->rb->t->desc;
  admisblock **rab;
  admisblock **cab;
  cdata.nG = enumerate_dh2matrix((pdh2matrix) G);

  rab = (admisblock **) allocmem(sizeof(admisblock *) * maxrname);
  for (i = 0; i < maxrname; i++) {
    rab[i] = (admisblock *) allocmem(sizeof(admisblock));
    rab[i]->next = NULL;
    rab[i]->length = 0;
  }

  cab = (admisblock **) allocmem(sizeof(admisblock *) * maxcname);
  for (i = 0; i < maxcname; i++) {
    cab[i] = (admisblock *) allocmem(sizeof(admisblock));
    cab[i]->next = NULL;
    cab[i]->length = 0;
  }
  cdata.rblock = rab;
  cdata.cblock = cab;

  create_lists(0, 0, 0, 0, 0, 0, 0, (void *) &cdata, &rab[0], &cab[0]);

  cdata.tm = tm;
  if (oro) {
    cdata.noro = enumerate_dclusteroperator(G->rb->t, oro);
  }
  else {
    cdata.noro = NULL;
  }
  if (oco) {
    cdata.noco = enumerate_dclusteroperator(G->cb->t, oco);
  }
  else {
    cdata.noco = NULL;
  }
  cdata.nco = enumerate_dclusteroperator(G->cb->t, yco);
  cdata.nro = enumerate_dclusteroperator(G->rb->t, yro);
  cdata.nrb = enumerate_dclusterbasis(G->rb->t, G->rb);
  cdata.ncb = enumerate_dclusterbasis(G->cb->t, G->cb);

  /* Compute weights */

  cdata.rows = true;
  iterate_dclusterbasis(G->rb, 0, max_pardepth, dweights, NULL,
			(void *) &cdata);

  cdata.rows = false;
  iterate_dclusterbasis(G->cb, 0, max_pardepth, dweights, NULL,
			(void *) &cdata);


  freemem(cdata.nco);
  freemem(cdata.noco);
  freemem(cdata.nro);
  freemem(cdata.noro);
  freemem(cdata.nG);

  for (i = 0; i < maxcname; i++) {
    freemem(cab[i]);
  }
  freemem(cab);

  for (i = 0; i < maxrname; i++) {
    freemem(rab[i]);
  }
  freemem(rab);

}


/* Truncation */

static void
truncate_pre(pdclusterbasis cb, uint tname, uint pardepth, void *data)
{

  pcompdata cdata = (pcompdata) data;
  pctruncmode tm = cdata->tm;
  preal     eps = cdata->eps;
  real      zeta_level;
  uint      i;
  uint      tname_next;

  /* Weights for truncation */
  zeta_level = (tm ? tm->zeta_level : 1.0);

  /* Compute accuracies for sons */
  tname_next = tname + 1;
  for (i = 0; i < cb->sons; i++) {
    eps[tname_next] = eps[tname] * zeta_level;

    tname_next += cb->t->son[i]->desc;
  }
  assert(tname_next == tname + cb->t->desc);
}

static void
truncate_post(pdclusterbasis cnew, uint tname, uint pardepth, void *data)
{

  pcompdata cdata = (pcompdata) data;
  pdclusteroperator co = cdata->nco[tname];
  pdclusteroperator bco = cdata->nbco[tname];
  pcdclusterbasis cold = cdata->ncb[tname];
  pctruncmode tm = cdata->tm;
  real      eps = cdata->eps[tname];

  uint      i, j;
  amatrix   tmp1, tmp2, tmp3, tmp4;
  pamatrix  Q, Qhat, U, Y, V, V1, X;
  avector   b1;
  pavector  tau;
  realavector br;
  prealavector sigma;
  uint      kmax, kold, kop, k;
  uint      off;


  /* Hopefully they have the right size */
  assert(cold->t == cnew->t);
  assert(cold->t == co->t);
  assert(cold->t == bco->t);

  for (i = 0; i < cold->directions; i++) {

    /* new k max vor every son */
    kmax = 0;
    kold = 0;
    if (cold->sons > 0) {
      for (j = 0; j < cold->sons; j++) {
	kmax += bco->son[j]->C[cnew->dirson[j][i]].rows;	/* Rows for every son */
	kold = (cold->E[j][i].cols > kold ? cold->E[j][i].cols : kold);	/* Maximal values for columns */
      }
      V = init_amatrix(&tmp1, kmax, kold);	/* Assembled cluster basis */
      clear_amatrix(V);
      off = 0;
      /* Compute R_t'c'* E_t'c for every son t' in direction c */
      for (j = 0; j < cold->sons; j++) {
	V1 =
	  init_sub_amatrix(&tmp2, V, bco->son[j]->C[cnew->dirson[j][i]].rows,
			   off, cold->E[j][i].cols, 0);
	addmul_amatrix(1.0, false, &bco->son[j]->C[cnew->dirson[j][i]], false,
		       &cold->E[j][i], V1);
	off += bco->son[j]->C[cnew->dirson[j][i]].rows;
	uninit_amatrix(V1);
      }
      assert(off == kmax);
    }
    else {
      kmax = cold->V[i].rows;
      kold = cold->V[i].cols;

      V = init_amatrix(&tmp1, kmax, kold);
      clear_amatrix(V);
      copy_amatrix(false, &cold->V[i], V);
    }

    if (V->rows > 0 && V->cols > 0) {
      if (V->rows > co->C[i].rows) {
	/* Set up X^c_t <- V_t * Z^*_t */
	X = init_amatrix(&tmp2, kmax, co->C[i].rows);
	clear_amatrix(X);
	addmul_amatrix(1.0, false, V, true, &co->C[i], X);

	/* One more time QR, so we have some small matrices for the svd */
	tau = init_avector(&b1, X->cols);
	qrdecomp_amatrix(X, tau);

	/* Y <- R and save Qhat from th qr decomposition */
	Y = init_amatrix(&tmp3, X->cols, X->cols);
	copy_upper_amatrix(X, false, Y);
	Qhat = init_amatrix(&tmp4, X->rows, Y->rows);
	qrexpand_amatrix(X, tau, Qhat);
	uninit_avector(tau);
	uninit_amatrix(X);

	k = Y->rows;
	U = init_amatrix(&tmp2, k, k);
	clear_amatrix(U);
	sigma = init_realavector(&br, k);
	svd_amatrix(Y, sigma, U, 0);
	uninit_amatrix(Y);

	Q = init_amatrix(&tmp3, Qhat->rows, U->cols);
	clear_amatrix(Q);
	addmul_amatrix(1.0, false, Qhat, false, U, Q);
	uninit_amatrix(U);
	uninit_amatrix(Qhat);
      }
      else {			/* For smaller matrix use the transpose */
	/* (X_t^c )^* = Z_t * V_t^*  */
	X = init_amatrix(&tmp2, co->C[i].rows, kmax);
	clear_amatrix(X);
	addmul_amatrix(1.0, false, &co->C[i], true, V, X);

	/* One more time QR, so we have some small matrices for the svd */
	tau = init_avector(&b1, X->cols);
	qrdecomp_amatrix(X, tau);

	/* Y <- R* from th qr decomposition */
	Qhat = init_amatrix(&tmp3, kmax, kmax);
	copy_upper_amatrix(X, false, Qhat);
	uninit_amatrix(X);
	/* Need R* for the svd */
	Y = init_amatrix(&tmp2, kmax, kmax);
	copy_amatrix(true, Qhat, Y);
	uninit_amatrix(Qhat);
	uninit_avector(tau);

	k = Y->rows;
	Q = init_amatrix(&tmp3, k, k);
	clear_amatrix(Q);
	sigma = init_realavector(&br, k);
	svd_amatrix(Y, sigma, Q, 0);
	uninit_amatrix(Y);
      }
      /* Find right size */
      kop = findrank_truncmode(tm, eps, sigma);
      /* Finally the right Q  =) */
      uninit_realavector(sigma);
      U = init_sub_amatrix(&tmp2, Q, Q->rows, 0, kop, 0);

      setrank_dclusterbasis(cnew, i, kop);
      resize_dclusteroperator(bco, kop, kold, i);
      /* Basis change matrix, saved in co ! */
      if (V->rows > 0 && V->cols > 0) {
	clear_amatrix(&bco->C[i]);
	addmul_amatrix(1.0, true, U, false, V, &bco->C[i]);
      }
      /* The wanted cluster basis */
      if (cold->sons > 0) {	/* Non-leaf */
	off = 0;
	for (j = 0; j < cold->sons; j++) {
	  Qhat =
	    init_sub_amatrix(&tmp4, U, cnew->son[j]->k[cold->dirson[j][i]],
			     off, kop, 0);
	  copy_amatrix(false, Qhat, &cnew->E[j][i]);
	  off += cnew->son[j]->k[cold->dirson[j][i]];
	  uninit_amatrix(Qhat);
	}
      }
      else {			/* Leaf */
	copy_amatrix(false, U, &cnew->V[i]);
      }
      uninit_amatrix(V);
      uninit_amatrix(U);
      uninit_amatrix(Q);
    }
    else {
      uninit_amatrix(V);
      setrank_dclusterbasis(cnew, i, 0);
      resize_dclusteroperator(bco, 0, kold, i);

    }
  }
  update_dclusterbasis(cnew);
}


void
truncate_dclusterbasis(pdclusterbasis cold, pdclusterbasis cnew,
		       pdclusteroperator co, pdclusteroperator bco,
		       pctruncmode tm, real eps)
{


  compdata  cdata;

  cdata.nco = enumerate_dclusteroperator(cold->t, co);
  cdata.nbco = enumerate_dclusteroperator(cold->t, bco);
  cdata.ncb = enumerate_dclusterbasis(cold->t, cold);
  cdata.tm = tm;
  cdata.eps = allocreal(cold->t->desc);
  cdata.eps[0] = eps;

  iterate_dclusterbasis(cnew, 0, max_pardepth, truncate_pre, truncate_post,
			(void *) &cdata);

  freemem(cdata.nco);
  freemem(cdata.nbco);
  freemem(cdata.ncb);
  freemem(cdata.eps);
}


pdh2matrix
compress_dh2matrix_dh2matrix(pcdh2matrix G, bool rbortho, bool cbortho,
			     pctruncmode tm, real eps)
{


  pdh2matrix New_G;
  pdclusteroperator oco, oro;	/* For orthogonalization if it is needed */
  pdclusteroperator yco, yro;	/* For the second orthogonalization (Y_t) */
  pdclusteroperator bcc, bcr;	/* To save basis change */
  pdclusterbasis cb, rb;

  if (rbortho == true) {
    oro = NULL;
  }
  else {
    oro = build_from_dclusterbasis_dclusteroperator(G->rb);
    weight_dclusterbasis_dclusteroperator(G->rb, oro);
  }
  if (cbortho == true) {
    oco = NULL;
  }
  else {
    oco = build_from_dclusterbasis_dclusteroperator(G->cb);
    weight_dclusterbasis_dclusteroperator(G->cb, oco);
  }

  yco = build_from_dclusterbasis_dclusteroperator(G->cb);
  yro = build_from_dclusterbasis_dclusteroperator(G->rb);

  compress_dweights(G, oro, oco, yro, yco, tm);

  /* First copy S^{*} / S or compute  R_s * S^{*} / R_t * S  
     lets orientate by the h2compression and call it local weigths */

  if (cbortho == false) {
    del_dclusteroperator(oco);
  }
  if (rbortho == false) {
    del_dclusteroperator(oro);
  }

  /* Truncation */
  /* Column */
  cb = clone_structure_dclusterbasis(G->cb);
  bcc = build_from_dclusterbasis_dclusteroperator(G->cb);
  truncate_dclusterbasis(G->cb, cb, yco, bcc, tm, eps);
  del_dclusteroperator(yco);

  /* Row */
  rb = clone_structure_dclusterbasis(G->rb);
  bcr = build_from_dclusterbasis_dclusteroperator(G->rb);
  truncate_dclusterbasis(G->rb, rb, yro, bcr, tm, eps);
  del_dclusteroperator(yro);

  /* Final projection */
  New_G = build_projected_dh2matrix(G, rb, cb, bcr, bcc);

  del_dclusteroperator(bcc);
  del_dclusteroperator(bcr);

  return New_G;
}
