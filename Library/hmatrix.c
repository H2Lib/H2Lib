/* ------------------------------------------------------------
 * This is the file "hmatrix.c" of the H2Lib package.
 * All rights reserved, Steffen Boerm 2010
 * ------------------------------------------------------------ */

#include <stdio.h>

#ifdef USE_NETCDF
#include <string.h>
#include <netcdf.h>
#endif

#include "hmatrix.h"
#include "basic.h"

/* ------------------------------------------------------------
 * Constructors and destructors
 * ------------------------------------------------------------ */

phmatrix
init_hmatrix(phmatrix hm, pccluster rc, pccluster cc)
{
  hm->rc = rc;
  hm->cc = cc;

  hm->r = NULL;
  hm->f = NULL;

  hm->son = NULL;
  hm->rsons = 0;
  hm->csons = 0;

  hm->refs = 0;
  hm->desc = 0;

  return hm;
}

void
uninit_hmatrix(phmatrix hm)
{
  uint      rsons = hm->rsons;
  uint      csons = hm->csons;
  uint      i, j;

  assert(hm->refs == 0);

  if (hm->son) {
    for (j = 0; j < csons; j++)
      for (i = 0; i < rsons; i++)
	unref_hmatrix(hm->son[i + j * rsons]);
    freemem(hm->son);
  }

  if (hm->f)
    del_amatrix(hm->f);

  if (hm->r)
    del_rkmatrix(hm->r);
}

phmatrix
new_hmatrix(pccluster rc, pccluster cc)
{
  phmatrix  hm;

  hm = allocmem(sizeof(hmatrix));

  init_hmatrix(hm, rc, cc);

  return hm;
}

phmatrix
new_rk_hmatrix(pccluster rc, pccluster cc, uint k)
{
  phmatrix  hm;

  hm = new_hmatrix(rc, cc);

  hm->r = new_rkmatrix(rc->size, cc->size, k);

  hm->desc = 1;

  return hm;
}

phmatrix
new_full_hmatrix(pccluster rc, pccluster cc)
{
  phmatrix  hm;

  hm = new_hmatrix(rc, cc);

  hm->f = new_amatrix(rc->size, cc->size);

  hm->desc = 1;

  return hm;
}

phmatrix
new_super_hmatrix(pccluster rc, pccluster cc, uint rsons, uint csons)
{
  phmatrix  hm;
  uint      i, j;

  hm = new_hmatrix(rc, cc);

  hm->rsons = rsons;
  hm->csons = csons;

  hm->son = (phmatrix *) allocmem((size_t) sizeof(phmatrix) * rsons * csons);
  for (j = 0; j < csons; j++)
    for (i = 0; i < rsons; i++)
      hm->son[i + j * rsons] = NULL;

  return hm;
}

phmatrix
clonestructure_hmatrix(pchmatrix src)
{
  const uint rsons = src->rsons;
  const uint csons = src->csons;

  phmatrix  hm1;
  uint      i, j;

  phmatrix  hm;

  if (src->son != NULL) {
    hm = new_super_hmatrix(src->rc, src->cc, rsons, csons);
    for (j = 0; j < csons; ++j) {
      for (i = 0; i < rsons; ++i) {
	hm1 = clonestructure_hmatrix(src->son[i + j * rsons]);
	ref_hmatrix(hm->son + i + j * rsons, hm1);
      }
    }
  }
  else if (src->r != NULL) {
    hm = new_rk_hmatrix(src->rc, src->cc, src->r->k);
  }
  else {
    assert(src->f != NULL);
    hm = new_full_hmatrix(src->rc, src->cc);
  }

  update_hmatrix(hm);

  return hm;
}

phmatrix
clone_hmatrix(pchmatrix src)
{
  const uint rsons = src->rsons;
  const uint csons = src->csons;

  phmatrix  hm1;
  uint      i, j;

  phmatrix  hm;

  if (src->son != NULL) {
    hm = new_super_hmatrix(src->rc, src->cc, rsons, csons);
    for (j = 0; j < csons; ++j) {
      for (i = 0; i < rsons; ++i) {
	hm1 = clone_hmatrix(src->son[i + j * rsons]);
	ref_hmatrix(hm->son + i + j * rsons, hm1);
      }
    }
    update_hmatrix(hm);
  }
  else if (src->r != NULL) {
    hm = new_rk_hmatrix(src->rc, src->cc, src->r->k);
    copy_amatrix(false, &src->r->A, &hm->r->A);
    copy_amatrix(false, &src->r->B, &hm->r->B);
  }
  else {
    assert(src->f != NULL);
    hm = new_full_hmatrix(src->rc, src->cc);
    copy_amatrix(false, src->f, hm->f);
  }

  update_hmatrix(hm);

  return hm;
}

void
update_hmatrix(phmatrix hm)
{
  uint      desc;
  uint      rsons, csons;
  uint      i, j;

  desc = 1;

  if (hm->son) {
    rsons = hm->rsons;
    csons = hm->csons;

    for (j = 0; j < csons; j++)
      for (i = 0; i < rsons; i++)
	desc += hm->son[i + j * rsons]->desc;
  }

  hm->desc = desc;
}

void
del_hmatrix(phmatrix hm)
{
  uninit_hmatrix(hm);

  freemem(hm);
}

/* ------------------------------------------------------------
 * Reference counting
 * ------------------------------------------------------------ */

void
ref_hmatrix(phmatrix *ptr, phmatrix hm)
{
  if (*ptr)
    unref_hmatrix(*ptr);

  *ptr = hm;

  if (hm)
    hm->refs++;
}

void
unref_hmatrix(phmatrix hm)
{
  assert(hm->refs > 0);

  hm->refs--;

  if (hm->refs == 0)
    del_hmatrix(hm);
}

/* ------------------------------------------------------------
 * Statistics
 * ------------------------------------------------------------ */

size_t
getsize_hmatrix(pchmatrix hm)
{
  size_t    sz;
  uint      rsons = hm->rsons;
  uint      csons = hm->csons;
  uint      i, j;

  sz = (size_t) sizeof(hmatrix);

  if (hm->r)
    sz += getsize_rkmatrix(hm->r);

  if (hm->f)
    sz += getsize_amatrix(hm->f);

  for (j = 0; j < csons; j++)
    for (i = 0; i < rsons; i++)
      sz += getsize_hmatrix(hm->son[i + j * rsons]);

  return sz;
}

size_t
getnearsize_hmatrix(pchmatrix hm)
{
  size_t    sz;
  uint      rsons = hm->rsons;
  uint      csons = hm->csons;
  uint      i, j;

  sz = 0;

  if (hm->f)
    sz += getsize_amatrix(hm->f);

  for (j = 0; j < csons; j++)
    for (i = 0; i < rsons; i++)
      sz += getnearsize_hmatrix(hm->son[i + j * rsons]);

  return sz;
}

size_t
getfarsize_hmatrix(pchmatrix hm)
{
  size_t    sz;
  uint      rsons = hm->rsons;
  uint      csons = hm->csons;
  uint      i, j;

  sz = 0;

  if (hm->r)
    sz += getsize_rkmatrix(hm->r);

  for (j = 0; j < csons; j++)
    for (i = 0; i < rsons; i++)
      sz += getfarsize_hmatrix(hm->son[i + j * rsons]);

  return sz;
}

/* ------------------------------------------------------------
 * Simple utility functions
 * ------------------------------------------------------------ */

void
clear_hmatrix(phmatrix hm)
{
  uint      rsons, csons;
  uint      i, j;

  if (hm->son) {
    rsons = hm->rsons;
    csons = hm->csons;

    for (j = 0; j < csons; j++) {
      for (i = 0; i < rsons; i++) {
	clear_hmatrix(hm->son[i + j * rsons]);
      }
    }
  }
  else if (hm->r) {
    setrank_rkmatrix(hm->r, 0);
  }
  else {
    assert(hm->f);
    clear_amatrix(hm->f);
  }
}

void
clear_upper_hmatrix(phmatrix hm, bool strict)
{
  uint      rsons, csons;
  uint      i, j;

  if (hm->son) {
    rsons = hm->rsons;
    csons = hm->csons;

    for (j = 0; j < csons; j++) {
      for (i = 0; i < UINT_MIN(j, rsons); i++) {
	clear_hmatrix(hm->son[i + j * rsons]);
      }
      if (j < rsons) {
	clear_upper_hmatrix(hm->son[j + j * rsons], strict);
      }
    }
  }
  else if (hm->r) {
    setrank_rkmatrix(hm->r, 0);
  }
  else {
    assert(hm->f);
    assert(hm->rc == hm->cc);
    clear_upper_amatrix(hm->f, strict);
  }
}

void
copy_hmatrix(pchmatrix src, phmatrix trg)
{
  const uint rsons = src->rsons;
  const uint csons = src->csons;

  uint      i, j;

  assert(rsons == trg->rsons);
  assert(csons == trg->csons);
  assert(src->rc == trg->rc);
  assert(src->cc == trg->cc);

  if (src->son != NULL) {
    for (j = 0; j < csons; ++j) {
      for (i = 0; i < rsons; ++i) {
	copy_hmatrix(src->son[i + j * rsons], trg->son[i + j * rsons]);
      }
    }
  }
  else if (src->r != NULL) {
    resize_rkmatrix(trg->r, src->rc->size, src->cc->size, src->r->k);
    copy_amatrix(false, &src->r->A, &trg->r->A);
    copy_amatrix(false, &src->r->B, &trg->r->B);
  }
  else if (src->f != NULL) {
    copy_amatrix(false, src->f, trg->f);
  }
}

void
identity_hmatrix(phmatrix hm)
{
  uint      rsons, csons;
  uint      i, j;

  if (hm->son) {
    rsons = hm->rsons;
    csons = hm->csons;

    for (j = 0; j < csons; j++) {
      for (i = 0; i < rsons; i++) {
	identity_hmatrix(hm->son[i + j * rsons]);
      }
    }
  }
  else if (hm->r) {
    setrank_rkmatrix(hm->r, 0);
  }
  else {
    assert(hm->f);
    if (hm->rc == hm->cc) {
      identity_amatrix(hm->f);
    }
    else {
      clear_amatrix(hm->f);
    }
  }
}

void
random_hmatrix(phmatrix hm, uint kmax)
{
  uint      rsons = hm->rsons;
  uint      csons = hm->csons;

  uint      i, j;

  if (hm->son != NULL) {
    for (j = 0; j < csons; ++j) {
      for (i = 0; i < rsons; ++i) {
	random_hmatrix(hm->son[i + j * rsons], kmax);
      }
    }
  }
  else if (hm->r != NULL) {
    resize_rkmatrix(hm->r, hm->rc->size, hm->cc->size, kmax);
    random_amatrix(&hm->r->A);
    random_amatrix(&hm->r->B);
  }
  else {
    assert(hm->f != NULL);
    random_amatrix(hm->f);
  }
}

/* ------------------------------------------------------------
 * Build H-matrix based on block tree
 * ------------------------------------------------------------ */

phmatrix
build_from_block_hmatrix(pcblock b, uint k)
{
  phmatrix  h, h1;
  pcblock   b1;
  int       rsons, csons;
  int       i, j;

  h = NULL;

  if (b->son) {
    rsons = b->rsons;
    csons = b->csons;

    h = new_super_hmatrix(b->rc, b->cc, rsons, csons);

    for (j = 0; j < csons; j++) {
      for (i = 0; i < rsons; i++) {
	b1 = b->son[i + j * rsons];

	h1 = build_from_block_hmatrix(b1, k);

	ref_hmatrix(h->son + i + j * rsons, h1);
      }
    }
  }
  else if (b->a > 0)
    h = new_rk_hmatrix(b->rc, b->cc, k);
  else
    h = new_full_hmatrix(b->rc, b->cc);

  update_hmatrix(h);

  return h;
}

/* ------------------------------------------------------------
 * Build block tree from H-matrix
 * ------------------------------------------------------------ */

pblock
build_from_hmatrix_block(pchmatrix G)
{
  uint      i, j, rsons, csons;
  pblock    b;

  b = 0;
  if (G->son) {
    rsons = G->rsons;
    csons = G->csons;

    b = new_block((pcluster) G->rc, (pcluster) G->cc, false, rsons, csons);

    for (j = 0; j < csons; j++) {
      for (i = 0; i < rsons; i++) {
	b->son[i + j * rsons] =
	  build_from_hmatrix_block(G->son[i + j * rsons]);
      }
    }
  }
  else if (G->r) {
    b = new_block((pcluster) G->rc, (pcluster) G->cc, true, 0, 0);
  }
  else {
    assert(G->f);
    b = new_block((pcluster) G->rc, (pcluster) G->cc, false, 0, 0);
  }

  update_block(b);

  return b;
}

/* ------------------------------------------------------------
 * Matrix-vector multiplication
 * ------------------------------------------------------------ */

void
mvm_hmatrix_avector(field alpha, bool atrans, pchmatrix a, pcavector x,
		    pavector y)
{
  if (atrans)
    addevaltrans_hmatrix_avector(alpha, a, x, y);
  else
    addeval_hmatrix_avector(alpha, a, x, y);
}

void
fastaddeval_hmatrix_avector(field alpha, pchmatrix hm, pcavector x,
			    pavector y)
{
  pavector  x1, y1;
  avector   xtmp, ytmp;
  uint      rsons, csons;
  uint      xoff, yoff, i, j;

  assert(x->dim == hm->cc->size);
  assert(y->dim == hm->rc->size);

  if (hm->r) {
    addeval_rkmatrix_avector(alpha, hm->r, x, y);
  }
  else if (hm->f) {
    mvm_amatrix_avector(alpha, false, hm->f, x, y);
  }
  else {
    rsons = hm->rsons;
    csons = hm->csons;

    xoff = 0;
    for (j = 0; j < csons; j++) {
      x1 = init_sub_avector(&xtmp, (pavector) x, hm->son[j * rsons]->cc->size,
			    xoff);

      yoff = 0;
      for (i = 0; i < rsons; i++) {
	y1 = init_sub_avector(&ytmp, y, hm->son[i]->rc->size, yoff);

	fastaddeval_hmatrix_avector(alpha, hm->son[i + j * rsons], x1, y1);

	uninit_avector(y1);

	yoff += hm->son[i]->rc->size;
      }
      assert(yoff == hm->rc->size);

      uninit_avector(x1);

      xoff += hm->son[j * rsons]->cc->size;
    }
    assert(xoff == hm->cc->size);
  }
}

void
addeval_hmatrix_avector(field alpha, pchmatrix hm, pcavector x, pavector y)
{
  pavector  xp, yp;
  avector   xtmp, ytmp;
  uint      i, ip;

  assert(x->dim == hm->cc->size);
  assert(y->dim == hm->rc->size);

  /* Permutation of x */
  xp = init_avector(&xtmp, x->dim);
  for (i = 0; i < xp->dim; i++) {
    ip = hm->cc->idx[i];
    assert(ip < x->dim);
    xp->v[i] = x->v[ip];
  }

  /* Permutation of y */
  yp = init_avector(&ytmp, y->dim);
  for (i = 0; i < yp->dim; i++) {
    ip = hm->rc->idx[i];
    assert(ip < y->dim);
    yp->v[i] = y->v[ip];
  }

  /* Matrix-vector multiplication */
  fastaddeval_hmatrix_avector(alpha, hm, xp, yp);

  /* Reverse permutation of y */
  for (i = 0; i < yp->dim; i++) {
    ip = hm->rc->idx[i];
    assert(ip < y->dim);
    y->v[ip] = yp->v[i];
  }

  uninit_avector(yp);
  uninit_avector(xp);
}

void
fastaddevaltrans_hmatrix_avector(field alpha, pchmatrix hm, pcavector x,
				 pavector y)
{
  pavector  x1, y1;
  avector   xtmp, ytmp;
  uint      rsons, csons;
  uint      xoff, yoff, i, j;

  assert(x->dim == hm->rc->size);
  assert(y->dim == hm->cc->size);

  if (hm->r) {
    addevaltrans_rkmatrix_avector(alpha, hm->r, x, y);
  }
  else if (hm->f) {
    mvm_amatrix_avector(alpha, true, hm->f, x, y);
  }
  else {
    rsons = hm->rsons;
    csons = hm->csons;

    yoff = 0;
    for (j = 0; j < csons; j++) {
      y1 = init_sub_avector(&ytmp, y, hm->son[j * rsons]->cc->size, yoff);

      xoff = 0;
      for (i = 0; i < rsons; i++) {
	x1 =
	  init_sub_avector(&xtmp, (pavector) x, hm->son[i]->rc->size, xoff);

	fastaddevaltrans_hmatrix_avector(alpha, hm->son[i + j * rsons], x1,
					 y1);

	uninit_avector(x1);

	xoff += hm->son[i]->rc->size;
      }
      assert(xoff == hm->rc->size);

      uninit_avector(y1);

      yoff += hm->son[j * rsons]->cc->size;
    }
    assert(yoff == hm->cc->size);
  }
}

void
addevaltrans_hmatrix_avector(field alpha, pchmatrix hm, pcavector x,
			     pavector y)
{
  pavector  xp, yp;
  avector   xtmp, ytmp;
  uint      i, ip;

  assert(x->dim == hm->rc->size);
  assert(y->dim == hm->cc->size);

  /* Permutation of x */
  xp = init_avector(&xtmp, x->dim);
  for (i = 0; i < xp->dim; i++) {
    ip = hm->rc->idx[i];
    assert(ip < x->dim);
    xp->v[i] = x->v[ip];
  }

  /* Permutation of y */
  yp = init_avector(&ytmp, y->dim);
  for (i = 0; i < yp->dim; i++) {
    ip = hm->cc->idx[i];
    assert(ip < y->dim);
    yp->v[i] = y->v[ip];
  }

  /* Matrix-vector multiplication */
  fastaddevaltrans_hmatrix_avector(alpha, hm, xp, yp);

  /* Reverse permutation of y */
  for (i = 0; i < yp->dim; i++) {
    ip = hm->cc->idx[i];
    assert(ip < y->dim);
    y->v[ip] = yp->v[i];
  }

  uninit_avector(yp);
  uninit_avector(xp);
}

static void
addevalsymm_offdiag(field alpha, pchmatrix hm, uint roff, uint coff,
		    pcavector xp, pavector yp)
{
  avector   tmp1, tmp2;
  pavector  xp1, yp1;
  uint      rsons, csons;
  uint      roff1, coff1;
  uint      i, j;

  if (hm->f) {
    xp1 = init_sub_avector(&tmp1, (pavector) xp, hm->cc->size, coff);
    yp1 = init_sub_avector(&tmp2, yp, hm->rc->size, roff);

    addeval_amatrix_avector(alpha, hm->f, xp1, yp1);

    uninit_avector(yp1);
    uninit_avector(xp1);

    xp1 = init_sub_avector(&tmp1, (pavector) xp, hm->rc->size, roff);
    yp1 = init_sub_avector(&tmp2, yp, hm->cc->size, coff);

    addevaltrans_amatrix_avector(alpha, hm->f, xp1, yp1);

    uninit_avector(yp1);
    uninit_avector(xp1);
  }
  else if (hm->r) {
    xp1 = init_sub_avector(&tmp1, (pavector) xp, hm->cc->size, coff);
    yp1 = init_sub_avector(&tmp2, yp, hm->rc->size, roff);

    addeval_rkmatrix_avector(alpha, hm->r, xp1, yp1);

    uninit_avector(yp1);
    uninit_avector(xp1);

    xp1 = init_sub_avector(&tmp1, (pavector) xp, hm->rc->size, roff);
    yp1 = init_sub_avector(&tmp2, yp, hm->cc->size, coff);

    addevaltrans_rkmatrix_avector(alpha, hm->r, xp1, yp1);

    uninit_avector(yp1);
    uninit_avector(xp1);
  }
  else {
    assert(hm->son != 0);

    rsons = hm->rsons;
    csons = hm->csons;

    coff1 = coff;
    for (j = 0; j < csons; j++) {
      roff1 = roff;

      for (i = 0; i < rsons; i++) {
	addevalsymm_offdiag(alpha, hm->son[i + j * rsons], roff1, coff1, xp,
			    yp);

	roff1 += hm->son[i]->rc->size;
      }
      assert(roff1 == roff + hm->rc->size);

      coff1 += hm->son[j * rsons]->cc->size;
    }
    assert(coff1 == coff + hm->cc->size);
  }
}

static void
addevalsymm_diag(field alpha, pchmatrix hm, uint off, pcavector xp,
		 pavector yp)
{
  avector   tmp1, tmp2;
  pavector  xp1, yp1;
  pfield    aa;
  uint      lda, sons;
  uint      roff, coff;
  uint      n;
  uint      i, j;

  assert(hm->rc == hm->cc);

  if (hm->f) {
    assert(hm->rc->size == hm->f->rows);
    assert(hm->cc->size == hm->f->cols);

    aa = hm->f->a;
    lda = hm->f->ld;

    n = hm->rc->size;
    xp1 = init_sub_avector(&tmp1, (pavector) xp, n, off);
    yp1 = init_sub_avector(&tmp2, yp, n, off);

    for (j = 0; j < n; j++) {
      yp1->v[j] += alpha * aa[j + j * lda] * xp1->v[j];
      for (i = j + 1; i < n; i++) {
	yp1->v[i] += alpha * aa[i + j * lda] * xp1->v[j];
	yp1->v[j] += alpha * CONJ(aa[i + j * lda]) * xp1->v[i];
      }
    }

    uninit_avector(yp1);
    uninit_avector(xp1);
  }
  else {
    assert(hm->son != 0);
    assert(hm->rsons == hm->csons);

    sons = hm->rsons;

    coff = off;
    for (j = 0; j < sons; j++) {
      roff = coff;

      addevalsymm_diag(alpha, hm->son[j + j * sons], coff, xp, yp);

      roff += hm->rc->son[j]->size;
      for (i = j + 1; i < sons; i++) {
	addevalsymm_offdiag(alpha, hm->son[i + j * sons], roff, coff, xp, yp);

	roff += hm->son[i]->rc->size;
      }
      assert(roff == off + hm->rc->size);

      coff += hm->son[j * sons]->cc->size;
    }
    assert(coff == off + hm->rc->size);
  }
}

void
fastaddevalsymm_hmatrix_avector(field alpha, pchmatrix hm, pcavector xp,
				pavector yp)
{
  assert(hm->rc == hm->cc);

  addevalsymm_diag(alpha, hm, 0, xp, yp);
}

void
addevalsymm_hmatrix_avector(field alpha, pchmatrix hm, pcavector x,
			    pavector y)
{
  pavector  xp, yp;
  avector   xtmp, ytmp;
  uint      n;
  const uint *idx;
  uint      i, ip;

  assert(hm->rc == hm->cc);
  assert(x->dim == hm->cc->size);
  assert(y->dim == hm->rc->size);

  n = hm->rc->size;
  idx = hm->rc->idx;

  /* Permutation of x */
  xp = init_avector(&xtmp, n);
  for (i = 0; i < n; i++) {
    ip = idx[i];
    assert(ip < n);
    xp->v[i] = x->v[ip];
  }

  /* Permutation of y */
  yp = init_avector(&ytmp, n);
  for (i = 0; i < n; i++) {
    ip = idx[i];
    assert(ip < n);
    yp->v[i] = y->v[ip];
  }

  /* Matrix-vector multiplication */
  fastaddevalsymm_hmatrix_avector(alpha, hm, xp, yp);

  /* Reverse permutation of y */
  for (i = 0; i < n; i++) {
    ip = idx[i];
    assert(ip < n);
    y->v[ip] = yp->v[i];
  }

  /* Clean up */
  uninit_avector(yp);
  uninit_avector(xp);
}

/* ------------------------------------------------------------
 * Enumeration
 * ------------------------------------------------------------ */

static void
enumerate(pcblock b, uint bname, phmatrix hm, phmatrix *hn)
{
  uint      bname1;
  uint      i, j;

  assert(hm->rc == b->rc);
  assert(hm->cc == b->cc);

  hn[bname] = hm;

  bname1 = bname + 1;

  if (hm == 0 || hm->son == 0)
    for (j = 0; j < b->csons; j++)
      for (i = 0; i < b->rsons; i++) {
	enumerate(b->son[i + j * b->rsons], bname1, 0, hn);

	bname1 += b->son[i + j * b->rsons]->desc;
      }
  else {
    assert(b->rsons == hm->rsons);
    assert(b->csons == hm->csons);

    for (j = 0; j < b->csons; j++)
      for (i = 0; i < b->rsons; i++) {
	enumerate(b->son[i + j * b->rsons], bname1, hm->son[i + j * b->rsons],
		  hn);

	bname1 += b->son[i + j * b->rsons]->desc;
      }
  }
  assert(bname1 == bname + b->desc);
}

phmatrix *
enumerate_hmatrix(pcblock b, phmatrix hm)
{
  phmatrix *hn;

  hn = (phmatrix *) allocmem((size_t) sizeof(phmatrix) * b->desc);

  enumerate(b, 0, hm, hn);

  return hn;
}

/* ------------------------------------------------------------
 * Simple utility functions
 * ------------------------------------------------------------ */

real
norm2_hmatrix(pchmatrix H)
{
  return norm2_matrix((mvm_t) mvm_hmatrix_avector, (void *) H, H->rc->size,
		      H->cc->size);
}

real
norm2diff_hmatrix(pchmatrix a, pchmatrix b)
{
  return norm2diff_matrix((mvm_t) mvm_hmatrix_avector, (void *) a,
			  (mvm_t) mvm_hmatrix_avector, (void *) b,
			  a->rc->size, a->cc->size);
}

/* ------------------------------------------------------------
 * File I/O
 * ------------------------------------------------------------ */

#ifdef USE_NETCDF
static void
write_count(pchmatrix G, size_t * blocks, size_t * coeffs)
{
  uint      rsons, csons;
  uint      i, j;

  /* Increase submatrix counter */
  (*blocks)++;

  if (G->f) {
    /* If it's an amatrix: add number of coefficients to coeffs */
    (*coeffs) += G->f->rows * G->f->cols;
  }
  else if (G->r) {
    /* If it's an rkmatrix: also add number of coefficients to coeffs */
    (*coeffs) +=
      (getrows_rkmatrix(G->r) +
       getcols_rkmatrix(G->r)) * getrank_rkmatrix(G->r);
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
write_cdf(pchmatrix G, size_t blocks, size_t coeffs,
	  size_t * blockidx, size_t * coeffidx,
	  int nc_file,
	  int nc_type, int nc_rows, int nc_cols, int nc_rank,
	  int nc_rsons, int nc_csons, int nc_coeff)
{
  size_t    start, count;
  ptrdiff_t stride;
  pamatrix  f;
  prkmatrix r;
  uint      rsons, csons;
  uint      rows, cols, k;
  uint      lda, ldb, ldf;
  int       val;
  uint      i, j;

  assert(*blockidx <= blocks);

  rows = getrows_hmatrix(G);
  cols = getcols_hmatrix(G);

  /* Write number of rows to nc_rows[*blockidx] */
  start = *blockidx;
  count = 1;
  stride = 1;
  val = rows;
  nc_put_vars(nc_file, nc_rows, &start, &count, &stride, &val);

  /* Write number of columns to nc_cols[*blockidx] */
  val = cols;
  nc_put_vars(nc_file, nc_cols, &start, &count, &stride, &val);

  if (G->f) {
    f = G->f;
    ldf = f->ld;
    assert(f->rows == rows);
    assert(f->cols == cols);

    /* Write type 1 to nc_type[*blockidx] */
    val = 1;
    nc_put_vars(nc_file, nc_type, &start, &count, &stride, &val);

    /* An amatrix has no rank, rsons, or csons */
    val = 0;
    nc_put_vars(nc_file, nc_rank, &start, &count, &stride, &val);
    nc_put_vars(nc_file, nc_rsons, &start, &count, &stride, &val);
    nc_put_vars(nc_file, nc_csons, &start, &count, &stride, &val);

    /* Write coefficients */
    start = *coeffidx;
    count = rows;
    assert(start + rows * cols <= coeffs);
    for (j = 0; j < cols; j++) {
      nc_put_vars(nc_file, nc_coeff, &start, &count, &stride, f->a + j * ldf);
      start += rows;
    }
    *coeffidx = start;

    /* Increase blocks counter */
    (*blockidx)++;
  }
  else if (G->r) {
    r = G->r;
    lda = r->A.ld;
    ldb = r->B.ld;
    k = r->A.cols;
    assert(r->B.cols == k);
    assert(r->A.rows == rows);
    assert(r->B.rows == cols);

    /* Write type 2 to nc_type[*blockidx] */
    val = 2;
    nc_put_vars(nc_file, nc_type, &start, &count, &stride, &val);

    /* Write rank to nc_rank[*blockidx] */
    val = k;
    nc_put_vars(nc_file, nc_rank, &start, &count, &stride, &val);

    /* An rkmatrix has no rsons, or csons */
    val = 0;
    nc_put_vars(nc_file, nc_rsons, &start, &count, &stride, &val);
    nc_put_vars(nc_file, nc_csons, &start, &count, &stride, &val);

    /* Write coefficients of A */
    start = *coeffidx;
    count = rows;
    assert(start + rows * k <= coeffs);
    for (j = 0; j < k; j++) {
      nc_put_vars(nc_file, nc_coeff, &start, &count, &stride,
		  r->A.a + j * lda);
      start += rows;
    }

    /* Write coefficients of B */
    count = cols;
    assert(start + cols * k <= coeffs);
    for (j = 0; j < k; j++) {
      nc_put_vars(nc_file, nc_coeff, &start, &count, &stride,
		  r->B.a + j * ldb);
      start += cols;
    }
    *coeffidx = start;

    /* Increase submatrix counter */
    (*blockidx)++;
  }
  else if (G->son) {
    rsons = G->rsons;
    csons = G->csons;

    /* Write type 3 to nc_type */
    val = 3;
    nc_put_vars(nc_file, nc_type, &start, &count, &stride, &val);

    /* Subdivided matrices have no rank */
    val = 0;
    nc_put_vars(nc_file, nc_rank, &start, &count, &stride, &val);

    /* Write rsons and csons.
     * rsons == 0 means that the sons use the same row cluster
     * as the father. The same convention is used for csons. */
    val = (G->rc == G->son[0]->rc ? 0 : rsons);
    nc_put_vars(nc_file, nc_rsons, &start, &count, &stride, &val);
    val = (G->cc == G->son[0]->cc ? 0 : csons);
    nc_put_vars(nc_file, nc_csons, &start, &count, &stride, &val);

    /* Increase submatrix counter */
    (*blockidx)++;

    /* Take care of sons */
    for (j = 0; j < csons; j++)
      for (i = 0; i < rsons; i++)
	write_cdf(G->son[i + j * rsons], blocks, coeffs,
		  blockidx, coeffidx,
		  nc_file,
		  nc_type, nc_rows, nc_cols, nc_rank,
		  nc_rsons, nc_csons, nc_coeff);
  }
}

void
write_cdf_hmatrix(pchmatrix G, const char *name)
{
  size_t    blocks, blockidx;
  size_t    coeffs, coeffidx;
  int       nc_file, nc_type, nc_rows, nc_cols, nc_rank;
  int       nc_rsons, nc_csons, nc_coeff, nc_blocks, nc_coeffs;
  int       result;

  /* Count number of blocks and coefficients */
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

  /* Define "rows" variable */
  result = nc_def_var(nc_file, "rows", NC_INT, 1, &nc_blocks, &nc_rows);
  assert(result == NC_NOERR);

  /* Define "cols" variable */
  result = nc_def_var(nc_file, "cols", NC_INT, 1, &nc_blocks, &nc_cols);
  assert(result == NC_NOERR);

  /* Define "rank" variable */
  result = nc_def_var(nc_file, "rank", NC_INT, 1, &nc_blocks, &nc_rank);
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
	    nc_file, nc_type, nc_rows, nc_cols, nc_rank,
	    nc_rsons, nc_csons, nc_coeff);
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
write_cdfpart_hmatrix(pchmatrix G, int nc_file, const char *prefix)
{
  size_t    blocks, blockidx;
  size_t    coeffs, coeffidx;
  char     *buf;
  int       bufsize;
  int       nc_type, nc_rows, nc_cols, nc_rank;
  int       nc_rsons, nc_csons, nc_coeff, nc_blocks, nc_coeffs;
  int       result;

  /* Prepare buffer for prefixed names */
  bufsize = strlen(prefix) + 16;
  buf = (char *) allocmem(sizeof(char) * bufsize);

  /* Count number of blocks and coefficients */
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

  /* Define "rows" variable */
  prefix_name(buf, bufsize, prefix, "rows");
  result = nc_def_var(nc_file, buf, NC_INT, 1, &nc_blocks, &nc_rows);
  assert(result == NC_NOERR);

  /* Define "cols" variable */
  prefix_name(buf, bufsize, prefix, "cols");
  result = nc_def_var(nc_file, buf, NC_INT, 1, &nc_blocks, &nc_cols);
  assert(result == NC_NOERR);

  /* Define "rank" variable */
  prefix_name(buf, bufsize, prefix, "rank");
  result = nc_def_var(nc_file, buf, NC_INT, 1, &nc_blocks, &nc_rank);
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
	    nc_file, nc_type, nc_rows, nc_cols, nc_rank,
	    nc_rsons, nc_csons, nc_coeff);
  assert(blockidx == blocks);
  assert(coeffidx == coeffs);

  /* Clean up */
  nc_sync(nc_file);
  freemem(buf);
}

static phmatrix
read_cdf_part(int nc_file,
	      int nc_type, int nc_rows, int nc_cols, int nc_rank,
	      int nc_rsons, int nc_csons, int nc_coeffs,
	      pcluster *rc, pcluster *cc,
	      uint * ridx, uint * cidx, size_t * blockidx, size_t * coeffidx)
{
  phmatrix  G, G1;
  prkmatrix r;
  pamatrix  f;
  pcluster *rc1, *cc1;
  uint     *ridx1, *cidx1;
  uint      rows, cols;
  uint      rsons, csons;
  uint      lda, ldb, ldf;
  uint      k;
  uint      i, j;
  size_t    start, count;
  ptrdiff_t stride;
  int       val, result;

  /* Get number of rows */
  start = *blockidx;
  count = 1;
  stride = 1;
  result = nc_get_vars(nc_file, nc_rows, &start, &count, &stride, &val);
  assert(result == NC_NOERR);
  rows = val;

  /* Get number of columns */
  result = nc_get_vars(nc_file, nc_cols, &start, &count, &stride, &val);
  assert(result == NC_NOERR);
  cols = val;

  /* If we have not inherited a row cluster, create one */
  if (*rc == 0) {
    /* If we have not inherited a row index, create one */
    if (ridx == 0) {
      ridx = (uint *) allocmem(sizeof(uint) * rows);
      for (i = 0; i < rows; i++)
	ridx[i] = i;
    }

    *rc = new_cluster(rows, ridx, 0, 1);
  }
  else {
    assert((*rc)->size == rows);

    ridx = (*rc)->idx;
  }

  /* If we have not inherited a column cluster, create one */
  if (*cc == 0) {
    /* If we have not inherited a column index, create one */
    if (cidx == 0) {
      cidx = (uint *) allocmem(sizeof(uint) * cols);
      for (j = 0; j < cols; j++)
	cidx[j] = j;
    }

    *cc = new_cluster(cols, cidx, 0, 1);
  }
  else {
    assert((*cc)->size == cols);

    cidx = (*cc)->idx;
  }

  /* Get matrix type */
  result = nc_get_vars(nc_file, nc_type, &start, &count, &stride, &val);
  assert(result == NC_NOERR);

  G = 0;
  switch (val) {
  case 1:
    /* Create amatrix */
    G = new_full_hmatrix(*rc, *cc);
    f = G->f;
    ldf = f->ld;

    /* Get coefficients */
    start = *coeffidx;
    count = rows;
    for (j = 0; j < cols; j++) {
      nc_get_vars(nc_file, nc_coeffs, &start, &count, &stride,
		  f->a + j * ldf);
      start += rows;
    }
    *coeffidx = start;

    /* Switch to next submatrix */
    (*blockidx)++;

    break;

  case 2:
    /* Get rank */
    result = nc_get_vars(nc_file, nc_rank, &start, &count, &stride, &val);
    assert(result == NC_NOERR);
    k = val;

    /* Create matrix */
    G = new_rk_hmatrix(*rc, *cc, k);
    r = G->r;
    lda = r->A.ld;
    ldb = r->B.ld;

    /* Get matrix A */
    start = *coeffidx;
    count = rows;
    for (j = 0; j < k; j++) {
      nc_get_vars(nc_file, nc_coeffs, &start, &count, &stride,
		  r->A.a + j * lda);
      start += rows;
    }

    /* Get matrix B */
    count = cols;
    for (j = 0; j < k; j++) {
      nc_get_vars(nc_file, nc_coeffs, &start, &count, &stride,
		  r->B.a + j * ldb);
      start += cols;
    }
    *coeffidx = start;

    /* Switch to next submatrix */
    (*blockidx)++;

    break;

  case 3:
    /* Get number of row sons */
    result = nc_get_vars(nc_file, nc_rsons, &start, &count, &stride, &val);
    assert(result == NC_NOERR);
    rsons = val;

    /* Get number of column sons */
    result = nc_get_vars(nc_file, nc_csons, &start, &count, &stride, &val);
    assert(result == NC_NOERR);
    csons = val;

    /* If we need the row sons, ensure that we have them */
    if (rsons > 0) {
      /* Update row cluster if it didn't have sons before */
      if ((*rc)->sons == 0)
	setsons_cluster(*rc, rsons);
      else
	assert((*rc)->sons == rsons);

      rc1 = (*rc)->son;
    }
    else {
      rsons = 1;
      rc1 = rc;
    }

    /* If we need the column sons, ensure that we have them */
    if (csons > 0) {
      /* Update column cluster if it didn't have sons before */
      if ((*cc)->sons == 0)
	setsons_cluster(*cc, csons);
      else
	assert((*cc)->sons == csons);

      cc1 = (*cc)->son;
    }
    else {
      csons = 1;
      cc1 = cc;
    }

    /* Create H-matrix */
    G = new_super_hmatrix(*rc, *cc, rsons, csons);

    /* Switch to next submatrix */
    (*blockidx)++;

    cidx1 = cidx;
    for (j = 0; j < csons; j++) {
      ridx1 = ridx;

      for (i = 0; i < rsons; i++) {
	G1 = read_cdf_part(nc_file,
			   nc_type, nc_rows, nc_cols, nc_rank,
			   nc_rsons, nc_csons, nc_coeffs,
			   rc1 + i, cc1 + j, ridx1, cidx1,
			   blockidx, coeffidx);
	ref_hmatrix(G->son + i + j * rsons, G1);

	ridx1 += G1->rc->size;
      }
      assert(ridx1 == ridx + rows);

      cidx1 += G->son[j * rsons]->cc->size;
    }
    assert(cidx1 == cidx + cols);

    break;

  default:
    (void) fprintf(stderr, "Unexpected matrix type %d\n", val);
    abort();
  }

  update_hmatrix(G);

  return G;
}

phmatrix
read_cdf_hmatrix(const char *name, pccluster rc, pccluster cc)
{
  phmatrix  G;
  pcluster  rc0, cc0;
  size_t    blocks, blockidx;
  size_t    coeffs, coeffidx;
  char      dimname[NC_MAX_NAME + 1];
  int       nc_file, nc_type, nc_rows, nc_cols, nc_rank;
  int       nc_rsons, nc_csons, nc_coeff, nc_blocks, nc_coeffs;
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

  /* Get "rows" variable */
  result = nc_inq_varid(nc_file, "rows", &nc_rows);
  assert(result == NC_NOERR);

  /* Get "cols" variable */
  result = nc_inq_varid(nc_file, "cols", &nc_cols);
  assert(result == NC_NOERR);

  /* Get "rank" variable */
  result = nc_inq_varid(nc_file, "rank", &nc_rank);
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
  rc0 = (pcluster) rc;
  cc0 = (pcluster) cc;
  blockidx = 0;
  coeffidx = 0;

  G = read_cdf_part(nc_file,
		    nc_type, nc_rows, nc_cols, nc_rank,
		    nc_rsons, nc_csons, nc_coeff,
		    &rc0, &cc0, 0, 0, &blockidx, &coeffidx);
  assert(rc0 == rc);
  assert(cc0 == cc);
  assert(blockidx == blocks);
  assert(coeffidx == coeffs);

  return G;
}

phmatrix
read_cdfpart_hmatrix(int nc_file, const char *prefix,
		     pccluster rc, pccluster cc)
{
  phmatrix  G;
  pcluster  rc0, cc0;
  size_t    blocks, blockidx;
  size_t    coeffs, coeffidx;
  char      dimname[NC_MAX_NAME + 1];
  char     *buf;
  int       bufsize;
  int       nc_type, nc_rows, nc_cols, nc_rank;
  int       nc_rsons, nc_csons, nc_coeff, nc_blocks, nc_coeffs;
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

  /* Get "rows" variable */
  prefix_name(buf, bufsize, prefix, "rows");
  result = nc_inq_varid(nc_file, buf, &nc_rows);
  assert(result == NC_NOERR);

  /* Get "cols" variable */
  prefix_name(buf, bufsize, prefix, "cols");
  result = nc_inq_varid(nc_file, buf, &nc_cols);
  assert(result == NC_NOERR);

  /* Get "rank" variable */
  prefix_name(buf, bufsize, prefix, "rank");
  result = nc_inq_varid(nc_file, buf, &nc_rank);
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
  rc0 = (pcluster) rc;
  cc0 = (pcluster) cc;
  blockidx = 0;
  coeffidx = 0;

  G = read_cdf_part(nc_file,
		    nc_type, nc_rows, nc_cols, nc_rank,
		    nc_rsons, nc_csons, nc_coeff,
		    &rc0, &cc0, 0, 0, &blockidx, &coeffidx);
  assert(rc == 0 || rc0 == rc);
  assert(cc == 0 || cc0 == cc);
  assert(blockidx == blocks);
  assert(coeffidx == coeffs);

  /* Clean up */
  freemem(buf);

  return G;
}

void
write_cdfcomplete_hmatrix(pchmatrix G, const char *name)
{
  int       nc_file;
  int       result;

  /* Create NetCDF file */
  result = nc_create(name, NC_64BIT_OFFSET, &nc_file);
  assert(result == NC_NOERR);

  /* Write row cluster tree */
  write_cdfpart_cluster(G->rc, nc_file, "rc");

  if (G->cc != G->rc) {
    /* Write column cluster tree */
    write_cdfpart_cluster(G->cc, nc_file, "cc");
  }

  /* Write hmatrix */
  write_cdfpart_hmatrix(G, nc_file, "ma");

  /* Close NetCDF file */
  nc_close(nc_file);
}

phmatrix
read_cdfcomplete_hmatrix(const char *name)
{
  pcluster  rc, cc;
  phmatrix  G;
  int       nc_file, nc_sons;
  int       result;

  /* Open NetCDF file */
  result = nc_open(name, NC_NOWRITE, &nc_file);
  assert(result == NC_NOERR);

  /* Read row cluster tree */
  rc = read_cdfpart_cluster(nc_file, "rc");

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

  /* Read hmatrix */
  G = read_cdfpart_hmatrix(nc_file, "ma", rc, cc);

  /* Close NetCDF file */
  nc_close(nc_file);

  return G;
}
#endif

static void
write_hlib_part(pchmatrix G, uint roff, uint coff, FILE * out)
{
  pcrkmatrix r;
  pcamatrix f;
  uint      rows, cols;
  uint      rsons, csons;
  uint      roff1, coff1;
  uint      lda, ldb, ldf;
  uint      i, j, k, l;

  if (G->r) {
    r = G->r;
    k = r->k;

    rows = r->A.rows;
    cols = r->B.rows;

    lda = r->A.ld;
    ldb = r->B.ld;

    assert(G->rc->size == rows);
    assert(G->cc->size == cols);
    assert(r->A.cols == k);
    assert(r->B.cols == k);

    fprintf(out, "Type 1\n"
	    "row_offset = %u\n"
	    "col_offset = %u\n"
	    "rows = %u\n"
	    "cols = %u\n"
	    "k = %u\n" "kt = %u\n", roff, coff, rows, cols, k, k);

    for (l = 0; l < k; l++)
      for (i = 0; i < rows; i++)
	fprintf(out, FIELD_CS(.16, e) "\n", FIELD_ARG(r->A.a[i + lda * l]));

    for (l = 0; l < k; l++)
      for (j = 0; j < cols; j++)
	fprintf(out, FIELD_CS(.16, e) "\n", FIELD_ARG(r->B.a[j + ldb * l]));
  }
  else if (G->f) {
    f = G->f;

    rows = f->rows;
    cols = f->cols;

    ldf = f->ld;

    fprintf(out, "Type 2\n"
	    "row_offset = %u\n"
	    "col_offset = %u\n"
	    "rows = %u\n" "cols = %u\n", roff, coff, rows, cols);

    for (j = 0; j < cols; j++)
      for (i = 0; i < rows; i++)
	fprintf(out, FIELD_CS(.16, e) "\n", FIELD_ARG(f->a[i + ldf * j]));
  }
  else {
    assert(G->son != 0);

    rows = G->rc->size;
    cols = G->cc->size;

    rsons = G->rsons;
    csons = G->csons;

    fprintf(out, "Type 3\n"
	    "row_offset = %u\n"
	    "col_offset = %u\n"
	    "rows = %u\n"
	    "cols = %u\n"
	    "block_rows = %u\n"
	    "block_cols = %u\n", roff, coff, rows, cols, rsons, csons);

    coff1 = coff;
    for (j = 0; j < csons; j++) {
      roff1 = roff;
      for (i = 0; i < rsons; i++) {
	write_hlib_part(G->son[i + j * rsons], roff1, coff1, out);

	roff1 += G->son[i]->rc->size;
      }
      assert(roff1 == roff + G->rc->size);

      coff1 += G->son[j * rsons]->cc->size;
    }
    assert(coff1 == coff + G->cc->size);
  }
}

void
write_hlib_hmatrix(pchmatrix G, const char *filename)
{
  FILE     *out;

  out = fopen(filename, "w");
  assert(out != 0);

  fprintf(out, "Beginn der Matrix\n");
  write_hlib_part(G, 0, 0, out);
  fprintf(out, "Ende der Matrix\n");
}

static phmatrix
read_hlib_part(FILE * in, uint roff, uint coff, pcluster *rc,
	       pcluster *cc, uint * ridx, uint * cidx, uint * lineno)
{
  char      buf[80];
  char     *line, *c;
  uint      fields;
  phmatrix  G, G1;
  prkmatrix r;
  pamatrix  f;
  pcluster *rc1, *cc1;
  uint     *ridx1, *cidx1;
  uint      rows, cols;
  uint      rsons, csons;
  uint      roff1, coff1;
  uint      lda, ldb, ldf;
  uint      mtype;
  uint      i, j, k, l;

  G = NULL;
  G1 = NULL;

  line = fgets(buf, 80, in);
  (*lineno)++;

  if (line == 0) {
    fprintf(stderr, "Unexpected end of file\n");

    return 0;
  }

  /* Determine matrix type */
  fields = sscanf(line, "Type %u", &mtype);
  if (fields != 1) {
    for (c = line; *c != '\0' && *c != '\n'; c++);
    if (*c == '\n')
      *c = '\0';

    fprintf(stderr, "Expected \"Type ...\", got \"%s\"\n", line);

    return 0;
  }

  /* Read offsets and compare to what we expect */
  line = fgets(buf, 80, in);
  (*lineno)++;
  fields = sscanf(line, "row_offset = %u", &roff1);
  assert(fields == 1);

  line = fgets(buf, 80, in);
  (*lineno)++;
  fields = sscanf(line, "col_offset = %u", &coff1);
  assert(fields == 1);

  assert(roff == roff1);
  assert(coff == coff1);

  /* Read number of rows and columns */
  line = fgets(buf, 80, in);
  (*lineno)++;
  fields = sscanf(line, "rows = %u", &rows);
  assert(fields == 1);

  line = fgets(buf, 80, in);
  (*lineno)++;
  fields = sscanf(line, "cols = %u", &cols);
  assert(fields == 1);

  /* If we have not inherited a row index, create one */
  if (ridx == 0) {
    ridx = (uint *) allocmem(sizeof(uint) * rows);
    for (i = 0; i < rows; i++)
      ridx[i] = i;
  }

  /* If we have not inherited a column index, create one */
  if (cidx == 0) {
    cidx = (uint *) allocmem(sizeof(uint) * cols);
    for (j = 0; j < cols; j++)
      cidx[j] = j;
  }

  /* If we have not inherited a row cluster, create one */
  if (*rc == 0)
    *rc = new_cluster(rows, ridx, 0, 1);

  /* If we have not inherited a column cluster, create one */
  if (*cc == 0)
    *cc = new_cluster(cols, cidx, 0, 1);

  switch (mtype) {
  case 1:			/* rkmatrix */
    /* Read rank */
    line = fgets(buf, 80, in);
    (*lineno)++;
    fields = sscanf(line, "k = %u", &k);
    assert(fields == 1);

    line = fgets(buf, 80, in);
    (*lineno)++;
    fields = sscanf(line, "kt = %u", &k);
    assert(fields == 1);

    /* Create matrix */
    G = new_rk_hmatrix(*rc, *cc, k);
    r = G->r;
    lda = r->A.ld;
    ldb = r->B.ld;

    /* Read A */
    for (l = 0; l < k; l++)
      for (i = 0; i < rows; i++) {
	line = fgets(buf, 80, in);
	(*lineno)++;
	fields = sscanf(line, FIELD_SCANF_CS(, e),
			FIELD_ADDR(r->A.a + i + l * lda));
	assert(fields == 1);
      }

    /* Read B */
    for (l = 0; l < k; l++)
      for (j = 0; j < cols; j++) {
	line = fgets(buf, 80, in);
	(*lineno)++;
	fields = sscanf(line, FIELD_SCANF_CS(, e),
			FIELD_ADDR(r->B.a + j + l * ldb));
	assert(fields == 1);
      }
    break;

  case 2:			/* amatrix */
    /* Create matrix */
    G = new_full_hmatrix(*rc, *cc);
    f = G->f;
    ldf = f->ld;

    /* Read coefficients */
    for (j = 0; j < cols; j++)
      for (i = 0; i < rows; i++) {
	line = fgets(buf, 80, in);
	(*lineno)++;
	fields = sscanf(line, FIELD_SCANF_CS(, e),
			FIELD_ADDR(f->a + i + j * ldf));
	assert(fields == 1);
      }
    break;

  case 3:			/* subdivided matrix */
    /* Read number of block rows and columns */
    line = fgets(buf, 80, in);
    (*lineno)++;
    fields = sscanf(line, "block_rows = %u", &rsons);
    assert(fields == 1);

    line = fgets(buf, 80, in);
    (*lineno)++;
    fields = sscanf(line, "block_cols = %u", &csons);
    assert(fields == 1);

    /* Update row cluster if it didn't have sons before */
    if ((*rc)->sons == 0)
      setsons_cluster(*rc, rsons);
    else
      assert((*rc)->sons == rsons);

    /* Update column cluster if it didn't have sons before */
    if ((*cc)->sons == 0)
      setsons_cluster(*cc, csons);
    else
      assert((*cc)->sons == csons);

    G = new_super_hmatrix(*rc, *cc, rsons, csons);

    rc1 = (*rc)->son;
    cc1 = (*cc)->son;

    coff1 = coff;
    cidx1 = cidx;
    for (j = 0; j < csons; j++) {
      roff1 = roff;
      ridx1 = ridx;

      for (i = 0; i < rsons; i++) {
	G1 = read_hlib_part(in, roff1, coff1, rc1 + i, cc1 + j, ridx1, cidx1,
			    lineno);
	ref_hmatrix(G->son + i + j * rsons, G1);

	roff1 += G1->rc->size;
	ridx1 += G1->rc->size;
      }
      assert(roff1 == roff + rows);

      coff1 += G1->cc->size;
      cidx1 += G1->cc->size;
    }
    assert(coff1 == coff + cols);
    break;

  default:
    ;
  }

  update_hmatrix(G);

  return G;
}

phmatrix
read_hlib_hmatrix(const char *filename)
{
  phmatrix  G;
  FILE     *in;
  char      buf[80], *res;
  pcluster  rc, cc;
  uint      lineno;

  in = fopen(filename, "r");
  assert(in != 0);

  res = fgets(buf, 80, in);
  assert(res != NULL);
  lineno = 1;

  rc = 0;
  cc = 0;
  G = read_hlib_part(in, 0, 0, &rc, &cc, 0, 0, &lineno);

  res = fgets(buf, 80, in);
  assert(res != NULL);
  lineno++;

  fclose(in);

  return G;
}

static phmatrix
read_hlibsymm_part(FILE * in, uint roff, pcluster *rc,
		   uint * ridx, uint * lineno)
{
  char      buf[80];
  char     *line, *c;
  uint      fields;
  phmatrix  G, G1;
  prkmatrix r;
  pamatrix  f;
  pcluster *rc1;
  uint     *ridx1, *cidx1;
  uint      rows, cols;
  uint      rsons, csons;
  uint      roff1, coff1;
  uint      lda, ldb, ldf;
  uint      mtype;
  uint      i, j, k, l;

  line = fgets(buf, 80, in);
  (*lineno)++;

  if (line == 0) {
    fprintf(stderr, "Unexpected end of file\n");

    return 0;
  }

  /* Determine matrix type */
  fields = sscanf(line, "Type %u", &mtype);
  if (fields != 1) {
    for (c = line; *c != '\0' && *c != '\n'; c++);
    if (*c == '\n')
      *c = '\0';

    fprintf(stderr, "Expected \"Type ...\", got \"%s\"\n", line);

    return 0;
  }

  /* Read offsets and compare to what we expect */
  line = fgets(buf, 80, in);
  (*lineno)++;
  fields = sscanf(line, "row_offset = %u", &roff1);
  assert(fields == 1);

  line = fgets(buf, 80, in);
  (*lineno)++;
  fields = sscanf(line, "col_offset = %u", &coff1);
  assert(fields == 1);

  assert(roff == roff1);
  assert(roff == coff1);

  /* Read number of rows and columns */
  line = fgets(buf, 80, in);
  (*lineno)++;
  fields = sscanf(line, "rows = %u", &rows);
  assert(fields == 1);

  line = fgets(buf, 80, in);
  (*lineno)++;
  fields = sscanf(line, "cols = %u", &cols);
  assert(fields == 1);

  assert(rows == cols);

  /* If we have not inherited a row index, create one */
  if (ridx == 0) {
    ridx = (uint *) allocmem(sizeof(uint) * rows);
    for (i = 0; i < rows; i++)
      ridx[i] = i;
  }

  /* If we have not inherited a row cluster, create one */
  if (*rc == 0)
    *rc = new_cluster(rows, ridx, 0, 1);

  G = 0;
  switch (mtype) {
  case 1:			/* rkmatrix */
    /* Read rank */
    line = fgets(buf, 80, in);
    (*lineno)++;
    fields = sscanf(line, "k = %u", &k);
    assert(fields == 1);

    line = fgets(buf, 80, in);
    (*lineno)++;
    fields = sscanf(line, "kt = %u", &k);
    assert(fields == 1);

    /* Create matrix */
    G = new_rk_hmatrix(*rc, *rc, k);
    r = G->r;
    lda = r->A.ld;
    ldb = r->B.ld;

    /* Read A */
    for (l = 0; l < k; l++)
      for (i = 0; i < rows; i++) {
	line = fgets(buf, 80, in);
	(*lineno)++;
	fields = sscanf(line, FIELD_SCANF_CS(, e),
			FIELD_ADDR(r->A.a + i + l * lda));
	assert(fields == 1);
      }

    /* Read B */
    for (l = 0; l < k; l++)
      for (j = 0; j < cols; j++) {
	line = fgets(buf, 80, in);
	(*lineno)++;
	fields = sscanf(line, FIELD_SCANF_CS(, e),
			FIELD_ADDR(r->B.a + j + l * ldb));
	assert(fields == 1);
      }
    break;

  case 2:			/* amatrix */
    /* Create matrix */
    G = new_full_hmatrix(*rc, *rc);
    f = G->f;
    ldf = f->ld;

    /* Read coefficients */
    for (j = 0; j < cols; j++)
      for (i = 0; i < rows; i++) {
	line = fgets(buf, 80, in);
	(*lineno)++;
	fields = sscanf(line, FIELD_SCANF_CS(, e),
			FIELD_ADDR(f->a + i + j * ldf));
	assert(fields == 1);
      }
    break;

  case 3:			/* subdivided matrix */
    /* Read number of block rows and columns */
    line = fgets(buf, 80, in);
    (*lineno)++;
    fields = sscanf(line, "block_rows = %u", &rsons);
    assert(fields == 1);

    line = fgets(buf, 80, in);
    (*lineno)++;
    fields = sscanf(line, "block_cols = %u", &csons);
    assert(fields == 1);
    assert(rsons == csons);

    /* Update row cluster if it didn't have sons before */
    if ((*rc)->sons == 0)
      setsons_cluster(*rc, rsons);
    else
      assert((*rc)->sons == rsons);

    G = new_super_hmatrix(*rc, *rc, rsons, csons);

    rc1 = (*rc)->son;

    coff1 = roff;
    cidx1 = ridx;
    for (j = 0; j < csons; j++) {
      roff1 = roff;
      ridx1 = ridx;

      for (i = 0; i < rsons; i++) {
	G1 = (i == j ?
	      read_hlibsymm_part(in, roff1, rc1 + i, ridx1, lineno) :
	      read_hlib_part(in, roff1, coff1, rc1 + i, rc1 + j, ridx1, ridx1,
			     lineno));
	ref_hmatrix(G->son + i + j * rsons, G1);

	roff1 += G1->rc->size;
	ridx1 += G1->rc->size;
      }
      assert(roff1 == roff + rows);

      coff1 += G->son[j * rsons]->cc->size;
      cidx1 += G->son[j * rsons]->cc->size;
    }
    assert(coff1 == roff + cols);
    break;

  default:
    ;
  }

  return G;
}

phmatrix
read_hlibsymm_hmatrix(const char *filename)
{
  phmatrix  G;
  FILE     *in;
  char      buf[80], *line;
  pcluster  rc;
  uint      lineno;

  in = fopen(filename, "r");
  assert(in != 0);

  line = fgets(buf, 80, in);
  assert(line != 0);
  lineno = 1;

  rc = 0;
  G = read_hlibsymm_part(in, 0, &rc, 0, &lineno);

  line = fgets(buf, 80, in);
  assert(line != 0);
  lineno++;

  fclose(in);

  return G;
}

/* ------------------------------------------------------------
 * Drawing
 * ------------------------------------------------------------ */

#ifdef USE_CAIRO
static void
cairodraw(cairo_t * cr, pchmatrix hm, bool storage, uint levels)
{
  cairo_text_extents_t extents;
  char      buf[5];
  uint      rsons, csons;
  uint      rsize, csize;
  uint      roff, coff;
  uint      i, j;

  if (hm->son && levels != 1) {
    rsons = hm->rsons;
    csons = hm->csons;

    coff = 0;
    for (j = 0; j < csons; j++) {
      roff = 0;
      for (i = 0; i < rsons; i++) {
	cairo_save(cr);
	cairo_translate(cr, coff, roff);
	cairodraw(cr, hm->son[i + j * rsons], storage, levels - 1);
	cairo_restore(cr);

	roff += hm->son[i + j * rsons]->rc->size;
      }
      assert(roff == hm->rc->size);

      coff += hm->son[j * rsons]->cc->size;
    }
    assert(coff == hm->cc->size);
  }
  else {
    rsize = hm->rc->size;
    csize = hm->cc->size;

    if (hm->son) {
      cairo_rectangle(cr, 0.0, 0.0, csize, rsize);
      cairo_save(cr);
      cairo_set_source_rgb(cr, 0.9, 0.9, 1.0);
      cairo_fill_preserve(cr);
      cairo_restore(cr);
      cairo_stroke(cr);
    }
    else if (hm->r) {
      if (storage) {
	cairo_save(cr);
	cairo_set_source_rgb(cr, 0.2, 1.0, 0.2);
	cairo_rectangle(cr, 0.0, 0.0, csize, UINT_MIN(rsize, hm->r->k));
	cairo_rectangle(cr, 0.0, 0.0, UINT_MIN(csize, hm->r->k), rsize);
	cairo_fill(cr);
	cairo_restore(cr);

	cairo_rectangle(cr, 0.0, 0.0, csize, rsize);
	cairo_stroke(cr);

	if (hm->r->k > 0) {
	  sprintf(buf, "%d", hm->r->k);
	  cairo_set_font_size(cr, UINT_MIN(rsize, csize) / 2.5);
	  cairo_text_extents(cr, buf, &extents);

	  cairo_save(cr);
	  cairo_translate(cr, (csize - extents.width) / 2.0,
			  (rsize + extents.height) / 2.0);
	  cairo_show_text(cr, buf);
	  cairo_restore(cr);
	}

      }
      else {
	cairo_rectangle(cr, 0.0, 0.0, csize, rsize);
	if (hm->r->k > 0) {
	  cairo_save(cr);
	  cairo_set_source_rgb(cr, 0.2, 1.0, 0.2);
	  cairo_fill_preserve(cr);
	  cairo_restore(cr);
	}
	if (hm->r->k == 0) {
	  cairo_save(cr);
	  cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);
	  cairo_fill_preserve(cr);
	  cairo_restore(cr);
	}

	cairo_stroke(cr);

	if (hm->r->k > 0) {
	  sprintf(buf, "%d", hm->r->k);
	  cairo_set_font_size(cr, UINT_MIN(rsize, csize) / 2.5);
	  cairo_text_extents(cr, buf, &extents);

	  cairo_save(cr);
	  cairo_translate(cr, (csize - extents.width) / 2.0,
			  (rsize + extents.height) / 2.0);
	  cairo_show_text(cr, buf);
	  cairo_restore(cr);
	}
      }
    }
    else if (hm->f) {
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
draw_cairo_hmatrix(cairo_t * cr, pchmatrix hm, bool storage, uint levels)
{
  double    sx, sy, ex, ey;
  uint      rsize, csize;
  double    scalex, scaley, scale;

  /* Save Cairo state */
  cairo_save(cr);

  /* Obtain size of block */
  rsize = hm->rc->size;
  csize = hm->cc->size;

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
  cairodraw(cr, hm, storage, levels);

  /* Restore Cairo state */
  cairo_restore(cr);
}
#endif

void
copy_sparsematrix_hmatrix(psparsematrix sp, phmatrix hm)
{

  uint      i, j, k, l, m, rsize, csize;
  uint     *ridx, *cidx, *row, *col;
  pfield    coeff;
  pamatrix  f;

  if (hm->son) {
    for (j = 0; j < hm->csons; j++) {
      for (i = 0; i < hm->rsons; i++) {
	copy_sparsematrix_hmatrix(sp, hm->son[i + j * (hm->rsons)]);
      }
    }
  }
  else if (hm->r) {		/*hm->r is defined, but hm->r->A = hm->r->B = NulL => all entries are 0 */

  }
  else if (hm->f) {
    rsize = hm->rc->size;
    csize = hm->cc->size;
    ridx = hm->rc->idx;
    cidx = hm->cc->idx;
    row = sp->row;
    col = sp->col;
    coeff = sp->coeff;
    f = hm->f;

    for (j = 0; j < f->cols; j++)
      for (i = 0; i < f->rows; i++)
	f->a[i + j * f->rows] = 0.0;

    for (k = 0; k < rsize; k++) {
      i = ridx[k];
      for (m = row[i]; m < row[i + 1]; m++) {
	for (l = 0; l < csize; l++) {

	  if (col[m] == cidx[l]) {

	    setentry_amatrix(hm->f, k, l, coeff[m]);
	  }
	}
      }
    }

  }
}
