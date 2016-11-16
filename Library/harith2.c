
/* ------------------------------------------------------------
 * This is the file "harith2.c" of the H2Lib package.
 * All rights reserved, Steffen Boerm 2015
 * ------------------------------------------------------------ */

#include "harith2.h"

/* ------------------------------------------------------------
 * Representation of structured matrix products
 * ------------------------------------------------------------ */

struct _hprodentry {
  field     alpha;
  bool      xtrans;
  pchmatrix x;
  bool      ytrans;
  pchmatrix y;

  phprodentry next;
};

static    phprodentry
new_hprodentry(field alpha, bool xtrans, pchmatrix x,
	       bool ytrans, pchmatrix y, phprodentry next)
{
  phprodentry he;

  he = (phprodentry) allocmem(sizeof(hprodentry));
  he->alpha = alpha;
  he->xtrans = xtrans;
  he->x = x;
  he->ytrans = ytrans;
  he->y = y;
  he->next = next;

  return he;
}

static    phprodentry
del_hprodentry(phprodentry he)
{
  phprodentry next;

  next = he->next;
  freemem(he);

  return next;
}

/* ------------------------------------------------------------
 * Constructor and destructor
 * ------------------------------------------------------------ */

phaccum
new_haccum(phmatrix z, pctruncmode tm, real eps)
{
  phaccum   ha;

  ha = (phaccum) allocmem(sizeof(haccum));
  ha->z = z;
  ha->tm = tm;
  ha->eps = eps;
  ha->r = new_rkmatrix(z->rc->size, z->cc->size, 0);
  ha->xy = 0;

  return ha;
}

void
del_haccum(phaccum ha)
{
  phprodentry xy;

  xy = ha->xy;
  while (xy)
    xy = del_hprodentry(xy);
  del_rkmatrix(ha->r);
  freemem(ha);
}

/* ------------------------------------------------------------
 * Add a product to an H-matrix accumulator
 * ------------------------------------------------------------ */

void
addproduct_haccum(field alpha, bool xtrans, pchmatrix x, bool ytrans,
		  pchmatrix y, phaccum za)
{
  rkmatrix  tmp1;
  amatrix   tmp2;
  phmatrix  z = za->z;
  prkmatrix r = za->r;
  pccluster rc = z->rc;
  pccluster cc = z->cc;
  pctruncmode tm = za->tm;
  real      eps = za->eps;
  pamatrix  xf, yf, zf, id;
  prkmatrix xy;
  uint      k;

  assert(xtrans ? (x->cc == z->rc) : (x->rc == z->rc));
  assert(ytrans ? (y->rc == z->cc) : (y->cc == z->cc));

  if (x->f) {
    xf = x->f;
    if (z->f)			/* Z = Z + X Y  <=>  Z^* = Z^* + Y^* X^* */
      addmul_hmatrix_amatrix_amatrix(alpha, !ytrans, y, !xtrans, xf, true,
				     z->f);
    else {
      if (xtrans ? (xf->cols > xf->rows) : (xf->rows > xf->cols)) {
	/* Compute rkmatrix X Y = X (Y^* I^*)^* */
	if (xtrans) {
	  k = xf->rows;
	  xy = init_rkmatrix(&tmp1, rc->size, cc->size, k);
	  copy_amatrix(true, xf, &xy->A);
	}
	else {
	  k = xf->cols;
	  xy = init_rkmatrix(&tmp1, rc->size, cc->size, k);
	  copy_amatrix(false, xf, &xy->A);
	}
	id = init_amatrix(&tmp2, k, k);
	identity_amatrix(id);
	clear_amatrix(&xy->B);
	addmul_hmatrix_amatrix_amatrix(alpha, !ytrans, y, true, id, false,
				       &xy->B);

	/* Add rkmatrix to accumulator */
	assert(z->f == 0);
	add_rkmatrix(1.0, xy, tm, eps, r);

	/* Clean up */
	uninit_amatrix(id);
	uninit_rkmatrix(xy);
      }
      else {
	/* Compute matrix Zf = X Y, Zf^* = Y^* X^* */
	zf = init_amatrix(&tmp2, rc->size, cc->size);
	clear_amatrix(zf);
	addmul_hmatrix_amatrix_amatrix(alpha, !ytrans, y, !xtrans, xf, true,
				       zf);

	/* Add Zf to Z */
	add_amatrix_destructive_hmatrix(1.0, false, zf, tm, eps, z);

	/* Clean up */
	uninit_amatrix(zf);
      }
    }
  }
  else if (x->r) {
    /* Compute rkmatrix X Y = A (Y^* B)^* */
    k = x->r->k;
    xy = init_rkmatrix(&tmp1, rc->size, cc->size, k);
    if (xtrans) {
      copy_amatrix(false, &x->r->B, &xy->A);
      clear_amatrix(&xy->B);
      addmul_hmatrix_amatrix_amatrix(alpha, !ytrans, y, false, &x->r->A,
				     false, &xy->B);
    }
    else {
      copy_amatrix(false, &x->r->A, &xy->A);
      clear_amatrix(&xy->B);
      addmul_hmatrix_amatrix_amatrix(alpha, !ytrans, y, false, &x->r->B,
				     false, &xy->B);
    }

    /* Add rkmatrix to Z or accumulator */
    if (z->f)
      add_rkmatrix_amatrix(1.0, false, xy, z->f);
    else
      add_rkmatrix(1.0, xy, tm, eps, r);

    /* Clean up */
    uninit_rkmatrix(xy);
  }
  else {
    assert(x->son);

    if (y->f) {
      yf = y->f;
      if (z->f)			/* Z = Z + X Y */
	addmul_hmatrix_amatrix_amatrix(alpha, xtrans, x, ytrans, yf, false,
				       z->f);
      else {
	if (ytrans ? (yf->rows > yf->cols) : (yf->cols > yf->rows)) {
	  /* Compute rkmatrix X Y = (X I) (Y^*)^* */
	  if (ytrans) {
	    k = yf->cols;
	    xy = init_rkmatrix(&tmp1, rc->size, cc->size, k);
	    copy_amatrix(false, yf, &xy->B);
	  }
	  else {
	    k = yf->rows;
	    xy = init_rkmatrix(&tmp1, rc->size, cc->size, k);
	    copy_amatrix(true, yf, &xy->B);
	  }
	  id = init_amatrix(&tmp2, k, k);
	  identity_amatrix(id);
	  clear_amatrix(&xy->A);
	  addmul_hmatrix_amatrix_amatrix(alpha, xtrans, x, false, id, false,
					 &xy->A);

	  /* Add rkmatrix to accumulator */
	  assert(z->f == 0);
	  add_rkmatrix(1.0, xy, tm, eps, r);

	  /* Clean up */
	  uninit_amatrix(id);
	  uninit_rkmatrix(xy);
	}
	else {
	  /* Compute matrix Zf = X Y */
	  zf = init_amatrix(&tmp2, rc->size, cc->size);
	  clear_amatrix(zf);
	  addmul_hmatrix_amatrix_amatrix(alpha, xtrans, x, ytrans, yf, false,
					 zf);

	  /* Add Zf to Z */
	  add_amatrix_destructive_hmatrix(1.0, false, zf, tm, eps, z);

	  /* Clean up */
	  uninit_amatrix(zf);
	}
      }
    }
    else if (y->r) {
      /* Compute rkmatrix X Y = (X A) B^* */
      k = y->r->k;
      xy = init_rkmatrix(&tmp1, rc->size, cc->size, k);
      if (ytrans) {
	copy_amatrix(false, &y->r->A, &xy->B);
	clear_amatrix(&xy->A);
	addmul_hmatrix_amatrix_amatrix(alpha, xtrans, x, false, &y->r->B,
				       false, &xy->A);
      }
      else {
	copy_amatrix(false, &y->r->B, &xy->B);
	clear_amatrix(&xy->A);
	addmul_hmatrix_amatrix_amatrix(alpha, xtrans, x, false, &y->r->A,
				       false, &xy->A);
      }

      /* Add rkmatrix to Z or accumulator */
      if (z->f)
	add_rkmatrix_amatrix(1.0, false, xy, z->f);
      else
	add_rkmatrix(1.0, xy, tm, eps, r);

      /* Clean up */
      uninit_rkmatrix(xy);
    }
    else {
      assert(y->son);

      if ((xtrans ? x->cc == x->son[0]->cc : x->rc == x->son[0]->rc)
	  && (ytrans ? y->rc == y->son[0]->rc : y->cc == y->son[0]->cc)) {
	/* X and Y have submatrices, but they all contribute
	 * directly to Z */
	if (xtrans) {
	  assert(x->csons == 1);
	  if (ytrans) {
	    assert(y->rsons == 1);
	    assert(x->rsons == y->csons);

	    for (k = 0; k < x->rsons; k++)
	      addproduct_haccum(alpha, xtrans, x->son[k], ytrans, y->son[k],
				za);
	  }
	  else {
	    assert(y->csons == 1);
	    assert(x->rsons == y->rsons);

	    for (k = 0; k < x->rsons; k++)
	      addproduct_haccum(alpha, xtrans, x->son[k], ytrans, y->son[k],
				za);
	  }
	}
	else {
	  assert(x->rsons == 1);
	  if (ytrans) {
	    assert(y->rsons == 1);
	    assert(x->csons == y->csons);

	    for (k = 0; k < x->csons; k++)
	      addproduct_haccum(alpha, xtrans, x->son[k], ytrans, y->son[k],
				za);
	  }
	  else {
	    assert(y->csons == 1);
	    assert(x->csons == y->rsons);

	    for (k = 0; k < x->csons; k++)
	      addproduct_haccum(alpha, xtrans, x->son[k], ytrans, y->son[k],
				za);
	  }
	}
      }
      else
	za->xy = new_hprodentry(alpha, xtrans, x, ytrans, y, za->xy);
    }
  }
}

/* ------------------------------------------------------------
 * Split an H-matrix accumulator
 * ------------------------------------------------------------ */

phaccum  *
split_haccum(phmatrix z, phaccum za)
{
  phaccum  *zason;
  phmatrix  rson;
  phprodentry he;
  uint      rsons, csons;
  field     alpha;
  bool      xtrans, ytrans;
  pchmatrix x, y;
  uint      i, j, k;

  assert(z->son);

  rsons = z->rsons;
  csons = z->csons;

  zason = (phaccum *) allocmem(sizeof(phaccum) * rsons * csons);
  rson = split_sub_rkmatrix(za->r, z->rc, z->cc, (z->son[0]->rc != z->rc),
			    (z->son[0]->cc != z->cc));

  for (j = 0; j < csons; j++)
    for (i = 0; i < rsons; i++) {
      zason[i + j * rsons] =
	new_haccum(z->son[i + j * rsons], za->tm, za->eps);
      copy_rkmatrix(false, rson->son[i + j * rsons]->r,
		    zason[i + j * rsons]->r);

      for (he = za->xy; he; he = he->next) {
	alpha = he->alpha;
	xtrans = he->xtrans;
	x = he->x;
	ytrans = he->ytrans;
	y = he->y;

	if (xtrans) {
	  assert(rsons == x->csons);
	  assert(z->rc == x->cc);

	  if (ytrans) {
	    assert(csons == y->rsons);
	    assert(z->cc == y->rc);

	    assert(x->rc == y->cc);
	    assert(x->rsons == y->csons);

	    for (k = 0; k < x->rsons; k++)
	      addproduct_haccum(alpha, xtrans, x->son[k + i * x->rsons],
				ytrans, y->son[j + k * y->rsons],
				zason[i + j * rsons]);
	  }
	  else {
	    assert(csons == y->csons);
	    assert(z->cc == y->cc);

	    assert(x->rc == y->rc);
	    assert(x->rsons == y->rsons);

	    for (k = 0; k < x->rsons; k++)
	      addproduct_haccum(alpha, xtrans, x->son[k + i * x->rsons],
				ytrans, y->son[k + j * y->rsons],
				zason[i + j * rsons]);
	  }
	}
	else {
	  assert(rsons == x->rsons);
	  assert(z->rc == x->rc);

	  if (ytrans) {
	    assert(csons == y->rsons);
	    assert(z->cc == y->rc);

	    assert(x->cc == y->cc);
	    assert(x->csons == y->csons);

	    for (k = 0; k < x->csons; k++)
	      addproduct_haccum(alpha, xtrans, x->son[i + k * x->rsons],
				ytrans, y->son[j + k * y->rsons],
				zason[i + j * rsons]);
	  }
	  else {
	    assert(csons == y->csons);
	    assert(z->cc == y->cc);

	    assert(x->cc == y->rc);
	    assert(x->csons == y->rsons);

	    for (k = 0; k < x->csons; k++)
	      addproduct_haccum(alpha, xtrans, x->son[i + k * x->rsons],
				ytrans, y->son[k + j * y->rsons],
				zason[i + j * rsons]);
	  }
	}
      }
    }

  del_hmatrix(rson);

  return zason;
}

/* ------------------------------------------------------------
 * Flush an H-matrix accumulator
 * ------------------------------------------------------------ */

void
flush_haccum(phaccum ha)
{
  phaccum  *hason;
  phmatrix  ztmp;
  prkmatrix rtmp;
  uint      rsons, csons;
  uint      i, j;

  if (ha->xy == 0) {
    /* No structured products remain, only the low-rank update
     * has to be handled. */

    if (ha->r->k > 0) {
      add_rkmatrix_hmatrix(1.0, ha->r, ha->tm, ha->eps, ha->z);

      setrank_rkmatrix(ha->r, 0);
    }
  }
  else {
    rsons = (ha->z->rc->sons > 0 ? ha->z->rc->sons : 1);
    csons = (ha->z->cc->sons > 0 ? ha->z->cc->sons : 1);

    if (ha->z->son) {
      /* Z has submatrices: split the accumulator to match the
       * submatrices and proceed by recursion. */

      hason = split_haccum(ha->z, ha);

      for (j = 0; j < csons; j++)
	for (i = 0; i < rsons; i++) {
	  flush_haccum(hason[i + j * rsons]);

	  del_haccum(hason[i + j * rsons]);
	}

      freemem(hason);
    }
    else if (ha->z->f) {
      /* Z is a dense matrix: split Z into submatrices, recursively
       * add to these submatrices. */

      ztmp = split_sub_amatrix(ha->z->f, ha->z->rc, ha->z->cc, true, true);

      hason = split_haccum(ztmp, ha);

      for (j = 0; j < csons; j++)
	for (i = 0; i < rsons; i++) {
	  flush_haccum(hason[i + j * rsons]);

	  del_haccum(hason[i + j * rsons]);
	}

      freemem(hason);

      del_hmatrix(ztmp);
    }
    else {
      assert(ha->z->r);
      /* Z is a low-rank matrix: split Z into submatrices, recursively
       * add to these submatrices, then merge the submatrices. */

      ztmp = split_rkmatrix(ha->z->r, ha->z->rc, ha->z->cc, true, true, true);

      hason = split_haccum(ztmp, ha);

      for (j = 0; j < csons; j++)
	for (i = 0; i < rsons; i++) {
	  flush_haccum(hason[i + j * rsons]);

	  del_haccum(hason[i + j * rsons]);
	}

      freemem(hason);

      rtmp = merge_hmatrix_rkmatrix(ztmp, ha->tm, ha->eps);

      copy_rkmatrix(false, rtmp, ha->z->r);

      del_rkmatrix(rtmp);
      del_hmatrix(ztmp);
    }
  }
}

void
add_amatrix_destructive_rkmatrix(field alpha, bool atrans, pamatrix a,
				 pctruncmode tm, real eps, prkmatrix b)
{
  amatrix   tmp1, tmp2;
  realavector tmp3;
  pamatrix  u, vt;
  prealavector sigma;
  uint      rows, cols, k, knew;

#ifdef HARITH_AMATRIX_QUICK_EXIT
  if (normfrob_amatrix(a) == 0.0)
    return;
#endif

  if (atrans) {
    assert(a->rows == b->B.rows);
    assert(a->cols == b->A.rows);
  }
  else {
    assert(a->rows == b->A.rows);
    assert(a->cols == b->B.rows);
  }

  rows = b->A.rows;
  cols = b->B.rows;

  addmul_amatrix(alpha, false, &b->A, true, &b->B, a);

  /* Find singular value decomposition of A */
  k = UINT_MIN(rows, cols);
  u = init_amatrix(&tmp1, rows, k);
  vt = init_amatrix(&tmp2, k, cols);
  sigma = init_realavector(&tmp3, k);
  svd_amatrix(a, sigma, u, vt);

  /* Determine rank */
  knew = findrank_truncmode(tm, eps, sigma);

  /* Set new rank */
  setrank_rkmatrix(b, knew);

  /* Scale singular vectors */
  diageval_realavector_amatrix(1.0, true, sigma, true, u);

  /* Copy singular vectors */
  clear_amatrix(&b->A);
  copy_sub_amatrix(false, u, &b->A);
  clear_amatrix(&b->B);
  copy_sub_amatrix(true, vt, &b->B);

  /* Clean up */
  uninit_realavector(sigma);
  uninit_amatrix(vt);
  uninit_amatrix(u);
}

void
add_amatrix_destructive_hmatrix(field alpha, bool atrans, pamatrix a,
				pctruncmode tm, real eps, phmatrix b)
{
  amatrix   tmp;
  phmatrix  b1;
  pamatrix  a1;
  uint      rsons, csons;
  uint      roff, coff;
  uint      i, j;

#ifdef HARITH_AMATRIX_QUICK_EXIT
  if (normfrob_amatrix(a) == 0.0)
    return;
#endif

  if (b->r)
    add_amatrix_destructive_rkmatrix(alpha, atrans, a, tm, eps, b->r);
  else if (b->f)
    add_amatrix(alpha, atrans, a, b->f);
  else {
    rsons = b->rsons;
    csons = b->csons;

    coff = 0;
    for (j = 0; j < csons; j++) {
      roff = 0;
      for (i = 0; i < rsons; i++) {
	b1 = b->son[i + j * rsons];
	a1 =
	  init_sub_amatrix(&tmp, a, b1->rc->size, roff, b1->cc->size, coff);

	add_amatrix_destructive_hmatrix(alpha, atrans, a1, tm, eps, b1);

	uninit_amatrix(a1);

	roff += b1->rc->size;
      }
      assert(roff == b->rc->size);

      coff += b->son[j * rsons]->cc->size;
    }
    assert(coff == b->cc->size);
  }
}

/* ------------------------------------------------------------
 * Approximate the product of two H-matrices
 * ------------------------------------------------------------ */

void
addmul2_hmatrix(field alpha, bool xtrans, pchmatrix x, bool ytrans,
		pchmatrix y, pctruncmode tm, real eps, phmatrix z)
{
  phaccum   za;

  za = new_haccum(z, tm, eps);

  addproduct_haccum(alpha, xtrans, x, ytrans, y, za);

  flush_haccum(za);

  del_haccum(za);
}

/* ------------------------------------------------------------
 * H-matrix lower triangular solve
 * ------------------------------------------------------------ */

static void
lowersolve_nn_haccum(bool aunit, pchmatrix a, phaccum xa)
{
  phmatrix  atmp;
  phmatrix  x = xa->z;
  phaccum  *xa1;
  uint      sons;
  uint      i, j, k;

  assert(a->rc == a->cc);
  assert(a->cc == x->rc);

  if (x->f) {
    flush_haccum(xa);
    triangularinvmul_hmatrix_amatrix(true, aunit, false, a, false, x->f);
  }
  else if (x->r) {
    flush_haccum(xa);
    triangularinvmul_hmatrix_amatrix(true, aunit, false, a, false, &x->r->A);
  }
  else {
    if (a->f) {
      atmp = split_sub_amatrix(a->f, a->rc, a->cc, (x->son[0]->rc != x->rc),
			       (x->son[0]->rc != x->rc));

      assert(x->rsons == atmp->rsons);

      sons = atmp->rsons;

      xa1 = split_haccum(x, xa);

      for (k = 0; k < x->csons; k++)
	for (i = 0; i < sons; i++) {
	  lowersolve_nn_haccum(aunit, atmp->son[i + i * sons],
			       xa1[i + k * x->rsons]);

	  for (j = i + 1; j < sons; j++)
	    addproduct_haccum(-1.0, false, atmp->son[j + i * sons], false,
			      x->son[i + k * x->rsons],
			      xa1[j + k * x->rsons]);
	}

      for (k = 0; k < x->csons; k++)
	for (i = 0; i < sons; i++)
	  del_haccum(xa1[i + k * x->rsons]);
      freemem(xa1);

      del_hmatrix(atmp);
    }
    else {
      assert(a->rsons == a->csons);
      assert(x->rsons == a->rsons);

      sons = a->rsons;

      xa1 = split_haccum(x, xa);

      for (k = 0; k < x->csons; k++)
	for (i = 0; i < sons; i++) {
	  lowersolve_nn_haccum(aunit, a->son[i + i * sons],
			       xa1[i + k * x->rsons]);

	  for (j = i + 1; j < sons; j++)
	    addproduct_haccum(-1.0, false, a->son[j + i * sons], false,
			      x->son[i + k * x->rsons],
			      xa1[j + k * x->rsons]);
	}

      for (k = 0; k < x->csons; k++)
	for (i = 0; i < sons; i++)
	  del_haccum(xa1[i + k * x->rsons]);
      freemem(xa1);
    }
  }
}

static void
lowersolve_nt_haccum(bool aunit, pchmatrix a, phaccum xa)
{
  phmatrix  atmp;
  phmatrix  x = xa->z;
  phaccum  *xa1;
  uint      sons;
  uint      i, j, k;

  assert(a->rc == a->cc);
  assert(a->cc == x->cc);

  if (x->f) {
    flush_haccum(xa);
    triangularinvmul_hmatrix_amatrix(true, aunit, false, a, true, x->f);
  }
  else if (x->r) {
    flush_haccum(xa);
    triangularinvmul_hmatrix_amatrix(true, aunit, false, a, false, &x->r->B);
  }
  else {
    if (a->f) {
      atmp = split_sub_amatrix(a->f, a->rc, a->cc, (x->son[0]->cc != x->cc),
			       (x->son[0]->cc != x->cc));

      assert(x->csons == atmp->rsons);

      sons = atmp->rsons;

      xa1 = split_haccum(x, xa);

      for (k = 0; k < x->rsons; k++)
	for (i = 0; i < sons; i++) {
	  lowersolve_nt_haccum(aunit, atmp->son[i + i * sons],
			       xa1[k + i * x->rsons]);

	  for (j = i + 1; j < sons; j++)
	    addproduct_haccum(-1.0, false, x->son[k + i * x->rsons], true,
			      atmp->son[j + i * sons], xa1[k + j * x->rsons]);
	}

      for (k = 0; k < x->rsons; k++)
	for (i = 0; i < sons; i++)
	  del_haccum(xa1[k + i * x->rsons]);
      freemem(xa1);

      del_hmatrix(atmp);
    }
    else {
      assert(a->rsons == a->csons);
      assert(x->csons == a->rsons);

      sons = a->rsons;

      xa1 = split_haccum(x, xa);

      for (k = 0; k < x->rsons; k++)
	for (i = 0; i < sons; i++) {
	  lowersolve_nt_haccum(aunit, a->son[i + i * sons],
			       xa1[k + i * x->rsons]);

	  for (j = i + 1; j < sons; j++)
	    addproduct_haccum(-1.0, false, x->son[k + i * x->rsons], true,
			      a->son[j + i * sons], xa1[k + j * x->rsons]);
	}

      for (k = 0; k < x->rsons; k++)
	for (i = 0; i < sons; i++)
	  del_haccum(xa1[k + i * x->rsons]);
      freemem(xa1);
    }
  }
}

static void
lowersolve_tn_haccum(bool aunit, pchmatrix a, phaccum xa)
{
  phmatrix  atmp;
  phmatrix  x = xa->z;
  phaccum  *xa1;
  uint      sons;
  uint      i, j, k;

  assert(a->rc == a->cc);
  assert(a->cc == x->rc);

  if (x->f) {
    flush_haccum(xa);
    triangularinvmul_hmatrix_amatrix(true, aunit, true, a, false, x->f);
  }
  else if (x->r) {
    flush_haccum(xa);
    triangularinvmul_hmatrix_amatrix(true, aunit, true, a, false, &x->r->A);
  }
  else {
    if (a->f) {
      atmp = split_sub_amatrix(a->f, a->rc, a->cc, (x->son[0]->rc != x->rc),
			       (x->son[0]->rc != x->rc));

      assert(x->rsons == atmp->rsons);

      sons = atmp->rsons;

      xa1 = split_haccum(x, xa);

      for (k = 0; k < x->csons; k++)
	for (i = sons; i-- > 0;) {
	  lowersolve_tn_haccum(aunit, atmp->son[i + i * sons],
			       xa1[i + k * x->rsons]);

	  for (j = i; j-- > 0;)
	    addproduct_haccum(-1.0, true, atmp->son[i + j * sons], false,
			      x->son[i + k * x->rsons],
			      xa1[j + k * x->rsons]);
	}

      for (k = 0; k < x->csons; k++)
	for (i = 0; i < sons; i++)
	  del_haccum(xa1[i + k * x->rsons]);
      freemem(xa1);

      del_hmatrix(atmp);
    }
    else {
      assert(a->rsons == a->csons);
      assert(x->rsons == a->rsons);

      sons = a->rsons;

      xa1 = split_haccum(x, xa);

      for (k = 0; k < x->csons; k++)
	for (i = sons; i-- > 0;) {
	  lowersolve_tn_haccum(aunit, a->son[i + i * sons],
			       xa1[i + k * x->rsons]);

	  for (j = i; j-- > 0;)
	    addproduct_haccum(-1.0, true, a->son[i + j * sons], false,
			      x->son[i + k * x->rsons],
			      xa1[j + k * x->rsons]);
	}

      for (k = 0; k < x->csons; k++)
	for (i = 0; i < sons; i++)
	  del_haccum(xa1[i + k * x->rsons]);
      freemem(xa1);
    }
  }
}

static void
lowersolve_tt_haccum(bool aunit, pchmatrix a, phaccum xa)
{
  phmatrix  atmp;
  phmatrix  x = xa->z;
  phaccum  *xa1;
  uint      sons;
  uint      i, j, k;

  assert(a->rc == a->cc);
  assert(a->cc == x->cc);

  if (x->f) {
    flush_haccum(xa);
    triangularinvmul_hmatrix_amatrix(true, aunit, true, a, true, x->f);
  }
  else if (x->r) {
    flush_haccum(xa);
    triangularinvmul_hmatrix_amatrix(true, aunit, true, a, false, &x->r->B);
  }
  else {
    if (a->f) {
      atmp = split_sub_amatrix(a->f, a->rc, a->cc, (x->son[0]->cc != x->cc),
			       (x->son[0]->cc != x->cc));

      assert(x->csons == atmp->rsons);

      sons = atmp->rsons;

      xa1 = split_haccum(x, xa);

      for (k = 0; k < x->rsons; k++)
	for (i = sons; i-- > 0;) {
	  lowersolve_tt_haccum(aunit, atmp->son[i + i * sons],
			       xa1[k + i * x->rsons]);

	  for (j = i; j-- > 0;)
	    addproduct_haccum(-1.0, false, x->son[k + i * x->rsons], false,
			      atmp->son[i + j * sons], xa1[k + j * x->rsons]);
	}

      for (k = 0; k < x->rsons; k++)
	for (i = 0; i < sons; i++)
	  del_haccum(xa1[k + i * x->rsons]);
      freemem(xa1);

      del_hmatrix(atmp);
    }
    else {
      assert(a->rsons == a->csons);
      assert(x->csons == a->rsons);

      sons = a->rsons;

      xa1 = split_haccum(x, xa);

      for (k = 0; k < x->rsons; k++)
	for (i = sons; i-- > 0;) {
	  lowersolve_tt_haccum(aunit, a->son[i + i * sons],
			       xa1[k + i * x->rsons]);

	  for (j = i; j-- > 0;)
	    addproduct_haccum(-1.0, false, x->son[k + i * x->rsons], false,
			      a->son[i + j * sons], xa1[k + j * x->rsons]);
	}

      for (k = 0; k < x->rsons; k++)
	for (i = 0; i < sons; i++)
	  del_haccum(xa1[k + i * x->rsons]);
      freemem(xa1);
    }
  }
}

void
lowersolve_haccum(bool aunit, bool atrans, pchmatrix a, bool xtrans,
		  phaccum xa)
{
  if (atrans) {
    if (xtrans)
      lowersolve_tt_haccum(aunit, a, xa);
    else
      lowersolve_tn_haccum(aunit, a, xa);
  }
  else {
    if (xtrans)
      lowersolve_nt_haccum(aunit, a, xa);
    else
      lowersolve_nn_haccum(aunit, a, xa);
  }
}

void
lowersolve2_hmatrix_hmatrix(bool aunit, bool atrans, pchmatrix a,
			    pctruncmode tm, real eps, bool xtrans, phmatrix x)
{
  phaccum   xa;

  xa = new_haccum(x, tm, eps);

  lowersolve_haccum(aunit, atrans, a, xtrans, xa);

  del_haccum(xa);
}

/* ------------------------------------------------------------
 * H-matrix upper triangular solve
 * ------------------------------------------------------------ */

static void
uppersolve_nn_haccum(bool aunit, pchmatrix a, phaccum xa)
{
  phmatrix  atmp;
  phmatrix  x = xa->z;
  phaccum  *xa1;
  uint      sons;
  uint      i, j, k;

  assert(a->rc == a->cc);
  assert(a->cc == x->rc);

  if (x->f) {
    flush_haccum(xa);
    triangularinvmul_hmatrix_amatrix(false, aunit, false, a, false, x->f);
  }
  else if (x->r) {
    flush_haccum(xa);
    triangularinvmul_hmatrix_amatrix(false, aunit, false, a, false, &x->r->A);
  }
  else {
    if (a->f) {
      atmp = split_sub_amatrix(a->f, a->rc, a->cc, (x->son[0]->rc != x->rc),
			       (x->son[0]->rc != x->rc));

      assert(x->rsons == atmp->rsons);

      sons = atmp->rsons;

      xa1 = split_haccum(x, xa);

      for (k = 0; k < x->csons; k++)
	for (i = sons; i-- > 0;) {
	  uppersolve_nn_haccum(aunit, atmp->son[i + i * sons],
			       xa1[i + k * x->rsons]);

	  for (j = i; j-- > 0;)
	    addproduct_haccum(-1.0, false, atmp->son[j + i * sons], false,
			      x->son[i + k * x->rsons],
			      xa1[j + k * x->rsons]);
	}

      for (k = 0; k < x->csons; k++)
	for (i = 0; i < sons; i++)
	  del_haccum(xa1[i + k * x->rsons]);
      freemem(xa1);

      del_hmatrix(atmp);
    }
    else {
      assert(a->rsons == a->csons);
      assert(x->rsons == a->rsons);

      sons = a->rsons;

      xa1 = split_haccum(x, xa);

      for (k = 0; k < x->csons; k++)
	for (i = sons; i-- > 0;) {
	  uppersolve_nn_haccum(aunit, a->son[i + i * sons],
			       xa1[i + k * x->rsons]);

	  for (j = i; j-- > 0;)
	    addproduct_haccum(-1.0, false, a->son[j + i * sons], false,
			      x->son[i + k * x->rsons],
			      xa1[j + k * x->rsons]);
	}

      for (k = 0; k < x->csons; k++)
	for (i = 0; i < sons; i++)
	  del_haccum(xa1[i + k * x->rsons]);
      freemem(xa1);
    }
  }
}

static void
uppersolve_nt_haccum(bool aunit, pchmatrix a, phaccum xa)
{
  phmatrix  atmp;
  phmatrix  x = xa->z;
  phaccum  *xa1;
  uint      sons;
  uint      i, j, k;

  assert(a->rc == a->cc);
  assert(a->cc == x->cc);

  if (x->f) {
    flush_haccum(xa);
    triangularinvmul_hmatrix_amatrix(false, aunit, false, a, true, x->f);
  }
  else if (x->r) {
    flush_haccum(xa);
    triangularinvmul_hmatrix_amatrix(false, aunit, false, a, false, &x->r->B);
  }
  else {
    if (a->f) {
      atmp = split_sub_amatrix(a->f, a->rc, a->cc, (x->son[0]->cc != x->cc),
			       (x->son[0]->cc != x->cc));

      assert(x->csons == atmp->rsons);

      sons = atmp->rsons;

      xa1 = split_haccum(x, xa);

      for (k = 0; k < x->rsons; k++)
	for (i = sons; i-- > 0;) {
	  uppersolve_nt_haccum(aunit, atmp->son[i + i * sons],
			       xa1[k + i * x->rsons]);

	  for (j = i; j-- > 0;)
	    addproduct_haccum(-1.0, false, x->son[k + i * x->rsons], true,
			      atmp->son[j + i * sons], xa1[k + j * x->rsons]);
	}

      for (k = 0; k < x->rsons; k++)
	for (i = 0; i < sons; i++)
	  del_haccum(xa1[k + i * x->rsons]);
      freemem(xa1);

      del_hmatrix(atmp);
    }
    else {
      assert(a->rsons == a->csons);
      assert(x->csons == a->rsons);

      sons = a->rsons;

      xa1 = split_haccum(x, xa);

      for (k = 0; k < x->rsons; k++)
	for (i = sons; i-- > 0;) {
	  uppersolve_nt_haccum(aunit, a->son[i + i * sons],
			       xa1[k + i * x->rsons]);

	  for (j = i; j-- > 0;)
	    addproduct_haccum(-1.0, false, x->son[k + i * x->rsons], true,
			      a->son[j + i * sons], xa1[k + j * x->rsons]);
	}

      for (k = 0; k < x->rsons; k++)
	for (i = 0; i < sons; i++)
	  del_haccum(xa1[k + i * x->rsons]);
      freemem(xa1);
    }
  }
}

static void
uppersolve_tn_haccum(bool aunit, pchmatrix a, phaccum xa)
{
  phmatrix  atmp;
  phmatrix  x = xa->z;
  phaccum  *xa1;
  uint      sons;
  uint      i, j, k;

  assert(a->rc == a->cc);
  assert(a->cc == x->rc);

  if (x->f) {
    flush_haccum(xa);
    triangularinvmul_hmatrix_amatrix(false, aunit, true, a, false, x->f);
  }
  else if (x->r) {
    flush_haccum(xa);
    triangularinvmul_hmatrix_amatrix(false, aunit, true, a, false, &x->r->A);
  }
  else {
    if (a->f) {
      atmp = split_sub_amatrix(a->f, a->rc, a->cc, (x->son[0]->rc != x->rc),
			       (x->son[0]->rc != x->rc));

      assert(x->rsons == atmp->rsons);

      sons = atmp->rsons;

      xa1 = split_haccum(x, xa);

      for (k = 0; k < x->csons; k++)
	for (i = 0; i < sons; i++) {
	  uppersolve_tn_haccum(aunit, atmp->son[i + i * sons],
			       xa1[i + k * x->rsons]);

	  for (j = i + 1; j < sons; j++)
	    addproduct_haccum(-1.0, true, atmp->son[i + j * sons], false,
			      x->son[i + k * x->rsons],
			      xa1[j + k * x->rsons]);
	}

      for (k = 0; k < x->csons; k++)
	for (i = 0; i < sons; i++)
	  del_haccum(xa1[i + k * x->rsons]);
      freemem(xa1);

      del_hmatrix(atmp);
    }
    else {
      assert(a->rsons == a->csons);
      assert(x->rsons == a->rsons);

      sons = a->rsons;

      xa1 = split_haccum(x, xa);

      for (k = 0; k < x->csons; k++)
	for (i = 0; i < sons; i++) {
	  uppersolve_tn_haccum(aunit, a->son[i + i * sons],
			       xa1[i + k * x->rsons]);

	  for (j = i + 1; j < sons; j++)
	    addproduct_haccum(-1.0, true, a->son[i + j * sons], false,
			      x->son[i + k * x->rsons],
			      xa1[j + k * x->rsons]);
	}

      for (k = 0; k < x->csons; k++)
	for (i = 0; i < sons; i++)
	  del_haccum(xa1[i + k * x->rsons]);
      freemem(xa1);
    }
  }
}

static void
uppersolve_tt_haccum(bool aunit, pchmatrix a, phaccum xa)
{
  phmatrix  atmp;
  phmatrix  x = xa->z;
  phaccum  *xa1;
  uint      sons;
  uint      i, j, k;

  assert(a->rc == a->cc);
  assert(a->cc == x->cc);

  if (x->f) {
    flush_haccum(xa);
    triangularinvmul_hmatrix_amatrix(false, aunit, true, a, true, x->f);
  }
  else if (x->r) {
    flush_haccum(xa);
    triangularinvmul_hmatrix_amatrix(false, aunit, true, a, false, &x->r->B);
  }
  else {
    if (a->f) {
      atmp = split_sub_amatrix(a->f, a->rc, a->cc, (x->son[0]->cc != x->cc),
			       (x->son[0]->cc != x->cc));

      assert(x->csons == atmp->rsons);

      sons = atmp->rsons;

      xa1 = split_haccum(x, xa);

      for (k = 0; k < x->rsons; k++)
	for (i = 0; i < sons; i++) {
	  uppersolve_tt_haccum(aunit, atmp->son[i + i * sons],
			       xa1[k + i * x->rsons]);

	  for (j = i + 1; j < sons; j++)
	    addproduct_haccum(-1.0, false, x->son[k + i * x->rsons], false,
			      atmp->son[i + j * sons], xa1[k + j * x->rsons]);
	}

      for (k = 0; k < x->rsons; k++)
	for (i = 0; i < sons; i++)
	  del_haccum(xa1[k + i * x->rsons]);
      freemem(xa1);

      del_hmatrix(atmp);
    }
    else {
      assert(a->rsons == a->csons);
      assert(x->csons == a->rsons);

      sons = a->rsons;

      xa1 = split_haccum(x, xa);

      for (k = 0; k < x->rsons; k++)
	for (i = 0; i < sons; i++) {
	  uppersolve_tt_haccum(aunit, a->son[i + i * sons],
			       xa1[k + i * x->rsons]);

	  for (j = i + 1; j < sons; j++)
	    addproduct_haccum(-1.0, false, x->son[k + i * x->rsons], false,
			      a->son[i + j * sons], xa1[k + j * x->rsons]);
	}

      for (k = 0; k < x->rsons; k++)
	for (i = 0; i < sons; i++)
	  del_haccum(xa1[k + i * x->rsons]);
      freemem(xa1);
    }
  }
}

void
uppersolve_haccum(bool aunit, bool atrans, pchmatrix a, bool xtrans,
		  phaccum xa)
{
  if (atrans) {
    if (xtrans)
      uppersolve_tt_haccum(aunit, a, xa);
    else
      uppersolve_tn_haccum(aunit, a, xa);
  }
  else {
    if (xtrans)
      uppersolve_nt_haccum(aunit, a, xa);
    else
      uppersolve_nn_haccum(aunit, a, xa);
  }
}

void
uppersolve2_hmatrix_hmatrix(bool aunit, bool atrans, pchmatrix a,
			    pctruncmode tm, real eps, bool xtrans, phmatrix x)
{
  phaccum   xa;

  xa = new_haccum(x, tm, eps);

  uppersolve_haccum(aunit, atrans, a, xtrans, xa);

  del_haccum(xa);
}

/* ------------------------------------------------------------
 * Triangular factorizations
 * ------------------------------------------------------------ */

void
lrdecomp_haccum(phaccum aa)
{
  phmatrix  a = aa->z;
  phaccum  *aa1;
  uint      sons;
  uint      i, j, k;
  uint      res;

  assert(a->rc == a->cc);

  if (a->f) {
    flush_haccum(aa);
    res = lrdecomp_amatrix(a->f);
    assert(res == 0);
  }
  else {
    assert(a->son != 0);
    assert(a->rsons == a->csons);

    sons = a->rsons;

    aa1 = split_haccum(a, aa);

    for (k = 0; k < sons; k++) {
      lrdecomp_haccum(aa1[k + k * sons]);

      for (j = k + 1; j < sons; j++)
	lowersolve_nn_haccum(true, a->son[k + k * sons], aa1[k + j * sons]);

      for (i = k + 1; i < sons; i++)
	uppersolve_tt_haccum(false, a->son[k + k * sons], aa1[i + k * sons]);

      for (j = k + 1; j < sons; j++)
	for (i = k + 1; i < sons; i++)
	  addproduct_haccum(-1.0, false, a->son[i + k * sons], false,
			    a->son[k + j * sons], aa1[i + j * sons]);
    }

    for (k = 0; k < sons; k++)
      for (i = 0; i < sons; i++)
	del_haccum(aa1[i + k * sons]);
    freemem(aa1);
  }
}

void
lrdecomp2_hmatrix(phmatrix a, pctruncmode tm, real eps)
{
  phaccum   aa;

  aa = new_haccum(a, tm, eps);

  lrdecomp_haccum(aa);

  del_haccum(aa);
}

void
choldecomp_haccum(phaccum aa)
{
  phmatrix  a = aa->z;
  phaccum  *aa1;
  uint      sons;
  uint      i, j, k;
  uint      res;

  assert(a->rc == a->cc);

  if (a->f) {
    flush_haccum(aa);
    res = choldecomp_amatrix(a->f);
    assert(res == 0);
  }
  else {
    assert(a->son != 0);
    assert(a->rsons == a->csons);

    sons = a->rsons;

    aa1 = split_haccum(a, aa);

    for (k = 0; k < sons; k++) {
      choldecomp_haccum(aa1[k + k * sons]);

      for (i = k + 1; i < sons; i++)
	lowersolve_nt_haccum(false, a->son[k + k * sons], aa1[i + k * sons]);

      for (j = k + 1; j < sons; j++)
	for (i = j; i < sons; i++)
	  addproduct_haccum(-1.0, false, a->son[i + k * sons], true,
			    a->son[j + k * sons], aa1[i + j * sons]);
    }

    for (k = 0; k < sons; k++)
      for (i = 0; i < sons; i++)
	del_haccum(aa1[i + k * sons]);
    freemem(aa1);
  }
}

void
choldecomp2_hmatrix(phmatrix a, pctruncmode tm, real eps)
{
  phaccum   aa;

  aa = new_haccum(a, tm, eps);

  choldecomp_haccum(aa);

  del_haccum(aa);
}
