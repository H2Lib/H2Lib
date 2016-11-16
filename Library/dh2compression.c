
/* ------------------------------------------------------------
 * This is the file "dh2compression.h" of the H2Lib package.
 * All rights reserved, Steffen Boerm 2015
 * ------------------------------------------------------------ */

#include "dh2compression.h"

#include "basic.h"
#include "eigensolvers.h"
#include "factorizations.h"

/* ------------------------------------------------------------
 * Active blocks,
 * i.e., blocks that have to be approximated by the current cluster
 * basis
 * ------------------------------------------------------------ */

typedef struct _compactive compactive;
typedef compactive *pcompactive;

struct _compactive {
  pcamatrix G;			/* Source matrix */
  pcdblock  b;			/* Block */
  uint      iota;		/* Direction */
  amatrix   A;			/* Permuted or compressed submatrix */
  real      weight;		/* Weight factor */
  pcompactive next;		/* Next block in list */
};

/* ------------------------------------------------------------
 * Passive blocks,
 * i.e., blocks that do not have to be approximated by the current
 * cluster basis, but may have to be approximated by its descendants
 * ------------------------------------------------------------ */

typedef struct _comppassive comppassive;
typedef comppassive *pcomppassive;

struct _comppassive {
  pcamatrix G;
  pcdblock  b;
  pcomppassive next;
};

/* ------------------------------------------------------------
 * Destructors
 * ------------------------------------------------------------ */

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

/* ------------------------------------------------------------
 * Add a row block to the lists for a given cluster.
 * ------------------------------------------------------------ */

static void
addrow_comp(pcdcluster rc, pcamatrix G, pcdblock b, pctruncmode tm,
	    pcompactive * active, pcomppassive * passive)
{
  pcompactive ca;
  pcomppassive cp;
  pamatrix  Ahat;
  const uint *ridx, *cidx;
  size_t    ldA, ldG;
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
  else if (b->adm) {
    assert(b->rc == rc);

    /* Create an active block */
    rsize = b->rc->size;
    csize = b->cc->size;
    ca = (pcompactive) allocmem(sizeof(compactive));
    ca->G = G;
    ca->b = b;
    ca->iota = b->rd;
    Ahat = init_amatrix(&ca->A, rsize, csize);
    ca->next = *active;
    *active = ca;

    /* Copy entries from original matrix G */
    ldG = G->ld;
    ldA = Ahat->ld;
    ridx = b->rc->idx;
    cidx = b->cc->idx;
    for (j = 0; j < csize; j++)
      for (i = 0; i < rsize; i++)
	Ahat->a[i + j * ldA] = G->a[ridx[i] + cidx[j] * ldG];

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

/* ------------------------------------------------------------
 * Add a column block to the lists for a given cluster.
 * ------------------------------------------------------------ */

static void
addcol_comp(pcdcluster cc, pcamatrix G, pcdblock b, pctruncmode tm,
	    pcompactive * active, pcomppassive * passive)
{
  pcompactive ca;
  pcomppassive cp;
  pamatrix  Bhat;
  const uint *ridx, *cidx;
  size_t    ldB, ldG;
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
  else if (b->adm) {
    assert(b->cc == cc);

    /* Create an active block */
    rsize = b->rc->size;
    csize = b->cc->size;
    ca = (pcompactive) allocmem(sizeof(compactive));
    ca->G = G;
    ca->b = b;
    ca->iota = b->cd;
    Bhat = init_amatrix(&ca->A, csize, rsize);
    ca->next = *active;
    *active = ca;

    /* Copy entries from original matrix G */
    ldG = G->ld;
    ldB = Bhat->ld;
    ridx = b->rc->idx;
    cidx = b->cc->idx;
    for (i = 0; i < csize; i++)
      for (j = 0; j < rsize; j++)
	Bhat->a[i + j * ldB] = CONJ(G->a[ridx[j] + cidx[i] * ldG]);

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

typedef struct {
  bool      colbasis;
  pctruncmode tm;
  preal     eps;
  pcompactive *active;
  pcomppassive *passive;
} buildbasis_data;

static void
buildbasis_pre(pdclusterbasis cb, uint tname, uint pardepth, void *data)
{
  buildbasis_data *bd = (buildbasis_data *) data;
  bool      colbasis = bd->colbasis;
  pctruncmode tm = bd->tm;
  real      eps = bd->eps[tname];
  pcompactive active = bd->active[tname];
  pcomppassive passive = bd->passive[tname];
  pcdcluster t = cb->t;
  pcompactive active1, ca, ca1;
  pcomppassive passive1, cp;
  real      zeta_age, zeta_level;
  uint      i, off, tname1;

  (void) pardepth;

  zeta_age = (tm ? tm->zeta_age : 1.0);
  zeta_level = (tm ? tm->zeta_level : 1.0);

  if (cb->sons > 0) {
    assert(cb->sons == t->sons);

    off = 0;
    tname1 = tname + 1;
    for (i = 0; i < t->sons; i++) {
      active1 = 0;
      passive1 = 0;

      /* Check for passive blocks that become active in the son */
      if (colbasis) {
	for (cp = passive; cp; cp = cp->next)
	  addcol_comp(t->son[i], cp->G, cp->b, tm, &active1, &passive1);
      }
      else {
	for (cp = passive; cp; cp = cp->next)
	  addrow_comp(t->son[i], cp->G, cp->b, tm, &active1, &passive1);
      }

      /* Add submatrices for already active blocks to the list */
      for (ca = active; ca; ca = ca->next) {
	ca1 = (pcompactive) allocmem(sizeof(compactive));
	ca1->G = ca->G;
	ca1->b = ca->b;
	ca1->iota = cb->dirson[i][ca->iota];
	init_sub_amatrix(&ca1->A, &ca->A,
			 t->son[i]->size, off, ca->A.cols, 0);
	ca1->weight = ca->weight * zeta_age;
	ca1->next = active1;
	active1 = ca1;
      }

      bd->eps[tname1] = eps * zeta_level;
      bd->active[tname1] = active1;
      bd->passive[tname1] = passive1;

      off += cb->son[i]->t->size;
      tname1 += cb->son[i]->t->desc;
    }
    assert(off == cb->t->size);
    assert(tname1 == tname + cb->t->desc);
  }
}

static void
buildbasis_post(pdclusterbasis cb, uint tname, uint pardepth, void *data)
{
  amatrix   tmp1, tmp2, tmp3, tmp4;
  realavector tmp5;
  buildbasis_data *bd = (buildbasis_data *) data;
  pctruncmode tm = bd->tm;
  real      eps = bd->eps[tname];
  pcompactive active = bd->active[tname];
  pcdcluster t = cb->t;
  pamatrix  Ahat, Ahat0, Ahat1;
  pamatrix  Q, Q1;
  prealavector sigma;
  pcompactive ca;
  uint      directions;
  uint      iota, iota1;
  uint     *md;
  uint      i, off, m, n, k, tname1;

  (void) pardepth;

  directions = cb->directions;

  md = (uint *) allocmem(sizeof(uint) * directions);

  if (cb->sons > 0) {
    assert(cb->sons == t->sons);

    tname1 = tname + 1;
    for (i = 0; i < t->sons; i++) {
      /* Clean up sons' block lists */
      del_compactive(bd->active[tname1]);
      del_comppassive(bd->passive[tname1]);

      tname1 += cb->son[i]->t->desc;
    }
    assert(tname1 == tname + cb->t->desc);

    /* Compute number of rows for each direction */
    for (iota = 0; iota < directions; iota++) {
      m = 0;
      for (i = 0; i < cb->sons; i++) {
	iota1 = cb->dirson[i][iota];
	m += cb->son[i]->k[iota1];
      }

      md[iota] = m;
    }

    for (ca = active; ca; ca = ca->next) {
      assert(ca->iota < directions);

      iota = ca->iota;
      Ahat = init_amatrix(&tmp1, md[iota], ca->A.cols);

      /* Combine projections of son's submatrices */
      off = 0;
      m = 0;
      for (i = 0; i < t->sons; i++) {
	iota1 = cb->dirson[i][iota];

	/* Matrix prepared by recursive call */
	Ahat1 = init_sub_amatrix(&tmp2, &ca->A,
				 cb->son[i]->k[iota1], off, Ahat->cols, 0);

	/* Appropriate submatrix in Ahat */
	Ahat0 = init_sub_amatrix(&tmp3, Ahat,
				 cb->son[i]->k[iota1], m, Ahat->cols, 0);

	/* Copy submatrix to its place */
	copy_amatrix(false, Ahat1, Ahat0);
	uninit_amatrix(Ahat1);
	uninit_amatrix(Ahat0);

	off += t->son[i]->size;
	m += cb->son[i]->k[iota1];
      }
      assert(off == t->size);
      assert(m == md[iota]);

      Ahat0 = init_sub_amatrix(&tmp3, &ca->A, m, 0, Ahat->cols, 0);
      copy_amatrix(false, Ahat, Ahat0);
      uninit_amatrix(Ahat0);
      uninit_amatrix(Ahat);
    }
  }
  else {
    /* Set number of rows for each direction */
    m = t->size;
    for (iota = 0; iota < directions; iota++)
      md[iota] = m;
  }

  for (iota = 0; iota < directions; iota++) {
    /* Get number of rows of SVD matrix */
    m = md[iota];

    /* Determine number of columns of SVD matrix */
    n = 0;
    for (ca = active; ca; ca = ca->next)
      if (ca->iota == iota)
	n += ca->A.cols;

    /* Quick exit if zero columns */
    if (n == 0) {
      setrank_dclusterbasis(cb, iota, 0);

      continue;
    }

    /* Combine submatrices */
    Ahat = init_amatrix(&tmp1, m, n);
    n = 0;
    for (ca = active; ca; ca = ca->next)
      if (ca->iota == iota) {
	Ahat1 = init_sub_amatrix(&tmp2, Ahat, m, 0, ca->A.cols, n);
	Ahat0 = init_sub_amatrix(&tmp3, &ca->A, m, 0, ca->A.cols, 0);

	copy_amatrix(false, Ahat0, Ahat1);
	if (ca->weight != 1.0)
	  scale_amatrix(ca->weight, Ahat1);

	uninit_amatrix(Ahat0);
	uninit_amatrix(Ahat1);

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
    setrank_dclusterbasis(cb, iota, k);

    if (cb->sons == 0) {
      assert(cb->V[iota].rows == Q->rows);
      assert(cb->V[iota].cols == k);

      /* Copy cluster basis matrix to new cluster basis */
      Q1 = init_sub_amatrix(&tmp1, Q, Q->rows, 0, k, 0);
      copy_amatrix(false, Q1, cb->V + iota);
      uninit_amatrix(Q1);
    }
    else {
      off = 0;
      for (i = 0; i < cb->sons; i++) {
	iota1 = cb->dirson[i][iota];

	assert(cb->E[i][iota].rows == cb->son[i]->k[iota1]);
	assert(cb->E[i][iota].cols == k);

	/* Copy transfer matrix to new cluster basis */
	Q1 = init_sub_amatrix(&tmp1, Q, cb->son[i]->k[iota1], off, k, 0);
	copy_amatrix(false, Q1, cb->E[i] + iota);
	uninit_amatrix(Q1);

	off += cb->son[i]->k[iota1];
      }
      assert(off == m);
    }

    /* Prepare projected matrices for father */
    Q1 = init_sub_amatrix(&tmp3, Q, m, 0, k, 0);
    for (ca = active; ca; ca = ca->next)
      if (ca->iota == iota) {
	/* Copy original matrix */
	Ahat = init_amatrix(&tmp1, m, ca->A.cols);
	Ahat0 = init_sub_amatrix(&tmp4, &ca->A, m, 0, ca->A.cols, 0);
	copy_amatrix(false, Ahat0, Ahat);
	uninit_amatrix(Ahat0);

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
  }

  freemem(md);

  /* Now we're done with this subtree */
  update_dclusterbasis(cb);
}

pdclusterbasis
buildrow_amatrix_dclusterbasis(pcamatrix G, pcdblock b,
			       pctruncmode tm, real eps)
{
  pdclusterbasis rb;
  buildbasis_data bd;
  pcdcluster rc = b->rc;
  pcompactive active;
  pcomppassive passive;

  active = 0;
  passive = 0;
  addrow_comp(rc, G, b, tm, &active, &passive);

  bd.colbasis = false;
  bd.tm = tm;
  bd.eps = (preal) allocmem(sizeof(real) * rc->desc);
  bd.active = (pcompactive *) allocmem(sizeof(pcompactive) * rc->desc);
  bd.passive = (pcomppassive *) allocmem(sizeof(pcomppassive) * rc->desc);

  bd.eps[0] = eps;
  bd.active[0] = active;
  bd.passive[0] = passive;

  rb = buildfromdcluster_dclusterbasis(rc);

  iterate_dclusterbasis(rb, 0, max_pardepth,
			buildbasis_pre, buildbasis_post, &bd);

  freemem(bd.passive);
  freemem(bd.active);
  freemem(bd.eps);

  del_compactive(active);
  del_comppassive(passive);

  return rb;
}

pdclusterbasis
buildcol_amatrix_dclusterbasis(pcamatrix G, pcdblock b,
			       pctruncmode tm, real eps)
{
  pdclusterbasis cb;
  buildbasis_data bd;
  pcdcluster cc = b->cc;
  pcompactive active;
  pcomppassive passive;

  active = 0;
  passive = 0;
  addcol_comp(cc, G, b, tm, &active, &passive);

  bd.colbasis = true;
  bd.tm = tm;
  bd.eps = (preal) allocmem(sizeof(real) * cc->desc);
  bd.active = (pcompactive *) allocmem(sizeof(pcompactive) * cc->desc);
  bd.passive = (pcomppassive *) allocmem(sizeof(pcomppassive) * cc->desc);

  bd.eps[0] = eps;
  bd.active[0] = active;
  bd.passive[0] = passive;

  cb = buildfromdcluster_dclusterbasis(cc);

  iterate_dclusterbasis(cb, 0, max_pardepth,
			buildbasis_pre, buildbasis_post, &bd);

  freemem(bd.passive);
  freemem(bd.active);
  freemem(bd.eps);

  del_compactive(active);
  del_comppassive(passive);

  return cb;
}

/* ------------------------------------------------------------
   Orthogonal projection
   ------------------------------------------------------------ */

void
collectdense_dh2matrix(pcamatrix a,
		       pcdclusterbasis rb, uint rd,
		       pcdclusterbasis cb, uint cd, pamatrix s)
{
  amatrix   tmp1, tmp2;
  pamatrix  s1, s2;
  longindex lda = a->ld;
  longindex lds;
  pcdcluster rc, cc;
  uint      rd1, cd1;
  uint      i, j;

  assert(s->rows == rb->k[rd]);
  assert(s->cols == cb->k[cd]);

  clear_amatrix(s);

  if (rb->sons > 0) {
    if (cb->sons > 0) {		/* rb has sons, cb has sons */
      for (j = 0; j < cb->sons; j++) {
	cd1 = cb->dirson[j][cd];

	s1 = init_amatrix(&tmp1, rb->k[rd], cb->son[j]->k[cd1]);
	clear_amatrix(s1);

	for (i = 0; i < rb->sons; i++) {
	  rd1 = rb->dirson[i][rd];

	  s2 = init_amatrix(&tmp2, rb->son[i]->k[rd1], cb->son[j]->k[cd1]);

	  collectdense_dh2matrix(a, rb->son[i], rd1, cb->son[j], cd1, s2);

	  addmul_amatrix(1.0, true, rb->E[i] + rd, false, s2, s1);

	  uninit_amatrix(s2);
	}
	addmul_amatrix(1.0, false, s1, false, cb->E[j] + cd, s);

	uninit_amatrix(s1);
      }
    }
    else {			/* rb has sons, cb is a leaf */
      for (i = 0; i < rb->sons; i++) {
	rd1 = rb->dirson[i][rd];

	s2 = init_amatrix(&tmp2, rb->son[i]->k[rd1], cb->k[cd]);

	collectdense_dh2matrix(a, rb->son[i], rd1, cb, cd, s2);

	addmul_amatrix(1.0, true, rb->E[i] + rd, false, s2, s);

	uninit_amatrix(s2);
      }
    }
  }
  else {
    if (cb->sons > 0) {		/* rb is a leaf, cb has sons */
      for (j = 0; j < cb->sons; j++) {
	cd1 = cb->dirson[j][cd];

	s1 = init_amatrix(&tmp1, rb->k[rd], cb->son[j]->k[cd1]);

	collectdense_dh2matrix(a, rb, rd, cb->son[j], cd1, s1);

	addmul_amatrix(1.0, false, s1, false, cb->E[j] + cd, s);

	uninit_amatrix(s1);
      }
    }
    else {			/* rb is a leaf, cb is a leaf */
      rc = rb->t;
      cc = cb->t;

      s1 = init_amatrix(&tmp1, rb->k[rd], cc->size);
      clear_amatrix(s1);

      s2 = init_amatrix(&tmp2, rc->size, cc->size);
      lds = s2->ld;

      for (j = 0; j < cc->size; j++)
	for (i = 0; i < rc->size; i++)
	  s2->a[i + j * lds] = a->a[rc->idx[i] + cc->idx[j] * lda];

      addmul_amatrix(1.0, true, rb->V + rd, false, s2, s1);

      addmul_amatrix(1.0, false, s1, false, cb->V + cd, s);

      uninit_amatrix(s2);
      uninit_amatrix(s1);
    }
  }
}

typedef struct {
  pcamatrix a;
} projectdense_data;

static void
projectdense_pre(pdh2matrix h2, uint mname, uint rname, uint cname,
		 uint pardepth, void *data)
{
  projectdense_data *pd = (projectdense_data *) data;
  pcamatrix a = pd->a;
  pamatrix  f;
  longindex lda = a->ld;
  longindex ldf;
  pcdcluster rc, cc;
  const uint *ridx, *cidx;
  uint      rsize, csize;
  uint      i, j;

  (void) mname;
  (void) rname;
  (void) cname;
  (void) pardepth;

  if (h2->f) {
    f = h2->f;
    ldf = f->ld;
    rc = h2->rb->t;
    cc = h2->cb->t;
    ridx = rc->idx;
    cidx = cc->idx;
    rsize = rc->size;
    csize = cc->size;

    assert(h2->f->rows == rc->size);
    assert(h2->f->cols == cc->size);

    for (j = 0; j < csize; j++)
      for (i = 0; i < rsize; i++)
	f->a[i + j * ldf] = a->a[ridx[i] + cidx[j] * lda];
  }
  else if (h2->u)
    collectdense_dh2matrix(a, h2->rb, h2->u->rd, h2->cb, h2->u->cd,
			   &h2->u->S);
}

void
projectdense_dh2matrix(pcamatrix a, pdh2matrix h2)
{
  projectdense_data pd;

  pd.a = a;

  iterate_dh2matrix(h2, 0, 0, 0, max_pardepth, projectdense_pre, 0, &pd);
}
