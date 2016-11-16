
/* ------------------------------------------------------------
 * This is the file "clusteroperator.c" of the H2Lib package.
 * All rights reserved, Steffen Boerm 2009
 * ------------------------------------------------------------ */

#include <stdio.h>

#include "clusteroperator.h"
#include "basic.h"

static uint active_clusteroperator = 0;

/* ------------------------------------------------------------
 * Constructors and destructors
 * ------------------------------------------------------------ */

pclusteroperator
init_clusteroperator(pclusteroperator co, pccluster t)
{
  uint      i, sons;

  (void) init_leaf_clusteroperator(co, t);

  sons = t->sons;
  if (sons > 0) {
    co->sons = sons;
    co->son =
      (pclusteroperator *) allocmem((size_t) sizeof(pclusteroperator) * sons);
    for (i = 0; i < sons; i++)
      co->son[i] = NULL;
  }

  return co;
}

pclusteroperator
init_leaf_clusteroperator(pclusteroperator co, pccluster t)
{
  assert(co != NULL);

  co->t = t;
  co->krow = 0;
  co->kcol = 0;
  init_amatrix(&co->C, 0, 0);
  co->sons = 0;
  co->son = NULL;
  co->refs = 0;

#ifdef USE_OPENMP
#pragma omp atomic
#endif
  active_clusteroperator++;

  return co;
}

void
uninit_clusteroperator(pclusteroperator co)
{
  uint      i;

  assert(co->refs == 0);

  if (co->sons > 0) {
    for (i = 0; i < co->sons; i++)
      unref_clusteroperator(co->son[i]);
    freemem(co->son);
  }

  uninit_amatrix(&co->C);

  assert(active_clusteroperator > 0);

#ifdef USE_OPENMP
#pragma omp atomic
#endif
  active_clusteroperator--;
}

pclusteroperator
new_clusteroperator(pccluster t)
{
  pclusteroperator co;

  co = (pclusteroperator) allocmem(sizeof(clusteroperator));

  init_clusteroperator(co, t);

  return co;
}

pclusteroperator
new_leaf_clusteroperator(pccluster t)
{
  pclusteroperator co;

  co = (pclusteroperator) allocmem(sizeof(clusteroperator));

  init_leaf_clusteroperator(co, t);

  return co;
}

void
del_clusteroperator(pclusteroperator co)
{
  uninit_clusteroperator(co);

  freemem(co);
}

void
removesons_clusteroperator(pclusteroperator co)
{
  uint      i;

  for (i = 0; i < co->sons; i++)
    ref_clusteroperator(co->son + i, 0);

  freemem(co->son);

  co->son = 0;
  co->sons = 0;

  update_clusteroperator(co);
}

/* ------------------------------------------------------------
 * Reference counting
 * ------------------------------------------------------------ */

void
ref_clusteroperator(pclusteroperator * ptr, pclusteroperator co)
{
  if (*ptr)
    unref_clusteroperator(*ptr);

  *ptr = co;

  if (co)
#ifdef USE_OPENMP
#pragma omp atomic
#endif
    co->refs++;
}

void
unref_clusteroperator(pclusteroperator co)
{
  assert(co->refs > 0);

#ifdef USE_OPENMP
#pragma omp atomic
#endif
  co->refs--;

  if (co->refs == 0)
    del_clusteroperator(co);
}

/* ------------------------------------------------------------
 * Low-level management
 * ------------------------------------------------------------ */

void
update_clusteroperator(pclusteroperator co)
{
  (void) co;			/* nothing at the moment */
}

void
resize_clusteroperator(pclusteroperator co, uint krow, uint kcol)
{
  if (krow != co->krow || kcol != co->kcol) {
    resize_amatrix(&co->C, krow, kcol);

    co->krow = krow;
    co->kcol = kcol;
  }

  update_clusteroperator(co);
}

pclusteroperator
identify_son_clusterweight_clusteroperator(pcclusteroperator cwf, pccluster t)
{
  pclusteroperator cw;

  uint      i;

  cw = (pclusteroperator) cwf;
  i = 0;
  while (cw->t != t && i < cwf->sons) {
    cw = cwf->son[i];
    i++;
  }
  assert(cw->t == t);

  if (cw->sons == 0)
    cw = (pclusteroperator) cwf;

  return cw;
}

/* ------------------------------------------------------------
 * Build empty clusteroperators
 * ------------------------------------------------------------ */

pclusteroperator
build_from_cluster_clusteroperator(pccluster t)
{
  pclusteroperator co, co1;
  uint      i;

  co = new_clusteroperator(t);

  for (i = 0; i < co->sons; i++) {
    co1 = build_from_cluster_clusteroperator(t->son[i]);

    ref_clusteroperator(co->son + i, co1);
  }

  update_clusteroperator(co);

  return co;
}

pclusteroperator
build_from_clusterbasis_clusteroperator(pcclusterbasis cb)
{
  pclusteroperator co, co1;
  uint      i;

  co = (cb->son ?
	new_clusteroperator(cb->t) : new_leaf_clusteroperator(cb->t));

  resize_clusteroperator(co, 0, cb->k);

  for (i = 0; i < co->sons; i++) {
    co1 = build_from_clusterbasis_clusteroperator(cb->son[i]);

    ref_clusteroperator(co->son + i, co1);
  }

  update_clusteroperator(co);

  return co;
}

/* ------------------------------------------------------------
 * Statistics
 * ------------------------------------------------------------ */

uint
getactives_clusteroperator()
{
  return active_clusteroperator;
}

size_t
getsize_clusteroperator(pcclusteroperator co)
{
  size_t    sz;
  uint      i;

  sz = (size_t) sizeof(clusteroperator);
  sz += getsize_heap_amatrix(&co->C);

  if (co->sons > 0) {
    sz += (size_t) sizeof(pclusteroperator) * co->sons;
    for (i = 0; i < co->sons; i++)
      sz += getsize_clusteroperator(co->son[i]);
  }

  return sz;
}

/* ------------------------------------------------------------
 * Debugging
 * ------------------------------------------------------------ */

static void
print_tree(pcclusteroperator co, uint level)
{
  uint      i;

  for (i = 0; i < level; i++)
    printf("  ");

  printf("%u %u %.4g\n", co->krow, co->kcol, norm2_amatrix(&(co->C)));

  for (i = 0; i < co->sons; i++)
    print_tree(co->son[i], level + 1);
}

void
print_tree_clusteroperator(pcclusteroperator co)
{
  print_tree(co, 0);
}

static void
norm2diff(pcclusteroperator co1, pcclusteroperator co2, uint level)
{
  uint      i;

  for (i = 0; i < level; i++)
    printf("  ");

  if (co1->sons != co2->sons)
    printf("Tree mismatch ");

  if (co1->krow != co2->krow || co1->kcol != co2->kcol)
    printf("Dimension mismatch\n");
  else
    printf("%g %g\n",
	   norm2_amatrix(&co1->C), norm2diff_amatrix(&co1->C, &co2->C));

  if (co1->sons == co2->sons)
    for (i = 0; i < co1->sons; i++)
      norm2diff(co1->son[i], co2->son[i], level + 1);
}

void
norm2diff_clusteroperator(pcclusteroperator co1, pcclusteroperator co2)
{
  norm2diff(co1, co2, 0);
}

static    real
compareweights(pcclusteroperator co1, pcclusteroperator co2, uint level)
{
  amatrix   tmp;
  pamatrix  D;
  real      norm, error, error1;
  uint      k;
  uint      i;

  k = co1->kcol;
  assert(k == co2->kcol);

  D = init_amatrix(&tmp, k, k);

  clear_amatrix(D);
  addmul_amatrix(1.0, true, &co1->C, false, &co1->C, D);
  norm = norm2_amatrix(D);
  addmul_amatrix(-1.0, true, &co2->C, false, &co2->C, D);
  error = norm2_amatrix(D);

  uninit_amatrix(D);

  if (norm > 0.0)
    error /= norm;

  if (co1->sons != co2->sons)
    printf("Tree mismatch");
  else
    for (i = 0; i < co1->sons; i++) {
      error1 = compareweights(co1->son[i], co2->son[i], level + 1);
      if (error1 > error)
	error = error1;
    }

  return error;
}

real
compareweights_clusteroperator(pcclusteroperator co1, pcclusteroperator co2)
{
  return compareweights(co1, co2, 0);
}

/* ------------------------------------------------------------
 * Enumeration
 * ------------------------------------------------------------ */

static void
enumerate(pccluster t, uint tname, pclusteroperator co,
	  pclusteroperator * con)
{
  uint      tname1;
  uint      i;

  assert(co == 0 || co->t == t);

  con[tname] = co;

  tname1 = tname + 1;

  if (co == 0 || co->son == 0)
    for (i = 0; i < t->sons; i++) {
      enumerate(t->son[i], tname1, 0, con);

      tname1 += t->son[i]->desc;
    }
  else {
    assert(t->sons == co->sons);

    for (i = 0; i < t->sons; i++) {
      enumerate(t->son[i], tname1, co->son[i], con);

      tname1 += t->son[i]->desc;
    }
  }
  assert(tname1 == tname + t->desc);
}

pclusteroperator *
enumerate_clusteroperator(pccluster t, pclusteroperator co)
{
  pclusteroperator *con;

  con = allocmem((size_t) sizeof(pclusteroperator) * co->t->desc);

  enumerate(t, 0, co, con);

  return con;
}

/* ------------------------------------------------------------
 * Cluster basis product
 * ------------------------------------------------------------ */

void
basisproduct_clusteroperator(pcclusterbasis cb1, pcclusterbasis cb2,
			     pclusteroperator pr)
{
  pamatrix  X, Xt;
  amatrix   tmp1, tmp2;
  uint      i;

  assert(cb1->t == pr->t);
  assert(cb2->t == pr->t);

  resize_clusteroperator(pr, cb1->k, cb2->k);

  if (pr->sons > 0) {
    assert(cb1->sons == pr->sons);
    assert(cb2->sons == pr->sons);

    for (i = 0; i < pr->sons; i++)
      basisproduct_clusteroperator(cb1->son[i], cb2->son[i], pr->son[i]);

    clear_amatrix(&pr->C);
    for (i = 0; i < pr->sons; i++) {
      X = init_amatrix(&tmp1, cb1->son[i]->k, cb2->k);
      clear_amatrix(X);
      addmul_amatrix(1.0, false, &pr->son[i]->C, false, &cb2->son[i]->E, X);
      addmul_amatrix(1.0, true, &cb1->son[i]->E, false, X, &pr->C);
      uninit_amatrix(X);
    }
  }
  else {
    if (cb1->sons == 0) {
      if (cb2->sons == 0) {
	clear_amatrix(&pr->C);
	addmul_amatrix(1.0, true, &cb1->V, false, &cb2->V, &pr->C);
      }
      else {
	Xt = init_amatrix(&tmp1, cb2->kbranch, cb1->k);
	compress_clusterbasis_amatrix(cb2, &cb1->V, Xt);
	X = init_sub_amatrix(&tmp2, Xt, cb2->k, 0, cb1->k, 0);
	copy_amatrix(true, X, &pr->C);
	uninit_amatrix(X);
	uninit_amatrix(Xt);
      }
    }
    else {
      assert(cb2->sons == 0);	/* Could be generalized */

      Xt = init_amatrix(&tmp1, cb1->kbranch, cb2->k);
      compress_clusterbasis_amatrix(cb1, &cb2->V, Xt);
      X = init_sub_amatrix(&tmp2, Xt, cb1->k, 0, cb2->k, 0);
      copy_amatrix(false, X, &pr->C);
      uninit_amatrix(X);
      uninit_amatrix(Xt);
    }
  }
}
