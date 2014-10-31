
/* ------------------------------------------------------------
   This is the file "cluster.c" of the H2Lib package.
   All rights reserved, Knut Reimer 2009
   ------------------------------------------------------------ */

#include <assert.h>
#include <math.h>

#include "basic.h"
#include "cluster.h"

/* ------------------------------------------------------------
   Constructors and destructors
   ------------------------------------------------------------ */

pcluster
new_cluster(uint size, uint * idx, uint sons, uint dim)
{
  pcluster  t;

  uint      i;

  t = (pcluster) allocmem((size_t) sizeof(cluster));
  t->size = size;
  t->dim = dim;
  t->idx = idx;
  t->bmin = allocreal(dim);
  t->bmax = allocreal(dim);
  t->sons = sons;
  if (sons > 0) {
    t->son = (pcluster *) allocmem((size_t) sons * sizeof(pcluster));
    for (i = 0; i < sons; i++)
      t->son[i] = 0;
  }
  else
    t->son = 0;

  return t;
}


void
del_cluster(pcluster t)
{
  uint      i;

  if (t->sons > 0) {
    for (i = 0; i < t->sons; i++)
      del_cluster(t->son[i]);
    freemem(t->son);
  }
  freemem(t->bmax);
  freemem(t->bmin);
  freemem(t);
}


void
update_cluster(pcluster t)
{
  uint      i, desc;

  desc = 1;
  for (i = 0; i < t->sons; i++)
    desc += t->son[i]->desc;

  t->desc = desc;
}


/* ------------------------------------------------------------
   Statistics
   ------------------------------------------------------------ */

uint
getdepth_cluster(pccluster t)
{
  uint      l, lt;
  uint      i;

  l = 0;
  if (t->son) {
    for (i = 0; i < t->sons; i++) {
      lt = getdepth_cluster(t->son[i]);
      l = UINT_MAX(lt, l);
    }
    l++;
  }

  return l;
}

uint
getmindepth_cluster(pccluster t)
{
  uint      l, lt;
  uint      i;

  l = 0;
  if (t->son) {
    l = getmindepth_cluster(t->son[0]);
    for (i = 1; i < t->sons; i++) {
      lt = getmindepth_cluster(t->son[i]);
      l = UINT_MIN(lt, l);
    }
    l++;
  }

  return l;
}

/* ------------------------------------------------------------
   Structure Adaptation
   ------------------------------------------------------------ */

void
extend_cluster(pcluster t, uint depth)
{
  uint      i;

  if ((t->son) && (depth > 0))
    for (i = 0; i < t->sons; i++)
      extend_cluster(t->son[i], depth - 1);

  else if ((!t->son) && (depth > 0)) {
    t->sons = 1;
    t->son = (pcluster *) allocmem(sizeof(pcluster));
    t->son[0] = new_cluster(t->size, t->idx, 0, t->dim);
    for (i = 0; i < t->dim; i++) {
      t->son[0]->bmin[i] = t->bmin[i];
      t->son[0]->bmax[i] = t->bmax[i];
    }
    extend_cluster(t->son[0], depth - 1);
  }
  update_cluster(t);
}

void
cut_cluster(pcluster t, uint depth)
{
  uint      i;

  if ((t->son) && (depth > 0))
    for (i = 0; i < t->sons; i++)
      cut_cluster(t->son[i], depth - 1);

  else if ((t->son) && (depth == 0)) {
    for (i = 0; i < t->sons; i++)
      del_cluster(t->son[i]);
    t->sons = 0;
    freemem(t->son);
    t->son = NULL;
  }
  update_cluster(t);
}

void
balance_cluster(pcluster t, uint depth)
{
  uint      i;

  if ((t->son) && (depth > 0))
    for (i = 0; i < t->sons; i++)
      balance_cluster(t->son[i], depth - 1);

  else if ((t->son) && (depth == 0))
    cut_cluster(t, depth);

  else if ((!t->son) && (depth > 0))
    extend_cluster(t, depth);

  update_cluster(t);
}

void
coarsen_cluster(pcluster t, uint minsize)
{
  bool      a = false;
  uint      i;

  for (i = 0; i < t->sons; i++) {
    if (t->son[i]->size < minsize) {
      a = true;
      break;
    }
  }

  if (a) {
    for (i = 0; i < t->sons; i++)
      del_cluster(t->son[i]);
    t->sons = 0;
    freemem(t->son);
    t->son = NULL;
  }
  else
    for (i = 0; i < t->sons; i++)
      coarsen_cluster(t->son[i], minsize);

  update_cluster(t);
}

void
setsons_cluster(pcluster t, uint sons)
{
  uint      i;

  assert(t->sons == 0);

  t->sons = sons;
  t->son = (pcluster *) allocmem(sizeof(pcluster) * sons);

  for (i = 0; i < sons; i++)
    t->son[i] = 0;
}

/* ------------------------------------------------------------
   Bounding box
   ------------------------------------------------------------ */

real
getdiam_2_cluster(pccluster t)
{
  real      diam2;
  uint      i;

  diam2 = 0.0;
  for (i = 0; i < t->dim; i++)
    diam2 += REAL_SQR(t->bmax[i] - t->bmin[i]);

  return REAL_SQRT(diam2);
}


real
getdist_2_cluster(pccluster t, pccluster s)
{
  real      dist2;
  uint      i;

  assert(t->dim == s->dim);

  dist2 = 0.0;
  for (i = 0; i < t->dim; i++)
    if (t->bmax[i] < s->bmin[i])	/* t to the "left" of s */
      dist2 += REAL_SQR(s->bmin[i] - t->bmax[i]);
    else if (s->bmax[i] < t->bmin[i])	/* t to the "right" of s */
      dist2 += REAL_SQR(t->bmin[i] - s->bmax[i]);

  return REAL_SQRT(dist2);
}

real
getdiam_max_cluster(pccluster t)
{

  real      diam_max, diam;
  uint      i;

  diam_max = 0.0;
  for (i = 0; i < t->dim; i++) {
    diam = t->bmax[i] - t->bmin[i];
    if (diam > diam_max) {
      diam_max = diam;
    }
  }

  return diam_max;
}

real
getdist_max_cluster(pccluster t, pccluster s)
{

  real      dist_max, dist;
  uint      i;

  assert(t->dim == s->dim);

  dist_max = 0.0;
  for (i = 0; i < t->dim; i++) {
    if (t->bmax[i] < s->bmin[i]) {
      dist = s->bmin[i] - t->bmax[i];
      if (dist > dist_max) {
	dist_max = dist;
      }
    }
    else if (s->bmax[i] < t->bmin[i]) {
      dist = t->bmin[i] - s->bmax[i];
      if (dist > dist_max) {
	dist_max = dist;
      }
    }
  }

  return dist_max;
}

/* ------------------------------------------------------------
   Hierarchical iterator
   ------------------------------------------------------------ */

void
iterate_cluster(pccluster t, uint tname,
		void (*pre) (pccluster t, uint tname, void *data),
		void (*post) (pccluster t, uint tname, void *data),
		void *data)
{
  uint      tname1;
  uint      i;

  if (pre)
    pre(t, tname, data);

  tname1 = tname + 1;
  for (i = 0; i < t->sons; i++) {
    iterate_cluster(t->son[i], tname1, pre, post, data);

    tname1 += t->son[i]->desc;
  }
  assert(tname1 == tname + t->desc);

  if (post)
    post(t, tname, data);
}

void
iterate_parallel_cluster(pccluster t, uint tname, uint pardepth,
			 void (*pre) (pccluster t, uint tname, void *data),
			 void (*post) (pccluster t, uint tname, void *data),
			 void *data)
{
  uint     *tname1, tname2;
#ifdef USE_OPENMP
  uint      nthreads;		/* HACK: Solaris workaround */
#endif
  uint      i;

  if (pre)
    pre(t, tname, data);

  if (t->sons > 0) {
    tname1 = (uint *) allocmem((size_t) sizeof(uint) * t->sons);

    tname2 = tname + 1;
    for (i = 0; i < t->sons; i++) {
      tname1[i] = tname2;
      tname2 += t->son[i]->desc;
    }
    assert(tname2 == tname + t->desc);

#ifdef USE_OPENMP
    nthreads = t->sons;
    (void) nthreads;
#pragma omp parallel for if(pardepth > 0), num_threads(nthreads)
#endif
    for (i = 0; i < t->sons; i++)
      iterate_parallel_cluster(t->son[i], tname1[i],
			       (pardepth > 0 ? pardepth - 1 : 0),
			       pre, post, data);

    freemem(tname1);
  }

  if (post)
    post(t, tname, data);
}

/* ------------------------------------------------------------
   Enumeration
   ------------------------------------------------------------ */

static void
enumerate(pcluster t, uint tname, pcluster *tn)
{
  uint      tname1;
  uint      i;

  tn[tname] = t;

  tname1 = tname + 1;
  for (i = 0; i < t->sons; i++) {
    enumerate(t->son[i], tname1, tn);

    tname1 += t->son[i]->desc;
  }
  assert(tname1 == tname + t->desc);
}

pcluster *
enumerate_cluster(pcluster t)
{
  pcluster *tn;

  tn = (pcluster *) allocmem((size_t) sizeof(pcluster) * t->desc);

  enumerate(t, 0, tn);

  return tn;
}
