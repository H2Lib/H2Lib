/* ------------------------------------------------------------
 This is the file "oclbem3d.c" of the H2Lib package.
 All rights reserved, Sven Christophersen 2015
 ------------------------------------------------------------ */

/**
 * @file oclbem3d.c
 * @author Sven Christophersen
 * @date 2015
 */

#ifdef USE_OPENMP
#ifdef USE_OPENCL

/* C STD LIBRARY */
/* CORE 0 */
/* CORE 1 */
/* CORE 2 */
/* CORE 3 */
/* SIMPLE */
/* PARTICLES */
/* BEM */
#include "oclbem3d.h"

bem3d_ocl ocl_bem3d = (bem3d_ocl) {
  .num_kernels = 0
};

ptaskgroup nf = NULL;
omp_lock_t nf_lock;
ptaskgroup nf_dist = NULL;
omp_lock_t nf_dist_lock;
ptaskgroup nf_vert = NULL;
omp_lock_t nf_vert_lock;
ptaskgroup nf_edge = NULL;
omp_lock_t nf_edge_lock;
ptaskgroup nf_iden = NULL;
omp_lock_t nf_iden_lock;

op_callbacks op_cb;

size_t
getsize_nf(void *data)
{
  nearfield_args *nf_args;

  nf_args = (nearfield_args *) data;

  return nf_args->N->rows * nf_args->N->cols * sizeof(field);
}

void
cleanup_nf_task(void *data)
{
  nearfield_args *nf_args = (nearfield_args *) data;

  if (nf_args->N->owner != NULL) {
    del_amatrix(nf_args->N);
  }
  freemem(data);
}

void
cleanup_nf_merge(void *data)
{
  merge_data_nf *mdata = (merge_data_nf *) data;

  freemem(mdata->N);
  freemem(mdata->ridx);
  freemem(mdata->cidx);
  if (mdata->addr != NULL) {
    freemem(mdata->addr);
  }
  freemem(mdata);
}

void
split_nf_task(void *data, void ***split, uint * n)
{
  nearfield_args *nf_args = (nearfield_args *) data;
  uint      r = nf_args->N->rows;
  uint      c = nf_args->N->cols;

  nearfield_args **nfs;

  nfs = (nearfield_args **) allocmem(4 * sizeof(nearfield_args *));
  *split = (void **) nfs;

  nfs[0] = (nearfield_args *) allocmem(sizeof(nearfield_args));
  nfs[0]->N = new_sub_amatrix(nf_args->N, r / 2, 0, c / 2, 0);
  nfs[0]->bem = nf_args->bem;
  nfs[0]->cidx = nf_args->cidx;
  nfs[0]->ridx = nf_args->ridx;
  nfs[0]->dist = nf_args->dist;
  nfs[0]->ntrans = nf_args->ntrans;

  nfs[1] = (nearfield_args *) allocmem(sizeof(nearfield_args));
  nfs[1]->N = new_sub_amatrix(nf_args->N, r - r / 2, r / 2, c / 2, 0);
  nfs[1]->bem = nf_args->bem;
  nfs[1]->cidx = nf_args->cidx;
  nfs[1]->ridx = nf_args->ridx + r / 2;
  nfs[1]->dist = nf_args->dist;
  nfs[1]->ntrans = nf_args->ntrans;

  nfs[2] = (nearfield_args *) allocmem(sizeof(nearfield_args));
  nfs[2]->N = new_sub_amatrix(nf_args->N, r / 2, 0, c - c / 2, c / 2);
  nfs[2]->bem = nf_args->bem;
  nfs[2]->cidx = nf_args->cidx + c / 2;
  nfs[2]->ridx = nf_args->ridx;
  nfs[2]->dist = nf_args->dist;
  nfs[2]->ntrans = nf_args->ntrans;

  nfs[3] = (nearfield_args *) allocmem(sizeof(nearfield_args));
  nfs[3]->N = new_sub_amatrix(nf_args->N, r - r / 2, r / 2, c - c / 2, c / 2);
  nfs[3]->bem = nf_args->bem;
  nfs[3]->cidx = nf_args->cidx + c / 2;
  nfs[3]->ridx = nf_args->ridx + r / 2;
  nfs[3]->dist = nf_args->dist;
  nfs[3]->ntrans = nf_args->ntrans;

  *n = 4;
}

size_t
getsize_nf_sing(void *data)
{
  merge_data_nf *nf_se = (merge_data_nf *) data;

  return nf_se->pos * sizeof(field);
}

void
cleanup_nf_sing_task(void *data)
{
  merge_data_nf *mdata = (merge_data_nf *) data;

  freemem(mdata->N);
  freemem(mdata->ridx);
  freemem(mdata->cidx);
  freemem(mdata->trans);
  freemem(mdata->addr);
  freemem(mdata);
}

void
cleanup_nf_sing_merge(void *data)
{
  (void) data;
}

void
split_nf_sing_task(void *data, void ***split, uint * n)
{
  (void) data;
  (void) split;
  (void) n;

  assert(0);
}

void
close_nf_dist(ptaskgroup tg)
{

  (void) tg;

  omp_set_lock(&nf_dist_lock);
  nf_dist = NULL;
  omp_unset_lock(&nf_dist_lock);
}

void
close_nf_vert(ptaskgroup tg)
{

  (void) tg;

  omp_set_lock(&nf_vert_lock);
  nf_vert = NULL;
  omp_unset_lock(&nf_vert_lock);
}

void
close_nf_edge(ptaskgroup tg)
{

  (void) tg;

  omp_set_lock(&nf_edge_lock);
  nf_edge = NULL;
  omp_unset_lock(&nf_edge_lock);
}

void
close_nf_iden(ptaskgroup tg)
{

  (void) tg;

  omp_set_lock(&nf_iden_lock);
  nf_iden = NULL;
  omp_unset_lock(&nf_iden_lock);
}

void
merge_nf_sing(ptaskgroup tg, void **data)
{
  assert(tg != NULL);
  assert(tg->tasks != NULL);
  assert(tg->tasks->next == NULL);

  *data = tg->tasks->data;
}

void
distribute_nf_sing(ptaskgroup tg, void *data)
{
  merge_data_nf *nf_se = (merge_data_nf *) data;
  uint      i;

  (void) tg;

  for (i = 0; i < nf_se->pos; ++i) {
    *(nf_se->addr[i]) = nf_se->trans[i] ? CONJ(nf_se->N[i]) : nf_se->N[i];
  }
}

void
add_nf_entry(pcbem3d bem, op_callbacks * op_cb, uint c, uint tt, uint ss,
	     bool trans, field * entry)
{
  merge_data_nf *mdata;
  ptaskgroup *tg;
  void      (*callback) (void *data);
  void      (*close_taskgroup) (ptaskgroup tg);
  void      (*cpu_wrapper) (void *data);
  omp_lock_t *lock;
  size_t    pos;
  real      factor;

  switch (c) {
  case 0:
    lock = &nf_dist_lock;
    break;
  case 1:
    lock = &nf_vert_lock;
    break;
  case 2:
    lock = &nf_edge_lock;
    break;
  case 3:
    lock = &nf_iden_lock;
    break;
  default:
    fprintf(stderr, "Unknown quadrature case!\n");
    abort();
    break;
  }

  omp_set_lock(lock);

  switch (c) {
  case 0:
    tg = &nf_dist;
    callback = op_cb->dist;
    close_taskgroup = close_nf_dist;
    break;
  case 1:
    tg = &nf_vert;
    callback = op_cb->vert;
    close_taskgroup = close_nf_vert;
    break;
  case 2:
    tg = &nf_edge;
    callback = op_cb->edge;
    close_taskgroup = close_nf_edge;
    break;
  case 3:
    tg = &nf_iden;
    callback = op_cb->iden;
    close_taskgroup = close_nf_iden;
    break;
  default:
    fprintf(stderr, "Unknown quadrature case!\n");
    abort();
    break;
  }

  cpu_wrapper = op_cb->sing_cpu;

  if (*tg == NULL) {
    *tg =
      new_ocltaskgroup(GPU_ONLY, NULL, merge_nf_sing, cleanup_nf_sing_merge,
		       distribute_nf_sing, close_taskgroup, cpu_wrapper,
		       callback, getsize_nf_sing, split_nf_sing_task,
		       cleanup_nf_sing_task, op_cb);

    mdata = (merge_data_nf *) allocmem(sizeof(merge_data_nf));
    mdata->N = (field *) allocmem(ocl_system.max_package_size);
    mdata->ridx =
      (uint *) allocmem(ocl_system.max_package_size / sizeof(field) *
			sizeof(uint));
    mdata->cidx =
      (uint *) allocmem(ocl_system.max_package_size / sizeof(field) *
			sizeof(uint));
    mdata->trans =
      (bool *) allocmem(ocl_system.max_package_size / sizeof(field) *
			sizeof(bool));
    mdata->addr =
      (field **) allocmem(ocl_system.max_package_size / sizeof(field) *
			  sizeof(field *));
    mdata->pos = 0;
    mdata->bem = bem;

    add_task_taskgroup(tg, (void *) mdata);
  }

  mdata = (merge_data_nf *) (*tg)->tasks->data;
  pos = mdata->pos;

  assert(entry != NULL);
  mdata->addr[pos] = entry;
  mdata->ridx[pos] = tt;
  mdata->cidx[pos] = ss;
  mdata->trans[pos] = trans;
  mdata->pos++;
  (*tg)->size += sizeof(field);

  factor = (real) bem->sq->q / (real) bem->sq->q2;
  factor *= factor;
  factor *= 0.5 * factor;

  if ((*tg)->size >= ocl_system.max_package_size * factor) {
    add_taskgroup_to_scheduler(*tg);
    *tg = NULL;
  }

  omp_unset_lock(lock);
}

void
close_nf(ptaskgroup tg)
{
  (void) tg;

  omp_set_lock(&nf_lock);
  nf = NULL;
  omp_unset_lock(&nf_lock);
}

void
merge_nf(ptaskgroup tg, void **data)
{
  ptask     t2;
  uint      r, c;
  nearfield_args *nf_args;
  merge_data_nf *mdata;
  uint      i, j, tt, ss;

  mdata = (merge_data_nf *) allocmem(sizeof(merge_data_nf));
  mdata->N = (field *) allocmem(ocl_system.max_package_size);
  mdata->ridx =
    (uint *) allocmem(ocl_system.max_package_size / sizeof(field) *
		      sizeof(uint));
  mdata->cidx =
    (uint *) allocmem(ocl_system.max_package_size / sizeof(field) *
		      sizeof(uint));
  mdata->addr = NULL;
  mdata->pos = 0;
  mdata->bem = ((nearfield_args *) tg->tasks->data)->bem;

  t2 = tg->tasks;
  while (t2 != NULL) {
    nf_args = (nearfield_args *) t2->data;

    if (nf_args->ntrans) {
      c = nf_args->N->rows;
      r = nf_args->N->cols;

      for (i = 0; i < r; ++i) {
	tt = nf_args->ridx[i];
	for (j = 0; j < c; ++j) {
	  ss = nf_args->cidx[j];

	  mdata->ridx[mdata->pos] = tt;
	  mdata->cidx[mdata->pos] = ss;
	  mdata->pos++;
	}
      }
    }
    else {
      r = nf_args->N->rows;
      c = nf_args->N->cols;

      for (j = 0; j < c; ++j) {
	ss = nf_args->cidx[j];
	for (i = 0; i < r; ++i) {
	  tt = nf_args->ridx[i];

	  mdata->ridx[mdata->pos] = tt;
	  mdata->cidx[mdata->pos] = ss;
	  mdata->pos++;
	}
      }

    }

    t2 = t2->next;
  }

  *data = mdata;
}

void
distribute_nf(ptaskgroup tg, void *data)
{
  merge_data_nf *mdata = (merge_data_nf *) data;

  uint(*geo_t)[3];
  ptask     t2;
  uint      r, c;
  nearfield_args *nf_args;
  field    *aa;
  longindex ld;
  uint      i, j, tt, ss, l;
  bool      trans;
  uint      pos;

  assert(tg->data != NULL);

  pos = 0;
  t2 = tg->tasks;
  while (t2 != NULL) {
    nf_args = (nearfield_args *) t2->data;
    aa = nf_args->N->a;
    ld = nf_args->N->ld;
    geo_t = nf_args->bem->gr->t;
    trans = nf_args->ntrans;
    r = nf_args->N->rows;
    c = nf_args->N->cols;

    for (j = 0; j < c; ++j) {
      memcpy(aa + j * ld, mdata->N + pos, r * sizeof(field));
      pos += r;
    }
    if (trans) {
      conjugate_amatrix(nf_args->N);
    }

    /****************************************************
     * Row and column indices
     ****************************************************/

    if (!nf_args->dist) {

      if (trans) {
	c = nf_args->N->rows;
	r = nf_args->N->cols;
      }

      for (j = 0; j < c; ++j) {
	for (i = 0; i < r; ++i) {
	  tt = nf_args->ridx[i];
	  ss = nf_args->cidx[j];

	  l = fast_select_quadrature(geo_t, tt, ss);

	  if (l > 0) {
	    add_nf_entry(nf_args->bem, (op_callbacks *) tg->data, l, tt, ss,
			 trans, trans ? aa + j + i * ld : aa + i + j * ld);
	  }
	}
      }
    }

    t2 = t2->next;
  }
}

#endif
#endif
