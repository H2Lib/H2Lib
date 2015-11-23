/* ------------------------------------------------------------
 This is the file "hcoarsen.c" of the H2Lib package.
 All rights reserved, Sven Christophersen 2014
 ------------------------------------------------------------ */

#include "hcoarsen.h"

void
coarsen_hmatrix(phmatrix G, ptruncmode tm, real eps, bool recursive)
{
  uint      rsons = G->rsons;
  uint      csons = G->csons;

  phmatrix  son;
  prkmatrix R;
  pamatrix  A, B;
  amatrix   T, S;
  uint      i, j, leafs, ranksum, rankoffset, rowoffset, coloffset, rank;
  size_t    sizeold, sizenew;

  leafs = 0;

  /* recursion */
  if (rsons * csons > 0) {
    for (j = 0; j < csons; ++j) {
      for (i = 0; i < rsons; ++i) {
	son = G->son[i + j * rsons];
	if (recursive == true) {
	  coarsen_hmatrix(son, tm, eps, recursive);
	}
	leafs += son->rsons * son->csons;
      }
    }

    update_hmatrix(G);
  }
  else {
    /* matrix is a leaf -> northing to do */
    return;
  }

  if (leafs > 0) {
    /* matrix has sons which are not leafs  -> nothing to do */
    return;
  }
  else {
    if (G->rc == G->cc) {
      return;
    }

    /* determine ranksum and size of sons */
    ranksum = 0;
    sizeold = 0;
    for (j = 0; j < csons; ++j) {
      for (i = 0; i < rsons; ++i) {
	son = G->son[i + j * rsons];
	if (son->r) {
	  ranksum += son->r->k;
	  sizeold += getsize_rkmatrix(son->r);
	}
	else {
	  assert(son->f != NULL);
	  ranksum += son->f->cols;
	  sizeold += getsize_amatrix(son->f);
	}
      }
    }

    /* new rank-k-matrix */
    R = new_rkmatrix(G->rc->size, G->cc->size, ranksum);
    A = &R->A;
    B = &R->B;
    clear_amatrix(A);
    clear_amatrix(B);

    /* copy sons into a big rank-k-matrix */
    rankoffset = 0;
    coloffset = 0;
    for (j = 0; j < csons; ++j) {
      rowoffset = 0;
      for (i = 0; i < rsons; ++i) {
	son = G->son[i + j * rsons];
	rank = son->r ? son->r->k : son->f->cols;

	init_sub_amatrix(&T, A, son->rc->size, rowoffset, rank, rankoffset);
	init_sub_amatrix(&S, B, son->cc->size, coloffset, rank, rankoffset);

	if (son->r) {
	  copy_amatrix(false, &(son->r->A), &T);
	  copy_amatrix(false, &(son->r->B), &S);
	}
	else {
	  copy_amatrix(false, son->f, &T);
	  identity_amatrix(&S);
	}

	rankoffset += rank;
	rowoffset += son->rc->size;
	uninit_amatrix(&T);
	uninit_amatrix(&S);
      }
      coloffset += G->son[j * rsons]->cc->size;
    }

    /* compression */
    trunc_rkmatrix(tm, eps, R);

    sizenew = getsize_rkmatrix(R);

    /* use new rank-k-matrix or discard */
    if (sizenew < sizeold) {
      for (j = 0; j < csons; ++j) {
	for (i = 0; i < rsons; ++i) {
	  unref_hmatrix(G->son[i + j * rsons]);
	}
      }

      G->rsons = 0;
      G->csons = 0;
      freemem(G->son);
      G->son = NULL;
      G->f = NULL;
      G->r = R;
      G->desc = 1;
    }
    else {
      del_rkmatrix(R);
    }

  }
}
