
/* ------------------------------------------------------------
 * This is the file "uniform.c" of the H2Lib package.
 * All rights reserved, Steffen Boerm 2015
 * ------------------------------------------------------------ */

#include <math.h>

#include "duniform.h"
#include "basic.h"

/* ------------------------------------------------------------
 * Constructors and destructors
 * ------------------------------------------------------------ */

pduniform
new_duniform(pdclusterbasis rb, uint rd, pdclusterbasis cb, uint cd)
{
  pduniform u;

  u = allocmem(sizeof(duniform));

  u->rb = rb;
  u->cb = cb;
  u->rd = rd;
  u->cd = cd;

  init_amatrix(&u->S, rb->k[rd], cb->k[cd]);

  return u;
}

void
del_duniform(pduniform u)
{
  assert(u != 0);

  uninit_amatrix(&u->S);

  freemem(u);
}

/* ------------------------------------------------------------
 * Statistics
 * ------------------------------------------------------------ */

size_t
getsize_duniform(pcduniform u)
{
  size_t    sz;

  sz = sizeof(duniform);
  sz += getsize_heap_amatrix(&u->S);

  return sz;
}

/* ------------------------------------------------------------
 * Simple utility functions
 * ------------------------------------------------------------ */

void
clear_duniform(pduniform u)
{
  clear_amatrix(&u->S);
}

void
copy_duniform(pcduniform src, pduniform trg)
{
  assert(src->rb == trg->rb);
  assert(src->rd == trg->rd);
  assert(src->cb == trg->cb);
  assert(src->cd == trg->cd);

  copy_amatrix(false, &src->S, &trg->S);
}

/* ------------------------------------------------------------
 * Matrix-vector multiplication
 * ------------------------------------------------------------ */

void
fastaddeval_duniform_avector(field alpha, pcduniform u,
			     pcavector xt, pavector yt)
{
  avector   tmp1, tmp2;
  pavector  xc, yc;
  pcdclusterbasis rb, cb;
  uint      rd, cd;

  assert(xt->dim >= u->cb->koff[u->cb->directions]);
  assert(yt->dim >= u->rb->koff[u->rb->directions]);

  /* Get row and colum basis */
  rb = u->rb;
  cb = u->cb;

  /* Get row and column direction */
  rd = u->rd;
  cd = u->cd;

  assert(u->S.rows == rb->k[rd]);
  assert(u->S.cols == cb->k[cd]);

  /* This part of xt corresponds to the direction cd of the
   * column cluster */
  xc = init_sub_avector(&tmp1, (pavector) xt, cb->k[cd], cb->koff[cd]);

  /* This part of yt corresponds to the direction rd of the
   * row cluster */
  yc = init_sub_avector(&tmp2, yt, rb->k[rd], rb->koff[rd]);

  /* Multiply by coupling matrix */
  addeval_amatrix_avector(alpha, &u->S, xc, yc);

  /* Clean up source sub-vector */
  uninit_avector(yc);

  /* Clean up target sub-vector */
  uninit_avector(xc);
}

void
fastaddevaltrans_duniform_avector(field alpha, pcduniform u,
				  pcavector xt, pavector yt)
{
  avector   tmp1, tmp2;
  pavector  xc, yc;
  pcdclusterbasis rb, cb;
  uint      rd, cd;

  assert(xt->dim >= u->rb->koff[u->rb->directions]);
  assert(yt->dim >= u->cb->koff[u->cb->directions]);

  /* Get row and colum basis */
  rb = u->rb;
  cb = u->cb;

  /* Get row and column direction */
  rd = u->rd;
  cd = u->cd;

  assert(u->S.rows == rb->k[rd]);
  assert(u->S.cols == cb->k[cd]);

  /* This part of xt corresponds to the direction rd of the
   * row cluster */
  xc = init_sub_avector(&tmp1, (pavector) xt, rb->k[rd], rb->koff[rd]);

  /* This part of yt corresponds to the direction cd of the
   * column cluster */
  yc = init_sub_avector(&tmp2, yt, cb->k[cd], cb->koff[cd]);

  /* Multiply by adjoint coupling matrix */
  addevaltrans_amatrix_avector(alpha, &u->S, xc, yc);

  /* Clean up source sub-vector */
  uninit_avector(yc);

  /* Clean up target sub-vector */
  uninit_avector(xc);
}

void
slowaddeval_duniform_avector(field alpha, pcduniform u,
			     pcavector x, pavector y)
{
  avector   tmp1, tmp2;
  pavector  xp, yp, xt, yt;
  pcdclusterbasis rb = u->rb;
  pcdclusterbasis cb = u->cb;
  const uint *ridx = rb->t->idx;
  const uint *cidx = cb->t->idx;
  uint      i;

  xp = init_avector(&tmp1, cb->t->size);
  for (i = 0; i < cb->t->size; i++)
    xp->v[i] = x->v[cidx[i]];

  xt = init_avector(&tmp2, cb->kbranch);
  compress_dclusterbasis(cb, xp, xt);

  uninit_avector(xp);

  yt = init_avector(&tmp1, rb->kbranch);
  clear_avector(yt);

  fastaddeval_duniform_avector(alpha, u, xt, yt);

  uninit_avector(xt);

  yp = init_avector(&tmp2, rb->t->size);
  clear_avector(yp);

  expand_dclusterbasis(1.0, rb, yt, yp);

  for (i = 0; i < rb->t->size; i++)
    y->v[ridx[i]] += yp->v[i];

  uninit_avector(yp);
  uninit_avector(yt);
}

void
slowaddevaltrans_duniform_avector(field alpha, pcduniform u,
				  pcavector x, pavector y)
{
  avector   tmp1, tmp2;
  pavector  xp, yp, xt, yt;
  pcdclusterbasis rb = u->rb;
  pcdclusterbasis cb = u->cb;
  const uint *ridx = rb->t->idx;
  const uint *cidx = cb->t->idx;
  uint      i;

  xp = init_avector(&tmp1, rb->t->size);
  for (i = 0; i < rb->t->size; i++)
    xp->v[i] = x->v[ridx[i]];

  xt = init_avector(&tmp2, rb->kbranch);
  compress_dclusterbasis(rb, xp, xt);

  uninit_avector(xp);

  yt = init_avector(&tmp1, cb->kbranch);
  clear_avector(yt);

  fastaddevaltrans_duniform_avector(alpha, u, xt, yt);

  uninit_avector(xt);

  yp = init_avector(&tmp2, cb->t->size);
  clear_avector(yp);

  expand_dclusterbasis(1.0, cb, yt, yp);

  for (i = 0; i < cb->t->size; i++)
    y->v[cidx[i]] += yp->v[i];

  uninit_avector(yp);
  uninit_avector(yt);
}

/* ------------------------------------------------------------
 * Conversion to a full matrix
 * ------------------------------------------------------------ */

static void
expand(pcdclusterbasis rb, uint rd,
       pcdclusterbasis cb, uint cd, pcamatrix S, bool permute, pamatrix G)
{
  amatrix   tmp1, tmp2, tmp3;
  pamatrix  Sc, S1, Gc, Gp, G1;
  longindex ldG, ldGp;
  const uint *ridx;
  const uint *cidx;
  uint      rd1, cd1;
  uint      roff, coff;
  uint      i, j, ii, jj;

  assert(S->rows == rb->k[rd]);
  assert(S->cols == cb->k[cd]);

  if (rb->sons > 0) {
    if (cb->sons > 0) {
      coff = 0;
      for (j = 0; j < cb->sons; j++) {
	cd1 = cb->dirson[j][cd];

	Sc = init_amatrix(&tmp1, rb->k[rd], cb->son[j]->k[cd1]);

	clear_amatrix(Sc);
	addmul_amatrix(1.0, false, S, true, cb->E[j] + cd, Sc);

	roff = 0;
	for (i = 0; i < rb->sons; i++) {
	  rd1 = rb->dirson[i][rd];

	  S1 = init_amatrix(&tmp2, rb->son[i]->k[rd1], cb->son[j]->k[cd1]);

	  clear_amatrix(S1);
	  addmul_amatrix(1.0, false, rb->E[i] + rd, false, Sc, S1);

	  if (permute)
	    expand(rb->son[i], rd1, cb->son[j], cd1, S1, true, G);
	  else {
	    G1 = init_sub_amatrix(&tmp3, G,
				  rb->son[i]->t->size, roff,
				  cb->son[j]->t->size, coff);
	    expand(rb->son[i], rd1, cb->son[j], cd1, S1, false, G1);
	    uninit_amatrix(G1);
	  }

	  uninit_amatrix(S1);

	  roff += rb->son[i]->t->size;
	}
	assert(roff == rb->t->size);

	uninit_amatrix(Sc);

	coff += cb->son[j]->t->size;
      }
      assert(coff == cb->t->size);
    }
    else {
      roff = 0;
      for (i = 0; i < rb->sons; i++) {
	rd1 = rb->dirson[i][rd];

	S1 = init_amatrix(&tmp2, rb->son[i]->k[rd1], cb->k[cd]);

	clear_amatrix(S1);
	addmul_amatrix(1.0, false, rb->E[i] + rd, false, S, S1);

	if (permute)
	  expand(rb->son[i], rd1, cb, cd, S1, true, G);
	else {
	  G1 = init_sub_amatrix(&tmp3, G,
				rb->son[i]->t->size, roff, cb->t->size, 0);
	  expand(rb->son[i], rd1, cb, cd, S1, false, G1);
	  uninit_amatrix(G1);
	}

	uninit_amatrix(S1);

	roff += rb->son[i]->t->size;
      }
      assert(roff == rb->t->size);
    }
  }
  else {
    if (cb->sons > 0) {
      coff = 0;
      for (j = 0; j < cb->sons; j++) {
	cd1 = cb->dirson[j][cd];

	S1 = init_amatrix(&tmp2, rb->k[rd], cb->son[j]->k[cd1]);

	clear_amatrix(S1);
	addmul_amatrix(1.0, false, S, true, cb->E[j] + cd, S1);

	if (permute)
	  expand(rb, rd, cb->son[j], cd1, S1, true, G);
	else {
	  G1 = init_sub_amatrix(&tmp3, G,
				rb->t->size, 0, cb->son[j]->t->size, coff);
	  expand(rb, rd, cb->son[j], cd1, S1, false, G1);
	  uninit_amatrix(G1);
	}

	uninit_amatrix(S1);

	coff += cb->son[j]->t->size;
      }
      assert(coff == cb->t->size);
    }
    else {
      if (permute) {
	Gc = init_amatrix(&tmp1, rb->t->size, cb->k[cd]);

	clear_amatrix(Gc);
	addmul_amatrix(1.0, false, rb->V + rd, false, S, Gc);

	Gp = init_amatrix(&tmp2, rb->t->size, cb->t->size);

	clear_amatrix(Gp);
	addmul_amatrix(1.0, false, Gc, true, cb->V + cd, Gp);

	ridx = rb->t->idx;
	cidx = cb->t->idx;
	ldG = G->ld;
	ldGp = Gp->ld;
	for (j = 0; j < Gp->cols; j++) {
	  jj = cidx[j];
	  for (i = 0; i < Gp->rows; i++) {
	    ii = ridx[i];
	    G->a[ii + jj * ldG] += Gp->a[i + j * ldGp];
	  }
	}

	uninit_amatrix(Gp);
	uninit_amatrix(Gc);
      }
      else {
	Gc = init_amatrix(&tmp1, rb->t->size, cb->k[cd]);

	clear_amatrix(Gc);
	addmul_amatrix(1.0, false, rb->V + rd, false, S, Gc);

	addmul_amatrix(1.0, false, Gc, true, cb->V + cd, G);

	uninit_amatrix(Gc);
      }
    }
  }
}

void
expand_duniform(field alpha, pcduniform u, pamatrix G)
{
  amatrix   tmp1;
  pamatrix  S;

  if (alpha == 1.0)
    expand(u->rb, u->rd, u->cb, u->cd, &u->S, true, G);
  else {
    S = init_amatrix(&tmp1, u->rb->k[u->rd], u->cb->k[u->cd]);
    clear_amatrix(S);
    add_amatrix(alpha, false, &u->S, S);

    expand(u->rb, u->rd, u->cb, u->cd, S, true, G);

    uninit_amatrix(S);
  }
}


real
norm2_fast_duniform(pcduniform u, pcdclusteroperator ro,
		    pcdclusteroperator co)
{

  real      norm;
  amatrix   tmp1, tmp2;
  pamatrix  A, B;
  uint      rd, cd;

  rd = u->rd;
  cd = u->cd;

  if (ro) {
    if (co) {

      A = init_amatrix(&tmp1, ro->C[rd].rows, u->S.cols);
      B = init_amatrix(&tmp2, ro->C[rd].rows, co->C[cd].rows);
      clear_amatrix(A);
      clear_amatrix(B);

      addmul_amatrix(1.0, false, &ro->C[rd], false, &u->S, A);
      addmul_amatrix(1.0, false, A, true, &co->C[cd], B);

      uninit_amatrix(A);
    }
    else {
      B = init_amatrix(&tmp2, ro->C[rd].rows, u->S.cols);
      clear_amatrix(B);

      addmul_amatrix(1.0, false, &ro->C[rd], false, &u->S, B);
    }
  }
  else {
    if (co) {

      B = init_amatrix(&tmp2, u->rb->k[rd], co->C[cd].rows);
      clear_amatrix(B);

      addmul_amatrix(1.0, false, &u->S, true, &co->C[cd], B);
    }
    else {
      B = (pamatrix) &u->S;
    }
  }

  norm = norm2_amatrix(B);

  if (B != &u->S) {
    uninit_amatrix(B);
  }

  return norm;
}

real
normfrob_fast_duniform(pcduniform u, pcdclusteroperator ro,
		       pcdclusteroperator co)
{

  real      norm;
  amatrix   tmp1, tmp2;
  pamatrix  A, B;
  uint      cd, rd;

  cd = u->cd;
  rd = u->rd;

  if (ro) {
    if (co) {

      A = init_amatrix(&tmp1, ro->C[rd].rows, u->S.cols);
      B = init_amatrix(&tmp2, ro->C[rd].rows, co->C[cd].rows);
      clear_amatrix(A);
      clear_amatrix(B);

      addmul_amatrix(1.0, false, &ro->C[rd], false, &u->S, A);
      addmul_amatrix(1.0, false, A, true, &co->C[cd], B);

      uninit_amatrix(A);
    }
    else {
      B = init_amatrix(&tmp2, ro->C[rd].rows, u->S.cols);
      clear_amatrix(B);

      addmul_amatrix(1.0, false, &ro->C[rd], false, &u->S, B);
    }
  }
  else {
    if (co) {

      B = init_amatrix(&tmp2, u->rb->k[rd], co->C[cd].rows);
      clear_amatrix(B);

      addmul_amatrix(1.0, false, &u->S, true, &co->C[cd], B);
    }
    else {
      B = (pamatrix) &u->S;
    }
  }

  norm = normfrob_amatrix(B);

  if (B != &u->S) {
    uninit_amatrix(B);
  }

  return norm;
}


void
add_projected_duniform(pcduniform u, pcdclusteroperator ro,
		       pcdclusteroperator co, pduniform unew)
{

  amatrix   tmp1;
  pamatrix  A;

  if (u->rb == unew->rb) {
    if (u->cb == unew->cb) {	/* Nothing changed */
      add_amatrix(1.0, false, &u->S, &unew->S);
    }
    else {			/* Different column cluster basis */
      assert(co->krow[u->cd] == unew->cb->k[u->cd]);
      assert(co->kcol[u->cd] == u->cb->k[u->cd]);

      addmul_amatrix(1.0, false, &u->S, true, &co->C[u->cd], &unew->S);
    }
  }
  else {
    if (u->cb == unew->cb) {	/* Different row cluster basis */
      assert(ro->krow[u->rd] == unew->rb->k[u->rd]);
      assert(ro->kcol[u->rd] == u->rb->k[u->rd]);

      addmul_amatrix(1.0, false, &ro->C[u->rd], false, &u->S, &unew->S);
    }
    else {			/* Both cluster basis are different */
      assert(ro->krow[u->rd] == unew->rb->k[u->rd]);
      assert(ro->kcol[u->rd] == u->rb->k[u->rd]);
      assert(co->krow[u->cd] == unew->cb->k[u->cd]);
      assert(co->kcol[u->cd] == u->cb->k[u->cd]);

      A = init_amatrix(&tmp1, unew->rb->k[u->rd], u->cb->k[u->cd]);
      clear_amatrix(A);
      addmul_amatrix(1.0, false, &ro->C[u->rd], false, &u->S, A);
      addmul_amatrix(1.0, false, A, true, &co->C[u->cd], &unew->S);
      uninit_amatrix(A);
    }
  }
}
