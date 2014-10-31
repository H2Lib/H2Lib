/* ------------------------------------------------------------
 This is the file "singquad1d.c" of the H2Lib package.
 All rights reserved, Sven Christophersen 2010
 ------------------------------------------------------------ */

#include "singquad1d.h"

#include "basic.h"

INLINE_PREFIX void
build_single_singquad1d(psingquad1d sq, uint q, preal x, preal w)
{

  uint      i;

  sq->x_single = allocreal(q);
  sq->w_single = allocreal(q);
  sq->base_single = 0.0;
  sq->n_single = q;

  for (i = 0; i < q; ++i) {
    sq->x_single[i] = 0.5 + 0.5 * x[i];
    sq->w_single[i] = 0.5 * w[i];
  }
}

psingquad1d
build_log_singquad1d(uint q, preal x, preal w)
{

  psingquad1d sq;

  uint      i, j, k, nq, nq2;

  /*
   * Allocate memory
   */

  sq = (psingquad1d) allocmem(sizeof(singquad1d));
  nq = q * q;
  nq2 = 2 * nq;

  sq->x_id = allocreal(nq2);
  sq->y_id = allocreal(nq2);
  sq->w_id = allocreal(nq2);
  sq->n_id = nq2;

  sq->x_vert = allocreal(nq2);
  sq->y_vert = allocreal(nq2);
  sq->w_vert = allocreal(nq2);
  sq->n_vert = nq2;

  sq->x_dist = allocreal(nq);
  sq->y_dist = allocreal(nq);
  sq->w_dist = allocreal(nq);
  sq->n_dist = nq;

  /*
   * Compute needed quadrature points and weights for different cases.
   */

  /*
   * Identical case
   */

  k = 0;
  sq->base_id = 4.0 * REAL_LOG(0.5) / 3.0;

  /*
   * [0.5,0.75]x[0.0,0.25]
   */

  for (i = 0; i < q; ++i) {
    for (j = 0; j < q; ++j) {
      sq->x_id[k] = 0.125 * x[i] + 0.625;
      sq->y_id[k] = 0.125 * x[j] + 0.125;
      sq->w_id[k] = w[i] * w[j] / 6.0;
      k++;
    }
  }

  /*
   * [0.75,1.0]x[0.0,0.25]
   */

  for (i = 0; i < q; ++i) {
    for (j = 0; j < q; ++j) {
      sq->x_id[k] = 0.125 * x[i] + 0.875;
      sq->y_id[k] = 0.125 * x[j] + 0.125;
      sq->w_id[k] = w[i] * w[j] / 12.0;
      k++;
    }
  }

  /*
   * Vertex case
   */

  k = 0;
  sq->base_vert = REAL_LOG(0.5) / 3.0;

  /*
   * [-1.0,-0.5]x[0.0,0.5]
   */

  for (i = 0; i < q; ++i) {
    for (j = 0; j < q; ++j) {
      sq->x_vert[k] = 0.25 * x[i] - 0.75;
      sq->y_vert[k] = 0.25 * x[j] + 0.25;
      sq->w_vert[k] = w[i] * w[j] / 6.0;
      k++;
    }
  }

  /*
   * [-1.0,-0.5]x[0.5,1.0]
   */

  for (i = 0; i < q; ++i) {
    for (j = 0; j < q; ++j) {
      sq->x_vert[k] = 0.25 * x[i] - 0.75;
      sq->y_vert[k] = 0.25 * x[j] + 0.75;
      sq->w_vert[k] = w[i] * w[j] / 12.0;
      k++;
    }
  }

  /*
   * Distant case
   */

  k = 0;
  sq->base_dist = 0.0;

  /*
   * [0.0,1.0]x[0.0,1.0]
   */

  for (i = 0; i < q; ++i) {
    for (j = 0; j < q; ++j) {
      sq->x_dist[k] = 0.5 * x[i] + 0.5;
      sq->y_dist[k] = 0.5 * x[j] + 0.5;
      sq->w_dist[k] = w[i] * w[j] / 4.0;
      k++;
    }
  }

  /*
   * Quadrature for single boundary element
   */

  build_single_singquad1d(sq, q, x, w);

  return sq;
}

psingquad1d
build_pow_singquad1d(uint q, preal x, preal w, real alpha)
{

  psingquad1d sq;

  uint      i, j, k, nq, nq2;
  real      c;

  /*
   * Allocate memory
   */

  sq = (psingquad1d) allocmem(sizeof(singquad1d));
  nq = q * q;
  nq2 = 2 * nq;
  c = 0.25 * pow(0.5, alpha);

  sq->x_id = NULL;
  sq->y_id = NULL;
  sq->w_id = NULL;
  sq->n_id = 0;

  sq->x_vert = allocreal(nq2 + nq);
  sq->y_vert = allocreal(nq2 + nq);
  sq->w_vert = allocreal(nq2 + nq);
  sq->n_vert = nq2 + nq;

  sq->x_dist = allocreal(nq);
  sq->y_dist = allocreal(nq);
  sq->w_dist = allocreal(nq);
  sq->n_dist = nq;

  /*
   * Compute needed quadrature points and weights for different cases.
   */

  /*
   * Identical case
   */

  sq->base_id = 0.0;

  /*
   * Vertex case
   */

  k = 0;
  sq->base_vert = 0.0;

  /*
   * [-1.0,-0.5]x[0.0,0.5]
   */

  for (i = 0; i < q; ++i) {
    for (j = 0; j < q; ++j) {
      sq->x_vert[k] = 0.25 * x[i] - 0.75;
      sq->y_vert[k] = 0.25 * x[j] + 0.25;
      sq->w_vert[k] = w[i] * w[j] / (1.0 - c) / 16.0;
      k++;
    }
  }

  /*
   * [-0.5,0.0]x[0.5,1.0]
   */

  for (i = 0; i < q; ++i) {
    for (j = 0; j < q; ++j) {
      sq->x_vert[k] = 0.25 * x[i] - 0.25;
      sq->y_vert[k] = 0.25 * x[j] + 0.75;
      sq->w_vert[k] = w[i] * w[j] / (1.0 - c) / 16.0;
      k++;
    }
  }

  /*
   * [-1.0,-0.5]x[0.5,1.0]
   */

  for (i = 0; i < q; ++i) {
    for (j = 0; j < q; ++j) {
      sq->x_vert[k] = 0.25 * x[i] - 0.75;
      sq->y_vert[k] = 0.25 * x[j] + 0.75;
      sq->w_vert[k] = w[i] * w[j] / (1.0 - c) / 16.0;
      k++;
    }
  }

  /*
   * Distant case
   */
  k = 0;
  sq->base_dist = 0.0;

  /*
   * [0.0,1.0]x[0.0,1.0]
   */

  for (i = 0; i < q; ++i) {
    for (j = 0; j < q; ++j) {
      sq->x_dist[k] = 0.5 * x[i] + 0.5;
      sq->y_dist[k] = 0.5 * x[j] + 0.5;
      sq->w_dist[k] = w[i] * w[j] / 4.0;
      k++;
    }
  }

  /*
   * Quadrature for single boundary element
   */

  build_single_singquad1d(sq, q, x, w);

  return sq;
}

void
del_singquad1d(psingquad1d sq)
{
  if (sq->x_id != NULL)
    freemem(sq->x_id);
  if (sq->y_id != NULL)
    freemem(sq->y_id);
  if (sq->w_id != NULL)
    freemem(sq->w_id);

  if (sq->x_vert != NULL)
    freemem(sq->x_vert);
  if (sq->y_vert != NULL)
    freemem(sq->y_vert);
  if (sq->w_vert != NULL)
    freemem(sq->w_vert);

  if (sq->x_dist != NULL)
    freemem(sq->x_dist);
  if (sq->y_dist != NULL)
    freemem(sq->y_dist);
  if (sq->w_dist != NULL)
    freemem(sq->w_dist);

  if (sq->x_single != NULL)
    freemem(sq->x_single);
  if (sq->w_single != NULL)
    freemem(sq->w_single);

  freemem(sq);
}

uint
select_quadrature_singquad1d(pcsingquad1d sq, const uint * tv,
			     const uint * sv, uint * tp, uint * sp, real ** x,
			     real ** y, real ** w, uint * n, real * base)
{

  uint      p;

  p =
    (tv[0] == sv[0]) + (tv[0] == sv[1]) + (tv[1] == sv[0]) + (tv[1] == sv[1]);

  switch (p) {
  case 0:			/* distant */
    tp[0] = 0, tp[1] = 1;
    sp[0] = 0, sp[1] = 1;

    *x = sq->x_dist;
    *y = sq->y_dist;
    *w = sq->w_dist;
    *n = sq->n_dist;
    *base = sq->base_dist;
    break;
  case 1:			/* vertex */
    if (tv[0] == sv[0]) {
      tp[0] = 0, tp[1] = 1;
      sp[0] = 0, sp[1] = 1;
    }
    else if (tv[0] == sv[1]) {
      tp[0] = 0, tp[1] = 1;
      sp[0] = 1, sp[1] = 0;
    }
    else if (tv[1] == sv[0]) {
      tp[0] = 1, tp[1] = 0;
      sp[0] = 0, sp[1] = 1;
    }
    else if (tv[1] == sv[1]) {
      tp[0] = 1, tp[1] = 0;
      sp[0] = 1, sp[1] = 0;
    }

    *x = sq->x_vert;
    *y = sq->y_vert;
    *w = sq->w_vert;
    *n = sq->n_vert;
    *base = sq->base_vert;
    break;
  case 2:			/* identical */
    tp[0] = 0, tp[1] = 1;
    sp[0] = 0, sp[1] = 1;

    *x = sq->x_id;
    *y = sq->y_id;
    *w = sq->w_id;
    *n = sq->n_id;
    *base = sq->base_id;
    break;
  default:
    printf("ERROR: Unknown quadrature case!\n");
    exit(0);
    break;
  }

  return p;
}
