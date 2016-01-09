/* ------------------------------------------------------------
 This is the file "curve2d.c" of the H2Lib package.
 All rights reserved, Steffen Boerm 2010
 ------------------------------------------------------------ */

#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "curve2d.h"
#include "basic.h"

/* ------------------------------------------------------------
 constructors / destructors
 ------------------------------------------------------------ */

pcurve2d
new_curve2d(uint vertices, uint edges)
{
  pcurve2d  gr;

  gr = (pcurve2d) allocmem(sizeof(curve2d));
  gr->x = (real(*)[2]) allocmem((size_t) sizeof(real[2]) * vertices);
  gr->e = (uint(*)[2]) allocmem((size_t) sizeof(uint[2]) * edges);
  gr->n = (real(*)[2]) allocmem((size_t) sizeof(real[2]) * edges);
  gr->g = (real *) allocmem((size_t) sizeof(real) * edges);

  gr->vertices = vertices;
  gr->edges = edges;

  return gr;
}

void
prepare_curve2d(pcurve2d gr)
{
  const     real(*x)[2] = (const real(*)[2]) gr->x;
  const     uint(*e)[2] = (const uint(*)[2]) gr->e;
  real(*n)[2] = gr->n;
  real     *g = gr->g;
  real      norm;
  uint      edges = gr->edges;
  uint      i;

  for (i = 0; i < edges; i++) {
    n[i][0] = x[e[i][1]][1] - x[e[i][0]][1];
    n[i][1] = x[e[i][0]][0] - x[e[i][1]][0];
    norm = REAL_SQRT(REAL_SQR(n[i][0]) + REAL_SQR(n[i][1]));
    g[i] = norm;
    n[i][0] /= norm;
    n[i][1] /= norm;
  }
}

void
del_curve2d(pcurve2d gr)
{
  freemem(gr->g);
  freemem(gr->n);
  freemem(gr->e);
  freemem(gr->x);
  freemem(gr);
}

/* ------------------------------------------------------------
 curve creation
 ------------------------------------------------------------ */

pcurve2d
new_circle_curve2d(uint edges, real r)
{
  real(*x)[2];
  uint(*e)[2];
  pcurve2d  gr;
  uint      i;

  assert(edges >= 3);

  gr = new_curve2d(edges, edges);
  x = gr->x;
  e = gr->e;
  for (i = 0; i < edges; i++) {
    x[i][0] = r * sin(2.0 * M_PI * i / edges);
    x[i][1] = r * cos(2.0 * M_PI * i / edges);
    e[i][0] = i;
    e[i][1] = (i + 1) % edges;
  }

  prepare_curve2d(gr);

  return gr;
}

pcurve2d
new_square_curve2d(uint edges, real a)
{
  real(*x)[2];
  uint(*e)[2];
  pcurve2d  gr;
  uint      i, top, left, bottom;

  assert(edges >= 4);

  a *= 0.5;

  top = edges / 4;
  left = edges / 2;
  bottom = 3 * edges / 4;

  gr = new_curve2d(edges, edges);
  x = gr->x;
  e = gr->e;
  for (i = 0; i < top; i++) {
    x[i][0] = a - 2.0 * a * i / top;
    x[i][1] = a;
  }
  for (; i < left; i++) {
    x[i][0] = -a;
    x[i][1] = a - 2.0 * a * (i - top) / (left - top);
  }
  for (; i < bottom; i++) {
    x[i][0] = -a + 2.0 * a * (i - left) / (bottom - left);
    x[i][1] = -a;
  }
  for (; i < edges; i++) {
    x[i][0] = a;
    x[i][1] = -a + 2.0 * a * (i - bottom) / (edges - bottom);
  }
  for (i = 0; i < edges; i++) {
    e[i][0] = i;
    e[i][1] = (i + 1) % edges;
  }

  prepare_curve2d(gr);

  return gr;
}

typedef enum _HORSESHOE {
  LEFT = 0, BOTTOM = 1, RIGHT = 2, TOP = 3
} HORSESHOE;

INLINE_PREFIX uint
hilbert_rec_curve2d(pcurve2d gr, uint n, real lengthx,
		    real lengthy, HORSESHOE s, uint i, real * a, real * b)
{
  real(*x)[2] = gr->x;
  uint(*e)[2] = gr->e;

  real      c[8], d[8];
  HORSESHOE cases[4] = {LEFT};
  uint      j = 0;
  real      phi;

  if (n == 1) {
    switch (s) {
    case LEFT:
      phi = 0;

      x[i][0] = (b[0] + 3.0 * a[0]) * 0.25;
      x[i][1] = (b[1] + 3.0 * a[1]) * 0.25;

      for (j = 1; j < 4; ++j) {
	x[i + j][0] = x[i + j - 1][0] + cos(phi) * lengthx;
	x[i + j][1] = x[i + j - 1][1] + sin(phi) * lengthy;

	e[i + j - 1][0] = i + j - 1;
	e[i + j - 1][1] = i + j;

	phi += M_PI * 0.5;
      }
      break;
    case BOTTOM:
      phi = M_PI * 0.5;

      x[i][0] = (b[0] + 3.0 * a[0]) * 0.25;
      x[i][1] = (b[1] + 3.0 * a[1]) * 0.25;

      for (j = 1; j < 4; ++j) {
	x[i + j][0] = x[i + j - 1][0] + cos(phi) * lengthx;
	x[i + j][1] = x[i + j - 1][1] + sin(phi) * lengthy;

	e[i + j - 1][0] = i + j - 1;
	e[i + j - 1][1] = i + j;

	phi -= M_PI * 0.5;
      }
      break;
    case RIGHT:
      phi = M_PI;

      x[i][0] = (3.0 * b[0] + a[0]) * 0.25;
      x[i][1] = (3.0 * b[1] + a[1]) * 0.25;

      for (j = 1; j < 4; ++j) {
	x[i + j][0] = x[i + j - 1][0] + cos(phi) * lengthx;
	x[i + j][1] = x[i + j - 1][1] + sin(phi) * lengthy;

	e[i + j - 1][0] = i + j - 1;
	e[i + j - 1][1] = i + j;

	phi += M_PI * 0.5;
      }
      break;
    case TOP:
      phi = 1.5 * M_PI;

      x[i][0] = (3.0 * b[0] + a[0]) * 0.25;
      x[i][1] = (3.0 * b[1] + a[1]) * 0.25;

      for (j = 1; j < 4; ++j) {
	x[i + j][0] = x[i + j - 1][0] + cos(phi) * lengthx;
	x[i + j][1] = x[i + j - 1][1] + sin(phi) * lengthy;

	e[i + j - 1][0] = i + j - 1;
	e[i + j - 1][1] = i + j;

	phi -= M_PI * 0.5;
      }
      break;
    }

    i = i + 4;

  }
  else {
    assert(n > 1);

    switch (s) {
    case LEFT:
      /* left bottom */
      c[0] = a[0];
      c[1] = a[1];
      d[0] = a[0] + (b[0] - a[0]) * 0.5;
      d[1] = a[1] + (b[1] - a[1]) * 0.5;
      /* right bottom */
      c[2] = a[0] + (b[0] - a[0]) * 0.5;
      c[3] = a[1];
      d[2] = b[0];
      d[3] = a[1] + (b[1] - a[1]) * 0.5;
      /* right top */
      c[4] = a[0] + (b[0] - a[0]) * 0.5;
      c[5] = a[1] + (b[1] - a[1]) * 0.5;
      d[4] = b[0];
      d[5] = b[1];
      /* left top */
      c[6] = a[0];
      c[7] = a[1] + (b[1] - a[1]) * 0.5;
      d[6] = a[0] + (b[0] - a[0]) * 0.5;
      d[7] = b[1];

      cases[0] = BOTTOM;
      cases[1] = LEFT;
      cases[2] = LEFT;
      cases[3] = TOP;
      break;
    case BOTTOM:
      /* left bottom */
      c[0] = a[0];
      c[1] = a[1];
      d[0] = a[0] + (b[0] - a[0]) * 0.5;
      d[1] = a[1] + (b[1] - a[1]) * 0.5;
      /* left top */
      c[2] = a[0];
      c[3] = a[1] + (b[1] - a[1]) * 0.5;
      d[2] = a[0] + (b[0] - a[0]) * 0.5;
      d[3] = b[1];
      /* right top */
      c[4] = a[0] + (b[0] - a[0]) * 0.5;
      c[5] = a[1] + (b[1] - a[1]) * 0.5;
      d[4] = b[0];
      d[5] = b[1];
      /* right bottom */
      c[6] = a[0] + (b[0] - a[0]) * 0.5;
      c[7] = a[1];
      d[6] = b[0];
      d[7] = a[1] + (b[1] - a[1]) * 0.5;

      cases[0] = LEFT;
      cases[1] = BOTTOM;
      cases[2] = BOTTOM;
      cases[3] = RIGHT;
      break;
    case RIGHT:
      /* right top */
      c[0] = a[0] + (b[0] - a[0]) * 0.5;
      c[1] = a[1] + (b[1] - a[1]) * 0.5;
      d[0] = b[0];
      d[1] = b[1];
      /* left top */
      c[2] = a[0];
      c[3] = a[1] + (b[1] - a[1]) * 0.5;
      d[2] = a[0] + (b[0] - a[0]) * 0.5;
      d[3] = b[1];
      /* left bottom */
      c[4] = a[0];
      c[5] = a[1];
      d[4] = a[0] + (b[0] - a[0]) * 0.5;
      d[5] = a[1] + (b[1] - a[1]) * 0.5;
      /* right bottom */
      c[6] = a[0] + (b[0] - a[0]) * 0.5;
      c[7] = a[1];
      d[6] = b[0];
      d[7] = a[1] + (b[1] - a[1]) * 0.5;

      cases[0] = TOP;
      cases[1] = RIGHT;
      cases[2] = RIGHT;
      cases[3] = BOTTOM;
      break;
    case TOP:
      /* right top */
      c[0] = a[0] + (b[0] - a[0]) * 0.5;
      c[1] = a[1] + (b[1] - a[1]) * 0.5;
      d[0] = b[0];
      d[1] = b[1];
      /* right bottom */
      c[2] = a[0] + (b[0] - a[0]) * 0.5;
      c[3] = a[1];
      d[2] = b[0];
      d[3] = a[1] + (b[1] - a[1]) * 0.5;
      /* left bottom */
      c[4] = a[0];
      c[5] = a[1];
      d[4] = a[0] + (b[0] - a[0]) * 0.5;
      d[5] = a[1] + (b[1] - a[1]) * 0.5;
      /* left top */
      c[6] = a[0];
      c[7] = a[1] + (b[1] - a[1]) * 0.5;
      d[6] = a[0] + (b[0] - a[0]) * 0.5;
      d[7] = b[1];

      cases[0] = RIGHT;
      cases[1] = TOP;
      cases[2] = TOP;
      cases[3] = LEFT;
      break;
    }

    for (j = 0; j < 4; ++j) {
      i = hilbert_rec_curve2d(gr, n - 1, lengthx, lengthy, cases[j], i,
			      &c[2 * j], &d[2 * j]);

      assert(i > 0);

      if (i < gr->edges) {
	e[i - 1][0] = i - 1;
	e[i - 1][1] = i;
      }
    }
  }

  return i;
}

pcurve2d
new_hilbert_curve2d(uint n, real l)
{
  pcurve2d  gr;

  uint      i, vertices, edges;
  real      lengthx, lengthy;
  real      a[2], b[2];

  vertices = pow(2, 2 * n);
  edges = vertices - 1;

  gr = new_curve2d(vertices, edges);

  a[0] = -l * 0.5;
  a[1] = -l * 0.5;
  b[0] = l * 0.5;
  b[1] = l * 0.5;

  lengthx = (b[0] - a[0]) / pow(2, n);
  lengthy = (b[1] - a[1]) / pow(2, n);

  i = hilbert_rec_curve2d(gr, n, lengthx, lengthy, LEFT, 0, a, b);

  assert(i == vertices);

  prepare_curve2d(gr);

  return gr;
}

pcurve2d
new_star_curve2d(uint edges, real r)
{
  real(*x)[2];
  uint(*e)[2];
  pcurve2d  gr;
  uint      i, j, d, nu, spikes;
  real      ri;

  spikes = 16;

  assert(edges >= spikes);
  assert(edges % spikes == 0);

  d = edges / spikes;
  assert(d * spikes == edges);

  ri = 0.25 * r;

  gr = new_curve2d(edges, edges);
  x = gr->x;
  e = gr->e;

  for (i = 0; i < spikes; i++) {
    x[d * i][0] = sin(2.0 * M_PI * i / spikes);
    x[d * i][1] = cos(2.0 * M_PI * i / spikes);
    if (i % 2 == 0) {
      x[d * i][0] *= r;
      x[d * i][1] *= r;
    }
    else {
      x[d * i][0] *= ri;
      x[d * i][1] *= ri;
    }
  }

  nu = 0;
  for (i = 0; i < spikes; ++i) {
    e[nu][0] = nu;
    e[nu][1] = (nu + 1) % edges;
    nu++;
    for (j = 1; j < d; ++j) {
      x[nu][0] = x[d * i][0]
	+ (x[d * ((i + 1) % spikes)][0] -
	   x[d * i][0]) * ((real) j / (real) d);
      x[nu][1] = x[d * i][1]
	+ (x[d * ((i + 1) % spikes)][1] -
	   x[d * i][1]) * ((real) j / (real) d);
      e[nu][0] = nu;
      e[nu][1] = (nu + 1) % edges;
      nu++;
    }
  }

  assert(nu == edges);

  prepare_curve2d(gr);

  return gr;
}

/* ------------------------------------------------------------
 output functions
 ------------------------------------------------------------ */

void
print_curve2d(pccurve2d gr)
{
  uint      vertices = gr->vertices;
  uint      edges = gr->edges;
  const     real(*x)[2] = (const real(*)[2]) gr->x;
  const     uint(*e)[2] = (const uint(*)[2]) gr->e;
  const     real(*n)[2] = (const real(*)[2]) gr->n;
  const real *g = (const real *) gr->g;
  uint      i;

  (void) printf("curve2d(%u,%u)\n", vertices, edges);
  for (i = 0; i < vertices; i++)
    (void) printf(" (% .5e % .5e)\n", x[i][0], x[i][1]);
  for (i = 0; i < edges; i++)
    (void) printf(" (%u %u   % .5e % .5e   %.5e)\n", e[i][0], e[i][1],
		  n[i][0], n[i][1], g[i]);
}

#ifdef USE_CAIRO
static void
compute_boundingbox_curve2d(pccurve2d gr, real * bmin, real * bmax)
{
  uint      i;
  uint      vertices = gr->vertices;
  real(*x)[2] = gr->x;

  bmin[0] = x[0][0];
  bmin[1] = x[0][1];
  bmax[0] = x[0][0];
  bmax[1] = x[0][1];

  for (i = 1; i < vertices; ++i) {
    if (x[i][0] < bmin[0]) {
      bmin[0] = x[i][0];
    }
    if (x[i][0] > bmax[0]) {
      bmax[0] = x[i][0];
    }

    if (x[i][1] < bmin[1]) {
      bmin[1] = x[i][1];
    }
    if (x[i][1] > bmax[1]) {
      bmax[1] = x[i][1];
    }
  }
}

void
draw_curve2d(pccurve2d gr, cairo_t * cr, real scale)
{
  uint      i, edges;
  real     *p, bmin[2], bmax[2];
  real      s, tx, ty, tol;

  edges = gr->edges;

  tol = 0.05;

  compute_boundingbox_curve2d(gr, bmin, bmax);

  if ((bmax[0] - bmin[0]) >= (bmax[1] - bmin[1])) {
    s = (1.0 - tol) / (bmax[0] - bmin[0]);
  }
  else {
    s = (1.0 - tol) / (bmax[1] - bmin[1]);
  }

  tx = -1.5 * bmin[0] + 0.5 * (bmax[0] + bmin[0]);
  ty = -1.5 * bmin[1] + 0.5 * (bmax[1] + bmin[1]);

  cairo_scale(cr, s, s);
  cairo_translate(cr, tx / s, ty / s);

  cairo_set_line_width(cr, cairo_get_line_width(cr) / scale);

  for (i = 0; i < edges; ++i) {
    p = gr->x[gr->e[i][0]];
    cairo_move_to(cr, p[0], p[1]);
    p = gr->x[gr->e[i][1]];
    cairo_line_to(cr, p[0], p[1]);
  }

  cairo_stroke(cr);
}
#endif
