#include "tri2d.h"

#include "basic.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>

#ifdef USE_CAIRO
#include <cairo.h>
#include <cairo-pdf.h>
#endif

#define BUFSIZE 80

ptri2d
new_tri2d(uint vertices, uint edges, uint triangles)
{
  ptri2d    t2;

  t2 = (ptri2d) allocmem(sizeof(tri2d));

  t2->x = (real(*)[2]) allocmem(sizeof(real[2]) * vertices);
  t2->e = (uint(*)[2]) allocmem(sizeof(uint[2]) * edges);
  t2->t = (uint(*)[3]) allocmem(sizeof(uint[3]) * triangles);

  t2->xb = (uint *) allocmem(sizeof(uint) * vertices);
  t2->eb = (uint *) allocmem(sizeof(uint) * edges);

  t2->vertices = vertices;
  t2->edges = edges;
  t2->triangles = triangles;

  return t2;
}

void
del_tri2d(ptri2d t2)
{
  freemem(t2->eb);
  freemem(t2->xb);
  freemem(t2->t);
  freemem(t2->e);
  freemem(t2->x);
  freemem(t2);
}

ptri2d
new_unitsquare_tri2d()
{
  ptri2d    t2;

  t2 = new_tri2d(4, 5, 2);

  t2->x[0][0] = -1.0;
  t2->x[0][1] = -1.0;
  t2->x[1][0] = 1.0;
  t2->x[1][1] = -1.0;
  t2->x[2][0] = -1.0;
  t2->x[2][1] = 1.0;
  t2->x[3][0] = 1.0;
  t2->x[3][1] = 1.0;

  t2->xb[0] = t2->xb[1] = t2->xb[2] = t2->xb[3] = 1;

  t2->e[0][0] = 0;
  t2->e[0][1] = 1;
  t2->e[1][0] = 0;
  t2->e[1][1] = 2;
  t2->e[2][0] = 1;
  t2->e[2][1] = 2;
  t2->e[3][0] = 1;
  t2->e[3][1] = 3;
  t2->e[4][0] = 2;
  t2->e[4][1] = 3;

  t2->eb[0] = t2->eb[1] = t2->eb[3] = t2->eb[4] = 1;
  t2->eb[2] = 0;

  t2->t[0][0] = 0;
  t2->t[0][1] = 1;
  t2->t[0][2] = 2;
  t2->t[1][0] = 2;
  t2->t[1][1] = 3;
  t2->t[1][2] = 4;

  t2->vertices = 4;
  t2->edges = 5;
  t2->triangles = 2;

  return t2;
}

ptri2d
new_unitcircle_tri2d()
{
  ptri2d    t2;

  t2 = new_tri2d(5, 8, 4);

  t2->x[0][0] = 0.0;
  t2->x[0][1] = 0.0;
  t2->x[1][0] = -1.0;
  t2->x[1][1] = 0.0;
  t2->x[2][0] = 0.0;
  t2->x[2][1] = 1.0;
  t2->x[3][0] = 1.0;
  t2->x[3][1] = 0.0;
  t2->x[4][0] = 0.0;
  t2->x[4][1] = -1.0;

  t2->xb[0] = 0;
  t2->xb[1] = t2->xb[2] = t2->xb[3] = t2->xb[4] = 1;

  t2->e[0][0] = 1;
  t2->e[0][1] = 2;
  t2->e[1][0] = 2;
  t2->e[1][1] = 3;
  t2->e[2][0] = 3;
  t2->e[2][1] = 4;
  t2->e[3][0] = 4;
  t2->e[3][1] = 1;
  t2->e[4][0] = 1;
  t2->e[4][1] = 0;
  t2->e[5][0] = 2;
  t2->e[5][1] = 0;
  t2->e[6][0] = 3;
  t2->e[6][1] = 0;
  t2->e[7][0] = 4;
  t2->e[7][1] = 0;

  t2->eb[0] = t2->eb[1] = t2->eb[2] = t2->eb[3] = 1;
  t2->eb[4] = t2->eb[5] = t2->eb[6] = t2->eb[7] = 0;

  t2->t[0][0] = 0;
  t2->t[0][1] = 5;
  t2->t[0][2] = 4;
  t2->t[1][0] = 5;
  t2->t[1][1] = 1;
  t2->t[1][2] = 6;
  t2->t[2][0] = 7;
  t2->t[2][1] = 6;
  t2->t[2][2] = 2;
  t2->t[3][0] = 3;
  t2->t[3][1] = 4;
  t2->t[3][2] = 7;

  t2->vertices = 5;
  t2->edges = 8;
  t2->triangles = 4;

  return t2;
}

ptri2d
new_lshape_tri2d()
{
  ptri2d    t2;

  t2 = new_tri2d(8, 13, 6);

  t2->x[0][0] = -1.0;
  t2->x[0][1] = -1.0;
  t2->x[1][0] = 0.0;
  t2->x[1][1] = -1.0;
  t2->x[2][0] = 1.0;
  t2->x[2][1] = -1.0;
  t2->x[3][0] = -1.0;
  t2->x[3][1] = 0.0;
  t2->x[4][0] = 0.0;
  t2->x[4][1] = 0.0;
  t2->x[5][0] = 1.0;
  t2->x[5][1] = 0.0;
  t2->x[6][0] = -1.0;
  t2->x[6][1] = 1.0;
  t2->x[7][0] = 0.0;
  t2->x[7][1] = 1.0;

  t2->xb[0] = t2->xb[1] = t2->xb[2] = t2->xb[3] = t2->xb[4] = t2->xb[5] =
    t2->xb[6] = t2->xb[7] = 1;

  t2->e[0][0] = 0;
  t2->e[0][1] = 1;
  t2->e[1][0] = 1;
  t2->e[1][1] = 2;
  t2->e[2][0] = 0;
  t2->e[2][1] = 3;
  t2->e[3][0] = 1;
  t2->e[3][1] = 3;
  t2->e[4][0] = 1;
  t2->e[4][1] = 4;
  t2->e[5][0] = 2;
  t2->e[5][1] = 4;
  t2->e[6][0] = 2;
  t2->e[6][1] = 5;
  t2->e[7][0] = 3;
  t2->e[7][1] = 4;
  t2->e[8][0] = 4;
  t2->e[8][1] = 5;
  t2->e[9][0] = 3;
  t2->e[9][1] = 6;
  t2->e[10][0] = 4;
  t2->e[10][1] = 6;
  t2->e[11][0] = 4;
  t2->e[11][1] = 7;
  t2->e[12][0] = 6;
  t2->e[12][1] = 7;

  t2->eb[0] = t2->eb[1] = t2->eb[2] = t2->eb[6] = t2->eb[8] = t2->eb[9] =
    t2->eb[11] = t2->eb[12] = 1;
  t2->eb[3] = t2->eb[4] = t2->eb[5] = t2->eb[7] = t2->eb[10] = 0;

  t2->t[0][0] = 3;
  t2->t[0][1] = 0;
  t2->t[0][2] = 2;
  t2->t[1][0] = 3;
  t2->t[1][1] = 4;
  t2->t[1][2] = 7;
  t2->t[2][0] = 5;
  t2->t[2][1] = 4;
  t2->t[2][2] = 1;
  t2->t[3][0] = 5;
  t2->t[3][1] = 8;
  t2->t[3][2] = 6;
  t2->t[4][0] = 10;
  t2->t[4][1] = 9;
  t2->t[4][2] = 7;
  t2->t[5][0] = 10;
  t2->t[5][1] = 12;
  t2->t[5][2] = 11;

  t2->vertices = 8;
  t2->edges = 13;
  t2->triangles = 6;

  return t2;
}

ptri2d
new_ushape_tri2d()
{
  ptri2d    t2;
  uint      i;

  t2 = new_tri2d(23, 46, 24);

  t2->x[0][0] = -1.0;
  t2->x[0][1] = -1.0;
  t2->x[1][0] = -0.5;
  t2->x[1][1] = -1.0;
  t2->x[2][0] = 0.0;
  t2->x[2][1] = -1.0;
  t2->x[3][0] = 0.5;
  t2->x[3][1] = -1.0;
  t2->x[4][0] = 1.0;
  t2->x[4][1] = -1.0;
  t2->x[5][0] = -1.0;
  t2->x[5][1] = -0.5;
  t2->x[6][0] = -0.5;
  t2->x[6][1] = -0.5;
  t2->x[7][0] = 0.0;
  t2->x[7][1] = -0.5;
  t2->x[8][0] = 0.5;
  t2->x[8][1] = -0.5;
  t2->x[9][0] = 1.0;
  t2->x[9][1] = -0.5;
  t2->x[10][0] = -1.0;
  t2->x[10][1] = 0.0;
  t2->x[11][0] = -0.5;
  t2->x[11][1] = 0.0;
  t2->x[12][0] = 0.0;
  t2->x[12][1] = 0.0;
  t2->x[13][0] = 0.5;
  t2->x[13][1] = 0.0;
  t2->x[14][0] = 1.0;
  t2->x[14][1] = 0.0;
  t2->x[15][0] = -1.0;
  t2->x[15][1] = 0.5;
  t2->x[16][0] = -0.5;
  t2->x[16][1] = 0.5;
  t2->x[17][0] = 0.5;
  t2->x[17][1] = 0.5;
  t2->x[18][0] = 1.0;
  t2->x[18][1] = 0.5;
  t2->x[19][0] = -1.0;
  t2->x[19][1] = 1.0;
  t2->x[20][0] = -0.5;
  t2->x[20][1] = 1.0;
  t2->x[21][0] = 0.5;
  t2->x[21][1] = 1.0;
  t2->x[22][0] = 1.0;
  t2->x[22][1] = 1.0;

  t2->xb[0] = t2->xb[1] = t2->xb[2] = t2->xb[3] = t2->xb[4] = t2->xb[5] = 1;
  t2->xb[9] = t2->xb[10] = t2->xb[11] = t2->xb[12] = t2->xb[13] = t2->xb[14] =
    1;
  t2->xb[15] = t2->xb[16] = t2->xb[17] = t2->xb[18] = 1;
  t2->xb[19] = t2->xb[20] = t2->xb[21] = t2->xb[22] = 1;
  t2->xb[6] = t2->xb[7] = t2->xb[8] = 0;

  t2->e[0][0] = 0;
  t2->e[0][1] = 1;
  t2->e[1][0] = 1;
  t2->e[1][1] = 2;
  t2->e[2][0] = 2;
  t2->e[2][1] = 3;
  t2->e[3][0] = 3;
  t2->e[3][1] = 4;
  t2->e[4][0] = 5;
  t2->e[4][1] = 0;
  t2->e[5][0] = 5;
  t2->e[5][1] = 1;
  t2->e[6][0] = 6;
  t2->e[6][1] = 1;
  t2->e[7][0] = 6;
  t2->e[7][1] = 2;
  t2->e[8][0] = 7;
  t2->e[8][1] = 2;
  t2->e[9][0] = 7;
  t2->e[9][1] = 3;
  t2->e[10][0] = 8;
  t2->e[10][1] = 3;
  t2->e[11][0] = 8;
  t2->e[11][1] = 4;
  t2->e[12][0] = 9;
  t2->e[12][1] = 4;
  t2->e[13][0] = 5;
  t2->e[13][1] = 6;
  t2->e[14][0] = 6;
  t2->e[14][1] = 7;
  t2->e[15][0] = 7;
  t2->e[15][1] = 8;
  t2->e[16][0] = 8;
  t2->e[16][1] = 9;
  t2->e[17][0] = 10;
  t2->e[17][1] = 5;
  t2->e[18][0] = 10;
  t2->e[18][1] = 6;
  t2->e[19][0] = 11;
  t2->e[19][1] = 6;
  t2->e[20][0] = 11;
  t2->e[20][1] = 7;
  t2->e[21][0] = 12;
  t2->e[21][1] = 7;
  t2->e[22][0] = 12;
  t2->e[22][1] = 8;
  t2->e[23][0] = 13;
  t2->e[23][1] = 8;
  t2->e[24][0] = 13;
  t2->e[24][1] = 9;
  t2->e[25][0] = 14;
  t2->e[25][1] = 9;
  t2->e[26][0] = 10;
  t2->e[26][1] = 11;
  t2->e[27][0] = 11;
  t2->e[27][1] = 12;
  t2->e[28][0] = 12;
  t2->e[28][1] = 13;
  t2->e[29][0] = 13;
  t2->e[29][1] = 14;
  t2->e[30][0] = 15;
  t2->e[30][1] = 10;
  t2->e[31][0] = 15;
  t2->e[31][1] = 11;
  t2->e[32][0] = 16;
  t2->e[32][1] = 11;
  t2->e[33][0] = 17;
  t2->e[33][1] = 13;
  t2->e[34][0] = 17;
  t2->e[34][1] = 14;
  t2->e[35][0] = 18;
  t2->e[35][1] = 14;
  t2->e[36][0] = 15;
  t2->e[36][1] = 16;
  t2->e[37][0] = 17;
  t2->e[37][1] = 18;
  t2->e[38][0] = 19;
  t2->e[38][1] = 15;
  t2->e[39][0] = 19;
  t2->e[39][1] = 16;
  t2->e[40][0] = 20;
  t2->e[40][1] = 16;
  t2->e[41][0] = 21;
  t2->e[41][1] = 17;
  t2->e[42][0] = 21;
  t2->e[42][1] = 18;
  t2->e[43][0] = 22;
  t2->e[43][1] = 18;
  t2->e[44][0] = 19;
  t2->e[44][1] = 20;
  t2->e[45][0] = 21;
  t2->e[45][1] = 22;

  for (i = 0; i < 46; i++)
    t2->eb[i] = 0;

  t2->eb[0] = t2->eb[1] = t2->eb[2] = t2->eb[3] = t2->eb[4] = 1;
  t2->eb[12] = t2->eb[17] = t2->eb[25] = t2->eb[27] = t2->eb[28] = 1;
  t2->eb[30] = t2->eb[32] = t2->eb[33] = t2->eb[35] = 1;
  t2->eb[38] = t2->eb[40] = t2->eb[41] = t2->eb[43] = t2->eb[44] =
    t2->eb[45] = 1;

  t2->t[0][0] = 0;
  t2->t[0][1] = 4;
  t2->t[0][2] = 5;
  t2->t[1][0] = 5;
  t2->t[1][1] = 6;
  t2->t[1][2] = 13;
  t2->t[2][0] = 6;
  t2->t[2][1] = 1;
  t2->t[2][2] = 7;
  t2->t[3][0] = 7;
  t2->t[3][1] = 8;
  t2->t[3][2] = 14;
  t2->t[4][0] = 8;
  t2->t[4][1] = 2;
  t2->t[4][2] = 9;
  t2->t[5][0] = 9;
  t2->t[5][1] = 10;
  t2->t[5][2] = 15;
  t2->t[6][0] = 10;
  t2->t[6][1] = 3;
  t2->t[6][2] = 11;
  t2->t[7][0] = 11;
  t2->t[7][1] = 12;
  t2->t[7][2] = 16;
  t2->t[8][0] = 17;
  t2->t[8][1] = 13;
  t2->t[8][2] = 18;
  t2->t[9][0] = 18;
  t2->t[9][1] = 19;
  t2->t[9][2] = 26;
  t2->t[10][0] = 19;
  t2->t[10][1] = 14;
  t2->t[10][2] = 20;
  t2->t[11][0] = 20;
  t2->t[11][1] = 21;
  t2->t[11][2] = 27;
  t2->t[12][0] = 21;
  t2->t[12][1] = 15;
  t2->t[12][2] = 22;
  t2->t[13][0] = 22;
  t2->t[13][1] = 23;
  t2->t[13][2] = 28;
  t2->t[14][0] = 23;
  t2->t[14][1] = 16;
  t2->t[14][2] = 24;
  t2->t[15][0] = 24;
  t2->t[15][1] = 25;
  t2->t[15][2] = 29;
  t2->t[16][0] = 30;
  t2->t[16][1] = 26;
  t2->t[16][2] = 31;
  t2->t[17][0] = 31;
  t2->t[17][1] = 32;
  t2->t[17][2] = 36;
  t2->t[18][0] = 33;
  t2->t[18][1] = 29;
  t2->t[18][2] = 34;
  t2->t[19][0] = 34;
  t2->t[19][1] = 35;
  t2->t[19][2] = 37;
  t2->t[20][0] = 38;
  t2->t[20][1] = 36;
  t2->t[20][2] = 39;
  t2->t[21][0] = 39;
  t2->t[21][1] = 40;
  t2->t[21][2] = 44;
  t2->t[22][0] = 41;
  t2->t[22][1] = 37;
  t2->t[22][2] = 42;
  t2->t[23][0] = 42;
  t2->t[23][1] = 43;
  t2->t[23][2] = 45;

  t2->vertices = 23;
  t2->edges = 46;
  t2->triangles = 24;

  return t2;
}

void
fixnormals_tri2d(ptri2d gr)
{
  const     real(*x)[2] = (const real(*)[2]) gr->x;
  uint(*e)[2] = gr->e;
  const     uint(*t)[3] = (const uint(*)[3]) gr->t;
  const uint *eb = gr->eb;
  uint      triangles = gr->triangles;
  real      n[2], xi[2];
  uint      k, i, ei, ej, v;

  for (k = 0; k < triangles; k++)
    for (i = 0; i < 3; i++) {
      ei = t[k][i];
      if (eb[ei]) {
	/* Compute clockwise normal of edge */
	n[0] = x[e[ei][1]][1] - x[e[ei][0]][1];
	n[1] = x[e[ei][0]][0] - x[e[ei][1]][0];

	/* Find point opposite the edge */
	ej = t[k][(i + 1) % 3];
	if (e[ej][0] == e[ei][0] || e[ej][0] == e[ei][1]) {
	  assert(e[ej][1] != e[ei][0] && e[ej][1] != e[ei][1]);

	  xi[0] = x[e[ej][1]][0] - x[e[ej][0]][0];
	  xi[1] = x[e[ej][1]][1] - x[e[ej][0]][1];
	}
	else {
	  xi[0] = x[e[ej][0]][0] - x[e[ej][1]][0];
	  xi[1] = x[e[ej][0]][1] - x[e[ej][1]][1];
	}

	/* Swap start and end vertex if inner product positive */
	if (xi[0] * n[0] + xi[1] * n[1] > 0) {
	  v = e[ei][0];
	  e[ei][0] = e[ei][1];
	  e[ei][1] = v;
	}
      }
    }
}

void
check_tri2d(pctri2d t2)
{
  uint      vertices = t2->vertices;
  uint      edges = t2->edges;
  uint      triangles = t2->triangles;
  const     uint(*e)[2] = (const uint(*)[2]) t2->e;
  const     uint(*t)[3] = (const uint(*)[3]) t2->t;
  const uint *eb = (const uint *) t2->eb;
  const uint *xb = (const uint *) t2->xb;
  uint     *en;
  uint      v[6];
  uint      i, j, k, l;
  uint      errors;

  errors = 0;

  /* Check if each edge contains only correct vertices */
  for (i = 0; i < edges; i++) {
    if (e[i][0] >= vertices) {
      (void) printf(" Edge %u contains illegal vertex %u\n", i, e[i][0]);
      errors++;
    }
    if (e[i][1] >= vertices) {
      (void) printf(" Edge %u contains illegal vertex %u\n", i, e[i][1]);
      errors++;
    }
  }

  /* Check if each boundary edge contains only boundary vertices */
  for (i = 0; i < edges; i++) {
    if (eb[i] == 1) {
      if (xb[e[i][0]] == 0) {
	(void) printf(" Boundary edge %u contains non-boundary vertex %u\n",
		      i, e[i][0]);
	errors++;
      }
      if (xb[e[i][1]] == 0) {
	(void) printf(" Boundary edge %u contains non-boundary vertex %u\n",
		      i, e[i][1]);
	errors++;
      }
    }
  }

  /* Check if each triangle consists of three vertices */
  for (i = 0; i < triangles; i++) {
    j = 0;
    for (k = 0; k < 3; k++) {
      if (t[i][k] >= edges) {
	(void) printf(" Trinagle %u contains illegal edge %u\n", i, t[i][k]);
	errors++;
      }
      else {
	for (l = 0; l < j && e[t[i][k]][0] != v[l]; l++);
	if (l == j) {
	  v[j] = e[t[i][k]][0];
	  j++;
	}

	for (l = 0; l < j && e[t[i][k]][1] != v[l]; l++);
	if (l == j) {
	  v[j] = e[t[i][k]][1];
	  j++;
	}
      }
    }
    if (j != 3) {
      (void) printf(" Triangle %u consists of %u vertices:", i, j);
      while (j-- > 0)
	(void) printf(" %u", v[j]);
      (void) printf("\n");
      errors++;
    }
  }

  /* Check if each non-boundary edge is adjacent to two triangles */
  en = (uint *) allocmem(sizeof(uint) * edges);
  for (i = 0; i < edges; i++)
    en[i] = 0;
  for (i = 0; i < triangles; i++) {
    en[t[i][0]]++;
    en[t[i][1]]++;
    en[t[i][2]]++;
  }
  for (i = 0; i < edges; i++) {
    if (eb[i] == 0 && en[i] != 2) {
      (void) printf(" Non-boundary edge %u adjacent to %u triangles\n", i,
		    en[i]);
      errors++;
    }
    if (eb[i] == 1 && en[i] != 1) {
      (void) printf(" Boundary edge %u adjacent to %u triangles\n", i, en[i]);
      errors++;
    }
  }
  freemem(en);

  (void) printf(" %u errors found\n", errors);
}

/* ------------------------------------------------------------
 Get geometrical information
 ------------------------------------------------------------ */

/* Auxiliary function:
 find a common vertex in e1 and e2 */
static    uint
common_vertex(const uint(*e)[2], uint e1, uint e2)
{
  if (e[e1][0] == e[e2][0] || e[e1][0] == e[e2][1])
    return 0;			/*return 0, if start vertex of e1 is the common vertex */

  assert(e[e1][1] == e[e2][0] || e[e1][1] == e[e2][1]);
  return 1;			/*return 1, if end vertex of e1 is the common vertex */
}

static    uint
common_vertex_global(const uint(*e)[2], uint e1, uint e2)
{
  return e[e1][common_vertex(e, e1, e2)];	/*return number of common vertex in e1 and e2 */
}

/*Auxiliary function
 find the vertexes in triangular tn of grid t2 and save them in v*/
void
getvertices_tri2d(pctri2d t2, uint tn, uint v[])
{
  const     uint(*e)[2] = (const uint(*)[2]) t2->e;
  const     uint(*t)[3] = (const uint(*)[3]) t2->t;

  v[2] = common_vertex_global(e, t[tn][0], t[tn][1]);	/*common vertex in edge 0 and 1 in triangle tn */
  v[0] = common_vertex_global(e, t[tn][1], t[tn][2]);	/*common vertex in edge 1 and 2 in triangle tn */
  v[1] = common_vertex_global(e, t[tn][2], t[tn][0]);	/*common vertex in edge 2 and 0 in triangle tn */
}

void
write_tri2d(pctri2d t2, const char *name)
{
  FILE     *out;
  const     real(*x)[2] = (const real(*)[2]) t2->x;
  const     uint(*e)[2] = (const uint(*)[2]) t2->e;
  const     uint(*t)[3] = (const uint(*)[3]) t2->t;
  const uint *xb = (const uint *) t2->xb;
  const uint *eb = (const uint *) t2->eb;
  uint      vertices = t2->vertices;
  uint      edges = t2->edges;
  uint      triangles = t2->triangles;
  uint      i;

  out = fopen(name, "w");
  if (!out) {
    (void) fprintf(stderr, "Could not open file \"%s\" for writing\n", name);
    return;
  }

  (void) fprintf(out, "# Triangular mesh description\n"
		 "# Vertices, edges and triangles\n"
		 "%u %u %u\n", vertices, edges, triangles);

  (void) fprintf(out,
		 "# List of vertices, given by coordinates and boundary flags\n");
  for (i = 0; i < vertices; i++)
    (void) fprintf(out, "%.12e \t %.12e \t  %u\n", x[i][0], x[i][1], xb[i]);
  (void) fprintf(out,
		 "# List of edges, given by vertex numbers and boundary flags\n");
  for (i = 0; i < edges; i++)
    (void) fprintf(out, "%u \t %u\t  %u\n", e[i][0], e[i][1], eb[i]);
  (void) fprintf(out, "# List of triangles, given by edge numbers\n");
  for (i = 0; i < triangles; i++)
    (void) fprintf(out, "%u \t %u \t%u\n", t[i][0], t[i][1], t[i][2]);

  (void) fclose(out);
}

ptri2d
read_tri2d(const char *name)
{

  FILE     *in;
  ptri2d    t2;
  real(*x)[2];
  uint(*e)[2];
  uint(*t)[3];
  uint     *xb;
  uint     *eb;
  uint      i;

  uint      vertices, edges, triangles;
  uint      items;
  char      buf[80], *res;

  in = fopen(name, "r");
  if (!in) {
    (void) fprintf(stderr, "Could not open file \"%s\" for reading\n", name);
    return 0;
  }

  res = fgets(buf, BUFSIZE, in);
  assert(res != NULL);
  while (!feof(in) && buf[0] == '#') {
    res = fgets(buf, 80, in);
    assert(res != NULL);
  }
  items = sscanf(buf, "%u %u %u", &vertices, &edges, &triangles);
  if (items != 3) {
    (void) fprintf(stderr, "Could not get sizes from file \"%s\"\n", name);
    (void) fclose(in);
    return 0;
  }

  t2 = new_tri2d(vertices, edges, triangles);
  x = t2->x;
  e = t2->e;
  t = t2->t;
  xb = t2->xb;
  eb = t2->eb;

  /*vertices */
  for (i = 0; i < vertices; i++) {
    res = fgets(buf, 80, in);
    assert(res != NULL);
    while (!feof(in) && buf[0] == '#') {
      res = fgets(buf, 80, in);
      assert(res != NULL);
    }
    items = sscanf(buf, "%" SCANF_PREFIX "f %" SCANF_PREFIX "f %u", x[i],
		   x[i] + 1, xb + i);
    if (items != 3) {
      (void) fprintf(stderr, "Could not read vertex %u from file \"%s\"\n", i,
		     name);
      del_tri2d(t2);
      (void) fclose(in);
      return 0;
    }
  }
  /*edges */
  for (i = 0; i < edges; i++) {
    res = fgets(buf, 80, in);
    assert(res != NULL);
    while (!feof(in) && buf[0] == '#') {
      res = fgets(buf, 80, in);
      assert(res != NULL);
    }
    items = sscanf(buf, "%u %u %u", e[i], e[i] + 1, eb + i);
    if (items != 3) {
      (void) fprintf(stderr, "Could not read edge %u from file \"%s\"\n", i,
		     name);
      del_tri2d(t2);
      (void) fclose(in);
      return 0;
    }
  }
  /*triangles */
  for (i = 0; i < triangles; i++) {
    res = fgets(buf, 80, in);
    assert(res != NULL);
    while (!feof(in) && buf[0] == '#') {
      res = fgets(buf, 80, in);
      assert(res != NULL);
    }
    items = sscanf(buf, "%u %u %u", t[i], t[i] + 1, t[i] + 2);
    if (items != 3) {
      (void) fprintf(stderr, "Could not read triangle %u from file \"%s\"\n",
		     i, name);
      del_tri2d(t2);
      (void) fclose(in);
      return 0;
    }
  }

  (void) fclose(in);

  return t2;
}

/* ------------------------------------------------------------
 Display grid
 ------------------------------------------------------------ */

void
draw_cairo_tri2d(pctri2d t2, const char *filename, bool mark_refedges,
		 int mark_triangle)
{
#ifdef USE_CAIRO
  cairo_surface_t *surface;
  cairo_t  *cr;
  const     real(*x)[2] = (const real(*)[2]) t2->x;
  const     uint(*e)[2] = (const uint(*)[2]) t2->e;
  const     uint(*t)[3] = (const uint(*)[3]) t2->t;
  uint      vertices = t2->vertices;
  uint      edges = t2->edges;
  uint      triangles = t2->triangles;
  double    xmin, xmax, ymin, ymax;
  double    hmin, hk, x0, y0, dx, dy, nx, ny;
  uint      k;

  xmin = xmax = x[0][0];
  ymin = ymax = x[0][1];
  for (k = 1; k < vertices; k++) {
    if (x[k][0] < xmin)
      xmin = x[k][0];
    else if (x[k][0] > xmax)
      xmax = x[k][0];
    if (x[k][1] < ymin)
      ymin = x[k][1];
    else if (x[k][1] > ymax)
      ymax = x[k][1];
  }

  surface = cairo_pdf_surface_create(filename, 100.0 * (xmax - xmin),
				     100.0 * (ymax - ymin));
  cr = cairo_create(surface);
  cairo_surface_destroy(surface);
  cairo_translate(cr, -100.0 * xmin, -100.0 * ymin);
  cairo_scale(cr, 95.0, 95.0);
  cairo_set_line_width(cr, cairo_get_line_width(cr) / 95.0);
  cairo_scale(cr, 1.0, -1.0);

  cairo_set_source_rgb(cr, 1.0, 0.0, 0.0);
  if (mark_triangle >= 0) {
    k = mark_triangle;
    cairo_move_to(cr, x[e[t[k][0]][0]][0], x[e[t[k][0]][0]][1]);
    cairo_line_to(cr, x[e[t[k][0]][1]][0], x[e[t[k][0]][1]][1]);
    if (e[t[k][1]][0] == e[t[k][0]][0] || e[t[k][1]][0] == e[t[k][0]][1])
      cairo_line_to(cr, x[e[t[k][1]][1]][0], x[e[t[k][1]][1]][1]);
    else
      cairo_line_to(cr, x[e[t[k][1]][0]][0], x[e[t[k][1]][0]][1]);
    cairo_close_path(cr);
    cairo_fill(cr);
  }

  if (mark_refedges) {
    hmin =
      REAL_SQRT(REAL_SQR(x[e[0][0]][0] - x[e[0][1]][0]) +
		REAL_SQR(x[e[0][0]][1] - x[e[0][1]][1]));
    for (k = 1; k < edges; k++) {
      hk =
	REAL_SQRT(REAL_SQR(x[e[k][0]][0] - x[e[k][1]][0]) +
		  REAL_SQR(x[e[k][0]][1] - x[e[k][1]][1]));
      if (hk < hmin)
	hmin = hk;
    }

    hk = xmax - xmin;
    if (ymax - ymin < hk)
      hk = ymax - ymin;
    if (hmin > 0.05 * hk)
      hmin = 0.05 * hk;

    cairo_set_source_rgb(cr, 0.0, 0.0, 1.0);
    for (k = 0; k < triangles; k++) {
      cairo_move_to(cr, x[e[t[k][0]][0]][0], x[e[t[k][0]][0]][1]);
      cairo_line_to(cr, x[e[t[k][0]][1]][0], x[e[t[k][0]][1]][1]);

      nx = x[e[t[k][0]][1]][1] - x[e[t[k][0]][0]][1];
      ny = x[e[t[k][0]][0]][0] - x[e[t[k][0]][1]][0];
      hk = REAL_SQRT(REAL_SQR(nx) + REAL_SQR(ny));
      nx /= hk;
      ny /= hk;

      if (e[t[k][1]][0] == e[t[k][2]][0] || e[t[k][1]][0] == e[t[k][2]][1]) {
	x0 = x[e[t[k][1]][0]][0];
	y0 = x[e[t[k][1]][0]][1];
      }
      else {
	assert(e[t[k][1]][1] == e[t[k][2]][0]
	       || e[t[k][1]][1] == e[t[k][2]][1]);
	x0 = x[e[t[k][1]][1]][0];
	y0 = x[e[t[k][1]][1]][1];
      }

      dx = x0 - x[e[t[k][0]][1]][0];
      dy = y0 - x[e[t[k][0]][1]][1];
      hk = REAL_ABS(dx * nx + dy * ny);
      dx /= hk;
      dy /= hk;
      cairo_line_to(cr, x[e[t[k][0]][1]][0] + hmin * 0.2 * dx,
		    x[e[t[k][0]][1]][1] + hmin * 0.2 * dy);

      dx = x0 - x[e[t[k][0]][0]][0];
      dy = y0 - x[e[t[k][0]][0]][1];
      hk = REAL_ABS(dx * nx + dy * ny);
      dx /= hk;
      dy /= hk;
      cairo_line_to(cr, x[e[t[k][0]][0]][0] + hmin * 0.2 * dx,
		    x[e[t[k][0]][0]][1] + hmin * 0.2 * dy);

      cairo_close_path(cr);
      cairo_fill(cr);
    }
  }

  cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
  for (k = 0; k < edges; k++) {
    cairo_move_to(cr, x[e[k][0]][0], x[e[k][0]][1]);
    cairo_line_to(cr, x[e[k][1]][0], x[e[k][1]][1]);
    cairo_stroke(cr);
  }

  cairo_destroy(cr);
#else
  fprintf(stderr, "CAIRO not supported.\n");
#endif
}

/* ------------------------------------------------------------
 Regular refinement
 ------------------------------------------------------------ */

/* Auxiliary function:
 find two edges in {e1,e1+1} and {e2,e2+1} that have a common vertex */
static void
intersecting_edges(const uint(*e)[2], uint e1, uint e2, uint * e1s,
		   uint * e2s)
{
  if (e[e1][0] == e[e2][0] || e[e1][0] == e[e2][1] || e[e1][1] == e[e2][0]
      || e[e1][1] == e[e2][1]) {
    *e1s = e1;
    *e2s = e2;
    return;
  }

  if (e[e1][0] == e[e2 + 1][0] || e[e1][0] == e[e2 + 1][1]
      || e[e1][1] == e[e2 + 1][0] || e[e1][1] == e[e2 + 1][1]) {
    *e1s = e1;
    *e2s = e2 + 1;
    return;
  }

  if (e[e1 + 1][0] == e[e2][0] || e[e1 + 1][0] == e[e2][1]
      || e[e1 + 1][1] == e[e2][0] || e[e1 + 1][1] == e[e2][1]) {
    *e1s = e1 + 1;
    *e2s = e2;
    return;
  }

  assert(e[e1 + 1][0] == e[e2 + 1][0] || e[e1 + 1][0] == e[e2 + 1][1]
	 || e[e1 + 1][1] == e[e2 + 1][0] || e[e1 + 1][1] == e[e2 + 1][1]);
  *e1s = e1 + 1;
  *e2s = e2 + 1;
}

ptri2d
refine_tri2d(pctri2d t2, ptri2dref * t2r)
{
  ptri2d    r;
  uint      vertices = t2->vertices;
  uint      edges = t2->edges;
  uint      triangles = t2->triangles;
  uint      edge_vertices;	/* First vertex on an edge */
  uint      triangle_edges;	/* First edge on a triangle */
  uint      i, j;
  const     uint(*re)[2];

  uint     *xf, *xt, *ef, *et, *tf;
  xf = xt = ef = et = tf = 0;

  r = new_tri2d(vertices + edges, 2 * edges + 3 * triangles, 4 * triangles);
  re = (const uint(*)[2]) r->e;

  if (t2r) {
    *t2r = (ptri2dref) allocmem(sizeof(tri2dref));
    xf = (*t2r)->xf = (uint *) allocmem(sizeof(uint) * r->vertices);
    xt = (*t2r)->xt = (uint *) allocmem(sizeof(uint) * r->vertices);
    ef = (*t2r)->ef = (uint *) allocmem(sizeof(uint) * r->edges);
    et = (*t2r)->et = (uint *) allocmem(sizeof(uint) * r->edges);
    tf = (*t2r)->tf = (uint *) allocmem(sizeof(uint) * r->triangles);
  }

  /* Create vertices by copying old vertices */
  i = 0;
  for (j = 0; j < vertices; j++) {
    r->x[i][0] = t2->x[j][0];
    r->x[i][1] = t2->x[j][1];
    r->xb[i] = t2->xb[j];
    if (t2r) {
      xf[i] = j;
      xt[i] = 0;
    }
    i++;
  }

  /* Create vertices within edges */
  edge_vertices = i;
  for (j = 0; j < edges; j++) {
    r->x[i][0] = 0.5 * (t2->x[t2->e[j][0]][0] + t2->x[t2->e[j][1]][0]);
    r->x[i][1] = 0.5 * (t2->x[t2->e[j][0]][1] + t2->x[t2->e[j][1]][1]);
    r->xb[i] = t2->eb[j];
    if (t2r) {
      xf[i] = j;
      xt[i] = 1;
    }
    i++;
  }
  assert(i == r->vertices);

  /* Create edges by splitting old edges */
  i = 0;
  for (j = 0; j < edges; j++) {
    r->e[i][0] = t2->e[j][0];
    r->e[i][1] = edge_vertices + j;
    r->eb[i] = t2->eb[j];
    if (t2r) {
      ef[i] = j;
      et[i] = 1;
    }
    i++;

    r->e[i][0] = edge_vertices + j;
    r->e[i][1] = t2->e[j][1];
    r->eb[i] = t2->eb[j];
    if (t2r) {
      ef[i] = j;
      et[i] = 1;
    }
    i++;
  }

  /* Create edges within triangles */
  triangle_edges = i;
  for (j = 0; j < triangles; j++) {
    r->e[i][0] = edge_vertices + t2->t[j][1];
    r->e[i][1] = edge_vertices + t2->t[j][2];
    r->eb[i] = 0;
    if (t2r) {
      ef[i] = j;
      et[i] = 2;
    }
    i++;

    r->e[i][0] = edge_vertices + t2->t[j][2];
    r->e[i][1] = edge_vertices + t2->t[j][0];
    r->eb[i] = 0;
    if (t2r) {
      ef[i] = j;
      et[i] = 2;
    }
    i++;

    r->e[i][0] = edge_vertices + t2->t[j][0];
    r->e[i][1] = edge_vertices + t2->t[j][1];
    r->eb[i] = 0;
    if (t2r) {
      ef[i] = j;
      et[i] = 2;
    }
    i++;
  }

  /* Create triangles by splitting old triangles */
  i = 0;
  for (j = 0; j < triangles; j++) {
    intersecting_edges(re, 2 * t2->t[j][1], 2 * t2->t[j][2], r->t[i] + 1,
		       r->t[i] + 2);
    r->t[i][0] = triangle_edges + 3 * j;
    if (t2r) {
      tf[i] = j;
    }
    i++;

    intersecting_edges(re, 2 * t2->t[j][0], 2 * t2->t[j][2], r->t[i],
		       r->t[i] + 2);
    r->t[i][1] = triangle_edges + 3 * j + 1;
    if (t2r) {
      tf[i] = j;
    }
    i++;

    intersecting_edges(re, 2 * t2->t[j][0], 2 * t2->t[j][1], r->t[i],
		       r->t[i] + 1);
    r->t[i][2] = triangle_edges + 3 * j + 2;
    if (t2r) {
      tf[i] = j;
    }
    i++;

    r->t[i][0] = triangle_edges + 3 * j;
    r->t[i][1] = triangle_edges + 3 * j + 1;
    r->t[i][2] = triangle_edges + 3 * j + 2;
    if (t2r) {
      tf[i] = j;
    }
    i++;
  }

  return r;
}

void
smooth_unitcircle_tri2d(ptri2d t2)
{
  uint      i;
  real      normL2, normL1;

  /* Scaling vertices */
  for (i = 0; i < t2->vertices; i++) {
    normL1 = REAL_ABS(t2->x[i][0]) + REAL_ABS(t2->x[i][1]);
    normL2 = REAL_SQRT(t2->x[i][0] * t2->x[i][0] + t2->x[i][1] * t2->x[i][1]);
    if (normL1 != 0.0) {
      t2->x[i][0] = t2->x[i][0] * normL1 / normL2;
      t2->x[i][1] = t2->x[i][1] * normL1 / normL2;
    }
  }
}

void
del_tri2dref(ptri2dref t2r)
{
  freemem(t2r->tf);
  freemem(t2r->et);
  freemem(t2r->ef);
  freemem(t2r->xt);
  freemem(t2r->xf);
  freemem(t2r);
}

/* ------------------------------------------------------------
 * Tool for constructing meshes
 * ------------------------------------------------------------ */

typedef struct _edgeentry edgeentry;
typedef edgeentry *pedgeentry;

typedef struct _trientry trientry;
typedef trientry *ptrientry;

struct _tri2dbuilder {
  real(*x)[2];

  pedgeentry *elist;
  ptrientry tlist;

  uint      vertices;
  uint      edges;
  uint      triangles;
};

struct _edgeentry {
  uint      name;

  uint      v[2];

  uint      neighbours;

  pedgeentry next;
};

struct _trientry {
  uint      name;

  pedgeentry e[3];

  ptrientry next;
};

static    pedgeentry
new_edgeentry(uint name, uint v0, uint v1, pedgeentry next)
{
  pedgeentry e;

  e = (pedgeentry) allocmem(sizeof(edgeentry));
  e->name = name;
  e->v[0] = v0;
  e->v[1] = v1;
  e->neighbours = 0;
  e->next = next;

  return e;
}

static    ptrientry
new_trientry(uint name, pedgeentry e0, pedgeentry e1,
	     pedgeentry e2, ptrientry next)
{
  ptrientry t;

  t = (ptrientry) allocmem(sizeof(trientry));
  t->name = name;
  t->e[0] = e0;
  t->e[1] = e1;
  t->e[2] = e2;
  t->next = next;

  e0->neighbours++;
  e1->neighbours++;
  e2->neighbours++;

  return t;
}

ptri2dbuilder
new_tri2dbuilder(uint vertices)
{
  ptri2dbuilder tb;
  uint      i;

  tb = (ptri2dbuilder) allocmem(sizeof(tri2dbuilder));
  tb->x = (real(*)[2]) allocmem(sizeof(real[2]) * vertices);
  tb->elist = (pedgeentry *) allocmem(sizeof(pedgeentry) * vertices);
  tb->tlist = 0;

  tb->vertices = vertices;
  tb->edges = 0;
  tb->triangles = 0;

  for (i = 0; i < vertices; i++)
    tb->elist[i] = 0;

  return tb;
}

static void
del_trientry(ptrientry t)
{
  ptrientry next;

  while (t) {
    next = t->next;

    freemem(t);

    t = next;
  }
}

static void
del_edgeentry(pedgeentry e)
{
  pedgeentry next;

  while (e) {
    next = e->next;

    freemem(e);

    e = next;
  }
}

void
del_tri2dbuilder(ptri2dbuilder tb)
{
  uint      i;

  for (i = 0; i < tb->vertices; i++)
    del_edgeentry(tb->elist[i]);

  del_trientry(tb->tlist);

  freemem(tb->elist);
  freemem(tb->x);
  freemem(tb);
}

real(*getx_tri2dbuilder(ptri2dbuilder tb))[2]
{
  return tb->x;
}

static    pedgeentry
find_edgeentry(ptri2dbuilder tb, uint v0, uint v1)
{
  pedgeentry e;
  uint      vx;

  assert(v0 < tb->vertices);
  assert(v1 < tb->vertices);

  /* Vertices in an edge are numbered in ascending order
   * (fixnormals_tri2d can be used later to ensure ccw order) */
  if (v0 > v1) {
    vx = v0;
    v0 = v1;
    v1 = vx;
  }

  /* Check whether this edge is already known */
  for (e = tb->elist[v0]; e; e = e->next)
    if (e->v[1] == v1)
      break;

  /* If not known, create it */
  if (e == 0) {
    tb->elist[v0] = e = new_edgeentry(tb->edges, v0, v1, tb->elist[v0]);
    tb->edges++;
  }

  return e;
}

void
addtriangle_tri2dbuilder(ptri2dbuilder tb, uint v0, uint v1, uint v2)
{
  pedgeentry e01, e12, e20;

  e01 = find_edgeentry(tb, v0, v1);
  e12 = find_edgeentry(tb, v1, v2);
  e20 = find_edgeentry(tb, v2, v0);

  tb->tlist = new_trientry(tb->triangles, e12, e20, e01, tb->tlist);
  tb->triangles++;
}

ptri2d
buildmesh_tri2dbuilder(ptri2dbuilder tb)
{
  ptri2d    gr;
  pedgeentry e;
  ptrientry t;
  uint      i;

  /* Create mesh */
  gr = new_tri2d(tb->vertices, tb->edges, tb->triangles);

  /* Reset vertex boundary flags */
  for (i = 0; i < tb->vertices; i++)
    gr->xb[i] = false;

  /* Copy triangles */
  for (t = tb->tlist; t; t = t->next) {
    gr->t[t->name][0] = t->e[0]->name;
    gr->t[t->name][1] = t->e[1]->name;
    gr->t[t->name][2] = t->e[2]->name;
  }

  /* Copy edges */
  for (i = 0; i < tb->vertices; i++)
    for (e = tb->elist[i]; e; e = e->next) {
      gr->e[e->name][0] = e->v[0];
      gr->e[e->name][1] = e->v[1];

      /* Determine whether this is a boundary edge or not */
      assert(e->neighbours == 1 || e->neighbours == 2);
      gr->eb[e->name] = (e->neighbours == 1);

      /* If it is a boundary edge, its vertices are boundary vertices */
      if (e->neighbours == 1) {
	gr->xb[e->v[0]] = true;
	gr->xb[e->v[1]] = true;
      }
    }

  /* Copy vertices */
  for (i = 0; i < tb->vertices; i++) {
    gr->x[i][0] = tb->x[i][0];
    gr->x[i][1] = tb->x[i][1];
  }

  /* Ensure that boundary edges are oriented counter-clockwise */
  fixnormals_tri2d(gr);

  return gr;
}
