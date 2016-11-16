/* ------------------------------------------------------------
 This is the file "macrosurface3d.c" of the H2Lib package.
 All rights reserved, Steffen Boerm 2010
 ------------------------------------------------------------ */

#include "macrosurface3d.h"

#include <assert.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#include "basic.h"

/* ------------------------------------------------------------
 Constructor and destructor
 ------------------------------------------------------------ */

pmacrosurface3d
new_macrosurface3d(uint vertices, uint edges, uint triangles)
{
  pmacrosurface3d mg;

  mg = (pmacrosurface3d) allocmem(sizeof(macrosurface3d));
  mg->x = (real(*)[3]) allocmem((size_t) sizeof(real[3]) * vertices);
  mg->e = (uint(*)[2]) allocmem((size_t) sizeof(uint[2]) * edges);
  mg->t = (uint(*)[3]) allocmem((size_t) sizeof(uint[3]) * triangles);
  mg->s = (uint(*)[3]) allocmem((size_t) sizeof(uint[3]) * triangles);

  mg->vertices = vertices;
  mg->edges = edges;
  mg->triangles = triangles;

  mg->phi = 0;
  mg->phidata = 0;

  return mg;
}

void
del_macrosurface3d(pmacrosurface3d mg)
{
  freemem(mg->s);
  freemem(mg->t);
  freemem(mg->e);
  freemem(mg->x);
  freemem(mg);
}

/* ------------------------------------------------------------
 Examples
 ------------------------------------------------------------ */

static void
sphere_parametrization(uint i, real xr1, real xr2, void *data, real xt[3])
{
  pcmacrosurface3d mg = (pcmacrosurface3d) data;
  const     real(*x)[3] = (const real(*)[3]) mg->x;
  const     uint(*t)[3] = (const uint(*)[3]) mg->t;
  real      norm;

  assert(i < mg->triangles);
  assert(t[i][0] < mg->vertices);
  assert(t[i][1] < mg->vertices);
  assert(t[i][2] < mg->vertices);

  xt[0] = (x[t[i][0]][0] * (1.0 - xr1 - xr2) + x[t[i][1]][0] * xr1
	   + x[t[i][2]][0] * xr2);
  xt[1] = (x[t[i][0]][1] * (1.0 - xr1 - xr2) + x[t[i][1]][1] * xr1
	   + x[t[i][2]][1] * xr2);
  xt[2] = (x[t[i][0]][2] * (1.0 - xr1 - xr2) + x[t[i][1]][2] * xr1
	   + x[t[i][2]][2] * xr2);

  norm = REAL_NORM3(xt[0], xt[1], xt[2]);
  xt[0] /= norm;
  xt[1] /= norm;
  xt[2] /= norm;
}

pmacrosurface3d
new_sphere_macrosurface3d()
{
  pmacrosurface3d mg;

  mg = new_macrosurface3d(6, 12, 8);

  /* Front */
  mg->x[0][0] = 0.0;
  mg->x[0][1] = 0.0;
  mg->x[0][2] = 1.0;
  /* Right */
  mg->x[1][0] = 1.0;
  mg->x[1][1] = 0.0;
  mg->x[1][2] = 0.0;
  /* Back */
  mg->x[2][0] = 0.0;
  mg->x[2][1] = 0.0;
  mg->x[2][2] = -1.0;
  /* Left */
  mg->x[3][0] = -1.0;
  mg->x[3][1] = 0.0;
  mg->x[3][2] = 0.0;
  /* Top */
  mg->x[4][0] = 0.0;
  mg->x[4][1] = 1.0;
  mg->x[4][2] = 0.0;
  /* Bottom */
  mg->x[5][0] = 0.0;
  mg->x[5][1] = -1.0;
  mg->x[5][2] = 0.0;

  /* Equator */
  mg->e[0][0] = 0;
  mg->e[0][1] = 1;
  mg->e[1][0] = 1;
  mg->e[1][1] = 2;
  mg->e[2][0] = 2;
  mg->e[2][1] = 3;
  mg->e[3][0] = 3;
  mg->e[3][1] = 0;
  /* Top */
  mg->e[4][0] = 0;
  mg->e[4][1] = 4;
  mg->e[5][0] = 1;
  mg->e[5][1] = 4;
  mg->e[6][0] = 2;
  mg->e[6][1] = 4;
  mg->e[7][0] = 3;
  mg->e[7][1] = 4;
  /* Bottom */
  mg->e[8][0] = 0;
  mg->e[8][1] = 5;
  mg->e[9][0] = 1;
  mg->e[9][1] = 5;
  mg->e[10][0] = 2;
  mg->e[10][1] = 5;
  mg->e[11][0] = 3;
  mg->e[11][1] = 5;

  /* Right front top */
  mg->t[0][0] = 0;
  mg->t[0][1] = 1;
  mg->t[0][2] = 4;
  mg->s[0][0] = 5;
  mg->s[0][1] = 4;
  mg->s[0][2] = 0;
  /* Right back top */
  mg->t[1][0] = 1;
  mg->t[1][1] = 2;
  mg->t[1][2] = 4;
  mg->s[1][0] = 6;
  mg->s[1][1] = 5;
  mg->s[1][2] = 1;
  /* Left back top */
  mg->t[2][0] = 2;
  mg->t[2][1] = 3;
  mg->t[2][2] = 4;
  mg->s[2][0] = 7;
  mg->s[2][1] = 6;
  mg->s[2][2] = 2;
  /* Left front top */
  mg->t[3][0] = 3;
  mg->t[3][1] = 0;
  mg->t[3][2] = 4;
  mg->s[3][0] = 4;
  mg->s[3][1] = 7;
  mg->s[3][2] = 3;
  /* Right front bottom */
  mg->t[4][0] = 0;
  mg->t[4][1] = 5;
  mg->t[4][2] = 1;
  mg->s[4][0] = 9;
  mg->s[4][1] = 0;
  mg->s[4][2] = 8;
  /* Right back bottom */
  mg->t[5][0] = 1;
  mg->t[5][1] = 5;
  mg->t[5][2] = 2;
  mg->s[5][0] = 10;
  mg->s[5][1] = 1;
  mg->s[5][2] = 9;
  /* Left back bottom */
  mg->t[6][0] = 2;
  mg->t[6][1] = 5;
  mg->t[6][2] = 3;
  mg->s[6][0] = 11;
  mg->s[6][1] = 2;
  mg->s[6][2] = 10;
  /* Left front bottom */
  mg->t[7][0] = 3;
  mg->t[7][1] = 5;
  mg->t[7][2] = 0;
  mg->s[7][0] = 8;
  mg->s[7][1] = 3;
  mg->s[7][2] = 11;

  mg->phi = sphere_parametrization;
  mg->phidata = mg;

  return mg;
}

static void
parabolic_mirror_parametrization(uint i, real xr1, real xr2,
				 void *data, real xt[3])
{
  pcmacrosurface3d mg = (pcmacrosurface3d) data;
  const     real(*x)[3] = (const real(*)[3]) mg->x;
  const     uint(*t)[3] = (const uint(*)[3]) mg->t;
  real      norm;

  assert(i < mg->triangles);
  assert(t[i][0] < mg->vertices);
  assert(t[i][1] < mg->vertices);
  assert(t[i][2] < mg->vertices);

  xt[0] = (x[t[i][0]][0] * (1.0 - xr1 - xr2) + x[t[i][1]][0] * xr1
	   + x[t[i][2]][0] * xr2);
  xt[1] = (x[t[i][0]][1] * (1.0 - xr1 - xr2) + x[t[i][1]][1] * xr1
	   + x[t[i][2]][1] * xr2);
  xt[2] = (x[t[i][0]][2] * (1.0 - xr1 - xr2) + x[t[i][1]][2] * xr1
	   + x[t[i][2]][2] * xr2);

  norm = REAL_NORM3(xt[0], xt[1], xt[2]);
  xt[0] /= norm;
  xt[0] = xt[0] < 0.0 ? xt[0] * xt[0] : xt[0];
  xt[1] /= norm;
  xt[2] /= norm;
}

pmacrosurface3d
new_parabolic_mirror_macrosurface3d()
{
  pmacrosurface3d mg;

  mg = new_sphere_macrosurface3d();
  mg->phi = parabolic_mirror_parametrization;

  return mg;
}

static void
cube_parametrization(uint i, real xr1, real xr2, void *data, real xt[3])
{
  pcmacrosurface3d mg = (pcmacrosurface3d) data;
  const     real(*x)[3] = (const real(*)[3]) mg->x;
  const     uint(*t)[3] = (const uint(*)[3]) mg->t;

  assert(i < mg->triangles);
  assert(t[i][0] < mg->vertices);
  assert(t[i][1] < mg->vertices);
  assert(t[i][2] < mg->vertices);

  xt[0] = (x[t[i][0]][0] * (1.0 - xr1 - xr2) + x[t[i][1]][0] * xr1
	   + x[t[i][2]][0] * xr2);
  xt[1] = (x[t[i][0]][1] * (1.0 - xr1 - xr2) + x[t[i][1]][1] * xr1
	   + x[t[i][2]][1] * xr2);
  xt[2] = (x[t[i][0]][2] * (1.0 - xr1 - xr2) + x[t[i][1]][2] * xr1
	   + x[t[i][2]][2] * xr2);
}

pmacrosurface3d
new_cuboid_macrosurface3d(real ax, real bx, real ay, real by,
			  real az, real bz)
{
  real(*x)[3];
  uint(*e)[2];
  uint(*t)[3];
  uint(*s)[3];
  uint      nx, ne, nt;

  pmacrosurface3d mg;

  mg = new_macrosurface3d(8, 18, 12);

  x = mg->x;
  e = mg->e;
  t = mg->t;
  s = mg->s;

  /****************************************************
   * vertices
   ****************************************************/

  nx = 0;

  /* left */
  x[nx][0] = ax, x[nx][1] = ay, x[nx++][2] = az;
  x[nx][0] = ax, x[nx][1] = by, x[nx++][2] = az;
  x[nx][0] = ax, x[nx][1] = ay, x[nx++][2] = bz;
  x[nx][0] = ax, x[nx][1] = by, x[nx++][2] = bz;

  /* right */
  x[nx][0] = bx, x[nx][1] = ay, x[nx++][2] = az;
  x[nx][0] = bx, x[nx][1] = by, x[nx++][2] = az;
  x[nx][0] = bx, x[nx][1] = ay, x[nx++][2] = bz;
  x[nx][0] = bx, x[nx][1] = by, x[nx++][2] = bz;

  /****************************************************
   * edges
   ****************************************************/

  ne = 0;

  /* left */
  e[ne][0] = 1, e[ne++][1] = 0;
  e[ne][0] = 0, e[ne++][1] = 3;
  e[ne][0] = 3, e[ne++][1] = 1;
  e[ne][0] = 2, e[ne++][1] = 0;
  e[ne][0] = 3, e[ne++][1] = 2;

  /* right */
  e[ne][0] = 5, e[ne++][1] = 7;
  e[ne][0] = 7, e[ne++][1] = 4;
  e[ne][0] = 4, e[ne++][1] = 5;
  e[ne][0] = 4, e[ne++][1] = 6;
  e[ne][0] = 6, e[ne++][1] = 7;

  /* front */
  e[ne][0] = 2, e[ne++][1] = 6;
  e[ne][0] = 6, e[ne++][1] = 3;
  e[ne][0] = 3, e[ne++][1] = 7;

  /* back */
  e[ne][0] = 5, e[ne++][1] = 1;
  e[ne][0] = 1, e[ne++][1] = 4;
  e[ne][0] = 4, e[ne++][1] = 0;

  /* top */
  e[ne][0] = 1, e[ne++][1] = 7;

  /* bottom */
  e[ne][0] = 0, e[ne++][1] = 6;

  /****************************************************
   * triangles
   ****************************************************/

  nt = 0;

  /* left */
  t[nt][0] = 1, t[nt][1] = 0, t[nt][2] = 3;
  s[nt][0] = 1, s[nt][1] = 2, s[nt++][2] = 0;

  t[nt][0] = 3, t[nt][1] = 0, t[nt][2] = 2;
  s[nt][0] = 3, s[nt][1] = 4, s[nt++][2] = 1;

  /* right */
  t[nt][0] = 4, t[nt][1] = 5, t[nt][2] = 7;
  s[nt][0] = 5, s[nt][1] = 6, s[nt++][2] = 7;

  t[nt][0] = 4, t[nt][1] = 7, t[nt][2] = 6;
  s[nt][0] = 9, s[nt][1] = 8, s[nt++][2] = 6;

  /* front */
  t[nt][0] = 3, t[nt][1] = 2, t[nt][2] = 6;
  s[nt][0] = 10, s[nt][1] = 11, s[nt++][2] = 4;

  t[nt][0] = 3, t[nt][1] = 6, t[nt][2] = 7;
  s[nt][0] = 9, s[nt][1] = 12, s[nt++][2] = 11;

  /* back */
  t[nt][0] = 4, t[nt][1] = 0, t[nt][2] = 1;
  s[nt][0] = 0, s[nt][1] = 14, s[nt++][2] = 15;

  t[nt][0] = 4, t[nt][1] = 1, t[nt][2] = 5;
  s[nt][0] = 13, s[nt][1] = 7, s[nt++][2] = 14;

  /* bottom */
  t[nt][0] = 0, t[nt][1] = 6, t[nt][2] = 2;
  s[nt][0] = 10, s[nt][1] = 3, s[nt++][2] = 17;

  t[nt][0] = 4, t[nt][1] = 6, t[nt][2] = 0;
  s[nt][0] = 17, s[nt][1] = 15, s[nt++][2] = 8;

  /* top */
  t[nt][0] = 1, t[nt][1] = 3, t[nt][2] = 7;
  s[nt][0] = 12, s[nt][1] = 16, s[nt++][2] = 2;

  t[nt][0] = 1, t[nt][1] = 7, t[nt][2] = 5;
  s[nt][0] = 5, s[nt][1] = 13, s[nt++][2] = 16;

  mg->phi = cube_parametrization;
  mg->phidata = mg;

  return mg;
}

pmacrosurface3d
new_cube_macrosurface3d()
{
  return new_cuboid_macrosurface3d(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0);
}

static void
cylinder_parametrization(uint i, real xr1, real xr2, void *data, real xt[3])
{
  pcmacrosurface3d mg = (pcmacrosurface3d) data;
  const     real(*x)[3] = (const real(*)[3]) mg->x;
  const     uint(*t)[3] = (const uint(*)[3]) mg->t;
  real      normL1, normL2;

  assert(i < mg->triangles);
  assert(t[i][0] < mg->vertices);
  assert(t[i][1] < mg->vertices);
  assert(t[i][2] < mg->vertices);

  xt[0] = (x[t[i][0]][0] * (1.0 - xr1 - xr2) + x[t[i][1]][0] * xr1
	   + x[t[i][2]][0] * xr2);
  xt[1] = (x[t[i][0]][1] * (1.0 - xr1 - xr2) + x[t[i][1]][1] * xr1
	   + x[t[i][2]][1] * xr2);
  xt[2] = (x[t[i][0]][2] * (1.0 - xr1 - xr2) + x[t[i][1]][2] * xr1
	   + x[t[i][2]][2] * xr2);

  if (REAL_ABS(xt[0]) != 4.0) {
    normL2 = REAL_NORM2(xt[1], xt[2]);
    xt[1] /= normL2;
    xt[2] /= normL2;
  }
  else {
    normL1 = REAL_ABS(xt[1]) + REAL_ABS(xt[2]);
    normL2 = REAL_NORM2(xt[1], xt[2]);
    if (normL2 > 0.0) {
      xt[1] = xt[1] * normL1 / normL2;
      xt[2] = xt[2] * normL1 / normL2;
    }
  }

}

pmacrosurface3d
new_cylinder_macrosurface3d()
{
  pmacrosurface3d mg;
  real(*x)[3];
  uint(*e)[2];
  uint(*t)[3];
  uint(*s)[3];

  mg = new_macrosurface3d(22, 60, 40);

  x = mg->x;
  e = mg->e;
  t = mg->t;
  s = mg->s;

  /****************************************************
   * vertices
   ****************************************************/

  /* left circle */
  x[0][0] = -4.0, x[0][1] = 0.0, x[0][2] = 1.0;
  x[1][0] = -4.0, x[1][1] = 1.0, x[1][2] = 0.0;
  x[2][0] = -4.0, x[2][1] = 0.0, x[2][2] = -1.0;
  x[3][0] = -4.0, x[3][1] = -1.0, x[3][2] = 0.0;
  x[4][0] = -4.0, x[4][1] = 0.0, x[4][2] = 0.0;

  /* 1st mid circle */
  x[5][0] = -2.0, x[5][1] = REAL_SQRT(0.5), x[5][2] = REAL_SQRT(0.5);
  x[6][0] = -2.0, x[6][1] = REAL_SQRT(0.5), x[6][2] = -REAL_SQRT(0.5);
  x[7][0] = -2.0, x[7][1] = -REAL_SQRT(0.5), x[7][2] = -REAL_SQRT(0.5);
  x[8][0] = -2.0, x[8][1] = -REAL_SQRT(0.5), x[8][2] = REAL_SQRT(0.5);

  /* 2nd mid circle */
  x[9][0] = 0.0, x[9][1] = 0.0, x[9][2] = 1.0;
  x[10][0] = 0.0, x[10][1] = 1.0, x[10][2] = 0.0;
  x[11][0] = 0.0, x[11][1] = 0.0, x[11][2] = -1.0;
  x[12][0] = 0.0, x[12][1] = -1.0, x[12][2] = 0.0;

  /* 3rd mid circle */
  x[13][0] = 2.0, x[13][1] = REAL_SQRT(0.5), x[13][2] = REAL_SQRT(0.5);
  x[14][0] = 2.0, x[14][1] = REAL_SQRT(0.5), x[14][2] = -REAL_SQRT(0.5);
  x[15][0] = 2.0, x[15][1] = -REAL_SQRT(0.5), x[15][2] = -REAL_SQRT(0.5);
  x[16][0] = 2.0, x[16][1] = -REAL_SQRT(0.5), x[16][2] = REAL_SQRT(0.5);

  /* right circle */
  x[17][0] = 4.0, x[17][1] = 0.0, x[17][2] = 1.0;
  x[18][0] = 4.0, x[18][1] = 1.0, x[18][2] = 0.0;
  x[19][0] = 4.0, x[19][1] = 0.0, x[19][2] = -1.0;
  x[20][0] = 4.0, x[20][1] = -1.0, x[20][2] = 0.0;
  x[21][0] = 4.0, x[21][1] = 0.0, x[21][2] = 0.0;

  /****************************************************
   * edges
   ****************************************************/

  /* left circle */
  e[0][0] = 0, e[0][1] = 1;
  e[1][0] = 1, e[1][1] = 4;
  e[2][0] = 4, e[2][1] = 0;

  e[3][0] = 4, e[3][1] = 2;
  e[4][0] = 2, e[4][1] = 1;

  e[5][0] = 2, e[5][1] = 3;
  e[6][0] = 3, e[6][1] = 4;

  e[7][0] = 0, e[7][1] = 3;

  /* left to left mid circle */
  e[8][0] = 0, e[8][1] = 5;
  e[9][0] = 0, e[9][1] = 8;

  e[10][0] = 1, e[10][1] = 5;
  e[11][0] = 1, e[11][1] = 6;

  e[12][0] = 2, e[12][1] = 6;
  e[13][0] = 2, e[13][1] = 7;

  e[14][0] = 3, e[14][1] = 7;
  e[15][0] = 3, e[15][1] = 8;

  /* left mid to mid circle */
  e[16][0] = 5, e[16][1] = 9;
  e[17][0] = 5, e[17][1] = 10;

  e[18][0] = 6, e[18][1] = 10;
  e[19][0] = 6, e[19][1] = 11;

  e[20][0] = 7, e[20][1] = 11;
  e[21][0] = 7, e[21][1] = 12;

  e[22][0] = 8, e[22][1] = 9;
  e[23][0] = 8, e[23][1] = 12;

  /* left to mid circle */
  e[24][0] = 0, e[24][1] = 9;
  e[25][0] = 1, e[25][1] = 10;
  e[26][0] = 2, e[26][1] = 11;
  e[27][0] = 3, e[27][1] = 12;

  /* right to mid circle */
  e[28][0] = 17, e[28][1] = 9;
  e[29][0] = 18, e[29][1] = 10;
  e[30][0] = 19, e[30][1] = 11;
  e[31][0] = 20, e[31][1] = 12;

  /* right mid to mid circle */
  e[32][0] = 13, e[32][1] = 9;
  e[33][0] = 13, e[33][1] = 10;

  e[34][0] = 14, e[34][1] = 10;
  e[35][0] = 14, e[35][1] = 11;

  e[36][0] = 15, e[36][1] = 11;
  e[37][0] = 15, e[37][1] = 12;

  e[38][0] = 16, e[38][1] = 9;
  e[39][0] = 16, e[39][1] = 12;

  /* right to right mid circle */
  e[40][0] = 17, e[40][1] = 13;
  e[41][0] = 17, e[41][1] = 16;

  e[42][0] = 18, e[42][1] = 13;
  e[43][0] = 18, e[43][1] = 14;

  e[44][0] = 19, e[44][1] = 14;
  e[45][0] = 19, e[45][1] = 15;

  e[46][0] = 20, e[46][1] = 15;
  e[47][0] = 20, e[47][1] = 16;

  /* right circle */
  e[48][0] = 17, e[48][1] = 18;
  e[49][0] = 18, e[49][1] = 21;
  e[50][0] = 21, e[50][1] = 17;

  e[51][0] = 21, e[51][1] = 19;
  e[52][0] = 19, e[52][1] = 18;

  e[53][0] = 19, e[53][1] = 20;
  e[54][0] = 20, e[54][1] = 21;

  e[55][0] = 17, e[55][1] = 20;

  /* left mid to right mid circle */
  e[56][0] = 5, e[56][1] = 13;
  e[57][0] = 6, e[57][1] = 14;
  e[58][0] = 7, e[58][1] = 15;
  e[59][0] = 8, e[59][1] = 16;

  /****************************************************
   * triangles
   ****************************************************/

  /* left circle */
  t[0][0] = 0, t[0][1] = 1, t[0][2] = 4;
  s[0][0] = 1, s[0][1] = 2, s[0][2] = 0;

  t[1][0] = 1, t[1][1] = 2, t[1][2] = 4;
  s[1][0] = 3, s[1][1] = 1, s[1][2] = 4;

  t[2][0] = 2, t[2][1] = 3, t[2][2] = 4;
  s[2][0] = 6, s[2][1] = 3, s[2][2] = 5;

  t[3][0] = 0, t[3][1] = 4, t[3][2] = 3;
  s[3][0] = 6, s[3][1] = 7, s[3][2] = 2;

  /* left to left mid circle */
  t[4][0] = 0, t[4][1] = 5, t[4][2] = 1;
  s[4][0] = 10, s[4][1] = 0, s[4][2] = 8;
  t[5][0] = 1, t[5][1] = 6, t[5][2] = 2;
  s[5][0] = 12, s[5][1] = 4, s[5][2] = 11;
  t[6][0] = 2, t[6][1] = 7, t[6][2] = 3;
  s[6][0] = 14, s[6][1] = 5, s[6][2] = 13;
  t[7][0] = 3, t[7][1] = 8, t[7][2] = 0;
  s[7][0] = 9, s[7][1] = 7, s[7][2] = 15;

  /* left mid to mid circle */
  t[8][0] = 0, t[8][1] = 9, t[8][2] = 5;
  s[8][0] = 16, s[8][1] = 8, s[8][2] = 24;
  t[9][0] = 0, t[9][1] = 8, t[9][2] = 9;
  s[9][0] = 22, s[9][1] = 24, s[9][2] = 9;

  t[10][0] = 1, t[10][1] = 5, t[10][2] = 10;
  s[10][0] = 17, s[10][1] = 25, s[10][2] = 10;
  t[11][0] = 1, t[11][1] = 10, t[11][2] = 6;
  s[11][0] = 18, s[11][1] = 11, s[11][2] = 25;

  t[12][0] = 2, t[12][1] = 6, t[12][2] = 11;
  s[12][0] = 19, s[12][1] = 26, s[12][2] = 12;
  t[13][0] = 2, t[13][1] = 11, t[13][2] = 7;
  s[13][0] = 20, s[13][1] = 13, s[13][2] = 26;

  t[14][0] = 3, t[14][1] = 7, t[14][2] = 12;
  s[14][0] = 21, s[14][1] = 27, s[14][2] = 14;
  t[15][0] = 3, t[15][1] = 12, t[15][2] = 8;
  s[15][0] = 23, s[15][1] = 15, s[15][2] = 27;

  /* left mid to right mid circle */
  t[16][0] = 5, t[16][1] = 9, t[16][2] = 13;
  s[16][0] = 32, s[16][1] = 56, s[16][2] = 16;
  t[17][0] = 5, t[17][1] = 13, t[17][2] = 10;
  s[17][0] = 33, s[17][1] = 17, s[17][2] = 56;

  t[18][0] = 6, t[18][1] = 10, t[18][2] = 14;
  s[18][0] = 34, s[18][1] = 57, s[18][2] = 18;
  t[19][0] = 6, t[19][1] = 14, t[19][2] = 11;
  s[19][0] = 35, s[19][1] = 19, s[19][2] = 57;

  t[20][0] = 7, t[20][1] = 11, t[20][2] = 15;
  s[20][0] = 36, s[20][1] = 58, s[20][2] = 20;
  t[21][0] = 7, t[21][1] = 15, t[21][2] = 12;
  s[21][0] = 37, s[21][1] = 21, s[21][2] = 58;

  t[22][0] = 8, t[22][1] = 12, t[22][2] = 16;
  s[22][0] = 39, s[22][1] = 59, s[22][2] = 23;
  t[23][0] = 8, t[23][1] = 16, t[23][2] = 9;
  s[23][0] = 38, s[23][1] = 22, s[23][2] = 59;

  /* mid to right circle */
  t[24][0] = 9, t[24][1] = 17, t[24][2] = 13;
  s[24][0] = 40, s[24][1] = 32, s[24][2] = 28;
  t[25][0] = 9, t[25][1] = 16, t[25][2] = 17;
  s[25][0] = 41, s[25][1] = 28, s[25][2] = 38;

  t[26][0] = 10, t[26][1] = 13, t[26][2] = 18;
  s[26][0] = 42, s[26][1] = 29, s[26][2] = 33;
  t[27][0] = 10, t[27][1] = 18, t[27][2] = 14;
  s[27][0] = 43, s[27][1] = 34, s[27][2] = 29;

  t[28][0] = 11, t[28][1] = 14, t[28][2] = 19;
  s[28][0] = 44, s[28][1] = 30, s[28][2] = 35;
  t[29][0] = 11, t[29][1] = 19, t[29][2] = 15;
  s[29][0] = 45, s[29][1] = 36, s[29][2] = 30;

  t[30][0] = 12, t[30][1] = 15, t[30][2] = 20;
  s[30][0] = 46, s[30][1] = 31, s[30][2] = 37;
  t[31][0] = 12, t[31][1] = 20, t[31][2] = 16;
  s[31][0] = 47, s[31][1] = 39, s[31][2] = 31;

  /* right mid to right circle */
  t[32][0] = 13, t[32][1] = 17, t[32][2] = 18;
  s[32][0] = 48, s[32][1] = 42, s[32][2] = 40;
  t[33][0] = 14, t[33][1] = 18, t[33][2] = 19;
  s[33][0] = 52, s[33][1] = 44, s[33][2] = 43;
  t[34][0] = 15, t[34][1] = 19, t[34][2] = 20;
  s[34][0] = 53, s[34][1] = 46, s[34][2] = 45;
  t[35][0] = 16, t[35][1] = 20, t[35][2] = 17;
  s[35][0] = 55, s[35][1] = 41, s[35][2] = 47;

  /* right circle */
  t[36][0] = 17, t[36][1] = 21, t[36][2] = 18;
  s[36][0] = 49, s[36][1] = 48, s[36][2] = 50;

  t[37][0] = 18, t[37][1] = 21, t[37][2] = 19;
  s[37][0] = 51, s[37][1] = 52, s[37][2] = 49;

  t[38][0] = 19, t[38][1] = 21, t[38][2] = 20;
  s[38][0] = 54, s[38][1] = 53, s[38][2] = 51;

  t[39][0] = 17, t[39][1] = 20, t[39][2] = 21;
  s[39][0] = 54, s[39][1] = 50, s[39][2] = 55;

  mg->phi = cylinder_parametrization;
  mg->phidata = mg;

  return mg;
}

/* ------------------------------------------------------------
 Polygonal approximation
 ------------------------------------------------------------ */

psurface3d
build_from_macrosurface3d_surface3d(pcmacrosurface3d mg, uint split)
{
  psurface3d gr;
  uint      mvertices = mg->vertices;
  uint      medges = mg->edges;
  uint      mtriangles = mg->triangles;
  const     uint(*me)[2] = (const uint(*)[2]) mg->e;
  const     uint(*mt)[3] = (const uint(*)[3]) mg->t;
  const     uint(*ms)[3] = (const uint(*)[3]) mg->s;
  real(*x)[3];
  uint(*e)[2];
  uint(*t)[3];
  uint(*s)[3];
  uint     *vdone, *edone;
  uint     *vv, *ve, *ee;
  uint     *vloc, *evloc, *ehloc, *edloc;
  real      alpha, beta;
  uint      vertices, edges, triangles;
  uint      vcount, ecount, tcount;
  uint      tij, sij;
  uint      i, j, k;

  vertices = (mvertices + medges * (split - 1)
	      + mtriangles * (split - 1) * (split - 2) / 2);
  edges = (medges * split + mtriangles * 3 * split * (split - 1) / 2);
  triangles = mtriangles * split * split;

  gr = new_surface3d(vertices, edges, triangles);
  x = gr->x;
  e = gr->e;
  t = gr->t;
  s = gr->s;

  /* Flags indicating whether a macro vertex or edge has already
     been assigned vertices and edges */
  vdone = (uint *) allocmem((size_t) sizeof(uint) * mvertices);
  edone = (uint *) allocmem((size_t) sizeof(uint) * medges);

  /* Numbers of vertices assigned to macro vertices */
  vv = (uint *) allocmem((size_t) sizeof(uint) * mvertices);

  /* Numbers of vertices assigned to macro edges */
  ve = (uint *) allocmem((size_t) sizeof(uint) * medges);

  /* Numbers of edges assigned to macro edges */
  ee = (uint *) allocmem((size_t) sizeof(uint) * medges);

  /* Reset "done" flags for macro vertices and edges */
  for (i = 0; i < mvertices; i++)
    vdone[i] = 0;
  for (i = 0; i < medges; i++)
    edone[i] = 0;

  /* Create vertices and edges for macro vertices and edges */
  vcount = 0;
  ecount = 0;
  tcount = 0;
  for (i = 0; i < mtriangles; i++) {
    for (j = 0; j < 3; j++)
      if (!vdone[mt[i][j]]) {	/* Create vertex */
	tij = mt[i][j];

	vv[tij] = vcount;
	assert(vcount < vertices);

	switch (j) {
	case 0:
	  mg->phi(i, 0.0, 0.0, mg->phidata, x[vcount]);
	  break;
	case 1:
	  mg->phi(i, 1.0, 0.0, mg->phidata, x[vcount]);
	  break;
	case 2:
	  mg->phi(i, 0.0, 1.0, mg->phidata, x[vcount]);
	  break;
	default:
	  abort();
	}
	vcount++;

	vdone[tij] = 1;
      }

    for (j = 0; j < 3; j++)
      if (!edone[ms[i][j]]) {
	sij = ms[i][j];

	ve[sij] = vcount;
	assert(vcount < vertices);

	for (k = 1; k < split; k++) {	/* Create interior vertices of macro edge */
	  alpha = (me[sij][0] == mt[i][(j + 1) % 3] ? (real) k / split :
		   (real) (split - k) / split);
	  switch (j) {
	  case 0:
	    mg->phi(i, 1.0 - alpha, alpha, mg->phidata, x[vcount]);
	    break;
	  case 1:
	    mg->phi(i, 0.0, 1.0 - alpha, mg->phidata, x[vcount]);
	    break;
	  case 2:
	    mg->phi(i, alpha, 0.0, mg->phidata, x[vcount]);
	    break;
	  default:
	    abort();
	  }
	  vcount++;
	}

	if (split < 2) {	/* Create interior edges of macro edge */
	  ee[sij] = ecount;

	  assert(ecount < edges);
	  e[ecount][0] = vv[me[sij][0]];
	  e[ecount][1] = vv[me[sij][1]];
	  ecount++;
	}
	else {
	  ee[sij] = ecount;

	  assert(ecount < edges);
	  e[ecount][0] = vv[me[sij][0]];
	  e[ecount][1] = ve[sij];
	  ecount++;

	  for (k = 1; k < split - 1; k++) {
	    assert(ecount < edges);
	    e[ecount][0] = ve[sij] + k - 1;
	    e[ecount][1] = ve[sij] + k;
	    ecount++;
	  }

	  assert(ecount < edges);
	  e[ecount][0] = ve[sij] + split - 2;
	  e[ecount][1] = vv[me[sij][1]];
	  ecount++;
	}

	edone[sij] = 1;
      }
  }

  /* Prepare buffers for vertices and
     vertical, horizontal and diagonal edges
     (I know, they should by triangular) */
  vloc = (uint *) allocmem((size_t) sizeof(uint) * (split + 1) * (split + 1));
  evloc = (uint *) allocmem((size_t) sizeof(uint) * split * split);
  ehloc = (uint *) allocmem((size_t) sizeof(uint) * split * split);
  edloc = (uint *) allocmem((size_t) sizeof(uint) * split * split);

  /* Subdivide all macro triangles */
  for (i = 0; i < mtriangles; i++) {
    /* Find vertices for macro vertices */
    vloc[0] = vv[mt[i][0]];
    vloc[split] = vv[mt[i][1]];
    vloc[split * (split + 1)] = vv[mt[i][2]];

    /* Find vertices and edges on macro edge 0 */
    sij = ms[i][0];
    if (me[sij][0] == mt[i][1]) {	/* Macro edge oriented correctly */
      assert(me[sij][1] == mt[i][2]);

      for (j = 1; j < split; j++)
	vloc[(split - j) + j * (split + 1)] = ve[sij] + (j - 1);
      for (j = 0; j < split; j++)
	edloc[(split - 1 - j) + j * split] = ee[sij] + j;
    }
    else {			/* Macro edge reversed */
      assert(me[sij][0] == mt[i][2]);
      assert(me[sij][1] == mt[i][1]);

      for (j = 1; j < split; j++)
	vloc[j + (split - j) * (split + 1)] = ve[sij] + (j - 1);
      for (j = 0; j < split; j++)
	edloc[j + (split - 1 - j) * split] = ee[sij] + j;
    }

    /* Find vertices and edges on macro edge 1 */
    sij = ms[i][1];
    if (me[sij][0] == mt[i][2]) {	/* Macro edge oriented correctly */
      assert(me[sij][1] == mt[i][0]);

      for (j = 1; j < split; j++)
	vloc[(split - j) * (split + 1)] = ve[sij] + (j - 1);
      for (j = 0; j < split; j++)
	evloc[(split - 1 - j) * split] = ee[sij] + j;
    }
    else {			/* Macro edge reversed */
      assert(me[sij][0] == mt[i][0]);
      assert(me[sij][1] == mt[i][2]);

      for (j = 1; j < split; j++)
	vloc[j * (split + 1)] = ve[sij] + (j - 1);
      for (j = 0; j < split; j++)
	evloc[j * split] = ee[sij] + j;
    }

    /* Find vertices and edges on macro edge 2 */
    sij = ms[i][2];
    if (me[sij][0] == mt[i][0]) {	/* Macro edge oriented correctly */
      assert(me[sij][1] == mt[i][1]);

      for (j = 1; j < split; j++)
	vloc[j] = ve[sij] + (j - 1);
      for (j = 0; j < split; j++)
	ehloc[j] = ee[sij] + j;
    }
    else {			/* Macro edge reversed */
      assert(me[sij][0] == mt[i][1]);
      assert(me[sij][1] == mt[i][0]);

      for (j = 1; j < split; j++)
	vloc[split - j] = ve[sij] + (j - 1);
      for (j = 0; j < split; j++)
	ehloc[split - 1 - j] = ee[sij] + j;
    }

    /* Create interior vertices */
    for (j = 1; j < split; j++) {
      beta = (real) j / split;
      for (k = 1; k < split - j; k++) {
	alpha = (real) k / split;

	vloc[k + j * (split + 1)] = vcount;
	assert(vcount < vertices);
	mg->phi(i, alpha, beta, mg->phidata, x[vcount]);
	vcount++;
      }
    }

    /* Create interior diagonal edges */
    for (j = 0; j < split; j++)
      for (k = 0; k < split - 1 - j; k++) {
	edloc[k + j * split] = ecount;
	assert(ecount < edges);
	e[ecount][0] = vloc[(k + 1) + j * (split + 1)];
	e[ecount][1] = vloc[k + (j + 1) * (split + 1)];
	ecount++;
      }

    /* Create interior vertical edges */
    for (j = 0; j < split; j++)
      for (k = 0; k < split - 1 - j; k++) {
	evloc[(k + 1) + j * split] = ecount;
	assert(ecount < edges);
	e[ecount][0] = vloc[(k + 1) + j * (split + 1)];
	e[ecount][1] = vloc[(k + 1) + (j + 1) * (split + 1)];
	ecount++;
      }

    /* Create interior horizontal edges */
    for (j = 0; j < split; j++)
      for (k = 0; k < split - 1 - j; k++) {
	ehloc[k + (j + 1) * split] = ecount;
	assert(ecount < edges);
	e[ecount][0] = vloc[k + (j + 1) * (split + 1)];
	e[ecount][1] = vloc[(k + 1) + (j + 1) * (split + 1)];
	ecount++;
      }

    /* Create triangles */
    for (j = 0; j < split; j++) {
      for (k = 0; k < split - j - 1; k++) {
	assert(tcount < triangles);
	t[tcount][0] = vloc[k + j * (split + 1)];
	t[tcount][1] = vloc[(k + 1) + j * (split + 1)];
	t[tcount][2] = vloc[k + (j + 1) * (split + 1)];
	s[tcount][0] = edloc[k + j * split];
	s[tcount][1] = evloc[k + j * split];
	s[tcount][2] = ehloc[k + j * split];
	tcount++;

	assert(tcount < triangles);
	t[tcount][0] = vloc[(k + 1) + (j + 1) * (split + 1)];
	t[tcount][1] = vloc[k + (j + 1) * (split + 1)];
	t[tcount][2] = vloc[(k + 1) + j * (split + 1)];
	s[tcount][0] = edloc[k + j * split];
	s[tcount][1] = evloc[(k + 1) + j * split];
	s[tcount][2] = ehloc[k + (j + 1) * split];
	tcount++;
      }
      assert(tcount < triangles);
      t[tcount][0] = vloc[k + j * (split + 1)];
      t[tcount][1] = vloc[(k + 1) + j * (split + 1)];
      t[tcount][2] = vloc[k + (j + 1) * (split + 1)];
      s[tcount][0] = edloc[k + j * split];
      s[tcount][1] = evloc[k + j * split];
      s[tcount][2] = ehloc[k + j * split];
      tcount++;
    }
  }

  prepare_surface3d(gr);

  assert(vcount == vertices);
  assert(ecount == edges);
  assert(tcount == triangles);

  /* Clean up */
  freemem(edone);
  freemem(vdone);
  freemem(vv);
  freemem(ve);
  freemem(ee);
  freemem(vloc);
  freemem(evloc);
  freemem(ehloc);
  freemem(edloc);

  return gr;
}

/* ------------------------------------------------------------
 Interactive setup of a surface3d object
 ------------------------------------------------------------ */

psurface3d
build_interactive_surface3d()
{
  psurface3d gr;
  pmacrosurface3d mg;
  char      buf[100], *c, *res;
  FILE     *grid;
  uint      r;

  gr = 0;

  do {
    (void) printf("Name of grid? (file name or built-in)\n");
    buf[0] = '\0';
    res = fgets(buf, 100, stdin);
    assert(res != NULL);

    for (c = buf; *c != '\n' && *c != '\r' && *c != '\0'; c++);
    *c = '\0';

    grid = fopen(buf, "r");
    if (grid) {
      fclose(grid);
      (void) printf("Reading file \"%s\"\n", buf);
      gr = read_surface3d(buf);
    }
    else {
      switch (buf[0]) {
      case 's':
	if (sscanf(buf + 1, "%u", &r) == 1) {
	  (void) printf("Creating sphere geometry with %u subdivisions\n", r);
	  mg = new_sphere_macrosurface3d();
	  gr = build_from_macrosurface3d_surface3d(mg, r);
	  del_macrosurface3d(mg);
	}
	break;
      case 'c':
	switch (buf[1]) {
	case 'u':
	  if (sscanf(buf + 2, "%u", &r) == 1) {
	    (void) printf("Creating cube geometry with %u subdivisions\n", r);
	    mg = new_cube_macrosurface3d();
	    gr = build_from_macrosurface3d_surface3d(mg, r);
	    del_macrosurface3d(mg);
	  }
	  break;
	case 'y':
	  if (sscanf(buf + 2, "%u", &r) == 1) {
	    (void) printf("Creating cylinder geometry with %u subdivisions\n",
			  r);
	    mg = new_cylinder_macrosurface3d();
	    gr = build_from_macrosurface3d_surface3d(mg, r);
	    del_macrosurface3d(mg);
	  }
	  break;
	default:
	  gr = NULL;
	  break;
	}

	break;
      case 'p':
	if (sscanf(buf + 1, "%u", &r) == 1) {
	  (void)
	    printf
	    ("Creating parabolic mirror geometry with %u subdivisions\n", r);
	  mg = new_parabolic_mirror_macrosurface3d();
	  gr = build_from_macrosurface3d_surface3d(mg, r);
	  del_macrosurface3d(mg);
	}
	break;
      default:
	gr = NULL;
	break;
      }
    }
  } while (gr == NULL);

  prepare_surface3d(gr);

  (void) printf("%u vertices, %u edges, %u triangles\n", gr->vertices,
		gr->edges, gr->triangles);

  return gr;
}
