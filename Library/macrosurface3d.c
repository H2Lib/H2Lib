/* ------------------------------------------------------------
 This is the file "macrosurface3d.c" of the H2Lib package.
 All rights reserved, Steffen Boerm 2010
 ------------------------------------------------------------ */

#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "macrosurface3d.h"
#include "surface3d.h"
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

  norm = REAL_SQRT(REAL_SQR(xt[0]) + REAL_SQR(xt[1]) + REAL_SQR(xt[2]));
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
  char      buf[100], *c;
  FILE     *grid;
  uint      r;

  gr = 0;

  do {
    (void) printf("Name of grid? (file name or built-in)\n");
    buf[0] = '\0';
    (void) fgets(buf, 100, stdin);

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
      default:
	gr = 0;
	break;
      }
    }
  } while (gr == 0);

  prepare_surface3d(gr);

  (void) printf("%u vertices, %u edges, %u triangles\n", gr->vertices,
		gr->edges, gr->triangles);

  return gr;
}
