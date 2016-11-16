/* ------------------------------------------------------------
 This is the file "surface3d.c" of the H2Lib package.
 All rights reserved, Steffen Boerm 2010
 ------------------------------------------------------------ */

#include "surface3d.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#ifdef USE_ZLIB
#include "zlib.h"
#endif

#ifdef USE_NETCDF
#include "netcdf.h"
#endif

#include "basic.h"

/* ------------------------------------------------------------
 Constructor and destructor
 ------------------------------------------------------------ */

psurface3d
new_surface3d(uint vertices, uint edges, uint triangles)
{
  psurface3d gr;

  gr = (psurface3d) allocmem(sizeof(surface3d));
  gr->x = (real(*)[3]) allocmem((size_t) sizeof(real[3]) * vertices);
  gr->e = (uint(*)[2]) allocmem((size_t) sizeof(uint[2]) * edges);
  gr->t = (uint(*)[3]) allocmem((size_t) sizeof(uint[3]) * triangles);
  gr->s = (uint(*)[3]) allocmem((size_t) sizeof(uint[3]) * triangles);
  gr->n = (real(*)[3]) allocmem((size_t) sizeof(real[3]) * triangles);
  gr->g = (real *) allocmem((size_t) sizeof(real) * triangles);

  gr->vertices = vertices;
  gr->edges = edges;
  gr->triangles = triangles;

  gr->hmin = 1e30;
  gr->hmax = 0.0;

  return gr;
}

void
prepare_surface3d(psurface3d gr)
{
  const     real(*x)[3] = (const real(*)[3]) gr->x;
  const     uint(*t)[3] = (const uint(*)[3]) gr->t;
  real(*n)[3] = gr->n;
  real     *g = gr->g;
  real      dx[3], dy[3], dz[3];
  real      norm, height, hmin, hmax, a, b, c;
  uint      triangles = gr->triangles;
  uint      i;

  hmin = 1e30;
  hmax = 0.0;

  for (i = 0; i < triangles; i++) {
    dx[0] = x[t[i][1]][0] - x[t[i][0]][0];
    dx[1] = x[t[i][1]][1] - x[t[i][0]][1];
    dx[2] = x[t[i][1]][2] - x[t[i][0]][2];

    dy[0] = x[t[i][2]][0] - x[t[i][0]][0];
    dy[1] = x[t[i][2]][1] - x[t[i][0]][1];
    dy[2] = x[t[i][2]][2] - x[t[i][0]][2];

    dz[0] = x[t[i][2]][0] - x[t[i][1]][0];
    dz[1] = x[t[i][2]][1] - x[t[i][1]][1];
    dz[2] = x[t[i][2]][2] - x[t[i][1]][2];

    n[i][0] = dx[1] * dy[2] - dx[2] * dy[1];
    n[i][1] = dx[2] * dy[0] - dx[0] * dy[2];
    n[i][2] = dx[0] * dy[1] - dx[1] * dy[0];
    norm =
      REAL_SQRT(REAL_SQR(n[i][0]) + REAL_SQR(n[i][1]) + REAL_SQR(n[i][2]));
    g[i] = norm;
    n[i][0] /= norm;
    n[i][1] /= norm;
    n[i][2] /= norm;

    a = REAL_SQRT(REAL_SQR(dx[0]) + REAL_SQR(dx[1]) + REAL_SQR(dx[2]));
    b = REAL_SQRT(REAL_SQR(dy[0]) + REAL_SQR(dy[1]) + REAL_SQR(dy[2]));
    c = REAL_SQRT(REAL_SQR(dz[0]) + REAL_SQR(dz[1]) + REAL_SQR(dz[2]));

    height = 2.0 * REAL_MIN3(a, b, c);

    a *= a;
    b *= b;
    c *= c;

    height = REAL_SQRT(2 * (a * b + b * c + a * c) - (a * a + b * b + c * c))
      / height;

    hmin = REAL_MIN(hmin, height);
    hmax = REAL_MAX(hmax, height);
  }

  gr->hmin = hmin;
  gr->hmax = hmax;
}

void
del_surface3d(psurface3d gr)
{
  freemem(gr->g);
  freemem(gr->n);
  freemem(gr->s);
  freemem(gr->t);
  freemem(gr->e);
  freemem(gr->x);
  freemem(gr);
}

/* ------------------------------------------------------------
 * Statistics
 * ------------------------------------------------------------ */

void
getproperties_surface3d(pcsurface3d gr, preal hmin, preal hmax,
			preal anglemin, preal angleedge)
{
  uint      edges = gr->edges;
  uint      triangles = gr->triangles;
  const     real(*x)[3] = (const real(*)[3]) gr->x;
  const     uint(*t)[3] = (const uint(*)[3]) gr->t;
  const     uint(*s)[3] = (const uint(*)[3]) gr->s;
  const     real(*n)[3] = (const real(*)[3]) gr->n;
  const real *g = (const real *) gr->g;
  uint     *neighbour;
  real      dt[3], ds[3], h[3], nn, angle;
  real      h_min, h_max, angle_min, angle_edge;
  uint      i, j;

  /* Just to avoid warnings */
  h_min = h_max = angle_min = 0.0;

  for (i = 0; i < triangles; i++) {
    dt[0] = x[t[i][1]][0] - x[t[i][0]][0];
    dt[1] = x[t[i][1]][1] - x[t[i][0]][1];
    dt[2] = x[t[i][1]][2] - x[t[i][0]][2];

    ds[0] = x[t[i][2]][0] - x[t[i][0]][0];
    ds[1] = x[t[i][2]][1] - x[t[i][0]][1];
    ds[2] = x[t[i][2]][2] - x[t[i][0]][2];

    h[0] = (REAL_SQR(dt[0]) + REAL_SQR(dt[1]) + REAL_SQR(dt[2]));
    h[1] = (REAL_SQR(ds[0]) + REAL_SQR(ds[1]) + REAL_SQR(ds[2]));
    h[2] = (REAL_SQR(dt[0] - ds[0]) + REAL_SQR(dt[1] - ds[1])
	    + REAL_SQR(dt[2] - ds[2]));

    if (i == 0 || h[0] < h_min)
      h_min = h[0];

    if (i == 0 || h[0] > h_max)
      h_max = h[0];

    if (h[1] < h_min)
      h_min = h[1];

    if (h[1] > h_max)
      h_max = h[1];

    if (h[2] < h_min)
      h_min = h[2];

    if (h[2] > h_max)
      h_max = h[2];

    angle = g[i] / REAL_SQRT(h[0] * h[1]);

    if (i == 0 || angle < angle_min)
      angle_min = angle;
  }

  if (hmin)
    *hmin = REAL_SQRT(h_min);

  if (hmax)
    *hmax = REAL_SQRT(h_max);

  if (anglemin)
    *anglemin = REAL_ASIN(angle_min);

  if (angleedge) {
    neighbour = (uint *) allocmem(sizeof(uint) * edges);

    for (i = 0; i < edges; i++)
      neighbour[i] = triangles;

    angle_edge = 1.0;

    for (i = 0; i < triangles; i++) {
      if (neighbour[s[i][0]] < triangles) {
	j = neighbour[s[i][0]];
	nn = (n[i][0] * n[j][0] + n[i][1] * n[j][1] + n[i][2] * n[j][2]);
	if (nn < angle_edge)
	  angle_edge = nn;
      }
      else
	neighbour[s[i][0]] = i;

      if (neighbour[s[i][1]] < triangles) {
	j = neighbour[s[i][1]];
	nn = (n[i][0] * n[j][0] + n[i][1] * n[j][1] + n[i][2] * n[j][2]);
	if (nn < angle_edge)
	  angle_edge = nn;
      }
      else
	neighbour[s[i][1]] = i;

      if (neighbour[s[i][2]] < triangles) {
	j = neighbour[s[i][2]];
	nn = (n[i][0] * n[j][0] + n[i][1] * n[j][1] + n[i][2] * n[j][2]);
	if (nn < angle_edge)
	  angle_edge = nn;
      }
      else
	neighbour[s[i][2]] = i;
    }

    *angleedge = REAL_ACOS(angle_edge);

    freemem(neighbour);
  }
}

/* ------------------------------------------------------------
 * Debugging
 * ------------------------------------------------------------ */

void
print_surface3d(pcsurface3d gr)
{
  uint      vertices = gr->vertices;
  uint      edges = gr->edges;
  uint      triangles = gr->triangles;
  const     real(*x)[3] = (const real(*)[3]) gr->x;
  const     uint(*e)[2] = (const uint(*)[2]) gr->e;
  const     uint(*t)[3] = (const uint(*)[3]) gr->t;
  const     uint(*s)[3] = (const uint(*)[3]) gr->s;
  const     real(*n)[3] = (const real(*)[3]) gr->n;
  const real *g = (const real *) gr->g;
  uint      i;

  (void) printf("surface3d(%u,%u,%u)\n", vertices, edges, triangles);
  for (i = 0; i < vertices; i++)
    (void) printf(" %d: (% .5e % .5e % .5e)\n", i, x[i][0], x[i][1], x[i][2]);
  for (i = 0; i < edges; i++)
    (void) printf(" %d: (%u %u)\n", i, e[i][0], e[i][1]);
  for (i = 0; i < triangles; i++)
    (void) printf(" %d: (%u %u %u  %u %u %u  % .5e % .5e % .5e  % .5e)\n", i,
		  t[i][0], t[i][1], t[i][2], s[i][0], s[i][1], s[i][2],
		  n[i][0], n[i][1], n[i][2], g[i]);
}

uint
check_surface3d(pcsurface3d gr)
{
  uint      vertices = gr->vertices;
  uint      edges = gr->edges;
  uint      triangles = gr->triangles;
  const     uint(*e)[2] = (const uint(*)[2]) gr->e;
  const     uint(*t)[3] = (const uint(*)[3]) gr->t;
  const     uint(*s)[3] = (const uint(*)[3]) gr->s;
  uint      i, j, problems;

  problems = 0;

  for (i = 0; i < edges; i++)
    for (j = 0; j < 2; j++)
      if (e[i][j] >= vertices) {
	(void) printf(" Vertex %u in edge %u out of bounds\n", j, i);
	problems++;
      }

  for (i = 0; i < triangles; i++)
    for (j = 0; j < 3; j++) {
      if (t[i][j] >= vertices) {
	(void) printf(" Vertex %u in triangle %u out of bounds\n", j, i);
	problems++;
      }

      if (s[i][j] >= edges) {
	(void) printf(" Edge %u in triangle %u out of bounds\n", j, i);
	problems++;
      }
      else {
	if (e[s[i][j]][0] == t[i][(j + 1) % 3]) {
	  if (e[s[i][j]][1] != t[i][(j + 2) % 3]) {
	    (void) printf(" Mismatched edge %u in triangle %u\n", j, i);
	    problems++;
	  }
	}
	else {
	  if (e[s[i][j]][0] != t[i][(j + 2) % 3] || e[s[i][j]][1]
	      != t[i][(j + 1) % 3]) {
	    (void) printf(" *Mismatched edge %u in triangle %u\n", j, i);
	    problems++;
	  }
	}
      }
    }

  return problems;
}

bool
isclosed_surface3d(pcsurface3d gr)
{
  uint      edges = gr->edges;
  uint      triangles = gr->triangles;
  const     uint(*s)[3] = (const uint(*)[3]) gr->s;
  int      *visits;
  uint      i, j, problems;

  visits = (int *) allocmem((size_t) sizeof(int) * edges);

  for (i = 0; i < edges; i++)
    visits[i] = 0;

  for (i = 0; i < triangles; i++)
    for (j = 0; j < 3; j++)
      visits[s[i][j]]++;

  problems = 0;
  for (i = 0; i < edges; i++)
    if (visits[i] != 2)
      problems++;

  freemem(visits);

  return (problems == 0);
}

bool
isoriented_surface3d(pcsurface3d gr)
{
  uint      edges = gr->edges;
  uint      triangles = gr->triangles;
  const     uint(*e)[2] = (const uint(*)[2]) gr->e;
  const     uint(*t)[3] = (const uint(*)[3]) gr->t;
  const     uint(*s)[3] = (const uint(*)[3]) gr->s;
  int      *visits, *svisits;
  uint      i, j, k, problems;

  visits = (int *) allocmem((size_t) sizeof(int) * edges);
  svisits = (int *) allocmem((size_t) sizeof(int) * edges);

  for (i = 0; i < edges; i++)
    visits[i] = svisits[i] = 0;

  for (i = 0; i < triangles; i++)
    for (j = 0; j < 3; j++) {
      k = s[i][j];
      visits[k]++;
      if (e[k][0] == t[i][(j + 1) % 3]) {
	assert(e[k][1] == t[i][(j + 2) % 3]);
	svisits[k]++;
      }
      else {
	assert(e[k][0] == t[i][(j + 2) % 3]);
	assert(e[k][1] == t[i][(j + 1) % 3]);
	svisits[k]--;
      }
    }

  problems = 0;
  for (i = 0; i < edges; i++)
    if (visits[i] > 1 && (svisits[i] < -1 || svisits[i] > 1))
      problems++;

  freemem(svisits);
  freemem(visits);

  return (problems == 0);
}

void
scale_surface3d(psurface3d gr, real * a, real * b)
{
  uint      vertices = gr->vertices;
  real(*x)[3] = gr->x;

  uint      i, j, d;
  real      min[3], max[3], diam[3], bdiam[3];
  real      alpha, beta;

  for (j = 0; j < 3; ++j) {
    min[j] = 1e30;
    max[j] = -1e30;
    for (i = 0; i < vertices; ++i) {
      min[j] = (x[i][j] < min[j]) ? x[i][j] : min[j];
      max[j] = (x[i][j] > max[j]) ? x[i][j] : max[j];
    }
  }

  for (j = 0; j < 3; ++j) {
    diam[j] = max[j] - min[j];
    bdiam[j] = b[j] - a[j];
    diam[j] = diam[j] / bdiam[j];
  }

  d = (diam[0] > diam[1]) ? (diam[0] > diam[2] ? 0 : 2) :
    (diam[1] > diam[2] ? 1 : 2);
  beta = (b[d] - a[d]) / (max[d] - min[d]);
  alpha = a[d] - beta * min[d];

  for (j = 0; j < 3; ++j) {
    for (i = 0; i < vertices; ++i) {
      x[i][j] = alpha + beta * x[i][j];
    }
  }

  prepare_surface3d(gr);
}

void
translate_surface3d(psurface3d gr, real * t)
{
  uint      vertices = gr->vertices;

  uint      i;

  for (i = 0; i < vertices; ++i) {
    gr->x[i][0] += t[0];
    gr->x[i][1] += t[1];
    gr->x[i][2] += t[2];
  }

  prepare_surface3d(gr);
}

psurface3d
merge_surface3d(pcsurface3d gr1, pcsurface3d gr2)
{
  uint      i;

  psurface3d gr;

  gr = new_surface3d(gr1->vertices + gr2->vertices, gr1->edges + gr2->edges,
		     gr1->triangles + gr2->triangles);

  for (i = 0; i < gr1->vertices; ++i) {
    gr->x[i][0] = gr1->x[i][0];
    gr->x[i][1] = gr1->x[i][1];
    gr->x[i][2] = gr1->x[i][2];
  }

  for (i = 0; i < gr2->vertices; ++i) {
    gr->x[i + gr1->vertices][0] = gr2->x[i][0];
    gr->x[i + gr1->vertices][1] = gr2->x[i][1];
    gr->x[i + gr1->vertices][2] = gr2->x[i][2];
  }

  for (i = 0; i < gr1->edges; ++i) {
    gr->e[i][0] = gr1->e[i][0];
    gr->e[i][1] = gr1->e[i][1];
    gr->e[i][2] = gr1->e[i][2];
  }

  for (i = 0; i < gr2->edges; ++i) {
    gr->e[i + gr1->edges][0] = gr2->e[i][0] + gr1->vertices;
    gr->e[i + gr1->edges][1] = gr2->e[i][1] + gr1->vertices;
    gr->e[i + gr1->edges][2] = gr2->e[i][2] + gr1->vertices;
  }

  for (i = 0; i < gr1->triangles; ++i) {
    gr->t[i][0] = gr1->t[i][0];
    gr->t[i][1] = gr1->t[i][1];
    gr->t[i][2] = gr1->t[i][2];
    gr->s[i][0] = gr1->s[i][0];
    gr->s[i][1] = gr1->s[i][1];
    gr->s[i][2] = gr1->s[i][2];
    gr->n[i][0] = gr1->n[i][0];
    gr->n[i][1] = gr1->n[i][1];
    gr->n[i][2] = gr1->n[i][2];
    gr->g[i] = gr1->g[i];
  }

  for (i = 0; i < gr2->triangles; ++i) {
    gr->t[i + gr1->triangles][0] = gr2->t[i][0] + gr1->vertices;
    gr->t[i + gr1->triangles][1] = gr2->t[i][1] + gr1->vertices;
    gr->t[i + gr1->triangles][2] = gr2->t[i][2] + gr1->vertices;
    gr->s[i + gr1->triangles][0] = gr2->s[i][0] + gr1->edges;
    gr->s[i + gr1->triangles][1] = gr2->s[i][1] + gr1->edges;
    gr->s[i + gr1->triangles][2] = gr2->s[i][2] + gr1->edges;
    gr->n[i + gr1->triangles][0] = gr2->n[i][0];
    gr->n[i + gr1->triangles][1] = gr2->n[i][1];
    gr->n[i + gr1->triangles][2] = gr2->n[i][2];
    gr->g[i + gr1->triangles] = gr2->g[i];
  }

  return gr;
}

/* ------------------------------------------------------------
 File I/O
 ------------------------------------------------------------ */

void
write_surface3d(pcsurface3d gr, const char *filename)
{
  FILE     *out;
  uint      vertices = gr->vertices;
  uint      edges = gr->edges;
  uint      triangles = gr->triangles;
  const     real(*x)[3] = (const real(*)[3]) gr->x;
  const     uint(*e)[2] = (const uint(*)[2]) gr->e;
  const     uint(*t)[3] = (const uint(*)[3]) gr->t;
  const     uint(*s)[3] = (const uint(*)[3]) gr->s;
  uint      i;

  out = fopen(filename, "w");
  if (out == 0) {
    (void) fprintf(stderr, "Could not open file \"%s\" for writing\n",
		   filename);
    return;
  }

  (void) fprintf(out, "%u %u %u\n", vertices, edges, triangles);

  for (i = 0; i < vertices; i++)
    (void) fprintf(out, "% .8g % .8g % .8g\n", x[i][0], x[i][1], x[i][2]);

  for (i = 0; i < edges; i++)
    (void) fprintf(out, "%u %u\n", e[i][0], e[i][1]);

  for (i = 0; i < triangles; i++)
    (void) fprintf(out, "%u %u %u  %u %u %u\n", t[i][0], t[i][1], t[i][2],
		   s[i][0], s[i][1], s[i][2]);

  (void) fclose(out);
}

#ifdef USE_ZLIB
static char *
readline(char *buf, int bufsize, gzFile in, uint * ln)
#else
static char *
readline(char *buf, int bufsize, FILE * in, uint * ln)
#endif
{
  char     *line;

  (void) bufsize;

  do {
#ifdef USE_ZLIB
    line = gzgets(in, buf, 80);
#else
    line = fgets(buf, 80, in);
#endif
    (*ln)++;
  } while (line && line[0] == '#');

  return line;
}

psurface3d
read_surface3d(const char *filename)
{
  psurface3d gr;
#ifdef USE_ZLIB
  gzFile    in;
#else
  FILE     *in;
#endif
  uint      vertices, edges, triangles;
  real(*x)[3];
  uint(*e)[2];
  uint(*t)[3];
  uint(*s)[3];
  char      buf[80];
  char     *line;
  uint      i, ln;

#ifdef USE_ZLIB
  in = gzopen(filename, "rb");
#else
  in = fopen(filename, "r");
#endif
  if (in == 0) {
    (void) fprintf(stderr, "Could not open file \"%s\" for reading\n",
		   filename);
    return 0;
  }

  ln = 0;
  line = readline(buf, 80, in, &ln);

  if (line == 0 || sscanf(line, "%u %u %u", &vertices, &edges, &triangles)
      != 3) {
    (void) fprintf(stderr, "Could not read first line of file \"%s\"\n",
		   filename);

#ifdef USE_ZLIB
    gzclose(in);
#else
    fclose(in);
#endif
    return 0;
  }

  gr = new_surface3d(vertices, edges, triangles);
  x = gr->x;
  e = gr->e;
  t = gr->t;
  s = gr->s;

  for (i = 0; i < vertices; i++) {
    line = readline(buf, 80, in, &ln);

    if (line == 0 || sscanf(line, "%" SCANF_PREFIX "f %" SCANF_PREFIX "f %"
			    SCANF_PREFIX "f", x[i], x[i] + 1, x[i] + 2)
	!= 3) {
      (void) fprintf(stderr,
		     "Could not read vertex %u in line %u of file \"%s\"\n",
		     i, ln, filename);
      del_surface3d(gr);

#ifdef USE_ZLIB
      gzclose(in);
#else
      fclose(in);
#endif
    }
  }

  for (i = 0; i < edges; i++) {
    line = readline(buf, 80, in, &ln);

    if (line == 0 || sscanf(line, "%u %u", e[i], e[i] + 1) != 2) {
      (void) fprintf(stderr,
		     "Could not read edge %u in line %u of file \"%s\"\n", i,
		     ln, filename);
      del_surface3d(gr);

#ifdef USE_ZLIB
      gzclose(in);
#else
      fclose(in);
#endif
    }
  }

  for (i = 0; i < triangles; i++) {
    line = readline(buf, 80, in, &ln);

    if (line == 0 || sscanf(line, "%u %u %u  %u %u %u", t[i], t[i] + 1,
			    t[i] + 2, s[i], s[i] + 1, s[i] + 2)
	!= 6) {
      (void) fprintf(stderr,
		     "Could not read triangle %u in line %u of file \"%s\"\n",
		     i, ln, filename);
      del_surface3d(gr);

#ifdef USE_ZLIB
      gzclose(in);
#else
      fclose(in);
#endif
    }
  }

#ifdef USE_ZLIB
  gzclose(in);
#else
  fclose(in);
#endif

  prepare_surface3d(gr);

  return gr;
}

#ifdef USE_NETCDF
static void
nc_handle_error(int res)
{
  if (res != NC_NOERR)
    (void) fprintf(stderr, "NetCDF error %d, \"%s\"\n",
		   res, nc_strerror(res));
}
#endif

void
write_nc_surface3d(pcsurface3d gr, const char *filename)
{
#ifdef USE_NETCDF
  const     real(*x)[3] = (const real(*)[3]) gr->x;
  const     uint(*e)[2] = (const uint(*)[2]) gr->e;
  const     uint(*t)[3] = (const uint(*)[3]) gr->t;
  const     uint(*s)[3] = (const uint(*)[3]) gr->s;
  uint      vertices = gr->vertices;
  uint      edges = gr->edges;
  uint      triangles = gr->triangles;
  int       res, nc_out;
  int       x_id, e_id, t_id, s_id;
  int       vertices_id, edges_id, triangles_id;
  int       n2_id, n3_id;
  int       dimid[2];

  /* Create CDF file */
  res = nc_create(filename, NC_NETCDF4, &nc_out);
  nc_handle_error(res);

  /* Define dimensions */
  res = nc_def_dim(nc_out, "vertices", vertices, &vertices_id);
  nc_handle_error(res);
  res = nc_def_dim(nc_out, "edges", edges, &edges_id);
  nc_handle_error(res);
  res = nc_def_dim(nc_out, "triangles", triangles, &triangles_id);
  nc_handle_error(res);
  res = nc_def_dim(nc_out, "two", 2, &n2_id);
  nc_handle_error(res);
  res = nc_def_dim(nc_out, "three", 3, &n3_id);
  nc_handle_error(res);

  /* Define x variable */
  dimid[0] = vertices_id;
  dimid[1] = n3_id;
  res = nc_def_var(nc_out, "x", NC_DOUBLE, 2, dimid, &x_id);
  nc_handle_error(res);
  res = nc_def_var_deflate(nc_out, x_id, 0, 1, NC_DEFLATE_LEVEL);
  nc_handle_error(res);

  /* Define e variable */
  dimid[0] = edges_id;
  dimid[1] = n2_id;
  res = nc_def_var(nc_out, "e", NC_UINT, 2, dimid, &e_id);
  nc_handle_error(res);
  res = nc_def_var_deflate(nc_out, e_id, 0, 1, NC_DEFLATE_LEVEL);
  nc_handle_error(res);

  /* Define t variable */
  dimid[0] = triangles_id;
  dimid[1] = n3_id;
  res = nc_def_var(nc_out, "t", NC_UINT, 2, dimid, &t_id);
  nc_handle_error(res);
  res = nc_def_var_deflate(nc_out, t_id, 0, 1, NC_DEFLATE_LEVEL);
  nc_handle_error(res);

  /* Define s variable */
  dimid[0] = triangles_id;
  dimid[1] = n3_id;
  res = nc_def_var(nc_out, "s", NC_UINT, 2, dimid, &s_id);
  nc_handle_error(res);
  res = nc_def_var_deflate(nc_out, s_id, 0, 1, NC_DEFLATE_LEVEL);
  nc_handle_error(res);

  /* Leave define mode */
  res = nc_enddef(nc_out);
  nc_handle_error(res);

  /* Write variables */
  res = nc_put_var(nc_out, x_id, x);
  nc_handle_error(res);
  res = nc_put_var(nc_out, e_id, e);
  nc_handle_error(res);
  res = nc_put_var(nc_out, t_id, t);
  nc_handle_error(res);
  res = nc_put_var(nc_out, s_id, s);
  nc_handle_error(res);

  /* Close file */
  res = nc_close(nc_out);
  nc_handle_error(res);

#else
  (void) gr;
  (void) filename;

  (void) printf("Sorry, no NetCDF support.\n");
#endif
}

psurface3d
read_nc_surface3d(const char *filename)
{
#ifdef USE_NETCDF
  psurface3d gr;
  real(*x)[3];
  uint(*e)[2];
  uint(*t)[3];
  uint(*s)[3];
  uint      vertices;
  uint      edges;
  uint      triangles;
  int       res, nc_in;
  int       x_id, e_id, t_id, s_id;
  int       vertices_id, edges_id, triangles_id;
  size_t    dim;

  /* Open CDF file */
  res = nc_open(filename, NC_NOWRITE, &nc_in);
  nc_handle_error(res);

  /* Obtain dimensions */
  res = nc_inq_dimid(nc_in, "vertices", &vertices_id);
  nc_handle_error(res);
  res = nc_inq_dimid(nc_in, "edges", &edges_id);
  nc_handle_error(res);
  res = nc_inq_dimid(nc_in, "triangles", &triangles_id);
  nc_handle_error(res);

  /* Get values of dimensions */
  res = nc_inq_dimlen(nc_in, vertices_id, &dim);
  nc_handle_error(res);
  vertices = dim;
  res = nc_inq_dimlen(nc_in, edges_id, &dim);
  nc_handle_error(res);
  edges = dim;
  res = nc_inq_dimlen(nc_in, triangles_id, &dim);
  nc_handle_error(res);
  triangles = dim;

  /* Create surface3d object */
  gr = new_surface3d(vertices, edges, triangles);
  x = gr->x;
  e = gr->e;
  t = gr->t;
  s = gr->s;

  /* Obtain variables */
  res = nc_inq_varid(nc_in, "x", &x_id);
  nc_handle_error(res);
  res = nc_inq_varid(nc_in, "e", &e_id);
  nc_handle_error(res);
  res = nc_inq_varid(nc_in, "t", &t_id);
  nc_handle_error(res);
  res = nc_inq_varid(nc_in, "s", &s_id);
  nc_handle_error(res);

  /* Read variables */
  res = nc_get_var(nc_in, x_id, x);
  res = nc_get_var(nc_in, e_id, e);
  res = nc_get_var(nc_in, t_id, t);
  res = nc_get_var(nc_in, s_id, s);

  /* Close file */
  res = nc_close(nc_in);
  nc_handle_error(res);

  return gr;

#else
  (void) filename;

  (void) printf("Sorry, no NetCDF support.\n");
  return 0;
#endif
}

struct edge_list {
  uint      x;
  struct edge_list *next;
};

static struct edge_list *
prepend_edge(struct edge_list *list, uint vert)
{
  struct edge_list *list2, *new;

  assert(list != NULL);

  list2 = list;

  while (list2->next != NULL) {
    if (list2->x == vert) {
      return list;
    }
    list2 = list2->next;
  }
  new = (struct edge_list *) allocmem(sizeof(struct edge_list));
  new->x = vert;
  new->next = list;

  return new;
}

static void
mark_edge(struct edge_list *list, uint vert)
{
  while (list->next != NULL) {
    if (list->x == vert) {
      list->x = (uint) - 1;
      return;
    }
    list = list->next;
  }
}

static void
prepare_edges(uint vertices, uint(*t)[3], uint triangles,
	      uint(*e)[2], uint edges)
{
  struct edge_list *(*list);
  struct edge_list *l;
  uint      i, j, k;

  list =
    (struct edge_list **) allocmem(vertices * sizeof(struct edge_list *));

  for (i = 0; i < vertices; ++i) {
    list[i] = (struct edge_list *) allocmem(sizeof(struct edge_list));
    list[i]->x = (uint) - 1;
    list[i]->next = NULL;
  }

  for (i = 0; i < triangles; ++i) {
    assert(t[i][0] < vertices);
    assert(t[i][1] < vertices);
    assert(t[i][2] < vertices);

    list[t[i][0]] = prepend_edge(list[t[i][0]], t[i][1]);
    list[t[i][0]] = prepend_edge(list[t[i][0]], t[i][2]);

    list[t[i][1]] = prepend_edge(list[t[i][1]], t[i][0]);
    list[t[i][1]] = prepend_edge(list[t[i][1]], t[i][2]);

    list[t[i][2]] = prepend_edge(list[t[i][2]], t[i][0]);
    list[t[i][2]] = prepend_edge(list[t[i][2]], t[i][1]);
  }

  //  *edges = 0;
  //  for (i = 0; i < vertices; ++i) {
  //    l = list[i];
  //    while (l->next != NULL) {
  //      l = l->next;
  //      (*edges)++;
  //    }
  //  }
  //  *edges = *edges / 2;

  //  *e = (uint (*)[2]) allocmem(*edges * sizeof(uint[2]));

  i = 0;
  for (j = 0; j < vertices; ++j) {
    l = list[j];
    while (l->next != NULL) {
      k = l->x;
      if (k != (uint) - 1) {
	//        (*e)[i][0] = j;
	//        (*e)[i][1] = k;
	e[i][0] = j;
	e[i][1] = k;
	mark_edge(list[j], k);
	mark_edge(list[k], j);
	i++;
      }
      l = l->next;
    }
  }
  assert(i == edges);

  for (i = 0; i < vertices; ++i) {
    freemem(list[i]);
  }
  freemem(list);
}

static void
prepare_arrays_s(uint(*t)[3], uint(*e)[2], uint triangles,
		 uint edges, uint(*s)[3])
{
  uint      i, j, k;

  for (i = 0; i < triangles; ++i) {
    for (k = 0; k < 3; ++k) {
      for (j = 0; j < edges; ++j) {
	if (e[j][0] == t[i][(k + 1) % 3] && e[j][1] == t[i][(k + 2) % 3]) {
	  s[i][k] = j;
	  break;
	}
	if (e[j][1] == t[i][(k + 1) % 3] && e[j][0] == t[i][(k + 2) % 3]) {
	  s[i][k] = j;
	  break;
	}
      }
    }
  }
}

psurface3d
read_netgen_surface3d(const char *filename)
{
#ifdef USE_ZLIB
  gzFile    in;
#else
  FILE     *in;
#endif
  uint      vertices, edges, triangles;
  real(*x)[3];
  uint(*e)[2];
  uint(*t)[3];
  uint(*s)[3];
  uint      tmp[5];
  char      buf[255];
  char     *line;
  uint      i, ln;

  psurface3d gr;

  x = NULL;
  e = NULL;
  t = NULL;
  s = NULL;

#ifdef USE_ZLIB
  in = gzopen(filename, "rb");
#else
  in = fopen(filename, "r");
#endif
  if (in == 0) {
    (void) fprintf(stderr, "Could not open file \"%s\" for reading\n",
		   filename);
    return NULL;
  }

  ln = 0;
  line = readline(buf, 255, in, &ln);
  while (line) {
    /* retrieve number of triangles */
    if (strcmp(line, "surfaceelements\n") == 0 || strcmp(line,
							 "surfaceelementsgi\n")
	== 0) {
      line = readline(buf, 255, in, &ln);
      if (sscanf(line, "%u", &triangles)) {
	t = (uint(*)[3]) allocmem(triangles * 3 * sizeof(uint));
	s = (uint(*)[3]) allocmem(triangles * 3 * sizeof(uint));
	line = readline(buf, 255, in, &ln);
	i = 0;
	while (line != NULL && strcmp(line, "\n") != 0) {
	  sscanf(line, "%u %u %u %u %u %u %u %u", &tmp[0], &tmp[1], &tmp[2],
		 &tmp[3], &tmp[4], &t[i][0], &t[i][1], &t[i][2]);
	  t[i][0]--;
	  t[i][1]--;
	  t[i][2]--;
	  if (tmp[3] == 1 && tmp[2] == 0) {
	    tmp[3] = t[i][0];
	    t[i][0] = t[i][2];
	    t[i][2] = tmp[3];
	  }
	  i++;
	  line = readline(buf, 255, in, &ln);
	}
	assert(i == triangles);
      }
      else {
	(void) fprintf(stderr, "Unknown file format.\n");
      }
    }

    /* retrieve number of vertices */
    if (strcmp(line, "points\n") == 0) {
      line = readline(buf, 255, in, &ln);
      if (sscanf(line, "%u", &vertices)) {
	x = (real(*)[3]) allocmem(vertices * 3 * sizeof(real));
	line = readline(buf, 255, in, &ln);
	i = 0;
	while (line != NULL && strcmp(line, "\n") != 0) {
	  sscanf(line,
		 "%" SCANF_PREFIX "f %" SCANF_PREFIX "f %" SCANF_PREFIX "f",
		 &x[i][0], &x[i][1], &x[i][2]);
	  i++;
	  line = readline(buf, 255, in, &ln);
	}
	assert(i == vertices);
      }
      else {
	(void) fprintf(stderr, "Unknown file format.\n");
      }
    }
    line = readline(buf, 255, in, &ln);
  }
#ifdef USE_ZLIB
  gzclose(in);
#else
  fclose(in);
#endif

  edges = 1.5 * triangles;
  e = (uint(*)[2]) allocuint(edges);

  prepare_edges(vertices, t, triangles, e, edges);
  prepare_arrays_s(t, e, triangles, edges, s);

  printf("geometry with %u vertices, %u edges, %u triangles read.\n",
	 vertices, edges, triangles);

  gr = new_surface3d(vertices, edges, triangles);

  freemem(gr->t);
  gr->t = t;
  freemem(gr->x);
  gr->x = x;
  freemem(gr->e);
  gr->e = e;
  freemem(gr->s);
  gr->s = s;

  return gr;
}

psurface3d
read_gmsh_surface3d(const char *filename)
{
#ifdef USE_ZLIB
  gzFile    file;
#else
  FILE     *file;
#endif
  char      buf[255];
  uint      line;
  uint      vertices, edges, triangles;
  uint      i, tmp, tmp2, tmp3, tmp4, tmp5, tmp6;
  int       c;

  psurface3d gr;

#ifdef USE_ZLIB
  file = gzopen(filename, "rb");
#else
  file = fopen(filename, "r");
#endif

  line = 0;

  /****************************************************
   * Find number of vertices ...
   ****************************************************/

  while (strcmp(buf, "$Nodes\n")) {
    readline(buf, 255, file, &line);
  }
  readline(buf, 255, file, &line);
  sscanf(buf, "%d", &vertices);

  i = 0;
  while (strcmp(buf, "$EndNodes\n")) {
    i++;
    readline(buf, 255, file, &line);
  }
  i--;
  assert(i == vertices);

  /****************************************************
   * Find number of triangles ...
   ****************************************************/

  while (strcmp(buf, "$Elements\n")) {
    readline(buf, 255, file, &line);
  }
  readline(buf, 255, file, &line);
  sscanf(buf, "%d", &triangles);

  i = 0;
  while (strcmp(buf, "$EndElements\n")) {
    c = sscanf(buf, "%d 2 2 %d %d %d %d %d", &tmp, &tmp2, &tmp3, &tmp4, &tmp5,
	       &tmp6);
    if (c != 6) {
      triangles--;
    }
    i++;
    readline(buf, 255, file, &line);
  }
  i--;
  triangles++;

  /****************************************************
   * Rewind file / create geometry structure
   ****************************************************/

#ifdef USE_ZLIB
  gzrewind(file);
#else
  rewind(file);
#endif
  line = 0;

  edges = 1.5 * triangles;
  gr = new_surface3d(vertices, edges, triangles);

  /****************************************************
   * Read vertices
   ****************************************************/

  while (strcmp(buf, "$Nodes\n")) {
    readline(buf, 255, file, &line);
  }
  readline(buf, 255, file, &line);
  i = 0;
  while (strcmp(buf, "$EndNodes\n")) {
    readline(buf, 255, file, &line);
    sscanf(buf, "%d %" SCANF_PREFIX "f %" SCANF_PREFIX "f %" SCANF_PREFIX "f",
	   &tmp, gr->x[i] + 0, gr->x[i] + 1, gr->x[i] + 2);
    i++;
  }
  i--;
  assert(i == vertices);

  /****************************************************
   * Read triangles
   ****************************************************/

  while (strcmp(buf, "$Elements\n")) {
    readline(buf, 255, file, &line);
  }
  readline(buf, 255, file, &line);
  i = 0;
  while (strcmp(buf, "$EndElements\n") && i <= triangles) {
    readline(buf, 255, file, &line);
    c = sscanf(buf, "%d 2 2 %d %d %d %d %d", &tmp, &tmp2, &tmp3, gr->t[i] + 0,
	       gr->t[i] + 1, gr->t[i] + 2);
    if (c == 6) {
      gr->t[i][0]--;
      gr->t[i][1]--;
      gr->t[i][2]--;
      i++;
    }
  }
  assert(i == triangles);

#ifdef USE_ZLIB
  gzclose(file);
#else
  fclose(file);
#endif

  /****************************************************
   * Generate Edges and arrays gr->s
   ****************************************************/

  prepare_edges(vertices, gr->t, triangles, gr->e, edges);
  prepare_arrays_s(gr->t, gr->e, triangles, edges, gr->s);

  prepare_surface3d(gr);

  return gr;
}

static void
prepare_dynamic_edges(uint vertices, uint(*t)[3], uint triangles,
		      uint(**e)[2], uint * edges)
{
  struct edge_list *(*list);
  struct edge_list *l;
  uint      i, j, k;

  list =
    (struct edge_list **) allocmem(vertices * sizeof(struct edge_list *));

  for (i = 0; i < vertices; ++i) {
    list[i] = (struct edge_list *) allocmem(sizeof(struct edge_list));
    list[i]->x = (uint) - 1;
    list[i]->next = NULL;
  }

  for (i = 0; i < triangles; ++i) {
    assert(t[i][0] < vertices);
    assert(t[i][1] < vertices);
    assert(t[i][2] < vertices);

    list[t[i][0]] = prepend_edge(list[t[i][0]], t[i][1]);
    list[t[i][0]] = prepend_edge(list[t[i][0]], t[i][2]);

    list[t[i][1]] = prepend_edge(list[t[i][1]], t[i][0]);
    list[t[i][1]] = prepend_edge(list[t[i][1]], t[i][2]);

    list[t[i][2]] = prepend_edge(list[t[i][2]], t[i][0]);
    list[t[i][2]] = prepend_edge(list[t[i][2]], t[i][1]);
  }

  *edges = 0;
  for (i = 0; i < vertices; ++i) {
    l = list[i];
    while (l->next != NULL) {
      l = l->next;
      (*edges)++;
    }
  }
  (*edges) = (*edges) / 2;

  (*e) = (uint(*)[2]) allocmem(*edges * sizeof(uint[2]));

  i = 0;
  for (j = 0; j < vertices; ++j) {
    l = list[j];
    while (l->next != NULL) {
      k = l->x;
      if (k != (uint) - 1) {
	(*e)[i][0] = j;
	(*e)[i][1] = k;
	mark_edge(list[j], k);
	mark_edge(list[k], j);
	i++;
      }
      l = l->next;
    }
  }
  assert(i == *edges);

  for (i = 0; i < vertices; ++i) {
    freemem(list[i]);
  }
  freemem(list);
}

psurface3d
read_unv_surface3d(char *filename)
{
#ifdef USE_ZLIB
  gzFile    file;
#else
  FILE     *file;
#endif
  uint      i, j, k, line, triangles, edges, vertices;
  char     *buf;
  uint     *x;
  uint     *t;
  uint      tmp1, tmp2, tmp3, tmp4, tmp5;

  psurface3d gr;

  buf = (char *) allocmem(255 * sizeof(char));
#ifdef USE_ZLIB
  file = gzopen(filename, "rb");
#else
  file = fopen(filename, "r");
#endif

  line = 0;
  readline(buf, 255, file, &line);

  /****************************************************
   * Find number of vertices ...
   ****************************************************/

  /* Scan fileheader until definition of nodes starts. */
  while (strncmp(buf, "  2411", 6)) {
    readline(buf, 255, file, &line);
  }

  i = 0;
  while (strncmp(buf, "    -1", 6)) {
    i++;
    readline(buf, 255, file, &line);
  }
  i--;
  /* Number of vertices found. */
  vertices = i / 2;

  /****************************************************
   * Find number of triangles ...
   ****************************************************/

  /* Scan fileheader until definition of triangles starts. */
  while (strncmp(buf, "  2412", 6)) {
    readline(buf, 255, file, &line);
  }

  i = 0;
  while (strncmp(buf, "    -1", 6)) {
    i++;
    readline(buf, 255, file, &line);
  }
  i--;
  /* Number of vertices found. */
  triangles = i / 2;

  edges = 1.5 * triangles;
  gr = new_surface3d(vertices, edges, triangles);
  x = allocuint(vertices);
  t = allocuint(triangles);

  printf("%d, %d, %d\n", vertices, edges, triangles);

#ifdef USE_ZLIB
  gzrewind(file);
#else
  rewind(file);
#endif

  line = 0;
  readline(buf, 255, file, &line);

  /****************************************************
   * Read vertices ...
   ****************************************************/

  /* Scan fileheader until definition of nodes starts. */
  while (strncmp(buf, "  2411", 6)) {
    readline(buf, 255, file, &line);
  }

  i = 0;
  readline(buf, 255, file, &line);
  while (strncmp(buf, "    -1", 6)) {
    assert(sscanf
	   (buf, "    %d         %d         %d         %d", x + i, &tmp1,
	    &tmp2, &tmp3)
	   == 4);
    readline(buf, 255, file, &line);
    assert(sscanf
	   (buf, "%" SCANF_PREFIX "f %" SCANF_PREFIX "f %" SCANF_PREFIX "f",
	    gr->x[i] + 0, gr->x[i] + 1, gr->x[i] + 2) == 3);
    readline(buf, 255, file, &line);
    i++;
  }

  /****************************************************
   * Read triangles ...
   ****************************************************/

  /* Scan fileheader until definition of triangles starts. */
  while (strncmp(buf, "  2412", 6)) {
    readline(buf, 255, file, &line);
  }

  i = 0;
  readline(buf, 255, file, &line);
  while (strncmp(buf, "    -1", 6)) {
    assert(sscanf(buf,
		  "    %d        %d         %d         %d         %d         %d",
		  t + i, &tmp1, &tmp2, &tmp3, &tmp4, &tmp5)
	   == 6);
    readline(buf, 255, file, &line);
    assert(sscanf(buf, "    %d    %d    %d", gr->t[i] + 0, gr->t[i] + 1,
		  gr->t[i] + 2)
	   == 3);
    readline(buf, 255, file, &line);
    i++;
  }

#ifdef USE_ZLIB
  gzclose(file);
#else
  fclose(file);
#endif

  /****************************************************
   * Assign correct vertex ordering
   ****************************************************/

  for (i = 0; i < triangles; ++i) {
    for (k = 0; k < 3; ++k) {
      for (j = 0; j < vertices; ++j) {
	if (gr->t[i][k] == x[j]) {
	  gr->t[i][k] = j;
	  break;
	}
      }
      assert(j < vertices);
    }
  }

  /****************************************************
   * Compute redundant edges/triangles information
   ****************************************************/

  freemem(gr->e);
  prepare_dynamic_edges(vertices, gr->t, triangles, &gr->e, &gr->edges);
  edges = gr->edges;
  freemem(gr->s);
  gr->s = (uint(*)[3]) allocuint(3 * edges);
  prepare_arrays_s(gr->t, gr->e, triangles, edges, gr->s);

  prepare_surface3d(gr);

  assert(check_surface3d(gr) == 0);
  //  assert(isclosed_surface3d(gr));
  assert(isoriented_surface3d(gr));

  freemem(buf);
  freemem(x);
  freemem(t);

  return gr;
}

psurface3d
refine_red_surface3d(psurface3d in)
{
  uint      triangles = in->triangles;
  uint      edges = in->edges;
  uint      vertices = in->vertices;

  psurface3d gr;
  uint      newtriangles, newedges, newvertices;
  uint      i, j, s, t, e, v;
  real      x[3];

  newtriangles = 4 * triangles;
  newedges = 2 * edges + 3 * triangles;
  newvertices = vertices + edges;

  gr = new_surface3d(newvertices, newedges, newtriangles);

  for (v = 0; v < vertices; ++v) {
    for (j = 0; j < 3; ++j) {
      gr->x[v][j] = in->x[v][j];
    }
  }

  for (e = 0; e < edges; ++e) {
    for (j = 0; j < 3; ++j) {
      x[j] = 0.5 * (in->x[in->e[e][0]][j] + in->x[in->e[e][1]][j]);
      gr->x[v][j] = x[j];
    }
    gr->e[2 * e][0] = in->e[e][0];
    gr->e[2 * e][1] = v;
    gr->e[2 * e + 1][0] = v;
    gr->e[2 * e + 1][1] = in->e[e][1];
    v++;
  }
  e *= 2;

  s = 0;
  for (t = 0; t < triangles; ++t) {
    for (j = 0; j < 3; ++j) {

      gr->e[e][0] = gr->e[2 * in->s[t][(j + 1) % 3] + 1][0];
      gr->e[e][1] = gr->e[2 * in->s[t][j] + 1][0];
      e++;

      if (in->e[in->s[t][(j + 1) % 3]][1] == in->e[in->s[t][j]][1]) {
	gr->s[s][0] = 2 * in->s[t][j] + 1;
	gr->s[s][1] = 2 * in->s[t][(j + 1) % 3] + 1;
	gr->s[s][2] = e - 1;

	gr->t[s][0] = gr->e[gr->s[s][1]][0];
	gr->t[s][1] = gr->e[gr->s[s][2]][1];
	gr->t[s][2] = gr->e[gr->s[s][0]][1];
      }
      else if (in->e[in->s[t][(j + 1) % 3]][0] == in->e[in->s[t][j]][1]) {
	gr->s[s][0] = 2 * in->s[t][j] + 1;
	gr->s[s][1] = 2 * in->s[t][(j + 1) % 3];
	gr->s[s][2] = e - 1;

	gr->t[s][0] = gr->e[gr->s[s][1]][1];
	gr->t[s][1] = gr->e[gr->s[s][2]][1];
	gr->t[s][2] = gr->e[gr->s[s][0]][1];
      }
      else if (in->e[in->s[t][(j + 1) % 3]][1] == in->e[in->s[t][j]][0]) {
	gr->s[s][0] = 2 * in->s[t][j];
	gr->s[s][1] = 2 * in->s[t][(j + 1) % 3] + 1;
	gr->s[s][2] = e - 1;

	gr->t[s][0] = gr->e[gr->s[s][1]][0];
	gr->t[s][1] = gr->e[gr->s[s][2]][1];
	gr->t[s][2] = gr->e[gr->s[s][0]][0];
      }
      else if (in->e[in->s[t][(j + 1) % 3]][0] == in->e[in->s[t][j]][0]) {
	gr->s[s][0] = 2 * in->s[t][j];
	gr->s[s][1] = 2 * in->s[t][(j + 1) % 3];
	gr->s[s][2] = e - 1;

	gr->t[s][0] = gr->e[gr->s[s][1]][1];
	gr->t[s][1] = gr->e[gr->s[s][2]][1];
	gr->t[s][2] = gr->e[gr->s[s][0]][0];
      }
      else {
	printf("ERROR!\n");
	abort();
      }

      s++;
    }
    for (i = 0; i < 3; ++i) {
      gr->s[s][i] = e - 3 + i;
      gr->t[s][i] = gr->e[e - 3 + ((i + 1) % 3)][0];
    }
    s++;

  }
  t *= 4;

  assert(t == newtriangles);
  assert(s == newtriangles);
  assert(e == newedges);
  assert(v == newvertices);

  prepare_surface3d(gr);

  return gr;
}
