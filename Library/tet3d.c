/* ------------------------------------------------------------
 * This is the file "tet3d.c" of the H2Lib package.
 * All rights reserved, Steffen Boerm 2015
 * ------------------------------------------------------------ */

#include "tet3d.h"

#include "basic.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>

ptet3d
new_tet3d(uint vertices, uint edges, uint faces, uint tetrahedra)
{
  ptet3d    t3;

  t3 = (ptet3d) allocmem(sizeof(tet3d));

  t3->x = (real(*)[3]) allocmem(sizeof(real[3]) * vertices);
  t3->e = (uint(*)[2]) allocmem(sizeof(uint[2]) * edges);
  t3->f = (uint(*)[3]) allocmem(sizeof(uint[3]) * faces);
  t3->t = (uint(*)[4]) allocmem(sizeof(uint[4]) * tetrahedra);

  t3->xb = (uint *) allocmem(sizeof(uint) * vertices);
  t3->eb = (uint *) allocmem(sizeof(uint) * edges);
  t3->fb = (uint *) allocmem(sizeof(uint) * faces);

  t3->vertices = vertices;
  t3->edges = edges;
  t3->faces = faces;
  t3->tetrahedra = tetrahedra;

  return t3;
}

void
del_tet3d(ptet3d t3)
{
  freemem(t3->fb);
  freemem(t3->eb);
  freemem(t3->xb);
  freemem(t3->t);
  freemem(t3->f);
  freemem(t3->e);
  freemem(t3->x);
  freemem(t3);
}

ptet3d
new_axis_tet3d()
{
  ptet3d    t3;

  t3 = new_tet3d(4, 6, 4, 1);

  t3->x[0][0] = 0.0;
  t3->x[0][1] = 0.0;
  t3->x[0][2] = 0.0;
  t3->x[1][0] = 1.0;
  t3->x[1][1] = 0.0;
  t3->x[1][2] = 0.0;
  t3->x[2][0] = 1.0;
  t3->x[2][1] = 1.0;
  t3->x[2][2] = 0.0;
  t3->x[3][0] = 1.0;
  t3->x[3][1] = 1.0;
  t3->x[3][2] = 1.0;

  t3->e[0][0] = 1;
  t3->e[0][1] = 2;
  t3->e[1][0] = 2;
  t3->e[1][1] = 3;
  t3->e[2][0] = 3;
  t3->e[2][1] = 1;
  t3->e[3][0] = 0;
  t3->e[3][1] = 1;
  t3->e[4][0] = 0;
  t3->e[4][1] = 2;
  t3->e[5][0] = 0;
  t3->e[5][1] = 3;

  t3->f[0][0] = 1;
  t3->f[0][1] = 2;
  t3->f[0][2] = 0;
  t3->f[1][0] = 1;
  t3->f[1][1] = 4;
  t3->f[1][2] = 5;
  t3->f[2][0] = 2;
  t3->f[2][1] = 5;
  t3->f[2][2] = 3;
  t3->f[3][0] = 0;
  t3->f[3][1] = 3;
  t3->f[3][2] = 4;

  t3->t[0][0] = 0;
  t3->t[0][1] = 1;
  t3->t[0][2] = 2;
  t3->t[0][3] = 3;

  t3->xb[0] = t3->xb[1] = t3->xb[2] = t3->xb[3] = 1;
  t3->eb[0] = t3->eb[1] = t3->eb[2] = t3->eb[3] = t3->eb[4] = t3->eb[5] = 1;
  t3->fb[0] = t3->fb[1] = t3->fb[2] = t3->fb[3] = 1;

  t3->vertices = 4;
  t3->edges = 6;
  t3->faces = 4;
  t3->tetrahedra = 1;

  return t3;
}

ptet3d
new_regular_tet3d()
{
  ptet3d    t3;
  real      r, h, c, y;

  r = 1.0;
  h = REAL_SQRT(2.0) - 1.0;
  c = REAL_SQRT(0.75);
  y = 1.0 - REAL_SQRT(0.5);

  t3 = new_tet3d(4, 6, 4, 1);

  t3->x[0][0] = 0.0;
  t3->x[0][1] = -h - y;
  t3->x[0][2] = -r;
  t3->x[1][0] = -c * r;
  t3->x[1][1] = -h - y;
  t3->x[1][2] = 0.5 * r;
  t3->x[2][0] = c * r;
  t3->x[2][1] = -h - y;
  t3->x[2][2] = 0.5 * r;
  t3->x[3][0] = 0.0;
  t3->x[3][1] = 1.0 - y;
  t3->x[3][2] = 0.0;

  t3->e[0][0] = 1;
  t3->e[0][1] = 2;
  t3->e[1][0] = 2;
  t3->e[1][1] = 3;
  t3->e[2][0] = 3;
  t3->e[2][1] = 1;
  t3->e[3][0] = 0;
  t3->e[3][1] = 1;
  t3->e[4][0] = 0;
  t3->e[4][1] = 2;
  t3->e[5][0] = 0;
  t3->e[5][1] = 3;

  t3->f[0][0] = 1;
  t3->f[0][1] = 2;
  t3->f[0][2] = 0;
  t3->f[1][0] = 1;
  t3->f[1][1] = 4;
  t3->f[1][2] = 5;
  t3->f[2][0] = 2;
  t3->f[2][1] = 5;
  t3->f[2][2] = 3;
  t3->f[3][0] = 0;
  t3->f[3][1] = 3;
  t3->f[3][2] = 4;

  t3->t[0][0] = 0;
  t3->t[0][1] = 1;
  t3->t[0][2] = 2;
  t3->t[0][3] = 3;

  t3->xb[0] = t3->xb[1] = t3->xb[2] = t3->xb[3] = 1;
  t3->eb[0] = t3->eb[1] = t3->eb[2] = t3->eb[3] = t3->eb[4] = t3->eb[5] = 1;
  t3->fb[0] = t3->fb[1] = t3->fb[2] = t3->fb[3] = 1;

  t3->vertices = 4;
  t3->edges = 6;
  t3->faces = 4;
  t3->tetrahedra = 1;

  return t3;
}

ptet3d
new_unitcube_tet3d()
{
  ptet3d    gr;
  ptet3dbuilder tb;
  real(*x)[3];

  tb = new_tet3dbuilder(8);

  x = getx_tet3dbuilder(tb);

  x[0][0] = 0.0;
  x[0][1] = 0.0;
  x[0][2] = 0.0;
  x[1][0] = 1.0;
  x[1][1] = 0.0;
  x[1][2] = 0.0;
  x[2][0] = 0.0;
  x[2][1] = 1.0;
  x[2][2] = 0.0;
  x[3][0] = 1.0;
  x[3][1] = 1.0;
  x[3][2] = 0.0;
  x[4][0] = 0.0;
  x[4][1] = 0.0;
  x[4][2] = 1.0;
  x[5][0] = 1.0;
  x[5][1] = 0.0;
  x[5][2] = 1.0;
  x[6][0] = 0.0;
  x[6][1] = 1.0;
  x[6][2] = 1.0;
  x[7][0] = 1.0;
  x[7][1] = 1.0;
  x[7][2] = 1.0;

  addtetrahedron_tet3dbuilder(tb, 0, 1, 3, 7);
  addtetrahedron_tet3dbuilder(tb, 0, 1, 5, 7);
  addtetrahedron_tet3dbuilder(tb, 0, 2, 3, 7);
  addtetrahedron_tet3dbuilder(tb, 0, 2, 6, 7);
  addtetrahedron_tet3dbuilder(tb, 0, 4, 5, 7);
  addtetrahedron_tet3dbuilder(tb, 0, 4, 6, 7);

  gr = buildmesh_tet3dbuilder(tb);

  del_tet3dbuilder(tb);

  return gr;
}

void
write_tet3d(pctet3d t3, const char *name)
{
  FILE     *out;
  const     real(*x)[3] = (const real(*)[3]) t3->x;
  const     uint(*e)[2] = (const uint(*)[2]) t3->e;
  const     uint(*f)[3] = (const uint(*)[3]) t3->f;
  const     uint(*t)[4] = (const uint(*)[4]) t3->t;
  const uint *xb = (const uint *) t3->xb;
  const uint *eb = (const uint *) t3->eb;
  const uint *fb = (const uint *) t3->fb;
  uint      vertices = t3->vertices;
  uint      edges = t3->edges;
  uint      faces = t3->faces;
  uint      tetrahedra = t3->tetrahedra;
  uint      i;

  out = fopen(name, "w");
  if (!out) {
    (void) fprintf(stderr, "Could not open file \"%s\" for writing\n", name);
    return;
  }

  (void) fprintf(out, "# Tetrahedral mesh description\n"
		 "# Vertices, edges, faces and tetrahedra\n"
		 "%u %u %u %u\n", vertices, edges, faces, tetrahedra);

  (void) fprintf(out,
		 "# List of vertices, given by coordinates and boundary flags\n");
  for (i = 0; i < vertices; i++)
    (void) fprintf(out, "%.12e %.12e %.12e  %u\n", x[i][0], x[i][1], x[i][2],
		   xb[i]);
  (void) fprintf(out,
		 "# List of edges, given by vertex numbers and boundary flags\n");
  for (i = 0; i < edges; i++)
    (void) fprintf(out, "%u %u  %u\n", e[i][0], e[i][1], eb[i]);
  (void) fprintf(out,
		 "# List of faces, given by edge numbers and boundary flags\n");
  for (i = 0; i < faces; i++)
    (void) fprintf(out, "%u %u %u  %u\n", f[i][0], f[i][1], f[i][2], fb[i]);
  (void) fprintf(out, "# List of tetrahedra, given by face numbers\n");
  for (i = 0; i < tetrahedra; i++)
    (void) fprintf(out, "%u %u %u %u\n", t[i][0], t[i][1], t[i][2], t[i][3]);

  (void) fclose(out);
}

ptet3d
read_tet3d(const char *name)
{
  FILE     *in;
  ptet3d    t3;
  real(*x)[3];
  uint(*e)[2];
  uint(*f)[3];
  uint(*t)[4];
  uint     *xb;
  uint     *eb;
  uint     *fb;
  uint      vertices;
  uint      edges;
  uint      faces;
  uint      tetrahedra;
  char      buf[80], *line;
  uint      i, items;

  in = fopen(name, "r");
  if (!in) {
    (void) fprintf(stderr, "Could not open file \"%s\" for reading\n", name);
    return 0;
  }

  line = fgets(buf, 80, in);
  while (line && line[0] == '#')
    line = fgets(buf, 80, in);
  items = sscanf(line, "%u %u %u %u", &vertices, &edges, &faces, &tetrahedra);
  if (items != 4) {
    (void) fprintf(stderr, "Could not get sizes from file \"%s\"\n", name);
    (void) fclose(in);
    return 0;
  }

  t3 = new_tet3d(vertices, edges, faces, tetrahedra);
  x = t3->x;
  e = t3->e;
  f = t3->f;
  t = t3->t;
  xb = t3->xb;
  eb = t3->eb;
  fb = t3->fb;

  for (i = 0; i < vertices; i++) {
    line = fgets(buf, 80, in);
    while (line && line[0] == '#')
      line = fgets(buf, 80, in);
    items = sscanf(line,
		   "%" SCANF_PREFIX "f %" SCANF_PREFIX "f %" SCANF_PREFIX
		   "f  %u", x[i], x[i] + 1, x[i] + 2, xb + i);
    if (items != 4) {
      (void) fprintf(stderr, "Could not read vertex %u from file \"%s\"\n", i,
		     name);
      del_tet3d(t3);
      (void) fclose(in);
      return 0;
    }
  }

  for (i = 0; i < edges; i++) {
    line = fgets(buf, 80, in);
    while (line && line[0] == '#')
      line = fgets(buf, 80, in);
    items = sscanf(line, "%u %u  %u", e[i], e[i] + 1, eb + i);
    if (items != 3) {
      (void) fprintf(stderr, "Could not read edge %u from file \"%s\"\n", i,
		     name);
      del_tet3d(t3);
      (void) fclose(in);
      return 0;
    }
  }

  for (i = 0; i < faces; i++) {
    line = fgets(buf, 80, in);
    while (line && line[0] == '#')
      line = fgets(buf, 80, in);
    items = sscanf(line, "%u %u %u  %u", f[i], f[i] + 1, f[i] + 2, fb + i);
    if (items != 4) {
      (void) fprintf(stderr, "Could not read face %u from file \"%s\"\n", i,
		     name);
      del_tet3d(t3);
      (void) fclose(in);
      return 0;
    }
  }

  for (i = 0; i < tetrahedra; i++) {
    line = fgets(buf, 80, in);
    while (line && line[0] == '#')
      line = fgets(buf, 80, in);
    items = sscanf(line, "%u %u %u %u", t[i], t[i] + 1, t[i] + 2, t[i] + 3);
    if (items != 4) {
      (void) fprintf(stderr,
		     "Could not read tetrahedron %u from file \"%s\"\n", i,
		     name);
      del_tet3d(t3);
      (void) fclose(in);
      return 0;
    }
  }

  (void) fclose(in);

  return t3;
}

/* Assuming that e1 and e2 share a vertex, find out whether this
 * vertex is the start or the end vertex of e1 */
static    uint
common_vertex(const uint(*e)[2], uint e1, uint e2)
{
  if (e[e1][0] == e[e2][0] || e[e1][0] == e[e2][1])
    return 0;

  assert(e[e1][1] == e[e2][0] || e[e1][1] == e[e2][1]);
  return 1;
}

/* Assuming that e1 and e2 share a vertex, find out its number */
static    uint
common_vertex_global(const uint(*e)[2], uint e1, uint e2)
{
  return e[e1][common_vertex(e, e1, e2)];
}

void
getvertices_byface_tet3d(pctet3d gr, uint tn, uint fl, uint v[])
{
  const     uint(*e)[2] = (const uint(*)[2]) gr->e;
  const     uint(*f)[3] = (const uint(*)[3]) gr->f;
  const     uint(*t)[4] = (const uint(*)[4]) gr->t;
  uint      f0, f1, ex;

  assert(tn < gr->tetrahedra);
  assert(fl < 4);

  /* Get global index of chosen face */
  f0 = t[tn][fl];

  /* Get global index of another face of the same tetrahedron */
  f1 = t[tn][(fl + 1) % 4];

  /* Vertex 1 is the intersection of edges 1 and 2 in face f0 */
  v[1] = common_vertex_global(e, f[f0][1], f[f0][2]);

  /* Vertex 2 is the intersection of edges 2 and 0 in face f0 */
  v[2] = common_vertex_global(e, f[f0][2], f[f0][0]);

  /* Vertex 3 is the intersection of edges 0 and 1 in face f0 */
  v[3] = common_vertex_global(e, f[f0][0], f[f0][1]);

  /* Find an edge that is not part of f0 */
  ex = f[f1][0];
  if (ex == f[f0][0] || ex == f[f0][1] || ex == f[f0][2]) {
    ex = f[f1][1];
    if (ex == f[f0][0] || ex == f[f0][1] || ex == f[f0][2]) {
      ex = f[f1][2];
      assert(ex != f[f0][0] && ex != f[f0][1] && ex != f[f0][2]);
    }
  }

  /* Find a vertex in ex that is not a vertex of f0 */
  if (e[ex][0] != v[1] && e[ex][0] != v[2] && e[ex][0] != v[3])
    v[0] = e[ex][0];
  else {
    assert(e[ex][1] != v[1] && e[ex][1] != v[2] && e[ex][1] != v[3]);
    v[0] = e[ex][1];
  }
}

/* Auxiliary function:
 find a common edge in f1 and f2 */
static    uint
common_edge(const uint(*f)[3], uint f1, uint f2)
{
  if (f[f1][0] == f[f2][0] || f[f1][0] == f[f2][1] || f[f1][0] == f[f2][2])
    return 0;

  if (f[f1][1] == f[f2][0] || f[f1][1] == f[f2][1] || f[f1][1] == f[f2][2])
    return 1;

  assert(f[f1][2] == f[f2][0] || f[f1][2] == f[f2][1]
	 || f[f1][2] == f[f2][2]);
  return 2;
}

static    uint
common_edge_global(const uint(*f)[3], uint f1, uint f2)
{
  return f[f1][common_edge(f, f1, f2)];
}

void
getvertices_tet3d(pctet3d t3, uint tn, uint v[])
{
  const     uint(*e)[2] = (const uint(*)[2]) t3->e;
  const     uint(*f)[3] = (const uint(*)[3]) t3->f;
  const     uint(*t)[4] = (const uint(*)[4]) t3->t;
  uint      e12, e23, e30, e01;

  e12 = common_edge_global(f, t[tn][1], t[tn][2]);
  e23 = common_edge_global(f, t[tn][2], t[tn][3]);

  v[0] = common_vertex_global(e, e12, e23);

  e30 = common_edge_global(f, t[tn][3], t[tn][0]);

  v[1] = common_vertex_global(e, e23, e30);

  e01 = common_edge_global(f, t[tn][0], t[tn][1]);
  v[2] = common_vertex_global(e, e30, e01);

  v[3] = common_vertex_global(e, e12, e01);
}

void
getedges_tet3d(pctet3d t3, uint tn, uint e[])
{
  const     uint(*f)[3] = (const uint(*)[3]) t3->f;
  const     uint(*t)[4] = (const uint(*)[4]) t3->t;

  e[0] = common_edge_global(f, t[tn][1], t[tn][2]);
  e[1] = common_edge_global(f, t[tn][2], t[tn][3]);
  e[2] = common_edge_global(f, t[tn][3], t[tn][1]);
  e[3] = common_edge_global(f, t[tn][0], t[tn][1]);
  e[4] = common_edge_global(f, t[tn][0], t[tn][2]);
  e[5] = common_edge_global(f, t[tn][0], t[tn][3]);
}

void 
getvertices_face_tet3d(pctet3d t3, uint nf, uint v[]) {

  uint e[3];

  /*Get edges of face nf*/
  e[0] = t3->f[nf][0];
  e[1] = t3->f[nf][1];
  e[2] = t3->f[nf][2];
  /*Get vertices */
  v[0] = t3->e[e[0]][0];
  v[1] = t3->e[e[0]][1];

  if (t3->e[e[1]][0] == v[0] || t3->e[e[1]][0] == v[1])
    v[2] = t3->e[e[1]][1];
  else {
    /*t3->e[e[1]][1] == v[0] || t3->e[e[1]][1] == v[1]*/
    v[2] = t3->e[e[1]][0];
  }
}

uint
fixnormals_tet3d(ptet3d gr)
{
  const     real(*x)[3] = (const real(*)[3]) gr->x;
  uint(*f)[3] = gr->f;
  const     uint(*t)[4] = (const uint(*)[4]) gr->t;
  const uint *fb = gr->fb;
  uint      tetrahedra = gr->tetrahedra;
  uint      v[4];
  real      n[3], xi[3], dx[3], dy[3];
  uint      fixed;
  uint      k, i, j, ex;

  fixed = 0;

  for (k = 0; k < tetrahedra; k++)
    for (i = 0; i < 4; i++) {
      j = t[k][i];
      if (fb[j]) {
	/* Get vertices, ensuring v[0] is opposite of face i */
	getvertices_byface_tet3d(gr, k, i, v);

	/* Compute spanning vectors of the face */
	dx[0] = x[v[2]][0] - x[v[1]][0];
	dx[1] = x[v[2]][1] - x[v[1]][1];
	dx[2] = x[v[2]][2] - x[v[1]][2];
	dy[0] = x[v[3]][0] - x[v[1]][0];
	dy[1] = x[v[3]][1] - x[v[1]][1];
	dy[2] = x[v[3]][2] - x[v[1]][2];

	/* Compute normal vector by cross product */
	n[0] = dx[1] * dy[2] - dx[2] * dy[1];
	n[1] = dx[2] * dy[0] - dx[0] * dy[2];
	n[2] = dx[0] * dy[1] - dx[1] * dy[0];

	/* Get coordinates of vertex 0 */
	xi[0] = x[v[0]][0] - x[v[1]][0];
	xi[1] = x[v[0]][1] - x[v[1]][1];
	xi[2] = x[v[0]][2] - x[v[1]][2];

	/* If n is an outer normal vector, vertex 0 should lie
	 * on the opposite side of the face */
	if (n[0] * xi[0] + n[1] * xi[1] + n[2] * xi[2] > 0) {
	  /* Swap edges 1 and 2 of face to ensure correct orientation */
	  ex = f[j][1];
	  f[j][1] = f[j][2];
	  f[j][2] = ex;

	  fixed++;
	}
      }
    }

  return fixed;
}

void
check_tet3d(pctet3d t3)
{
  uint      vertices = t3->vertices;
  uint      edges = t3->edges;
  uint      faces = t3->faces;
  uint      tetrahedra = t3->tetrahedra;
  const     uint(*e)[2] = (const uint(*)[2]) t3->e;
  const     uint(*f)[3] = (const uint(*)[3]) t3->f;
  const     uint(*t)[4] = (const uint(*)[4]) t3->t;
  const uint *fb = (const uint *) t3->fb;
  uint     *fn;
  uint      v[24];
  uint      i, j, k, l, m;
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

  /* Check if each face consists of three vertices */
  for (i = 0; i < faces; i++) {
    j = 0;
    for (k = 0; k < 3; k++) {
      if (f[i][k] >= edges) {
	(void) printf(" Face %u contains illegal edge %u\n", i, f[i][k]);
	errors++;
      }
      else {
	for (l = 0; l < j && e[f[i][k]][0] != v[l]; l++);
	if (l == j) {
	  v[j] = e[f[i][k]][0];
	  j++;
	}

	for (l = 0; l < j && e[f[i][k]][1] != v[l]; l++);
	if (l == j) {
	  v[j] = e[f[i][k]][1];
	  j++;
	}
      }
    }
    if (j != 3) {
      (void) printf(" Face %u consists of %u vertices:", i, j);
      while (j-- > 0)
	(void) printf(" %u", v[j]);
      (void) printf("\n");
      errors++;
    }
  }

  /* Check if each tetrahedron consists of six edges */
  for (i = 0; i < tetrahedra; i++) {
    j = 0;
    for (k = 0; k < 4; k++) {
      if (t[i][k] >= faces) {
	(void) printf(" Tetrahedron %u contains illegal face %u\n", i,
		      t[i][k]);
	errors++;
      }
      else {
	for (l = 0; l < 3; l++) {
	  for (m = 0; m < j && f[t[i][k]][l] != v[m]; m++);
	  if (m == j) {
	    v[j] = f[t[i][k]][l];
	    j++;
	  }
	}
      }
    }
    if (j != 6) {
      (void) printf(" Tetrahedron %u consists of %u edges:", i, j);
      while (j-- > 0)
	(void) printf(" %u", v[j]);
      (void) printf("\n");
      errors++;
    }
  }

  /* Check if each tetrahedron consists of four vertices */
  for (i = 0; i < tetrahedra; i++) {
    j = 0;
    for (k = 0; k < 4; k++) {
      for (l = 0; l < 3; l++) {
	for (m = 0; m < j && e[f[t[i][k]][l]][0] != v[m]; m++);
	if (m == j) {
	  v[j] = e[f[t[i][k]][l]][0];
	  j++;
	}

	for (m = 0; m < j && e[f[t[i][k]][l]][1] != v[m]; m++);
	if (m == j) {
	  v[j] = e[f[t[i][k]][l]][1];
	  j++;
	}
      }
    }
    if (j != 4) {
      (void) printf(" Tetrahedron %u consists of %u vertices:", i, j);
      while (j-- > 0)
	(void) printf(" %u", v[j]);
      (void) printf("\n");
      errors++;
    }
  }

  /* Check if each non-boundary face is adjacent to two tetrahedra */
  fn = (uint *) allocmem(sizeof(uint) * faces);
  for (i = 0; i < faces; i++)
    fn[i] = 0;
  for (i = 0; i < tetrahedra; i++) {
    fn[t[i][0]]++;
    fn[t[i][1]]++;
    fn[t[i][2]]++;
    fn[t[i][3]]++;
  }
  for (i = 0; i < faces; i++)
    if (!fb[i] && fn[i] != 2) {
      (void) printf(" Face %u adjacent to %u tetrahedra\n", i, fn[i]);
      errors++;
    }
  freemem(fn);

  (void) printf(" %u errors found\n", errors);

}

void
statistics_tet3d(pctet3d t3, preal hmin, preal hmax, preal volmin,
		 preal volmax, preal relvolmin, preal relvolmax)
{
  const     real(*x)[3] = (const real(*)[3]) t3->x;
  uint      tetrahedra = t3->tetrahedra;
  real      xl[4][3];
  real      xd[6][3];
  real      h, norm, vol;
  uint      v[4];
  uint      i, j;

  for (i = 0; i < tetrahedra; i++) {
    getvertices_tet3d(t3, i, v);

    /* Get coordinates of vertices */
    for (j = 0; j < 4; j++) {
      xl[j][0] = x[v[j]][0];
      xl[j][1] = x[v[j]][1];
      xl[j][2] = x[v[j]][2];
    }

    /* Compute edge vectors */
    xd[0][0] = xl[1][0] - xl[0][0];
    xd[0][1] = xl[1][1] - xl[0][1];
    xd[0][2] = xl[1][2] - xl[0][2];

    xd[1][0] = xl[2][0] - xl[0][0];
    xd[1][1] = xl[2][1] - xl[0][1];
    xd[1][2] = xl[2][2] - xl[0][2];

    xd[2][0] = xl[3][0] - xl[0][0];
    xd[2][1] = xl[3][1] - xl[0][1];
    xd[2][2] = xl[3][2] - xl[0][2];

    xd[3][0] = xl[2][0] - xl[1][0];
    xd[3][1] = xl[2][1] - xl[1][1];
    xd[3][2] = xl[2][2] - xl[1][2];

    xd[4][0] = xl[3][0] - xl[1][0];
    xd[4][1] = xl[3][1] - xl[1][1];
    xd[4][2] = xl[3][2] - xl[1][2];

    xd[5][0] = xl[3][0] - xl[2][0];
    xd[5][1] = xl[3][1] - xl[2][1];
    xd[5][2] = xl[3][2] - xl[2][2];

    /* Update h, hmin and hmax */
    h = 0.0;
    for (j = 0; j < 6; j++) {
      norm =
	REAL_SQRT(REAL_SQR(xd[j][0]) + REAL_SQR(xd[j][1]) +
		  REAL_SQR(xd[j][2]));
      if (i == 0 || norm < *hmin)
	*hmin = norm;
      if (i == 0 || norm > *hmax)
	*hmax = norm;
      if (norm > h)
	h = norm;
    }

    /* Compute volume */
    vol = fabs(xd[0][0] * xd[1][1] * xd[2][2] + xd[1][0] * xd[2][1] * xd[0][2]
	       + xd[2][0] * xd[0][1] * xd[1][2]
	       - xd[0][0] * xd[2][1] * xd[1][2]
	       - xd[1][0] * xd[0][1] * xd[2][2]
	       - xd[2][0] * xd[1][1] * xd[0][2])
      / 6.0;
    if (i == 0 || vol < *volmin)
      *volmin = vol;
    if (i == 0 || vol > *volmax)
      *volmax = vol;

    /* Compute relative volume */
    vol /= h * h * h;
    if (i == 0 || vol < *relvolmin)
      *relvolmin = vol;
    if (i == 0 || vol > *relvolmax)
      *relvolmax = vol;
  }
}

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

ptet3d
refine_tet3d(pctet3d t3, ptet3dref * t3r)
{
  ptet3d    r;
  uint      vertices = t3->vertices;
  uint      edges = t3->edges;
  uint      faces = t3->faces;
  uint      tetrahedra = t3->tetrahedra;
  const     real(*x)[3] = (const real(*)[3]) t3->x;
  const     uint(*e)[2] = (const uint(*)[2]) t3->e;
  const     uint(*f)[3] = (const uint(*)[3]) t3->f;
  const     uint(*t)[4] = (const uint(*)[4]) t3->t;
  const uint *xb = (const uint *) t3->xb;
  const uint *eb = (const uint *) t3->eb;
  const uint *fb = (const uint *) t3->fb;
  const     uint(*re)[2];
  uint     *xf, *xt, *ef, *et, *ff, *ft, *tf;
  uint      edge_vertices;	/* First vertex on an edge */
  uint      face_edges;		/* First edge on a face */
  uint      tetrahedron_edges;	/* First edge in a tetrahedron */
  uint      vertex_faces;	/* First vertex face in a tetrahedron */
  uint      center_faces;	/* First center face in a tetrahedron */
  uint      i, j;

  r = new_tet3d(vertices + edges, 2 * edges + 3 * faces + tetrahedra,
		4 * faces + 8 * tetrahedra, 8 * tetrahedra);
  re = (const uint(*)[2]) r->e;

  xf = xt = ef = et = ff = ft = tf = 0;
  if (t3r) {
    *t3r = (ptet3dref) allocmem(sizeof(tet3dref));
    xf = (*t3r)->xf = (uint *) allocmem(sizeof(uint) * r->vertices);
    xt = (*t3r)->xt = (uint *) allocmem(sizeof(uint) * r->vertices);
    ef = (*t3r)->ef = (uint *) allocmem(sizeof(uint) * r->edges);
    et = (*t3r)->et = (uint *) allocmem(sizeof(uint) * r->edges);
    ff = (*t3r)->ff = (uint *) allocmem(sizeof(uint) * r->faces);
    ft = (*t3r)->ft = (uint *) allocmem(sizeof(uint) * r->faces);
    tf = (*t3r)->tf = (uint *) allocmem(sizeof(uint) * r->tetrahedra);
  }

  /* Create vertices by copying old vertices */
  i = 0;
  for (j = 0; j < vertices; j++) {
    r->x[i][0] = x[j][0];
    r->x[i][1] = x[j][1];
    r->x[i][2] = x[j][2];
    r->xb[i] = xb[j];
    if (t3r) {
      xf[i] = j;
      xt[i] = 0;
    }
    i++;
  }

  /* Create vertices within edges */
  edge_vertices = i;
  for (j = 0; j < edges; j++) {
    r->x[i][0] = 0.5 * (x[e[j][0]][0] + x[e[j][1]][0]);
    r->x[i][1] = 0.5 * (x[e[j][0]][1] + x[e[j][1]][1]);
    r->x[i][2] = 0.5 * (x[e[j][0]][2] + x[e[j][1]][2]);
    r->xb[i] = eb[j];
    if (t3r) {
      xf[i] = j;
      xt[i] = 1;
    }
    i++;
  }
  assert(i == r->vertices);

  /* Create edges by splitting old edges */
  i = 0;
  for (j = 0; j < edges; j++) {
    r->e[i][0] = e[j][0];
    r->e[i][1] = edge_vertices + j;
    r->eb[i] = eb[j];
    if (t3r) {
      ef[i] = j;
      et[i] = 1;
    }
    i++;

    r->e[i][0] = edge_vertices + j;
    r->e[i][1] = e[j][1];
    r->eb[i] = eb[j];
    if (t3r) {
      ef[i] = j;
      et[i] = 1;
    }
    i++;
  }

  /* Create edges within faces */
  face_edges = i;
  for (j = 0; j < faces; j++) {
    r->e[i][0] = edge_vertices + f[j][1];
    r->e[i][1] = edge_vertices + f[j][2];
    r->eb[i] = fb[j];
    if (t3r) {
      ef[i] = j;
      et[i] = 2;
    }
    i++;

    r->e[i][0] = edge_vertices + f[j][2];
    r->e[i][1] = edge_vertices + f[j][0];
    r->eb[i] = fb[j];
    if (t3r) {
      ef[i] = j;
      et[i] = 2;
    }
    i++;

    r->e[i][0] = edge_vertices + f[j][0];
    r->e[i][1] = edge_vertices + f[j][1];
    r->eb[i] = fb[j];
    if (t3r) {
      ef[i] = j;
      et[i] = 2;
    }
    i++;
  }

  /* Create edges within tetrahedra */
  tetrahedron_edges = i;
  for (j = 0; j < tetrahedra; j++) {
    r->e[i][0] = edge_vertices + common_edge_global(f, t[j][3], t[j][1]);
    r->e[i][1] = edge_vertices + common_edge_global(f, t[j][2], t[j][0]);
    r->eb[i] = 0;
    if (t3r) {
      ef[i] = j;
      et[i] = 3;
    }
    i++;
  }
  assert(i == r->edges);

  /* Create faces by splitting old faces */
  i = 0;
  for (j = 0; j < faces; j++) {
    intersecting_edges(re, 2 * f[j][1], 2 * f[j][2], r->f[i] + 1,
		       r->f[i] + 2);
    r->f[i][0] = face_edges + 3 * j;
    r->fb[i] = fb[j];
    if (t3r) {
      ff[i] = j;
      ft[i] = 2;
    }
    i++;

    intersecting_edges(re, 2 * f[j][0], 2 * f[j][2], r->f[i], r->f[i] + 2);
    r->f[i][1] = face_edges + 3 * j + 1;
    r->fb[i] = fb[j];
    if (t3r) {
      ff[i] = j;
      ft[i] = 2;
    }
    i++;

    intersecting_edges(re, 2 * f[j][0], 2 * f[j][1], r->f[i], r->f[i] + 1);
    r->f[i][2] = face_edges + 3 * j + 2;
    r->fb[i] = fb[j];
    if (t3r) {
      ff[i] = j;
      ft[i] = 2;
    }
    i++;

    r->f[i][0] = face_edges + 3 * j;
    r->f[i][1] = face_edges + 3 * j + 1;
    r->f[i][2] = face_edges + 3 * j + 2;
    r->fb[i] = fb[j];
    if (t3r) {
      ff[i] = j;
      ft[i] = 2;
    }
    i++;
  }

  /* Create faces closest to vertices within tetrahedra */
  vertex_faces = i;
  for (j = 0; j < tetrahedra; j++) {
    r->f[i][0] = face_edges + 3 * t[j][1] + common_edge(f, t[j][1], t[j][0]);
    r->f[i][1] = face_edges + 3 * t[j][2] + common_edge(f, t[j][2], t[j][0]);
    r->f[i][2] = face_edges + 3 * t[j][3] + common_edge(f, t[j][3], t[j][0]);
    r->fb[i] = 0;
    if (t3r) {
      ff[i] = j;
      ft[i] = 3;
    }
    i++;

    r->f[i][0] = face_edges + 3 * t[j][2] + common_edge(f, t[j][2], t[j][1]);
    r->f[i][1] = face_edges + 3 * t[j][3] + common_edge(f, t[j][3], t[j][1]);
    r->f[i][2] = face_edges + 3 * t[j][0] + common_edge(f, t[j][0], t[j][1]);
    r->fb[i] = 0;
    if (t3r) {
      ff[i] = j;
      ft[i] = 3;
    }
    i++;

    r->f[i][0] = face_edges + 3 * t[j][3] + common_edge(f, t[j][3], t[j][2]);
    r->f[i][1] = face_edges + 3 * t[j][0] + common_edge(f, t[j][0], t[j][2]);
    r->f[i][2] = face_edges + 3 * t[j][1] + common_edge(f, t[j][1], t[j][2]);
    r->fb[i] = 0;
    if (t3r) {
      ff[i] = j;
      ft[i] = 3;
    }
    i++;

    r->f[i][0] = face_edges + 3 * t[j][0] + common_edge(f, t[j][0], t[j][3]);
    r->f[i][1] = face_edges + 3 * t[j][1] + common_edge(f, t[j][1], t[j][3]);
    r->f[i][2] = face_edges + 3 * t[j][2] + common_edge(f, t[j][2], t[j][3]);
    r->fb[i] = 0;
    if (t3r) {
      ff[i] = j;
      ft[i] = 3;
    }
    i++;
  }

  /* Create faces touching the central diagonal within tetrahedra */
  center_faces = i;
  for (j = 0; j < tetrahedra; j++) {
    r->f[i][0] = tetrahedron_edges + j;
    r->f[i][1] = face_edges + 3 * t[j][0] + common_edge(f, t[j][0], t[j][1]);
    r->f[i][2] = face_edges + 3 * t[j][3] + common_edge(f, t[j][3], t[j][2]);
    r->fb[i] = 0;
    if (t3r) {
      ff[i] = j;
      ft[i] = 3;
    }
    i++;

    r->f[i][0] = tetrahedron_edges + j;
    r->f[i][1] = face_edges + 3 * t[j][0] + common_edge(f, t[j][0], t[j][3]);
    r->f[i][2] = face_edges + 3 * t[j][1] + common_edge(f, t[j][1], t[j][2]);
    r->fb[i] = 0;
    if (t3r) {
      ff[i] = j;
      ft[i] = 3;
    }
    i++;

    r->f[i][0] = tetrahedron_edges + j;
    r->f[i][1] = face_edges + 3 * t[j][2] + common_edge(f, t[j][2], t[j][3]);
    r->f[i][2] = face_edges + 3 * t[j][1] + common_edge(f, t[j][1], t[j][0]);
    r->fb[i] = 0;
    if (t3r) {
      ff[i] = j;
      ft[i] = 3;
    }
    i++;

    r->f[i][0] = tetrahedron_edges + j;
    r->f[i][1] = face_edges + 3 * t[j][2] + common_edge(f, t[j][2], t[j][1]);
    r->f[i][2] = face_edges + 3 * t[j][3] + common_edge(f, t[j][3], t[j][0]);
    if (t3r) {
      ff[i] = j;
      ft[i] = 3;
    }
    r->fb[i] = 0;
    i++;
  }
  assert(i == r->faces);

  /* Create vertex tetrahedra */
  i = 0;
  for (j = 0; j < tetrahedra; j++) {
    r->t[i][0] = vertex_faces + 4 * j;
    r->t[i][1] = 4 * t[j][1] + common_edge(f, t[j][1], t[j][0]);
    r->t[i][2] = 4 * t[j][2] + common_edge(f, t[j][2], t[j][0]);
    r->t[i][3] = 4 * t[j][3] + common_edge(f, t[j][3], t[j][0]);
    if (t3r)
      tf[i] = j;
    i++;

    r->t[i][0] = 4 * t[j][0] + common_edge(f, t[j][0], t[j][1]);
    r->t[i][1] = vertex_faces + 4 * j + 1;
    r->t[i][2] = 4 * t[j][2] + common_edge(f, t[j][2], t[j][1]);
    r->t[i][3] = 4 * t[j][3] + common_edge(f, t[j][3], t[j][1]);
    if (t3r)
      tf[i] = j;
    i++;

    r->t[i][0] = 4 * t[j][0] + common_edge(f, t[j][0], t[j][2]);
    r->t[i][1] = 4 * t[j][1] + common_edge(f, t[j][1], t[j][2]);
    r->t[i][2] = vertex_faces + 4 * j + 2;
    r->t[i][3] = 4 * t[j][3] + common_edge(f, t[j][3], t[j][2]);
    if (t3r)
      tf[i] = j;
    i++;

    r->t[i][0] = 4 * t[j][0] + common_edge(f, t[j][0], t[j][3]);
    r->t[i][1] = 4 * t[j][1] + common_edge(f, t[j][1], t[j][3]);
    r->t[i][2] = 4 * t[j][2] + common_edge(f, t[j][2], t[j][3]);
    r->t[i][3] = vertex_faces + 4 * j + 3;
    if (t3r)
      tf[i] = j;
    i++;
  }

  /* Create interior tetrahedra */
  for (j = 0; j < tetrahedra; j++) {
    r->t[i][0] = center_faces + 4 * j + 2;
    r->t[i][1] = 4 * t[j][2] + 3;
    r->t[i][2] = center_faces + 4 * j + 3;
    r->t[i][3] = vertex_faces + 4 * j;
    if (t3r)
      tf[i] = j;
    i++;

    r->t[i][0] = center_faces + 4 * j;
    r->t[i][1] = vertex_faces + 4 * j + 1;
    r->t[i][2] = center_faces + 4 * j + 3;
    r->t[i][3] = 4 * t[j][3] + 3;
    if (t3r)
      tf[i] = j;
    i++;

    r->t[i][0] = vertex_faces + 4 * j + 3;
    r->t[i][1] = center_faces + 4 * j + 1;
    r->t[i][2] = 4 * t[j][1] + 3;
    r->t[i][3] = center_faces + 4 * j + 2;
    if (t3r)
      tf[i] = j;
    i++;

    r->t[i][0] = 4 * t[j][0] + 3;
    r->t[i][1] = center_faces + 4 * j + 1;
    r->t[i][2] = vertex_faces + 4 * j + 2;
    r->t[i][3] = center_faces + 4 * j;
    if (t3r)
      tf[i] = j;
    i++;
  }
  assert(i == r->tetrahedra);

  return r;
}

void
del_tet3dref(ptet3dref t3r)
{
  freemem(t3r->tf);
  freemem(t3r->ft);
  freemem(t3r->ff);
  freemem(t3r->et);
  freemem(t3r->ef);
  freemem(t3r->xt);
  freemem(t3r->xf);
  freemem(t3r);
}

/* ------------------------------------------------------------
 * Tool for constructing meshes
 * ------------------------------------------------------------ */

typedef struct _edgeentry edgeentry;
typedef edgeentry *pedgeentry;

typedef struct _faceentry faceentry;
typedef faceentry *pfaceentry;

typedef struct _tetentry tetentry;
typedef tetentry *ptetentry;

struct _tet3dbuilder {
  real(*x)[3];

  pedgeentry *elist;
  ptetentry tlist;

  uint      vertices;
  uint      edges;
  uint      faces;
  uint      tetrahedra;
};

struct _edgeentry {
  uint      name;

  uint      v[2];

  pfaceentry flist;

  bool      boundary;

  pedgeentry next;
};

struct _faceentry {
  uint      name;

  pedgeentry e[3];

  uint      neighbours;

  pfaceentry next;
};

struct _tetentry {
  uint      name;

  pfaceentry f[4];

  ptetentry next;
};

static    pedgeentry
new_edgeentry(uint name, uint v0, uint v1, pedgeentry next)
{
  pedgeentry e;

  e = (pedgeentry) allocmem(sizeof(edgeentry));
  e->name = name;
  e->v[0] = v0;
  e->v[1] = v1;
  e->flist = 0;
  e->boundary = false;
  e->next = next;

  return e;
}

static    pfaceentry
new_faceentry(uint name, pedgeentry e0, pedgeentry e1,
	      pedgeentry e2, pfaceentry next)
{
  pfaceentry f;

  f = (pfaceentry) allocmem(sizeof(faceentry));
  f->name = name;
  f->e[0] = e0;
  f->e[1] = e1;
  f->e[2] = e2;
  f->neighbours = 0;
  f->next = next;

  return f;
}

static    ptetentry
new_tetentry(uint name, pfaceentry f0, pfaceentry f1,
	     pfaceentry f2, pfaceentry f3, ptetentry next)
{
  ptetentry t;

  t = (ptetentry) allocmem(sizeof(tetentry));
  t->name = name;
  t->f[0] = f0;
  t->f[1] = f1;
  t->f[2] = f2;
  t->f[3] = f3;
  t->next = next;

  f0->neighbours++;
  f1->neighbours++;
  f2->neighbours++;
  f3->neighbours++;

  return t;
}

ptet3dbuilder
new_tet3dbuilder(uint vertices)
{
  ptet3dbuilder tb;
  uint      i;

  tb = (ptet3dbuilder) allocmem(sizeof(tet3dbuilder));
  tb->x = (real(*)[3]) allocmem(sizeof(real[3]) * vertices);
  tb->elist = (pedgeentry *) allocmem(sizeof(pedgeentry) * vertices);
  tb->tlist = 0;

  tb->vertices = vertices;
  tb->edges = 0;
  tb->faces = 0;
  tb->tetrahedra = 0;

  for (i = 0; i < vertices; i++)
    tb->elist[i] = 0;

  return tb;
}

static void
del_tetentry(ptetentry t)
{
  ptetentry next;

  while (t) {
    next = t->next;

    freemem(t);

    t = next;
  }
}

static void
del_faceentry(pfaceentry f)
{
  pfaceentry next;

  while (f) {
    next = f->next;

    freemem(f);

    f = next;
  }
}

static void
del_edgeentry(pedgeentry e)
{
  pedgeentry next;

  while (e) {
    next = e->next;

    del_faceentry(e->flist);

    freemem(e);

    e = next;
  }
}

void
del_tet3dbuilder(ptet3dbuilder tb)
{
  uint      i;

  for (i = 0; i < tb->vertices; i++)
    del_edgeentry(tb->elist[i]);

  del_tetentry(tb->tlist);

  freemem(tb->elist);
  freemem(tb->x);
  freemem(tb);
}

real(*getx_tet3dbuilder(ptet3dbuilder tb))[3]
{
  return tb->x;
}

static    pedgeentry
find_edgeentry(ptet3dbuilder tb, uint v0, uint v1)
{
  pedgeentry e;
  uint      vx;

  assert(v0 < tb->vertices);
  assert(v1 < tb->vertices);

  /* Vertices in an edge are numbered in ascending order */
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

static    pfaceentry
find_faceentry(ptet3dbuilder tb, pedgeentry e0, pedgeentry e1, pedgeentry e2)
{
  pfaceentry f;
  pedgeentry ex;

  /* Edges in a face are numbered in ascending order */
  if (e0->name > e1->name) {
    ex = e0;
    e0 = e1;
    e1 = ex;
  }
  if (e1->name > e2->name) {
    ex = e1;
    e1 = e2;
    e2 = ex;
  }
  if (e0->name > e1->name) {
    ex = e0;
    e0 = e1;
    e1 = ex;
  }

  /* Check if this face is already known */
  for (f = e0->flist; f; f = f->next)
    if (f->e[1] == e1 && f->e[2] == e2)
      break;

  /* If not known, create it */
  if (f == 0) {
    e0->flist = f = new_faceentry(tb->faces, e0, e1, e2, e0->flist);
    tb->faces++;
  }

  return f;
}

void
addtetrahedron_tet3dbuilder(ptet3dbuilder tb, uint v0, uint v1,
			    uint v2, uint v3)
{
  pedgeentry e01, e02, e03, e12, e13, e23;
  pfaceentry f0, f1, f2, f3;

  e01 = find_edgeentry(tb, v0, v1);
  e02 = find_edgeentry(tb, v0, v2);
  e03 = find_edgeentry(tb, v0, v3);
  e12 = find_edgeentry(tb, v1, v2);
  e13 = find_edgeentry(tb, v1, v3);
  e23 = find_edgeentry(tb, v2, v3);

  f0 = find_faceentry(tb, e12, e23, e13);
  f1 = find_faceentry(tb, e02, e23, e03);
  f2 = find_faceentry(tb, e01, e13, e03);
  f3 = find_faceentry(tb, e01, e12, e02);

  tb->tlist = new_tetentry(tb->tetrahedra, f0, f1, f2, f3, tb->tlist);
  tb->tetrahedra++;
}

ptet3d
buildmesh_tet3dbuilder(ptet3dbuilder tb)
{
  ptet3d    gr;
  pedgeentry e;
  pfaceentry f;
  ptetentry t;
  uint      i;

  /* Create mesh */
  gr = new_tet3d(tb->vertices, tb->edges, tb->faces, tb->tetrahedra);

  /* Reset vertex boundary flags */
  for (i = 0; i < tb->vertices; i++)
    gr->xb[i] = false;

  /* Copy tetrahedra */
  for (t = tb->tlist; t; t = t->next) {
    gr->t[t->name][0] = t->f[0]->name;
    gr->t[t->name][1] = t->f[1]->name;
    gr->t[t->name][2] = t->f[2]->name;
    gr->t[t->name][3] = t->f[3]->name;
  }

  /* Copy faces */
  for (i = 0; i < tb->vertices; i++)
    for (e = tb->elist[i]; e; e = e->next)
      for (f = e->flist; f; f = f->next) {
	gr->f[f->name][0] = f->e[0]->name;
	gr->f[f->name][1] = f->e[1]->name;
	gr->f[f->name][2] = f->e[2]->name;

	/* Determine whether this a boundary face or not */
	assert(f->neighbours == 1 || f->neighbours == 2);
	gr->fb[f->name] = (f->neighbours == 1);

	/* If it is a boundary face, its edges are boundary edges */
	if (f->neighbours == 1) {
	  f->e[0]->boundary = true;
	  f->e[1]->boundary = true;
	  f->e[2]->boundary = true;
	}
      }

  /* Copy edges */
  for (i = 0; i < tb->vertices; i++)
    for (e = tb->elist[i]; e; e = e->next) {
      gr->e[e->name][0] = e->v[0];
      gr->e[e->name][1] = e->v[1];
      gr->eb[e->name] = e->boundary;

      /* If this is a boundary edge, its vertices are boundary vertices */
      if (e->boundary) {
	gr->xb[e->v[0]] = true;
	gr->xb[e->v[1]] = true;
      }
    }

  /* Copy vertices */
  for (i = 0; i < tb->vertices; i++) {
    gr->x[i][0] = tb->x[i][0];
    gr->x[i][1] = tb->x[i][1];
    gr->x[i][2] = tb->x[i][2];
  }

  /* Ensure that boundary faces are oriented counter-clockwise */
  fixnormals_tet3d(gr);

  return gr;
}
