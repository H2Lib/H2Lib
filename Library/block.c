
/* ------------------------------------------------------------
 * This is the file "block.c" of the H2Lib package.
 * All rights reserved, Knut Reimer 2009
 * ------------------------------------------------------------ */

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifdef USE_CAIRO
#include <cairo/cairo.h>
#endif
#ifdef USE_FREEGLUT
#include <GL/freeglut.h>
#include <GL/glu.h>
#include <GL/gl.h>
#endif

#include "basic.h"
#include "cluster.h"
#include "block.h"

/* ------------------------------------------------------------
 * Admissibility conditions
 * ------------------------------------------------------------ */

bool
admissible_2_cluster(pcluster rc, pcluster cc, void *data)
{
  real      eta = *(real *) data;
  real      diamt, diams, dist, a;
  uint      i;

  diamt = 0.0;
  diams = 0.0;
  dist = 0.0;

  for (i = 0; i < rc->dim; i++) {
    a = rc->bmax[i] - rc->bmin[i];
    diamt += a * a;
    a = cc->bmax[i] - cc->bmin[i];
    diams += a * a;
    a = REAL_MAX3(0.0, rc->bmin[i] - cc->bmax[i], cc->bmin[i] - rc->bmax[i]);
    dist += a * a;
  }

  return (diamt <= eta * eta * dist && diams <= eta * eta * dist);
}

bool
admissible_max_cluster(pcluster rc, pcluster cc, void *data)
{
  real      eta = *(real *) data;
  real      diamt, diams, dist;

  diamt = getdiam_max_cluster(rc);
  diams = getdiam_max_cluster(cc);
  dist = getdist_max_cluster(rc, cc);

  return (diamt <= eta * dist && diams <= eta * dist);
}

bool
admissible_sphere_cluster(pcluster rc, pcluster cc, void *data)
{
  real      eta = *(real *) data;
  uint      dim = rc->dim;
  real      diamt, diams, dist;
  uint      i;

  assert(cc->dim == dim);

  diamt = getdiam_2_cluster(rc);
  diams = getdiam_2_cluster(cc);

  dist = 0.0;
  for (i = 0; i < dim; ++i) {
    dist +=
      REAL_SQR(0.5 *
	       ((rc->bmin[i] + rc->bmax[i]) - (cc->bmin[i] + cc->bmax[i])));
  }
  dist = REAL_SQRT(dist) - 0.5 * (diamt + diams);

  return (diamt <= eta * dist && diams <= eta * dist);
}

bool
admissible_2_min_cluster(pcluster rc, pcluster cc, void *data)
{
  real      eta = *(real *) data;
  real      diamt, diams, dist;

  diamt = getdiam_2_cluster(rc);
  diams = getdiam_2_cluster(cc);
  dist = getdist_2_cluster(rc, cc);

  return (diamt <= eta * dist || diams <= eta * dist);
}

/* ------------------------------------------------------------
 Constructors and destructors
 ------------------------------------------------------------ */

pblock
new_block(pcluster rc, pcluster cc, bool a, uint rsons, uint csons)
{
  pblock    b;

  uint      i;

  b = (pblock) allocmem((size_t) sizeof(block));
  b->rc = rc;
  b->cc = cc;
  b->a = a;
  b->rsons = rsons;
  b->csons = csons;
  b->son = NULL;
  if (rsons > 0 && csons > 0)
    b->son = (pblock *) allocmem((size_t) rsons * csons * sizeof(pblock));
  for (i = 0; i < rsons * csons; i++) {
    b->son[i] = NULL;
  }

  return b;
}

void
del_block(pblock b)
{
  uint      i, j;

  for (j = 0; j < b->csons; j++)
    for (i = 0; i < b->rsons; i++)
      del_block(b->son[i + j * b->rsons]);
  freemem(b->son);
  freemem(b);
}

/* ------------------------------------------------------------
 Block clustering strategies
 ------------------------------------------------------------ */

void
update_block(pblock b)
{
  uint      rsons = b->rsons;
  uint      csons = b->csons;

  uint      desc;

  uint      i, j;

  desc = 1;
  for (j = 0; j < csons; j++)
    for (i = 0; i < rsons; i++)
      desc += b->son[i + j * rsons]->desc;

  b->desc = desc;
}

pblock
build_nonstrict_block(pcluster rc, pcluster cc, void *eta, admissible admis)
{
  pblock    b;

  bool      a;
  uint      rsons, csons;
  uint      i, j;

  a = admis(rc, cc, eta);

  if (a == false) {
    /* inadmissible leaf */
    if (rc->sons * cc->sons == 0) {
      rsons = 0;
      csons = 0;
    }
    /* no leaf */
    else {
      rsons = rc->sons;
      csons = cc->sons;
    }
  }
  /* admissible leaf */
  else {
    assert(a == true);
    rsons = 0;
    csons = 0;
  }

  b = new_block(rc, cc, a, rsons, csons);

  for (i = 0; i < rsons; i++) {
    for (j = 0; j < csons; j++) {
      b->son[i + j * rsons] =
	build_nonstrict_block(rc->son[i], cc->son[j], eta, admis);
    }
  }

  update_block(b);

  return b;
}

pblock
build_nonstrict_lower_block(pcluster rc, pcluster cc, void *eta,
			    admissible admis)
{
  pblock    b;

  bool      a;
  uint      rsons, csons;
  uint      i, j;

  a = admis(rc, cc, eta);

  if (a == false) {
    /* inadmissible leaf */
    if (rc->sons * cc->sons == 0) {
      rsons = 0;
      csons = 0;
    }
    /* no leaf */
    else {
      rsons = rc->sons;
      csons = cc->sons;
    }
  }
  /* admissible leaf */
  else {
    assert(a == true);
    rsons = 0;
    csons = 0;
  }

  b = new_block(rc, cc, a, rsons, csons);

  for (j = 0; j < csons; j++) {
    for (i = 0; i < j; ++i) {
      b->son[i + j * rsons] = new_block(rc->son[i], cc->son[j], true, 0, 0);
      update_block(b->son[i + j * rsons]);
    }
    b->son[i + j * rsons] =
      build_nonstrict_lower_block(rc->son[i], cc->son[j], eta, admis);
    for (i = j + 1; i < rsons; i++) {
      b->son[i + j * rsons] =
	build_nonstrict_block(rc->son[i], cc->son[j], eta, admis);
    }
  }

  update_block(b);

  return b;
}

pblock
build_strict_block(pcluster rc, pcluster cc, void *eta, admissible admis)
{
  pblock    b;

  bool      a;
  uint      rsons, csons;
  uint      i, j;

  a = admis(rc, cc, eta);

  if (a == false) {
    if (rc->sons == 0) {
      /* inadmissible leaf */
      if (cc->sons == 0) {
	rsons = 0;
	csons = 0;
	b = new_block(rc, cc, a, rsons, csons);
      }
      /* no leaf */
      else {
	rsons = 1;
	csons = cc->sons;
	b = new_block(rc, cc, a, rsons, csons);
	for (j = 0; j < csons; j++) {
	  b->son[j] = build_strict_block(rc, cc->son[j], eta, admis);
	}
      }
    }
    else {
      /* no leaf */
      if (cc->sons == 0) {
	rsons = rc->sons;
	csons = 1;
	b = new_block(rc, cc, a, rsons, csons);
	for (i = 0; i < rsons; i++) {
	  b->son[i] = build_strict_block(rc->son[i], cc, eta, admis);
	}
      }
      /* no leaf */
      else {
	rsons = rc->sons;
	csons = cc->sons;
	b = new_block(rc, cc, a, rsons, csons);
	for (i = 0; i < rsons; i++) {
	  for (j = 0; j < csons; j++) {
	    b->son[i + j * rsons] = build_strict_block(rc->son[i], cc->son[j],
						       eta, admis);
	  }
	}
      }
    }
  }
  /* admissible leaf */
  else {
    assert(a == true);
    rsons = 0;
    csons = 0;
    b = new_block(rc, cc, a, rsons, csons);
  }

  update_block(b);

  return b;
}

pblock
build_strict_lower_block(pcluster rc, pcluster cc, void *eta,
			 admissible admis)
{
  pblock    b;

  bool      a;
  uint      rsons, csons;
  uint      i, j;

  a = admis(rc, cc, eta);

  if (a == false) {
    if (rc->sons == 0) {
      /* inadmissible leaf */
      if (cc->sons == 0) {
	rsons = 0;
	csons = 0;
	b = new_block(rc, cc, a, rsons, csons);
      }
      /* no leaf */
      else {
	rsons = 1;
	csons = cc->sons;
	b = new_block(rc, cc, a, rsons, csons);
	for (j = 0; j < csons; j++) {
	  b->son[j] = build_strict_block(rc, cc->son[j], eta, admis);
	}
      }
    }
    else {
      /* no leaf */
      if (cc->sons == 0) {
	rsons = rc->sons;
	csons = 1;
	b = new_block(rc, cc, a, rsons, csons);
	for (i = 0; i < rsons; i++) {
	  b->son[i] = build_strict_block(rc->son[i], cc, eta, admis);
	}
      }
      /* no leaf */
      else {
	rsons = rc->sons;
	csons = cc->sons;
	b = new_block(rc, cc, a, rsons, csons);
	for (j = 0; j < csons; j++) {
	  for (i = 0; i < j; ++i) {
	    b->son[i + j * rsons] = new_block(rc->son[i], cc->son[j], true, 0,
					      0);
	    update_block(b->son[i + j * rsons]);
	  }
	  b->son[i + j * rsons] = build_strict_lower_block(rc->son[i],
							   cc->son[j], eta,
							   admis);
	  for (i = j + 1; i < rsons; i++) {
	    b->son[i + j * rsons] = build_strict_block(rc->son[i], cc->son[j],
						       eta, admis);
	  }
	}
      }
    }
  }
  /* admissible leaf */
  else {
    assert(a == true);
    rsons = 0;
    csons = 0;
    b = new_block(rc, cc, a, rsons, csons);
  }

  update_block(b);

  return b;
}

/* ------------------------------------------------------------
 Drawing block cluster trees
 ------------------------------------------------------------ */

#ifdef USE_CAIRO
static void
draw_cairo_subblock(cairo_t * cr, pcblock b, int levels)
{
  uint      rsons, csons;
  uint      rsize, csize;
  uint      roff, coff;
  uint      i, j;

  if (b->son && levels != 1) {
    rsons = b->rsons;
    csons = b->csons;

    coff = 0;
    for (j = 0; j < csons; j++) {
      roff = 0;
      for (i = 0; i < rsons; i++) {
	cairo_save(cr);
	cairo_translate(cr, coff, roff);
	draw_cairo_subblock(cr, b->son[i + j * rsons], levels - 1);
	cairo_restore(cr);

	roff += b->son[i + j * rsons]->rc->size;
      }
      assert(roff == b->rc->size);

      coff += b->son[j * rsons]->cc->size;
    }
    assert(coff == b->cc->size);
  }
  else {
    rsize = b->rc->size;
    csize = b->cc->size;

    if (b->son) {
      cairo_new_path(cr);
      cairo_set_source_rgba(cr, 0.0, 0.0, 0.0, 1.0);
      cairo_rectangle(cr, 0.0, 0.0, csize, rsize);
      cairo_stroke(cr);
    }
    else if (b->a > 0) {
      cairo_new_path(cr);
      cairo_set_source_rgba(cr, 0.0, 0.0, 0.0, 1.0);
      cairo_rectangle(cr, 0.0, 0.0, csize, rsize);
      cairo_save(cr);
      cairo_set_source_rgba(cr, 0.2, 0.2, 1.0, 1.0);
      cairo_fill_preserve(cr);
      cairo_restore(cr);
      cairo_stroke(cr);
    }
    else {
      cairo_new_path(cr);
      cairo_set_source_rgba(cr, 0.0, 0.0, 0.0, 1.0);
      cairo_rectangle(cr, 0.0, 0.0, csize, rsize);
      cairo_save(cr);
      cairo_set_source_rgba(cr, 1.0, 0.0, 0.0, 1.0);
      cairo_fill_preserve(cr);
      cairo_restore(cr);
      cairo_stroke(cr);
    }
  }
}

void
draw_cairo_block(cairo_t * cr, pcblock b, int levels)
{
  double    sx, sy, ex, ey;
  uint      rsize, csize;
  double    scale, scalex, scaley;

  cairo_clip_extents(cr, &sx, &sy, &ex, &ey);

  rsize = b->rc->size;
  csize = b->cc->size;

  scalex = (ex - sx) / rsize;
  scaley = (ey - sy) / csize;

  scale = (scalex < scaley ? scalex : scaley);

  cairo_translate(cr, 0.5 * (ex + sx - scale * rsize),
		  0.5 * (ey + sy - scale * csize));
  cairo_scale(cr, scale, scale);
  cairo_set_line_width(cr, 0.5 * cairo_get_line_width(cr) / scale);

  draw_cairo_subblock(cr, b, levels);
}

#endif

/* ------------------------------------------------------------
 Interactive visualization
 ------------------------------------------------------------ */

#ifdef USE_FREEGLUT
static int blocklist = -1;
static float angle_x = 0.0;
static float angle_y = 0.0;
static float start_x = 0.0;
static float start_y = 0.0;
static int base_x = 0;
static int base_y = 0;

static void
display_glut_block()
{
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glRotatef(angle_y, 1.0, 0.0, 0.0);
  glRotatef(angle_x, 0.0, 0.0, 1.0);

  glDisable(GL_LIGHTING);
  glColor3f(1.0, 1.0, 1.0);
  glBegin(GL_LINE_LOOP);
  glVertex3f(-1.0, -1.0, 0.0);
  glVertex3f(1.0, -1.0, 0.0);
  glVertex3f(1.0, 1.0, 0.0);
  glVertex3f(-1.0, 1.0, 0.0);
  glEnd();

  if (blocklist != -1)
    glCallList(blocklist);

  glPopMatrix();
  glFlush();
  glutSwapBuffers();
}

static void
reshape_glut_block(int width, int height)
{
  glViewport(0, 0, width, height);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  gluLookAt(0.0, 0.0, 40.0,	/* My position */
	    0.0, 0.0, 0.0,	/* Where I am looking */
	    0.0, 1.0, 0.0);	/* Which direction is up */
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  if (width > height) {
    glFrustum(-1.5 * width / height, 1.5 * width / height, -1.5, 1.5, 38.5,
	      41.5);
  }
  else {
    glFrustum(-1.5, 1.5, -1.5 * height / width, 1.5 * height / width, 38.5,
	      41.5);
  }
  glMatrixMode(GL_MODELVIEW);
  glEnable(GL_LIGHT0);
  glEnable(GL_DEPTH_TEST);
  glCullFace(GL_BACK);
  glEnable(GL_CULL_FACE);
  glutPostRedisplay();
}

static void
mouse_glut_block(int button, int state, int x, int y)
{
  if (state == GLUT_DOWN) {
    start_x = angle_x;
    start_y = angle_y;
    base_x = x;
    base_y = y;
  }
}

static void
motion_glut_block(int x, int y)
{
  angle_x = start_x + (x - base_x);
  angle_y = start_y + (y - base_y);

  glutPostRedisplay();
}

static void
keyboard_glut_block(unsigned char key, int x, int y)
{
  switch (key) {
  case 'q':
    glutLeaveMainLoop();
    break;
  }
}

static void
displaylist_block(pcblock b, uint roff, uint coff, uint level,
		  pcblock left, pcblock right, pcblock down, pcblock up)
{
  GLfloat   mat_wall[] = {
    0.3, 0.3, 0.3, 1.0
  };
  GLfloat   mat_adm[] = {
    0.0, 1.0, 0.0, 1.0
  };
  GLfloat   mat_inadm[] = {
    1.0, 0.0, 0.0, 1.0
  };
  pcblock   up1, down1, left1, right1;
  uint      rows = b->rc->size;
  uint      cols = b->cc->size;
  uint      roff1, coff1;
  uint      i, j;

  if (level > 0) {
    glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_wall);
    if (!left) {
      glBegin(GL_TRIANGLE_STRIP);
      glNormal3f(1.0, 0.0, 0.0);
      glVertex3f(coff, roff, -0.05 * level);
      glVertex3f(coff, roff + rows, -0.05 * level);
      glVertex3f(coff, roff, -0.05 * (level - 1));
      glVertex3f(coff, roff + rows, -0.05 * (level - 1));
      glEnd();
    }
    if (!right) {
      glBegin(GL_TRIANGLE_STRIP);
      glNormal3f(-1.0, 0.0, 0.0);
      glVertex3f(coff + cols, roff + rows, -0.05 * level);
      glVertex3f(coff + cols, roff, -0.05 * level);
      glVertex3f(coff + cols, roff + rows, -0.05 * (level - 1));
      glVertex3f(coff + cols, roff, -0.05 * (level - 1));
      glEnd();
    }
    if (!down) {
      glBegin(GL_TRIANGLE_STRIP);
      glNormal3f(0.0, 1.0, 0.0);
      glVertex3f(coff + cols, roff, -0.05 * level);
      glVertex3f(coff, roff, -0.05 * level);
      glVertex3f(coff + cols, roff, -0.05 * (level - 1));
      glVertex3f(coff, roff, -0.05 * (level - 1));
      glEnd();
    }
    if (!up) {
      glBegin(GL_TRIANGLE_STRIP);
      glNormal3f(0.0, -1.0, 0.0);
      glVertex3f(coff, roff + rows, -0.05 * level);
      glVertex3f(coff + cols, roff + rows, -0.05 * level);
      glVertex3f(coff, roff + rows, -0.05 * (level - 1));
      glVertex3f(coff + cols, roff + rows, -0.05 * (level - 1));
      glEnd();
    }
  }

  if (b->son) {
    coff1 = coff;
    for (j = 0; j < b->csons; j++) {
      roff1 = roff;
      for (i = 0; i < b->rsons; i++) {
	left1 = 0;
	if (j == 0) {
	  if (left && left->son) {
	    assert(left->rsons == b->rsons);
	    left1 = left->son[i + (left->csons - 1) * left->rsons];
	  }
	}
	else
	  left1 = b->son[i + (j - 1) * b->rsons];

	right1 = 0;
	if (j + 1 < b->csons)
	  right1 = b->son[i + (j + 1) * b->rsons];
	else {
	  if (right && right->son) {
	    assert(right->rsons == b->rsons);
	    right1 = right->son[i];
	  }
	}

	down1 = 0;
	if (i == 0) {
	  if (down && down->son) {
	    assert(down->csons == b->csons);
	    down1 = down->son[(down->rsons - 1) + j * down->rsons];
	  }
	}
	else
	  down1 = b->son[(i - 1) + j * b->rsons];

	up1 = 0;
	if (i + 1 < b->rsons)
	  up1 = b->son[(i + 1) + j * b->rsons];
	else {
	  if (up && up->son) {
	    assert(up->csons == b->csons);
	    up1 = up->son[j * up->rsons];
	  }
	}

	displaylist_block(b->son[i + j * b->rsons], roff1, coff1, level + 1,
			  left1, right1, down1, up1);
	roff1 += b->son[i + j * b->rsons]->rc->size;
      }
      assert(roff1 == roff + b->rc->size);
      coff1 += b->son[j * b->rsons]->cc->size;
    }
    assert(coff1 == coff + b->cc->size);
  }
  else {
    if (b->a)
      glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_adm);
    else
      glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_inadm);
    glNormal3f(0.0, 0.0, 1.0);
    glBegin(GL_TRIANGLE_STRIP);
    glVertex3f(coff, roff, -0.05 * level);
    glVertex3f(coff + cols, roff, -0.05 * level);
    glVertex3f(coff, roff + rows, -0.05 * level);
    glVertex3f(coff + cols, roff + rows, -0.05 * level);
    glEnd();
  }
}

void
view_block(pcblock b)
{
  uint      rows = b->rc->size;
  uint      cols = b->cc->size;

  glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE, 1);

  glutInitDisplayMode(GLUT_DEPTH | GLUT_RGBA | GLUT_DOUBLE);
  glutInitWindowSize(800, 600);
  glutCreateWindow("View Block Tree");
  glClearColor(0.0, 0.0, 0.0, 0.0);

  if (blocklist == -1)
    blocklist = glGenLists(1);

  glNewList(blocklist, GL_COMPILE);
  glEnable(GL_LIGHTING);
  glEnable(GL_NORMALIZE);
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  if (rows > cols) {
    glTranslatef(-1.0, -1.0 * cols / rows, 0.0);
    glScalef(2.0 / rows, 2.0 / rows, 1.0);
  }
  else {
    glTranslatef(-1.0 * rows / cols, -1.0, 0.0);
    glScalef(2.0 / cols, 2.0 / cols, 1.0);
  }
  displaylist_block(b, 0, 0, 0, 0, 0, 0, 0);
  glPopMatrix();
  glEndList();

  glutDisplayFunc(display_glut_block);
  glutReshapeFunc(reshape_glut_block);
  glutMouseFunc(mouse_glut_block);
  glutMotionFunc(motion_glut_block);
  glutKeyboardFunc(keyboard_glut_block);

  glutMainLoop();
}
#else
void
view_block(pcblock b)
{
  (void) b;
  (void) fprintf(stderr, "Sorry, GLUT is not available.\n");
}
#endif

/* ------------------------------------------------------------
 Hierarchical iterators
 ------------------------------------------------------------ */

void
iterate_block(pcblock b, uint bname, uint rname, uint cname,
	      void (*pre) (pcblock b, uint bname, uint rname, uint cname,
			   uint pardepth, void *data),
	      void (*post) (pcblock b, uint bname, uint rname, uint cname,
			    uint pardepth, void *data), void *data)
{
  pcblock   b1;
  uint      bname1, rname1, cname1;
  uint      i, j;

  if (pre)
    pre(b, bname, rname, cname, 0, data);

  if (b->son) {
    bname1 = bname + 1;
    cname1 = (b->son[0]->cc == b->cc ? cname : cname + 1);

    for (j = 0; j < b->csons; j++) {
      rname1 = (b->son[0]->rc == b->rc ? rname : rname + 1);

      for (i = 0; i < b->rsons; i++) {
	b1 = b->son[i + j * b->rsons];

	iterate_block(b1, bname1, rname1, cname1, pre, post, data);

	bname1 += b1->desc;
	rname1 += b1->rc->desc;
      }
      assert(rname1 == rname + b->rc->desc);

      cname1 += b->son[j * b->rsons]->cc->desc;
    }
    assert(cname1 == cname + b->cc->desc);
    assert(bname1 == bname + b->desc);
  }

  if (post)
    post(b, bname, rname, cname, 0, data);
}

static void
del_blockentry(pblockentry be)
{
  pblockentry next;

  while (be) {
    next = be->next;
    freemem(be);
    be = next;
  }
}

static    pblockentry
addrow(pcblock b, uint bname, uint rname, uint cname,
       pcblockentry father, pblockentry next)
{
  pblockentry pb;
  pcblockentry pbf;
  pccluster cc;
  uint      bname1, cname1;
  uint      j;

  pb = (pblockentry) allocmem(sizeof(blockentry));
  pb->b = b;
  pb->bname = bname;
  pb->rname = rname;
  pb->cname = cname;
  pb->father = father;
  pb->next = next;

  pbf = pb;

  if (b->son && b->son[0]->rc == b->rc) {
    cc = b->cc;

    assert(b->csons == cc->sons);

    bname1 = bname + 1;
    cname1 = cname + 1;

    for (j = 0; j < b->csons; j++) {
      pb = addrow(b->son[j * b->rsons], bname1, rname, cname1, pbf, pb);

      bname1 += b->son[j]->desc;
      cname1 += cc->son[j]->desc;
    }
    assert(bname1 == bname + b->desc);
    assert(cname1 == cname + cc->desc);
  }

  return pb;
}

static void
iterate_rowlist(pccluster rc, uint rname, pblockentry pb,
		uint pardepth, void (*pre) (pcblockentry pb, uint pardepth,
					    void *data),
		void (*post) (pcblockentry pb, uint pardepth, void *data),
		void *data)
{
  pblockentry pb0, *pb1;
  pcblock   b;
  uint      bname1, rname1, cname1;
  uint     *rnames;
#ifdef USE_OPENMP
  uint      nthreads;		/* HACK: Solaris workaround */
#endif
  uint      i, j;

  /* Quick exit */
  if (pb == 0)
    return;

  /* Use "pre" callback, if one is given */
  if (pre)
    pre(pb, pardepth, data);

  /* Consider sons */
  if (rc->son) {
    /* Initialize block lists for all sons */
    pb1 = (pblockentry *) allocmem((size_t) sizeof(pblockentry) * rc->sons);
    for (i = 0; i < rc->sons; i++)
      pb1[i] = 0;

    /* Initialize names for all sons */
    rnames = (uint *) allocmem((size_t) sizeof(uint) * rc->sons);
    rname1 = rname + 1;
    for (i = 0; i < rc->sons; i++) {
      rnames[i] = rname1;
      rname1 += rc->son[i]->desc;
    }
    assert(rname1 == rname + rc->desc);

    /* Fill block lists for all sons */
    for (pb0 = pb; pb0; pb0 = pb0->next) {
      b = pb0->b;

      assert(b->rc == rc);

      /* Consider block sons */
      if (b->son && b->son[0]->rc != rc) {
	assert(b->rsons == rc->sons);
	assert(pb0->rname == rname);

	bname1 = pb0->bname + 1;

	cname1 = (b->son[0]->cc == b->cc ? pb0->cname : pb0->cname + 1);
	for (j = 0; j < b->csons; j++) {
	  rname1 = rname + 1;
	  for (i = 0; i < b->rsons; i++) {
	    assert(b->son[i + j * b->rsons]->rc == rc->son[i]);

	    pb1[i] = addrow(b->son[i + j * b->rsons], bname1, rname1, cname1,
			    pb0, pb1[i]);

	    bname1 += b->son[i + j * b->rsons]->desc;
	    rname1 += rc->son[i]->desc;
	  }
	  assert(rname1 == rname + rc->desc);

	  cname1 += b->son[j * b->rsons]->cc->desc;
	}
	assert(bname1 == pb0->bname + b->desc);
	assert(cname1 == pb0->cname + b->cc->desc);
      }
    }

    /* Recursively handle sons */
#ifdef USE_OPENMP
    nthreads = rc->sons;
    (void) nthreads;
#pragma omp parallel for if(pardepth > 0), num_threads(nthreads)
#endif
    for (i = 0; i < rc->sons; i++)
      iterate_rowlist(rc->son[i], rnames[i], pb1[i],
		      (pardepth > 0 ? pardepth - 1 : 0), pre, post, data);

    /* Clean up */
    for (i = 0; i < rc->sons; i++)
      del_blockentry(pb1[i]);
    freemem(pb1);
    freemem(rnames);
  }

  /* Use "post" callback, if one is given */
  if (post)
    post(pb, pardepth, data);
}

void
iterate_rowlist_block(pcblock b, uint bname, uint rname, uint cname,
		      uint pardepth, void (*pre) (pcblockentry pb,
						  uint pardepth, void *data),
		      void (*post) (pcblockentry pb, uint pardepth,
				    void *data), void *data)
{
  pblockentry pb;

  pb = addrow(b, bname, rname, cname, 0, 0);

  iterate_rowlist(b->rc, rname, pb, pardepth, pre, post, data);

  del_blockentry(pb);
}

static    pblockentry
addcol(pcblock b, uint bname, uint rname, uint cname,
       pcblockentry father, pblockentry next)
{
  pblockentry pb;
  pcblockentry pbf;
  pccluster rc;
  uint      bname1, rname1;
  uint      i;

  pb = (pblockentry) allocmem(sizeof(blockentry));
  pb->b = b;
  pb->bname = bname;
  pb->rname = rname;
  pb->cname = cname;
  pb->father = father;
  pb->next = next;

  pbf = pb;

  if (b->son && b->son[0]->cc == b->cc) {
    rc = b->rc;

    assert(b->rsons == rc->sons);

    bname1 = bname + 1;
    rname1 = rname + 1;

    for (i = 0; i < b->rsons; i++) {
      pb = addcol(b->son[i], bname1, rname1, cname, pbf, pb);

      bname1 += b->son[i]->desc;
      rname1 += rc->son[i]->desc;
    }
    assert(bname1 == bname + b->desc);
    assert(rname1 == rname + rc->desc);
  }

  return pb;
}

static void
iterate_collist(pccluster cc, uint cname, pblockentry pb,
		uint pardepth, void (*pre) (pcblockentry pb, uint pardepth,
					    void *data),
		void (*post) (pcblockentry pb, uint pardepth, void *data),
		void *data)
{
  pblockentry pb0, *pb1;
  pcblock   b;
  uint      bname1, rname1, cname1;
  uint     *cnames;
#ifdef USE_OPENMP
  uint      nthreads;		/* HACK: Solaris workaround */
#endif
  uint      i, j;

  /* Quick exit */
  if (pb == 0)
    return;

  /* Use "pre" callback, if one is given */
  if (pre)
    pre(pb, pardepth, data);

  /* Consider sons */
  if (cc->son) {
    /* Initialize block lists for all sons */
    pb1 = (pblockentry *) allocmem((size_t) sizeof(pblockentry) * cc->sons);
    for (j = 0; j < cc->sons; j++)
      pb1[j] = 0;

    /* Initialize names for all sons */
    cnames = (uint *) allocmem((size_t) sizeof(uint) * cc->sons);
    cname1 = cname + 1;
    for (j = 0; j < cc->sons; j++) {
      cnames[j] = cname1;
      cname1 += cc->son[j]->desc;
    }
    assert(cname1 == cname + cc->desc);

    /* Fill block lists for all sons */
    for (pb0 = pb; pb0; pb0 = pb0->next) {
      b = pb0->b;

      assert(b->cc == cc);

      /* Consider block sons */
      if (b->son && b->son[0]->cc != cc) {
	assert(b->csons == cc->sons);
	assert(pb0->cname == cname);

	bname1 = pb0->bname + 1;

	cname1 = cname + 1;
	for (j = 0; j < b->csons; j++) {
	  rname1 = (b->son[0]->rc == b->rc ? pb0->rname : pb0->rname + 1);

	  for (i = 0; i < b->rsons; i++) {
	    assert(b->son[i + j * b->rsons]->cc == cc->son[j]);

	    pb1[j] = addcol(b->son[i + j * b->rsons], bname1, rname1, cname1,
			    pb0, pb1[j]);

	    bname1 += b->son[i + j * b->rsons]->desc;
	    rname1 += b->son[i + j * b->rsons]->rc->desc;
	  }
	  assert(rname1 == pb0->rname + b->rc->desc);

	  cname1 += b->son[j * b->rsons]->cc->desc;
	}
	assert(bname1 == pb0->bname + b->desc);
	assert(cname1 == cname + b->cc->desc);
      }
    }

    /* Recursively handle sons */
#ifdef USE_OPENMP
    nthreads = cc->sons;
    (void) nthreads;
#pragma omp parallel for if(pardepth > 0), num_threads(nthreads)
#endif
    for (j = 0; j < cc->sons; j++)
      iterate_collist(cc->son[j], cnames[j], pb1[j],
		      (pardepth > 0 ? pardepth - 1 : 0), pre, post, data);

    /* Clean up */
    for (j = 0; j < cc->sons; j++)
      del_blockentry(pb1[j]);
    freemem(pb1);
    freemem(cnames);
  }

  /* Use "post" callback, if one is given */
  if (post)
    post(pb, pardepth, data);
}

void
iterate_collist_block(pcblock b, uint bname, uint rname, uint cname,
		      uint pardepth, void (*pre) (pcblockentry pb,
						  uint pardepth, void *data),
		      void (*post) (pcblockentry pb, uint pardepth,
				    void *data), void *data)
{
  pblockentry pb;

  pb = addcol(b, bname, rname, cname, 0, 0);

  iterate_collist(b->cc, cname, pb, pardepth, pre, post, data);

  del_blockentry(pb);
}

struct _listdata {
  void      (*pre) (pcblock b, uint bname, uint rname, uint cname,
		    uint pardepth, void *data);
  void      (*post) (pcblock b, uint bname, uint rname, uint cname,
		     uint pardepth, void *data);
  void     *data;
};

static void
call_reversed_blockentry(pcblockentry pb, uint pardepth,
			 void (*pre) (pcblock b, uint bname, uint rname,
				      uint cname, uint pardepth, void *data),
			 void *data)
{
  if (pb->next)
    call_reversed_blockentry(pb->next, pardepth, pre, data);

  pre(pb->b, pb->bname, pb->rname, pb->cname, pardepth, data);
}

static void
pre_blocklist(pcblockentry pb, uint pardepth, void *data)
{
  struct _listdata *ld = (struct _listdata *) data;

  if (ld->pre)
    call_reversed_blockentry(pb, pardepth, ld->pre, ld->data);
}

static void
call_inorder_blockentry(pcblockentry pb, uint pardepth,
			void (*post) (pcblock b, uint bname, uint rname,
				      uint cname, uint pardepth, void *data),
			void *data)
{
  post(pb->b, pb->bname, pb->rname, pb->cname, pardepth, data);

  if (pb->next)
    call_inorder_blockentry(pb->next, pardepth, post, data);
}

static void
post_blocklist(pcblockentry pb, uint pardepth, void *data)
{
  struct _listdata *ld = (struct _listdata *) data;

  if (ld->post)
    call_inorder_blockentry(pb, pardepth, ld->post, ld->data);
}

void
iterate_byrow_block(pcblock b, uint bname, uint rname, uint cname,
		    uint pardepth,
		    void (*pre) (pcblock b, uint bname, uint rname,
				 uint cname, uint pardepth, void *data),
		    void (*post) (pcblock b, uint bname, uint rname,
				  uint cname, uint pardepth, void *data),
		    void *data)
{
  struct _listdata ld;
  pblockentry pb;

  ld.pre = pre;
  ld.post = post;
  ld.data = data;

  pb = addrow(b, bname, rname, cname, 0, 0);

  iterate_rowlist(b->rc, rname, pb, pardepth, pre_blocklist, post_blocklist,
		  &ld);

  del_blockentry(pb);
}

void
iterate_bycol_block(pcblock b, uint bname, uint rname, uint cname,
		    uint pardepth,
		    void (*pre) (pcblock b, uint bname, uint rname,
				 uint cname, uint pardepth, void *data),
		    void (*post) (pcblock b, uint bname, uint rname,
				  uint cname, uint pardepth, void *data),
		    void *data)
{
  struct _listdata ld;
  pblockentry pb;

  ld.pre = pre;
  ld.post = post;
  ld.data = data;

  pb = addcol(b, bname, rname, cname, 0, 0);

  iterate_collist(b->cc, cname, pb, pardepth, pre_blocklist, post_blocklist,
		  &ld);

  del_blockentry(pb);
}

/* ------------------------------------------------------------
 Enumeration
 ------------------------------------------------------------ */

static void
enumerate(pblock b, uint bname, pblock *bn)
{
  uint      bname1;
  uint      i, j;

  bn[bname] = b;

  bname1 = bname + 1;
  for (j = 0; j < b->csons; j++)
    for (i = 0; i < b->rsons; i++) {
      enumerate(b->son[i + j * b->rsons], bname1, bn);

      bname1 += b->son[i + j * b->rsons]->desc;
    }
  assert(bname1 == bname + b->desc);
}

pblock   *
enumerate_block(pblock b)
{
  pblock   *bn;

  bn = (pblock *) allocmem((size_t) sizeof(pblock) * b->desc);

  enumerate(b, 0, bn);

  return bn;
}

static void
enumerate_level(pblock b, uint bname, uint * bln, uint p)
{
  uint      bname1;
  uint      i, j;

  bln[bname] = p + 1;

  bname1 = bname + 1;
  for (j = 0; j < b->csons; j++)
    for (i = 0; i < b->rsons; i++) {
      enumerate_level(b->son[i + j * b->rsons], bname1, bln, p + 1);

      bname1 += b->son[i + j * b->rsons]->desc;
    }
  assert(bname1 == bname + b->desc);
}

uint     *
enumerate_level_block(pblock b)
{
  uint     *bln;

  bln = allocuint(b->desc);

  enumerate_level(b, 0, bln, 0);

  return bln;
}

uint
getdepth_block(pcblock b)
{
  uint      rsons = b->rsons;
  uint      csons = b->csons;

  uint      p, i, j;

  p = 0;

  if (b->son != NULL) {
    for (j = 0; j < csons; ++j) {
      for (i = 0; i < rsons; ++i) {
	p = UINT_MAX(p, getdepth_block(b->son[i + j * rsons]));
      }
    }
    p++;
  }

  return p;
}

static void
pre_compute_csp_block(pcblock b, uint bname, uint rname, uint cname,
		      uint pardepth, void *data)
{
  uint     *row, *col;
  uint    **cspdata;

  (void) b;
  (void) bname;
  (void) pardepth;

  cspdata = (uint **) data;
  row = cspdata[0];
  col = cspdata[1];

  row[rname]++;
  col[cname]++;
}

uint
compute_csp_block(pcblock b)
{
  uint      i, csp;
  uint     *row, *col;
  uint    **cspdata;

  row = (uint *) allocmem((size_t) sizeof(uint) * b->rc->desc);
  col = (uint *) allocmem((size_t) sizeof(uint) * b->cc->desc);
  cspdata = (uint **) allocmem((size_t) sizeof(uint *) * 2);
  cspdata[0] = row;
  cspdata[1] = col;

  for (i = 0; i < b->rc->desc; i++)
    row[i] = 0;
  for (i = 0; i < b->cc->desc; i++)
    col[i] = 0;

  iterate_block(b, 0, 0, 0, pre_compute_csp_block, NULL, cspdata);

  csp = 0;
  for (i = 0; i < b->rc->desc; i++)
    csp = (row[i] > csp ? row[i] : csp);
  for (i = 0; i < b->cc->desc; i++)
    csp = (col[i] > csp ? col[i] : csp);

  freemem(row);
  freemem(col);
  freemem(cspdata);
  return csp;
}
