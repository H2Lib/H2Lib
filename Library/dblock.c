
/* ------------------------------------------------------------
 * This is the file "dblock.c" of the H2Lib package.
 * All rights reserved, Steffen Boerm 2015
 * ------------------------------------------------------------ */

#include "dblock.h"

#include "basic.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>

static uint active_dblock = 0;

/* ------------------------------------------------------------
 * Constructors and destructors
 * ------------------------------------------------------------ */

pdblock
new_dblock(pdcluster rc, pdcluster cc, uint rd, uint cd,
	   uint rsons, uint csons)
{
  pdblock   b;
  uint      i, j;

  b = (pdblock) allocmem(sizeof(dblock));
  b->rc = rc;
  b->cc = cc;
  b->rd = rd;
  b->cd = cd;
  b->adm = false;
  b->rsons = rsons;
  b->csons = csons;
  b->desc = 0;

  b->son = NULL;
  if (rsons > 0 && csons > 0) {
    b->son = (pdblock *) allocmem((size_t) sizeof(dblock) * rsons * csons);
    for (j = 0; j < csons; j++)
      for (i = 0; i < rsons; i++)
	b->son[i + j * rsons] = NULL;
  }

#ifdef USE_OPENMP
#pragma omp atomic
#endif
  active_dblock++;

  return b;
}

void
update_dblock(pdblock b)
{
  uint      desc;
  uint      i, j;

  desc = 1;
  for (j = 0; j < b->csons; j++)
    for (i = 0; i < b->rsons; i++)
      desc += b->son[i + j * b->rsons]->desc;

  b->desc = desc;
}

void
del_dblock(pdblock b)
{
  uint      i, j;

  if (b->son) {
    for (j = 0; j < b->csons; j++)
      for (i = 0; i < b->rsons; i++)
	del_dblock(b->son[i + j * b->rsons]);
    freemem(b->son);
  }

  freemem(b);

  assert(active_dblock > 0);

#ifdef USE_OPENMP
#pragma omp atomic
#endif
  active_dblock--;
}

/* ------------------------------------------------------------
 * Statistics
 * ------------------------------------------------------------ */

uint
getactives_dblock()
{
  return active_dblock;
}

size_t
getsize_dblock(pcdblock b)
{
  size_t    sz;
  uint      i, j;

  sz = (size_t) sizeof(dblock);
  if (b->son) {
    sz += (size_t) sizeof(pdblock) * b->rsons * b->csons;
    for (j = 0; j < b->csons; j++)
      for (i = 0; i < b->rsons; i++)
	sz += getsize_dblock(b->son[i + j * b->rsons]);
  }

  return sz;
}

uint
getdepth_dblock(pcdblock b)
{
  uint      depth, depth1;
  uint      i, j;

  depth = 0;
  if (b->son) {
    for (j = 0; j < b->csons; j++)
      for (i = 0; i < b->rsons; i++) {
	depth1 = getdepth_dblock(b->son[i + j * b->rsons]) + 1;
	if (depth1 > depth)
	  depth = depth1;
      }
  }
  return depth;
}

/* ------------------------------------------------------------
 * Drawing
 * ------------------------------------------------------------ */

#ifdef USE_CAIRO
static void
cairodraw_subblock(cairo_t * cr, pcdblock b, int levels)
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
	cairodraw_subblock(cr, b->son[i + j * rsons], levels - 1);
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
      cairo_rectangle(cr, 0.0, 0.0, csize, rsize);
      cairo_stroke(cr);
    }
    else if (b->adm) {
      cairo_rectangle(cr, 0.0, 0.0, csize, rsize);
      cairo_save(cr);
      cairo_set_source_rgb(cr, 0.2, 0.2, 1.0);
      cairo_fill_preserve(cr);
      cairo_restore(cr);
      cairo_stroke(cr);
    }
    else {
      cairo_rectangle(cr, 0.0, 0.0, csize, rsize);
      cairo_save(cr);
      cairo_set_source_rgb(cr, 1.0, 0.0, 0.0);
      cairo_fill_preserve(cr);
      cairo_restore(cr);
      cairo_stroke(cr);
    }
  }
}

void
cairodraw_dblock(cairo_t * cr, pcdblock b, int levels)
{
  double    sx, sy, ex, ey;
  uint      rsize, csize;
  double    scalex, scaley, scale;

  /* Save Cairo state */
  cairo_save(cr);

  /* Obtain size of block */
  rsize = b->rc->size;
  csize = b->cc->size;

  /* Obtain size of current Cairo bounding box */
  cairo_clip_extents(cr, &sx, &sy, &ex, &ey);

  /* Compute scaling factor */
  scalex = (ex - sx) / rsize;
  scaley = (ey - sy) / csize;
  scale = (scalex < scaley ? scalex : scaley);

  /* Center block in bounding box */
  cairo_translate(cr,
		  0.5 * (ex - sx - scale * rsize),
		  0.5 * (ey - sy - scale * csize));

  /* Scale coordinates */
  cairo_scale(cr, scale, scale);
  cairo_set_line_width(cr, cairo_get_line_width(cr) / scale);

  /* Start drawing */
  cairodraw_subblock(cr, b, levels);

  /* Restore Cairo state */
  cairo_restore(cr);
}
#endif

/* ------------------------------------------------------------
 * Building standard block trees
 * ------------------------------------------------------------ */

pdblock
build_dblock(pdcluster rc, pdcluster cc, uint l,
	     bool(*admissible) (pdcluster rc, pdcluster cc, uint l,
				uint * rd, uint * cd, void *data), void *data)
{
  pdblock   b;
  uint      rd, cd;
  uint      i, j;

  b = 0;
  if (admissible(rc, cc, l, &rd, &cd, data)) {
    b = new_dblock(rc, cc, rd, cd, 0, 0);
    b->adm = true;
  }
  else if (rc->sons == 0 || cc->sons == 0) {
    b = new_dblock(rc, cc, 0, 0, 0, 0);
  }
  else {
    assert(rc->sons > 0);
    assert(cc->sons > 0);

    b = new_dblock(rc, cc, 0, 0, rc->sons, cc->sons);
    for (j = 0; j < cc->sons; j++)
      for (i = 0; i < rc->sons; i++)
	b->son[i + j * rc->sons] = build_dblock(rc->son[i], cc->son[j], l + 1,
						admissible, data);
  }

  update_dblock(b);

  return b;
}

bool
parabolic_admissibility(pdcluster rc, pdcluster cc, uint l,
			uint * rd, uint * cd, void *data)
{
  pcdiradmdata dad = (pcdiradmdata) data;
  real      eta2 = dad->eta2;
  real      wave_k = dad->wave_k;
  preal     xy = dad->xy;
  real      rdiam, cdiam, dist, dist2, invdist;
  uint      dim;
  uint      i, iota;

  dim = rc->dim;
  assert(cc->dim == dim);

  /* Compute cluster diameters */
  rdiam = diam_dcluster(rc);
  cdiam = diam_dcluster(cc);

  /* Compute distance */
  dist = dist_dcluster(rc, cc);

  /* Check parabolic admissibility condition */
  if (wave_k * rdiam * rdiam <= eta2 * dist &&
      wave_k * cdiam * cdiam <= eta2 * dist &&
      rdiam <= eta2 * dist && cdiam <= eta2 * dist) {

    /* Compute distance of midpoints */
    dist2 = 0.0;
    for (i = 0; i < dim; i++) {
      xy[i] = 0.5 * (rc->bmax[i] + rc->bmin[i]
		     - cc->bmax[i] - cc->bmin[i]);
      dist2 += REAL_SQR(xy[i]);
    }
    invdist = 1.0 / REAL_SQRT(dist2);

    /* Scale distance vector... */
    for (i = 0; i < dim; i++)
      xy[i] *= invdist;

    /* ... and choose a direction that is sufficiently close */
    iota = finddirection_leveldir(dad->ld, l, 1.0, xy);
    *rd = iota;
    *cd = iota;

    return true;
  }

  *rd = 0;
  *cd = 0;

  return false;
}

bool
standard_admissibility(pdcluster rc, pdcluster cc, uint l,
		       uint * rd, uint * cd, void *data)
{
  pcdiradmdata dad = (pcdiradmdata) data;
  real      eta2 = dad->eta2;
  preal     xy = dad->xy;
  real      rdiam, cdiam, dist, dist2, invdist;
  uint      dim;
  uint      i, iota;

  dim = rc->dim;
  assert(cc->dim == dim);

  /* Compute cluster diameters */
  rdiam = diam_dcluster(rc);
  cdiam = diam_dcluster(cc);

  /* Compute distance */
  dist = dist_dcluster(rc, cc);

  /* Check standard admissibility condition */
  if (rdiam <= eta2 * dist && cdiam <= eta2 * dist) {

    /* Compute distance of midpoints */
    dist2 = 0.0;
    for (i = 0; i < dim; i++) {
      xy[i] = 0.5 * (rc->bmax[i] + rc->bmin[i]
		     - cc->bmax[i] - cc->bmin[i]);
      dist2 += REAL_SQR(xy[i]);
    }
    invdist = 1.0 / REAL_SQRT(dist2);

    /* Scale distance vector... */
    for (i = 0; i < dim; i++)
      xy[i] *= invdist;

    /* ... and choose a direction that is sufficiently close */
    iota = finddirection_leveldir(dad->ld, l, 1.0, xy);
    *rd = iota;
    *cd = iota;

    return true;
  }

  *rd = 0;
  *cd = 0;

  return false;
}

real
getmaxeta_dblock(pcdblock b)
{
  real      eta, eta1, rdiam, cdiam, dist;
  uint      rsons, csons;
  uint      i, j;

  eta = 0.0;

  if (b->son) {
    rsons = b->rsons;
    csons = b->csons;

    for (j = 0; j < csons; j++)
      for (i = 0; i < rsons; i++) {
	eta1 = getmaxeta_dblock(b->son[i + j * rsons]);

	if (eta1 > eta)
	  eta = eta1;
      }
  }
  else if (b->adm) {
    rdiam = diam_dcluster(b->rc);
    cdiam = diam_dcluster(b->cc);
    dist = dist_dcluster(b->rc, b->cc);
    eta = REAL_MAX(rdiam, cdiam) / dist;
  }

  return eta;
}

static void
scan_tree(pdblock b, uint ** used, uint ** dirson, uint * fill, uint count)
{

  uint      rsons, csons;
  uint      i, j;
  uint      directions;


  rsons = b->rsons;
  csons = b->csons;

  directions = b->rc->directions;
  assert(directions == b->cc->directions);

  if ((rsons + csons) > 0) {	/* We have sons */
    if (fill[count] == 0) {	/* Copy important information */
      fill[count] = 1;

      if (b->rc->sons > 0) {
	for (i = 0; i < directions; i++) {
	  dirson[count][i] = b->rc->dirson[0][i];
	}
      }
      else {
	for (i = 0; i < directions; i++) {
	  dirson[count][i] = b->cc->dirson[0][i];
	}
      }
    }
    if (directions > 1) {	/* We always have a symmetric direction object -> more than one direction */
      used[count][b->rd] = 1;
      used[count][b->cd] = 1;
    }
    count++;

    /* Recursiv call */
    if (rsons > 0) {
      if (csons > 0) {
	for (i = 0; i < csons; i++) {
	  for (j = 0; j < rsons; j++) {
	    scan_tree(b->son[j + i * rsons], used, dirson, fill, count);
	  }
	}
      }
      else {
	for (j = 0; j < rsons; j++) {
	  scan_tree(b->son[j], used, dirson, fill, count);
	}
      }
    }
    else {
      for (i = 0; i < rsons; i++) {
	scan_tree(b->son[i], used, dirson, fill, count);
      }
    }
  }
  else {			/* Leaf case */
    /* Save used directions */
    if (directions > 1) {
      used[count][b->rd] = 1;
      used[count][b->cd] = 1;
    }
  }

}


static void
check_directions(uint ** used, uint ** tmp_dirson, pleveldir ld)
{

  uint      i, j;
  uint      l;

  i = 0;

  /* Starts testing with root if for every direction the compatible son direction is inside
     and if not add it! */
  while (i < ld->depth) {
    for (j = 0; j < ld->directions[i]; j++) {
      if (used[i][j] > 0) {
	l = tmp_dirson[i][j];
	/* If it is one of the used, everything is all right, if not it will be fixed */
	if (used[i + 1][l] == 0) {
	  //printf("correct \n");          
	  used[i + 1][l] = 1;
	}
      }
    }
    i++;
  }

}


static void
copy_direction(uint ** used, uint ** idx, pleveldir ldir, pleveldir lold)
{

  uint      i, j, l;
  uint      entrys;

  i = 0;
  while (i < ldir->depth + 1) {
    entrys = 0;
    for (j = 0; j < lold->directions[i]; j++) {
      if (used[i][j] > 0) {
	idx[i][j] = entrys;
	for (l = 0; l < ldir->dim; l++) {
	  ldir->dir[i][entrys][l] = lold->dir[i][j][l];
	}
	entrys++;
      }
    }
    assert(entrys == ldir->directions[i]);
    i++;
  }

}


static void
change_dcluster(pdcluster t, uint ** used, uint ** idx, uint ** tmp_dirson,
		pleveldir ldir, uint count)
{

  uint      i, j, m;
  uint      entrys;

  if (t->sons) {

    count++;
    for (i = 0; i < t->sons; i++) {
      entrys = 0;
      change_dcluster(t->son[i], used, idx, tmp_dirson, ldir, count);
      if (t->dirson) {
	freemem(t->dirson[i]);
	t->dirson[i] =
	  (uint *) allocmem(sizeof(uint) * ldir->directions[count - 1]);
	for (j = 0; j < t->directions; j++) {
	  if (used[count - 1][j] > 0) {
	    if (ldir->directions[count] > 0) {
	      m = tmp_dirson[count - 1][j];
	      t->dirson[i][entrys] = idx[count][m];
	      entrys++;
	    }
	    else {		/* Son hasn't directions */
	      t->dirson[i][entrys] = 0;
	      entrys++;
	    }
	  }
	}

	assert(entrys == ldir->directions[count - 1]);
      }
    }
    count--;
    t->directions = ldir->directions[count];
    t->dir = (pcreal *) ldir->dir[count];

  }
  else {			/* Leaf case */
    t->directions = ldir->directions[count];
    t->dir = (pcreal *) ldir->dir[count];
  }
}


static void
change_dblock(pdblock b, uint ** idx, uint count)
{

  uint      rsons, csons;
  uint      i, j;

  rsons = b->rsons;
  csons = b->csons;

  if ((rsons + csons) > 0) {	/* We have sons */
    count++;
    for (j = 0; j < b->csons; j++) {
      for (i = 0; i < b->rsons; i++) {
	change_dblock(b->son[i + j * rsons], idx, count);
      }
    }
    count--;
    if (b->rc->directions > 0) {
      b->rd = idx[count][b->rd];
      b->cd = idx[count][b->cd];
    }
    else {
      b->rd = 0;
      b->cd = 0;
    }
  }
  else {			/* Leaf case */
    if (b->rc->directions > 0) {
      b->rd = idx[count][b->rd];
      b->cd = idx[count][b->cd];
    }
    else {
      b->rd = 0;
      b->cd = 0;
    }
  }

}


pleveldir
remove_unused_direction(pdblock b, pdcluster t, pleveldir lold)
{

  uint      i, j;
  uint      entrys, dirs;
  preal     dirmem;
  pleveldir ldir;

  if (t->directions == 0) {
    return lold;
  }
  ldir = new_leveldir(getdepth_dcluster(t), t->dim);

  uint      depth = getdepth_dblock(b);
  uint    **idx = (uint **) allocmem(sizeof(uint *) * (depth + 1));
  uint    **used = (uint **) allocmem(sizeof(uint *) * (depth + 1));
  uint    **tmp_dirson = (uint **) allocmem(sizeof(uint *) * (depth));
  uint     *fill = (uint *) allocmem(sizeof(uint) * (depth + 1));

  /* Copy Information that won't change */
  for (i = 0; i < ldir->depth + 1; i++) {
    ldir->maxdiam[i] = lold->maxdiam[i];
    ldir->splits[i] = lold->splits[i];
  }

  /* Set up used array */
  for (i = 0; i < depth; i++) {
    used[i] = (uint *) allocmem(sizeof(uint) * lold->directions[i]);
    tmp_dirson[i] = (uint *) allocmem(sizeof(uint) * lold->directions[i]);
    fill[i] = 0;
    for (j = 0; j < lold->directions[i]; j++) {
      used[i][j] = 0;
    }

  }
  used[depth] = (uint *) allocmem(sizeof(uint) * lold->directions[depth]);
  fill[depth] = 0;
  for (j = 0; j < lold->directions[depth]; j++) {
    used[depth][j] = 0;
  }

  /* Scan directional block tree for used directions */
  scan_tree(b, used, tmp_dirson, fill, 0);
  freemem(fill);

  /* Check if all needed directions are inside and copy important directions */
  check_directions(used, tmp_dirson, lold);

  /* Set up memory for new leveldir object and array for new numeration */
  dirs = 0;
  for (i = 0; i < depth + 1; i++) {
    entrys = 0;
    for (j = 0; j < lold->directions[i]; j++) {	/* Count new directions */
      if (used[i][j] > 0) {
	entrys++;
      }
    }
    idx[i] = (uint *) allocmem(sizeof(uint) * lold->directions[i]);
    ldir->directions[i] = entrys;
    ldir->dir[i] = (preal *) allocmem(sizeof(preal) * entrys);
    dirs += entrys;
  }

  ldir->dirmem = dirmem = allocreal(dirs * ldir->dim);
  for (i = 0; i < depth + 1; i++) {
    for (j = 0; j < ldir->directions[i]; j++) {
      ldir->dir[i][j] = dirmem;
      dirmem += ldir->dim;
    }
  }

  /* Fill new leveldir object */
  copy_direction(used, idx, ldir, lold);

  /* Last part is changing directional cluster and block */
  /* Since we only have one directional cluster for both rows and columns we only need to call it once */
  change_dcluster(t, used, idx, tmp_dirson, ldir, 0);

  for (i = 0; i < depth; i++) {
    freemem(tmp_dirson[i]);
    freemem(used[i]);
  }
  freemem(used[depth]);
  free(tmp_dirson);
  free(used);

  /* Now the block */
  change_dblock(b, idx, 0);

  for (i = 0; i < depth + 1; i++) {
    freemem(idx[i]);
  }
  freemem(idx);
  del_leveldir(lold);

  return ldir;
}

/* ------------------------------------------------------------
 * Enumeration
 * ------------------------------------------------------------ */


static void
enumerate(pdblock b, uint bname, pdblock *bn)
{

  uint      bname1;
  uint      i, j;

  bn[bname] = b;

  bname1 = bname + 1;
  for (j = 0; j < b->csons; j++) {
    for (i = 0; i < b->rsons; i++) {

      enumerate(b->son[i + j * b->rsons], bname1, bn);
      bname1 += b->son[i + j * b->rsons]->desc;
    }
  }
  assert(bname1 == bname + b->desc);
}


pdblock  *
enumerate_dblock(pdblock b)
{

  pdblock  *bn;

  bn = (pdblock *) allocmem((size_t) sizeof(pdblock) * b->desc);

  enumerate(b, 0, bn);

  return bn;
}
