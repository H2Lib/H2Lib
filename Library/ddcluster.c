
#include "basic.h"
#include "ddcluster.h"

#include <stdio.h>

pcluster
build_regular_interface_cluster(pclustergeometry cg, uint size, uint * idx,
				uint clf, uint dim, uint direction,
				uint levelint)
{

  pcluster  c;			/*new cluster */
  uint      size0, size1;	/*sizes of sons */
  uint      i, j, newd;		/*Indices, new direction */
  real      a, b, m;		/*left, right and midpoint of the bounding box in direction */

  size0 = 0;
  size1 = 0;

  /*Select new direction */
  if (direction < dim - 1) {
    newd = direction + 1;
  }
  else {
    newd = 0;
  }

  if (size > clf) {

    if (levelint % dim) {
      /*build two sons */

      levelint++;
      m = cg->hmax[direction] - cg->hmin[direction];

      if (m > 0.0) {
	/*Split idx */
	m = (cg->hmax[direction] + cg->hmin[direction]) / 2.0;
	for (i = 0; i < size; i++) {
	  if (cg->x[idx[i]][direction] < m) {
	    j = idx[i];
	    idx[i] = idx[size0];
	    idx[size0] = j;
	    size0++;
	  }
	  else {
	    size1++;
	  }
	}

	/* build sons */
	if (size0 > 0) {
	  if (size1 > 0) {
	    /* both sons are not empty */
	    assert(size0 > 0);
	    assert(size1 > 0);
	    c = new_cluster(size, idx, 2, cg->dim);

	    a = cg->hmin[direction];
	    b = cg->hmax[direction];
	    cg->hmax[direction] = m;
	    c->son[0] =
	      build_regular_interface_cluster(cg, size0, idx, clf, dim, newd,
					      levelint);

	    cg->hmax[direction] = b;
	    cg->hmin[direction] = m;
	    c->son[1] =
	      build_regular_interface_cluster(cg, size1, idx + size0, clf,
					      dim, newd, levelint);

	    cg->hmin[direction] = a;
	    update_bbox_cluster(c);
	  }
	  else {
	    /* the first son is not empty , the second son is empty */
	    assert(size0 > 0);
	    assert(size1 == 0);
	    c = new_cluster(size, idx, 1, cg->dim);

	    b = cg->hmax[direction];
	    cg->hmax[direction] = m;
	    c->son[0] =
	      build_regular_interface_cluster(cg, size, idx, clf, dim, newd,
					      levelint);

	    cg->hmax[direction] = b;
	    update_bbox_cluster(c);
	  }
	}
	else {
	  /* the first son is empty, the second son is not empty */
	  assert(size1 > 0);
	  assert(size0 == 0);
	  c = new_cluster(size, idx, 1, cg->dim);

	  a = cg->hmin[direction];
	  cg->hmin[direction] = m;
	  c->son[0] =
	    build_regular_interface_cluster(cg, size, idx, clf, dim, newd,
					    levelint);

	  cg->hmin[direction] = a;
	  update_bbox_cluster(c);
	}
      }
      else {
	assert(m == 0);
	c = new_cluster(size, idx, 1, cg->dim);

	c->son[0] =
	  build_regular_interface_cluster(cg, size, idx, clf, dim, newd,
					  levelint);

	update_bbox_cluster(c);
      }
    }
    else {			/*build one son */
      levelint++;
      c = new_cluster(size, idx, 1, cg->dim);

      c->son[0] =
	build_regular_interface_cluster(cg, size, idx, clf, dim, newd,
					levelint);
      update_bbox_cluster(c);
    }

  }
  else {
    c = new_cluster(size, idx, 0, cg->dim);
    update_support_bbox_cluster(cg, c);
  }
  c->type = 2;
  update_cluster(c);
  return c;
}

pcluster
build_regular_dd_cluster(pclustergeometry cg, uint size, uint * idx, uint clf,
			 psparsematrix sp, uint dim, uint direction,
			 uint * flag)
{

  pcluster  c;
  uint      i, j, newd, tmp, tmp2;
  uint      size0, size1, size2;	/*sizes of sons */
  real      a, b, m;		/*left, right and midpoint of the bounding box in direction direction */
  real      h, diam_max, diam;

  if (size > clf) {

    size0 = 0;
    size1 = 0;
    size2 = 0;

    /*Find new direction */
    if (direction < dim - 1) {
      newd = direction + 1;
    }
    else {
      newd = 0;
    }

    m = cg->hmax[direction] - cg->hmin[direction];

    if (m > 0.0) {

      /*Split idx in two pieces */
      m = (cg->hmax[direction] + cg->hmin[direction]) / 2.0;
      for (i = 0; i < size; i++) {
	if (cg->x[idx[i]][direction] < m) {
	  j = idx[i];
	  idx[i] = idx[size0];
	  idx[size0] = j;
	  size0++;
	}
	else {
	  size1++;
	}
      }

      for (i = 0; i < size0; i++)
	flag[idx[i]] = 1;
      tmp2 = 0;
      i = size0;
      while (i < size0 + size1) {
	tmp2 = 0;
	for (j = sp->row[idx[i]]; j < sp->row[idx[i] + 1]; j++) {
	  if (flag[sp->col[j]] == 1) {
	    tmp = idx[i];
	    idx[i] = idx[size - size2 - 1];
	    idx[size - size2 - 1] = tmp;
	    size2++;
	    size1--;
	  }
	  if (flag[sp->col[j]] == 1) {
	    tmp2 = 1;
	    break;
	  }
	}
	if (tmp2 == 0) {
	  i++;
	}
      }

      for (i = 0; i < size0; i++)
	flag[idx[i]] = 0;

      /*Build sons */
      if (size0 > 0) {

	if (size1 > 0) {

	  if (size2 > 0) {
	    /* all three sons are not empty */
	    assert(size0 > 0);
	    assert(size1 > 0);
	    assert(size2 > 0);

	    c = new_cluster(size, idx, 3, cg->dim);

	    a = cg->hmin[direction];
	    b = cg->hmax[direction];
	    cg->hmax[direction] = m;
	    c->son[0] =
	      build_regular_dd_cluster(cg, size0, idx, clf, sp, dim, newd,
				       flag);

	    cg->hmax[direction] = b;
	    cg->hmin[direction] = m;
	    c->son[1] =
	      build_regular_dd_cluster(cg, size1, idx + size0, clf, sp, dim,
				       newd, flag);

	    diam_max = 0.0;
	    for (j = size1; j < size; j++) {
	      for (i = 0; i < c->dim; i++) {
		diam = cg->smax[idx[j]][i] - cg->smin[idx[j]][i];
		if (diam > diam_max) {
		  diam_max = diam;
		}
	      }
	    }


	    h = diam_max;
	    cg->hmin[direction] = (a + b) / 2.0 - h;
	    cg->hmax[direction] = (a + b) / 2.0 + h;
	    c->son[2] =
	      build_regular_interface_cluster(cg, size2, idx + size0 + size1,
					      clf, dim, newd, 1);

	    cg->hmin[direction] = a;
	    cg->hmax[direction] = b;
	    update_bbox_cluster(c);
	  }
	  else {
	    /* first and second son are not empty, the third son is empty */
	    assert(size2 == 0);
	    assert(size0 > 0);
	    assert(size1 > 0);

	    c = new_cluster(size, idx, 2, cg->dim);

	    a = cg->hmin[direction];
	    b = cg->hmax[direction];
	    cg->hmax[direction] = m;
	    c->son[0] =
	      build_regular_dd_cluster(cg, size0, idx, clf, sp, dim, newd,
				       flag);

	    cg->hmax[direction] = b;
	    cg->hmin[direction] = m;
	    c->son[1] =
	      build_regular_dd_cluster(cg, size1, idx + size0, clf, sp, dim,
				       newd, flag);

	    cg->hmin[direction] = a;
	    cg->hmax[direction] = b;
	    update_bbox_cluster(c);
	  }
	}
	else {
	  if (size2 > 0) {
	    /*Dieser Fall tritt nie ein!!!??? */
	    /* the first ist not empty, the second son is empty, the third is not */
	    assert(size0 > 0);
	    assert(size1 == 0);
	    assert(size2 > 0);

	    assert(0);
	  }
	  else {
	    /*first son is not empty,second and third son are empty */
	    assert(size2 == 0);
	    assert(size1 == 0);
	    assert(size0 > 0);
	    c = new_cluster(size, idx, 1, cg->dim);

	    b = cg->hmax[direction];
	    cg->hmax[direction] = m;
	    c->son[0] =
	      build_regular_dd_cluster(cg, size0, idx, clf, sp, dim, newd,
				       flag);

	    cg->hmax[direction] = b;
	    update_bbox_cluster(c);
	  }
	}
      }

      else {
	if (size1 > 0) {
	  /*the first in empty, the second son is not empty, the third is empty */
	  assert(size0 == 0);
	  assert(size1 > 0), assert(size2 == 0);

	  c = new_cluster(size, idx, 1, cg->dim);
	  a = cg->hmin[direction];
	  cg->hmin[direction] = m;

	  c->son[0] =
	    build_regular_dd_cluster(cg, size1, idx + size0, clf, sp, dim,
				     newd, flag);
	  cg->hmin[direction] = a;
	  update_bbox_cluster(c);
	}
	else {
	  assert(0);
	}
      }


    }
    else {
      c = new_cluster(size, idx, 1, cg->dim);

      c->son[0] =
	build_regular_dd_cluster(cg, size, idx, clf, sp, dim, newd, flag);

      update_bbox_cluster(c);
    }
  }
  else {
    c = new_cluster(size, idx, 0, cg->dim);
    update_support_bbox_cluster(cg, c);
  }
  c->type = 1;
  update_cluster(c);

  return c;
}


bool
admissible_dd_cluster(pcluster s, pcluster t, void *data)
{

  bool      b, a;
// real eta;

// eta = *(real *) data;

  a = admissible_2_min_cluster(s, t, data);

  if (a == true)
    b = true;
  else {
    a = (s == t);
    if ((s->type == 1 && s->type == t->type && a == false /*&& s != t */ )) {
      b = true;
    }

    else
      b = false;
  }
  return b;
}

pcluster
build_adaptive_interface_cluster(pclustergeometry cg, uint size, uint * idx,
				 uint clf, uint dim, uint levelint)
{
  pcluster  c;

  uint      size0, size1;
  uint      i, j, direction;
  real      a, m;

  size0 = 0;
  size1 = 0;

  if (size > clf) {
    if (levelint % dim) {
      levelint++;

      update_point_bbox_clustergeometry(cg, size, idx);

      /* compute the direction of partition */
      direction = 0;
      a = cg->hmax[0] - cg->hmin[0];
      for (j = 1; j < cg->dim; j++) {
	m = cg->hmax[j] - cg->hmin[j];
	if (a < m) {
	  a = m;
	  direction = j;
	}
      }

      /* build sons */
      if (a > 0.0) {
	m = (cg->hmax[direction] + cg->hmin[direction]) / 2.0;
	size0 = 0;
	size1 = 0;

	for (i = 0; i < size; i++) {
	  if (cg->x[idx[i]][direction] < m) {
	    j = idx[i];
	    idx[i] = idx[size0];
	    idx[size0] = j;
	    size0++;
	  }
	  else {
	    size1++;
	  }
	}
	if (size0 > 0) {
	  if (size1 > 0) {
	    c = new_cluster(size, idx, 2, cg->dim);

	    c->son[0] =
	      build_adaptive_interface_cluster(cg, size0, idx, clf, dim,
					       levelint);
	    c->son[1] =
	      build_adaptive_interface_cluster(cg, size1, idx + size0, clf,
					       dim, levelint);
	    update_bbox_cluster(c);
	  }
	  else {
	    assert(size0 > 0);
	    assert(size1 == 0);
	    c = new_cluster(size, idx, 1, cg->dim);
	    c->son[0] =
	      build_adaptive_interface_cluster(cg, size, idx, clf, dim,
					       levelint);
	    update_bbox_cluster(c);
	  }
	}
	else {
	  assert(size0 == 0);
	  assert(size1 > 0);
	  c = new_cluster(size, idx, 1, cg->dim);
	  c->son[0] =
	    build_adaptive_interface_cluster(cg, size, idx, clf, dim,
					     levelint);
	  update_bbox_cluster(c);
	}

      }
      else {
	assert(a == 0.0);
	c = new_cluster(size, idx, 0, cg->dim);
	update_support_bbox_cluster(cg, c);
      }
    }

    else {
      levelint++;
      c = new_cluster(size, idx, 1, cg->dim);
      c->son[0] =
	build_adaptive_interface_cluster(cg, size, idx, clf, dim, levelint);
      update_bbox_cluster(c);
    }


  }

  else {
    /* size <= clf */
    c = new_cluster(size, idx, 0, cg->dim);
    update_support_bbox_cluster(cg, c);
  }

  c->type = 2;
  update_cluster(c);

  return c;
}



pcluster
build_adaptive_dd_cluster(pclustergeometry cg, uint size, uint * idx,
			  uint clf, psparsematrix sp, uint dim, uint * flag)
{
  pcluster  c;

  uint      i, j, direction, tmp;
  bool      inter;
  uint      size0, size1, size2;
  real      a, m;

  if (size > clf) {
    size0 = 0;
    size1 = 0;
    size2 = 0;

    update_point_bbox_clustergeometry(cg, size, idx);

    /* compute the direction of partition */
    direction = 0;
    a = cg->hmax[0] - cg->hmin[0];

    for (j = 1; j < cg->dim; j++) {
      m = cg->hmax[j] - cg->hmin[j];
      if (a < m) {
	a = m;
	direction = j;
      }
    }

    /* build sons */
    if (a > 0.0) {

      m = (cg->hmax[direction] + cg->hmin[direction]) / 2.0;

      for (i = 0; i < size; i++) {

	if (cg->x[idx[i]][direction] < m) {
	  j = idx[i];
	  idx[i] = idx[size0];
	  idx[size0] = j;
	  size0++;
	}
	else
	  size1++;
      }

      for (i = 0; i < size0; i++)
	flag[idx[i]] = 1;


      i = size0;
      while (i < size0 + size1) {
	inter = false;
	for (j = sp->row[idx[i]]; j < sp->row[idx[i] + 1]; j++) {

	  if (flag[sp->col[j]] == 1) {
	    tmp = idx[i];
	    idx[i] = idx[size - size2 - 1];
	    idx[size - size2 - 1] = tmp;
	    size2++;
	    size1--;
	    inter = true;
	    break;
	  }
	}
	if (inter == false)
	  i++;
      }

      for (i = 0; i < size0; i++)
	flag[idx[i]] = 0;

      /*Build sons */

      if (size0 > 0) {
	if (size1 > 0) {
	  if (size2 > 0) {
	    /* all three sons are not empty */
	    c = new_cluster(size, idx, 3, cg->dim);

	    c->son[0] =
	      build_adaptive_dd_cluster(cg, size0, idx, clf, sp, dim, flag);
	    c->son[1] =
	      build_adaptive_dd_cluster(cg, size1, idx + size0, clf, sp, dim,
					flag);
	    c->son[2] =
	      build_adaptive_interface_cluster(cg, size2, idx + size0 + size1,
					       clf, dim, 1);

	    update_bbox_cluster(c);
	  }
	  else {
	    assert(size2 == 0);
	    /* first and second son are not empty, the third son is empty */
	    c = new_cluster(size, idx, 2, cg->dim);

	    c->son[0] =
	      build_adaptive_dd_cluster(cg, size0, idx, clf, sp, dim, flag);
	    c->son[1] =
	      build_adaptive_dd_cluster(cg, size1, idx + size0, clf, sp, dim,
					flag);

	    update_bbox_cluster(c);
	  }
	}
	else {
	  /* size1 == 0 */
	  assert(0);		/*case not possible */
	  assert(size1 == 0);
	  assert(size2 == 0);
	  assert(size0 > 0);
	  /* the first ist not empty, the second son is empty, the third is not */
	  c = new_cluster(size, idx, 1, cg->dim);

	  c->son[0] =
	    build_adaptive_dd_cluster(cg, size0, idx, clf, sp, dim, flag);

	  update_bbox_cluster(c);
	}

      }
      else {
	assert(0);		/*case not possible */
	assert(size0 == 0);
	if (size1 > 0) {
	  c = new_cluster(size, idx, 1, cg->dim);

	  c->son[0] =
	    build_adaptive_dd_cluster(cg, size0, idx, clf, sp, dim, flag);

	  update_bbox_cluster(c);
	}
	else {
	  assert(0);
	}
      }

    }

    else {
      /* a == 0 */
      c = new_cluster(size, idx, 0, cg->dim);
      update_support_bbox_cluster(cg, c);
    }
  }

  else {
    /* size <= clf */
    c = new_cluster(size, idx, 0, cg->dim);
    update_support_bbox_cluster(cg, c);
  }
  c->type = 1;
  update_cluster(c);

  return c;
}

/**************************************************************
 * Clusterstrategy based on Clusterstrategy for the triangles
 **************************************************************/

pcluster
build_adaptive_interface_fromcluster_cluster(pcluster c1, pcluster c2,
					     pclustergeometry cg, uint size,
					     uint * idx, uint clf, uint dim,
					     uint levelint)
{
  pcluster  c;

  uint      size0, size1;
  uint      i, j, direction;
  real      a, m;

  size0 = 0;
  size1 = 0;

  if (size > clf) {
    if (levelint % dim) {
      levelint++;

      update_point_bbox_clustergeometry(cg, size, idx);

      /* compute the direction of partition */
      direction = 0;
      a = cg->hmax[0] - cg->hmin[0];
      for (j = 1; j < cg->dim; j++) {
	m = cg->hmax[j] - cg->hmin[j];
	if (a < m) {
	  a = m;
	  direction = j;
	}
      }

      /* build sons */
      if (a > 0.0) {
	m = (cg->hmax[direction] + cg->hmin[direction]) / 2.0;
	size0 = 0;
	size1 = 0;

	for (i = 0; i < size; i++) {
	  if (cg->x[idx[i]][direction] < m) {
	    j = idx[i];
	    idx[i] = idx[size0];
	    idx[size0] = j;
	    size0++;
	  }
	  else {
	    size1++;
	  }
	}
	c = new_cluster(size, idx, 2, cg->dim);

	c->son[0] =
	  build_adaptive_interface_fromcluster_cluster(c1, c2, cg, size0, idx,
						       clf, dim, levelint);
	c->son[1] =
	  build_adaptive_interface_fromcluster_cluster(c1, c2, cg, size1,
						       idx + size0, clf, dim,
						       levelint);

	update_bbox_cluster(c);
      }
      else {
	assert(a == 0.0);
	c = new_cluster(size, idx, 0, cg->dim);
	update_support_bbox_cluster(cg, c);
      }
    }

    else {
      levelint++;
      c = new_cluster(size, idx, 1, cg->dim);
      c->son[0] =
	build_adaptive_interface_fromcluster_cluster(c1, c2, cg, size, idx,
						     clf, dim, levelint);
      update_bbox_cluster(c);
    }


  }

  else {
    /* size <= clf */
    c = new_cluster(size, idx, 0, cg->dim);
    update_support_bbox_cluster(cg, c);
  }

  c->type = 2;
  update_cluster(c);

  return c;
}

pcluster
build_adaptive_fromcluster_cluster(pcluster ctri, pclustergeometry cg,
				   uint size, uint * idx, uint clf, uint dim)
{
  pcluster  c;
  pcluster  son1, son2;
  uint      size1, size2, size3;
  uint     *flag, *flag2;
  uint      i, j, count, k;

  if (ctri->sons > 0) {		/*Cluster for triangles/tetrahedra  has sons */
    assert(ctri->sons == 2);
    update_point_bbox_clustergeometry(cg, size, idx);
    /*Build sons */
    son1 = ctri->son[0];
    son2 = ctri->son[1];

    flag = allocuint(size);
    for (i = 0; i < size; i++)
      flag[i] = 0;

    /*Mark edges/faces of son1 */
    for (i = 0; i < size; i++) {
      count = 0;
      for (j = 0; j < ctri->dim; j++) {
	if ((son1->bmin[j] <= cg->x[idx[i]][j])
	    && (cg->x[idx[i]][j] <= son1->bmax[j]))
	  count++;
      }
      if (count == ctri->dim)
	flag[i] = 1;		/*Belongs to son1 */
    }

    /* Mark edges of son 2 und edges of both sons */
    for (i = 0; i < size; i++) {
      count = 0;
      for (j = 0; j < ctri->dim; j++) {
	if ((son2->bmin[j] <= cg->x[idx[i]][j])
	    && (cg->x[idx[i]][j] <= son2->bmax[j]))
	  count++;
      }
      if (count == ctri->dim) {
	if (flag[i] == 1)
	  flag[i] = 3;		/*edges of both sons */
	else if (flag[i] == 0)
	  flag[i] = 2;		/*edge of son2 */
	else
	  printf("error Kante übrig\t");
      }
    }
    flag2 = allocuint(size);	/*flag2 is copy of flag */
    for (i = 0; i < size; i++)
      flag2[i] = flag[i];

    /*Sort indices */
    size1 = 0;
    size2 = 0;
    size3 = 0;
    for (i = 0; i < size; i++) {
      if (flag[i] == 1) {	/*son1 */
	j = idx[i];
	idx[i] = idx[size1];
	idx[size1] = j;
	k = flag2[i];
	flag2[i] = flag2[size1];
	flag2[size1] = k;
	size1++;
      }
    }

    for (i = size1; i < size; i++) {
      if (flag2[i] == 2) {	/*son2 */
	j = idx[i];
	idx[i] = idx[size1 + size2];
	idx[size1 + size2] = j;
	k = flag2[i];
	flag2[i] = flag2[size1 + size2];
	flag2[size1 + size2] = k;
	size2++;
      }
    }

    for (i = size1 + size2; i < size; i++) {	/*Interface */
      if (flag2[i] == 3)
	size3++;
      else
	printf("error size\t");
    }
    //printf("cg->nidx = %u size1 = %u size2 = %u size3 = %u size = %u\n", cg->nidx, size1, size2, size3, size);
    if (size1 > 0) {
      if (size2 > 0) {
	if (size3 > 0) {
	  c = new_cluster(size, idx, 3, ctri->dim);
	  c->son[0] =
	    build_adaptive_fromcluster_cluster(son1, cg, size1, idx, clf,
					       dim);
	  c->son[1] =
	    build_adaptive_fromcluster_cluster(son2, cg, size2, idx + size1,
					       clf, dim);
	  /*adaptive clustering for interface clusters */
	  c->son[2] =
	    build_adaptive_interface_fromcluster_cluster(son1, son2, cg,
							 size3,
							 idx + size1 + size2,
							 clf, dim, 1);
	}
	else {
	  assert(size1 > 0);
	  assert(size2 > 0);
	  assert(size3 == 0);
	  c = new_cluster(size, idx, 2, ctri->dim);
	  c->son[0] =
	    build_adaptive_fromcluster_cluster(son1, cg, size1, idx, clf,
					       dim);
	  c->son[1] =
	    build_adaptive_fromcluster_cluster(son2, cg, size2, idx + size1,
					       clf, dim);
	}
      }
      else {
	assert(size1 > 0);
	assert(size2 == 0);
	if (size3 > 0) {
	  c = new_cluster(size, idx, 2, ctri->dim);
	  c->son[0] =
	    build_adaptive_fromcluster_cluster(son1, cg, size1, idx, clf,
					       dim);
	  /*adaptive clustering for interface clusters */
	  c->son[1] =
	    build_adaptive_interface_fromcluster_cluster(son1, son2, cg,
							 size3,
							 idx + size1 + size2,
							 clf, dim, 1);
	}
	else {
	  assert(size3 == 0);
	  c = new_cluster(size, idx, 1, ctri->dim);
	  c->son[0] =
	    build_adaptive_fromcluster_cluster(son1, cg, size1, idx, clf,
					       dim);
	}
      }
    }
    else {
      assert(size1 == 0);
      if (size2 > 0) {
	if (size3 > 0) {
	  c = new_cluster(size, idx, 2, ctri->dim);
	  c->son[0] =
	    build_adaptive_fromcluster_cluster(son1, cg, size2, idx, clf,
					       dim);
	  /*adaptive clustering for interface clusters */
	  c->son[1] =
	    build_adaptive_interface_fromcluster_cluster(son1, son2, cg,
							 size3,
							 idx + size1 + size2,
							 clf, dim, 1);
	}
	else {
	  assert(size2 > 0);
	  assert(size3 == 0);
	  c = new_cluster(size, idx, 1, ctri->dim);
	  c->son[0] =
	    build_adaptive_fromcluster_cluster(son1, cg, size2, idx, clf,
					       dim);
	}
      }
      else {
	assert(size1 == 0);
	assert(size2 == 0);
	if (size3 > 0) {
	  c = new_cluster(size, idx, 1, ctri->dim);
	  /*adaptive clustering for interface clusters */
	  c->son[0] =
	    build_adaptive_interface_fromcluster_cluster(son1, son2, cg,
							 size3,
							 idx + size1 + size2,
							 clf, dim, 1);
	}
	else {
	  printf("error\n");
	  c = new_cluster(size, idx, 0, ctri->dim);
	}

      }

    }
    update_bbox_cluster(c);
    freemem(flag);
    freemem(flag2);
  }
  else {
    //printf("size = %u\n", size);
    /*Cluster has no sons */
    c = new_cluster(size, idx, 0, ctri->dim);
    update_support_bbox_cluster(cg, c);
  }
  c->type = 1;
  update_cluster(c);

  return c;
}
