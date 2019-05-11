
#include "tri2dp1.h"

#include <stdio.h>
#include <assert.h>
#include <math.h>

#include "basic.h"

ptri2dp1
new_tri2dp1(pctri2d t2)
{
  ptri2dp1  dc;

  bool     *is_dof;
  uint     *idx2dof;
  uint      ndof, nfix;
  uint      i;

  dc = (ptri2dp1) allocmem(sizeof(tri2dp1));
  is_dof = dc->is_dof =
    (bool *) allocmem((size_t) sizeof(bool) * t2->vertices);
  idx2dof = dc->idx2dof = (uint *) allocuint(t2->vertices);
  dc->t2 = t2;

  for (i = 0; i < t2->vertices; i++)	/* Find d.o.f.s */
    is_dof[i] = (t2->xb[i] == 0);

  ndof = 0;
  nfix = 0;
  for (i = 0; i < t2->vertices; i++) {
    if (is_dof[i])
      idx2dof[i] = ndof++;
    else
      idx2dof[i] = nfix++;
  }

  dc->ndof = ndof;
  dc->nfix = nfix;

  return dc;
}

void
del_tri2dp1(ptri2dp1 dc)
{
  freemem(dc->is_dof);
  freemem(dc->idx2dof);
  freemem(dc);
}

psparsematrix
build_tri2dp1_sparsematrix(pctri2dp1 dc)
{
  pctri2d   t2 = dc->t2;
  const bool *is_dof = dc->is_dof;
  const uint *idx2dof = dc->idx2dof;
  uint      ndof = dc->ndof;
  uint      triangles = dc->t2->triangles;
  psparsepattern sp;
  psparsematrix A;
  uint      i, j, ii, jj, d, xt[3];

  sp = new_sparsepattern(ndof, ndof);

  for (d = 0; d < triangles; d++) {
    /* Get vertices */
    getvertices_tri2d(t2, d, xt);
    /* Build sparsepattern */
    for (i = 0; i < 3; i++) {
      if (is_dof[xt[i]]) {
	ii = idx2dof[xt[i]];
	for (j = 0; j < 3; j++) {
	  if (is_dof[xt[j]]) {
	    jj = idx2dof[xt[j]];
	    addnz_sparsepattern(sp, ii, jj);
	  }
	}
      }
    }
  }

  /* Build sparsematrix from sparsepattern */
  A = new_zero_sparsematrix(sp);

  del_sparsepattern(sp);

  return A;
}

psparsematrix
build_tri2dp1_interaction_sparsematrix(pctri2dp1 dc)
{
  pctri2d   t2 = dc->t2;
  const bool *is_dof = dc->is_dof;
  const uint *idx2dof = dc->idx2dof;
  uint      ndof = dc->ndof;
  uint      nfix = dc->nfix;
  psparsepattern sp;
  psparsematrix Af;
  uint      i, j, ii, jj, t, xt[3];

  sp = new_sparsepattern(ndof, nfix);

  for (t = 0; t < t2->triangles; t++) {
    /* Get vertices */
    getvertices_tri2d(t2, t, xt);
    /* Build sparsepattern */
    for (i = 0; i < 3; i++) {
      if (is_dof[xt[i]]) {
	ii = idx2dof[xt[i]];
	assert(ii < ndof);

	for (j = 0; j < 3; j++) {
	  if (!is_dof[xt[j]]) {
	    jj = idx2dof[xt[j]];
	    assert(jj < nfix);

	    addnz_sparsepattern(sp, ii, jj);
	  }
	}
      }
    }
  }

  /* Build sparsematrix from sparsepattern */
  Af = new_zero_sparsematrix(sp);

  del_sparsepattern(sp);

  return Af;
}

psparsematrix
build_tri2dp1_prolongation_sparsematrix(pctri2dp1 dfine, pctri2dp1 dcoarse,
					pctri2dref rf)
{
  pctri2d   gr = dfine->t2;
  const     uint(*e)[2] = (const uint(*)[2]) dcoarse->t2->e;
  const bool *fine_dof = dfine->is_dof;
  const uint *fine_idx2dof = dfine->idx2dof;
  const bool *coarse_dof = dcoarse->is_dof;
  const uint *coarse_idx2dof = dcoarse->idx2dof;
  uint      nfine = dfine->ndof;
  uint      ncoarse = dcoarse->ndof;
  psparsepattern sp;
  psparsematrix P;
  uint      i, j, k, ii, jj;

  sp = new_sparsepattern(nfine, ncoarse);

  for (i = 0; i < gr->vertices; i++)
    if (fine_dof[i]) {
      ii = fine_idx2dof[i];
      switch (rf->xt[i]) {
      case 0:
	j = rf->xf[i];
	if (coarse_dof[j]) {
	  jj = coarse_idx2dof[j];
	  addnz_sparsepattern(sp, ii, jj);
	}
	break;
      case 1:
	k = rf->xf[i];
	j = e[k][0];
	if (coarse_dof[j]) {
	  jj = coarse_idx2dof[j];
	  addnz_sparsepattern(sp, ii, jj);
	}
	j = e[k][1];
	if (coarse_dof[j]) {
	  jj = coarse_idx2dof[j];
	  addnz_sparsepattern(sp, ii, jj);
	}
	break;
      default:
	(void) fprintf(stderr,
		       "Unknown father type %u of vertex %u\n", rf->xf[i], i);
	del_sparsepattern(sp);
	return 0;
      }
    }

  P = new_zero_sparsematrix(sp);

  del_sparsepattern(sp);

  for (i = 0; i < gr->vertices; i++)
    if (fine_dof[i]) {
      ii = fine_idx2dof[i];
      switch (rf->xt[i]) {
      case 0:
	j = rf->xf[i];
	if (coarse_dof[j]) {
	  jj = coarse_idx2dof[j];
	  addentry_sparsematrix(P, ii, jj, 1.0);
	}
	break;
      case 1:
	k = rf->xf[i];
	j = e[k][0];
	if (coarse_dof[j]) {
	  jj = coarse_idx2dof[j];
	  addentry_sparsematrix(P, ii, jj, 0.5);
	}
	j = e[k][1];
	if (coarse_dof[j]) {
	  jj = coarse_idx2dof[j];
	  addentry_sparsematrix(P, ii, jj, 0.5);
	}
	break;
      default:
	(void) fprintf(stderr,
		       "Unknown father type %u of vertex %u\n", rf->xf[i], i);
	del_sparsematrix(P);
	return 0;
      }
    }

  return P;
}

void
assemble_tri2dp1_laplace_sparsematrix(pctri2dp1 dc, psparsematrix A,
				      psparsematrix Af)
{
  const     real(*x)[2] = (const real(*)[2]) dc->t2->x;
  const     uint(*e)[2] = (const uint(*)[2]) dc->t2->e;
  const     uint(*t)[3] = (const uint(*)[3]) dc->t2->t;
  const bool *is_dof = dc->is_dof;
  const uint *idx2dof = dc->idx2dof;
  uint      triangles = dc->t2->triangles;
  uint      ndof = dc->ndof;
  uint      nfix = dc->nfix;
  uint      i, j, ii, jj, d;
  real      det;
  real      gr[3][2];
  real      At[3][3];
  uint      xt[3];

  for (d = 0; d < triangles; d++) {
    /* Get vertices */
    xt[0] = e[t[d][0]][0];
    xt[1] = e[t[d][0]][1];
    if (e[t[d][1]][0] == e[t[d][0]][0] || e[t[d][1]][0] == e[t[d][0]][1])
      xt[2] = e[t[d][1]][1];
    else
      xt[2] = e[t[d][1]][0];

    assert(xt[0] != xt[1] && xt[0] != xt[2] && xt[1] != xt[2]);

    det = (x[xt[0]][0] - x[xt[1]][0]) * (x[xt[2]][1] - x[xt[1]][1])
      - (x[xt[0]][1] - x[xt[1]][1]) * (x[xt[2]][0] - x[xt[1]][0]);

    gr[0][0] = (x[xt[2]][1] - x[xt[1]][1]);
    gr[0][1] = (x[xt[1]][0] - x[xt[2]][0]);

    gr[1][0] = (x[xt[0]][1] - x[xt[2]][1]);
    gr[1][1] = (x[xt[2]][0] - x[xt[0]][0]);

    gr[2][0] = (x[xt[1]][1] - x[xt[0]][1]);
    gr[2][1] = (x[xt[0]][0] - x[xt[1]][0]);

    /* Compute element matrix */
    for (i = 0; i < 3; i++)
      for (j = 0; j < 3; j++)
	At[i][j] =
	  (gr[i][0] * gr[j][0] + gr[i][1] * gr[j][1]) / fabs(det) / 2.0;

    /* Add to system matrix */
    for (i = 0; i < 3; i++) {
      if (is_dof[xt[i]]) {
	ii = idx2dof[xt[i]];
	for (j = 0; j < 3; j++) {
	  if (is_dof[xt[j]]) {
	    jj = idx2dof[xt[j]];
	    addentry_sparsematrix(A, ii, jj, At[i][j]);
	  }
	}
      }
    }

    /* Add to interaction matrix */
    if (Af) {
      for (i = 0; i < 3; i++) {
	if (is_dof[xt[i]]) {
	  ii = idx2dof[xt[i]];
	  assert(ii < ndof);

	  for (j = 0; j < 3; j++) {
	    if (!is_dof[xt[j]]) {
	      jj = idx2dof[xt[j]];
	      assert(jj < nfix);
	      addentry_sparsematrix(Af, ii, jj, At[i][j]);
	    }
	  }
	}
      }
    }

  }
}

void
assemble_tri2dp1_mass_sparsematrix(pctri2dp1 dc, psparsematrix M,
				   psparsematrix Mf)
{
  const     real(*x)[2] = (const real(*)[2]) dc->t2->x;
  const     uint(*e)[2] = (const uint(*)[2]) dc->t2->e;
  const     uint(*t)[3] = (const uint(*)[3]) dc->t2->t;
  const bool *is_dof = dc->is_dof;
  const uint *idx2dof = dc->idx2dof;
  uint      triangles = dc->t2->triangles;
  uint      ndof = dc->ndof;
  uint      nfix = dc->nfix;
  uint      i, j, ii, jj, d;
  real      det;
  real      Mt[3][3];
  uint      xt[3];

  for (d = 0; d < triangles; d++) {
    /* Get vertices */
    xt[0] = e[t[d][0]][0];
    xt[1] = e[t[d][0]][1];
    if (e[t[d][1]][0] == e[t[d][0]][0] || e[t[d][1]][0] == e[t[d][0]][1])
      xt[2] = e[t[d][1]][1];
    else
      xt[2] = e[t[d][1]][0];

    assert(xt[0] != xt[1] && xt[0] != xt[2] && xt[1] != xt[2]);

    det = (x[xt[1]][0] - x[xt[0]][0]) * (x[xt[2]][1] - x[xt[0]][1])
      - (x[xt[1]][1] - x[xt[0]][1]) * (x[xt[2]][0] - x[xt[0]][0]);

    /* Compute element matrix */
    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
	if (i == j)
	  Mt[i][j] = fabs(det) / 12.0;
	else
	  Mt[i][j] = fabs(det) / 24.0;
      }
    }

    /* Add to system matrix */
    for (i = 0; i < 3; i++) {
      if (is_dof[xt[i]]) {
	ii = idx2dof[xt[i]];
	for (j = 0; j < 3; j++) {
	  if (is_dof[xt[j]]) {
	    jj = idx2dof[xt[j]];
	    addentry_sparsematrix(M, ii, jj, Mt[i][j]);
	  }
	}
      }
    }

    /* Add to interaction matrix */
    if (Mf) {
      for (i = 0; i < 3; i++) {
	if (is_dof[xt[i]]) {
	  ii = idx2dof[xt[i]];
	  assert(ii < ndof);

	  for (j = 0; j < 3; j++) {
	    if (!is_dof[xt[j]]) {
	      jj = idx2dof[xt[j]];
	      assert(jj < nfix);

	      addentry_sparsematrix(Mf, ii, jj, Mt[i][j]);
	    }
	  }
	}
      }
    }

  }
}

void
assemble_tri2dp1_dirichlet_avector(pctri2dp1 dc,
				   field(*d) (const real * x, void *fdata),
				   void *fdata, pavector dv)
{
  const     real(*x)[2] = (const real(*)[2]) dc->t2->x;
  const bool *is_dof = dc->is_dof;
  const uint *idx2dof = dc->idx2dof;
  uint      vertices = dc->t2->vertices;
  uint      nfix = dc->nfix;
  uint      i, ii;

  assert(dv->dim == nfix);

  for (i = 0; i < vertices; i++) {
    if (!is_dof[i]) {
      ii = idx2dof[i];
      assert(ii < nfix);

      dv->v[ii] = (d ? d(x[i], fdata) : 0.0);
    }
  }
}

void
assemble_tri2dp1_functional_avector(ptri2dp1 dc,
				    field(*f) (const real * x, void *fdata),
				    void *fdata, pavector fv)
{
  pctri2d   t2 = dc->t2;
  const     real(*x)[2] = (const real(*)[2]) t2->x;
  const     uint(*e)[2] = (const uint(*)[2]) t2->e;
  const     uint(*t)[3] = (const uint(*)[3]) t2->t;
  const bool *is_dof = dc->is_dof;
  const uint *idx2dof = dc->idx2dof;
  uint      triangles = dc->t2->triangles;
  uint      ndof = dc->ndof;
  real      xt[3][2], xd[2][2], xm[2], vol;
  field     f01, f02, f12;
  uint      i, ii, d, v[3];

  assert(fv->dim == ndof);

  clear_avector(fv);

  for (d = 0; d < triangles; d++) {
    /* Get vertices */
    v[0] = e[t[d][0]][0];
    v[1] = e[t[d][0]][1];
    if (e[t[d][1]][0] == e[t[d][0]][0] || e[t[d][1]][0] == e[t[d][0]][1])
      v[2] = e[t[d][1]][1];
    else
      v[2] = e[t[d][1]][0];

    /* Get vertex coordinates */
    for (i = 0; i < 3; i++) {
      xt[i][0] = x[v[i]][0];
      xt[i][1] = x[v[i]][1];
    }

    /* Compute difference vectors */
    for (i = 0; i < 2; i++) {
      xd[i][0] = xt[i + 1][0] - xt[0][0];
      xd[i][1] = xt[i + 1][1] - xt[0][1];
    }

    /* Compute determinant */
    vol = fabs(xd[0][0] * xd[1][1] - xd[0][1] * xd[1][0]);

    /* Evaluate function in midpoints */
    xm[0] = 0.5 * (xt[0][0] + xt[1][0]);
    xm[1] = 0.5 * (xt[0][1] + xt[1][1]);
    f01 = f(xm, fdata);
    xm[0] = 0.5 * (xt[0][0] + xt[2][0]);
    xm[1] = 0.5 * (xt[0][1] + xt[2][1]);
    f02 = f(xm, fdata);
    xm[0] = 0.5 * (xt[1][0] + xt[2][0]);
    xm[1] = 0.5 * (xt[1][1] + xt[2][1]);
    f12 = f(xm, fdata);

    /* Approximate integrals by midpoint rule */
    if (is_dof[v[0]]) {
      ii = idx2dof[v[0]];
      assert(ii < ndof);
      fv->v[ii] += 0.5 * (f01 + f02) * vol / 6.0;
    }
    if (is_dof[v[1]]) {
      ii = idx2dof[v[1]];
      assert(ii < ndof);
      fv->v[ii] += 0.5 * (f01 + f12) * vol / 6.0;
    }
    if (is_dof[v[2]]) {
      ii = idx2dof[v[2]];
      assert(ii < ndof);
      fv->v[ii] += 0.5 * (f02 + f12) * vol / 6.0;
    }
  }
}

real
normmax_tri2dp1(pctri2dp1 dc,
		field(*u) (const real * x, void *fdata), void *fdata,
		pavector xs)
{
  const     real(*x)[2] = (const real(*)[2]) dc->t2->x;
  const bool *is_dof = dc->is_dof;
  const uint *idx2dof = dc->idx2dof;
  uint      vertices = dc->t2->vertices;
  uint      ndof = dc->ndof;
  uint      i, ii;
  real      val, norm;

  assert(xs->dim == ndof);

  norm = 0.0;
  for (i = 0; i < vertices; i++)
    if (is_dof[i]) {
      ii = idx2dof[i];
      assert(ii < ndof);

      val = fabs(u(x[i], fdata) - xs->v[ii]);
      if (val > norm)
	norm = val;
    }

  return norm;
}

real
norml2_tri2dp1(pctri2dp1 dc,
	       field(*f) (const real * x, void *fdata), void *fdata,
	       pcavector b, pcavector bf)
{
  const     real(*x)[2] = (const real(*)[2]) dc->t2->x;
  const     uint(*e)[2] = (const uint(*)[2]) dc->t2->e;
  const     uint(*t)[3] = (const uint(*)[3]) dc->t2->t;
  const bool *is_dof = dc->is_dof;
  const uint *idx2dof = dc->idx2dof;
  uint      triangles = dc->t2->triangles;
  uint      ndof = dc->ndof;
  uint      nfix = dc->nfix;
  real      xt[3][2], xm[2], area;
  real      val, norm, tnorm;
  uint      vl[3];
  uint      i, ii, j, k;

  assert(b == 0 || b->dim == ndof);
  assert(bf == 0 || bf->dim == nfix);

  norm = 0.0;

  for (k = 0; k < triangles; k++) {
    /* Get vertices (to compute Jacobi determinant) */
    getvertices_tri2d(dc->t2, k, vl);

    /* Get vertex coordinates */
    for (i = 0; i < 3; i++) {
      xt[i][0] = x[vl[i]][0];
      xt[i][1] = x[vl[i]][1];
    }

    /* Compute area */
    area = fabs((xt[1][0] - xt[0][0]) * (xt[2][1] - xt[0][1]) -
		(xt[1][1] - xt[0][1]) * (xt[2][0] - xt[0][0])) * 0.5;

    /* Approximate integral by midpoint quadrature */
    tnorm = 0.0;
    for (j = 0; j < 3; j++) {
      val = 0.0;

      /* Evaluate in start point */
      i = e[t[k][j]][0];
      if (is_dof[i]) {
	if (b) {
	  ii = idx2dof[i];
	  assert(ii < ndof);

	  val += 0.5 * b->v[ii];
	}
      }
      else {
	if (bf) {
	  ii = idx2dof[i];
	  assert(ii < nfix);

	  val += 0.5 * bf->v[ii];
	}
      }

      /* Contribution to midpoint coordinates */
      xm[0] = 0.5 * x[i][0];
      xm[1] = 0.5 * x[i][1];

      /* Evaluate in end point */
      i = e[t[k][j]][1];
      if (is_dof[i]) {
	if (b) {
	  ii = idx2dof[i];
	  assert(ii < ndof);

	  val += 0.5 * b->v[ii];
	}
      }
      else {
	if (bf) {
	  ii = idx2dof[i];
	  assert(ii < nfix);

	  val += 0.5 * bf->v[ii];
	}
      }

      /* Find midpoint coordinates */
      xm[0] += 0.5 * x[i][0];
      xm[1] += 0.5 * x[i][1];

      /* Subtract reference function */
      if (f)
	val -= f(xm, fdata);

      tnorm += ABSSQR(val);
    }

    norm += tnorm * area / 3.0;
  }

  return REAL_SQRT(norm);
}


pclustergeometry
build_tri2dp1_clustergeometry(pctri2dp1 p1, uint * idx)
{

  pclustergeometry cg;
  uint      i, j, k, c;
  uint      ndof = p1->ndof;
  pctri2d   t2 = p1->t2;
  uint      vertices = t2->vertices;
  uint      triangles = t2->triangles;
  const bool *is_dof = p1->is_dof;
  const uint *idx2dof = p1->idx2dof;
  const     real(*x)[2] = (const real(*)[2]) t2->x;
  uint      v[3];
  real      min[2], max[2];

  cg = new_clustergeometry(2, ndof);

  /* Copy characteristic vertices into clustergeometry structure and
   * compute initial values of smin and smax */
  for (i = 0; i < vertices; i++) {
    if (is_dof[i]) {
      c = idx2dof[i];

      cg->x[c][0] = x[i][0];
      cg->x[c][1] = x[i][1];

      cg->smin[c][0] = x[i][0];
      cg->smin[c][1] = x[i][1];
      cg->smax[c][0] = x[i][0];
      cg->smax[c][1] = x[i][1];

      idx[c] = c;
    }
  }

  for (i = 0; i < triangles; i++) {	/* all triangles */
    getvertices_tri2d(t2, i, v);	/* get all vertices of triangle i */

    min[0] = REAL_MIN3(x[v[0]][0], x[v[1]][0], x[v[2]][0]);
    min[1] = REAL_MIN3(x[v[0]][1], x[v[1]][1], x[v[2]][1]);
    max[0] = REAL_MAX3(x[v[0]][0], x[v[1]][0], x[v[2]][0]);
    max[1] = REAL_MAX3(x[v[0]][1], x[v[1]][1], x[v[2]][1]);

    for (j = 0; j < 3; j++) {	/* all vertices */
      if (is_dof[v[j]]) {	/* vertex is degree of freedom */
	k = idx2dof[v[j]];

	cg->smin[k][0] = REAL_MIN(cg->smin[k][0], min[0]);
	cg->smin[k][1] = REAL_MIN(cg->smin[k][1], min[1]);
	cg->smax[k][0] = REAL_MAX(cg->smax[k][0], max[0]);
	cg->smax[k][1] = REAL_MAX(cg->smax[k][1], max[1]);
      }
    }
  }

  update_point_bbox_clustergeometry(cg, ndof, idx);

  return cg;
}
