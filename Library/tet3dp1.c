
/* ------------------------------------------------------------
 * This is the file "tet3dp1.c" of the H2Lib package.
 * All rights reserved, Steffen Boerm 2015
 * ------------------------------------------------------------ */

#include "tet3dp1.h"

#include <stdio.h>

#include "basic.h"
#include "sparsepattern.h"

ptet3dp1
new_tet3dp1(pctet3d gr)
{
  ptet3dp1  dc;
  bool     *is_dof;
  uint     *idx2dof;
  uint      ndof, nfix;
  uint      i;

  dc = (ptet3dp1) allocmem(sizeof(tet3dp1));
  dc->is_dof = is_dof =
    (bool *) allocmem((size_t) sizeof(bool) * gr->vertices);
  dc->idx2dof = idx2dof = (uint *) allocuint(gr->vertices);
  dc->gr = gr;

  for (i = 0; i < gr->vertices; i++)	/* Find d.o.f.s */
    is_dof[i] = (gr->xb[i] == 0);

  ndof = 0;
  nfix = 0;
  for (i = 0; i < gr->vertices; i++)
    if (is_dof[i])
      idx2dof[i] = ndof++;
    else
      idx2dof[i] = nfix++;

  dc->ndof = ndof;
  dc->nfix = nfix;

  return dc;
}

void
del_tet3dp1(ptet3dp1 dc)
{
  freemem(dc->idx2dof);
  freemem(dc->is_dof);
  freemem(dc);
}

psparsematrix
build_tet3dp1_sparsematrix(pctet3dp1 dc)
{
  pctet3d   gr = dc->gr;
  const bool *is_dof = dc->is_dof;
  const uint *idx2dof = dc->idx2dof;
  uint      ndof = dc->ndof;
  psparsepattern sp;
  psparsematrix A;
  uint      i, j, ii, jj, t, v[4];

  sp = new_sparsepattern(ndof, ndof);

  for (t = 0; t < gr->tetrahedra; t++) {
    getvertices_tet3d(gr, t, v);
    for (i = 0; i < 4; i++)
      if (is_dof[v[i]]) {
	ii = idx2dof[v[i]];

	for (j = 0; j < 4; j++)
	  if (is_dof[v[j]]) {
	    jj = idx2dof[v[j]];

	    addnz_sparsepattern(sp, ii, jj);
	  }
      }
  }

  A = new_zero_sparsematrix(sp);

  del_sparsepattern(sp);

  return A;
}

psparsematrix
build_tet3dp1_interaction_sparsematrix(pctet3dp1 dc)
{
  pctet3d   gr = dc->gr;
  const bool *is_dof = dc->is_dof;
  const uint *idx2dof = dc->idx2dof;
  uint      ndof = dc->ndof;
  uint      nfix = dc->nfix;
  psparsepattern sp;
  psparsematrix Af;
  uint      i, j, ii, jj, t, v[4];

  sp = new_sparsepattern(ndof, nfix);

  for (t = 0; t < gr->tetrahedra; t++) {
    getvertices_tet3d(gr, t, v);
    for (i = 0; i < 4; i++)
      if (is_dof[v[i]]) {
	ii = idx2dof[v[i]];
	assert(ii < ndof);

	for (j = 0; j < 4; j++)
	  if (!is_dof[v[j]]) {
	    jj = idx2dof[v[j]];
	    assert(jj < nfix);

	    addnz_sparsepattern(sp, ii, jj);
	  }
      }
  }

  Af = new_zero_sparsematrix(sp);

  del_sparsepattern(sp);

  return Af;
}

psparsematrix
build_tet3dp1_prolongation_sparsematrix(pctet3dp1 dfine, pctet3dp1 dcoarse,
					pctet3dref rf)
{
  pctet3d   gr = dfine->gr;
  const     uint(*e)[2] = (const uint(*)[2]) dcoarse->gr->e;
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
assemble_tet3dp1_laplace_sparsematrix(pctet3dp1 dc, psparsematrix A,
				      psparsematrix Af)
{
  pctet3d   gr = dc->gr;
  const     real(*x)[3] = (const real(*)[3]) gr->x;
  uint      tetrahedra = gr->tetrahedra;
  const bool *is_dof = dc->is_dof;
  const uint *idx2dof = dc->idx2dof;
  uint      ndof = dc->ndof;
  uint      nfix = dc->nfix;
  real      xt[4][3];
  real      g[4][3];
  real      At[4][4];
  real      det;
  uint      v[4];
  uint      t, i, j, ii, jj;

  for (t = 0; t < tetrahedra; t++) {
    /* Get vertices */
    getvertices_tet3d(gr, t, v);

    /* Get vertex coordinates */
    for (i = 0; i < 4; i++) {
      xt[i][0] = x[v[i]][0];
      xt[i][1] = x[v[i]][1];
      xt[i][2] = x[v[i]][2];
    }

    /* Compute gradients and Jacobi determinant */
    for (i = 0; i < 4; i++) {
      g[i][0] = ((xt[(i + 2) % 4][1] - xt[(i + 1) % 4][1])
		 * (xt[(i + 3) % 4][2] - xt[(i + 1) % 4][2])
		 - (xt[(i + 2) % 4][2] - xt[(i + 1) % 4][2])
		 * (xt[(i + 3) % 4][1] - xt[(i + 1) % 4][1]));
      g[i][1] = ((xt[(i + 2) % 4][2] - xt[(i + 1) % 4][2])
		 * (xt[(i + 3) % 4][0] - xt[(i + 1) % 4][0])
		 - (xt[(i + 2) % 4][0] - xt[(i + 1) % 4][0])
		 * (xt[(i + 3) % 4][2] - xt[(i + 1) % 4][2]));
      g[i][2] = ((xt[(i + 2) % 4][0] - xt[(i + 1) % 4][0])
		 * (xt[(i + 3) % 4][1] - xt[(i + 1) % 4][1])
		 - (xt[(i + 2) % 4][1] - xt[(i + 1) % 4][1])
		 * (xt[(i + 3) % 4][0] - xt[(i + 1) % 4][0]));

      det = ((xt[i][0] - xt[(i + 1) % 4][0]) * g[i][0]
	     + (xt[i][1] - xt[(i + 1) % 4][1]) * g[i][1]
	     + (xt[i][2] - xt[(i + 1) % 4][2]) * g[i][2]);

      g[i][0] /= det;
      g[i][1] /= det;
      g[i][2] /= det;
    }

    /* Compute element matrix */
    for (i = 0; i < 4; i++)
      for (j = 0; j < 4; j++)
	At[i][j] = (g[i][0] * g[j][0]
		    + g[i][1] * g[j][1]
		    + g[i][2] * g[j][2]) * fabs(det) / 6.0;

    /* Add to system matrix */
    if (A)
      for (i = 0; i < 4; i++)
	if (is_dof[v[i]]) {
	  ii = idx2dof[v[i]];

	  for (j = 0; j < 4; j++)
	    if (is_dof[v[j]]) {
	      jj = idx2dof[v[j]];
	      addentry_sparsematrix(A, ii, jj, At[i][j]);
	    }
	}

    /* Add to interaction matrix */
    if (Af)
      for (i = 0; i < 4; i++)
	if (is_dof[v[i]]) {
	  ii = idx2dof[v[i]];
	  assert(ii < ndof);

	  for (j = 0; j < 4; j++)
	    if (!is_dof[v[j]]) {
	      jj = idx2dof[v[j]];
	      assert(jj < nfix);

	      addentry_sparsematrix(Af, ii, jj, At[i][j]);
	    }
	}
  }
}

void
assemble_tet3dp1_mass_sparsematrix(pctet3dp1 dc, psparsematrix M,
				   psparsematrix Mf)
{
  const real Mt[4][4] = { {1.0 / 60.0, 1.0 / 120.0, 1.0 / 120.0, 1.0 / 120.0},
  {1.0 / 120.0, 1.0 / 60.0, 1.0 / 120.0, 1.0 / 120.0},
  {1.0 / 120.0, 1.0 / 120.0, 1.0 / 60.0, 1.0 / 120.0},
  {1.0 / 120.0, 1.0 / 120.0, 1.0 / 120.0, 1.0 / 60.0}
  };
  pctet3d   gr = dc->gr;
  const     real(*x)[3] = (const real(*)[3]) gr->x;
  uint      tetrahedra = gr->tetrahedra;
  const bool *is_dof = dc->is_dof;
  const uint *idx2dof = dc->idx2dof;
  uint      ndof = dc->ndof;
  uint      nfix = dc->nfix;
  real      xt[4][3];
  real      g0[3];
  real      adet;
  uint      v[4];
  uint      t, i, j, ii, jj;

  for (t = 0; t < tetrahedra; t++) {
    /* Get vertices */
    getvertices_tet3d(gr, t, v);

    /* Get vertex coordinates */
    for (i = 0; i < 4; i++) {
      xt[i][0] = x[v[i]][0];
      xt[i][1] = x[v[i]][1];
      xt[i][2] = x[v[i]][2];
    }

    g0[0] = ((xt[2][1] - xt[1][1]) * (xt[3][2] - xt[1][2])
	     - (xt[2][2] - xt[1][2]) * (xt[3][1] - xt[1][1]));
    g0[1] = ((xt[2][2] - xt[1][2]) * (xt[3][0] - xt[1][0])
	     - (xt[2][0] - xt[1][0]) * (xt[3][2] - xt[1][2]));
    g0[2] = ((xt[2][0] - xt[1][0]) * (xt[3][1] - xt[1][1])
	     - (xt[2][1] - xt[1][1]) * (xt[3][0] - xt[1][0]));

    adet = fabs((xt[0][0] - xt[1][0]) * g0[0]
		+ (xt[0][1] - xt[1][1]) * g0[1]
		+ (xt[0][2] - xt[1][2]) * g0[2]);

    /* Add to system matrix */
    if (M)
      for (i = 0; i < 4; i++)
	if (is_dof[v[i]]) {
	  ii = idx2dof[v[i]];

	  for (j = 0; j < 4; j++)
	    if (is_dof[v[j]]) {
	      jj = idx2dof[v[j]];
	      addentry_sparsematrix(M, ii, jj, adet * Mt[i][j]);
	    }
	}

    /* Add to interaction matrix */
    if (Mf)
      for (i = 0; i < 4; i++)
	if (is_dof[v[i]]) {
	  ii = idx2dof[v[i]];
	  assert(ii < ndof);

	  for (j = 0; j < 4; j++)
	    if (!is_dof[v[j]]) {
	      jj = idx2dof[v[j]];
	      assert(jj < nfix);

	      addentry_sparsematrix(Mf, ii, jj, adet * Mt[i][j]);
	    }
	}
  }
}

void
assemble_tet3dp1_dirichlet_avector(pctet3dp1 dc,
				   field(*d) (const real * x, void *fdata),
				   void *fdata, pavector dv)
{
  const     real(*x)[3] = (const real(*)[3]) dc->gr->x;
  const bool *is_dof = dc->is_dof;
  const uint *idx2dof = dc->idx2dof;
  uint      vertices = dc->gr->vertices;
  uint      nfix = dc->nfix;
  uint      i, ii;

  assert(dv->dim == nfix);

  for (i = 0; i < vertices; i++)
    if (!is_dof[i]) {
      ii = idx2dof[i];
      assert(ii < nfix);

      dv->v[ii] = (d ? d(x[i], fdata) : 0.0);
    }
}

void
assemble_tet3dp1_functional_avector(pctet3dp1 dc,
				    field(*f) (const real * x, void *fdata),
				    void *fdata, pavector fv)
{
  pctet3d   gr = dc->gr;
  const     real(*x)[3] = (const real(*)[3]) gr->x;
  const bool *is_dof = dc->is_dof;
  const uint *idx2dof = dc->idx2dof;
  uint      tetrahedra = dc->gr->tetrahedra;
  uint      ndof = dc->ndof;
  real      xt[4][3], xd[3][3], xm[3], vol;
  field     f01, f02, f03, f12, f13, f23;
  uint      i, ii, t, v[4];

  assert(fv->dim == ndof);

  clear_avector(fv);

  for (t = 0; t < tetrahedra; t++) {
    /* Get vertices */
    getvertices_tet3d(gr, t, v);

    /* Get vertex coordinates */
    for (i = 0; i < 4; i++) {
      xt[i][0] = x[v[i]][0];
      xt[i][1] = x[v[i]][1];
      xt[i][2] = x[v[i]][2];
    }

    /* Compute difference vectors */
    for (i = 0; i < 3; i++) {
      xd[i][0] = xt[i + 1][0] - xt[0][0];
      xd[i][1] = xt[i + 1][1] - xt[0][1];
      xd[i][2] = xt[i + 1][2] - xt[0][2];
    }

    /* Compute volume by Sarrus' rule */
    vol = fabs(xd[0][0] * xd[1][1] * xd[2][2]
	       + xd[1][0] * xd[2][1] * xd[0][2]
	       + xd[2][0] * xd[0][1] * xd[1][2]
	       - xd[0][0] * xd[2][1] * xd[1][2]
	       - xd[1][0] * xd[0][1] * xd[2][2]
	       - xd[2][0] * xd[1][1] * xd[0][2]) / 6.0;

    /* Evaluate function in midpoints */
    xm[0] = 0.5 * (xt[0][0] + xt[1][0]);
    xm[1] = 0.5 * (xt[0][1] + xt[1][1]);
    xm[2] = 0.5 * (xt[0][2] + xt[1][2]);
    f01 = f(xm, fdata);
    xm[0] = 0.5 * (xt[0][0] + xt[2][0]);
    xm[1] = 0.5 * (xt[0][1] + xt[2][1]);
    xm[2] = 0.5 * (xt[0][2] + xt[2][2]);
    f02 = f(xm, fdata);
    xm[0] = 0.5 * (xt[0][0] + xt[3][0]);
    xm[1] = 0.5 * (xt[0][1] + xt[3][1]);
    xm[2] = 0.5 * (xt[0][2] + xt[3][2]);
    f03 = f(xm, fdata);
    xm[0] = 0.5 * (xt[1][0] + xt[2][0]);
    xm[1] = 0.5 * (xt[1][1] + xt[2][1]);
    xm[2] = 0.5 * (xt[1][2] + xt[2][2]);
    f12 = f(xm, fdata);
    xm[0] = 0.5 * (xt[1][0] + xt[3][0]);
    xm[1] = 0.5 * (xt[1][1] + xt[3][1]);
    xm[2] = 0.5 * (xt[1][2] + xt[3][2]);
    f13 = f(xm, fdata);
    xm[0] = 0.5 * (xt[2][0] + xt[3][0]);
    xm[1] = 0.5 * (xt[2][1] + xt[3][1]);
    xm[2] = 0.5 * (xt[2][2] + xt[3][2]);
    f23 = f(xm, fdata);

    /* Approximate integrals by midpoint rule */
    if (is_dof[v[0]]) {
      ii = idx2dof[v[0]];
      assert(ii < ndof);

      fv->v[ii] += 0.5 * (f01 + f02 + f03) * vol / 6.0;
    }
    if (is_dof[v[1]]) {
      ii = idx2dof[v[1]];
      assert(ii < ndof);

      fv->v[ii] += 0.5 * (f01 + f12 + f13) * vol / 6.0;
    }
    if (is_dof[v[2]]) {
      ii = idx2dof[v[2]];
      assert(ii < ndof);

      fv->v[ii] += 0.5 * (f02 + f12 + f23) * vol / 6.0;
    }
    if (is_dof[v[3]]) {
      ii = idx2dof[v[3]];
      assert(ii < ndof);

      fv->v[ii] += 0.5 * (f03 + f13 + f23) * vol / 6.0;
    }
  }
}

real
normmax_tet3dp1(pctet3dp1 dc,
		field(*u) (const real * x, void *fdata), void *fdata,
		pavector xs)
{
  const     real(*x)[3] = (const real(*)[3]) dc->gr->x;
  const bool *is_dof = dc->is_dof;
  const uint *idx2dof = dc->idx2dof;
  uint      vertices = dc->gr->vertices;
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
norml2_tet3dp1(pctet3dp1 dc,
	       field(*f) (const real * x, void *fdata), void *fdata,
	       pcavector xs, pcavector xf)
{
  const     real(*x)[3] = (const real(*)[3]) dc->gr->x;
  const     uint(*e)[2] = (const uint(*)[2]) dc->gr->e;
  const bool *is_dof = dc->is_dof;
  const uint *idx2dof = dc->idx2dof;
  uint      tetrahedra = dc->gr->tetrahedra;
  uint      ndof = dc->ndof;
  uint      nfix = dc->nfix;
  real      xt[4][3], xd[3][3], xm[3], vol;
  real      val, norm, tnorm;
  uint      vl[4], el[6];
  uint      t, i, ii, j;

  assert(xs == 0 || xs->dim == ndof);
  assert(xf == 0 || xf->dim == nfix);

  norm = 0.0;

  for (t = 0; t < tetrahedra; t++) {
    /* Get vertices (to compute Jacobi determinant) */
    getvertices_tet3d(dc->gr, t, vl);

    /* Get vertex coordinates */
    for (i = 0; i < 4; i++) {
      xt[i][0] = x[vl[i]][0];
      xt[i][1] = x[vl[i]][1];
      xt[i][2] = x[vl[i]][2];
    }

    /* Compute difference vectors */
    for (i = 0; i < 3; i++) {
      xd[i][0] = xt[i + 1][0] - xt[0][0];
      xd[i][1] = xt[i + 1][1] - xt[0][1];
      xd[i][2] = xt[i + 1][2] - xt[0][2];
    }

    /* Compute volume by Sarrus' rule */
    vol = fabs(xd[0][0] * xd[1][1] * xd[2][2]
	       + xd[1][0] * xd[2][1] * xd[0][2]
	       + xd[2][0] * xd[0][1] * xd[1][2]
	       - xd[0][0] * xd[2][1] * xd[1][2]
	       - xd[1][0] * xd[0][1] * xd[2][2]
	       - xd[2][0] * xd[1][1] * xd[0][2]) / 6.0;

    /* Approximate integral by midpoint quadrature */
    getedges_tet3d(dc->gr, t, el);
    tnorm = 0.0;
    for (j = 0; j < 6; j++) {
      val = 0.0;

      /* Evaluate in start point */
      i = e[el[j]][0];
      if (is_dof[i]) {
	if (xs) {
	  ii = idx2dof[i];
	  assert(ii < ndof);

	  val += 0.5 * xs->v[ii];
	}
      }
      else {
	if (xf) {
	  ii = idx2dof[i];
	  assert(ii < nfix);

	  val += 0.5 * xf->v[ii];
	}
      }

      /* Contribution to midpoint coordinates */
      xm[0] = 0.5 * x[i][0];
      xm[1] = 0.5 * x[i][1];
      xm[2] = 0.5 * x[i][2];

      /* Evaluate in end point */
      i = e[el[j]][1];
      if (is_dof[i]) {
	if (xs) {
	  ii = idx2dof[i];
	  assert(ii < ndof);

	  val += 0.5 * xs->v[ii];
	}
      }
      else {
	if (xf) {
	  ii = idx2dof[i];
	  assert(ii < nfix);

	  val += 0.5 * xf->v[ii];
	}
      }

      /* Find midpoint coordinates */
      xm[0] += 0.5 * x[i][0];
      xm[1] += 0.5 * x[i][1];
      xm[2] += 0.5 * x[i][2];

      /* Subtract reference function */
      if (f)
	val -= f(xm, fdata);

      tnorm += ABSSQR(val);
    }

    norm += tnorm * vol / 6.0;
  }

  return REAL_SQRT(norm);
}

pclustergeometry
build_tet3dp1_clustergeometry(pctet3dp1 p1, uint * idx)
{

  pclustergeometry cg;
  uint      i, j, k, c;
  uint      ndof = p1->ndof;
  pctet3d   t3 = p1->gr;
  uint      vertices = t3->vertices;
  uint      tetrahedra = t3->tetrahedra;
  const bool *is_dof = p1->is_dof;
  const uint *idx2dof = p1->idx2dof;
  const     real(*x)[3] = (const real(*)[3]) p1->gr->x;
  uint      v[4];
  real      min[3], max[3];

  cg = new_clustergeometry(3, ndof);

  /* Copying characteristic vertices into clustergeometry structure and 
   *computing initial values of smin and smax */
  for (i = 0; i < vertices; i++) {
    if (is_dof[i]) {
      c = idx2dof[i];

      cg->x[c][0] = x[i][0];
      cg->x[c][1] = x[i][1];
      cg->x[c][2] = x[i][2];

      cg->smin[c][0] = x[i][0];
      cg->smin[c][1] = x[i][1];
      cg->smin[c][2] = x[i][2];
      cg->smax[c][0] = x[i][0];
      cg->smax[c][1] = x[i][1];
      cg->smax[c][2] = x[i][2];

      idx[c] = c;
    }
  }

  for (i = 0; i < tetrahedra; i++) {	/* all tetrahedra */
    getvertices_tet3d(t3, i, v);	/* get all vertices of tetrahedron i */

    min[0] = REAL_MIN(REAL_MIN(x[v[0]][0], x[v[1]][0]),
		      REAL_MIN(x[v[2]][0], x[v[3]][0]));
    min[1] = REAL_MIN(REAL_MIN(x[v[0]][1], x[v[1]][1]),
		      REAL_MIN(x[v[2]][1], x[v[3]][1]));
    min[2] = REAL_MIN(REAL_MIN(x[v[0]][2], x[v[1]][2]),
		      REAL_MIN(x[v[2]][2], x[v[3]][2]));
    max[0] = REAL_MAX(REAL_MAX(x[v[0]][0], x[v[1]][0]),
		      REAL_MAX(x[v[2]][0], x[v[3]][0]));
    max[1] = REAL_MAX(REAL_MAX(x[v[0]][1], x[v[1]][1]),
		      REAL_MAX(x[v[2]][1], x[v[3]][1]));
    max[2] = REAL_MAX(REAL_MAX(x[v[0]][2], x[v[1]][2]),
		      REAL_MAX(x[v[2]][2], x[v[3]][2]));

    for (j = 0; j < 4; j++) {	/* all vertices */
      if (is_dof[v[j]]) {	/* vertex is degree of freedom */
	k = idx2dof[v[j]];

	cg->smin[k][0] = REAL_MIN(cg->smin[k][0], min[0]);
	cg->smin[k][1] = REAL_MIN(cg->smin[k][1], min[1]);
	cg->smin[k][2] = REAL_MIN(cg->smin[k][2], min[2]);
	cg->smax[k][0] = REAL_MAX(cg->smax[k][0], max[0]);
	cg->smax[k][1] = REAL_MAX(cg->smax[k][1], max[1]);
	cg->smax[k][2] = REAL_MAX(cg->smax[k][2], max[2]);
      }
    }
  }

  update_point_bbox_clustergeometry(cg, ndof, idx);

  return cg;
}
