#include "tet3drt0.h"

#include <stdio.h>
#include <assert.h>
#include <math.h>

#include "basic.h"

ptet3drt0 new_tet3drt0(pctet3d t3) {
  ptet3drt0 dc;

  uint *is_dof;
  uint *idx2dof;
  uint ndof, nfix;
  uint i;

  dc = (ptet3drt0) allocmem(sizeof(tet3drt0));
  is_dof = dc->is_dof = (uint *) allocmem((size_t ) sizeof(uint) * t3->faces);
  idx2dof = dc->idx2dof = (uint *) allocuint(t3->faces);
  dc->t3 = t3;

  for (i = 0; i < t3->faces; i++) { /* Find d.o.f.s */
    if (t3->fb[i] == 0)
      is_dof[i] = 0;
    else
      is_dof[i] = 1;
  }

  ndof = 0;
  nfix = 0;
  for (i = 0; i < t3->faces; i++) {
    if (is_dof[i] == 0)
      idx2dof[i] = ndof++;
    else
      idx2dof[i] = nfix++;
  }

  dc->ndof = ndof;
  dc->nfix = nfix;

  return dc;
}

void update_tet3drt0(ptet3drt0 rt0) {

  pctet3d t3 = rt0->t3;
  uint *is_dof = rt0->is_dof;
  uint *idx2dof = rt0->idx2dof;
  uint i, nfix, ndof;

  nfix = 0;
  ndof = 0;
  for (i = 0; i < t3->faces; i++) {
    if (is_dof[i] == 0 || is_dof[i] == 1) {
      idx2dof[i] = ndof++; // printf("0/1 \t");
    }
    //else if (is_dof[i] == 1){
    //  idx2dof[i] = ndof++; printf("1 \t");
    //}
    else {
      assert(is_dof[i] == 2);
      idx2dof[i] = nfix++;   //printf("2 \t");
    }
  }
  rt0->ndof = ndof;
  rt0->nfix = nfix;
}

void del_tet3drt0(ptet3drt0 dc) {
  freemem(dc->is_dof);
  freemem(dc->idx2dof);
  freemem(dc);
}

//psparsematrix build_tet3drt0_A_sparsematrix(pctet3drt0 dc) {
//  pctet3d t3 = dc->t3;
//  const uint *is_dof = dc->is_dof;
//  const uint *idx2dof = dc->idx2dof;
//  uint ndof = dc->ndof;
//  uint tetrahedra = dc->t3->tetrahedra;
//  uint faces = dc->t3->faces;
//  psparsepattern sp;
//  psparsematrix A;
//  uint i, j, ii, jj, d, f, ft[4];
//
//  sp = new_sparsepattern(ndof, ndof);
//
//  for (d = 0; d < tetrahedra; d++) {
//    /* Get faces */
//    ft[0] = t3->t[d][0];
//    ft[1] = t3->t[d][1];
//    ft[2] = t3->t[d][2];
//    ft[3] = t3->t[d][3];
//    /* Build sparsepattern */
//    for (i = 0; i < 4; i++) {
//      if (is_dof[ft[i]] == 0 || is_dof[ft[i]] == 1) { /*inner or Dirichlet edge*/
//        ii = idx2dof[ft[i]];
//        for (j = i + 1; j < 4; j++) {
//          if (is_dof[ft[j]] == 0 || is_dof[ft[j]] == 1) { /*inner or Dirichlet*/
//            jj = idx2dof[ft[j]];
//            addnz_sparsepattern(sp, ii, jj);
//            addnz_sparsepattern(sp, jj, ii);
//          }
//        }
//      }
//    }
//  }
//
//  for (f = 0; f < faces; ++f) {
//    if (is_dof[f] == 0 || is_dof[f] == 1) { /*inner or Dirichlet edge*/
//      ii = idx2dof[f];
//      addnz_sparsepattern(sp, ii, ii);
//    }
//  }
//
//  /* Build sparsematrix from sparsepattern */
//  A = new_zero_sparsematrix(sp);
//
//  del_sparsepattern(sp);
//
//  return A;
//}

psparsematrix build_tet3drt0_A_sparsematrix(pctet3drt0 dc) {
  pctet3d t3 = dc->t3;
  const uint *is_dof = dc->is_dof;
  const uint *idx2dof = dc->idx2dof;
  uint ndof = dc->ndof;
  uint tetrahedra = dc->t3->tetrahedra;
  uint faces = dc->t3->faces;
  psparsematrix A;
  uint *row_nnz;
  uint i, j, ii, jj, d, f, ft[4], nnz, nnz2, is_dof_i, is_dof_j;

  A = new_raw_sparsematrix(ndof, ndof, 0);
  row_nnz = allocuint(ndof);

  for (f = 0; f < ndof; ++f) {
    row_nnz[f] = 0;
  }
  nnz = 0;

  /****************************************************
   * count non-zeroes in off-diagonal part
   ****************************************************/

  for (d = 0; d < tetrahedra; d++) {
    /* Get faces */
    ft[0] = t3->t[d][0];
    ft[1] = t3->t[d][1];
    ft[2] = t3->t[d][2];
    ft[3] = t3->t[d][3];
    for (i = 0; i < 4; i++) {
      is_dof_i = is_dof[ft[i]];
      if (is_dof_i < 2) { /*inner or Dirichlet edge*/
        ii = idx2dof[ft[i]];
        for (j = i + 1; j < 4; j++) {
          is_dof_j = is_dof[ft[j]];
          if (is_dof_j < 2) { /*inner or Dirichlet*/
            jj = idx2dof[ft[j]];
            row_nnz[ii]++;
            row_nnz[jj]++;
            nnz += 2;
          }
        }
      }
    }
  }

  /****************************************************
   * count non-zeroes in diagonal part
   ****************************************************/

  for (f = 0; f < faces; ++f) {
    is_dof_i = is_dof[f];
    if (is_dof_i < 2) { /*inner or Dirichlet edge*/
      ii = idx2dof[f];
      row_nnz[ii]++;
      nnz++;
    }
  }

  /****************************************************
   * update A->row array
   ****************************************************/

  A->row[0] = 0;
  nnz2 = 0;
  for (f = 1; f <= ndof; ++f) {
    nnz2 += row_nnz[f - 1];
    A->row[f] = nnz2;
  }
  assert(nnz2 == nnz);

  for (f = 0; f < ndof; ++f) {
    row_nnz[f] = A->row[f];
  }

  /****************************************************
   * Allocate memory for remaining fields of A
   ****************************************************/

  freemem(A->coeff); freemem(A->col);
  
  A->coeff = allocfield(nnz);
  for (f = 0; f < nnz; ++f) {
    A->coeff[f] = 0.0;
  }
  A->col = allocuint(nnz);
  A->nz = nnz;

  /****************************************************
   * Insert entries in off-diagonal part
   ****************************************************/

  for (d = 0; d < tetrahedra; d++) {
    /* Get faces */
    ft[0] = t3->t[d][0];
    ft[1] = t3->t[d][1];
    ft[2] = t3->t[d][2];
    ft[3] = t3->t[d][3];
    for (i = 0; i < 4; i++) {
      is_dof_i = is_dof[ft[i]];
      if (is_dof_i < 2) { /*inner or Dirichlet edge*/
        ii = idx2dof[ft[i]];
        for (j = i + 1; j < 4; j++) {
          is_dof_j = is_dof[ft[j]];
          if (is_dof_j < 2) { /*inner or Dirichlet*/
            jj = idx2dof[ft[j]];
            A->col[row_nnz[ii]++] = jj;
            A->col[row_nnz[jj]++] = ii;
          }
        }
      }
    }
  }

  /****************************************************
   * Insert entries in off-diagonal part
   ****************************************************/

  for (f = 0; f < faces; ++f) {
    is_dof_i = is_dof[f];
    if (is_dof_i < 2) { /*inner or Dirichlet edge*/
      ii = idx2dof[f];
      A->col[row_nnz[ii]++] = ii;
    }
  }

  sort_sparsematrix(A);

  freemem(row_nnz);

  return A;
}

//psparsematrix build_tet3drt0_A_interaction_sparsematrix(pctet3drt0 dc) {
//  pctet3d t3 = dc->t3;
//  const uint *is_dof = dc->is_dof;
//  const uint *idx2dof = dc->idx2dof;
//  uint ndof = dc->ndof;
//  uint nfix = dc->nfix;
//  psparsepattern sp;
//  psparsematrix Af;
//  uint i, j, ii, jj, t, ft[4];
//
//  sp = new_sparsepattern(ndof, nfix);
//
//  for (t = 0; t < t3->tetrahedra; t++) {
//    /* Get faces */
//    ft[0] = t3->t[t][0];
//    ft[1] = t3->t[t][1];
//    ft[2] = t3->t[t][2];
//    ft[3] = t3->t[t][3];
//    /* Build sparsepattern */
//    for (i = 0; i < 4; i++) {
//      if (is_dof[ft[i]] == 0 || is_dof[ft[i]] == 1) { /*inner or Dirichlet*/
//        ii = idx2dof[ft[i]];
//        assert(ii < ndof);
//        for (j = 0; j < 4; j++) {
//          if (is_dof[ft[j]] == 2) { /*Neumann*/
//            jj = idx2dof[ft[j]];
//            assert(jj < nfix);
//            addnz_sparsepattern(sp, ii, jj);
//          }
//        }
//      }
//    }
//  }
//
//  /* Build sparsematrix from sparsepattern */
//  Af = new_zero_sparsematrix(sp);
//
//  del_sparsepattern(sp);
//
//  return Af;
//}

psparsematrix build_tet3drt0_A_interaction_sparsematrix(pctet3drt0 dc) {
  pctet3d t3 = dc->t3;
  const uint *is_dof = dc->is_dof;
  const uint *idx2dof = dc->idx2dof;
  uint ndof = dc->ndof;
  uint nfix = dc->nfix;
  psparsematrix Af;
  uint *row_nnz;
  uint i, j, ii, jj, t, f, ft[4], nnz, nnz2, is_dof_i, is_dof_j;

  Af = new_raw_sparsematrix(ndof, nfix, 0);
  row_nnz = allocuint(ndof);

  for (f = 0; f < ndof; ++f) {
    row_nnz[f] = 0;
  }
  nnz = 0;

  /****************************************************
   * count non-zeroes
   ****************************************************/

  for (t = 0; t < t3->tetrahedra; t++) {
    /* Get faces */
    ft[0] = t3->t[t][0];
    ft[1] = t3->t[t][1];
    ft[2] = t3->t[t][2];
    ft[3] = t3->t[t][3];
    for (i = 0; i < 4; i++) {
      is_dof_i = is_dof[ft[i]];
      if (is_dof_i < 2) { /*inner or Dirichlet*/
        ii = idx2dof[ft[i]];
        assert(ii < ndof);
        for (j = 0; j < 4; j++) {
          is_dof_j = is_dof[ft[j]];
          if (is_dof_j == 2) { /*Neumann*/
            jj = idx2dof[ft[j]];
            assert(jj < nfix);
            row_nnz[ii]++;
            nnz++;
          }
        }
      }
    }
  }

  /****************************************************
   * update Af->row array
   ****************************************************/

  Af->row[0] = 0;
  nnz2 = 0;
  for (f = 1; f <= ndof; ++f) {
    nnz2 += row_nnz[f - 1];
    Af->row[f] = nnz2;
  }
  assert(nnz2 == nnz);

  for (f = 0; f < ndof; ++f) {
    row_nnz[f] = Af->row[f];
  }

  /****************************************************
   * Allocate memory for remaining fields of Af
   ****************************************************/

  freemem(Af->coeff); freemem(Af->col);
  
  Af->coeff = allocfield(nnz);
  for (f = 0; f < nnz; ++f) {
    Af->coeff[f] = 0.0;
  }
  Af->col = allocuint(nnz);
  Af->nz = nnz;

  /****************************************************
   * Insert entries
   ****************************************************/

  for (t = 0; t < t3->tetrahedra; t++) {
    /* Get faces */
    ft[0] = t3->t[t][0];
    ft[1] = t3->t[t][1];
    ft[2] = t3->t[t][2];
    ft[3] = t3->t[t][3];
    for (i = 0; i < 4; i++) {
      is_dof_i = is_dof[ft[i]];
      if (is_dof_i < 2) { /*inner or Dirichlet*/
        ii = idx2dof[ft[i]];
        assert(ii < ndof);
        for (j = 0; j < 4; j++) {
          is_dof_j = is_dof[ft[j]];
          if (is_dof_j == 2) { /*Neumann*/
            jj = idx2dof[ft[j]];
            assert(jj < nfix);
            Af->col[row_nnz[ii]++] = jj;
          }
        }
      }
    }
  }

  sort_sparsematrix(Af);

  freemem(row_nnz);

  return Af;
}

//psparsematrix build_tet3drt0_B_sparsematrix(pctet3drt0 dc) {
//  pctet3d t3 = dc->t3;
//  const uint *is_dof = dc->is_dof;
//  const uint *idx2dof = dc->idx2dof;
//  uint ndof = dc->ndof;
//  uint tetrahedra = dc->t3->tetrahedra;
//  psparsepattern sp;
//  psparsematrix A;
//  uint i, ii, d, ft[4], nt;
//
//  nt = dc->t3->tetrahedra;
//  sp = new_sparsepattern(nt, ndof);
//
//  for (d = 0; d < tetrahedra; d++) {
//    /* Get faces */
//    ft[0] = t3->t[d][0];
//    ft[1] = t3->t[d][1];
//    ft[2] = t3->t[d][2];
//    ft[3] = t3->t[d][3];
//    /* Build sparsepattern */
//    for (i = 0; i < 4; i++) {
//      if (is_dof[ft[i]] == 0 || is_dof[ft[i]] == 1) { /*inner or Dirichlet*/
//        ii = idx2dof[ft[i]];
//        addnz_sparsepattern(sp, d, ii);
//      }
//    }
//  }
//
//  /* Build sparsematrix from sparsepattern */
//  A = new_zero_sparsematrix(sp);
//
//  del_sparsepattern(sp);
//
//  return A;
//}

psparsematrix build_tet3drt0_B_sparsematrix(pctet3drt0 dc) {
  pctet3d t3 = dc->t3;
  const uint *is_dof = dc->is_dof;
  const uint *idx2dof = dc->idx2dof;
  uint ndof = dc->ndof;
  uint tetrahedra = dc->t3->tetrahedra;
  psparsematrix A;
  uint *row_nnz;
  uint i, ii, d, ft[4], f, nt, nnz, nnz2, is_dof_i;

  nt = dc->t3->tetrahedra;

  A = new_raw_sparsematrix(nt, ndof, 0);
  row_nnz = allocuint(nt);

  for (f = 0; f < nt; ++f) {
    row_nnz[f] = 0;
  }
  nnz = 0;

  /****************************************************
   * count non-zeroes
   ****************************************************/

  for (d = 0; d < tetrahedra; d++) {
    /* Get faces */
    ft[0] = t3->t[d][0];
    ft[1] = t3->t[d][1];
    ft[2] = t3->t[d][2];
    ft[3] = t3->t[d][3];
    for (i = 0; i < 4; i++) {
      is_dof_i = is_dof[ft[i]];
      if (is_dof_i < 2) { /*inner or Dirichlet*/
//        ii = idx2dof[ft[i]];
        row_nnz[d]++;
        nnz++;
      }
    }
  }

  /****************************************************
   * update A->row array
   ****************************************************/

  A->row[0] = 0;
  nnz2 = 0;
  for (f = 1; f <= nt; ++f) {
    nnz2 += row_nnz[f - 1];
    A->row[f] = nnz2;
  }
  assert(nnz2 == nnz);

  for (f = 0; f < nt; ++f) {
    row_nnz[f] = A->row[f];
  }

  /****************************************************
   * Allocate memory for remaining fields of A
   ****************************************************/

  freemem(A->coeff); freemem(A->col);
  
  A->coeff = allocfield(nnz);
  for (f = 0; f < nnz; ++f) {
    A->coeff[f] = 0.0;
  }
  A->col = allocuint(nnz);
  A->nz = nnz;

  /****************************************************
   * Insert entries
   ****************************************************/

  for (d = 0; d < tetrahedra; d++) {
    /* Get faces */
    ft[0] = t3->t[d][0];
    ft[1] = t3->t[d][1];
    ft[2] = t3->t[d][2];
    ft[3] = t3->t[d][3];
    for (i = 0; i < 4; i++) {
      is_dof_i = is_dof[ft[i]];
      if (is_dof_i < 2) { /*inner or Dirichlet*/
        ii = idx2dof[ft[i]];
        A->col[row_nnz[d]++] = ii;
      }
    }
  }

  sort_sparsematrix(A);

  freemem(row_nnz);

  return A;
}

//psparsematrix build_tet3drt0_B_interaction_sparsematrix(pctet3drt0 dc) {
//  pctet3d t3 = dc->t3;
//  const uint *is_dof = dc->is_dof;
//  const uint *idx2dof = dc->idx2dof;
//  // uint ndof = dc->ndof;
//  uint nfix = dc->nfix;
//  uint tetrahedra = dc->t3->tetrahedra;
//  psparsepattern sp;
//  psparsematrix A;
//  uint i, ii, d, ft[4], nt;
//
//  nt = dc->t3->tetrahedra;
//  sp = new_sparsepattern(nt, nfix);
//
//  for (d = 0; d < tetrahedra; d++) {
//    /* Get faces */
//    ft[0] = t3->t[d][0];
//    ft[1] = t3->t[d][1];
//    ft[2] = t3->t[d][2];
//    ft[3] = t3->t[d][3];
//    /* Build sparsepattern */
//    for (i = 0; i < 4; i++) {
//      if (is_dof[ft[i]] == 2) { /*Neumann*/
//        ii = idx2dof[ft[i]];
//        addnz_sparsepattern(sp, d, ii);
//      }
//    }
//  }
//
//  /* Build sparsematrix from sparsepattern */
//  A = new_zero_sparsematrix(sp);
//
//  del_sparsepattern(sp);
//
//  return A;
//}

psparsematrix build_tet3drt0_B_interaction_sparsematrix(pctet3drt0 dc) {
  pctet3d t3 = dc->t3;
  const uint *is_dof = dc->is_dof;
  const uint *idx2dof = dc->idx2dof;
  // uint ndof = dc->ndof;
  uint nfix = dc->nfix;
  uint tetrahedra = dc->t3->tetrahedra;
  psparsematrix Af;
  uint *row_nnz;
  uint i, ii, d, ft[4], f, nt, nnz, nnz2, is_dof_i;

  nt = dc->t3->tetrahedra;
  Af = new_raw_sparsematrix(nt, nfix, 0);
  row_nnz = allocuint(nt);

  for (f = 0; f < nt; ++f) {
    row_nnz[f] = 0;
  }
  nnz = 0;

  /****************************************************
   * count non-zeroes
   ****************************************************/

  for (d = 0; d < tetrahedra; d++) {
    /* Get faces */
    ft[0] = t3->t[d][0];
    ft[1] = t3->t[d][1];
    ft[2] = t3->t[d][2];
    ft[3] = t3->t[d][3];
    for (i = 0; i < 4; i++) {
      is_dof_i = is_dof[ft[i]];
      if (is_dof_i == 2) { /*Neumann*/
//        ii = idx2dof[ft[i]];
        row_nnz[d]++;
        nnz++;
      }
    }
  }

  /****************************************************
   * update Af->row array
   ****************************************************/

  Af->row[0] = 0;
  nnz2 = 0;
  for (f = 1; f <= nt; ++f) {
    nnz2 += row_nnz[f - 1];
    Af->row[f] = nnz2;
  }
  assert(nnz2 == nnz);

  for (f = 0; f < nt; ++f) {
    row_nnz[f] = Af->row[f];
  }

  /****************************************************
   * Allocate memory for remaining fields of Af
   ****************************************************/

  freemem(Af->coeff); freemem(Af->col);
  
  Af->coeff = allocfield(nnz);
  for (f = 0; f < nnz; ++f) {
    Af->coeff[f] = 0.0;
  }
  Af->col = allocuint(nnz);
  Af->nz = nnz;

  /****************************************************
   * Insert entries
   ****************************************************/

  for (d = 0; d < tetrahedra; d++) {
    /* Get faces */
    ft[0] = t3->t[d][0];
    ft[1] = t3->t[d][1];
    ft[2] = t3->t[d][2];
    ft[3] = t3->t[d][3];
    for (i = 0; i < 4; i++) {
      is_dof_i = is_dof[ft[i]];
      if (is_dof_i == 2) { /*Neumann*/
        ii = idx2dof[ft[i]];
        Af->col[row_nnz[d]++] = ii;
      }
    }
  }

  sort_sparsematrix(Af);

  freemem(row_nnz);

  return Af;
}

#if 0

void get_edges_and_opposite_vertices_tri2d(pctri2d t2, uint d, uint et[], uint v[])
{
  const uint (*t)[3] = (const uint (*)[3]) t2->t;
  const uint (*e)[2] = (const uint (*)[2]) t2->e;
  uint xt[3];

  /*Get edges*/
  et[0] = t[d][0];
  et[1] = t[d][1];
  et[2] = t[d][2];
  /*Get vertices*/
  getvertices_tri2d(t2, d, xt);

  if( (e[et[0]][0] == xt[0] && e[et[0]][1] == xt[1]) || (e[et[0]][1] == xt[0] && e[et[0]][0] == xt[1]) )
  v[0] = xt[2];
  else if( (e[et[0]][0] == xt[0] && e[et[0]][1] == xt[2]) || (e[et[0]][0] == xt[2] && e[et[0]][1] == xt[0]) )
  v[0]=xt[1];
  else
  v[1] = xt[0];

  if( (e[et[1]][0] == xt[0] && e[et[1]][1] == xt[1]) || (e[et[1]][1] == xt[0] && e[et[1]][0] == xt[1]) )
  v[0] = xt[2];
  else if( (e[et[1]][0] == xt[0] && e[et[1]][1] == xt[2]) || (e[et[1]][0] == xt[2] && e[et[1]][1] == xt[0]) )
  v[0]=xt[1];
  else
  v[1] = xt[0];

  if( (e[et[2]][0] == xt[0] && e[et[2]][1] == xt[1]) || (e[et[2]][1] == xt[0] && e[et[2]][0] == xt[1]) )
  v[0] = xt[2];
  else if( (e[et[2]][0] == xt[0] && e[et[2]][1] == xt[2]) || (e[et[2]][0] == xt[2] && e[et[2]][1] == xt[0]) )
  v[0]=xt[1];
  else
  v[1] = xt[0];

}
#endif

static real
scalar(real x[3], real y[3]) {
  real s;
  s = x[0] * y[0] + x[1] * y[1] + x[2] * y[2];

  return s;
}

/*if clockwise normal is outer normal at face i in tetrahedra t: return  1
 else   -1 */
real compute_type_of_face_tet3d(pctet3d t3, uint d, uint i) {

  const real (*x)[3] = (const real (*)[3]) t3->x;
  real signum;
  real nr[3], a[3], b[3], c[3];
  uint v[4];

  /*Get vertices, v[0] is opposite of face i */
  getvertices_byface_tet3d(t3, d, i, v);
  /*Compute vectors spanning the plane, which contains face i*/
  a[0] = x[v[2]][0] - x[v[1]][0];
  a[1] = x[v[2]][1] - x[v[1]][1];
  a[2] = x[v[2]][2] - x[v[1]][2];
  b[0] = x[v[3]][0] - x[v[1]][0];
  b[1] = x[v[3]][1] - x[v[1]][1];
  b[2] = x[v[3]][2] - x[v[1]][2];
  /*Compute normal vector by cross product*/
  nr[0] = a[1] * b[2] - a[2] * b[1];
  nr[1] = a[2] * b[0] - a[0] * b[2];
  nr[2] = a[0] * b[1] - a[1] * b[0];

  c[0] = x[v[0]][0] - x[v[1]][0];
  c[1] = x[v[0]][1] - x[v[1]][1];
  c[2] = x[v[0]][2] - x[v[1]][2];

  if (nr[0] * c[0] + nr[1] * c[1] + nr[2] * c[2] > 0)
    signum = -1.0; /*c and nr lie on the same side of face i*/
  else
    signum = 1.0; /*c and nr lie on different sides of face i*/
  return signum;
}

void assemble_tet3drt0_darcy_A_sparsematrix(pctet3drt0 dc, psparsematrix A,
    psparsematrix Af, pavector K) {
  const real (*x)[3] = (const real (*)[3]) dc->t3->x;
  // const uint (*e)[2] = (const uint (*)[2]) dc->t3->e;
  const uint (*t)[4] = (const uint (*)[4]) dc->t3->t;
  const uint *is_dof = dc->is_dof;
  const uint *idx2dof = dc->idx2dof;
  uint tetrahedra = dc->t3->tetrahedra;
  uint ndof = dc->ndof;
  uint nfix = dc->nfix;
  uint i, j, ii, jj, d, k, l;

  real At[4][4];
  uint v[4];
  uint ft[4];
  // int sigma[3];
  real T, signum_i, signum_j;
  //real E[3],
  real x_ij[3], x_jk[3], /*x_kl[3],*/x_li[3], x_ik[3], x_jl[3];
  //real a,b;
  //int eti_end, eti_start;

  //tetrahedra =1;
  for (d = 0; d < tetrahedra; d++) {
    // printf("tetrahedra %d\n", d);
    /*Get faces of tetrahedra d*/
    ft[0] = t[d][0];
    ft[1] = t[d][1];
    ft[2] = t[d][2];
    ft[3] = t[d][3];
    // printf("faces %d %d %d %d\n", ft[0], ft[1], ft[2], ft[3]);
    /*Get vertices of tetrahedra d*/
    getvertices_tet3d(dc->t3, d, v);
    //     printf("vertices\n");
    // printf("%d %f %f %f\n", xt[0], x[xt[0]][0], x[xt[0]][1], x[xt[0]][2]);
    // printf("%d %f %f %f\n", xt[1], x[xt[1]][0], x[xt[1]][1], x[xt[1]][2]);
    // printf("%d %f %f %f\n", xt[2], x[xt[2]][0], x[xt[2]][1], x[xt[2]][2]);
    // printf("%d %f %f %f\n", xt[3], x[xt[3]][0], x[xt[3]][1], x[xt[3]][2]);

    /*volume of tetrahedra d*/
    T = (1.0 / 6.0) * fabs(
            (x[v[1]][0] - x[v[0]][0]) * ((x[v[2]][1] - x[v[0]][1])
                * (x[v[3]][2] - x[v[0]][2])
                                         - (x[v[2]][2] - x[v[0]][2]) * (x[v[3]][1]
                                             - x[v[0]][1]))

            - (x[v[1]][1] - x[v[0]][1]) * ((x[v[2]][0] - x[v[0]][0])
                * (x[v[3]][2] - x[v[0]][2])
                                           - (x[v[2]][2] - x[v[0]][2]) * (x[v[3]][0]
                                               - x[v[0]][0]))

            + (x[v[1]][2] - x[v[0]][2]) * ((x[v[2]][0] - x[v[0]][0])
                * (x[v[3]][1] - x[v[0]][1])
                                           - (x[v[2]][1] - x[v[0]][1]) * (x[v[3]][0]
                                               - x[v[0]][0])));
    //printf("T = %f \n",T);
    /*algebraic sign*/
    //sigma[0]=1; sigma[1]= 1;sigma[2]=1; /*VZ noch richtig bestimmen!!!!!!!!!!!!*/
    // for(i=0;i<3;i++){
    //  eti_start = e[et[i]][0]; eti_end = e[et[i]][1];
    // a = x[eti_start][0] - x[eti_end][0];
    // b = x[eti_start][1] - x[eti_end][1];
    /*length of edge i*/
    // E[i] = sqrt(a*a+b*b);
    //}
    /*Initialise At with zero*/
    for (i = 0; i < 4; i++)
      for (j = 0; j < 4; j++)
        At[i][j] = 0.0;

    /* Compute element matrix */
    for (i = 0; i < 4; i++)
      for (j = 0; j < 4; j++) {
        /*non diagonal entry*/
        if (i != j) {
          if ((i == 0 && j == 1) || (i == 1 && j == 0)) {
            k = 2;
            l = 3;
          } else if ((i == 0 && j == 2) || (i == 2 && j == 0)) {
            k = 1;
            l = 3;
          } else if ((i == 0 && j == 3) || (i == 3 && j == 0)) {
            k = 1;
            l = 2;
          } else if ((i == 1 && j == 2) || (i == 2 && j == 1)) {
            k = 0;
            l = 3;
          } else if ((i == 1 && j == 3) || (i == 3 && j == 1)) {
            k = 0;
            l = 2;
          } else {
            k = 0;
            l = 1;
          }
          // printf("i=%d, j=%d, k=%d, l=%d \n",i,j,k,l);

          x_ij[0] = x[v[i]][0] - x[v[j]][0];
          x_ij[1] = x[v[i]][1] - x[v[j]][1];
          x_ij[2] = x[v[i]][2] - x[v[j]][2];
          x_jk[0] = x[v[j]][0] - x[v[k]][0];
          x_jk[1] = x[v[j]][1] - x[v[k]][1];
          x_jk[2] = x[v[j]][2] - x[v[k]][2];
          //x_kl[0] = x[v[k]][0] - x[v[l]][0]; x_kl[1] = x[v[k]][1] - x[v[l]][1]; x_kl[2] = x[v[k]][2] - x[v[l]][2];
          x_li[0] = x[v[l]][0] - x[v[i]][0];
          x_li[1] = x[v[l]][1] - x[v[i]][1];
          x_li[2] = x[v[l]][2] - x[v[i]][2];
          x_ik[0] = x[v[i]][0] - x[v[k]][0];
          x_ik[1] = x[v[i]][1] - x[v[k]][1];
          x_ik[2] = x[v[i]][2] - x[v[k]][2];
          x_jl[0] = x[v[j]][0] - x[v[l]][0];
          x_jl[1] = x[v[j]][1] - x[v[l]][1];
          x_jl[2] = x[v[j]][2] - x[v[l]][2];

          At[i][j] = -2 * scalar(x_ij, x_ij) - 2 * scalar(x_ij, x_ik)
              + 2 * scalar(x_jk, x_ij) + 3 * scalar(x_jk, x_ik)
                     - 2 * scalar(x_jk, x_li)
                     + 2 * scalar(x_jl, x_ij) + 2 * scalar(x_jl, x_ik)
                     - 3 * scalar(x_jl, x_li)
                     + 2 * scalar(x_ij, x_li);
          signum_i = compute_type_of_face_tet3d(dc->t3, d, i);
          signum_j = compute_type_of_face_tet3d(dc->t3, d, j);
          At[i][j] = (signum_i * signum_j * At[i][j] * K->v[d]) / (324 * T); 
          //printf("At[%u][%u] = %f signum[%u] = %f signum[%u] = %f\n",i, j,  At[i][j], i, signum_i, j, signum_j);
        }

        else { /*diagonal entry*/
          assert(i == j);
          i = ((j + 1) % 4);
          k = ((j + 2) % 4);
          l = ((j + 3) % 4);

          x_ij[0] = x[v[i]][0] - x[v[j]][0];
          x_ij[1] = x[v[i]][1] - x[v[j]][1];
          x_ij[2] = x[v[i]][2] - x[v[j]][2];
          x_jk[0] = x[v[j]][0] - x[v[k]][0];
          x_jk[1] = x[v[j]][1] - x[v[k]][1];
          x_jk[2] = x[v[j]][2] - x[v[k]][2];
          //x_kl[0] = x[v[k]][0] - x[v[l]][0]; x_kl[1] = x[v[k]][1] - x[v[l]][1]; x_kl[2] = x[v[k]][2] - x[v[l]][2];
          x_li[0] = x[v[l]][0] - x[v[i]][0];
          x_li[1] = x[v[l]][1] - x[v[i]][1];
          x_li[2] = x[v[l]][2] - x[v[i]][2];
          x_ik[0] = x[v[i]][0] - x[v[k]][0];
          x_ik[1] = x[v[i]][1] - x[v[k]][1];
          x_ik[2] = x[v[i]][2] - x[v[k]][2];
          x_jl[0] = x[v[j]][0] - x[v[l]][0];
          x_jl[1] = x[v[j]][1] - x[v[l]][1];
          x_jl[2] = x[v[j]][2] - x[v[l]][2];

          At[j][j] = 3 * scalar(x_ij, x_ij) + 3 * scalar(x_jk, x_jk)
                     + 3 * scalar(x_jl, x_jl)
                     - 4 * scalar(x_ij, x_jk)
                     + 4 * scalar(x_jl, x_jk)
                     - 4 * scalar(x_jl, x_ij);
          At[j][j] = (At[j][j] * K->v[d]) / (324 * T); 

          i = j;
        }

      }
    /* Add to system matrix */
    for (i = 0; i < 4; i++) {
      if (is_dof[ft[i]] == 0 || is_dof[ft[i]] == 1) {
        ii = idx2dof[ft[i]];
        for (j = 0; j < 4; j++) {
          if (is_dof[ft[j]] == 0 || is_dof[ft[j]] == 1) {
            jj = idx2dof[ft[j]];
            addentry_sparsematrix(A, ii, jj, At[i][j]);
          }
        }
      }
    }

    /* Add to interaction matrix */
    if (Af) {
      for (i = 0; i < 4; i++) {
        if (is_dof[ft[i]] == 0 || is_dof[ft[i]] == 1) { /*inner or Dirichlet*/
          ii = idx2dof[ft[i]];
          assert(ii < ndof);
          for (j = 0; j < 4; j++) {
            if (is_dof[ft[j]] == 2) { /*Neumann*/
              jj = idx2dof[ft[j]];
              assert(jj < nfix);
              addentry_sparsematrix(Af, ii, jj, At[i][j]);
            }
          }
        }
      }
    }

  }
}

void assemble_tet3drt0_darcy_B_sparsematrix(pctet3drt0 dc, psparsematrix A,
    psparsematrix Af) {
  const uint (*t)[4] = (const uint (*)[4]) dc->t3->t;
  const uint *is_dof = dc->is_dof;
  const uint *idx2dof = dc->idx2dof;
  uint tetrahedra = dc->t3->tetrahedra;
  uint i, ii, d; //E_i;
  real At;
  uint ft[4];
  //int eti_start, eti_end;
  //real x1[2], x2[2];
  //real a,b, sigma_i;
  real signum;

  for (d = 0; d < tetrahedra; d++) {
    /*Get faces of tetrahedra d*/
    ft[0] = t[d][0];
    ft[1] = t[d][1];
    ft[2] = t[d][2];
    ft[3] = t[d][3];
    /* Compute element matrix 4 * 1*/
    for (i = 0; i < 4; i++) {
      /*Find direction of normal at face i*/
      signum = compute_type_of_face_tet3d(dc->t3, d, i);
      At = signum * -1.0;
      /*Add to system matrix*/
      if (is_dof[ft[i]] == 0 || is_dof[ft[i]] == 1) { /*face i is degree of freedom; inner or dirichlet*/
        ii = idx2dof[ft[i]];
        addentry_sparsematrix(A, d, ii, At);
      }
      /*Add to interaction matrix*/
      if (Af) {
        if (is_dof[ft[i]] == 2) { /*Neumann edge?*/
          ii = idx2dof[ft[i]];
          addentry_sparsematrix(Af, d, ii, At);
        }
      }
    }
  }
}


void assemble_tet3drt0_b_D_avector(pctet3drt0 dc,
    field (*d)(const real *e, void *fdata), void *fdata, pavector dv) {
  // const uint (*f)[3] = (const uint(*)[3]) dc->t3->f;
  pctet3d t3 = dc->t3;
  // const real (*x)[2] = (const real(*)[2]) dc->t2->x;
  const uint *is_dof = dc->is_dof;
  const uint *idx2dof = dc->idx2dof;
  uint faces = dc->t3->faces;
  //uint nfix = dc->nfix;
  uint i, ii;  //, sv, ev;
  //real start[2], end[2], m[2];
  uint v[3];
  real s[3], alpha;
  // assert(dv->dim == nfix);

  clear_avector(dv);
  for (i = 0; i < faces; i++) {
    if (is_dof[i] == 1) { /*Dirichlet*/
      ii = idx2dof[i];
      getvertices_face_tet3d(dc->t3, i, v);
      alpha = 1.0L / 3.0L;
      s[0] = alpha * (t3->x[v[0]][0] + t3->x[v[1]][0] + t3->x[v[2]][0]);
      s[1] = alpha * (t3->x[v[0]][1] + t3->x[v[1]][1] + t3->x[v[2]][1]);
      s[2] = alpha * (t3->x[v[0]][2] + t3->x[v[1]][2] + t3->x[v[2]][2]);
      dv->v[ii] = (d ? -d(s, fdata) : 0.0);
    }
  }
}

void assemble_tet3drt0_b_f_avector(ptet3drt0 dc,
    field (*f)(const real *x, void *fdata), void *fdata, pavector fv) {
  pctet3d t3 = dc->t3;
  const real (*x)[3] = (const real (*)[3]) t3->x;
  uint tetrahedra = dc->t3->tetrahedra;
  // uint ndof = dc->ndof;
  real xt[4][3], xm[3], vol, alpha;
  field f012, f123, f230, f301;
  uint i, d, v[4];

  //assert(fv->dim == ndof);

  clear_avector(fv);

  for (d = 0; d < tetrahedra; d++) {
    //printf("tetrahedra %u\t", d);
    /*Get vertices of tetrahedra d*/
    // printf("Get vertices\t");
    getvertices_tet3d(dc->t3, d, v);
    //     printf("vertices\n");
    // printf("%d %f %f %f\n", xt[0], x[xt[0]][0], x[xt[0]][1], x[xt[0]][2]);
    // printf("%d %f %f %f\n", xt[1], x[xt[1]][0], x[xt[1]][1], x[xt[1]][2]);
    // printf("%d %f %f %f\n", xt[2], x[xt[2]][0], x[xt[2]][1], x[xt[2]][2]);
    // printf("%d %f %f %f\n", xt[3], x[xt[3]][0], x[xt[3]][1], x[xt[3]][2]);
    /* Get vertex coordinates */
    // printf("Get vertex coordinates\t");
    for (i = 0; i < 4; i++) {
      xt[i][0] = x[v[i]][0];
      xt[i][1] = x[v[i]][1];
      xt[i][2] = x[v[i]][2];
    }
    /*volume of tetrahedra d*/
    // printf("vol\t");
    vol = (1.0 / 6.0)
        * fabs(
            (xt[1][0] - xt[0][0]) * ((xt[2][1] - xt[0][1])
                * (xt[3][2] - xt[0][2])
                                     - (xt[2][2] - xt[0][2]) * (xt[3][1]
                                         - xt[0][1]))

            - (xt[1][1] - xt[0][1]) * ((xt[2][0] - xt[0][0])
                * (xt[3][2] - xt[0][2])
                                       - (xt[2][2] - xt[0][2]) * (xt[3][0]
                                           - xt[0][0]))

            + (xt[1][2] - xt[0][2]) * ((xt[2][0] - xt[0][0])
                * (xt[3][1] - xt[0][1])
                                       - (xt[2][1] - xt[0][1]) * (xt[3][0]
                                           - xt[0][0])));
    /* Evaluate function in midpoints */
    //printf("Eval f\t");
    alpha = 1.0 / 3.0;
    xm[0] = alpha * (xt[0][0] + xt[1][0] + xt[2][0]);
    xm[1] = alpha * (xt[0][1] + xt[1][1] + xt[2][1]);
    xm[2] = alpha * (xt[0][2] + xt[1][2] + xt[2][2]);
    f012 = f(xm, fdata);
    xm[0] = alpha * (xt[1][0] + xt[2][0] + xt[3][0]);
    xm[1] = alpha * (xt[1][1] + xt[2][1] + xt[3][1]);
    xm[2] = alpha * (xt[1][2] + xt[2][2] + xt[3][2]);
    f123 = f(xm, fdata);
    xm[0] = alpha * (xt[2][0] + xt[3][0] + xt[0][0]);
    xm[1] = alpha * (xt[2][1] + xt[3][1] + xt[0][1]);
    xm[2] = alpha * (xt[2][2] + xt[3][2] + xt[0][2]);
    f230 = f(xm, fdata);
    xm[0] = alpha * (xt[3][0] + xt[0][0] + xt[1][0]);
    xm[1] = alpha * (xt[3][1] + xt[0][1] + xt[1][1]);
    xm[2] = alpha * (xt[3][2] + xt[0][2] + xt[1][2]);
    f301 = f(xm, fdata);
    //printf("Fill fv\n");
    fv->v[d] = -(vol / 4.0) * (f012 + f123 + f230 + f301);
  }
}

void assemble_tet3drt0_g_N_avector(ptet3drt0 dc,
    field (*f)(const uint *e, void *data), void *data, pavector g) {

  const uint (*e)[2] = (const uint (*)[2]) dc->t3->e;
  const real (*x)[3] = (const real (*)[3]) dc->t3->x;
  uint faces = dc->t3->faces;
  const uint *is_dof = dc->is_dof;
  uint i, d, j;
  real xt[3][3], xd[2][3], c[3], area;
  uint v[3];

  d = 0;
  for (i = 0; i < faces; i++) {
    // printf("face %u\t", i);

    if (is_dof[i] == 2) { /*Neumann face*/
      // printf("face %u\t", i);

      /*area of face i*/
      /*Get vertices of face i*/
      // printf("Get vertives of face\t");
      getvertices_face_tet3d(dc->t3, i, v);
      /*Get coordinates of the vertices of face i*/
      //printf("Get vertex coordinates\t");
      for (j = 0; j < 3; j++) {
        xt[j][0] = x[v[j]][0];
        xt[j][1] = x[v[j]][1];
        xt[j][2] = x[v[j]][2];
      }
      /*Compute difference vectors*/
      for (j = 0; j < 2; j++) {
        xd[j][0] = xt[j + 1][0] - xt[0][0];
        xd[j][1] = xt[j + 1][1] - xt[0][1];
        xd[j][2] = xt[j + 1][2] - xt[0][2];
      }
      c[0] = xd[0][1] * xd[1][2] - xd[0][2] * xd[1][1];
      c[1] = xd[0][2] * xd[1][0] - xd[0][0] * xd[1][2];
      c[2] = xd[0][0] * xd[1][1] - xd[0][1] * xd[1][0];
      area = REAL_SQRT(c[0] * c[0] + c[1] * c[1] + c[2] * c[2]);
      area = area / 2;
      // printf("Fill g\n");
      g->v[d] = (f ? (area * f(e[i], data)) : 0.0);
      d = d + 1;
    }
  }
}

real getvolume_tetrahedra_tet3d(pctet3d t3, uint t) {

  uint v[4], i;
  real vol, xt[4][3];
  const real (*x)[3] = (const real (*)[3]) t3->x;

  /*Get vertices of tetrahedra t*/
  getvertices_tet3d(t3, t, v);
  /* Get vertex coordinates */
  for (i = 0; i < 4; i++) {
    xt[i][0] = x[v[i]][0];
    xt[i][1] = x[v[i]][1];
    xt[i][2] = x[v[i]][2];
  }
  /*volume of tetrahedra t*/
  vol = (1.0 / 6.0)
      * fabs(
          (xt[1][0] - xt[0][0]) * ((xt[2][1] - xt[0][1]) * (xt[3][2] - xt[0][2])
              - (xt[2][2] - xt[0][2]) * (xt[3][1] - xt[0][1]))

          - (xt[1][1] - xt[0][1]) * ((xt[2][0] - xt[0][0])
              * (xt[3][2] - xt[0][2])
                                     - (xt[2][2] - xt[0][2]) * (xt[3][0]
                                         - xt[0][0]))

          + (xt[1][2] - xt[0][2]) * ((xt[2][0] - xt[0][0])
              * (xt[3][1] - xt[0][1])
                                     - (xt[2][1] - xt[0][1]) * (xt[3][0]
                                         - xt[0][0])));
  return vol;
}

/*Schwerpunktsregel*/
real norml2_pressure_centroid_tet3drt0(pctet3drt0 dc,
    field (*p)(const real * x, void *fdata), void *fdata, pcavector x2) {
  const real (*x)[3] = (const real (*)[3]) dc->t3->x;
  uint tetrahedra = dc->t3->tetrahedra;
  real s[3], xt[4][3];
  real val, norm, vol, alpha;
  uint v[4];
  uint i, k;

  norm = 0.0;

  for (k = 0; k < tetrahedra; k++) {

    /*Get vertices*/
    getvertices_tet3d(dc->t3, k, v);
    /* Get vertex coordinates */
    for (i = 0; i < 4; i++) {
      xt[i][0] = x[v[i]][0];
      xt[i][1] = x[v[i]][1];
      xt[i][2] = x[v[i]][2];
    }
    /*area of triangle k*/
    vol = getvolume_tetrahedra_tet3d(dc->t3, k);
    /* Approximate integral */
    alpha = 1.0 / 4.0;
    s[0] = alpha * (xt[0][0] + xt[1][0] + xt[2][0] + xt[3][0]);
    s[1] = alpha * (xt[0][1] + xt[1][1] + xt[2][1] + xt[3][1]);
    s[2] = alpha * (xt[0][2] + xt[1][2] + xt[2][2] + xt[3][2]);
    val = p(s, fdata) - x2->v[k];
    //printf("p ? %f x2 = %f val = %f\n", (double) p(s, fdata), (double)x2->v[k], val);
    val = val * val;

    norm = norm + vol * val;

  }

  return REAL_SQRT(norm);
}

/*Kantenmittelpunkt*/
real norml2_pressure_facemidpoint_tet3drt0(pctet3drt0 dc,
    field (*p)(const real * x, void *fdata), void *fdata, pcavector x2) {
  const real (*x)[3] = (const real (*)[3]) dc->t3->x;
  uint tetrahedra = dc->t3->tetrahedra;
  //const     uint(*t)[4] = (const uint(*)[4]) dc->t3->t;
  //const     uint *is_dof = dc->is_dof;
  real s[3];
  real val, norm, vol, alpha;
  //uint      fd[3];
  uint v[4];
  uint k, f;
  //uint l;
  real tnorm;

  norm = 0.0;

  for (k = 0; k < tetrahedra; k++) {

    /*volume of tetrahedra k*/
    vol = getvolume_tetrahedra_tet3d(dc->t3, k);

    /*Get faces of tetrahedra k*/
    //fd[0] = t[k][0];
    //fd[1] = t[k][1];
    //fd[2] = t[k][2];
    //fd[3] = t[k][3];
    tnorm = 0.0;
    /*Approximate integral*/
    for (f = 0; f <= 3; f++) {
      /*Get vertices of tetrahedra k, v[1], v[2], v[3] are the vertices of face f*/
      getvertices_byface_tet3d(dc->t3, k, f, v);
      /*Compute centroid of face f in T_k*/
      alpha = 1.0 / 3.0;
      s[0] = alpha * (x[v[1]][0] + x[v[2]][0] + x[v[3]][0]);
      s[1] = alpha * (x[v[1]][1] + x[v[2]][1] + x[v[3]][1]);
      s[2] = alpha * (x[v[1]][2] + x[v[2]][2] + x[v[3]][2]);
      val = p(s, fdata) - x2->v[k]; //printf("p ? %f x2 = %f val = %f\n", (double) p(s, fdata), (double)x2->v[k], val);
      //if(is_dof[fd[f]] == 0){/*inner edge*/
      //  for(l=0;l<tetrahedra;l++){
      //    if((fd[f] == t[l][0] ||fd[f] == t[l][1] || fd[f] == t[l][2]) && l!=k){
      //    assert(k!=l);
      //    val = val - x2->v[l];
      //      }
      //    }
      //  }
      tnorm = tnorm + val * val;
    }
    norm = norm + (vol / 4.0) * tnorm;
  }
  return REAL_SQRT(norm);
}

real norml2_flux_centroid_tet3drt0(pctet3drt0 dc,
    void (*q)(const real * x, void *fdata, pavector v), void *fdata,
    pcavector x1, pcavector g) {
  const real (*x)[3] = (const real (*)[3]) dc->t3->x;
  const uint (*t)[4] = (const uint (*)[4]) dc->t3->t;
  const uint *is_dof = dc->is_dof;
  const uint *idx2dof = dc->idx2dof;
  uint tetrahedra = dc->t3->tetrahedra;
  uint faces = dc->t3->faces;
  uint ndof = dc->ndof;
  uint nfix = dc->nfix;
  real xt[4][3], vol, s[3], alpha, p[3];
  real val[3], norm, tnorm, signum;
  uint fd[4];
  uint i, j, k, jj, num;
  uint v[4];
  pavector v_flux;

  assert(x1 == 0 || x1->dim == ndof);
  assert(g == 0 || g->dim == nfix);

  norm = 0.0;
  v_flux = new_avector(3);

  for (k = 0; k < tetrahedra; k++) {
    /*Area of triangle k*/    //printf("triangle %u\t", k);
    vol = getvolume_tetrahedra_tet3d(dc->t3, k);
    /*faces of tetrahedra k*/
    fd[0] = t[k][0];
    fd[1] = t[k][1];
    fd[2] = t[k][2];
    fd[3] = t[k][3];
    /*Get vertices*/
    getvertices_tet3d(dc->t3, k, v);
    /* Get vertex coordinates */
    for (i = 0; i < 4; i++) {
      xt[i][0] = x[v[i]][0];
      xt[i][1] = x[v[i]][1];
      xt[i][2] = x[v[i]][2];
    }
    /*Compute centroid*/
    alpha = 1.0 / 4.0;
    s[0] = alpha * (xt[0][0] + xt[1][0] + xt[2][0] + xt[3][0]);
    s[1] = alpha * (xt[0][1] + xt[1][1] + xt[2][1] + xt[3][1]);
    s[2] = alpha * (xt[0][2] + xt[1][2] + xt[2][1] + xt[3][2]);

    val[0] = 0.0;
    val[1] = 0.0;
    val[2] = 0.0;
    /* Approximate integral*/    //printf("Approx int\t");
    for (j = 0; j < faces; j++) {
      if (fd[0] == j || fd[1] == j || fd[2] == j || fd[3] == j) { /*face j belongs to tetrahedra k*/
        /*Find vertex opposite of face j*/
        if (j == fd[0])
          num = 0;
        else if (j == fd[1])
          num = 1;
        else if (j == fd[2])
          num = 2;
        else
          num = 3;
        p[0] = x[v[num]][0];
        p[1] = x[v[num]][1];
        p[2] = x[v[num]][2];
        signum = compute_type_of_face_tet3d(dc->t3, k, num);

        if (is_dof[j] == 2) { /*Neumann face*/
          jj = idx2dof[j];
          alpha = g->v[jj] / (3.0 * vol);

          val[0] = val[0] + alpha * (s[0] - p[0]);
          val[1] = val[1] + alpha * (s[1] - p[1]);
          val[2] = val[2] + alpha * (s[2] - p[2]);
        }

        if (is_dof[j] == 0 || is_dof[j] == 1) { /*Inner or Dirichlet face*/
          jj = idx2dof[j];
          alpha = (signum * x1->v[jj]) / (3.0 * vol);

          val[0] = val[0] + alpha * (s[0] - p[0]);
          val[1] = val[1] + alpha * (s[1] - p[1]);
          val[2] = val[2] + alpha * (s[2] - p[2]);
        }

      }
    }
    // printf("assemble\n");
    q(s, fdata, v_flux); /*q exact*/
    val[0] = val[0] - v_flux->v[0];
    val[1] = val[1] - v_flux->v[1];
    val[2] = val[2] - v_flux->v[2];
    tnorm = val[0] * val[0] + val[1] * val[1] + val[2] * val[2];

    norm = norm + tnorm * vol;
  }

  del_avector(v_flux);
  return REAL_SQRT(norm);
}

real norml2_flux_facemidpoint_tet3drt0(pctet3drt0 dc,
    void (*q)(const real * x, void *fdata, pavector v), void *fdata,
    pcavector x1, pcavector g) {
  const real (*x)[3] = (const real (*)[3]) dc->t3->x;
  const uint (*t)[4] = (const uint (*)[4]) dc->t3->t;
  const uint *is_dof = dc->is_dof;
  const uint *idx2dof = dc->idx2dof;
  uint tetrahedra = dc->t3->tetrahedra;
  uint faces = dc->t3->faces;
  uint ndof = dc->ndof;
  uint nfix = dc->nfix;
  real vol, alpha, p[3], s[3];
  real val[3], norm, tnorm, value, signum;
  uint fd[4];
  uint j, k, l, ll, num;
  uint v[4];
  pavector v_flux;

  assert(x1 == 0 || x1->dim == ndof);
  assert(g == 0 || g->dim == nfix);

  norm = 0.0;
  tnorm = 0.0;
  v_flux = new_avector(3);

  for (k = 0; k < tetrahedra; k++) {
    /*volume of tetrahedra k*/    //printf("triangle %u\n", k);
    vol = getvolume_tetrahedra_tet3d(dc->t3, k);
    /*faces of tetrahedra k*/
    fd[0] = t[k][0];
    fd[1] = t[k][1];
    fd[2] = t[k][2];
    fd[3] = t[k][3];
    // printf("f: %u %u %u %u\n", fd[0], fd[1], fd[2], fd[3]);
    /*Get vertices*/
    getvertices_tet3d(dc->t3, k, v);

    tnorm = 0.0;
    for (j = 0; j <= 3; j++) {/*faces of tetrahedra k, for sum calculation*/
      /*centroid of face j in tetrahedra k*/
      getvertices_byface_tet3d(dc->t3, k, j, v);
      //  printf("face %u\n", j);
      // printf("v: %u %u %u %u\n", v[0], v[1], v[2], v[3]);
      alpha = 1.0 / 3.0;
      s[0] =
          alpha * (x[v[(j + 1) % 4]][0] + x[v[(j + 2) % 4]][0]
                   + x[v[(j + 3) % 4]][0]);
      s[1] =
          alpha * (x[v[(j + 1) % 4]][1] + x[v[(j + 2) % 4]][1]
                   + x[v[(j + 3) % 4]][1]);
      s[2] =
          alpha * (x[v[(j + 1) % 4]][2] + x[v[(j + 2) % 4]][2]
                   + x[v[(j + 3) % 4]][2]);
      //printf("s: %f %f %f\n", s[0], s[1], s[2]);
      val[0] = 0.0;
      val[1] = 0.0;
      val[2] = 0.0;
      for (l = 0; l <= faces; l++) { /*all faces*/
        if (fd[0] == l || fd[1] == l || fd[2] == l || fd[3] == l) { /*face l belongs to tetrahedra k*/
          getvertices_tet3d(dc->t3, k, v);
          // printf("v: %u %u %u %u\n", v[0], v[1], v[2], v[3]);
          /*Find vertex opposite of face j*/
          if (l == fd[0])
            num = 0;
          else if (l == fd[1])
            num = 1;
          else if (l == fd[2])
            num = 2;
          else
            num = 3;
          //printf("num = %u\n",num);
          p[0] = x[v[num]][0];
          p[1] = x[v[num]][1];
          p[2] = x[v[num]][2];
          /*vertex opposite of face j*/
          //p[0] = x[v[j]][0];
          //p[1] = x[v[j]][1];
          //p[2] = x[v[j]][2];
          signum = compute_type_of_face_tet3d(dc->t3, k, num);

          if (is_dof[l] == 2) { /*Neumann edge*/
            ll = idx2dof[l];
            alpha = g->v[ll] / (3.0 * vol);

            val[0] = val[0] + alpha * (s[0] - p[0]);
            val[1] = val[1] + alpha * (s[1] - p[1]);
            val[2] = val[2] + alpha * (s[2] - p[2]);
          }

          if (is_dof[l] == 0 || is_dof[l] == 1) { /*Inner or Dirichlet edge*/
            ll = idx2dof[l];
            alpha = (signum * x1->v[ll]) / (3.0 * vol);

            val[0] = val[0] + alpha * (s[0] - p[0]);
            val[1] = val[1] + alpha * (s[1] - p[1]);
            val[2] = val[2] + alpha * (s[2] - p[2]);
          }
        }
      }

      q(s, fdata, v_flux);
      val[0] = val[0] - v_flux->v[0];
      val[1] = val[1] - v_flux->v[1];
      val[2] = val[2] - v_flux->v[2];
      value = val[0] * val[0] + val[1] * val[1] + val[2] * val[2];
      tnorm = tnorm + value;
    }
    norm = norm + (vol / 4.0) * tnorm;
    //printf("tnorm = %f, norm = %f\n", tnorm, norm);
  }
  del_avector(v_flux);
  return REAL_SQRT(norm);
}

pclustergeometry build_tet3drt0_A_clustergeometry(pctet3drt0 rt,
    uint *idx) {

  pclustergeometry cg;
  uint i, j, k, l, m, c;
  uint nidx = rt->ndof;
  pctet3d t3 = rt->t3;
  uint faces = t3->faces;
  uint tetrahedra = t3->tetrahedra;
  const uint *is_dof = rt->is_dof;
  const uint *idx2dof = rt->idx2dof;
  const real (*x)[3] = (const real (*)[3]) rt->t3->x;
  const uint (*f)[3] = (const uint (*)[3]) rt->t3->f;
  const uint (*e)[2] = (const uint (*)[2]) rt->t3->e;
  cg = new_clustergeometry(3, nidx);

  uint e1, sv1, ev1, e2, sv2, ev2;
  real v1[3], v2[3], v3[3], mid[3];
  uint vt[4], v;

  /* Copying characteristic vertices into clustergeometry structure and
   *computing initial values of smin and smax */
  //c = 0;
  for (i = 0; i < faces; i++) {
    if (is_dof[i] == 0 || is_dof[i] == 1) { /*inner or Dirichlet edge*/
      c = idx2dof[i];
      /*Compute centroid of face (triangle) i*/
      e1 = f[i][0]; /*first edge of face i*/
      sv1 = e[e1][0];
      ev1 = e[e1][1]; /*start and end vertex of edge e1*/
      e2 = f[i][1]; /*second edge of face i*/
      sv2 = e[e2][0];
      ev2 = e[e2][1]; /*start and end vertex of edge e2*/
      if ((sv2 == sv1) || (sv2 == ev1)) {
        v3[0] = x[ev2][0];
        v3[1] = x[ev2][1];
        v3[2] = x[ev2][2];
      } else { /* (ev2 == sv1) || (ev2 == ev1)*/
        v3[0] = x[sv2][0];
        v3[1] = x[sv2][1];
        v3[2] = x[sv2][2];
      }
      v1[0] = x[sv1][0];
      v1[1] = x[sv1][1];
      v1[2] = x[sv1][2];
      v2[0] = x[ev1][0];
      v2[1] = x[ev1][1];
      v2[2] = x[ev1][2];

      mid[0] = (v1[0] + v2[0] + v3[0]) / 3.0;
      mid[1] = (v1[1] + v2[1] + v3[1]) / 3.0;
      mid[2] = (v1[2] + v2[2] + v3[2]) / 3.0;

      /*Write centroid of face i in cg*/
      cg->x[c][0] = mid[0];
      cg->x[c][1] = mid[1];
      cg->x[c][2] = mid[2];
      cg->smin[c][0] = mid[0];
      cg->smin[c][1] = mid[1];
      cg->smin[c][2] = mid[2];
      cg->smax[c][0] = mid[0];
      cg->smax[c][1] = mid[1];
      cg->smax[c][2] = mid[2];
      //c++;
      
      idx[c] = c;
    }
  }

  /* Updating values of smin and smax */
  for (i = 0; i < tetrahedra; i++) { /*all tetrahedra*/
    getvertices_tet3d(t3, i, vt); /*all vertices of tetrahedra i*/
    for (j = 0; j <= 3; j++) { /*all faces of i*/
      k = t3->t[i][j]; /*face k*/
      if (is_dof[k] == 0 || is_dof[k] == 1) { /*inner or Dirichlet face?*/
        m = idx2dof[k]; /*number of k in idx2dof*/
        for (l = 0; l <= 3; l++) {
          v = vt[l];
          if (cg->smin[m][0] > x[v][0])
            cg->smin[m][0] = x[v][0];
          if (cg->smax[m][0] < x[v][0])
            cg->smax[m][0] = x[v][0];
          if (cg->smin[m][1] > x[v][1])
            cg->smin[m][1] = x[v][1];
          if (cg->smax[m][1] < x[v][1])
            cg->smax[m][1] = x[v][1];
          if (cg->smin[m][2] > x[v][2])
            cg->smin[m][2] = x[v][2];
          if (cg->smax[m][2] < x[v][2])
            cg->smax[m][2] = x[v][2];
        }
      }
    }
  }

  update_point_bbox_clustergeometry(cg, nidx, idx);

  return cg;
}

pclustergeometry build_tet3drt0_B_clustergeometry(pctet3drt0 rt,
    uint *idx) {

  pclustergeometry cg;
  uint i, l, c;
  pctet3d t3 = rt->t3;
  uint tetrahedra = t3->tetrahedra;
  const real (*x)[3] = (const real (*)[3]) rt->t3->x;
  cg = new_clustergeometry(3, tetrahedra);
  real cen[3];
  uint v[4], w;

  /* Copying characteristic vertices into clustergeometry structure and
   *computing initial values of smin and smax */
  c = 0;
  for (i = 0; i < tetrahedra; i++) {
   // if(is_dof[i]){
     // c = idx2dof[i];
    getvertices_tet3d(t3, i, v);
    cen[0] = (x[v[0]][0] + x[v[1]][0] + x[v[2]][0] + x[v[3]][0]) / 4.0;
    cen[1] = (x[v[0]][1] + x[v[1]][1] + x[v[2]][1] + x[v[3]][1]) / 4.0;
    cen[2] = (x[v[0]][2] + x[v[1]][2] + x[v[2]][2] + x[v[3]][2]) / 4.0;
    cg->x[c][0] = cen[0];
    cg->x[c][1] = cen[1];
    cg->x[c][2] = cen[2];
    cg->smin[c][0] = cg->x[c][0];
    cg->smin[c][1] = cg->x[c][1];
    cg->smin[c][2] = cg->x[c][2];
    cg->smax[c][0] = cg->x[c][0];
    cg->smax[c][1] = cg->x[c][1];
    cg->smax[c][2] = cg->x[c][2];
    c++;
   // idx[c] = c;
  //}
  }
  /* Updating values of smin and smax */
  for (i = 0; i < tetrahedra; i++) { /*all tetrahedra*/
    getvertices_tet3d(t3, i, v); /*all vertices of tetrahedra i*/
    for (l = 0; l <= 3; l++) {/*All vertices of i*/
      w = v[l];
      if (cg->smin[i][0] > x[w][0])
        cg->smin[i][0] = x[w][0];
      if (cg->smax[i][0] < x[w][0])
        cg->smax[i][0] = x[w][0];
      if (cg->smin[i][1] > x[w][1])
        cg->smin[i][1] = x[w][1];
      if (cg->smax[i][1] < x[w][1])
        cg->smax[i][1] = x[w][1];
      if (cg->smin[i][2] > x[w][2])
        cg->smin[i][2] = x[w][2];
      if (cg->smax[i][2] < x[w][2])
        cg->smax[i][2] = x[w][2];
    }
  }

  update_point_bbox_clustergeometry(cg, tetrahedra, idx);
  return cg;
}

