
#include "tri2drt0.h"

#include <stdio.h>
#include <assert.h>
#include <math.h>

#include "basic.h"


ptri2drt0
new_tri2drt0(pctri2d t2)
{
  ptri2drt0 dc;
  
  uint *is_dof;
  uint *idx2dof;
  uint ndof, nfix;
  uint i;
  
  dc      = (ptri2drt0) allocmem(sizeof(tri2drt0));
  is_dof = dc->is_dof = (uint *) allocmem((size_t) sizeof(uint) * t2->edges);
  idx2dof = dc->idx2dof = (uint *) allocuint(t2->edges);
  dc->t2  = t2;
 
  for(i=0;i<t2->edges;i++){
    if(t2->eb[i] == 0) is_dof[i] = 0;
    else is_dof[i] = 1;
  } 
  
  ndof = 0;
  nfix = 0;
  for(i=0; i<t2->edges; i++){
    if(is_dof[i] == 0)
      idx2dof[i] = ndof++;
    else
      idx2dof[i] = nfix++;
  }

  dc->ndof = ndof;
  dc->nfix = nfix;

  return dc;  
}


void update_tri2drt0(ptri2drt0 rt0){ 
  
 pctri2d t2 = rt0->t2;
 uint *is_dof = rt0->is_dof;
 uint *idx2dof = rt0->idx2dof;
 uint i, nfix, ndof;
 
 nfix = 0; ndof = 0;
 for(i=0; i<t2->edges;i++){
   if(is_dof[i] == 0 || is_dof[i] == 1) {
     idx2dof[i] = ndof++;// printf("0/1 \t");
   }
   //else if (is_dof[i] == 1){
   //  idx2dof[i] = ndof++; printf("1 \t");
   //}
   else{
     assert(is_dof[i] == 2);
     idx2dof[i] = nfix++;//printf("2 \t");
   }
  }
  rt0->ndof = ndof;
  rt0->nfix = nfix;
}

void
del_tri2drt0(ptri2drt0 dc)
{
  freemem(dc->is_dof);
  freemem(dc->idx2dof);
  freemem(dc);
}


psparsematrix
build_tri2drt0_A_sparsematrix(pctri2drt0 dc)
{
  pctri2d t2 = dc->t2;
  const uint *is_dof = dc->is_dof;
  const uint *idx2dof = dc->idx2dof;
  uint ndof = dc->ndof;
  uint triangles  = dc->t2->triangles;
  psparsepattern sp;
  psparsematrix A;
  uint i, j, ii, jj, d, et[3];
  
  sp = new_sparsepattern(ndof, ndof);

  for(d=0; d<triangles; d++){
    /* Get edges */
    et[0] = t2->t[d][0];
    et[1] = t2->t[d][1];
    et[2] = t2->t[d][2];
    /* Build sparsepattern */
    for(i=0; i<3; i++){
      if((is_dof[et[i]] == 0) || (is_dof[et[i]]==1)){ /*inner edge or Dirichlet edge*/
	ii = idx2dof[et[i]];
	for(j=0; j<3; j++){
	  if((is_dof[et[j]] == 0) || (is_dof[et[j]] == 1)){ /*inner or Dirichlet*/
	    jj = idx2dof[et[j]];
	    addnz_sparsepattern(sp,ii,jj);
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
build_tri2drt0_A_interaction_sparsematrix(pctri2drt0 dc)
{
  pctri2d   t2 = dc->t2;
  const uint *is_dof = dc->is_dof;
  const uint *idx2dof = dc->idx2dof;
  uint      ndof = dc->ndof;
  uint      nfix = dc->nfix;
  psparsepattern sp;
  psparsematrix Af;
  uint      i, j, ii, jj, t, et[3];

  sp = new_sparsepattern(ndof, nfix);

  for (t = 0; t < t2->triangles; t++) {
    /* Get edges */
    et[0] = t2->t[t][0];
    et[1] = t2->t[t][1];
    et[2] = t2->t[t][2];
    /* Build sparsepattern */
    for (i = 0; i < 3; i++) {
      if ((is_dof[et[i]] == 0) || (is_dof[et[i]] == 1)) { /*inner or Dirichlet*/
	ii = idx2dof[et[i]];
	assert(ii < ndof);

	for (j = 0; j < 3; j++) {
	  if (is_dof[et[j]] == 2) { /*Neumann*/
	    jj = idx2dof[et[j]];
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
build_tri2drt0_B_sparsematrix(pctri2drt0 dc)
{
  pctri2d t2 = dc->t2;
  const uint *is_dof = dc->is_dof;
  const uint *idx2dof = dc->idx2dof;
  uint ndof = dc->ndof;
  uint triangles  = dc->t2->triangles;
  psparsepattern sp;
  psparsematrix A;
  uint i, ii, d, et[3], nt;
  
  nt = dc->t2->triangles;
  sp = new_sparsepattern(nt, ndof);

  for(d=0; d<triangles; d++){
    /* Get edges */
    et[0] = t2->t[d][0];
    et[1] = t2->t[d][1];
    et[2] = t2->t[d][2];
    /* Build sparsepattern */
    for(i=0; i<3; i++){
      if((is_dof[et[i]] == 0) || (is_dof[et[i]] == 1)){ /*Dirichlet or inner*/
	ii = idx2dof[et[i]];
	addnz_sparsepattern(sp,d,ii);
      }
    }
  }
  
  /* Build sparsematrix from sparsepattern */
  A = new_zero_sparsematrix(sp);
  
  del_sparsepattern(sp);
  
  return A;
}

psparsematrix
build_tri2drt0_B_interaction_sparsematrix(pctri2drt0 dc)
{
  pctri2d t2 = dc->t2;
  const uint *is_dof = dc->is_dof;
  const uint *idx2dof = dc->idx2dof;
 // uint ndof = dc->ndof;
  uint nfix = dc->nfix;
  uint triangles  = dc->t2->triangles;
  psparsepattern sp;
  psparsematrix A;
  uint i, ii, d, et[3], nt;
  
  nt = dc->t2->triangles;
  sp = new_sparsepattern(nt, nfix);

  for(d=0; d<triangles; d++){
    /* Get edges */
    et[0] = t2->t[d][0];
    et[1] = t2->t[d][1];
    et[2] = t2->t[d][2];
    /* Build sparsepattern */
    for(i=0; i<3; i++){
      if(is_dof[et[i]]==2){
	ii = idx2dof[et[i]];
	addnz_sparsepattern(sp,d,ii);
      }
    }
  }
  
  /* Build sparsematrix from sparsepattern */
  A = new_zero_sparsematrix(sp);
  
  del_sparsepattern(sp);
  
  return A;
}
#if  0
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
      v[0] = xt[0];
  
    if( (e[et[1]][0] == xt[0] && e[et[1]][1] == xt[1]) || (e[et[1]][1] == xt[0] && e[et[1]][0] == xt[1]) )
    v[1] = xt[2];
    else if( (e[et[1]][0] == xt[0] && e[et[1]][1] == xt[2]) || (e[et[1]][0] == xt[2] && e[et[1]][1] == xt[0]) )
      v[1]=xt[1];
    else 
      v[1] = xt[0];
    
    if( (e[et[2]][0] == xt[0] && e[et[2]][1] == xt[1]) || (e[et[2]][1] == xt[0] && e[et[2]][0] == xt[1]) )
    v[2] = xt[2];
    else if( (e[et[2]][0] == xt[0] && e[et[2]][1] == xt[2]) || (e[et[2]][0] == xt[2] && e[et[2]][1] == xt[0]) )
      v[2]=xt[1];
    else 
      v[2] = xt[0];
    
    
}
#endif

static real
scalar(real x[2], real y[2])
{
  real s;
  s = x[0] * y[0] + x[1] * y[1];
  
 return s; 
}

/*if clockwise normal is outer normal at edge i in triangle t: return  1
                                                               else   -1 */
real
compute_type_of_edge_tri2d(pctri2d t2, uint d, uint i){
  
  const uint(*t)[3] = (const uint(*)[3]) t2->t;
  const uint(*e)[2] = (const uint(*)[2]) t2->e;
  const real(*x)[2] = (const real(*)[2]) t2->x;
  real signum, nr[2], v[2];
  uint ei, ej;
  
  /*Find global numbers of edges*/
  ei = t[d][i];  
  ej = t[d][(i+1)%3];
  /*Compute clockwise normal of edge ei*/
  nr[0] = x[e[ei][1]][1] - x[e[ei][0]][1];
  nr[1] = x[e[ei][0]][0] - x[e[ei][1]][0];

  if(e[ej][0] == e[ei][0] || e[ej][0] == e[ei][1]){ /*startvertex of ej belongs to ei*/
    v[0] = x[e[ej][1]][0] - x[e[ej][0]][0];
    v[1] = x[e[ej][1]][1] - x[e[ej][0]][1];
  }
  else{ /*endvertex of ej belongs to ei*/
    v[0] = x[e[ej][0]][0] - x[e[ej][1]][0];
    v[1] = x[e[ej][0]][1] - x[e[ej][1]][1];
  }
  /*Return -1, if innerproduct of v and nr is > 0, else +1*/  
  if(nr[0] * v[0] + nr[1] * v[1] > 0) /*nr is inner normal at ei and triangle d*/
    signum = -1.0;
  else /*nr is outer normal at ei and triangle d*/
    signum = 1.0;
    
  return signum;
}


void
assemble_tri2drt0_darcy_A_sparsematrix(pctri2drt0 dc, psparsematrix A, psparsematrix Af, pavector K)
{
  const real (*x)[2] = (const real (*)[2]) dc->t2->x;
 // const uint (*e)[2] = (const uint (*)[2]) dc->t2->e;
  const uint (*t)[3] = (const uint (*)[3]) dc->t2->t;
  const uint *is_dof = dc->is_dof;
  const uint *idx2dof = dc->idx2dof;
  uint triangles  = dc->t2->triangles;
  
  uint ndof = dc->ndof;
  uint nfix = dc->nfix;
  uint i, j, ii, jj, d, k;
  real At[3][3];
  uint et[3],v[3];
  real T, signum_i, signum_j;
  real x_ij[2], x_jk[2], x_ki[2];

  //triangles = 1;
  
  for(d=0; d<triangles; d++){
     /*Get edges of triangle d*/
    et[0] = t[d][0];
    et[1] = t[d][1];
    et[2] = t[d][2];
    /*Get vertices of triangle d*/
    getvertices_tri2d(dc->t2, d, v); //printf("vertices %u %u %u \n", v[0], v[1], v[2]);
    /*get area of T*/
    T = (1.0/2)* fabs( ((x[v[2]][0]-x[v[0]][0])*(x[v[1]][1]-x[v[0]][1]))
		      -((x[v[2]][1]-x[v[0]][1])*(x[v[1]][0]-x[v[0]][0])) ); 
    /*Initialise At with zero*/
    for(i=0;i<3;i++)
      for(j=0;j<3;j++)
	At[i][j] = 0.0;
    /* Compute element matrix */
    for(i=0;i<3;i++){
      for(j=0;j<3;j++){
	/*non diagonal entry*/
	if (i!=j){
	   if((i==0 && j==1) || (i==1 && j==0)) k=2;
	   else if((i==0 && j==2) || (i==2 && j==0)) k=1;
           else k=0;
	   x_ij[0] = x[v[i]][0] - x[v[j]][0]; x_ij[1] = x[v[i]][1] - x[v[j]][1]; 
	   x_jk[0] = x[v[j]][0] - x[v[k]][0]; x_jk[1] = x[v[j]][1] - x[v[k]][1]; 
	   x_ki[0] = x[v[k]][0] - x[v[i]][0]; x_ki[1] = x[v[k]][1] - x[v[i]][1];
	   At[i][j] = - scalar(x_ij,x_ij) + scalar(x_jk,x_ij) 
	              - 2 * scalar(x_jk,x_ki) + scalar(x_ij,x_ki);
           //printf("At[%u][%u] = %f \n", i, j, At[i][j]);
	   signum_i = compute_type_of_edge_tri2d(dc->t2, d, i);
	   signum_j = compute_type_of_edge_tri2d(dc->t2, d, j);
	   At[i][j] = (signum_i * signum_j * At[i][j] * K->v[d]) / (48 * T); /*Faktor k fehlt noch*/
        }
	else {/*i=j*/ /*j is the vertex opposite of edge j*/
	  i = ((j+1)%3);
	  k = ((j+2)%3);
	 // printf("i = %u, j = %u, k = %u\n", i, j, k);
	 // printf("v[i] = %u, v[j] = %u, v[k] = %u\n", v[i], v[j], v[k]);
	  x_ij[0] = x[v[i]][0] - x[v[j]][0]; x_ij[1] = x[v[i]][1] - x[v[j]][1]; 
	  x_jk[0] = x[v[j]][0] - x[v[k]][0]; x_jk[1] = x[v[j]][1] - x[v[k]][1]; 
	  x_ki[0] = x[v[k]][0] - x[v[i]][0]; x_ki[1] = x[v[k]][1] - x[v[i]][1];
	  At[j][j] = scalar(x_ij, x_ij) + scalar(x_jk,x_jk) - scalar(x_jk,x_ij);
	  //printf("At[%u][%u] = %f \n", j, j, At[j][j]);
	  At[j][j] = (At[j][j] *K->v[d])/ (24*T); /*Faktor k fehlt noch*/
	  i = j;
        }
      } 
    }

   // printf("At = \n");
   // for(i=0;i<3;i++){
   //   for(j=0;j<3;j++)
   //     printf("At[%u][%u] = %f \t",i, j,  At[i][j]);
   //  printf("\n"); 
  //  }
      
    /* Add to system matrix */
    for(i=0;i<3;i++){
      if(is_dof[et[i]] == 0 || is_dof[et[i]] == 1){
	ii = idx2dof[et[i]];
	for(j=0;j<3;j++){
	  if(is_dof[et[j]] == 0 || is_dof[et[j]] == 1){
	  jj = idx2dof[et[j]];
	// printf("ii = %u jj = %u\n", ii, jj);
	  addentry_sparsematrix(A, ii, jj, At[i][j]);
	 // print_sparsematrix(A);
	  }
        }
      }
    }
    /* Add to interaction matrix */
    if(Af){
      for(i=0; i<3; i++){
	if(is_dof[et[i]] == 0 || is_dof[et[i]] == 1){
	  ii = idx2dof[et[i]];
	  assert(ii < ndof);
	   for(j=0; j<3; j++){
            if(is_dof[et[j]] == 2){
	      jj = idx2dof[et[j]];
	      assert(jj < nfix);
	      addentry_sparsematrix(Af, ii, jj, At[i][j]);
	    }
	  }
        }
      }
    }
    
    //printf("triangle %u\n", d);
    //printf("with edges %u %u %u\n", et[0], et[1], et[2]);
   // print_sparsematrix(A);
   // print_sparsematrix(Af);	

  }

}


void
assemble_tri2drt0_darcy_B_sparsematrix(pctri2drt0 dc, psparsematrix A, psparsematrix Af)
{
  //const real (*x)[2] = (const real (*)[2]) dc->t2->x;
  //const uint (*e)[2] = (const uint (*)[2]) dc->t2->e;
  const uint (*t)[3] = (const uint (*)[3]) dc->t2->t;
  const uint *is_dof = dc->is_dof;
  const uint *idx2dof = dc->idx2dof;
  uint triangles  = dc->t2->triangles;
  uint i, ii, d; //E_i;
  real signum, At;
  uint et[3];
 // int eti_start, eti_end;
 // real x1[2], x2[2];
  //real a,b, sigma_i;
  
  for(d=0; d<triangles; d++){
    /*Get edges of triangle d*/
    et[0] = t[d][0];
    et[1] = t[d][1];
    et[2] = t[d][2];
    /* Compute part of matrices computed by edges of triangle d*/
    for(i=0;i<3;i++){
      /*Find direction of normal at edge i*/
      signum = compute_type_of_edge_tri2d(dc->t2, d, i);
      At = signum * -1;
      /*Add to system matrix*/
      if(is_dof[et[i]] == 0 || is_dof[et[i]] == 1){ /*inner or Dirichlet edge?*/
        ii = idx2dof[et[i]];
	addentry_sparsematrix(A, d, ii, At);
      }
      /*Add to interaction matrix*/
      if(Af){
        if(is_dof[et[i]] == 2){ /*Neumann edge?*/
	  ii = idx2dof[et[i]];
	  addentry_sparsematrix(Af, d, ii, At);
        }
      }  
    }
   // printf("triangle %d\n", d);
    //print_sparsematrix(A);
  }
}


void
assemble_tri2drt0_b_D_avector(pctri2drt0 dc,
		  field (*d)(const real *e, void *fdata), void *fdata,
		  pavector dv)
{
  const uint (*e)[2] = (const uint(*)[2]) dc->t2->e;
  const real (*x)[2] = (const real(*)[2]) dc->t2->x;
  const uint *is_dof = dc->is_dof;
  const uint *idx2dof = dc->idx2dof;
  uint edges = dc->t2->edges;
  //uint nfix = dc->nfix;
  uint i, ii, sv, ev;
  real start[2], end[2], m[2];

 // assert(dv->dim == nfix);

  clear_avector(dv);
  for(i=0; i<edges; i++){
    if(is_dof[i] == 1) {
      ii = idx2dof[i];
      //assert(ii < nfix);
      sv = e[i][0]; ev = e[i][1];
      start[0] = x[sv][0]; start[1] = x[sv][1];
      end[0]   = x[ev][0]; end[1]   = x[ev][1];
      m[0] = (start[0] + end[0]) / 2.0;
      m[1] = (start[1] + end[1]) / 2.0;
      //printf("edge %u\n", i);
      //printf("m[0] = %f, m[1] = %f\n", m[0], m[1]);
      dv->v[ii] = (d ? - d(m,fdata) : 0.0);
    }
    
  }
}


void
assemble_tri2drt0_b_f_avector(ptri2drt0 dc,
                   field (*f)(const real *x, void *fdata), void *fdata,
                   pavector fv)
{
  pctri2d t2 = dc->t2;
  const real (*x)[2] = (const real (*)[2]) t2->x;
  uint triangles = dc->t2->triangles;
 // uint ndof = dc->ndof;
  real xt[3][2], xd[2][2], xm[2], vol;
  field f01, f02, f12;
  uint i, d, v[3];
  
  //assert(fv->dim == ndof);
  
  clear_avector(fv);
  
  for(d=0; d<triangles; d++)
  {
    /* Get vertices of triangle d*/
    getvertices_tri2d(t2, d, v);
    
    /* Get vertex coordinates */
    for(i=0; i<3;i++)
    {
      xt[i][0] = x[v[i]][0];
      xt[i][1] = x[v[i]][1];
    }
    
    /* Compute difference vectors */
    for(i=0; i<2; i++)
         
    {
      xd[i][0] = xt[i+1][0] - xt[0][0];
      xd[i][1] = xt[i+1][1] - xt[0][1];
    }
    
    /* Compute determinant */
    vol = fabs(xd[0][0]*xd[1][1] - xd[0][1]*xd[1][0]);
    
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
    
    fv->v[d] = - (vol/2.0) * (f01 + f02 +f12) * (1.0/3.0);
    }
  
}

void assemble_tri2drt0_g_N_avector(ptri2drt0 dc, 
				   field (*f)(const uint *e, void *data), void *data,
				   pavector g)
{
  
  const uint (*e)[2] = (const uint(*)[2]) dc->t2->e;
  const real (*x)[2] = (const real(*)[2]) dc->t2->x;
  const uint *idx2dof = dc->idx2dof;
  uint edges = dc->t2->edges;
  const uint *is_dof = dc->is_dof;
  uint i, sv, ev, ii;
  real start[2], end[2], z, y, e_i;
 
  for(i=0;i<edges;i++){
    if(is_dof[i] == 2){ /*Neumann edges*/
     ii = idx2dof[i];
      /*leghth of edges*/
     sv = e[i][0]; ev = e[i][1]; 
     start[0] = x[sv][0]; start[1] = x[sv][1];
     end[0]   = x[ev][0]; end[1]   = x[ev][1];
     z = end[0] - start[0];
     y = end[1] - start[1];
     e_i = sqrt(z*z + y*y);
     // m[0] = (start[0] + end[0]) / 2.0;
     // m[1] = (start[1] + end[1]) / 2.0;
     //printf("e_i = %f\n", e_i);
     g->v[ii] = (f ? (e_i * f(e[i], data)) : 0.0);
    }
  }
}


real
getarea_triangle_tri2d(pctri2d t2, uint t){
  
  uint v[3], i;
  real area, xt[3][2];
  const real(*x)[2] = (const real(*)[2]) t2->x;
  
  /*Get vertices*/
  getvertices_tri2d(t2, t, v);
  /* Get vertex coordinates */
  for (i = 0; i < 3; i++) {
    xt[i][0] = x[v[i]][0];
    xt[i][1] = x[v[i]][1];
  }
  /* Compute area */
  area = fabs((xt[1][0] - xt[0][0]) * (xt[2][1] - xt[0][1]) -
              (xt[1][1] - xt[0][1]) * (xt[2][0] - xt[0][0])) * 0.5;
	      
  return area;  
}

/*Schwerpunktsregel*/
real
norml2_pressure_centroid_tri2drt0(pctri2drt0 dc,
	       field(*p) (const real * x, void *fdata), void *fdata,
	       pcavector x2)
{
  const     real(*x)[2] = (const real(*)[2]) dc->t2->x;
  uint      triangles = dc->t2->triangles;
  real      s[2], xt[3][2];
  real      val, norm, area, alpha;
  uint      v[3];
  uint      i, k;
  
  norm = 0.0;
  
  for (k = 0; k < triangles; k++) {

    /*Get vertices*/
    getvertices_tri2d(dc->t2, k, v);
    /* Get vertex coordinates */
    for (i = 0; i < 3; i++) {
      xt[i][0] = x[v[i]][0];
      xt[i][1] = x[v[i]][1];
    }
    /*area of triangle k*/		
    area = getarea_triangle_tri2d(dc->t2, k);
    /* Approximate integral */
    alpha = 1.0/3.0;
    s[0] = alpha * (xt[0][0] + xt[1][0] + xt[2][0]);
    s[1] = alpha * (xt[0][1] + xt[1][1] + xt[2][1]);
    val = p(s, fdata) - x2->v[k];
    //printf("p ? %f x2 = %f val = %f\n", (double) p(s, fdata), (double)x2->v[k], val);
    val = val * val;
  
    norm = norm + area * val;
   
  }

  return REAL_SQRT(norm);
}

/*Kantenmittelpunkt*/
real
norml2_pressure_edgemidpoint_tri2drt0(pctri2drt0 dc,
	       field(*p) (const real * x, void *fdata), void *fdata,
	       pcavector x2)
{
  const     real(*x)[2] = (const real(*)[2]) dc->t2->x;
  const     uint(*e)[2] = (const uint(*)[2]) dc->t2->e;
  uint      triangles = dc->t2->triangles;
  const     uint(*t)[3] = (const uint(*)[3]) dc->t2->t;
  //const     uint *is_dof = dc->is_dof;
  real       s[2];
  real      val, norm, area, alpha;
  uint      ed[3];
  uint       k, edge, start, end;
//  uint l;
  real tnorm;
  
  norm = 0.0;
  
  for (k = 0; k < triangles; k++) {

    /*area of triangle k*/		
    area = getarea_triangle_tri2d(dc->t2, k);
    
    /*Get edges of triangle k*/
    ed[0] = t[k][0];
    ed[1] = t[k][1];
    ed[2] = t[k][2];
    tnorm = 0.0;
    /*Approximate integral*/
    for(edge=0;edge<=2;edge++){
      start = e[ed[edge]][0];
      end   = e[ed[edge]][1];
      /*Compute midpoint of edge in T_k*/
      alpha = 1.0/2.0;
      s[0] = alpha * (x[start][0] + x[end][0]);
      s[1] = alpha * (x[start][1] + x[end][1]);
      val = p(s,fdata) - x2->v[k];//printf("p ? %f x2 = %f val = %f\n", (double) p(s, fdata), (double)x2->v[k], val);
#if 0 
     if(is_dof[ed[edge]] == 0){/*inner edge*/
      for(l=0;l<triangles;l++){
          if((ed[edge] == t[l][0] ||ed[edge] == t[l][1] || ed[edge] == t[l][2]) && l!=k){
	    assert(k!=l);
	    val = val - x2->v[l];
         }
        } 
     }
#endif     
      tnorm = tnorm + val * val;
    }
    norm = norm + (area / 3.0) * tnorm;
  }
  return REAL_SQRT(norm);
}



real
norml2_flux_centroid_tri2drt0(pctri2drt0 dc,
	       void (*q) (const real * x, void *fdata, pavector v), void *fdata,
	       pcavector x1, pcavector g)
{
  const     real(*x)[2] = (const real(*)[2]) dc->t2->x;
  const     uint(*t)[3] = (const uint(*)[3]) dc->t2->t;
  const uint *is_dof = dc->is_dof;
  const uint *idx2dof = dc->idx2dof;
  uint      triangles = dc->t2->triangles;
  uint 	    edges = dc->t2->edges;
  uint      ndof = dc->ndof;
  uint      nfix = dc->nfix;
  real      xt[3][2], area, s[2], alpha, p[2];
  real      val[2], norm, tnorm, signum;
  uint      ed[3];
  uint      i, j, k, jj, num;
  uint v[3];
  pavector v_flux;

  assert(x1 == 0 || x1->dim == ndof);
  assert(g == 0 || g->dim == nfix);

  norm = 0.0;
  v_flux = new_avector(2);
  signum = 0.0;

  for (k = 0; k < triangles; k++) {
    /*Area of triangle k*///printf("triangle %u\t", k);
    area = getarea_triangle_tri2d(dc->t2, k);
    /*edges of triangle k*/
    ed[0] = t[k][0];
    ed[1] = t[k][1];
    ed[2] = t[k][2];
    /*Get vertices*/
    getvertices_tri2d(dc->t2, k, v);
    /* Get vertex coordinates */
    for(i = 0; i < 3; i++) {
      xt[i][0] = x[v[i]][0];
      xt[i][1] = x[v[i]][1];
    }
    /*Compute centroid*/
    alpha = 1.0/3.0;
    s[0]  = alpha * (xt[0][0] + xt[1][0] + xt[2][0]);
    s[1]  = alpha * (xt[0][1] + xt[1][1] + xt[2][1]);
    
    val[0] = 0.0; val[1] = 0.0;
    /* Approximate integral*///printf("Approx int\t");
    for(j=0;j<edges;j++){
      if(ed[0] == j || ed[1] == j || ed[2] == j){ /*edge j belongs to triangle k*/
        /*Find vertex opposite of edge j*/
	 if(j==ed[0]) num = 0;
	 else if(j == ed[1]) num = 1;
	 else num = 2;
	 p[0] = x[v[num]][0];
	 p[1] = x[v[num]][1];
	 signum = compute_type_of_edge_tri2d(dc->t2, k, num); //printf("signum = %f \n", signum);
	 
	if(is_dof[j] == 2){ /*Neumann edge*/
	  jj =  idx2dof[j];
	  alpha = g->v[jj] / (2.0 * area);
	  
	  val[0] = val[0] + alpha * (s[0] - p[0]);
	  val[1] = val[1] + alpha * (s[1] - p[1]);
	}
	
      	if(is_dof[j] == 0 || is_dof[j] == 1){ /*Inner or Dirichlet edge*/
	  jj =  idx2dof[j];
	  alpha = (signum * x1->v[jj]) / (2.0 * area);
	  
	  val[0] = val[0] + alpha * (s[0] - p[0]);
	  val[1] = val[1] + alpha * (s[1] - p[1]);
	}
	
      }
    }
   // printf("assemble\n");
    q(s, fdata, v_flux);
    val[0] = val[0] - v_flux->v[0];
    val[1] = val[1] - v_flux->v[1];
    tnorm = val[0] * val[0] + val[1] * val[1];
    

    norm = norm + tnorm * area;
  }

 del_avector(v_flux);
  return REAL_SQRT(norm);
}

real
norml2_flux_edgemidpoint_tri2drt0(pctri2drt0 dc,
	       void (*q) (const real * x, void *fdata, pavector v), void *fdata,
	       pcavector x1, pcavector g)
{
  const     real(*x)[2] = (const real(*)[2]) dc->t2->x;
  const     uint(*e)[2] = (const uint(*)[2]) dc->t2->e;
  const     uint(*t)[3] = (const uint(*)[3]) dc->t2->t;
  const uint *is_dof = dc->is_dof;
  const uint *idx2dof = dc->idx2dof;
  uint      triangles = dc->t2->triangles;
  uint 	    edges = dc->t2->edges;
  uint      ndof = dc->ndof;
  uint      nfix = dc->nfix;
  real      area, alpha, p[2], m[2];
  real      val[2], norm, tnorm, value, signum;
  uint      ed[3];
  uint      j, k, num, start, end, l, ll;
  uint v[3];
  pavector v_flux;

  assert(x1 == 0 || x1->dim == ndof);
  assert(g == 0 || g->dim == nfix);

  norm = 0.0; 
  tnorm = 0.0;
  v_flux = new_avector(2);

  for (k = 0; k < triangles; k++) {
    /*Area of triangle k*///printf("triangle %u\t", k);
    area = getarea_triangle_tri2d(dc->t2, k);
    /*edges of triangle k*/
    ed[0] = t[k][0];
    ed[1] = t[k][1];
    ed[2] = t[k][2];
    /*Get vertices*/
    getvertices_tri2d(dc->t2, k, v);
    
    tnorm = 0.0;
    for(j=0;j<=2;j++){/*edges of triangle k, for sum calculation*/
      /*midpoint of edge j in triangle k*/
      start = e[ed[j]][0];
      end   = e[ed[j]][1];
      alpha = 1.0/2.0;
      m[0] = alpha * (x[start][0] + x[end][0]);
      m[1] = alpha * (x[start][1] + x[end][1]);
      val[0] = 0.0; val[1] = 0.0;
      for(l=0;l<=edges;l++){ /*all edges*/
	if(ed[0] == l || ed[1] == l || ed[2] == l){ /*edge l belongs to triangle k*/
          /*Find vertex opposite of edge l*/
	  if(l==ed[0]) num = 0;
	  else if(l == ed[1]) num = 1;
	  else num = 2;
	  p[0] = x[v[num]][0];
	  p[1] = x[v[num]][1];
	  signum = compute_type_of_edge_tri2d(dc->t2, k, num);
	 
	 if(is_dof[l] == 2){ /*Neumann edge*/
	   ll =  idx2dof[l];
	   alpha = g->v[ll] / (2.0 * area);
	  
	   val[0] = val[0] + alpha * (m[0] - p[0]);
	   val[1] = val[1] + alpha * (m[1] - p[1]);
	 }
	
      	 if(is_dof[l] == 0 || is_dof[l] == 1){ /*Inner or Dirichlet edge*/
	   ll =  idx2dof[l];
	   alpha = (signum * x1->v[ll]) / (2.0 * area);
	  
	   val[0] = val[0] + alpha * (m[0] - p[0]);
	   val[1] = val[1] + alpha * (m[1] - p[1]);
	 }
        }
      }
      
      q(m, fdata, v_flux);
      val[0] = val[0] - v_flux->v[0];
      val[1] = val[1] - v_flux->v[1];
      value = val[0] * val[0] + val[1] * val[1];
      tnorm = tnorm + value;
    }
    norm = norm + (area/3.0) * tnorm;  
    //printf("tnorm = %f, norm = %f\n", tnorm, norm);
  }  
  del_avector(v_flux);
  return REAL_SQRT(norm);
}

pclustergeometry build_tri2drt0_A_clustergeometry(pctri2drt0 rt, uint *idx)
{
  
  pclustergeometry cg;
  uint i, j, k, l, m,c;
  uint nidx = rt->ndof;
  pctri2d t2 = rt->t2;
  uint edges = t2->edges;
  uint triangles = t2->triangles;
  const uint *is_dof = rt->is_dof;
  const uint *idx2dof = rt->idx2dof;
  const real (*x)[2] = (const real (*)[2])rt->t2->x;
  const uint (*e)[2] = (const uint (*)[2])rt->t2->e;
  cg = new_clustergeometry(2, nidx);
  uint sv, ev, v;
  real start[2], end[2], mid[2];
  uint vt[3];
	
  /* Copying characteristic vertices into clustergeometry structure and 
   *computing initial values of smin and smax */
  c = 0;
  for(i=0; i<edges; i++){
    if(is_dof[i] == 0 || is_dof[i] == 1){ /*inner or Dirichlet edge*/
      sv = e[i][0]; ev = e[i][1];
      start[0] = x[sv][0]; start[1] = x[sv][1];
        end[0] = x[ev][0];   end[1] = x[ev][1];
      mid[0] = (start[0] + end[0]) / 2.0;
      mid[1] = (start[1] + end[1]) / 2.0;
      cg->x[c][0] = mid[0];
      cg->x[c][1] = mid[1];
      cg->smin[c][0] = mid[0];
      cg->smin[c][1] = mid[1];
      cg->smax[c][0] = mid[0];
      cg->smax[c][1] = mid[1];
      c++;
    }
  }
  
  /* Updating values of smin and smax */
  for(i=0; i<triangles; i++){ /*all triangles*/
    getvertices_tri2d(t2, i, vt); /*all vertices of triangle i*/
    for(j=0;j<=2;j++){ /*all edges of triangle i*/
      k = t2->t[i][j]; /*edge k*/
      if(is_dof[k]== 0 || is_dof[k] == 1){
        m = idx2dof[k]; /*number of k in idx2dof*/
        for(l = 0;l<=2;l++){/*all vertices of triangle i*/
	  v = vt[l];
          if(cg->smin[m][0] > x[v][0])
	    cg->smin[m][0] = x[v][0];
	  if(cg->smax[m][0] < x[v][0])
	    cg->smax[m][0] = x[v][0];
	  if(cg->smin[m][1] > x[v][1])
	    cg->smin[m][1] = x[v][1];
	  if(cg->smax[m][1] < x[v][1])
	    cg->smax[m][1] = x[v][1];
         }
      }
    }
  }
  
  update_point_bbox_clustergeometry(cg, nidx, idx);

  return cg;
}

pclustergeometry build_tri2drt0_B_clustergeometry(pctri2drt0 rt, uint *idx)
{
  
  pclustergeometry cg;
  uint i,l,c; 
  pctri2d t2 = rt->t2;
  uint triangles = t2->triangles;
  const real (*x)[2] = (const real (*)[2])rt->t2->x;
  cg = new_clustergeometry(2, triangles);
  real cen[2];
  uint v[3], w;
	
  /* Copying characteristic vertices into clustergeometry structure and 
   *computing initial values of smin and smax */
  c = 0;
  for(i=0; i<triangles; i++){
      getvertices_tri2d(t2, i, v);
      cen[0] = (x[v[0]][0] + x[v[1]][0] + x[v[2]][0]) / 3.0;
      cen[1] = (x[v[0]][1] + x[v[1]][1] + x[v[2]][1]) / 3.0;
      cg->x[c][0] = cen[0];
      cg->x[c][1] = cen[1];
      cg->smin[c][0] = cg->x[c][0];
      cg->smin[c][1] = cg->x[c][1];
      cg->smax[c][0] = cg->x[c][0];
      cg->smax[c][1] = cg->x[c][1]; 
      c++;
  }
 
  /* Updating values of smin and smax */
  for(i=0; i<triangles; i++){ /*all triangles*/
    getvertices_tri2d(t2, i, v); /*all vertices of triangle i*/
      for(l = 0;l<=2;l++){/*All vertices*/
	w = v[l];
        if(cg->smin[i][0] > x[w][0])
	  cg->smin[i][0] = x[w][0];
        if(cg->smax[i][0] < x[w][0])
         cg->smax[i][0] = x[w][0];
	if(cg->smin[i][1] > x[w][1])
	  cg->smin[i][1] = x[w][1];
	if(cg->smax[i][1] < x[w][1])
	  cg->smax[i][1] = x[w][1];
      }
   }

  update_point_bbox_clustergeometry(cg, triangles, idx);

  return cg;
}

