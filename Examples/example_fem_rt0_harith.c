/*----------------------------------------------------------------------------*/
/* Working with lowest order Raviart-Thomas elements and hierarchical matrices*/ 
/* for Darcy's Law                                                            */
/*----------------------------------------------------------------------------*/ 
#include <stdio.h>

#include "parameters.h"	/* Read parameters interactively */

#include "tri2d.h"	/* 2-dimensional mesh*/
#include "tri2drt0.h"	/* discretisation with Raviart-Thomas elements based on tri2d.h*/

#include "tet3d.h"	/* 3-dimensional mesh */
#include "tet3drt0.h"	/* discretisation with Raviart-Thomas elements based on tet3d.h*/

#include "ddcluster.h"	/* Domain decomposition clustering*/
#include "hmatrix.h"	/* Hierarchical Matrices*/
#include "harith.h"	/* Arithmetic for hierarchical matrices*/

#include "truncation.h" /* Auxiliary function for truncation*/
#include "hcoarsen.h"	/* Coarsening of hierarchical matrices */
#include "matrixnorms.h"	/* Computing difference norms for two matrices */

real
norm2lu_sparsematrix(pchmatrix LU, pcsparsematrix sp)
{
  avector tmp1, tmp2;
  uint rows = LU->rc->size;
  uint cols = LU->cc->size;
  pavector x, y;
  real norm;
  uint j;

  assert(sp->rows == rows);
  assert(sp->cols == cols);

  x = init_avector(&tmp1, cols);
  y = init_avector(&tmp2, rows);

  random_avector(x);
  norm = norm2_avector(x);
  scale_avector(1.0 / norm, x);
  
  for(j=0; j < NORM_STEPS ; j++) {
   // printf("norm = %g \n", sqrt( norm));
    clear_avector(y);
    mvm_sparsematrix_avector(1.0, false, sp ,x, y);
    triangularsolve_hmatrix_avector(true ,true, false, LU, y);
    triangularsolve_hmatrix_avector(false,false, false, LU, y);
    add_avector(-1.0, y, x);
    copy_avector(x, y);
    triangularsolve_hmatrix_avector(false,false, true, LU, y);
    triangularsolve_hmatrix_avector(true ,true , true, LU, y);
    mvm_sparsematrix_avector(-1.0, true, sp ,y, x);
    norm = norm2_avector(x);
    scale_avector(1.0/norm, x);
  }
  
  uninit_avector(y);
  uninit_avector(x);

  return REAL_SQRT(norm);
}

real
norm2chol_sparsematrix(pchmatrix LU, pcsparsematrix sp)
{
  avector tmp1, tmp2;
  uint rows = LU->rc->size;
  uint cols = LU->cc->size;
  pavector x, y;
  real norm;
  uint j;

  assert(sp->rows == rows);
  assert(sp->cols == cols);

  x = init_avector(&tmp1, cols);
  y = init_avector(&tmp2, rows);

  random_avector(x);
  norm = norm2_avector(x);
  scale_avector(1.0 / norm, x);
  
  for(j=0; j < NORM_STEPS ; j++) {
   // printf("norm = %g \n", sqrt( norm));
    clear_avector(y);
    mvm_sparsematrix_avector(1.0, false, sp ,x, y);
    triangularsolve_hmatrix_avector(true, false, false, LU, y);
    triangularsolve_hmatrix_avector(true, false, true,  LU, y);
    add_avector(-1.0, y, x);
    copy_avector(x, y);
    triangularsolve_hmatrix_avector(true, false, false, LU, y);
    triangularsolve_hmatrix_avector(true, false , true, LU, y);
    mvm_sparsematrix_avector(-1.0, true, sp ,y, x);
    norm = norm2_avector(x);
    scale_avector(1.0/norm, x);
  }
  
  uninit_avector(y);
  uninit_avector(x);

  return REAL_SQRT(norm);
}
  
/*****************************************
 *Boundary conditions for domain
 ******************************************/
void set_boundary_unitsquare_tri2drt0(ptri2drt0 rt0){
  
  pctri2d t2 = rt0->t2;
  real(*x)[2] = t2->x ;
  uint edges = t2->edges;
  uint (*e)[2] = t2->e;
  uint xt[2];
  
  uint i;
  for(i=0; i<edges;i++){
    /*start and end vertex of edge i*/
    xt[0] = e[i][0]; xt[1] = e[i][1];
    if(x[xt[0]][0] == -1 && x[xt[1]][0] == -1){
    //  printf("left Dirichlet \n");
     /*e in left boundary, Dirichlet*/ 
    }
    else if(x[xt[0]][0] == 1 && x[xt[1]][0] == 1){
      /*e in right boundary, Dirichlet*/
      //printf("right Dirichlet\n");
    }
    else if(x[xt[0]][1] == -1 && x[xt[1]][1] == -1){
     /*e in lower boundary, Neumann*/ 
      rt0->is_dof[i] = 2;
     // printf(" lower Neumann\n");
    }
    else if(x[xt[0]][1] == 1 && x[xt[1]][1] == 1){
      /*e in upper boundary, Neumann*/
      rt0->is_dof[i] = 2;
      //printf("upper Neumann\n");
    }
   
  }
 
}
void set_boundary_unitcube_tet3drt0(ptet3drt0 rt0){
  
  pctet3d t3 = rt0->t3;
  uint faces = t3->faces;
  real (*x)[3] = t3->x;
  real xt[3][3];
  uint v[3];
  uint i, j;
  
  for(i=0; i<faces;i++){
   // printf("face %u\t", i);
    /*Get vertices of face i*/
    getvertices_face_tet3d(t3, i, v);
   /*Get coordinates of the vertices of face i*/
    for(j=0; j<3;j++){
      xt[j][0] = x[v[j]][0]; 
      xt[j][1] = x[v[j]][1];
      xt[j][2] = x[v[j]][2];
      //printf("vertex %u: %f %f %f\t", j, xt[j][0], xt[j][1], xt[j][2]);
     }
    /*start and end vertex of edge i*/
    /*Boden, z=0*/
    if(xt[0][2]==0 && xt[1][2]==0 && xt[2][2]==0)
      rt0->is_dof[i] = 2;
    /*Dach, z=1*/
    else if(xt[0][2]==1 && xt[1][2]==1 && xt[2][2]==1)
      rt0->is_dof[i] = 2;
    /*vorne, y=0*/
    else if(xt[0][1]==0 && xt[1][1]==0 && xt[2][1]==0)
      rt0->is_dof[i] = 2;
    /*hinten, y=1*/
    else if(xt[0][1]==1 && xt[1][1]==1 && xt[2][1]==1)
      rt0->is_dof[i] = 2;
    
    //printf("%u\n", rt0->is_dof[i]);
  }
}
/*************************************************
 *Sets values for the permeability
 *************************************************/
void function_permeability_2d(ptri2d t2, pavector k)
{
  uint i;
  
  for(i=0;i<k->dim;i++)
    k->v[i] = 1.0;//1e-2;
}
void function_permeability_3d(ptet3d t3, pavector k)
{
  uint i;
  
  for(i=0;i<k->dim;i++)
    k->v[i] = 1.0;//1e2;
}

int main (int argc, char **argv) {
  
  uint prob;		/* Used for selecting the problem: 2d or 3d? */
  uint L;		/* Number of grid refinements */
  uint clf;		/* Leafsize */
  uint clustermode;	/* Used for selecting the cluster strategy */ 
  uint decomp;		/* Used for selecting the type of the decomposition */
  uint i, j;		/* Auxiliary variables for loops */
  ptri2d *gr_2d;	/* 2d mesh hierarchy */
  ptri2drt0 rt0_2d;	/* Raviart_Thomas functions in 2d */	
  ptet3d *gr_3d;	/* 3d mesh hierarchy */
  ptet3drt0 rt0_3d;	/* Raviart-Thomas functions in 3d */
  psparsematrix sp_A, sp_B;	/* Sparsematrix object*/
  pavector k;		/* Vector for permeabilities */
  uint ndof;		/* Number of degree of freedom */
  uint *idx_B_rows, *idx_B_cols;	/* Array of indices for the clustergeometry */
  pclustergeometry cg_B_rows, cg_B_cols;/* Clustergeometry object */
  uint dim;		/* Dimension for domain decomposition clustering */
  pcluster root_B_cols, root_B_rows;	/* Cluster tree object */
  pblock broot_A, broot_B;		/* Block cluster tree object */
  real eta;		/* Admissibilty parameter*/
  phmatrix hm_A, hm_B;	/* Hierarchical matrices*/
  real tol_decomp;	/* Tolerance for decompositions*/
  real tol_coarsen;	/* Tolerance for coarsening of the hierarchical matrix */
  phmatrix l;		/* Hierarchical matrices for storing the cholesky factor*/
  ptruncmode tm;	/* Truncmode for truncations within the arithmetic */
  real error;		/* Auxiliary variable for testing*/
  pstopwatch sw;	/* Stopwatch for time measuring */
  real time;		/* Variable for time measurement */
  size_t sz, sz_rk, sz_full; /* Variables for memory footprint */
  uint elements;	/*Auxiliary variable: numberof elements*/
  
  /* First initialise the library */
  init_h2lib(&argc, &argv);	
  
  /*Initialise variables*/  
  eta = 2.0;
  tol_coarsen = 1e-16;
  tm = new_releucl_truncmode();
  sw = new_stopwatch();
  
  /* Set problme parameters */
  printf("Select problem\n");
  prob = askforint("1 = fem2d,\t 2 = fem3d\t \n","h2lib_fem", 1);
  L = askforint("Refinement?\n", "h2lib_L", 4);
  clf = askforint("Leafsize?\n", "h2lib_clf", 32);
  clustermode = askforint("1 = adaptive, \t 2 = adaptive from cluster", "h2lib_clusterstrategy", 2);
  decomp = askforint("1 = LU-decomposition of hm_A,\t, 2 = cholesky-decomposition of hm_A\n","h2lib_decomp",1);
  tol_decomp = askforreal("tol_decomp=?\n", "h2lib_tol_decomp", 1e-13);
  
  /* Build geomtery and discretisation for Darcy's equation */
  if(prob==1){
    printf("========================================\n"
	   "  Create and fill fem2d sparsematrix\n");
    /* Mesh hierarchy */
    gr_2d = (ptri2d *) allocmem((size_t) sizeof(ptri2d) * (L+1)); 
    gr_2d[0] = new_unitsquare_tri2d(); /* Set domain */
    //gr_2d[0]=new_unitcircle_tri2d();
    //gr_2d[0] = new_lshape_tri2d();
    for(i=0;i<L;i++){ /* Mesh refinements */
     gr_2d[i+1] = refine_tri2d(gr_2d[i], NULL); 
     fixnormals_tri2d(gr_2d[i]);
    }
    fixnormals_tri2d(gr_2d[i]);
    check_tri2d(gr_2d[L]); /*Check mesh for inconsistencies */
    rt0_2d = new_tri2drt0(gr_2d[L]); /*Build discretisation */
    set_boundary_unitsquare_tri2drt0(rt0_2d);  
    update_tri2drt0(rt0_2d);
    sp_A = build_tri2drt0_A_sparsematrix(rt0_2d); /*Build corresponding sparsematrices */
    sp_B = build_tri2drt0_B_sparsematrix(rt0_2d);
    k = new_avector(gr_2d[L]->triangles); /*Build and fill vector for the permeabilities */
    function_permeability_2d(gr_2d[L], k);
    assemble_tri2drt0_darcy_A_sparsematrix(rt0_2d, sp_A, 0, k); /*Fill the sparsematrices */
    assemble_tri2drt0_darcy_B_sparsematrix(rt0_2d, sp_B, 0);
    ndof = rt0_2d->ndof;printf("ndof = %u\n", ndof);
    /*Initialise index arrays for the cluster geometries */
    idx_B_rows = allocuint(rt0_2d->t2->triangles); 
    for(i=0;i<rt0_2d->t2->triangles;i++)
      idx_B_rows[i]  = i;
    idx_B_cols = allocuint(rt0_2d->ndof);
    for(i=0;i<rt0_2d->ndof;i++)
      idx_B_cols[i] = i;
    /* Build clustergeomtries for the problem */
    cg_B_rows = build_tri2drt0_B_clustergeometry(rt0_2d, idx_B_rows);
    cg_B_cols = build_tri2drt0_A_clustergeometry(rt0_2d, idx_B_cols);
    elements = gr_2d[L]->triangles;
    del_tri2drt0(rt0_2d);
    for(i=0;i<=L;i++){
      j =  L-i;
      del_tri2d(gr_2d[j]);
    }
    freemem(gr_2d);
    del_avector(k);
    dim = 2; /* Set dimenison for domain decomposition clustering */
  }
  else { assert(prob==2); 
    printf("========================================\n"
	   "  Create and fill fem3d sparsematrix\n");
    /* Mesh hierarchy */
    gr_3d = (ptet3d *) allocmem((size_t) sizeof(ptet3d) * (L+1));
    gr_3d[0] = new_unitcube_tet3d(); /* Set domain */
    for(i=0;i<L;i++){
      gr_3d[i+1] = refine_tet3d(gr_3d[i],NULL); 
      fixnormals_tet3d(gr_3d[i]);
    }
    fixnormals_tet3d(gr_3d[i]);
    check_tet3d(gr_3d[L]); /*Check mesh for inconsistencies */
    rt0_3d = new_tet3drt0(gr_3d[L]); /*build discretisation */
    set_boundary_unitcube_tet3drt0(rt0_3d);  
    update_tet3drt0(rt0_3d);
    sp_A = build_tet3drt0_A_sparsematrix(rt0_3d); /*Build corresponding sparsematrices*/
    sp_B = build_tet3drt0_B_sparsematrix(rt0_3d);
    k = new_avector(gr_3d[L]->tetrahedra); /* Build and fill vector for the permeabilities */
    function_permeability_3d(gr_3d[L], k);
    assemble_tet3drt0_darcy_A_sparsematrix(rt0_3d, sp_A, 0, k); /*Fill the sparsematrices*/
    assemble_tet3drt0_darcy_B_sparsematrix(rt0_3d, sp_B, 0);
    ndof = rt0_3d->ndof;printf("ndof = %u\n", ndof);
    /*Initialise index arrays for the cluster geometries */
    idx_B_rows = allocuint(rt0_3d->t3->tetrahedra); 
    for(i=0;i<rt0_3d->t3->tetrahedra;i++)
      idx_B_rows[i]  = i;
    idx_B_cols = allocuint(rt0_3d->ndof);
    for(i=0;i<rt0_3d->ndof;i++)
      idx_B_cols[i] = i;
    /* Build clustergeomtries for the problem */
    cg_B_rows = build_tet3drt0_B_clustergeometry(rt0_3d, idx_B_rows);
    cg_B_cols = build_tet3drt0_A_clustergeometry(rt0_3d, idx_B_cols);
    elements = gr_3d[L]->tetrahedra;
    del_tet3drt0(rt0_3d);
    for(i=0;i<=L;i++){
      j = L -i;
      del_tet3d(gr_3d[j]);
    }
    freemem(gr_3d);
    del_avector(k);
    dim = 3;
  }
  
  /* Build and fill the corresponding hierarchical matrices*/
  printf("========================================\n"
	 "  Create and fill hierarchical matrix\n");
  if(clustermode ==1){
    /* Build cluster trees based on adaptive clustering*/
    root_B_rows = build_cluster(cg_B_rows, elements, idx_B_rows, clf, H2_ADAPTIVE); 
    root_B_cols = build_cluster(cg_B_cols, ndof,     idx_B_cols, clf, H2_ADAPTIVE);
    /*Build block cluster trees using the euclidian (maximum) admissibilty condition*/
    broot_B = build_nonstrict_block(root_B_rows, root_B_cols, &eta, admissible_2_cluster);
    broot_A = build_nonstrict_block(root_B_cols, root_B_cols, &eta, admissible_2_cluster);
    /*Build hierarchical matrices from block*/
    hm_A = build_from_block_hmatrix(broot_A, 0);
    hm_B = build_from_block_hmatrix(broot_B, 0);
    /*Fill hierarchical matrices with sparsematrix entries*/
    copy_sparsematrix_hmatrix(sp_A, hm_A);
    copy_sparsematrix_hmatrix(sp_B, hm_B);
    /* Compute error */
    printf("sp_A: rows %u cols: %u\n", sp_A->rows, sp_A->cols);
    printf("sp_B: rows %u cols: %u\n", sp_B->rows, sp_B->cols);
    printf("hm_A: rows %u cols: %u\n", hm_A->rc->size, hm_A->cc->size);
    printf("hm_B: rows %u cols: %u\n", hm_B->rc->size, hm_B->cc->size);
    error = norm2diff_sparsematrix_hmatrix(hm_A, sp_A);
    printf("error A = %g\n", error);
    error = norm2diff_sparsematrix_hmatrix(hm_B, sp_B);
    printf("error B = %g\n", error);
    del_block(broot_A);
    del_block(broot_B);
  }
  else{ /*clustermode == 2*/
    /* Build cluster tree based on adaptive clustering */
    root_B_rows = build_cluster(cg_B_rows, elements, idx_B_rows, clf, H2_ADAPTIVE); 
    /* Build cluster tree based adaptive doamin decomposition clustering */
    root_B_cols = build_adaptive_fromcluster_cluster(root_B_rows, cg_B_cols, ndof, idx_B_cols, clf, dim);
   // freemem(flag);
    /* Build block cluster trees using the domain decompostion admissibilty (rt0) condition*/
    broot_B = build_nonstrict_block(root_B_rows, root_B_cols, &eta, admissible_2_cluster);
    //broot_B = build_nonstrict_block(root_B_rows, root_B_cols, &eta, admissible_dd_rt0_cluster);
    broot_A = build_nonstrict_block(root_B_cols, root_B_cols, &eta, admissible_dd_cluster);
    /* Build hierarchical matrices from block cluster tree*/
    hm_A = build_from_block_hmatrix(broot_A, 0);
    hm_B = build_from_block_hmatrix(broot_B, 0);
    /* Fill hierarchical matrices with sparsematrix entries*/
    copy_sparsematrix_hmatrix(sp_A, hm_A);
    copy_sparsematrix_hmatrix(sp_B, hm_B);
    /* Compute error */
    printf("sp_A: rows %u cols: %u\n", sp_A->rows, sp_A->cols);
    printf("sp_B: rows %u cols: %u\n", sp_B->rows, sp_B->cols);
    printf("hm_A: rows %u cols: %u\n", hm_A->rc->size, hm_A->cc->size);
    printf("hm_B: rows %u cols: %u\n", hm_B->rc->size, hm_B->cc->size);
    error = norm2diff_sparsematrix_hmatrix(hm_A, sp_A);
    printf("error A = %g\n", error);
    error = norm2diff_sparsematrix_hmatrix(hm_B, sp_B);
    printf("error B = %g\n", error);
    del_block(broot_A);
    del_block(broot_B);
  }
    
    /* Draw hierarchical matrix with cairo */
#if 0
    #ifdef USE_CAIRO
    cairo_t *cr;
    printf("Draw hmatrix to \"hm_p1.pdf\"\n");
    cr = new_cairopdf("../hm_p1.pdf", 1024.0, 1024.0);
    draw_cairo_hmatrix(cr, hm, false, 0);
    cairo_destroy(cr);
    #endif
#endif

  /*Compute size of the hierarchical matrix A*/
  sz = getsize_hmatrix(hm_A);
  sz_rk = getfarsize_hmatrix(hm_A);
  sz_full = getnearsize_hmatrix(hm_A);
  printf("matrix A:\t \t %.2f MB \t %.2f KB \t %.2f KB/DoF\n",
	 sz/1024.0/1024.0, sz/1024.0, sz/1024.0/ndof);
  printf("size rk: \t %.2f MB\t size full: \t %.2f MB \n", sz_rk/1024.0/1024.0, sz_full
            /1024.0/1024.0);
  printf("Coarsen hmatrix\n");
  
  /*Coarsening of the hierarchical matrix A to safe storage */
  coarsen_hmatrix(hm_A, tm, tol_coarsen, true);
  
  /* Compute size of the hierarchical matrix A after coarsening*/
  sz = getsize_hmatrix(hm_A);
  sz_rk = getfarsize_hmatrix(hm_A);
  sz_full = getnearsize_hmatrix(hm_A);
  printf("matrix A:\t \t %.2f MB \t %.2f KB \t %.2f KB/DoF\n",
	   sz/1024.0/1024.0, sz/1024.0, sz/1024.0/ndof);
  printf("size rk: \t %.2f MB\t size full: \t %.2f MB \n", sz_rk/1024.0/1024.0, sz_full
            /1024.0/1024.0);
  /* Compute error after coarseing */  
  error = norm2diff_sparsematrix_hmatrix(hm_A, sp_A);
  printf("error = %g\n", error);
  
  /*Compute size of the hierarchical matrix B */
  sz = getsize_hmatrix(hm_B);
  sz_rk = getfarsize_hmatrix(hm_B);
  sz_full = getnearsize_hmatrix(hm_B);
  printf("matrix B:\t \t %.2f MB \t %.2f KB \t %.2f KB/DoF\n",
	 sz/1024.0/1024.0, sz/1024.0, sz/1024.0/ndof);
  printf("size rk: \t %.2f MB\t size full: \t %.2f MB \n", sz_rk/1024.0/1024.0, sz_full
            /1024.0/1024.0);
  printf("Coarsen hmatrix\n");
  
  /*Coarsening of the hierarchical matrix B to safe storage */
  coarsen_hmatrix(hm_B, tm, tol_coarsen, true);
  
  /* Compute size of the hierarchical matrix B after coarsening*/
  sz = getsize_hmatrix(hm_B);
  sz_rk = getfarsize_hmatrix(hm_B);
  sz_full = getnearsize_hmatrix(hm_B);
  printf("matrix B:\t \t %.2f MB \t %.2f KB \t %.2f KB/DoF\n",
	   sz/1024.0/1024.0, sz/1024.0, sz/1024.0/ndof);
  printf("size rk: \t %.2f MB\t size full: \t %.2f MB \n", sz_rk/1024.0/1024.0, sz_full
            /1024.0/1024.0);
  /* Compute error after coarseing */  
  error = norm2diff_sparsematrix_hmatrix(hm_B, sp_B);
  printf("error = %g\n", error);
  
  /* Compute decomposition of the hierarchical matrix A */
  if(decomp == 1){
    printf("========================================\n"
	 "  Test lu-decomposition (A)\n");
    /*Compute size of the hierarchical matrix*/
    sz = getsize_hmatrix(hm_A);
    printf("size original matrix:\t \t %.2f MB \t %.2f KB \t %.2f KB/DoF\n",
	   sz/1024.0/1024.0, sz/1024.0, sz/1024.0/ndof);
    /*Compute lu-decomposition of the hierarchical matrix*/
    start_stopwatch(sw);
    lrdecomp_hmatrix(hm_A, 0, tol_decomp); 
    time = stop_stopwatch(sw);
    /*Compute size of the lu-decomposition*/
    sz = getsize_hmatrix(hm_A);
    sz_rk = getfarsize_hmatrix(hm_A);
    sz_full = getnearsize_hmatrix(hm_A);
    printf("size lu-decomposition (lu):\t %.2f MB \t %.2f KB \t %.2f KB/DoF\n", 
	   sz/1024.0/1024.0, sz/1024.0, sz/1024.0/ndof);
    printf("size rk: \t %.2f MB\t size full: \t %.2f MB \n", sz_rk/1024.0/1024.0, sz_full
            /1024.0/1024.0);
    printf("time = %f s\n", time);

    /*Compute error*/
    error = norm2lu_sparsematrix(hm_A, sp_A);
    printf("||I-(lu)^-1 sp||_2 ");
    printf("error = %g\n", error);
  }
  else{ /*decomp == 2*/
    printf("========================================\n"
	 "  Test cholesky-decomposition (A)\n");
    /*Compute size of the hierarchical matrix*/
    sz = getsize_hmatrix(hm_A);
    printf("size original matrix:\t %.2f MB \t %.2f KB \t %.2f KB/DoF\n",
	   sz/1024.0/1024.0, sz/1024.0, sz/1024.0/ndof);
    /* Compute cholesky-decomposition*/
    start_stopwatch(sw);
    choldecomp_hmatrix(hm_A, 0, tol_decomp);
    time = stop_stopwatch(sw);
    /* Compute size of the cholesky-factor*/
    l = clone_lower_hmatrix(false, hm_A);
    sz = getsize_hmatrix(l);
    sz_rk = getfarsize_hmatrix(l);
    sz_full = getnearsize_hmatrix(l);
    printf("size chol-decomposition (l):\t %.2f MB \t %.2f KB \t %.2f KB/DoF\n", 
	   sz/1024.0/1024.0, sz/1024.0, sz/1024.0/ndof);
    printf("size rk: \t %.2f MB\t size full: \t %.2f MB \n", sz_rk/1024.0/1024.0, sz_full
            /1024.0/1024.0);
    printf("time = %f s\n", time);

    /* Compute error */
    error = norm2chol_sparsematrix(hm_A, sp_A);
    printf("||I-(ll*)^-1 sp||_2 ");
    printf("error = %g\n", error);

    del_hmatrix(l);
  }
  
  /*Cleaning up*/
  del_truncmode(tm);
  del_stopwatch(sw);
  del_hmatrix(hm_A); del_hmatrix(hm_B);
  
  
  del_cluster(root_B_cols); del_cluster(root_B_rows);
  del_clustergeometry(cg_B_rows); del_clustergeometry(cg_B_cols);
  del_sparsematrix(sp_A); del_sparsematrix(sp_B);
  freemem(idx_B_rows); freemem(idx_B_cols);
  (void) printf("  %u matrices and\n"
    "  %u vectors still active\n",
    getactives_amatrix(),
    getactives_avector());
  uninit_h2lib();
  printf("The end\n");
  return 0;
}
