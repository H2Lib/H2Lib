/*------------------------------------------------------------------------*/
/*  How to create an DH²-matrix approximation for the Helmholtz equation  */
/*------------------------------------------------------------------------*/

#ifdef USE_COMPLEX
/* Important, for Helmholtz problems complex numbers are necessary 
(To change, use the "options.inc" / DUSE_COMPLEX) */

#include <stdio.h>
#include <stdlib.h>

#include "basic.h"		/* Initialize library and for basic definitions */
#include "surface3d.h"		/* For surface triangulation */
#include "macrosurface3d.h"	/* For surface prototype */
#include "bem3d.h"		/* For boundary element method objects */
#include "helmholtzbem3d.h"	/* Helmholtz depending BEM */
#include "clustergeometry.h"	/* For geometry objects (point sets) corresponding to clusters */
#include "dblock.h"		/* For directional blocks */
#include "dclusterbasis.h"	/* For directional cluster basis */

int
main(int argc, char **argv)
{

  psurface3d gr;		/* Grid for BEM */
  pmacrosurface3d mg;		/* Used to create grid */
  real      k_wave = 2.0;	/* Wave number */
  pbem3d    bem;		/* BEM object */
  uint      nq = 3;		/* Number of quadrature points */
  uint      ni = 4;		/* Number of interpolation points */
  pclustergeometry cgeo;	/* Object for the cluster geometry */
  uint     *idx;		/* Array of indices for the cluster geometry */
  uint      res = 32;		/* Resolution for the cluster tree */
  pcluster  root;		/* Cluster tree object */
  real      eta1 = 1.0;		/* Directional admissibility parameter */
  pleveldir ld;			/* Object for level-wise set of directions */
  pdcluster droot;		/* Directional cluster tree object */
  real      eta2 = 1.0;		/* Parabolic admissibility parameter */
  diradmdata dad;		/* Set of parameters for admissibility condition */
  pdblock   dbroot;		/* Directional block tree */
  pdclusterbasis rb, cb;	/* Directional row and column cluster basis */
  pdh2matrix Gh;		/* DH²-matrix approximation */


  init_h2lib(&argc, &argv);	/* First initialize the libaray */
  /* Creating surface */
  printf("  Start building DH^2 - matrix\n");
  mg = new_sphere_macrosurface3d();	/* Choose the unit sphere as surface */
  gr = build_from_macrosurface3d_surface3d(mg, 10);	/* Split the triangles '10' times */
  del_macrosurface3d(mg);	/* Delete redundant macrosurface */
  prepare_surface3d(gr);	/* Prepare grid (evaluate normal vector etc.) */
  /* Creating BEM object */
  bem = new_slp_helmholtz_bem3d((field) k_wave, gr, nq, nq + 2, BASIS_CONSTANT_BEM3D, BASIS_CONSTANT_BEM3D);	/* Set up BEM object (this time single layer potential, orthers are possible -> 'helmholtzbem3d.h') */
  /* Creating cluster tree */
  cgeo = build_bem3d_clustergeometry(bem, &idx, BASIS_CONSTANT_BEM3D);	/* Create cluster geometry for the BEM objekt with respect to the chosen basis functions and prepare it for clustering */
  update_point_bbox_clustergeometry(cgeo, cgeo->nidx, idx);
  root = build_adaptive_cluster(cgeo, cgeo->nidx, idx, res);	/* Build a cluster tree object from the cluster geometry and indices */
  del_clustergeometry(cgeo);	/* Delete redundant cluster geometry */
  /* Creating directional cluster tree */
  droot = buildfromcluster_dcluster(root);
  ld = builddirections_box_dcluster(droot, eta1 / k_wave);	/* Build set of directions */
  /* Creating directional block tree */
  dad.eta1 = eta1;
  dad.eta2 = eta2;
  dad.wave_k = k_wave;		/* Fill set of parameters */
  dad.xy = allocreal(3);
  dad.ld = ld;
  dbroot = build_dblock(droot, droot, 0, parabolic_admissibility, &dad);	/* Create directional block tree */
  freemem(dad.xy);
  /* Creating directional cluster basis */
  rb = buildfromdcluster_dclusterbasis(droot);	/* Build directional cluster basis */
  cb = buildfromdcluster_dclusterbasis(droot);
  findranks_dclusterbasis(ni * ni * ni, dbroot, rb, cb);	/* Only if a matrix in direction 'i' is used, its rank is set to 'ni*ni*ni' otherwise its rank is set to zero */
  initmatrices_dclusterbasis(rb);	/* Initialize cluster basis matrices */
  initmatrices_dclusterbasis(cb);
  /* Creating DH²-matrix and filling directional cluster basis */
  setup_dh2matrix_aprx_inter_bem3d(bem, rb, cb, dbroot, ni);	/* Simplest set up routine, more effective ones are available (take a look at 'bem3d.h') */
  Gh = buildfromblock_dh2matrix(dbroot, rb, cb);	/* Build DH²-matrix */
  assemble_bem3d_dh2matrix_row_dclusterbasis(bem, rb);	/* Fill directional cluster basis */
  assemble_bem3d_dh2matrix_col_dclusterbasis(bem, cb);
  /* Filling DH²-matrix */
  assemble_bem3d_farfield_dh2matrix(bem, Gh);
  assemble_bem3d_nearfield_dh2matrix(bem, Gh);
  printf("  Building DH^2 - matrix complete!\n");

  /* Cleaning up */
  del_helmholtz_bem3d(bem);
  del_dh2matrix(Gh);
  del_dclusterbasis(rb);
  del_dclusterbasis(cb);
  del_leveldir(ld);
  del_dblock(dbroot);
  del_dcluster(droot);
  del_cluster(root);
  freemem(idx);
  del_surface3d(gr);

  /* Test, if all matrices and vectors are deleted */
  printf
    ("\n  %u directional cluster bases,\n  %u cluster operators, \n  %u matrices and\n"
     "  %u vectors still active\n", getactives_dclusterbasis(),
     getactives_dclusteroperator(), getactives_amatrix(),
     getactives_avector());
  /* Last thing, uninitialize the library */
  uninit_h2lib();

  return EXIT_SUCCESS;
}

#else

#include <stdio.h>

int
main(int argc, char **argv)
{

  fprintf(stderr,
	  "Need complex numbers, activate with '-DUSE_COMPLEX' (options.inc) \n");

  return 0;
}

#endif
