/*    Simple test to demonstrate the application of functions from the visualize file   */

#ifdef USE_COMPLEX
/* Important, for Helmholtz problems complex numbers are necessary 
(To change, use the "options.inc" / DUSE_COMPLEX) */

#include <stdlib.h>
#include <stdio.h>
#include "basic.h"
#include "dblock.h"
#include "dclusterbasis.h"
#include "surface3d.h"
#include "macrosurface3d.h"
#include "clustergeometry.h"
#include "parameters.h"
#include "helmholtzbem3d.h"
#include "visualize.h"

/*******************************************************************/
/*  Use test_type to change between objects to visualize           */
/*                                                                 */
/* 0 = surface as example for tri2d, tet3d and surface objects     */
/* 1 = cluster as example for cluster and dcluster                 */
/* 2 = dblock as example for block and dblock                      */
/* 3 = boundary values on surface                                  */
/* 4 = solution of the Helmholtz equation or animated solution     */
/*******************************************************************/




int
main(int argc, char **argv)
{

  psurface3d gr;		/* BEM grid */
  pclustergeometry cgeo;	/* Description of cluster geometry */
  uint     *idx;		/* Index permutation */
  pcluster  root;		/* Root of cluster tree */
  pdcluster droot;		/* Root of directional cluster tree */
  pleveldir ld;			/* Levelwise directions */
  uint      res;		/* Cluster resolution */
  real      wave_k;		/* Wave number */
  real      eta1;		/* Directional admissibility parameter */
  real      eta2;		/* Parabolic admissibility parameter */
  real      eta_eff;		/* Effective admissibility parameter */
  uint      nq;			/* Quadrature order */
  uint      ni;			/* Interpolation order */
  diradmdata dad;		/* Parameters for admissibility condition */
  pdblock   broot;		/* Root of block tree */
  pbem3d    bem;		/* BEM object with parameters and callback functions */
  pavector  val;		/* Vector with values for the BEM object on the surface */

  uint      test_type;		/* Variable used to change visualize function for the test */
  uint      dirs;
  uint      i;




  test_type = 2;		/* Choose test_type value */

  init_h2lib(&argc, &argv);	/* Initilize the library */

  printf("----------------------------------------\n");
  printf("Preparing surface\n ");	/* Surface: typing s6 e.g. means sphere with 6 splits */
  gr = build_interactive_surface3d();
  prepare_surface3d(gr);
  printf("  %u vertices, %u edges, %u triangles\n", gr->vertices,
	 gr->edges, gr->triangles);
  fflush(stdout);

  if (test_type == 0) {		/* Visualize surface: two possible functions for a surface object */
    //visualize_surface3d(gr, argc, argv);                                        /* -> shows the whole surface */
    visualize_certain_surface_triangle(gr, argc, argv);	/* -> made to choose one triangle out of all surface triangles */
  }

  /* Set up a BEM object with the following parameters */
  printf("Preparing Helmholtz BEM object\n");

  wave_k = askforreal("Wave number?", "h2lib_wavek", 8.0);
  eta1 =
    askforreal("Directional admissibility parameter?", "h2lib_eta1", 5.0);
  eta2 = askforreal("Parabolic admissibility parameter?", "h2lib_eta2", 1.0);
  ni = 4;			/* Order of interpolation */
  res = askforint("Cluster resolution?", "h2lib_leafsize", 2 * ni * ni * ni);

  nq = 3;			/* Order of quadrature */


  printf("Creating dbem3d object\n");	/* Prepare wave vector */
  bem = new_slp_helmholtz_bem3d((field) wave_k, gr, nq, nq + 2, BASIS_CONSTANT_BEM3D, BASIS_CONSTANT_BEM3D);	/* Set up new single layer bem object for the helmholtz equation */
  printf("  Wave number %g\n"
	 "  Interpolation order %u\n"
	 "  %u degrees of freedom\n", REAL(bem->k), ni, gr->triangles);
  fflush(stdout);

  printf("Creating cluster tree with resolution %u\n", res);

  cgeo = build_bem3d_clustergeometry(bem, &idx, BASIS_CONSTANT_BEM3D);
  root = build_adaptive_cluster(cgeo, cgeo->nidx, idx, res);
  printf("  %u clusters\n", root->desc);
  fflush(stdout);

  if (test_type == 1) {
    visualize_cluster_bbox(root, cgeo, true, argc, argv);	/* Visualize the cluster 'root' with corresponding clustergeometry cgeo */
  }

  printf("Creating directional cluster tree\n");	/* Compute the directional cluster tree for 'root' */
  droot = buildfromcluster_dcluster(root);
  ld = builddirections_box_dcluster(droot, eta1 / wave_k);
  dirs = getalldirections_dcluster(droot);
  for (i = 0; i <= ld->depth; i++)
    printf("  Level %2u: maxdiam %g, %u directions\n", i, ld->maxdiam[i],
	   ld->directions[i]);

  printf("  %u directions (%.1f per cluster)\n", dirs,
	 (real) dirs / root->desc);
  fflush(stdout);

  printf("Creating directional block tree\n");	/* And now the directional block tree for 'root' times 'root' */
  dad.eta1 = eta1;
  dad.eta2 = eta2;
  dad.wave_k = wave_k;
  dad.xy = allocreal(3);
  dad.ld = ld;
  broot = build_dblock(droot, droot, 0, parabolic_admissibility, &dad);
  eta_eff = getmaxeta_dblock(broot);
  freemem(dad.xy);
  printf("  %u blocks (%.1f per cluster)\n"
	 "  Effective eta %.3f (predicted rho %.3f)\n", broot->desc,
	 (real) broot->desc / root->desc, eta_eff,
	 1.0 / (REAL_SQRT(1.0 / REAL_SQR(eta_eff) + 1.0) + 1.0 / eta_eff));
  fflush(stdout);

  if (test_type == 2) {		/* This time we have three possible functions */
    visualize_dblock_certain_bbox(broot, cgeo, cgeo, argc, argv);	/* -> made to pick out one dblock and show the corresponding dclusters */
    //visualize_dblock_level_certain_bbox(broot, cgeo, cgeo, 1, getdepth_dblock(broot), argc, argv);    /* -> similar to the first one, but only for previously determined level */
    //visualize_dblock_bbox(broot, cgeo, cgeo, true, argc, argv);               /* -> made to show the whole block tree (levelwise) */

  }

  /* More complicate versions, only for the helmholtz equation and their solution */
  if (test_type == 3) {		/* Visualize boundary values on surface, for the sake of simplicity this time random values.... */
    val = new_avector(gr->triangles);	/* Vector for the random values */
    random_avector(val);	/* Fill the vector */
    visualize_boundaryvalue_surface_triangle(gr, val, argc, argv);	/* Show the surface, the color gradient is chosen accordingly to the entries of the vector val */
    del_avector(val);
  }

  if (test_type == 4) {		/* This time the solution will be visualized, both functions needs an additional file with values of the 
				   solution, a simple example file 'example_solutiondata.dat' is given  */
    //visualize_helmholtz_solution_surface(gr, argc, argv);                     /* The simple version only shows the solution */
    animate_helmholtz_solution(gr, argc, argv);	/* This one shows an animated version for the solution */
  }

  (void) printf("Cleaning up!\n");
  del_clustergeometry(cgeo);
  del_helmholtz_bem3d(bem);
  del_dblock(broot);
  del_dcluster(droot);
  del_cluster(root);
  freemem(idx);
  del_surface3d(gr);

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
