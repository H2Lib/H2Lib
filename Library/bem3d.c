/* ------------------------------------------------------------
 This is the file "bem3d.c" of the H2Lib package.
 All rights reserved, Sven Christophersen 2011
 ------------------------------------------------------------ */

/* C STD LIBRARY */
#include <string.h>
/* CORE 0 */
#ifdef USE_SIMD
#include "simd.h"
#endif
#include "basic.h"
/* CORE 1 */
/* CORE 2 */
/* CORE 3 */
/* SIMPLE */
/* PARTICLES */
/* BEM */
#include "bem3d.h"

/* ------------------------------------------------------------
 Structs and typedefs
 ------------------------------------------------------------ */

/*
 * This constant defines the minimal width of an interval, that is used for
 * interpolation. If an interval @f$[a,b]@f$ is smaller than \ref INTERPOLATION_EPS_BEM3D
 * it is expanded to @f$[a - 0.5 \cdot \texttt{INTERPOLATION\_EPS\_BEM3D},
 * b + 0.5 \cdot \texttt{INTERPOLATION\_EPS\_BEM3D} ]@f$ .
 */
#define INTERPOLATION_EPS_BEM3D 5.0e-3

/*
 * Just an abbreviation for the struct _greencluster3d .
 */
typedef struct _greencluster3d greencluster3d;
/*
 * Pointer to a @ref greencluster3d object.
 */
typedef greencluster3d *pgreencluster3d;
/*
 * Pointer to a constant @ref greencluster3d object.
 */
typedef const greencluster3d *pcgreencluster3d;

/*
 * Just an abbreviation for the struct _greencluster3d .
 */
typedef struct _greenclusterbasis3d greenclusterbasis3d;
/*
 * Pointer to a @ref greenclusterbasis3d object.
 */
typedef greenclusterbasis3d *pgreenclusterbasis3d;
/*
 * Pointer to a constant @ref greenclusterbasis3d object.
 */
typedef const greenclusterbasis3d *pcgreenclusterbasis3d;

/*
 * @brief Substructure used for approximating @ref _hmatrix "h-", @ref
 * _uniformhmatrix "uniformh-" and @ref _h2matrix "h2matrices".
 *
 * This struct contains many needed parameters for the different approximation
 * techniques such as interpolation or green based methods.
 */
struct _aprxbem3d {

  /*
   * @brief Number of Tschebyscheff-interpolation points in each spatial dimension.
   */
  uint      m_inter;

  /*
   * @brief One dimensional Tschebyscheff-interpolation points in [-1,1].
   */
  preal     x_inter;

  /*
   * @brief Rank for interpolation based methods: @f$ k = m^2 @f$ .
   */
  uint      k_inter;

  /*
   * @brief Number of interval segments for green-quadrature.
   */
  uint      l_green;

  /*
   * @brief Number of gaussian quadrature points for green based methods.
   */
  uint      m_green;

  /*
   * @brief Product of @f$ m_{green} @f$ and @f$ l_{green} @f$.
   */
  uint      ml_green;

  /*
   * @brief Rank of green based methods.
   *
   * Depending on the utilized parameterization
   * the rank can vary. In case of cube parameterization @ref
   * build_bem3d_cube_quadpoints the rank applies to
   * @f$ k = 12 \, m^2 \cdot l^2 @f$ .
   */
  uint      k_green;

  /*
   * @brief @f$ \delta @f$ controls the distance between the bounding box of the current
   * cluster and the parameterization.
   *
   * This parameter is used as a relative
   * value. In fact the distance will be computed as @f$ \delta = \widetilde
   * \delta \cdot \operatorname{diam\_max\_cluser(t)} @f$ if cube parameterization
   * @ref build_bem3d_cube_quadpoints is used.
   */
  real      delta_green;

  /*
   * @brief One dimensional gaussian quadrature points in [-1,1] for green based
   * methods.
   */
  preal     t_green;

  /*
   * @brief One dimensional gaussian quadrature weight in [-1,1] for green based
   * methods.
   */

  preal     w_green;

  /*
   * @brief This is a callback function of type @ref quadpoints3d.
   *
   * It defines which
   * parameterization should be used for green based approximation techniques.
   */
  quadpoints3d quadpoints;

  /*
   * @brief Additional Information for greenhybrid methods.
   *
   * When using the greenhybrid methods, one needs to save the pivot elements
   * and the matrix @f$ V_t @f$ for each cluster. There we clone the structure of a
   * clustertree in a new struct _greencluster3d "greencluster3d" and add
   * these information to the clusters. <tt>grc_green</tt> is responsible for
   * the row clustertree.
   */
  pgreencluster3d grc_green;

  /*
   * @brief Additional Information for greenhybrid methods.
   *
   * When using the greenhybrid methods, one needs to save the pivot elements
   * and the matrix @f$ W_t @f$ for each cluster. There we clone the structure of a
   * clustertree in a new struct _greencluster3d "greencluster3d" and add
   * these information to the clusters. <tt>gcc_green</tt> is responsible for
   * the column clustertree.
   */
  pgreencluster3d gcc_green;

  /*
   * @brief Additional Information for greenhybrid methods.
   *
   * When using the greenhybrid methods, one needs to save the pivot elements
   * and the matrix @f$ V_t @f$ for each cluster. There we clone the structure of a
   * clusterbasis in a new struct _greenclusterbasis3d "greenclusterbasis3d" and add
   * these information to the clusters. <tt>grb_green</tt> is responsible for
   * the row clusterbasis. The matrix @f$ V_t @f$ is stored in the usual way inside
   * the clusterbasis itself.
   */
  pgreenclusterbasis3d grb_green;

  /*
   * @brief Additional Information for greenhybrid methods.
   *
   * When using the greenhybrid methods, one needs to save the pivot elements
   * and the matrix @f$ W_t @f$ for each cluster. There we clone the structure of a
   * clusterbasis in a new struct _greenclusterbasis3d "greenclusterbasis3d" and add
   * these information to the clusters. <tt>gcb_green</tt> is responsible for
   * the column clusterbasis. The matrix @f$ W_t @f$ is stored in the usual way inside
   * the clusterbasis itself.
   */
  pgreenclusterbasis3d gcb_green;

  /*
   * @brief Accuracy for ACA based algorithms
   */
  real      accur_aca;

  /*
   * @brief This flag indicated if blockwise recompression technique should be used or
   * not.
   *
   * Default value is <tt>recomp = false</tt>
   */
  bool      recomp;

  /*
   * @brief Accuracy for blockwise recompression.
   *
   * If <tt>recomp == true</tt> then <tt>accur_recomp</tt> controls the minimum
   * accuracy for each block within the @ref _hmatrix "hmatrix".
   */
  real      accur_recomp;

  /*
   * @brief This flag indicates whether \"coarsening\" should be used while constructing
   * @ref _hmatrix "hmatrices" or not.
   */
  bool      coarsen;

  /*
   * @brief Accuracy for coarsening.
   * If <tt>coarsen == true</tt> then <tt>accur_coarsen</tt> controls the minimum
   * accuracy for each coarsening step.
   */
  real      accur_coarsen;

  /*
   * @brief This flag indicated whether \"hierarchical compression\" should be used to
   * construct @ref _h2matrix "h2matrices"
   */
  bool      hiercomp;

  /*
   * @brief Accuracy for hierarchical recompression.
   *
   * If <tt>hiercomp == true</tt> then <tt>accur_hiercomp</tt> controls the minimum
   * accuracy for hierarchical compression technique.
   */
  real      accur_hiercomp;

  /*
   * @brief Additional information used by truncation-routines.
   */
  ptruncmode tm;
};

struct _parbem3d {
  /*
   * special members for different H- and H2-matrix approximations
   */

  phmatrix *hn;			/* temporary enumerated list of hmatrices */
  ph2matrix *h2n;		/* temporary enumerated list of h2matrices */
  pdh2matrix *dh2n;		/* temporary enumerated list of h2matrices */
  pclusterbasis *rbn;
  pclusterbasis *cbn;
  pdclusterbasis *drbn;
  pdclusterbasis *dcbn;
  pclusteroperator *rwn;	/* temporary enumerated list of clusteroperators for row-cluster */
  pclusteroperator *cwn;	/* temporary enumerated list of clusteroperators for col-cluster */
  pdclusteroperator *ron;	/* temporary enumerated list of dclusteroperators for row-cluster */
  pdclusteroperator *con;	/* temporary enumerated list of dclusteroperators for col-cluster */
  uint     *leveln;		/* temporary list of levelnumber for each block in blocktree. */
  pgreencluster3d *grcn;
  uint      grcnn;
  pgreencluster3d *gccn;
  uint      gccnn;
  pgreenclusterbasis3d *grbn;
  uint      grbnn;
  pgreenclusterbasis3d *gcbn;
  uint      gcbnn;
};

struct _greencluster3d {
  uint     *xi;
  /** local indices of pivot elements */
  uint     *xihat;
  /** global indices of pivot elements */
  pamatrix  V;
  /** Interpolation operator */
  pccluster t; /** corresponding cluster */
  uint      sons;
/** number of sons for current cluster */
};

struct _greenclusterbasis3d {
  uint     *xi;
  /** local indices of pivot elements */
  uint     *xihat;
  /** global indices of pivot elements */
  pamatrix  Qinv;
  /** Triangular factor of QR-decomposition */
  pcclusterbasis cb; /** corresponding clusterbasis */
  uint      sons;
  /** number of sons for current clusterbasis */
  uint      m;
};

static void
uninit_interpolation_bem3d(paprxbem3d aprx)
{
  if (aprx->x_inter != NULL) {
    freemem(aprx->x_inter);
    aprx->x_inter = NULL;
  }
  aprx->m_inter = 0;
  aprx->k_inter = 0;
}

static void
uninit_green_bem3d(paprxbem3d aprx)
{
  if (aprx->t_green != NULL) {
    freemem(aprx->t_green);
    aprx->t_green = NULL;
  }
  if (aprx->w_green != NULL) {
    freemem(aprx->w_green);
    aprx->w_green = NULL;
  }
  aprx->m_green = 0;
  aprx->l_green = 0;
  aprx->delta_green = 0.0;
  aprx->ml_green = 0;
  aprx->k_green = 0;

  aprx->grc_green = NULL;
  aprx->gcc_green = NULL;
  aprx->grb_green = NULL;
  aprx->gcb_green = NULL;

}

static void
uninit_aca_bem3d(paprxbem3d aprx)
{
  aprx->accur_aca = 0.0;
}

static void
uninit_recompression_bem3d(paprxbem3d aprx)
{
  aprx->recomp = false;
  aprx->accur_recomp = 0.0;
  aprx->coarsen = false;
  aprx->accur_coarsen = 0.0;
  aprx->hiercomp = false;
  aprx->accur_hiercomp = 0.0;
  if (aprx->tm != NULL) {
    aprx->tm = NULL;
  }
}

/* ------------------------------------------------------------
 Constructors and destructors
 ------------------------------------------------------------ */

static    pgreencluster3d
new_greencluster3d(pccluster c)
{
  pgreencluster3d gc;
  uint      sons = c->sons;

  gc = (pgreencluster3d) allocmem(sizeof(greencluster3d));

  gc->V = new_amatrix(0, 0);
  gc->xi = NULL;
  gc->xihat = NULL;

  gc->t = c;
  gc->sons = sons;

  return gc;
}

static void
del_greencluster3d(pgreencluster3d gc)
{
  if (gc->xi != NULL) {
    freemem(gc->xi);
  }

  if (gc->xihat != NULL) {
    freemem(gc->xihat);
  }

  if (gc->V != NULL) {
    del_amatrix(gc->V);
  }

  freemem(gc);
}

static    pgreenclusterbasis3d
new_greenclusterbasis3d(pcclusterbasis cb)
{
  pgreenclusterbasis3d gcb;
  uint      sons = cb->sons;

  gcb = (pgreenclusterbasis3d) allocmem(sizeof(greenclusterbasis3d));

  gcb->xi = NULL;
  gcb->xihat = NULL;
  gcb->Qinv = NULL;

  gcb->cb = cb;
  gcb->sons = sons;
  gcb->m = 0;

  return gcb;
}

static void
del_greenclusterbasis3d(pgreenclusterbasis3d gcb)
{
  if (gcb->xi != NULL) {
    freemem(gcb->xi);
  }
  if (gcb->xihat != NULL) {
    freemem(gcb->xihat);
  }

  if (gcb->Qinv != NULL) {
    del_amatrix(gcb->Qinv);
    gcb->Qinv = NULL;
  }

  freemem(gcb);
}

static    plistnode
new_listnode(uint data, plistnode next)
{
  plistnode ln;

  ln = (plistnode) allocmem(sizeof(listnode));
  ln->data = data;
  ln->next = next;

  return ln;
}

static void
del_listnode(plistnode ln)
{
  plistnode ln1;

  if (ln) {
    ln1 = ln;
    while (ln1->next) {
      ln = ln1;
      ln1 = ln1->next;
      freemem(ln);
    }
    freemem(ln1);
  }

}

static    paprxbem3d
new_aprxbem3d()
{
  paprxbem3d aprx;

  aprx = (paprxbem3d) allocmem(sizeof(aprxbem3d));

  /* Interpolation */
  aprx->x_inter = NULL;
  aprx->m_inter = 0;
  aprx->k_inter = 0;

  /* Green */
  aprx->m_green = 0;
  aprx->l_green = 0;
  aprx->ml_green = 0;
  aprx->k_green = 0;
  aprx->delta_green = 0.0;
  aprx->t_green = NULL;
  aprx->w_green = NULL;
  aprx->quadpoints = NULL;
  aprx->grc_green = NULL;
  aprx->gcc_green = NULL;
  aprx->grb_green = NULL;
  aprx->gcb_green = NULL;

  /* ACA */
  aprx->accur_aca = 0.0;

  /* Recompression */
  aprx->recomp = false;
  aprx->accur_recomp = 0.0;
  aprx->coarsen = false;
  aprx->accur_coarsen = 0.0;
  aprx->hiercomp = false;
  aprx->accur_hiercomp = 0.0;
  aprx->tm = NULL;

  return aprx;
}

static void
del_aprxbem3d(paprxbem3d aprx)
{
  uninit_interpolation_bem3d(aprx);
  uninit_green_bem3d(aprx);
  uninit_aca_bem3d(aprx);
  uninit_recompression_bem3d(aprx);

  freemem(aprx);
}

static    pkernelbem3d
new_kernelbem3d()
{
  pkernelbem3d kernels;

  kernels = allocmem(sizeof(kernelbem3d));

  kernels->kernel_row = NULL;
  kernels->kernel_col = NULL;
  kernels->dnz_kernel_row = NULL;
  kernels->dnz_kernel_col = NULL;
  kernels->fundamental_row = NULL;
  kernels->fundamental_col = NULL;
  kernels->dnz_fundamental_row = NULL;
  kernels->dnz_fundamental_col = NULL;
  kernels->lagrange_row = NULL;
  kernels->lagrange_wave_row = NULL;
  kernels->lagrange_col = NULL;
  kernels->lagrange_wave_col = NULL;

  return kernels;
}

static void
del_kernelbem3d(pkernelbem3d kernels)
{
  freemem(kernels);
}

static    pparbem3d
new_parbem3d()
{
  pparbem3d par;

  par = (pparbem3d) allocmem(sizeof(parbem3d));

  par->rbn = NULL;
  par->cbn = NULL;
  par->drbn = NULL;
  par->dcbn = NULL;
  par->hn = NULL;
  par->h2n = NULL;
  par->dh2n = NULL;
  par->rwn = NULL;
  par->cwn = NULL;
  par->ron = NULL;
  par->con = NULL;
  par->leveln = NULL;
  par->grcn = NULL;
  par->grcnn = 0;
  par->gccn = NULL;
  par->gccnn = 0;
  par->grbn = NULL;
  par->grbnn = 0;
  par->gcbn = NULL;
  par->gcbnn = 0;

  return par;
}

static void
del_parbem3d(pparbem3d par)
{
  uint      n, i;

  if (par->hn != NULL) {
    freemem(par->hn);
  }

  if (par->h2n != NULL) {
    freemem(par->h2n);
  }

  if (par->rwn != NULL) {
    freemem(par->rwn);
  }

  if (par->cwn != NULL) {
    freemem(par->cwn);
  }

  if (par->ron != NULL) {
    freemem(par->ron);
  }

  if (par->con != NULL) {
    freemem(par->con);
  }

  if (par->leveln != NULL) {
    freemem(par->leveln);
  }

  /*
   * greencluster
   */

  n = par->grcnn;
  if (par->grcn != NULL && n != 0) {
    for (i = 0; i < n; ++i) {
      if (par->grcn[i] != NULL) {
	del_greencluster3d(par->grcn[i]);
      }
    }
    freemem(par->grcn);
  }

  n = par->gccnn;
  if (par->gccn != NULL && n != 0) {
    for (i = 0; i < n; ++i) {
      if (par->gccn[i] != NULL) {
	del_greencluster3d(par->gccn[i]);
      }
    }
    freemem(par->gccn);
  }

  /*
   * greenclusterbasis
   */

  n = par->grbnn;
  if (par->grbn != NULL && n != 0) {
    for (i = 0; i < n; ++i) {
      if (par->grbn[i] != NULL) {
	del_greenclusterbasis3d(par->grbn[i]);
      }
    }
    freemem(par->grbn);
  }

  n = par->gcbnn;
  if (par->gcbn != NULL && n != 0) {
    for (i = 0; i < n; ++i) {
      if (par->gcbn[i] != NULL) {
	del_greenclusterbasis3d(par->gcbn[i]);
      }
    }
    freemem(par->gcbn);
  }

  freemem(par);
}

pbem3d
new_bem3d(pcsurface3d gr, basisfunctionbem3d row_basis,
	  basisfunctionbem3d col_basis)
{
  pbem3d    bem;

  bem = (pbem3d) allocmem(sizeof(bem3d));

  bem->gr = gr;

  bem->mass = NULL;
  bem->v2t = NULL;
  bem->alpha = 0.0;

  bem->row_basis = row_basis;
  bem->col_basis = col_basis;

  bem->nearfield = NULL;
  bem->farfield_rk = NULL;
  bem->farfield_u = NULL;
  bem->farfield_wave_u = NULL;
  bem->leaf_row = NULL;
  bem->leaf_wave_row = NULL;
  bem->leaf_col = NULL;
  bem->leaf_wave_col = NULL;
  bem->transfer_row = NULL;
  bem->transfer_wave_row = NULL;
  bem->transfer_wave_wave_row = NULL;
  bem->transfer_col = NULL;
  bem->transfer_wave_col = NULL;
  bem->transfer_wave_wave_col = NULL;

  bem->aprx = new_aprxbem3d();
  bem->kernels = new_kernelbem3d();
  bem->par = new_parbem3d();


  if (row_basis == BASIS_LINEAR_BEM3D || col_basis == BASIS_LINEAR_BEM3D) {
    setup_vertex_to_triangle_map_bem3d(bem);
  }

  return bem;
}

void
del_bem3d(pbem3d bem)
{
  uint      i, n;

  if (bem->mass != NULL) {
    freemem(bem->mass);
  }

  if (bem->sq != NULL) {
    del_singquad2d(bem->sq);
  }

  del_aprxbem3d(bem->aprx);
  del_kernelbem3d(bem->kernels);
  del_parbem3d(bem->par);

  n = bem->gr->vertices;
  if (bem->v2t != NULL) {
    for (i = 0; i < n; ++i) {
      if (bem->v2t[i] != NULL) {
	del_listnode(bem->v2t[i]);
      }
    }
    freemem(bem->v2t);
  }

  freemem(bem);
}

pvert_list
new_vert_list(pvert_list next)
{
  pvert_list vl;

  vl = (pvert_list) allocmem(sizeof(vert_list));
  vl->v = 0;
  vl->next = next;

  return vl;
}

void
del_vert_list(pvert_list vl)
{
  pvert_list vl1;

  assert(vl != NULL);

  vl1 = vl;
  while (vl1) {
    vl = vl1->next;
    freemem(vl1);
    vl1 = vl;
  }
}

ptri_list
new_tri_list(ptri_list next)
{
  ptri_list tl;

  tl = (ptri_list) allocmem(sizeof(tri_list));
  tl->t = 0;
  tl->vl = NULL;
  tl->next = next;

  return tl;
}

void
del_tri_list(ptri_list tl)
{
  ptri_list tl1;

  assert(tl != NULL);

  tl1 = tl;
  while (tl1) {
    tl = tl1->next;
    del_vert_list(tl1->vl);
    freemem(tl1);
    tl1 = tl;
  }
}

void
setup_vertex_to_triangle_map_bem3d(pbem3d bem)
{
  pcsurface3d gr = bem->gr;
  const uint vertices = gr->vertices;
  const uint triangles = gr->triangles;
  uint(*tri)[3] = gr->t;

  plistnode *v2t = allocmem(vertices * sizeof(plistnode));
  uint      i, *t;

  /* TODO no sentinel node */
  for (i = 0; i < vertices; ++i) {
    v2t[i] = new_listnode(0, NULL);
  }

  for (i = 0; i < triangles; ++i) {
    t = tri[i];

    v2t[t[0]] = new_listnode(i, v2t[t[0]]);
    v2t[t[1]] = new_listnode(i, v2t[t[1]]);
    v2t[t[2]] = new_listnode(i, v2t[t[2]]);
  }

  bem->v2t = v2t;
}

/* ------------------------------------------------------------
 Methods to build clustertrees
 ------------------------------------------------------------ */

pclustergeometry
build_bem3d_const_clustergeometry(pcbem3d bem, uint ** idx)
{
  pcsurface3d gr = bem->gr;
  const     real(*x)[3] = (const real(*)[3]) gr->x;
  const     uint(*t)[3] = (const uint(*)[3]) gr->t;
  const real *g = (const real *) gr->g;
  uint      triangles = gr->triangles;

  pclustergeometry cg;
  uint      i;

  cg = new_clustergeometry(3, triangles);
  *idx = allocuint(triangles);

  for (i = 0; i < triangles; i++) {
    (*idx)[i] = i;

    /* Center of gravity as characteristic point */
    cg->x[i][0] = (x[t[i][0]][0] + x[t[i][1]][0] + x[t[i][2]][0]) / 3.0;
    cg->x[i][1] = (x[t[i][0]][1] + x[t[i][1]][1] + x[t[i][2]][1]) / 3.0;
    cg->x[i][2] = (x[t[i][0]][2] + x[t[i][1]][2] + x[t[i][2]][2]) / 3.0;

    /* Lower front left corner of bounding box */
    cg->smin[i][0] = REAL_MIN3(x[t[i][0]][0], x[t[i][1]][0], x[t[i][2]][0]);
    cg->smin[i][1] = REAL_MIN3(x[t[i][0]][1], x[t[i][1]][1], x[t[i][2]][1]);
    cg->smin[i][2] = REAL_MIN3(x[t[i][0]][2], x[t[i][1]][2], x[t[i][2]][2]);

    /* Upper back right corner of bounding box */
    cg->smax[i][0] = REAL_MAX3(x[t[i][0]][0], x[t[i][1]][0], x[t[i][2]][0]);
    cg->smax[i][1] = REAL_MAX3(x[t[i][0]][1], x[t[i][1]][1], x[t[i][2]][1]);
    cg->smax[i][2] = REAL_MAX3(x[t[i][0]][2], x[t[i][1]][2], x[t[i][2]][2]);

    cg->w[i] = g[i];
  }

  return cg;
}

pclustergeometry
build_bem3d_linear_clustergeometry(pcbem3d bem, uint ** idx)
{
  pcsurface3d gr = bem->gr;
  const     real(*x)[3] = (const real(*)[3]) gr->x;
  const     uint(*t)[3] = (const uint(*)[3]) gr->t;
  uint      triangles = gr->triangles;
  uint      vertices = gr->vertices;

  pclustergeometry cg;
  uint      tv0, tv1, tv2, i, j;
  real      tmin[3], tmax[3];

  cg = new_clustergeometry(3, vertices);
  *idx = allocuint(vertices);

  for (i = 0; i < vertices; ++i) {
    (*idx)[i] = i;

    /* Vertices as characteristic points */
    cg->x[i][0] = x[i][0];
    cg->x[i][1] = x[i][1];
    cg->x[i][2] = x[i][2];
    cg->smin[i][0] = x[i][0];
    cg->smin[i][1] = x[i][1];
    cg->smin[i][2] = x[i][2];
    cg->smax[i][0] = x[i][0];
    cg->smax[i][1] = x[i][1];
    cg->smax[i][2] = x[i][2];
  }

  for (i = 0; i < triangles; i++) {
    tv0 = t[i][0];
    tv1 = t[i][1];
    tv2 = t[i][2];

    /* Lower front left corner of bounding box for current triangle */
    tmin[0] = REAL_MIN3(x[tv0][0], x[tv1][0], x[tv2][0]);
    tmin[1] = REAL_MIN3(x[tv0][1], x[tv1][1], x[tv2][1]);
    tmin[2] = REAL_MIN3(x[tv0][2], x[tv1][2], x[tv2][2]);

    /* Upper back right corner of bounding box for current triangle */
    tmax[0] = REAL_MAX3(x[tv0][0], x[tv1][0], x[tv2][0]);
    tmax[1] = REAL_MAX3(x[tv0][1], x[tv1][1], x[tv2][1]);
    tmax[2] = REAL_MAX3(x[tv0][2], x[tv1][2], x[tv2][2]);

    /* update the bounding box for every vertex of the current triangle */
    for (j = 0; j < 3; ++j) {
      tv0 = t[i][j];

      /* Lower front left corner of bounding box */
      cg->smin[tv0][0] = REAL_MIN(cg->smin[tv0][0], tmin[0]);
      cg->smin[tv0][1] = REAL_MIN(cg->smin[tv0][1], tmin[1]);
      cg->smin[tv0][2] = REAL_MIN(cg->smin[tv0][2], tmin[2]);

      /* Upper back right corner of bounding box */
      cg->smax[tv0][0] = REAL_MAX(cg->smax[tv0][0], tmax[0]);
      cg->smax[tv0][1] = REAL_MAX(cg->smax[tv0][1], tmax[1]);
      cg->smax[tv0][2] = REAL_MAX(cg->smax[tv0][2], tmax[2]);
    }
  }

  return cg;
}

pclustergeometry
build_bem3d_clustergeometry(pcbem3d bem, uint ** idx,
			    basisfunctionbem3d basis)
{
  pclustergeometry cg;

  if (basis == BASIS_CONSTANT_BEM3D) {
    cg = build_bem3d_const_clustergeometry(bem, idx);
  }
  else {
    assert(basis == BASIS_LINEAR_BEM3D);
    cg = build_bem3d_linear_clustergeometry(bem, idx);
  }

  return cg;
}

pcluster
build_bem3d_cluster(pcbem3d bem, uint clf, basisfunctionbem3d basis)
{
  pclustergeometry cg;
  pcluster  c;
  uint     *idx;
  uint      n;

  if (basis == BASIS_CONSTANT_BEM3D) {
    cg = build_bem3d_const_clustergeometry(bem, &idx);
    n = bem->gr->triangles;
  }
  else {
    assert(basis == BASIS_LINEAR_BEM3D);
    cg = build_bem3d_linear_clustergeometry(bem, &idx);
    n = bem->gr->vertices;
  }

  c = build_adaptive_cluster(cg, n, idx, clf);

  del_clustergeometry(cg);

  return c;
}

/****************************************************
 * Nearfield integration routines
 ****************************************************/

#ifdef USE_SIMD
void
assemble_cc_simd_near_bem3d(const uint * ridx, const uint * cidx,
			    pcbem3d bem, bool ntrans, pamatrix N,
			    kernel_simd_func3d kernel)
{
  const pcsurface3d gr = bem->gr;
  const     real(*gr_x)[3] = (const real(*)[3]) gr->x;
  const     uint(*gr_t)[3] = (const uint(*)[3]) gr->t;
  const     real(*gr_n)[3] = (const real(*)[3]) gr->n;
  const real *gr_g = (const real *) gr->g;
  field    *aa = N->a;
  uint      rows = ntrans ? N->cols : N->rows;
  uint      cols = ntrans ? N->rows : N->cols;
  longindex ld = N->ld;

#ifdef USE_OPENMP
#pragma omp parallel if(!omp_in_parallel() && cols >= 256) num_threads(1 << max_pardepth)
  {
#endif
    vreal     vt[3][3], vs[3][3], nx[3], ny[3];
    const uint *tri_t, *tri_s;
    real     *xq, *yq, *wq;
    uint      tp[3], sp[3];
    real      factor, factor2, base;
    vreal     ct[3], cs[3], tx, sx, ty, sy, w, x[3], y[3], sum_r, sum_i,
      eval_r, eval_i, vcount, cmp;
    uint      q, nq, vnq, remainder, ss, tt, s, t, i, j;
#ifdef USE_TRIQUADPOINTS
    real     *tri_tx, *tri_ty, *tri_tz, *tri_sx, *tri_sy, *tri_sz;
    uint      c;
    uint      q2;
    uint      nq2 = bem->sq->n_single;
    uint      vnq2 = ROUNDUP(nq2, VREAL);
#endif

    vreal     c_one = vset1(1.0);

#ifdef USE_OPENMP
#pragma omp for
#endif
    for (s = 0; s < cols; ++s) {
      ss = (cidx == NULL ? s : cidx[s]);
      tri_s = gr_t[ss];
      factor = gr_g[ss] * bem->kernel_const;
      for (i = 0; i < 3; ++i) {
	ny[i] = vload1(gr_n[ss] + i);
      }
#ifdef USE_TRIQUADPOINTS
      tri_sx = bem->sq->tri_x + ss * vnq2;
      tri_sy = bem->sq->tri_y + ss * vnq2;
      tri_sz = bem->sq->tri_z + ss * vnq2;
#endif

      for (t = 0; t < rows; ++t) {
	tt = (ridx == NULL ? t : ridx[t]);
	tri_t = gr_t[tt];
	factor2 = factor * gr_g[tt];
	for (i = 0; i < 3; ++i) {
	  nx[i] = vload1(gr_n[tt] + i);
	}
#ifdef USE_TRIQUADPOINTS
	tri_tx = bem->sq->tri_x + tt * vnq2;
	tri_ty = bem->sq->tri_y + tt * vnq2;
	tri_tz = bem->sq->tri_z + tt * vnq2;
#endif

#ifdef USE_TRIQUADPOINTS
	c = select_quadrature_singquad2d(bem->sq, tri_t, tri_s, tp, sp, &xq,
					 &yq, &wq, &nq, &base);
	if (nq % VREAL || nq2 % VREAL) {
	  for (q = 0; q < VREAL; ++q) {
	    vcount[q] = q;
	  }
	}
#else
	(void) select_quadrature_singquad2d(bem->sq, tri_t, tri_s, tp, sp,
					    &xq, &yq, &wq, &nq, &base);
	if (nq % VREAL) {
	  for (q = 0; q < VREAL; ++q) {
	    vcount[q] = q;
	  }
	}
#endif
	vnq = ROUNDUP(nq, VREAL);

	wq += 9 * vnq;

	sum_r = vsetzero();
	sum_i = vsetzero();

#ifdef USE_TRIQUADPOINTS
	if (c == 0) {
	  remainder = vnq2 - VREAL;

	  for (q = 0; q < nq2; q++) {
	    x[0] = vload1(tri_tx + q);
	    x[1] = vload1(tri_ty + q);
	    x[2] = vload1(tri_tz + q);
	    for (q2 = 0; q2 < vnq2; q2 += VREAL) {
	      y[0] = vload(tri_sx + q2);
	      y[1] = vload(tri_sy + q2);
	      y[2] = vload(tri_sz + q2);
	      w = vloadu(wq + q2 + q * nq2);

	      kernel(x, y, nx, ny, (void *) bem, &eval_r, &eval_i);

	      if (nq2 % VREAL && q2 >= remainder) {
		cmp =
		  vcmplt(vadd(vcount, vset1((real) q2)), vset1((real) nq2));
		eval_r = vand(cmp, eval_r);
		eval_i = vand(cmp, eval_i);
	      }

	      sum_r = vfmadd(w, eval_r, sum_r);
	      sum_i = vfmadd(w, eval_i, sum_i);
	    }
	  }
	}
	else {
#endif
	  remainder = vnq - VREAL;

	  for (i = 0; i < 3; ++i) {	// x-, y-, z- component
	    for (j = 0; j < 3; ++j) {	// vertex A, B, C
	      vt[i][j] = vload1(gr_x[tri_t[tp[j]]] + i);
	      vs[i][j] = vload1(gr_x[tri_s[sp[j]]] + i);
	    }
	  }

	  for (q = 0; q < vnq; q += VREAL) {
	    tx = vload(xq + q);
	    sx = vload(xq + q + vnq);
	    ty = vload(yq + q);
	    sy = vload(yq + q + vnq);
	    w = vload(wq + q);

	    ct[0] = vsub(c_one, tx);
	    ct[1] = vsub(tx, sx);
	    ct[2] = sx;
	    cs[0] = vsub(c_one, ty);
	    cs[1] = vsub(ty, sy);
	    cs[2] = sy;

	    for (i = 0; i < 3; ++i) {
	      x[i] = vdot3(vt[i], ct);
	      y[i] = vdot3(vs[i], cs);
	    }

	    kernel(x, y, nx, ny, (void *) bem, &eval_r, &eval_i);

	    if (nq % VREAL && q >= remainder) {
	      cmp = vcmplt(vadd(vcount, vset1((real) q)), vset1((real) nq));
	      eval_r = vand(cmp, eval_r);
	      eval_i = vand(cmp, eval_i);
	    }

	    sum_r = vfmadd(w, eval_r, sum_r);
	    sum_i = vfmadd(w, eval_i, sum_i);
	  }
#ifdef USE_TRIQUADPOINTS
	}
#endif

	if (ntrans) {
	  aa[s + t * ld] = ((vreduce(sum_r) + base) * factor2)
	    - (vreduce(sum_i) * factor2) * I;
	}
	else {
	  aa[t + s * ld] = ((vreduce(sum_r) + base) * factor2)
	    + (vreduce(sum_i) * factor2) * I;
	}

	if (bem->alpha != 0.0 && tt == ss) {
	  if (ntrans) {
	    aa[t + t * ld] += 0.5 * CONJ(bem->alpha) * gr_g[tt];
	  }
	  else {
	    aa[t + t * ld] += 0.5 * bem->alpha * gr_g[tt];
	  }
	}
      }
    }
#ifdef USE_OPENMP
  }
#endif
}

#endif

void
assemble_cc_near_bem3d(const uint * ridx, const uint * cidx, pcbem3d bem,
		       bool ntrans, pamatrix N, kernel_func3d kernel)
{
  const pcsurface3d gr = bem->gr;
  const     real(*gr_x)[3] = (const real(*)[3]) gr->x;
  const     uint(*gr_t)[3] = (const uint(*)[3]) gr->t;
  const     real(*gr_n)[3] = (const real(*)[3]) gr->n;
  const real *gr_g = (const real *) gr->g;
  field    *aa = N->a;
  uint      rows = ntrans ? N->cols : N->rows;
  uint      cols = ntrans ? N->rows : N->cols;
  longindex ld = N->ld;

#ifdef USE_OPENMP
#pragma omp parallel if(!omp_in_parallel() && cols >= 256) num_threads(1 << max_pardepth)
  {
#endif
    const real *A_t, *B_t, *C_t, *A_s, *B_s, *C_s, *nx, *ny;
    const uint *tri_t, *tri_s;
    real     *xq, *yq, *wq;
    uint      tp[3], sp[3];
    real      Ax, Bx, Cx, Ay, By, Cy, tx, sx, ty, sy, x[3], y[3], factor,
      factor2, base;
    field     sum;
    uint      q, nq, vnq, ss, tt, s, t;

#ifdef USE_OPENMP
#pragma omp for
#endif
    for (s = 0; s < cols; ++s) {
      ss = (cidx == NULL ? s : cidx[s]);
      tri_s = gr_t[ss];
      factor = gr_g[ss] * bem->kernel_const;
      ny = gr_n[ss];

      for (t = 0; t < rows; ++t) {
	tt = (ridx == NULL ? t : ridx[t]);
	tri_t = gr_t[tt];
	factor2 = factor * gr_g[tt];
	nx = gr_n[tt];

	select_quadrature_singquad2d(bem->sq, tri_t, tri_s, tp, sp, &xq, &yq,
				     &wq, &nq, &base);
	vnq = ROUNDUP(nq, VREAL);
	wq += 9 * vnq;

	A_t = gr_x[tri_t[tp[0]]];
	B_t = gr_x[tri_t[tp[1]]];
	C_t = gr_x[tri_t[tp[2]]];
	A_s = gr_x[tri_s[sp[0]]];
	B_s = gr_x[tri_s[sp[1]]];
	C_s = gr_x[tri_s[sp[2]]];

	sum = base;

	for (q = 0; q < nq; ++q) {
	  tx = xq[q];
	  sx = xq[q + vnq];
	  ty = yq[q];
	  sy = yq[q + vnq];
	  Ax = 1.0 - tx;
	  Bx = tx - sx;
	  Cx = sx;
	  Ay = 1.0 - ty;
	  By = ty - sy;
	  Cy = sy;

	  x[0] = A_t[0] * Ax + B_t[0] * Bx + C_t[0] * Cx;
	  x[1] = A_t[1] * Ax + B_t[1] * Bx + C_t[1] * Cx;
	  x[2] = A_t[2] * Ax + B_t[2] * Bx + C_t[2] * Cx;
	  y[0] = A_s[0] * Ay + B_s[0] * By + C_s[0] * Cy;
	  y[1] = A_s[1] * Ay + B_s[1] * By + C_s[1] * Cy;
	  y[2] = A_s[2] * Ay + B_s[2] * By + C_s[2] * Cy;

	  sum += wq[q] * kernel(x, y, nx, ny, (void *) bem);
	}
	if (ntrans) {
	  aa[s + t * ld] = CONJ(sum) * factor2;
	}
	else {
	  aa[t + s * ld] = sum * factor2;
	}

	if (bem->alpha != 0.0 && tt == ss) {
	  if (ntrans) {
	    aa[t + t * ld] += 0.5 * CONJ(bem->alpha) * gr_g[tt];
	  }
	  else {
	    aa[t + t * ld] += 0.5 * bem->alpha * gr_g[tt];
	  }
	}
      }
    }
#ifdef USE_OPENMP
  }
#endif
}

#ifdef USE_SIMD
void
assemble_cc_simd_far_bem3d(const uint * ridx, const uint * cidx,
			   pcbem3d bem, bool ntrans, pamatrix N,
			   kernel_simd_func3d kernel)
{
  const pcsurface3d gr = bem->gr;
  const     real(*gr_n)[3] = (const real(*)[3]) gr->n;
  const real *gr_g = (const real *) gr->g;
  field    *aa = N->a;
  uint      rows = ntrans ? N->cols : N->rows;
  uint      cols = ntrans ? N->rows : N->cols;
  longindex ld = N->ld;

#ifdef USE_OPENMP
#pragma omp parallel if(!omp_in_parallel() && cols >= 256) num_threads(1 << max_pardepth)
  {
#endif
    vreal     nx[3], ny[3];
    real     *wq;
    real      factor, factor2, base;
    vreal     w, x[3], y[3], sum_r, sum_i, eval_r, eval_i, vcount, cmp;
    uint      q, nq, vnq, remainder, ss, tt, s, t, i;
#ifdef USE_TRIQUADPOINTS
    real     *tri_tx, *tri_ty, *tri_tz, *tri_sx, *tri_sy, *tri_sz;
    uint      q2;
    uint      nq2 = bem->sq->n_single;
    uint      vnq2 = ROUNDUP(nq2, VREAL);
#else
    const     real(*gr_x)[3] = (const real(*)[3]) gr->x;
    const     uint(*gr_t)[3] = (const uint(*)[3]) gr->t;

    uint      j;
    vreal     c_one = vset1(1.0);
    vreal     vt[3][3], vs[3][3];
    vreal     ct[3], cs[3], tx, sx, ty, sy;
    real     *xq = bem->sq->x_dist;
    real     *yq = bem->sq->y_dist;
    const uint *tri_t, *tri_s;
#endif

    nq = bem->sq->n_dist;
    vnq = ROUNDUP(nq, VREAL);
    wq = bem->sq->w_dist + 9 * vnq;
    base = bem->sq->base_dist;

#ifdef USE_OPENMP
#pragma omp for
#endif
    for (s = 0; s < cols; ++s) {
      ss = (cidx == NULL ? s : cidx[s]);

      factor = gr_g[ss] * bem->kernel_const;
      for (i = 0; i < 3; ++i) {
	ny[i] = vload1(gr_n[ss] + i);
      }
#ifdef USE_TRIQUADPOINTS
      tri_sx = bem->sq->tri_x + ss * vnq2;
      tri_sy = bem->sq->tri_y + ss * vnq2;
      tri_sz = bem->sq->tri_z + ss * vnq2;
#else
      tri_s = gr_t[ss];
#endif

      for (t = 0; t < rows; ++t) {
	tt = (ridx == NULL ? t : ridx[t]);

	factor2 = factor * gr_g[tt];
	for (i = 0; i < 3; ++i) {
	  nx[i] = vload1(gr_n[tt] + i);
	}
#ifdef USE_TRIQUADPOINTS
	tri_tx = bem->sq->tri_x + tt * vnq2;
	tri_ty = bem->sq->tri_y + tt * vnq2;
	tri_tz = bem->sq->tri_z + tt * vnq2;
#else
	tri_t = gr_t[tt];
#endif

#ifdef USE_TRIQUADPOINTS
	if (nq % VREAL || nq2 % VREAL) {
	  for (q = 0; q < VREAL; ++q) {
	    vcount[q] = q;
	  }
	}
#else
	if (nq % VREAL) {
	  for (q = 0; q < VREAL; ++q) {
	    vcount[q] = q;
	  }
	}
#endif

	sum_r = vsetzero();
	sum_i = vsetzero();

#ifdef USE_TRIQUADPOINTS
	remainder = vnq2 - VREAL;

	for (q = 0; q < nq2; q++) {
	  x[0] = vload1(tri_tx + q);
	  x[1] = vload1(tri_ty + q);
	  x[2] = vload1(tri_tz + q);
	  for (q2 = 0; q2 < vnq2; q2 += VREAL) {
	    y[0] = vload(tri_sx + q2);
	    y[1] = vload(tri_sy + q2);
	    y[2] = vload(tri_sz + q2);
	    w = vloadu(wq + q2 + q * nq2);

	    kernel(x, y, nx, ny, (void *) bem, &eval_r, &eval_i);

	    if (nq2 % VREAL && q2 >= remainder) {
	      cmp = vcmplt(vadd(vcount, vset1((real) q2)), vset1((real) nq2));
	      eval_r = vand(cmp, eval_r);
	      eval_i = vand(cmp, eval_i);
	    }

	    sum_r = vfmadd(w, eval_r, sum_r);
	    sum_i = vfmadd(w, eval_i, sum_i);
	  }
	}
#else
	remainder = vnq - VREAL;

	for (i = 0; i < 3; ++i) {	// x-, y-, z- component
	  for (j = 0; j < 3; ++j) {	// vertex A, B, C
	    vt[i][j] = vload1(gr_x[tri_t[j]] + i);
	    vs[i][j] = vload1(gr_x[tri_s[j]] + i);
	  }
	}

	for (q = 0; q < vnq; q += VREAL) {
	  tx = vload(xq + q);
	  sx = vload(xq + q + vnq);
	  ty = vload(yq + q);
	  sy = vload(yq + q + vnq);
	  w = vload(wq + q);

	  ct[0] = vsub(c_one, tx);
	  ct[1] = vsub(tx, sx);
	  ct[2] = sx;
	  cs[0] = vsub(c_one, ty);
	  cs[1] = vsub(ty, sy);
	  cs[2] = sy;

	  for (i = 0; i < 3; ++i) {
	    x[i] = vdot3(vt[i], ct);
	    y[i] = vdot3(vs[i], cs);
	  }

	  kernel(x, y, nx, ny, (void *) bem, &eval_r, &eval_i);

	  if (nq % VREAL && q >= remainder) {
	    cmp = vcmplt(vadd(vcount, vset1((real) q)), vset1((real) nq));
	    eval_r = vand(cmp, eval_r);
	    eval_i = vand(cmp, eval_i);
	  }

	  sum_r = vfmadd(w, eval_r, sum_r);
	  sum_i = vfmadd(w, eval_i, sum_i);
	}
#endif

	if (ntrans) {
	  aa[s + t * ld] = ((vreduce(sum_r) + base) * factor2)
	    - (vreduce(sum_i) * factor2) * I;
	}
	else {
	  aa[t + s * ld] = ((vreduce(sum_r) + base) * factor2)
	    + (vreduce(sum_i) * factor2) * I;
	}

	if (bem->alpha != 0.0 && tt == ss) {
	  if (ntrans) {
	    aa[t + t * ld] += 0.5 * CONJ(bem->alpha) * gr_g[tt];
	  }
	  else {
	    aa[t + t * ld] += 0.5 * bem->alpha * gr_g[tt];
	  }
	}
      }
    }
#ifdef USE_OPENMP
  }
#endif
}

#endif

void
assemble_cc_far_bem3d(const uint * ridx, const uint * cidx, pcbem3d bem,
		      bool ntrans, pamatrix N, kernel_func3d kernel)
{
  const pcsurface3d gr = bem->gr;
  const     real(*gr_x)[3] = (const real(*)[3]) gr->x;
  const     uint(*gr_t)[3] = (const uint(*)[3]) gr->t;
  const     real(*gr_n)[3] = (const real(*)[3]) gr->n;
  const real *gr_g = (const real *) gr->g;
  field    *aa = N->a;
  uint      rows = ntrans ? N->cols : N->rows;
  uint      cols = ntrans ? N->rows : N->cols;
  longindex ld = N->ld;

#ifdef USE_OPENMP
#pragma omp parallel if(!omp_in_parallel() && cols >= 256) num_threads(1 << max_pardepth)
  {
#endif
    const real *A_t, *B_t, *C_t, *A_s, *B_s, *C_s, *nx, *ny;
    const uint *tri_t, *tri_s;
    real     *xq, *yq, *wq;
    real      Ax, Bx, Cx, Ay, By, Cy, tx, sx, ty, sy, x[3], y[3], factor,
      factor2;
    field     sum;
    uint      q, nq, vnq, ss, tt, s, t;

    xq = bem->sq->x_dist;
    yq = bem->sq->y_dist;
    nq = bem->sq->n_dist;
    vnq = ROUNDUP(nq, VREAL);
    wq = bem->sq->w_dist + 9 * vnq;

#ifdef USE_OPENMP
#pragma omp for
#endif
    for (s = 0; s < cols; ++s) {
      ss = (cidx == NULL ? s : cidx[s]);
      tri_s = gr_t[ss];
      factor = gr_g[ss] * bem->kernel_const;

      A_s = gr_x[tri_s[0]];
      B_s = gr_x[tri_s[1]];
      C_s = gr_x[tri_s[2]];
      ny = gr_n[ss];

      for (t = 0; t < rows; ++t) {
	tt = (ridx == NULL ? t : ridx[t]);
	tri_t = gr_t[tt];
	factor2 = factor * gr_g[tt];

	A_t = gr_x[tri_t[0]];
	B_t = gr_x[tri_t[1]];
	C_t = gr_x[tri_t[2]];
	nx = gr_n[tt];

	sum = bem->sq->base_dist;

	for (q = 0; q < nq; ++q) {
	  tx = xq[q];
	  sx = xq[q + vnq];
	  ty = yq[q];
	  sy = yq[q + vnq];
	  Ax = 1.0 - tx;
	  Bx = tx - sx;
	  Cx = sx;
	  Ay = 1.0 - ty;
	  By = ty - sy;
	  Cy = sy;

	  x[0] = A_t[0] * Ax + B_t[0] * Bx + C_t[0] * Cx;
	  x[1] = A_t[1] * Ax + B_t[1] * Bx + C_t[1] * Cx;
	  x[2] = A_t[2] * Ax + B_t[2] * Bx + C_t[2] * Cx;

	  y[0] = A_s[0] * Ay + B_s[0] * By + C_s[0] * Cy;
	  y[1] = A_s[1] * Ay + B_s[1] * By + C_s[1] * Cy;
	  y[2] = A_s[2] * Ay + B_s[2] * By + C_s[2] * Cy;

	  sum += wq[q] * kernel(x, y, nx, ny, (void *) bem);
	}
	if (ntrans) {
	  aa[s + t * ld] = CONJ(sum) * factor2;
	}
	else {
	  aa[t + s * ld] = sum * factor2;
	}
      }
    }
#ifdef USE_OPENMP
  }
#endif
}

#ifdef USE_SIMD
void
assemble_cl_simd_near_bem3d(const uint * ridx, const uint * cidx,
			    pcbem3d bem, bool ntrans, pamatrix N,
			    kernel_simd_func3d kernel)
{
  const pcsurface3d gr = bem->gr;
  const     real(*gr_x)[3] = (const real(*)[3]) gr->x;
  const     uint(*gr_t)[3] = (const uint(*)[3]) gr->t;
  const preal gr_g = (const preal) gr->g;
  const     real(*gr_n)[3] = (const real(*)[3]) gr->n;
  const uint triangles = gr->triangles;
  plistnode *v2t = bem->v2t;

  field    *aa = N->a;
  uint      rows = ntrans ? N->cols : N->rows;
  uint      cols = ntrans ? N->rows : N->cols;
  longindex ld = N->ld;

  vreal     vt[3][3], vs[3][3], nx[3], ny[3];
  real     *quad_r, *quad_i;
  ptri_list tl, tl1;
  pvert_list vl;
  const uint *tri_t, *tri_s;
  plistnode v;
  real     *xq, *yq, *wq, *mass;
  uint      tp[3], sp[3], tri_tp[3], tri_sp[3];
  vreal     ct[3], cs[3], tx, sx, ty, sy, w, x[3], y[3], sum_r, sum_i, eval_r,
    eval_i, vcount, cmp;
  real      base, factor, factor2;
  uint      i, j, t, s, q, nq, vnq, remainder, cj;
  uint      ii, jj, tt, ss, vv;
#ifdef USE_TRIQUADPOINTS
  real     *tri_tx, *tri_ty, *tri_tz, *tri_sx, *tri_sy, *tri_sz;
  uint      c;
  uint      q2;
  uint      nq2 = bem->sq->n_single;
  uint      vnq2 = ROUNDUP(nq2, VREAL);
#endif

  clear_amatrix(N);

  quad_r = allocreal(ROUNDUP(bem->sq->nmax, VREAL));
  quad_i = allocreal(ROUNDUP(bem->sq->nmax, VREAL));
  vreal     c_one = vset1(1.0);

  tl = NULL;

  cj = 0;
  for (i = 0; i < cols; ++i) {
    ii = (cidx == NULL ? i : cidx[i]);
    for (v = v2t[ii], vv = v->data; v->next != NULL;
	 v = v->next, vv = v->data) {

      tl1 = tl;
      while (tl1 && tl1->t != vv) {
	tl1 = tl1->next;
      }

      if (tl1 == NULL) {
	tl1 = tl = new_tri_list(tl);
	tl->t = vv;
	cj++;
      }

      tl1->vl = new_vert_list(tl1->vl);
      tl1->vl->v = i;
    }
  }

  for (s = 0, tl1 = tl; s < cj; s++, tl1 = tl1->next) {
    ss = tl1->t;
    assert(ss < triangles);
    factor = gr_g[ss] * bem->kernel_const;
    tri_s = gr_t[ss];
    for (i = 0; i < 3; ++i) {
      ny[i] = vload1(gr_n[ss] + i);
    }
#ifdef USE_TRIQUADPOINTS
    tri_sx = bem->sq->tri_x + ss * vnq2;
    tri_sy = bem->sq->tri_y + ss * vnq2;
    tri_sz = bem->sq->tri_z + ss * vnq2;
#endif
    for (t = 0; t < rows; t++) {
      tt = (ridx == NULL ? t : ridx[t]);
      assert(tt < triangles);
      factor2 = factor * gr_g[tt];
      tri_t = gr_t[tt];
      for (i = 0; i < 3; ++i) {
	nx[i] = vload1(gr_n[tt] + i);
      }
#ifdef USE_TRIQUADPOINTS
      tri_tx = bem->sq->tri_x + tt * vnq2;
      tri_ty = bem->sq->tri_y + tt * vnq2;
      tri_tz = bem->sq->tri_z + tt * vnq2;
#endif

#ifdef USE_TRIQUADPOINTS
      c =
	select_quadrature_singquad2d(bem->sq, tri_t, tri_s, tp, sp, &xq, &yq,
				     &wq, &nq, &base);
      if (nq % VREAL || nq2 % VREAL) {
	for (q = 0; q < VREAL; ++q) {
	  vcount[q] = q;
	}
      }

#else
      (void) select_quadrature_singquad2d(bem->sq, tri_t, tri_s, tp, sp, &xq,
					  &yq, &wq, &nq, &base);
      if (nq % VREAL) {
	for (q = 0; q < VREAL; ++q) {
	  vcount[q] = q;
	}
      }

#endif
      vnq = ROUNDUP(nq, VREAL);

#ifdef USE_TRIQUADPOINTS
      if (c == 0) {
	remainder = vnq2 - VREAL;

	for (i = 0; i < 3; ++i) {
	  tri_tp[i] = tri_t[i];
	  tri_sp[i] = tri_s[i];
	}

	for (q = 0; q < nq2; q++) {
	  x[0] = vload1(tri_tx + q);
	  x[1] = vload1(tri_ty + q);
	  x[2] = vload1(tri_tz + q);
	  for (q2 = 0; q2 < vnq2; q2 += VREAL) {
	    y[0] = vload(tri_sx + q2);
	    y[1] = vload(tri_sy + q2);
	    y[2] = vload(tri_sz + q2);

	    kernel(x, y, nx, ny, (void *) bem, &eval_r, &eval_i);

	    if (nq2 % VREAL && q2 >= remainder) {
	      cmp = vcmplt(vadd(vcount, vset1((real) q2)), vset1((real) nq2));
	      eval_r = vand(cmp, eval_r);
	      eval_i = vand(cmp, eval_i);
	    }

	    vstoreu(quad_r + q2 + q * nq2, eval_r);
	    vstoreu(quad_i + q2 + q * nq2, eval_i);
	  }
	}
      }
      else {
#endif
	remainder = vnq - VREAL;

	for (i = 0; i < 3; ++i) {
	  tri_tp[i] = tri_t[tp[i]];
	  tri_sp[i] = tri_s[sp[i]];
	}

	for (i = 0; i < 3; ++i) {	// x-, y-, z- component
	  for (j = 0; j < 3; ++j) {	// vertex A, B, C
	    vt[i][j] = vload1(gr_x[tri_tp[j]] + i);
	    vs[i][j] = vload1(gr_x[tri_sp[j]] + i);
	  }
	}

	for (q = 0; q < vnq; q += VREAL) {
	  tx = vload(xq + q);
	  sx = vload(xq + q + vnq);
	  ty = vload(yq + q);
	  sy = vload(yq + q + vnq);

	  ct[0] = vsub(c_one, tx);
	  ct[1] = vsub(tx, sx);
	  ct[2] = sx;
	  cs[0] = vsub(c_one, ty);
	  cs[1] = vsub(ty, sy);
	  cs[2] = sy;

	  for (i = 0; i < 3; ++i) {
	    x[i] = vdot3(vt[i], ct);
	    y[i] = vdot3(vs[i], cs);
	  }

	  kernel(x, y, nx, ny, (void *) bem, &eval_r, &eval_i);

	  if (nq % VREAL && q >= remainder) {
	    cmp = vcmplt(vadd(vcount, vset1((real) q)), vset1((real) nq));
	    eval_r = vand(cmp, eval_r);
	    eval_i = vand(cmp, eval_i);
	  }

	  vstore(quad_r + q, eval_r);
	  vstore(quad_i + q, eval_i);
	}
#ifdef USE_TRIQUADPOINTS
      }
#endif
      vl = tl1->vl;
      while (vl) {
	j = vl->v;
	if (j < cols) {
	  jj = cidx == NULL ? j : cidx[j];
	  for (i = 0; i < 3; ++i) {
	    if (jj == tri_sp[i]) {
	      sum_r = vsetzero();
	      sum_i = vsetzero();

	      for (q = 0; q < vnq; q += VREAL) {
		eval_r = vload(quad_r + q);
		eval_i = vload(quad_i + q);
		w = vload(wq + q);
		sum_r = vfmadd(w, eval_r, sum_r);
		sum_i = vfmadd(w, eval_i, sum_i);
	      }

	      if (ntrans) {
		aa[j + t * ld] += ((vreduce(sum_r) + base) * factor2)
		  - (vreduce(sum_i) * factor2) * I;
	      }
	      else {
		aa[t + j * ld] += ((vreduce(sum_r) + base) * factor2)
		  + (vreduce(sum_i) * factor2) * I;
	      }
	    }
	    wq += vnq;
	  }
	  wq -= 3 * vnq;
	}
	vl = vl->next;
      }

      if (bem->alpha != 0.0 && tt == ss) {

	for (i = 0; i < 3; ++i) {
	  tri_sp[i] = tri_s[i];
	}

	mass = bem->mass;
	factor2 = bem->alpha * gr_g[tt];

	vl = tl1->vl;
	while (vl) {
	  j = vl->v;
	  if (j < cols) {
	    jj = cidx == NULL ? j : cidx[j];
	    for (i = 0; i < 3; ++i) {
	      if (jj == tri_sp[i]) {
		if (ntrans) {
		  aa[j + t * ld] += factor2 * *mass;
		}
		else {
		  aa[t + j * ld] += factor2 * *mass;
		}
	      }
	      mass++;
	    }
	    mass = bem->mass;
	  }
	  vl = vl->next;
	}
      }
    }
  }
}
#endif

void
assemble_cl_near_bem3d(const uint * ridx, const uint * cidx, pcbem3d bem,
		       bool ntrans, pamatrix N, kernel_func3d kernel)
{
  const pcsurface3d gr = bem->gr;
  const     real(*gr_x)[3] = (const real(*)[3]) gr->x;
  const     uint(*gr_t)[3] = (const uint(*)[3]) gr->t;
  const preal gr_g = (const preal) gr->g;
  const     real(*gr_n)[3] = (const real(*)[3]) gr->n;
  const uint triangles = gr->triangles;
  plistnode *v2t = bem->v2t;

  field    *aa = N->a;
  uint      rows = ntrans ? N->cols : N->rows;
  uint      cols = ntrans ? N->rows : N->cols;
  longindex ld = N->ld;

  field    *quad;
  ptri_list tl, tl1;
  pvert_list vl;
  const real *A_t, *B_t, *C_t, *A_s, *B_s, *C_s, *ns, *nt;
  const uint *tri_t, *tri_s;
  plistnode v;
  real     *xq, *yq, *wq, *mass;
  uint      tp[3], sp[3], tri_tp[3], tri_sp[3];
  real      Ax, Bx, Cx, Ay, By, Cy, tx, sx, ty, sy, x[3], y[3], factor,
    factor2, base;
  field     res;
  uint      i, j, t, s, q, nq, vnq, cj;
  uint      ii, jj, tt, ss, vv;

  clear_amatrix(N);

  quad = allocfield(bem->sq->nmax);

  tl = NULL;

  cj = 0;
  for (i = 0; i < cols; ++i) {
    ii = (cidx == NULL ? i : cidx[i]);
    for (v = v2t[ii], vv = v->data; v->next != NULL;
	 v = v->next, vv = v->data) {

      tl1 = tl;
      while (tl1 && tl1->t != vv) {
	tl1 = tl1->next;
      }

      if (tl1 == NULL) {
	tl1 = tl = new_tri_list(tl);
	tl->t = vv;
	cj++;
      }

      tl1->vl = new_vert_list(tl1->vl);
      tl1->vl->v = i;
    }
  }

  for (s = 0, tl1 = tl; s < cj; s++, tl1 = tl1->next) {
    ss = tl1->t;
    assert(ss < triangles);
    factor = gr_g[ss] * bem->kernel_const;
    tri_s = gr_t[ss];
    ns = gr_n[ss];
    for (t = 0; t < rows; t++) {
      tt = (ridx == NULL ? t : ridx[t]);
      assert(tt < triangles);
      tri_t = gr_t[tt];
      nt = gr_n[tt];

      factor2 = factor * gr_g[tt];

      select_quadrature_singquad2d(bem->sq, tri_t, tri_s, tp, sp, &xq, &yq,
				   &wq, &nq, &base);
      vnq = ROUNDUP(nq, VREAL);

      for (i = 0; i < 3; ++i) {
	tri_tp[i] = tri_t[tp[i]];
	tri_sp[i] = tri_s[sp[i]];
      }

      A_t = gr_x[tri_tp[0]];
      B_t = gr_x[tri_tp[1]];
      C_t = gr_x[tri_tp[2]];
      A_s = gr_x[tri_sp[0]];
      B_s = gr_x[tri_sp[1]];
      C_s = gr_x[tri_sp[2]];

      for (q = 0; q < nq; ++q) {
	tx = xq[q];
	sx = xq[q + vnq];
	ty = yq[q];
	sy = yq[q + vnq];
	Ax = 1.0 - tx;
	Bx = tx - sx;
	Cx = sx;
	Ay = 1.0 - ty;
	By = ty - sy;
	Cy = sy;

	x[0] = A_t[0] * Ax + B_t[0] * Bx + C_t[0] * Cx;
	x[1] = A_t[1] * Ax + B_t[1] * Bx + C_t[1] * Cx;
	x[2] = A_t[2] * Ax + B_t[2] * Bx + C_t[2] * Cx;

	y[0] = A_s[0] * Ay + B_s[0] * By + C_s[0] * Cy;
	y[1] = A_s[1] * Ay + B_s[1] * By + C_s[1] * Cy;
	y[2] = A_s[2] * Ay + B_s[2] * By + C_s[2] * Cy;

	quad[q] = kernel(x, y, nt, ns, (void *) bem);
      }

      vl = tl1->vl;
      while (vl) {
	j = vl->v;
	if (j < cols) {
	  jj = cidx == NULL ? j : cidx[j];
	  for (i = 0; i < 3; ++i) {
	    if (jj == tri_sp[i]) {
	      res = base;

	      for (q = 0; q < nq; ++q) {
		res += wq[q] * quad[q];
	      }

	      if (ntrans) {
		aa[j + t * ld] += res * factor2;
	      }
	      else {
		aa[t + j * ld] += res * factor2;
	      }
	    }
	    wq += vnq;
	  }
	  wq -= 3 * vnq;
	}
	vl = vl->next;
      }

      if (bem->alpha != 0.0 && tt == ss) {

	for (i = 0; i < 3; ++i) {
	  tri_sp[i] = tri_s[i];
	}

	mass = bem->mass;
	factor2 = bem->alpha * gr_g[tt];

	vl = tl1->vl;
	while (vl) {
	  j = vl->v;
	  if (j < cols) {
	    jj = cidx == NULL ? j : cidx[j];
	    for (i = 0; i < 3; ++i) {
	      if (jj == tri_sp[i]) {
		if (ntrans) {
		  aa[j + t * ld] += factor2 * *mass;
		}
		else {
		  aa[t + j * ld] += factor2 * *mass;
		}
	      }
	      mass++;
	    }
	    mass = bem->mass;
	  }
	  vl = vl->next;
	}
      }
    }
  }

  del_tri_list(tl);
  freemem(quad);
}

#ifdef USE_SIMD
void
assemble_cl_simd_far_bem3d(const uint * ridx, const uint * cidx,
			   pcbem3d bem, bool ntrans, pamatrix N,
			   kernel_simd_func3d kernel)
{
  //TODO implement nice 'far' version
  assemble_cl_simd_near_bem3d(ridx, cidx, bem, ntrans, N, kernel);
}
#endif

void
assemble_cl_far_bem3d(const uint * ridx, const uint * cidx, pcbem3d bem,
		      bool ntrans, pamatrix N, kernel_func3d kernel)
{
  const pcsurface3d gr = bem->gr;
  const     real(*gr_x)[3] = (const real(*)[3]) gr->x;
  const     uint(*gr_t)[3] = (const uint(*)[3]) gr->t;
  const preal gr_g = (const preal) gr->g;
  const     real(*gr_n)[3] = (const real(*)[3]) gr->n;
  const uint triangles = gr->triangles;
  plistnode *v2t = bem->v2t;
  field    *quad;
  field    *aa = N->a;
  uint      rows = ntrans ? N->cols : N->rows;
  uint      cols = ntrans ? N->rows : N->cols;
  longindex ld = N->ld;

  ptri_list tl, tl1;
  pvert_list vl;
  const real *A_t, *B_t, *C_t, *A_s, *B_s, *C_s, *ns, *nt;
  const uint *tri_t, *tri_s;
  plistnode v;
  real     *xq, *yq, *wq;
  real      Ax, Bx, Cx, Ay, By, Cy, tx, sx, ty, sy, x[3], y[3], factor,
    factor2, base;
  field     res;
  uint      i, j, t, s, q, nq, vnq, cj;
  uint      ii, jj, tt, ss, vv;

  clear_amatrix(N);

  quad = allocfield(bem->sq->nmax);

  tl = NULL;

  xq = bem->sq->x_dist;
  yq = bem->sq->y_dist;
  wq = bem->sq->w_dist;
  nq = bem->sq->n_dist;
  vnq = ROUNDUP(nq, VREAL);
  base = bem->sq->base_dist;

  cj = 0;
  for (i = 0; i < cols; ++i) {
    ii = (cidx == NULL ? i : cidx[i]);
    for (v = v2t[ii], vv = v->data; v->next != NULL;
	 v = v->next, vv = v->data) {

      tl1 = tl;
      while (tl1 && tl1->t != vv) {
	tl1 = tl1->next;
      }

      if (tl1 == NULL) {
	tl1 = tl = new_tri_list(tl);
	tl->t = vv;
	cj++;
      }

      tl1->vl = new_vert_list(tl1->vl);
      tl1->vl->v = i;
    }
  }

  for (s = 0, tl1 = tl; s < cj; s++, tl1 = tl1->next) {
    ss = tl1->t;
    assert(ss < triangles);
    factor = gr_g[ss] * bem->kernel_const;
    tri_s = gr_t[ss];
    ns = gr_n[ss];

    A_s = gr_x[tri_s[0]];
    B_s = gr_x[tri_s[1]];
    C_s = gr_x[tri_s[2]];

    for (t = 0; t < rows; t++) {
      tt = (ridx == NULL ? t : ridx[t]);
      assert(tt < triangles);
      tri_t = gr_t[tt];
      factor2 = factor * gr_g[tt];
      nt = gr_n[tt];

      A_t = gr_x[tri_t[0]];
      B_t = gr_x[tri_t[1]];
      C_t = gr_x[tri_t[2]];

      for (q = 0; q < nq; ++q) {
	tx = xq[q];
	sx = xq[q + vnq];
	ty = yq[q];
	sy = yq[q + vnq];
	Ax = 1.0 - tx;
	Bx = tx - sx;
	Cx = sx;
	Ay = 1.0 - ty;
	By = ty - sy;
	Cy = sy;

	x[0] = A_t[0] * Ax + B_t[0] * Bx + C_t[0] * Cx;
	x[1] = A_t[1] * Ax + B_t[1] * Bx + C_t[1] * Cx;
	x[2] = A_t[2] * Ax + B_t[2] * Bx + C_t[2] * Cx;

	y[0] = A_s[0] * Ay + B_s[0] * By + C_s[0] * Cy;
	y[1] = A_s[1] * Ay + B_s[1] * By + C_s[1] * Cy;
	y[2] = A_s[2] * Ay + B_s[2] * By + C_s[2] * Cy;

	quad[q] = kernel(x, y, nt, ns, (void *) bem);
      }

      vl = tl1->vl;
      while (vl) {
	j = vl->v;
	if (j < cols) {
	  jj = cidx == NULL ? j : cidx[j];
	  for (i = 0; i < 3; ++i) {
	    if (jj == tri_s[i]) {
	      res = base;

	      for (q = 0; q < nq; ++q) {
		res += wq[q] * quad[q];
	      }

	      if (ntrans) {
		aa[j + t * ld] += res * factor2;
	      }
	      else {
		aa[t + j * ld] += res * factor2;
	      }
	    }
	    wq += vnq;
	  }
	  wq -= 3 * vnq;
	}
	vl = vl->next;
      }
    }
  }

  del_tri_list(tl);
  freemem(quad);
}

#ifdef USE_SIMD
void
assemble_lc_simd_near_bem3d(const uint * ridx, const uint * cidx,
			    pcbem3d bem, bool ntrans, pamatrix N,
			    kernel_simd_func3d kernel)
{
  const pcsurface3d gr = bem->gr;
  const     real(*gr_x)[3] = (const real(*)[3]) gr->x;
  const     uint(*gr_t)[3] = (const uint(*)[3]) gr->t;
  const preal gr_g = (const preal) gr->g;
  const     real(*gr_n)[3] = (const real(*)[3]) gr->n;
  const uint triangles = gr->triangles;
  plistnode *v2t = bem->v2t;

  field    *aa = N->a;
  uint      rows = ntrans ? N->cols : N->rows;
  uint      cols = ntrans ? N->rows : N->cols;
  longindex ld = N->ld;

  vreal     vt[3][3], vs[3][3], nx[3], ny[3];
  real     *quad_r, *quad_i;
  ptri_list tl, tl1;
  pvert_list vl;
  const uint *tri_t, *tri_s;
  plistnode v;
  real     *xq, *yq, *wq, *mass;
  uint      tp[3], sp[3], tri_tp[3], tri_sp[3];
  vreal     ct[3], cs[3], tx, sx, ty, sy, w, x[3], y[3], sum_r, sum_i, eval_r,
    eval_i, vcount, cmp;
  real      base, factor, factor2;
  uint      i, j, t, s, q, nq, vnq, remainder, cj;
  uint      ii, tt, ss, vv;
#ifdef USE_TRIQUADPOINTS
  real     *tri_tx, *tri_ty, *tri_tz, *tri_sx, *tri_sy, *tri_sz;
  uint      c;
  uint      q2;
  uint      nq2 = bem->sq->n_single;
  uint      vnq2 = ROUNDUP(nq2, VREAL);
#endif

  clear_amatrix(N);

  quad_r = allocreal(ROUNDUP(bem->sq->nmax, VREAL));
  quad_i = allocreal(ROUNDUP(bem->sq->nmax, VREAL));
  vreal     c_one = vset1(1.0);

  tl = NULL;

  cj = 0;
  for (j = 0; j < rows; ++j) {
    ii = (ridx == NULL ? j : ridx[j]);
    for (v = v2t[ii], vv = v->data; v->next != NULL;
	 v = v->next, vv = v->data) {

      tl1 = tl;
      while (tl1 && tl1->t != vv) {
	tl1 = tl1->next;
      }

      if (tl1 == NULL) {
	tl1 = tl = new_tri_list(tl);
	tl->t = vv;
	cj++;
      }

      tl1->vl = new_vert_list(tl1->vl);
      tl1->vl->v = j;
    }
  }

  for (t = 0, tl1 = tl; t < cj; t++, tl1 = tl1->next) {
    tt = tl1->t;
    assert(tt < triangles);
    factor = gr_g[tt] * bem->kernel_const;
    tri_s = gr_t[tt];
    for (i = 0; i < 3; ++i) {
      nx[i] = vload1(gr_n[tt] + i);
    }
#ifdef USE_TRIQUADPOINTS
    tri_tx = bem->sq->tri_x + tt * vnq2;
    tri_ty = bem->sq->tri_y + tt * vnq2;
    tri_tz = bem->sq->tri_z + tt * vnq2;
#endif
    for (s = 0; s < cols; s++) {
      ss = (cidx == NULL ? s : cidx[s]);
      assert(ss < triangles);
      factor2 = factor * gr_g[ss];
      tri_t = gr_t[ss];
      for (i = 0; i < 3; ++i) {
	ny[i] = vload1(gr_n[ss] + i);
      }
#ifdef USE_TRIQUADPOINTS
      tri_sx = bem->sq->tri_x + ss * vnq2;
      tri_sy = bem->sq->tri_y + ss * vnq2;
      tri_sz = bem->sq->tri_z + ss * vnq2;
#endif

#ifdef USE_TRIQUADPOINTS
      c =
	select_quadrature_singquad2d(bem->sq, tri_t, tri_s, tp, sp, &xq, &yq,
				     &wq, &nq, &base);
      if (nq % VREAL || nq2 % VREAL) {
	for (q = 0; q < VREAL; ++q) {
	  vcount[q] = q;
	}
      }

#else
      (void) select_quadrature_singquad2d(bem->sq, tri_t, tri_s, tp, sp, &xq,
					  &yq, &wq, &nq, &base);
      if (nq % VREAL) {
	for (q = 0; q < VREAL; ++q) {
	  vcount[q] = q;
	}
      }

#endif
      vnq = ROUNDUP(nq, VREAL);

#ifdef USE_TRIQUADPOINTS
      if (c == 0) {
	remainder = vnq2 - VREAL;

	for (i = 0; i < 3; ++i) {
	  tri_tp[i] = tri_t[i];
	  tri_sp[i] = tri_s[i];
	}

	for (q = 0; q < nq2; q++) {
	  x[0] = vload1(tri_tx + q);
	  x[1] = vload1(tri_ty + q);
	  x[2] = vload1(tri_tz + q);
	  for (q2 = 0; q2 < vnq2; q2 += VREAL) {
	    y[0] = vload(tri_sx + q2);
	    y[1] = vload(tri_sy + q2);
	    y[2] = vload(tri_sz + q2);

	    kernel(x, y, nx, ny, (void *) bem, &eval_r, &eval_i);

	    if (nq2 % VREAL && q2 >= remainder) {
	      cmp = vcmplt(vadd(vcount, vset1((real) q2)), vset1((real) nq2));
	      eval_r = vand(cmp, eval_r);
	      eval_i = vand(cmp, eval_i);
	    }

	    vstoreu(quad_r + q2 + q * nq2, eval_r);
	    vstoreu(quad_i + q2 + q * nq2, eval_i);
	  }
	}
      }
      else {
#endif
	remainder = vnq - VREAL;

	for (i = 0; i < 3; ++i) {
	  tri_tp[i] = tri_t[tp[i]];
	  tri_sp[i] = tri_s[sp[i]];
	}

	for (i = 0; i < 3; ++i) {	// x-, y-, z- component
	  for (j = 0; j < 3; ++j) {	// vertex A, B, C
	    vt[i][j] = vload1(gr_x[tri_tp[j]] + i);
	    vs[i][j] = vload1(gr_x[tri_sp[j]] + i);
	  }
	}

	for (q = 0; q < vnq; q += VREAL) {
	  tx = vload(xq + q);
	  sx = vload(xq + q + vnq);
	  ty = vload(yq + q);
	  sy = vload(yq + q + vnq);

	  ct[0] = vsub(c_one, tx);
	  ct[1] = vsub(tx, sx);
	  ct[2] = sx;
	  cs[0] = vsub(c_one, ty);
	  cs[1] = vsub(ty, sy);
	  cs[2] = sy;

	  for (i = 0; i < 3; ++i) {
	    x[i] = vdot3(vt[i], ct);
	    y[i] = vdot3(vs[i], cs);
	  }

	  kernel(x, y, nx, ny, (void *) bem, &eval_r, &eval_i);

	  if (nq % VREAL && q >= remainder) {
	    cmp = vcmplt(vadd(vcount, vset1((real) q)), vset1((real) nq));
	    eval_r = vand(cmp, eval_r);
	    eval_i = vand(cmp, eval_i);
	  }

	  vstore(quad_r + q, eval_r);
	  vstore(quad_i + q, eval_i);
	}
#ifdef USE_TRIQUADPOINTS
      }
#endif
      vl = tl1->vl;
      while (vl) {
	i = vl->v;
	if (i < rows) {
	  ii = ridx == NULL ? i : ridx[i];
	  for (j = 0; j < 3; ++j) {
	    if (ii == tri_tp[j]) {
	      sum_r = vsetzero();
	      sum_i = vsetzero();

	      for (q = 0; q < vnq; q += VREAL) {
		eval_r = vload(quad_r + q);
		eval_i = vload(quad_i + q);
		w = vload(wq + q);
		sum_r = vfmadd(w, eval_r, sum_r);
		sum_i = vfmadd(w, eval_i, sum_i);
	      }

	      if (ntrans) {
		aa[s + i * ld] += ((vreduce(sum_r) + base) * factor2)
		  - (vreduce(sum_i) * factor2) * I;
	      }
	      else {
		aa[i + s * ld] += ((vreduce(sum_r) + base) * factor2)
		  + (vreduce(sum_i) * factor2) * I;
	      }
	    }
	    wq += vnq;
	  }
	  wq -= 3 * vnq;
	}
	vl = vl->next;
      }

      if (bem->alpha != 0.0 && tt == ss) {

	for (i = 0; i < 3; ++i) {
	  tri_tp[i] = tri_t[i];
	}

	mass = bem->mass;
	factor2 = bem->alpha * gr_g[tt];

	vl = tl1->vl;
	while (vl) {
	  i = vl->v;
	  if (i < rows) {
	    ii = ridx == NULL ? i : ridx[i];
	    for (j = 0; j < 3; ++j) {
	      if (ii == tri_tp[j]) {
		if (ntrans) {
		  aa[s + i * ld] += factor2 * *mass;
		}
		else {
		  aa[i + s * ld] += factor2 * *mass;
		}
	      }
	      mass++;
	    }
	    mass = bem->mass;
	  }
	  vl = vl->next;
	}
      }
    }
  }
}
#endif

void
assemble_lc_near_bem3d(const uint * ridx, const uint * cidx, pcbem3d bem,
		       bool ntrans, pamatrix N, kernel_func3d kernel)
{
  const pcsurface3d gr = bem->gr;
  const     real(*gr_x)[3] = (const real(*)[3]) gr->x;
  const     uint(*gr_t)[3] = (const uint(*)[3]) gr->t;
  const preal gr_g = (const preal) gr->g;
  const     real(*gr_n)[3] = (const real(*)[3]) gr->n;
  const uint triangles = gr->triangles;
  plistnode *v2t = bem->v2t;

  field    *aa = N->a;
  uint      rows = ntrans ? N->cols : N->rows;
  uint      cols = ntrans ? N->rows : N->cols;
  longindex ld = N->ld;

  field    *quad;
  ptri_list tl, tl1;
  pvert_list vl;
  const real *A_t, *B_t, *C_t, *A_s, *B_s, *C_s, *ns, *nt;
  const uint *tri_t, *tri_s;
  plistnode v;
  real     *xq, *yq, *wq, *mass;
  uint      tp[3], sp[3], tri_tp[3], tri_sp[3];
  real      Ax, Bx, Cx, Ay, By, Cy, tx, sx, ty, sy, x[3], y[3], factor,
    factor2, base;
  field     res;
  uint      i, j, t, s, q, nq, vnq, cj;
  uint      ii, tt, ss, vv;

  clear_amatrix(N);

  quad = allocfield(bem->sq->nmax);

  tl = NULL;

  cj = 0;
  for (j = 0; j < rows; ++j) {
    ii = (ridx == NULL ? j : ridx[j]);
    for (v = v2t[ii], vv = v->data; v->next != NULL;
	 v = v->next, vv = v->data) {

      tl1 = tl;
      while (tl1 && tl1->t != vv) {
	tl1 = tl1->next;
      }

      if (tl1 == NULL) {
	tl1 = tl = new_tri_list(tl);
	tl->t = vv;
	cj++;
      }

      tl1->vl = new_vert_list(tl1->vl);
      tl1->vl->v = j;
    }
  }

  for (t = 0, tl1 = tl; t < cj; t++, tl1 = tl1->next) {
    tt = tl1->t;
    assert(tt < triangles);
    factor = gr_g[tt] * bem->kernel_const;
    tri_t = gr_t[tt];
    nt = gr_n[tt];
    for (s = 0; s < cols; s++) {
      ss = (cidx == NULL ? s : cidx[s]);
      assert(ss < triangles);
      tri_s = gr_t[ss];
      ns = gr_n[ss];

      factor2 = factor * gr_g[ss];

      select_quadrature_singquad2d(bem->sq, tri_t, tri_s, tp, sp, &xq, &yq,
				   &wq, &nq, &base);
      vnq = ROUNDUP(nq, VREAL);

      for (i = 0; i < 3; ++i) {
	tri_tp[i] = tri_t[tp[i]];
	tri_sp[i] = tri_s[sp[i]];
      }

      A_t = gr_x[tri_tp[0]];
      B_t = gr_x[tri_tp[1]];
      C_t = gr_x[tri_tp[2]];
      A_s = gr_x[tri_sp[0]];
      B_s = gr_x[tri_sp[1]];
      C_s = gr_x[tri_sp[2]];

      for (q = 0; q < nq; ++q) {
	tx = xq[q];
	sx = xq[q + vnq];
	ty = yq[q];
	sy = yq[q + vnq];
	Ax = 1.0 - tx;
	Bx = tx - sx;
	Cx = sx;
	Ay = 1.0 - ty;
	By = ty - sy;
	Cy = sy;

	x[0] = A_t[0] * Ax + B_t[0] * Bx + C_t[0] * Cx;
	x[1] = A_t[1] * Ax + B_t[1] * Bx + C_t[1] * Cx;
	x[2] = A_t[2] * Ax + B_t[2] * Bx + C_t[2] * Cx;

	y[0] = A_s[0] * Ay + B_s[0] * By + C_s[0] * Cy;
	y[1] = A_s[1] * Ay + B_s[1] * By + C_s[1] * Cy;
	y[2] = A_s[2] * Ay + B_s[2] * By + C_s[2] * Cy;

	quad[q] = kernel(x, y, nt, ns, (void *) bem);
      }

      vl = tl1->vl;
      while (vl) {
	i = vl->v;
	if (i < rows) {
	  ii = ridx == NULL ? i : ridx[i];
	  for (j = 0; j < 3; ++j) {
	    if (ii == tri_tp[j]) {
	      res = base;

	      for (q = 0; q < nq; ++q) {
		res += wq[q] * quad[q];
	      }

	      if (ntrans) {
		aa[s + i * ld] += res * factor2;
	      }
	      else {
		aa[i + s * ld] += res * factor2;
	      }
	    }
	    wq += vnq;
	  }
	  wq -= 3 * vnq;
	}
	vl = vl->next;
      }

      if (bem->alpha != 0.0 && tt == ss) {

	for (j = 0; j < 3; ++j) {
	  tri_tp[j] = tri_t[j];
	}

	mass = bem->mass;
	factor2 = bem->alpha * gr_g[tt];

	vl = tl1->vl;
	while (vl) {
	  i = vl->v;
	  if (j < cols) {
	    ii = ridx == NULL ? i : ridx[i];
	    for (j = 0; j < 3; ++j) {
	      if (ii == tri_tp[j]) {
		if (ntrans) {
		  aa[s + i * ld] += factor2 * *mass;
		}
		else {
		  aa[i + s * ld] += factor2 * *mass;
		}
	      }
	      mass++;
	    }
	    mass = bem->mass;
	  }
	  vl = vl->next;
	}
      }
    }
  }

  del_tri_list(tl);
  freemem(quad);
}

#ifdef USE_SIMD
void
assemble_lc_simd_far_bem3d(const uint * ridx, const uint * cidx,
			   pcbem3d bem, bool ntrans, pamatrix N,
			   kernel_simd_func3d kernel)
{
  //TODO implement nice 'far' version
  assemble_lc_simd_near_bem3d(ridx, cidx, bem, ntrans, N, kernel);
}
#endif

void
assemble_lc_far_bem3d(const uint * ridx, const uint * cidx, pcbem3d bem,
		      bool ntrans, pamatrix N, kernel_func3d kernel)
{
  const pcsurface3d gr = bem->gr;
  const     real(*gr_x)[3] = (const real(*)[3]) gr->x;
  const     uint(*gr_t)[3] = (const uint(*)[3]) gr->t;
  const preal gr_g = (const preal) gr->g;
  const     real(*gr_n)[3] = (const real(*)[3]) gr->n;
  const uint triangles = gr->triangles;
  plistnode *v2t = bem->v2t;
  field    *quad;
  field    *aa = N->a;
  uint      rows = ntrans ? N->cols : N->rows;
  uint      cols = ntrans ? N->rows : N->cols;
  longindex ld = N->ld;

  ptri_list tl, tl1;
  pvert_list vl;
  const real *A_t, *B_t, *C_t, *A_s, *B_s, *C_s, *ns, *nt;
  const uint *tri_t, *tri_s;
  plistnode v;
  real     *xq, *yq, *wq;
  real      Ax, Bx, Cx, Ay, By, Cy, tx, sx, ty, sy, x[3], y[3], factor,
    factor2, base;
  field     res;
  uint      i, j, t, s, q, nq, vnq, cj;
  uint      ii, tt, ss, vv;

  clear_amatrix(N);

  quad = allocfield(bem->sq->nmax);

  tl = NULL;

  xq = bem->sq->x_dist;
  yq = bem->sq->y_dist;
  wq = bem->sq->w_dist;
  nq = bem->sq->n_dist;
  vnq = ROUNDUP(nq, VREAL);
  base = bem->sq->base_dist;

  cj = 0;
  for (j = 0; j < rows; ++j) {
    ii = (ridx == NULL ? j : ridx[j]);
    for (v = v2t[ii], vv = v->data; v->next != NULL;
	 v = v->next, vv = v->data) {

      tl1 = tl;
      while (tl1 && tl1->t != vv) {
	tl1 = tl1->next;
      }

      if (tl1 == NULL) {
	tl1 = tl = new_tri_list(tl);
	tl->t = vv;
	cj++;
      }

      tl1->vl = new_vert_list(tl1->vl);
      tl1->vl->v = j;
    }
  }

  for (t = 0, tl1 = tl; t < cj; t++, tl1 = tl1->next) {
    tt = tl1->t;
    assert(tt < triangles);
    factor = gr_g[tt] * bem->kernel_const;
    tri_t = gr_t[tt];
    nt = gr_n[tt];

    A_s = gr_x[tri_t[0]];
    B_s = gr_x[tri_t[1]];
    C_s = gr_x[tri_t[2]];

    for (s = 0; s < cols; s++) {
      ss = (cidx == NULL ? s : cidx[s]);
      assert(ss < triangles);
      tri_s = gr_t[ss];
      factor2 = factor * gr_g[ss];
      ns = gr_n[ss];

      A_t = gr_x[tri_s[0]];
      B_t = gr_x[tri_s[1]];
      C_t = gr_x[tri_s[2]];

      for (q = 0; q < nq; ++q) {
	tx = xq[q];
	sx = xq[q + vnq];
	ty = yq[q];
	sy = yq[q + vnq];
	Ax = 1.0 - tx;
	Bx = tx - sx;
	Cx = sx;
	Ay = 1.0 - ty;
	By = ty - sy;
	Cy = sy;

	x[0] = A_t[0] * Ax + B_t[0] * Bx + C_t[0] * Cx;
	x[1] = A_t[1] * Ax + B_t[1] * Bx + C_t[1] * Cx;
	x[2] = A_t[2] * Ax + B_t[2] * Bx + C_t[2] * Cx;

	y[0] = A_s[0] * Ay + B_s[0] * By + C_s[0] * Cy;
	y[1] = A_s[1] * Ay + B_s[1] * By + C_s[1] * Cy;
	y[2] = A_s[2] * Ay + B_s[2] * By + C_s[2] * Cy;

	quad[q] = kernel(x, y, nt, ns, (void *) bem);
      }

      vl = tl1->vl;
      while (vl) {
	i = vl->v;
	if (i < rows) {
	  ii = ridx == NULL ? i : ridx[i];
	  for (j = 0; j < 3; ++j) {
	    if (ii == tri_t[j]) {
	      res = base;

	      for (q = 0; q < nq; ++q) {
		res += wq[q] * quad[q];
	      }

	      if (ntrans) {
		aa[s + i * ld] += res * factor2;
	      }
	      else {
		aa[i + s * ld] += res * factor2;
	      }
	    }
	    wq += vnq;
	  }
	  wq -= 3 * vnq;
	}
	vl = vl->next;
      }
    }
  }

  del_tri_list(tl);
  freemem(quad);
}

#ifdef USE_SIMD
void
assemble_ll_simd_near_bem3d(const uint * ridx, const uint * cidx,
			    pcbem3d bem, bool ntrans, pamatrix N,
			    kernel_simd_func3d kernel)
{
  const pcsurface3d gr = bem->gr;
  const     real(*gr_x)[3] = (const real(*)[3]) gr->x;
  const     uint(*gr_t)[3] = (const uint(*)[3]) gr->t;
  const     real(*gr_n)[3] = (const real(*)[3]) gr->n;
  const preal gr_g = (const preal) gr->g;
  const uint triangles = gr->triangles;
  plistnode *v2t = bem->v2t;
  field    *aa = N->a;
  uint      rows = ntrans ? N->cols : N->rows;
  uint      cols = ntrans ? N->rows : N->cols;
  longindex ld = N->ld;

  ptri_list tl_r, tl1_r, tl_c, tl1_c;
  pvert_list vl_r, vl_c;
  vreal     vt[3][3], vs[3][3], nx[3], ny[3];
  real     *quad_r, *quad_i;
  const uint *tri_t, *tri_s;
  plistnode v;
  real     *xq, *yq, *wq, *ww;
  uint      tp[3], sp[3], tri_tp[3], tri_sp[3];
  vreal     ct[3], cs[3], tx, sx, ty, sy, w, x[3], y[3], sum_r, sum_i, eval_r,
    eval_i, vcount, cmp;
  real      base, factor, factor2;
  real     *mass;
  uint      i, j, t, s, q, nq, vnq, remainder, rj, cj;
  uint      ii, jj, tt, ss, vv, k, l;
#ifdef USE_TRIQUADPOINTS
  real     *tri_tx, *tri_ty, *tri_tz, *tri_sx, *tri_sy, *tri_sz;
  uint      c;
  uint      q2;
  uint      nq2 = bem->sq->n_single;
  uint      vnq2 = ROUNDUP(nq2, VREAL);
#endif

  quad_r = allocreal(ROUNDUP(bem->sq->nmax, VREAL));
  quad_i = allocreal(ROUNDUP(bem->sq->nmax, VREAL));
  vreal     c_one = vset1(1.0);

  clear_amatrix(N);

  tl_r = NULL;
  tl_c = NULL;

  rj = 0;
  for (i = 0; i < rows; ++i) {
    ii = (ridx == NULL ? i : ridx[i]);
    for (v = v2t[ii], vv = v->data; v->next != NULL;
	 v = v->next, vv = v->data) {

      tl1_r = tl_r;
      while (tl1_r && tl1_r->t != vv) {
	tl1_r = tl1_r->next;
      }

      if (tl1_r == NULL) {
	tl1_r = tl_r = new_tri_list(tl_r);
	tl_r->t = vv;
	rj++;
      }

      tl1_r->vl = new_vert_list(tl1_r->vl);
      tl1_r->vl->v = i;
    }
  }

  cj = 0;
  for (i = 0; i < cols; ++i) {
    ii = (cidx == NULL ? i : cidx[i]);
    for (v = v2t[ii], vv = v->data; v->next != NULL;
	 v = v->next, vv = v->data) {

      tl1_c = tl_c;
      while (tl1_c && tl1_c->t != vv) {
	tl1_c = tl1_c->next;
      }

      if (tl1_c == NULL) {
	tl1_c = tl_c = new_tri_list(tl_c);
	tl_c->t = vv;
	cj++;
      }

      tl1_c->vl = new_vert_list(tl1_c->vl);
      tl1_c->vl->v = i;
    }
  }

  for (s = 0, tl1_c = tl_c; s < cj; s++, tl1_c = tl1_c->next) {
    ss = tl1_c->t;
    assert(ss < triangles);
    factor = gr_g[ss] * bem->kernel_const;
    tri_s = gr_t[ss];
    for (i = 0; i < 3; ++i) {
      ny[i] = vload1(gr_n[ss] + i);
    }
#ifdef USE_TRIQUADPOINTS
    tri_sx = bem->sq->tri_x + ss * vnq2;
    tri_sy = bem->sq->tri_y + ss * vnq2;
    tri_sz = bem->sq->tri_z + ss * vnq2;
#endif
    for (t = 0, tl1_r = tl_r; t < rj; t++, tl1_r = tl1_r->next) {
      tt = tl1_r->t;
      assert(tt < triangles);
      factor2 = factor * gr_g[tt];
      tri_t = gr_t[tt];
      for (i = 0; i < 3; ++i) {
	nx[i] = vload1(gr_n[tt] + i);
      }
#ifdef USE_TRIQUADPOINTS
      tri_tx = bem->sq->tri_x + tt * vnq2;
      tri_ty = bem->sq->tri_y + tt * vnq2;
      tri_tz = bem->sq->tri_z + tt * vnq2;
#endif

#ifdef USE_TRIQUADPOINTS
      c =
	select_quadrature_singquad2d(bem->sq, tri_t, tri_s, tp, sp, &xq, &yq,
				     &wq, &nq, &base);
      if (nq % VREAL || nq2 % VREAL) {
	for (q = 0; q < VREAL; ++q) {
	  vcount[q] = q;
	}
      }

#else
      (void) select_quadrature_singquad2d(bem->sq, tri_t, tri_s, tp, sp, &xq,
					  &yq, &wq, &nq, &base);
      if (nq % VREAL) {
	for (q = 0; q < VREAL; ++q) {
	  vcount[q] = q;
	}
      }

#endif
      vnq = ROUNDUP(nq, VREAL);

#ifdef USE_TRIQUADPOINTS
      if (c == 0) {
	remainder = vnq2 - VREAL;

	for (i = 0; i < 3; ++i) {
	  tri_tp[i] = tri_t[i];
	  tri_sp[i] = tri_s[i];
	}

	for (q = 0; q < nq2; q++) {
	  x[0] = vload1(tri_tx + q);
	  x[1] = vload1(tri_ty + q);
	  x[2] = vload1(tri_tz + q);
	  for (q2 = 0; q2 < vnq2; q2 += VREAL) {
	    y[0] = vload(tri_sx + q2);
	    y[1] = vload(tri_sy + q2);
	    y[2] = vload(tri_sz + q2);

	    kernel(x, y, nx, ny, (void *) bem, &eval_r, &eval_i);

	    if (nq2 % VREAL && q2 >= remainder) {
	      cmp = vcmplt(vadd(vcount, vset1((real) q2)), vset1((real) nq2));
	      eval_r = vand(cmp, eval_r);
	      eval_i = vand(cmp, eval_i);
	    }

	    vstoreu(quad_r + q2 + q * nq2, eval_r);
	    vstoreu(quad_i + q2 + q * nq2, eval_i);
	  }
	}
      }
      else {
#endif
	remainder = vnq - VREAL;

	for (i = 0; i < 3; ++i) {
	  tri_tp[i] = tri_t[tp[i]];
	  tri_sp[i] = tri_s[sp[i]];
	}

	for (i = 0; i < 3; ++i) {	// x-, y-, z- component
	  for (j = 0; j < 3; ++j) {	// vertex A, B, C
	    vt[i][j] = vload1(gr_x[tri_tp[j]] + i);
	    vs[i][j] = vload1(gr_x[tri_sp[j]] + i);
	  }
	}

	for (q = 0; q < vnq; q += VREAL) {
	  tx = vload(xq + q);
	  sx = vload(xq + q + vnq);
	  ty = vload(yq + q);
	  sy = vload(yq + q + vnq);

	  ct[0] = vsub(c_one, tx);
	  ct[1] = vsub(tx, sx);
	  ct[2] = sx;
	  cs[0] = vsub(c_one, ty);
	  cs[1] = vsub(ty, sy);
	  cs[2] = sy;

	  for (i = 0; i < 3; ++i) {
	    x[i] = vdot3(vt[i], ct);
	    y[i] = vdot3(vs[i], cs);
	  }

	  kernel(x, y, nx, ny, (void *) bem, &eval_r, &eval_i);

	  if (nq % VREAL && q >= remainder) {
	    cmp = vcmplt(vadd(vcount, vset1((real) q)), vset1((real) nq));
	    eval_r = vand(cmp, eval_r);
	    eval_i = vand(cmp, eval_i);
	  }

	  vstore(quad_r + q, eval_r);
	  vstore(quad_i + q, eval_i);
	}
#ifdef USE_TRIQUADPOINTS
      }
#endif

      vl_c = tl1_c->vl;
      while (vl_c) {
	j = vl_c->v;
	assert(j < cols);
	jj = ((cidx == NULL) ? j : cidx[j]);
	for (k = 0; k < 3; ++k) {
	  if (jj == tri_sp[k]) {
	    vl_r = tl1_r->vl;
	    while (vl_r) {
	      i = vl_r->v;
	      assert(i < rows);
	      ii = ((ridx == NULL) ? i : ridx[i]);
	      for (l = 0; l < 3; ++l) {
		if (ii == tri_tp[l]) {
		  sum_r = vsetzero();
		  sum_i = vsetzero();
		  ww = wq + (l + k * 3) * vnq;

		  for (q = 0; q < vnq; q += VREAL) {
		    eval_r = vload(quad_r + q);
		    eval_i = vload(quad_i + q);
		    w = vload(ww + q);
		    sum_r = vfmadd(w, eval_r, sum_r);
		    sum_i = vfmadd(w, eval_i, sum_i);
		  }

		  if (ntrans) {
		    aa[j + i * ld] += ((vreduce(sum_r) + base) * factor2)
		      - (vreduce(sum_i) * factor2) * I;
		  }
		  else {
		    aa[i + j * ld] += ((vreduce(sum_r) + base) * factor2)
		      + (vreduce(sum_i) * factor2) * I;
		  }

		}
	      }
	      vl_r = vl_r->next;
	    }
	  }
	}
	vl_c = vl_c->next;
      }

      if (bem->alpha != 0.0 && tt == ss) {
	for (i = 0; i < 3; ++i) {
	  tri_tp[i] = tri_t[i];
	  tri_sp[i] = tri_s[i];
	}

	mass = bem->mass;
	factor2 = bem->alpha * gr_g[tt];

	vl_c = tl1_c->vl;
	while (vl_c) {
	  j = vl_c->v;
	  assert(j < cols);
	  jj = ((cidx == NULL) ? j : cidx[j]);
	  for (k = 0; k < 3; ++k) {
	    if (jj == tri_sp[k]) {
	      vl_r = tl1_r->vl;
	      while (vl_r) {
		i = vl_r->v;
		assert(i < rows);
		ii = ((ridx == NULL) ? i : ridx[i]);
		for (l = 0; l < 3; ++l) {
		  if (ii == tri_tp[l]) {
		    if (ntrans) {
		      aa[j + i * ld] += CONJ(mass[l + k * 3] * factor2);
		    }
		    else {
		      aa[i + j * ld] += mass[l + k * 3] * factor2;
		    }
		  }
		}
		vl_r = vl_r->next;
	      }
	    }
	  }
	  vl_c = vl_c->next;
	}
      }
    }
  }

  del_tri_list(tl_r);
  del_tri_list(tl_c);

  freemem(quad_r);
  freemem(quad_i);
}
#endif

void
assemble_ll_near_bem3d(const uint * ridx, const uint * cidx, pcbem3d bem,
		       bool ntrans, pamatrix N, kernel_func3d kernel)
{
  const pcsurface3d gr = bem->gr;
  const     real(*gr_x)[3] = (const real(*)[3]) gr->x;
  const     uint(*gr_t)[3] = (const uint(*)[3]) gr->t;
  const     real(*gr_n)[3] = (const real(*)[3]) gr->n;
  const preal gr_g = (const preal) gr->g;
  const uint triangles = gr->triangles;
  plistnode *v2t = bem->v2t;
  field    *quad;
  field    *aa = N->a;
  uint      rows = ntrans ? N->cols : N->rows;
  uint      cols = ntrans ? N->rows : N->cols;
  longindex ld = N->ld;

  ptri_list tl_r, tl1_r, tl_c, tl1_c;
  pvert_list vl_r, vl_c;
  const real *A_t, *B_t, *C_t, *A_s, *B_s, *C_s, *nt, *ns;
  const uint *tri_t, *tri_s;
  plistnode v;
  real     *xq, *yq, *wq, *ww;
  uint      tp[3], sp[3], tri_tp[3], tri_sp[3];
  real      Ax, Bx, Cx, Ay, By, Cy, tx, sx, ty, sy, x[3], y[3], factor,
    factor2, base;
  field     res;
  real     *mass;
  uint      i, j, t, s, k, l, rj, cj, tt, ss, q, nq, vnq, ii, jj, vv;

  quad = allocfield(bem->sq->nmax);

  clear_amatrix(N);

  tl_r = NULL;
  tl_c = NULL;

  rj = 0;
  for (i = 0; i < rows; ++i) {
    ii = (ridx == NULL ? i : ridx[i]);
    for (v = v2t[ii], vv = v->data; v->next != NULL;
	 v = v->next, vv = v->data) {

      tl1_r = tl_r;
      while (tl1_r && tl1_r->t != vv) {
	tl1_r = tl1_r->next;
      }

      if (tl1_r == NULL) {
	tl1_r = tl_r = new_tri_list(tl_r);
	tl_r->t = vv;
	rj++;
      }

      tl1_r->vl = new_vert_list(tl1_r->vl);
      tl1_r->vl->v = i;
    }
  }

  cj = 0;
  for (i = 0; i < cols; ++i) {
    ii = (cidx == NULL ? i : cidx[i]);
    for (v = v2t[ii], vv = v->data; v->next != NULL;
	 v = v->next, vv = v->data) {

      tl1_c = tl_c;
      while (tl1_c && tl1_c->t != vv) {
	tl1_c = tl1_c->next;
      }

      if (tl1_c == NULL) {
	tl1_c = tl_c = new_tri_list(tl_c);
	tl_c->t = vv;
	cj++;
      }

      tl1_c->vl = new_vert_list(tl1_c->vl);
      tl1_c->vl->v = i;
    }
  }

  for (s = 0, tl1_c = tl_c; s < cj; s++, tl1_c = tl1_c->next) {
    ss = tl1_c->t;
    assert(ss < triangles);
    factor = gr_g[ss] * bem->kernel_const;
    tri_s = gr_t[ss];
    ns = gr_n[ss];
    for (t = 0, tl1_r = tl_r; t < rj; t++, tl1_r = tl1_r->next) {
      tt = tl1_r->t;
      assert(tt < triangles);
      factor2 = factor * gr_g[tt];
      tri_t = gr_t[tt];
      nt = gr_n[tt];

      select_quadrature_singquad2d(bem->sq, tri_t, tri_s, tp, sp, &xq, &yq,
				   &wq, &nq, &base);
      vnq = ROUNDUP(nq, VREAL);

      for (i = 0; i < 3; ++i) {
	tri_tp[i] = tri_t[tp[i]];
	tri_sp[i] = tri_s[sp[i]];
      }

      A_t = gr_x[tri_tp[0]];
      B_t = gr_x[tri_tp[1]];
      C_t = gr_x[tri_tp[2]];
      A_s = gr_x[tri_sp[0]];
      B_s = gr_x[tri_sp[1]];
      C_s = gr_x[tri_sp[2]];

      for (q = 0; q < nq; ++q) {
	tx = xq[q];
	sx = xq[q + vnq];
	ty = yq[q];
	sy = yq[q + vnq];
	Ax = 1.0 - tx;
	Bx = tx - sx;
	Cx = sx;
	Ay = 1.0 - ty;
	By = ty - sy;
	Cy = sy;

	x[0] = A_t[0] * Ax + B_t[0] * Bx + C_t[0] * Cx;
	x[1] = A_t[1] * Ax + B_t[1] * Bx + C_t[1] * Cx;
	x[2] = A_t[2] * Ax + B_t[2] * Bx + C_t[2] * Cx;
	y[0] = A_s[0] * Ay + B_s[0] * By + C_s[0] * Cy;
	y[1] = A_s[1] * Ay + B_s[1] * By + C_s[1] * Cy;
	y[2] = A_s[2] * Ay + B_s[2] * By + C_s[2] * Cy;

	quad[q] = kernel(x, y, nt, ns, (void *) bem);
      }

      vl_c = tl1_c->vl;
      while (vl_c) {
	j = vl_c->v;
	assert(j < cols);
	jj = ((cidx == NULL) ? j : cidx[j]);
	for (k = 0; k < 3; ++k) {
	  if (jj == tri_sp[k]) {
	    vl_r = tl1_r->vl;
	    while (vl_r) {
	      i = vl_r->v;
	      assert(i < rows);
	      ii = ((ridx == NULL) ? i : ridx[i]);
	      for (l = 0; l < 3; ++l) {
		if (ii == tri_tp[l]) {
		  res = base;

		  ww = wq + (l + k * 3) * vnq;
		  for (q = 0; q < nq; ++q) {
		    res += ww[q] * quad[q];
		  }

		  if (ntrans) {
		    aa[j + i * ld] += CONJ(res * factor2);
		  }
		  else {
		    aa[i + j * ld] += res * factor2;
		  }
		}
	      }
	      vl_r = vl_r->next;
	    }
	  }
	}
	vl_c = vl_c->next;
      }

      if (bem->alpha != 0.0 && tt == ss) {
	for (i = 0; i < 3; ++i) {
	  tri_tp[i] = tri_t[i];
	  tri_sp[i] = tri_s[i];
	}

	mass = bem->mass;
	factor2 = bem->alpha * gr_g[tt];

	vl_c = tl1_c->vl;
	while (vl_c) {
	  j = vl_c->v;
	  assert(j < cols);
	  jj = ((cidx == NULL) ? j : cidx[j]);
	  for (k = 0; k < 3; ++k) {
	    if (jj == tri_sp[k]) {
	      vl_r = tl1_r->vl;
	      while (vl_r) {
		i = vl_r->v;
		assert(i < rows);
		ii = ((ridx == NULL) ? i : ridx[i]);
		for (l = 0; l < 3; ++l) {
		  if (ii == tri_tp[l]) {
		    if (ntrans) {
		      aa[j + i * ld] += CONJ(mass[l + k * 3] * factor2);
		    }
		    else {
		      aa[i + j * ld] += mass[l + k * 3] * factor2;
		    }
		  }
		}
		vl_r = vl_r->next;
	      }
	    }
	  }
	  vl_c = vl_c->next;
	}
      }
    }
  }

  del_tri_list(tl_r);
  del_tri_list(tl_c);

  freemem(quad);
}

#ifdef USE_SIMD
void
assemble_ll_simd_far_bem3d(const uint * ridx, const uint * cidx,
			   pcbem3d bem, bool ntrans, pamatrix N,
			   kernel_simd_func3d kernel)
{
  //TODO implement nice 'far' version
  assemble_ll_simd_near_bem3d(ridx, cidx, bem, ntrans, N, kernel);
}
#endif

void
assemble_ll_far_bem3d(const uint * ridx, const uint * cidx, pcbem3d bem,
		      bool ntrans, pamatrix N, kernel_func3d kernel)
{
  const pcsurface3d gr = bem->gr;
  const     real(*gr_x)[3] = (const real(*)[3]) gr->x;
  const     uint(*gr_t)[3] = (const uint(*)[3]) gr->t;
  const     real(*gr_n)[3] = (const real(*)[3]) gr->n;
  const preal gr_g = (const preal) gr->g;
  const uint triangles = gr->triangles;
  plistnode *v2t = bem->v2t;
  field    *quad;
  field    *aa = N->a;
  uint      rows = ntrans ? N->cols : N->rows;
  uint      cols = ntrans ? N->rows : N->cols;
  longindex ld = N->ld;

  ptri_list tl_r, tl1_r, tl_c, tl1_c;
  pvert_list vl_r, vl_c;
  const real *A_t, *B_t, *C_t, *A_s, *B_s, *C_s, *nt, *ns;
  const uint *tri_t, *tri_s;
  plistnode v;
  real     *xq, *yq, *wq, *ww;
  real      Ax, Bx, Cx, Ay, By, Cy, tx, sx, ty, sy, x[3], y[3], factor,
    factor2, base;
  field     res;
  uint      i, j, t, s, k, l, rj, cj, tt, ss, q, nq, vnq, ii, jj, vv;

  xq = bem->sq->x_dist;
  yq = bem->sq->y_dist;
  nq = bem->sq->n_dist;
  vnq = ROUNDUP(nq, VREAL);
  wq = bem->sq->w_dist;

  quad = allocfield(bem->sq->nmax);

  clear_amatrix(N);

  tl_r = NULL;
  tl_c = NULL;

  rj = 0;
  for (i = 0; i < rows; ++i) {
    ii = (ridx == NULL ? i : ridx[i]);
    for (v = v2t[ii], vv = v->data; v->next != NULL;
	 v = v->next, vv = v->data) {

      tl1_r = tl_r;
      while (tl1_r && tl1_r->t != vv) {
	tl1_r = tl1_r->next;
      }

      if (tl1_r == NULL) {
	tl1_r = tl_r = new_tri_list(tl_r);
	tl_r->t = vv;
	rj++;
      }

      tl1_r->vl = new_vert_list(tl1_r->vl);
      tl1_r->vl->v = i;
    }
  }

  cj = 0;
  for (i = 0; i < cols; ++i) {
    ii = (cidx == NULL ? i : cidx[i]);
    for (v = v2t[ii], vv = v->data; v->next != NULL;
	 v = v->next, vv = v->data) {

      tl1_c = tl_c;
      while (tl1_c && tl1_c->t != vv) {
	tl1_c = tl1_c->next;
      }

      if (tl1_c == NULL) {
	tl1_c = tl_c = new_tri_list(tl_c);
	tl_c->t = vv;
	cj++;
      }

      tl1_c->vl = new_vert_list(tl1_c->vl);
      tl1_c->vl->v = i;
    }
  }

  for (s = 0, tl1_c = tl_c; s < cj; s++, tl1_c = tl1_c->next) {
    ss = tl1_c->t;
    assert(ss < triangles);
    factor = gr_g[ss] * bem->kernel_const;
    tri_s = gr_t[ss];
    ns = gr_n[ss];

    A_s = gr_x[tri_s[0]];
    B_s = gr_x[tri_s[1]];
    C_s = gr_x[tri_s[2]];

    for (t = 0, tl1_r = tl_r; t < rj; t++, tl1_r = tl1_r->next) {
      tt = tl1_r->t;
      assert(tt < triangles);
      factor2 = factor * gr_g[tt];
      tri_t = gr_t[tt];
      nt = gr_n[tt];

      A_t = gr_x[tri_t[0]];
      B_t = gr_x[tri_t[1]];
      C_t = gr_x[tri_t[2]];

      base = bem->sq->base_dist;

      for (q = 0; q < nq; ++q) {
	tx = xq[q];
	sx = xq[q + vnq];
	ty = yq[q];
	sy = yq[q + vnq];
	Ax = 1.0 - tx;
	Bx = tx - sx;
	Cx = sx;
	Ay = 1.0 - ty;
	By = ty - sy;
	Cy = sy;

	x[0] = A_t[0] * Ax + B_t[0] * Bx + C_t[0] * Cx;
	x[1] = A_t[1] * Ax + B_t[1] * Bx + C_t[1] * Cx;
	x[2] = A_t[2] * Ax + B_t[2] * Bx + C_t[2] * Cx;
	y[0] = A_s[0] * Ay + B_s[0] * By + C_s[0] * Cy;
	y[1] = A_s[1] * Ay + B_s[1] * By + C_s[1] * Cy;
	y[2] = A_s[2] * Ay + B_s[2] * By + C_s[2] * Cy;

	quad[q] = kernel(x, y, nt, ns, (void *) bem);
      }

      vl_c = tl1_c->vl;
      while (vl_c) {
	j = vl_c->v;
	assert(j < cols);
	jj = ((cidx == NULL) ? j : cidx[j]);
	for (k = 0; k < 3; ++k) {
	  if (jj == tri_s[k]) {
	    vl_r = tl1_r->vl;
	    while (vl_r) {
	      i = vl_r->v;
	      assert(i < rows);
	      ii = ((ridx == NULL) ? i : ridx[i]);
	      for (l = 0; l < 3; ++l) {
		if (ii == tri_t[l]) {
		  res = base;

		  ww = wq + (l + k * 3) * vnq;
		  for (q = 0; q < nq; ++q) {
		    res += ww[q] * quad[q];
		  }

		  if (ntrans) {
		    aa[j + i * ld] += CONJ(res * factor2);
		  }
		  else {
		    aa[i + j * ld] += res * factor2;
		  }
		}
	      }
	      vl_r = vl_r->next;
	    }
	  }
	}
	vl_c = vl_c->next;
      }
    }
  }

  del_tri_list(tl_r);
  del_tri_list(tl_c);

  freemem(quad);
}

void
fill_bem3d(pcbem3d bem, const real(*X)[3], const real(*Y)[3],
	   const real(*NX)[3], const real(*NY)[3], pamatrix V,
	   kernel_func3d kernel)
{
  uint      rows = V->rows;
  uint      cols = V->cols;
  longindex ld = V->ld;
  field     c = bem->kernel_const;

  uint      i, j;

  if (NX == NULL && NY == NULL) {
    for (j = 0; j < cols; ++j) {
      for (i = 0; i < rows; ++i) {
	V->a[i + j * ld] = c * kernel(X[i], Y[j], NULL, NULL, (void *) bem);
      }
    }
    return;
  }

  if (NX != NULL && NY == NULL) {
    for (j = 0; j < cols; ++j) {
      for (i = 0; i < rows; ++i) {
	V->a[i + j * ld] = c * kernel(X[i], Y[j], NX[i], NULL, (void *) bem);
      }
    }
    return;
  }

  if (NX == NULL && NY != NULL) {
    for (j = 0; j < cols; ++j) {
      for (i = 0; i < rows; ++i) {
	V->a[i + j * ld] = c * kernel(X[i], Y[j], NULL, NY[j], (void *) bem);
      }
    }
    return;
  }

  if (NX != NULL && NY != NULL) {
    for (j = 0; j < cols; ++j) {
      for (i = 0; i < rows; ++i) {
	V->a[i + j * ld] = c * kernel(X[i], Y[j], NX[i], NY[j], (void *) bem);
      }
    }
    return;
  }

}

void
fill_wave_bem3d(pcbem3d bem, const real(*X)[3], const real(*Y)[3],
		const real(*NX)[3], const real(*NY)[3], pamatrix V,
		pcreal dir, kernel_wave_func3d kernel)
{
  uint      rows = V->rows;
  uint      cols = V->cols;
  longindex ld = V->ld;
  field     c = bem->kernel_const;

  uint      i, j;

  if (NX == NULL && NY == NULL) {
    for (j = 0; j < cols; ++j) {
      for (i = 0; i < rows; ++i) {
	V->a[i + j * ld] = c
	  * kernel(X[i], Y[j], NULL, NULL, dir, (void *) bem);
      }
    }
    return;
  }

  if (NX != NULL && NY == NULL) {
    for (j = 0; j < cols; ++j) {
      for (i = 0; i < rows; ++i) {
	V->a[i + j * ld] = c
	  * kernel(X[i], Y[j], NX[i], NULL, dir, (void *) bem);
      }
    }
    return;
  }

  if (NX == NULL && NY != NULL) {
    for (j = 0; j < cols; ++j) {
      for (i = 0; i < rows; ++i) {
	V->a[i + j * ld] = c
	  * kernel(X[i], Y[j], NULL, NY[j], dir, (void *) bem);
      }
    }
    return;
  }

  if (NX != NULL && NY != NULL) {
    for (j = 0; j < cols; ++j) {
      for (i = 0; i < rows; ++i) {
	V->a[i + j * ld] = c
	  * kernel(X[i], Y[j], NX[i], NY[j], dir, (void *) bem);
      }
    }
    return;
  }
}

void
fill_row_c_bem3d(const uint * idx, const real(*Z)[3], pcbem3d bem,
		 pamatrix V, kernel_func3d kernel)
{
  pcsurface3d gr = bem->gr;
  const     real(*gr_x)[3] = (const real(*)[3]) gr->x;
  const     uint(*gr_t)[3] = (const uint(*)[3]) gr->t;
  const     real(*gr_n)[3] = (const real(*)[3]) gr->n;
  const preal gr_g = (const preal) gr->g;
  uint      rows = V->rows;
  uint      cols = V->cols;
  longindex ld = V->ld;

  uint      nq = bem->sq->n_single;
  uint      vnq = ROUNDUP(nq, VREAL);
  real     *xx = bem->sq->x_single;
  real     *yy = bem->sq->y_single;
  real     *ww = bem->sq->w_single + 3 * vnq;

  const real *A, *B, *C, *n;
  uint      t, tt, i, q;
  real      gt_fac, x[3], tx, sx, Ax, Bx, Cx;
  field     sum;

  for (t = 0; t < rows; ++t) {
    tt = (idx == NULL ? t : idx[t]);
    gt_fac = gr_g[tt] * bem->kernel_const;
    A = gr_x[gr_t[tt][0]];
    B = gr_x[gr_t[tt][1]];
    C = gr_x[gr_t[tt][2]];
    n = gr_n[tt];

    for (i = 0; i < cols; ++i) {

      sum = 0.0;

      for (q = 0; q < nq; ++q) {
	tx = xx[q];
	sx = yy[q];
	Ax = 1.0 - tx;
	Bx = tx - sx;
	Cx = sx;

	x[0] = A[0] * Ax + B[0] * Bx + C[0] * Cx;
	x[1] = A[1] * Ax + B[1] * Bx + C[1] * Cx;
	x[2] = A[2] * Ax + B[2] * Bx + C[2] * Cx;

	sum += ww[q] * kernel(x, Z[i], n, NULL, (void *) bem);
      }

      V->a[t + i * ld] = sum * gt_fac;
    }
  }
}

#ifdef USE_SIMD
#ifdef USE_TRIQUADPOINTS
void
fill_row_simd_c_bem3d(const uint * idx, const real(*Z)[3], pcbem3d bem,
		      pamatrix V, kernel_simd_func3d kernel)
{
  pcsurface3d gr = bem->gr;
  const     real(*gr_n)[3] = (const real(*)[3]) gr->n;
  const preal gr_g = (const preal) gr->g;
  uint      rows = V->rows;
  uint      cols = V->cols;
  longindex ld = V->ld;

  uint      nq = bem->sq->n_single;
  uint      vnq = ROUNDUP(nq, VREAL);
  real     *wq = bem->sq->w_single + 3 * vnq;
  real      base = bem->sq->base_single;

  const real *n;
  uint      t, tt, i, j, q, remainder;
  real      gt_fac;
  field     sum;
  real     *tri_tx, *tri_ty, *tri_tz;

  vreal     x[3], z[3], nx[3];
  vreal     w, eval_r, eval_i, sum_r, sum_i, cmp, vcount;

  remainder = vnq - VREAL;

  if (nq % VREAL) {
    for (q = 0; q < VREAL; ++q) {
      vcount[q] = q;
    }
  }

  for (t = 0; t < rows; ++t) {
    tt = (idx == NULL ? t : idx[t]);
    gt_fac = gr_g[tt] * bem->kernel_const;
    n = gr_n[tt];
    tri_tx = bem->sq->tri_x + tt * vnq;
    tri_ty = bem->sq->tri_y + tt * vnq;
    tri_tz = bem->sq->tri_z + tt * vnq;

    for (j = 0; j < 3; ++j) {
      nx[j] = vset1(n[j]);
    }

    for (i = 0; i < cols; ++i) {
      sum = 0.0;

      sum_r = vsetzero();
      sum_i = vsetzero();

      for (j = 0; j < 3; ++j) {
	z[j] = vset1(Z[i][j]);
      }

      for (q = 0; q < vnq; q += VREAL) {
	x[0] = vload(tri_tx + q);
	x[1] = vload(tri_ty + q);
	x[2] = vload(tri_tz + q);
	w = vload(wq + q);

	kernel(x, z, nx, NULL, (void *) bem, &eval_r, &eval_i);

	if (nq % VREAL && q >= remainder) {
	  cmp = vcmplt(vadd(vcount, vset1((real) q)), vset1((real) nq));
	  eval_r = vand(cmp, eval_r);
	  eval_i = vand(cmp, eval_i);
	}

	sum_r = vfmadd(w, eval_r, sum_r);
	sum_i = vfmadd(w, eval_i, sum_i);
      }

      sum =
	((vreduce(sum_r) + base) * gt_fac) + (vreduce(sum_i) * gt_fac) * I;

      V->a[t + i * ld] = sum;
    }
  }
}
#else
void
fill_row_simd_c_bem3d(const uint * idx, const real(*Z)[3], pcbem3d bem,
		      pamatrix V, kernel_simd_func3d kernel)
{
  pcsurface3d gr = bem->gr;
  const     real(*gr_x)[3] = (const real(*)[3]) gr->x;
  const     uint(*gr_t)[3] = (const uint(*)[3]) gr->t;
  const     real(*gr_n)[3] = (const real(*)[3]) gr->n;
  const preal gr_g = (const preal) gr->g;
  uint      rows = V->rows;
  uint      cols = V->cols;
  longindex ld = V->ld;

  uint      nq = bem->sq->n_single;
  uint      vnq = ROUNDUP(nq, VREAL);
  real     *xq = bem->sq->x_single;
  real     *yq = bem->sq->y_single;
  real     *wq = bem->sq->w_single + 3 * vnq;
  real      base = bem->sq->base_single;

  const real *n;
  const uint *tri_t;
  uint      t, tt, i, j, q, remainder;
  real      gt_fac;
  field     sum;

  vreal     vt[3][3], ct[3], x[3], z[3], nx[3];
  vreal     tx, sx, w, c_one, eval_r, eval_i, sum_r, sum_i, cmp, vcount;

  c_one = vset1(1.0);

  remainder = vnq - VREAL;

  if (nq % VREAL) {
    for (q = 0; q < VREAL; ++q) {
      vcount[q] = q;
    }
  }

  for (t = 0; t < rows; ++t) {
    tt = (idx == NULL ? t : idx[t]);
    gt_fac = gr_g[tt] * bem->kernel_const;
    tri_t = gr_t[tt];
    n = gr_n[tt];

    for (j = 0; j < 3; ++j) {
      nx[j] = vset1(n[j]);
    }

    for (i = 0; i < 3; ++i) {	// x-, y-, z- component
      for (j = 0; j < 3; ++j) {	// vertex A, B, C
	vt[i][j] = vload1(gr_x[tri_t[j]] + i);
      }
    }

    for (i = 0; i < cols; ++i) {

      sum = 0.0;

      sum_r = vsetzero();
      sum_i = vsetzero();

      for (j = 0; j < 3; ++j) {
	z[j] = vset1(Z[i][j]);
      }

      for (q = 0; q < vnq; q += VREAL) {
	tx = vload(xq + q);
	sx = vload(yq + q);
	w = vload(wq + q);

	ct[0] = vsub(c_one, tx);
	ct[1] = vsub(tx, sx);
	ct[2] = sx;

	for (j = 0; j < 3; ++j) {
	  x[j] = vdot3(vt[j], ct);
	}

	kernel(x, z, nx, NULL, (void *) bem, &eval_r, &eval_i);

	if (nq % VREAL && q >= remainder) {
	  cmp = vcmplt(vadd(vcount, vset1((real) q)), vset1((real) nq));
	  eval_r = vand(cmp, eval_r);
	  eval_i = vand(cmp, eval_i);
	}

	sum_r = vfmadd(w, eval_r, sum_r);
	sum_i = vfmadd(w, eval_i, sum_i);
      }

      sum =
	((vreduce(sum_r) + base) * gt_fac) + (vreduce(sum_i) * gt_fac) * I;

      V->a[t + i * ld] = sum;
    }
  }
}
#endif
#endif

void
fill_col_c_bem3d(const uint * idx, const real(*Z)[3], pcbem3d bem,
		 pamatrix V, kernel_func3d kernel)
{
  pcsurface3d gr = bem->gr;
  const     real(*gr_x)[3] = (const real(*)[3]) gr->x;
  const     uint(*gr_t)[3] = (const uint(*)[3]) gr->t;
  const     real(*gr_n)[3] = (const real(*)[3]) gr->n;
  const preal gr_g = (const preal) gr->g;
  uint      rows = V->rows;
  uint      cols = V->cols;
  longindex ld = V->ld;

  uint      nq = bem->sq->n_single;
  uint      vnq = ROUNDUP(nq, VREAL);
  real     *xx = bem->sq->x_single;
  real     *yy = bem->sq->y_single;
  real     *ww = bem->sq->w_single + 3 * vnq;

  const real *A, *B, *C, *n;
  uint      s, ss, i, q;
  real      gs_fac, x[3], tx, sx, Ax, Bx, Cx;
  field     sum;

  for (s = 0; s < rows; ++s) {
    ss = (idx == NULL ? s : idx[s]);
    gs_fac = gr_g[ss] * bem->kernel_const;
    A = gr_x[gr_t[ss][0]];
    B = gr_x[gr_t[ss][1]];
    C = gr_x[gr_t[ss][2]];
    n = gr_n[ss];

    for (i = 0; i < cols; ++i) {

      sum = 0.0;

      for (q = 0; q < nq; ++q) {
	tx = xx[q];
	sx = yy[q];
	Ax = 1.0 - tx;
	Bx = tx - sx;
	Cx = sx;

	x[0] = A[0] * Ax + B[0] * Bx + C[0] * Cx;
	x[1] = A[1] * Ax + B[1] * Bx + C[1] * Cx;
	x[2] = A[2] * Ax + B[2] * Bx + C[2] * Cx;

	sum += ww[q] * kernel(Z[i], x, NULL, n, (void *) bem);
      }

      V->a[s + i * ld] = sum * gs_fac;
    }
  }
}

#ifdef USE_SIMD
#ifdef USE_TRIQUADPOINTS
void
fill_col_simd_c_bem3d(const uint * idx, const real(*Z)[3], pcbem3d bem,
		      pamatrix V, kernel_simd_func3d kernel)
{
  pcsurface3d gr = bem->gr;
  const     real(*gr_n)[3] = (const real(*)[3]) gr->n;
  const preal gr_g = (const preal) gr->g;
  uint      rows = V->rows;
  uint      cols = V->cols;
  longindex ld = V->ld;

  uint      nq = bem->sq->n_single;
  uint      vnq = ROUNDUP(nq, VREAL);
  real     *wq = bem->sq->w_single + 3 * vnq;
  real      base = bem->sq->base_single;

  const real *n;
  uint      s, ss, i, j, q, remainder;
  real      gs_fac;
  field     sum;
  real     *tri_sx, *tri_sy, *tri_sz;

  vreal     y[3], z[3], ny[3];
  vreal     w, eval_r, eval_i, sum_r, sum_i, cmp, vcount;

  remainder = vnq - VREAL;

  if (nq % VREAL) {
    for (q = 0; q < VREAL; ++q) {
      vcount[q] = q;
    }
  }

  for (s = 0; s < rows; ++s) {
    ss = (idx == NULL ? s : idx[s]);
    gs_fac = gr_g[ss] * bem->kernel_const;
    n = gr_n[ss];
    tri_sx = bem->sq->tri_x + ss * vnq;
    tri_sy = bem->sq->tri_y + ss * vnq;
    tri_sz = bem->sq->tri_z + ss * vnq;

    for (j = 0; j < 3; ++j) {
      ny[j] = vset1(n[j]);
    }

    for (i = 0; i < cols; ++i) {
      sum = 0.0;

      sum_r = vsetzero();
      sum_i = vsetzero();

      for (j = 0; j < 3; ++j) {
	z[j] = vset1(Z[i][j]);
      }

      for (q = 0; q < vnq; q += VREAL) {
	y[0] = vload(tri_sx + q);
	y[1] = vload(tri_sy + q);
	y[2] = vload(tri_sz + q);
	w = vload(wq + q);

	kernel(z, y, NULL, ny, (void *) bem, &eval_r, &eval_i);

	if (nq % VREAL && q >= remainder) {
	  cmp = vcmplt(vadd(vcount, vset1((real) q)), vset1((real) nq));
	  eval_r = vand(cmp, eval_r);
	  eval_i = vand(cmp, eval_i);
	}

	sum_r = vfmadd(w, eval_r, sum_r);
	sum_i = vfmadd(w, eval_i, sum_i);
      }

      sum =
	((vreduce(sum_r) + base) * gs_fac) + (vreduce(sum_i) * gs_fac) * I;

      V->a[s + i * ld] = sum;
    }
  }
}
#else
void
fill_col_simd_c_bem3d(const uint * idx, const real(*Z)[3], pcbem3d bem,
		      pamatrix V, kernel_simd_func3d kernel)
{
  pcsurface3d gr = bem->gr;
  const     real(*gr_x)[3] = (const real(*)[3]) gr->x;
  const     uint(*gr_t)[3] = (const uint(*)[3]) gr->t;
  const     real(*gr_n)[3] = (const real(*)[3]) gr->n;
  const preal gr_g = (const preal) gr->g;
  uint      rows = V->rows;
  uint      cols = V->cols;
  longindex ld = V->ld;

  uint      nq = bem->sq->n_single;
  uint      vnq = ROUNDUP(nq, VREAL);
  real     *xq = bem->sq->x_single;
  real     *yq = bem->sq->y_single;
  real     *wq = bem->sq->w_single + 3 * vnq;
  real      base = bem->sq->base_single;

  const real *n;
  const uint *tri_s;
  uint      s, ss, i, j, q, remainder;
  real      gs_fac;
  field     sum;

  vreal     vs[3][3], cs[3], y[3], z[3], ny[3];
  vreal     ty, sy, w, c_one, eval_r, eval_i, sum_r, sum_i, cmp, vcount;

  c_one = vset1(1.0);

  remainder = vnq - VREAL;

  if (nq % VREAL) {
    for (q = 0; q < VREAL; ++q) {
      vcount[q] = q;
    }
  }

  for (s = 0; s < rows; ++s) {
    ss = (idx == NULL ? s : idx[s]);
    gs_fac = gr_g[ss] * bem->kernel_const;
    tri_s = gr_t[ss];
    n = gr_n[ss];

    for (j = 0; j < 3; ++j) {
      ny[j] = vset1(n[j]);
    }

    for (i = 0; i < 3; ++i) {	// x-, y-, z- component
      for (j = 0; j < 3; ++j) {	// vertex A, B, C
	vs[i][j] = vload1(gr_x[tri_s[j]] + i);
      }
    }

    for (i = 0; i < cols; ++i) {

      sum = 0.0;

      sum_r = vsetzero();
      sum_i = vsetzero();

      for (j = 0; j < 3; ++j) {
	z[j] = vset1(Z[i][j]);
      }

      for (q = 0; q < vnq; q += VREAL) {
	ty = vload(xq + q);
	sy = vload(yq + q);
	w = vload(wq + q);

	cs[0] = vsub(c_one, ty);
	cs[1] = vsub(ty, sy);
	cs[2] = sy;

	for (j = 0; j < 3; ++j) {
	  y[j] = vdot3(vs[j], cs);
	}

	kernel(z, y, NULL, ny, (void *) bem, &eval_r, &eval_i);

	if (nq % VREAL && q >= remainder) {
	  cmp = vcmplt(vadd(vcount, vset1((real) q)), vset1((real) nq));
	  eval_r = vand(cmp, eval_r);
	  eval_i = vand(cmp, eval_i);
	}

	sum_r = vfmadd(w, eval_r, sum_r);
	sum_i = vfmadd(w, eval_i, sum_i);
      }

      sum =
	((vreduce(sum_r) + base) * gs_fac) + (vreduce(sum_i) * gs_fac) * I;

      V->a[s + i * ld] = sum;
    }
  }
}
#endif
#endif

void
fill_row_l_bem3d(const uint * idx, const real(*Z)[3], pcbem3d bem,
		 pamatrix V, kernel_func3d kernel)
{
  pcsurface3d gr = bem->gr;
  const     real(*gr_x)[3] = (const real(*)[3]) gr->x;
  const     uint(*gr_t)[3] = (const uint(*)[3]) gr->t;
  const     real(*gr_n)[3] = (const real(*)[3]) gr->n;
  const preal gr_g = (const preal) gr->g;
  plistnode *v2t = bem->v2t;
  uint      rows = V->rows;
  uint      cols = V->cols;
  field    *aa = V->a;
  longindex ld = V->ld;
  uint      nq = bem->sq->n_single;
  uint      vnq = ROUNDUP(nq, VREAL);
  real     *xx = bem->sq->x_single;
  real     *yy = bem->sq->y_single;
  real     *ww = bem->sq->w_single;
  real      base = bem->sq->base_single;
  field    *quad;

  ptri_list tl, tl1;
  pvert_list vl;
  const real *A, *B, *C, *n;
  const uint *tri_t;
  plistnode v;
  uint      t, i, j, k, q, rj;
  real      gt_fac, x[3], tx, sx, Ax, Bx, Cx;
  longindex ii, tt, vv;
  field     sum;

  quad = allocfield(vnq);

  clear_amatrix(V);

  tl = NULL;

  rj = 0;
  for (i = 0; i < rows; ++i) {
    ii = (idx == NULL ? i : idx[i]);
    for (v = v2t[ii], vv = v->data; v->next != NULL;
	 v = v->next, vv = v->data) {

      tl1 = tl;
      while (tl1 && tl1->t != vv) {
	tl1 = tl1->next;
      }

      if (tl1 == NULL) {
	tl1 = tl = new_tri_list(tl);
	tl->t = vv;
	rj++;
      }

      tl1->vl = new_vert_list(tl1->vl);
      tl1->vl->v = i;
    }
  }

  for (t = 0, tl1 = tl; t < rj; t++, tl1 = tl1->next) {
    tt = tl1->t;
    tri_t = gr_t[tt];
    gt_fac = gr_g[tt] * bem->kernel_const;
    A = gr_x[gr_t[tt][0]];
    B = gr_x[gr_t[tt][1]];
    C = gr_x[gr_t[tt][2]];
    n = gr_n[tt];

    for (j = 0; j < cols; ++j) {

      for (q = 0; q < nq; ++q) {
	tx = xx[q];
	sx = yy[q];
	Ax = 1.0 - tx;
	Bx = tx - sx;
	Cx = sx;

	x[0] = A[0] * Ax + B[0] * Bx + C[0] * Cx;
	x[1] = A[1] * Ax + B[1] * Bx + C[1] * Cx;
	x[2] = A[2] * Ax + B[2] * Bx + C[2] * Cx;

	quad[q] = kernel(x, Z[j], n, NULL, (void *) bem);
      }

      ww = bem->sq->w_single;
      vl = tl1->vl;
      while (vl) {
	k = vl->v;
	if (k < rows) {
	  ii = idx == NULL ? k : idx[k];
	  for (i = 0; i < 3; ++i) {
	    if (ii == tri_t[i]) {
	      sum = base;

	      for (q = 0; q < nq; ++q) {
		sum += ww[q] * quad[q];
	      }

	      aa[k + j * ld] += sum * gt_fac;
	    }
	    ww += vnq;
	  }
	  ww = bem->sq->w_single;
	}
	vl = vl->next;
      }

    }
  }

  del_tri_list(tl);

  freemem(quad);
}

void
fill_col_l_bem3d(const uint * idx, const real(*Z)[3], pcbem3d bem,
		 pamatrix V, kernel_func3d kernel)
{
  pcsurface3d gr = bem->gr;
  const     real(*gr_x)[3] = (const real(*)[3]) gr->x;
  const     uint(*gr_t)[3] = (const uint(*)[3]) gr->t;
  const     real(*gr_n)[3] = (const real(*)[3]) gr->n;
  const preal gr_g = (const preal) gr->g;
  plistnode *v2t = bem->v2t;
  uint      rows = V->rows;
  uint      cols = V->cols;
  field    *aa = V->a;
  longindex ld = V->ld;
  uint      nq = bem->sq->n_single;
  uint      vnq = ROUNDUP(nq, VREAL);
  real     *xx = bem->sq->x_single;
  real     *yy = bem->sq->y_single;
  real     *ww = bem->sq->w_single;
  real      base = bem->sq->base_single;
  field    *quad;

  ptri_list tl, tl1;
  pvert_list vl;
  const real *A, *B, *C, *n;
  const uint *tri_s;
  plistnode v;
  uint      s, i, j, k, q, rj;
  real      gs_fac, x[3], tx, sx, Ax, Bx, Cx;
  longindex ii, ss, vv;
  field     sum;

  quad = allocfield(vnq);

  clear_amatrix(V);

  tl = NULL;

  rj = 0;
  for (i = 0; i < rows; ++i) {
    ii = (idx == NULL ? i : idx[i]);
    for (v = v2t[ii], vv = v->data; v->next != NULL;
	 v = v->next, vv = v->data) {

      tl1 = tl;
      while (tl1 && tl1->t != vv) {
	tl1 = tl1->next;
      }

      if (tl1 == NULL) {
	tl1 = tl = new_tri_list(tl);
	tl->t = vv;
	rj++;
      }

      tl1->vl = new_vert_list(tl1->vl);
      tl1->vl->v = i;
    }
  }

  for (s = 0, tl1 = tl; s < rj; s++, tl1 = tl1->next) {
    ss = tl1->t;
    tri_s = gr_t[ss];
    gs_fac = gr_g[ss] * bem->kernel_const;
    A = gr_x[gr_t[ss][0]];
    B = gr_x[gr_t[ss][1]];
    C = gr_x[gr_t[ss][2]];
    n = gr_n[ss];

    for (j = 0; j < cols; ++j) {

      for (q = 0; q < nq; ++q) {
	tx = xx[q];
	sx = yy[q];
	Ax = 1.0 - tx;
	Bx = tx - sx;
	Cx = sx;

	x[0] = A[0] * Ax + B[0] * Bx + C[0] * Cx;
	x[1] = A[1] * Ax + B[1] * Bx + C[1] * Cx;
	x[2] = A[2] * Ax + B[2] * Bx + C[2] * Cx;

	quad[q] = kernel(Z[j], x, NULL, n, (void *) bem);
      }

      ww = bem->sq->w_single;
      vl = tl1->vl;
      while (vl) {
	k = vl->v;
	if (k < rows) {
	  ii = idx == NULL ? k : idx[k];
	  for (i = 0; i < 3; ++i) {
	    if (ii == tri_s[i]) {
	      sum = base;

	      for (q = 0; q < nq; ++q) {
		sum += ww[q] * quad[q];
	      }

	      aa[k + j * ld] += sum * gs_fac;
	    }
	    ww += vnq;
	  }
	  ww = bem->sq->w_single;
	}
	vl = vl->next;
      }

    }
  }

  del_tri_list(tl);

  freemem(quad);
}

void
fill_dnz_row_c_bem3d(const uint * idx, const real(*Z)[3],
		     const real(*N)[3], pcbem3d bem, pamatrix V,
		     kernel_func3d kernel)
{
  pcsurface3d gr = bem->gr;
  const     real(*gr_x)[3] = (const real(*)[3]) gr->x;
  const     uint(*gr_t)[3] = (const uint(*)[3]) gr->t;
  const     real(*gr_n)[3] = (const real(*)[3]) gr->n;
  const preal gr_g = (const preal) gr->g;
  uint      rows = V->rows;
  uint      cols = V->cols;
  longindex ld = V->ld;

  uint      nq = bem->sq->n_single;
  uint      vnq = ROUNDUP(nq, VREAL);
  real     *xx = bem->sq->x_single;
  real     *yy = bem->sq->y_single;
  real     *ww = bem->sq->w_single + 3 * vnq;

  const real *A, *B, *C, *n;
  uint      t, tt, i, q;
  real      gt_fac, x[3], tx, sx, Ax, Bx, Cx;
  field     sum;

  for (t = 0; t < rows; ++t) {
    tt = (idx == NULL ? t : idx[t]);
    gt_fac = gr_g[tt] * bem->kernel_const;
    A = gr_x[gr_t[tt][0]];
    B = gr_x[gr_t[tt][1]];
    C = gr_x[gr_t[tt][2]];
    n = gr_n[tt];

    for (i = 0; i < cols; ++i) {

      sum = 0.0;

      for (q = 0; q < nq; ++q) {
	tx = xx[q];
	sx = yy[q];
	Ax = 1.0 - tx;
	Bx = tx - sx;
	Cx = sx;

	x[0] = A[0] * Ax + B[0] * Bx + C[0] * Cx;
	x[1] = A[1] * Ax + B[1] * Bx + C[1] * Cx;
	x[2] = A[2] * Ax + B[2] * Bx + C[2] * Cx;

	sum += ww[q] * kernel(x, Z[i], n, N[i], (void *) bem);
      }

      V->a[t + i * ld] = sum * gt_fac;
    }
  }
}

#ifdef USE_SIMD
#ifdef USE_TRIQUADPOINTS
void
fill_dnz_row_simd_c_bem3d(const uint * idx, const real(*Z)[3],
			  const real(*N)[3], pcbem3d bem, pamatrix V,
			  kernel_simd_func3d kernel)
{
  pcsurface3d gr = bem->gr;
  const     real(*gr_n)[3] = (const real(*)[3]) gr->n;
  const preal gr_g = (const preal) gr->g;
  uint      rows = V->rows;
  uint      cols = V->cols;
  longindex ld = V->ld;

  uint      nq = bem->sq->n_single;
  uint      vnq = ROUNDUP(nq, VREAL);
  real     *wq = bem->sq->w_single + 3 * vnq;
  real      base = bem->sq->base_single;

  const real *n;
  uint      t, tt, i, j, q, remainder;
  real      gt_fac;
  field     sum;
  real     *tri_tx, *tri_ty, *tri_tz;

  vreal     x[3], z[3], nx[3], nz[3];
  vreal     w, eval_r, eval_i, sum_r, sum_i, cmp, vcount;

  remainder = vnq - VREAL;

  if (nq % VREAL) {
    for (q = 0; q < VREAL; ++q) {
      vcount[q] = q;
    }
  }

  for (t = 0; t < rows; ++t) {
    tt = (idx == NULL ? t : idx[t]);
    gt_fac = gr_g[tt] * bem->kernel_const;
    n = gr_n[tt];
    tri_tx = bem->sq->tri_x + tt * vnq;
    tri_ty = bem->sq->tri_y + tt * vnq;
    tri_tz = bem->sq->tri_z + tt * vnq;

    for (j = 0; j < 3; ++j) {
      nx[j] = vset1(n[j]);
    }

    for (i = 0; i < cols; ++i) {
      sum = 0.0;

      sum_r = vsetzero();
      sum_i = vsetzero();

      for (j = 0; j < 3; ++j) {
	z[j] = vset1(Z[i][j]);
	nz[j] = vset1(N[i][j]);
      }

      for (q = 0; q < vnq; q += VREAL) {
	x[0] = vload(tri_tx + q);
	x[1] = vload(tri_ty + q);
	x[2] = vload(tri_tz + q);
	w = vload(wq + q);

	kernel(x, z, nx, nz, (void *) bem, &eval_r, &eval_i);

	if (nq % VREAL && q >= remainder) {
	  cmp = vcmplt(vadd(vcount, vset1((real) q)), vset1((real) nq));
	  eval_r = vand(cmp, eval_r);
	  eval_i = vand(cmp, eval_i);
	}

	sum_r = vfmadd(w, eval_r, sum_r);
	sum_i = vfmadd(w, eval_i, sum_i);
      }

      sum =
	((vreduce(sum_r) + base) * gt_fac) + (vreduce(sum_i) * gt_fac) * I;

      V->a[t + i * ld] = sum;
    }
  }
}
#else
void
fill_dnz_row_simd_c_bem3d(const uint * idx, const real(*Z)[3],
			  const real(*N)[3], pcbem3d bem, pamatrix V,
			  kernel_simd_func3d kernel)
{
  pcsurface3d gr = bem->gr;
  const     real(*gr_x)[3] = (const real(*)[3]) gr->x;
  const     uint(*gr_t)[3] = (const uint(*)[3]) gr->t;
  const     real(*gr_n)[3] = (const real(*)[3]) gr->n;
  const preal gr_g = (const preal) gr->g;
  uint      rows = V->rows;
  uint      cols = V->cols;
  longindex ld = V->ld;

  uint      nq = bem->sq->n_single;
  uint      vnq = ROUNDUP(nq, VREAL);
  real     *xq = bem->sq->x_single;
  real     *yq = bem->sq->y_single;
  real     *wq = bem->sq->w_single + 3 * vnq;
  real      base = bem->sq->base_single;

  const real *n;
  const uint *tri_t;
  uint      t, tt, i, j, q, remainder;
  real      gt_fac;
  field     sum;

  vreal     vt[3][3], ct[3], x[3], z[3], nx[3], nz[3];
  vreal     tx, sx, w, c_one, eval_r, eval_i, sum_r, sum_i, cmp, vcount;

  c_one = vset1(1.0);

  remainder = vnq - VREAL;

  if (nq % VREAL) {
    for (q = 0; q < VREAL; ++q) {
      vcount[q] = q;
    }
  }

  for (t = 0; t < rows; ++t) {
    tt = (idx == NULL ? t : idx[t]);
    gt_fac = gr_g[tt] * bem->kernel_const;
    tri_t = gr_t[tt];
    n = gr_n[tt];

    for (j = 0; j < 3; ++j) {
      nx[j] = vset1(n[j]);
    }

    for (i = 0; i < 3; ++i) {	// x-, y-, z- component
      for (j = 0; j < 3; ++j) {	// vertex A, B, C
	vt[i][j] = vload1(gr_x[tri_t[j]] + i);
      }
    }

    for (i = 0; i < cols; ++i) {

      sum = 0.0;

      sum_r = vsetzero();
      sum_i = vsetzero();

      for (j = 0; j < 3; ++j) {
	z[j] = vset1(Z[i][j]);
	nz[j] = vset1(N[i][j]);
      }

      for (q = 0; q < vnq; q += VREAL) {
	tx = vload(xq + q);
	sx = vload(yq + q);
	w = vload(wq + q);

	ct[0] = vsub(c_one, tx);
	ct[1] = vsub(tx, sx);
	ct[2] = sx;

	for (j = 0; j < 3; ++j) {
	  x[j] = vdot3(vt[j], ct);
	}

	kernel(x, z, nx, nz, (void *) bem, &eval_r, &eval_i);

	if (nq % VREAL && q >= remainder) {
	  cmp = vcmplt(vadd(vcount, vset1((real) q)), vset1((real) nq));
	  eval_r = vand(cmp, eval_r);
	  eval_i = vand(cmp, eval_i);
	}

	sum_r = vfmadd(w, eval_r, sum_r);
	sum_i = vfmadd(w, eval_i, sum_i);
      }

      sum =
	((vreduce(sum_r) + base) * gt_fac) + (vreduce(sum_i) * gt_fac) * I;

      V->a[t + i * ld] = sum;
    }
  }
}
#endif
#endif

void
fill_dnz_col_c_bem3d(const uint * idx, const real(*Z)[3],
		     const real(*N)[3], pcbem3d bem, pamatrix V,
		     kernel_func3d kernel)
{
  pcsurface3d gr = bem->gr;
  const     real(*gr_x)[3] = (const real(*)[3]) gr->x;
  const     uint(*gr_t)[3] = (const uint(*)[3]) gr->t;
  const     real(*gr_n)[3] = (const real(*)[3]) gr->n;
  const preal gr_g = (const preal) gr->g;
  uint      rows = V->rows;
  uint      cols = V->cols;
  longindex ld = V->ld;

  uint      nq = bem->sq->n_single;
  uint      vnq = ROUNDUP(nq, VREAL);
  real     *xx = bem->sq->x_single;
  real     *yy = bem->sq->y_single;
  real     *ww = bem->sq->w_single + 3 * vnq;

  const real *A, *B, *C, *n;
  uint      s, ss, i, q;
  real      gs_fac, x[3], tx, sx, Ax, Bx, Cx;
  field     sum;

  for (s = 0; s < rows; ++s) {
    ss = (idx == NULL ? s : idx[s]);
    gs_fac = gr_g[ss] * bem->kernel_const;
    A = gr_x[gr_t[ss][0]];
    B = gr_x[gr_t[ss][1]];
    C = gr_x[gr_t[ss][2]];
    n = gr_n[ss];

    for (i = 0; i < cols; ++i) {

      sum = 0.0;

      for (q = 0; q < nq; ++q) {
	tx = xx[q];
	sx = yy[q];
	Ax = 1.0 - tx;
	Bx = tx - sx;
	Cx = sx;

	x[0] = A[0] * Ax + B[0] * Bx + C[0] * Cx;
	x[1] = A[1] * Ax + B[1] * Bx + C[1] * Cx;
	x[2] = A[2] * Ax + B[2] * Bx + C[2] * Cx;

	sum += ww[q] * kernel(Z[i], x, N[i], n, (void *) bem);
      }

      V->a[s + i * ld] = sum * gs_fac;
    }
  }
}

#ifdef USE_SIMD
#ifdef USE_TRIQUADPOINTS
void
fill_dnz_col_simd_c_bem3d(const uint * idx, const real(*Z)[3],
			  const real(*N)[3], pcbem3d bem, pamatrix V,
			  kernel_simd_func3d kernel)
{
  pcsurface3d gr = bem->gr;
  const     real(*gr_n)[3] = (const real(*)[3]) gr->n;
  const preal gr_g = (const preal) gr->g;
  uint      rows = V->rows;
  uint      cols = V->cols;
  longindex ld = V->ld;

  uint      nq = bem->sq->n_single;
  uint      vnq = ROUNDUP(nq, VREAL);
  real     *wq = bem->sq->w_single + 3 * vnq;
  real      base = bem->sq->base_single;

  const real *n;
  uint      s, ss, i, j, q, remainder;
  real      gs_fac;
  field     sum;
  real     *tri_sx, *tri_sy, *tri_sz;

  vreal     y[3], z[3], ny[3], nz[3];
  vreal     w, eval_r, eval_i, sum_r, sum_i, cmp, vcount;

  remainder = vnq - VREAL;

  if (nq % VREAL) {
    for (q = 0; q < VREAL; ++q) {
      vcount[q] = q;
    }
  }

  for (s = 0; s < rows; ++s) {
    ss = (idx == NULL ? s : idx[s]);
    gs_fac = gr_g[ss] * bem->kernel_const;
    n = gr_n[ss];
    tri_sx = bem->sq->tri_x + ss * vnq;
    tri_sy = bem->sq->tri_y + ss * vnq;
    tri_sz = bem->sq->tri_z + ss * vnq;

    for (j = 0; j < 3; ++j) {
      ny[j] = vset1(n[j]);
    }

    for (i = 0; i < cols; ++i) {
      sum = 0.0;

      sum_r = vsetzero();
      sum_i = vsetzero();

      for (j = 0; j < 3; ++j) {
	z[j] = vset1(Z[i][j]);
	nz[j] = vset1(N[i][j]);
      }

      for (q = 0; q < vnq; q += VREAL) {
	y[0] = vload(tri_sx + q);
	y[1] = vload(tri_sy + q);
	y[2] = vload(tri_sz + q);
	w = vload(wq + q);

	kernel(z, y, nz, ny, (void *) bem, &eval_r, &eval_i);

	if (nq % VREAL && q >= remainder) {
	  cmp = vcmplt(vadd(vcount, vset1((real) q)), vset1((real) nq));
	  eval_r = vand(cmp, eval_r);
	  eval_i = vand(cmp, eval_i);
	}

	sum_r = vfmadd(w, eval_r, sum_r);
	sum_i = vfmadd(w, eval_i, sum_i);
      }

      sum =
	((vreduce(sum_r) + base) * gs_fac) + (vreduce(sum_i) * gs_fac) * I;

      V->a[s + i * ld] = sum;
    }
  }
}
#else
void
fill_dnz_col_simd_c_bem3d(const uint * idx, const real(*Z)[3],
			  const real(*N)[3], pcbem3d bem, pamatrix V,
			  kernel_simd_func3d kernel)
{
  pcsurface3d gr = bem->gr;
  const     real(*gr_x)[3] = (const real(*)[3]) gr->x;
  const     uint(*gr_t)[3] = (const uint(*)[3]) gr->t;
  const     real(*gr_n)[3] = (const real(*)[3]) gr->n;
  const preal gr_g = (const preal) gr->g;
  uint      rows = V->rows;
  uint      cols = V->cols;
  longindex ld = V->ld;

  uint      nq = bem->sq->n_single;
  uint      vnq = ROUNDUP(nq, VREAL);
  real     *xq = bem->sq->x_single;
  real     *yq = bem->sq->y_single;
  real     *wq = bem->sq->w_single + 3 * vnq;
  real      base = bem->sq->base_single;

  const real *n;
  const uint *tri_s;
  uint      s, ss, i, j, q, remainder;
  real      gs_fac;
  field     sum;

  vreal     vs[3][3], cs[3], y[3], z[3], ny[3], nz[3];
  vreal     ty, sy, w, c_one, eval_r, eval_i, sum_r, sum_i, cmp, vcount;

  c_one = vset1(1.0);

  remainder = vnq - VREAL;

  if (nq % VREAL) {
    for (q = 0; q < VREAL; ++q) {
      vcount[q] = q;
    }
  }

  for (s = 0; s < rows; ++s) {
    ss = (idx == NULL ? s : idx[s]);
    gs_fac = gr_g[ss] * bem->kernel_const;
    tri_s = gr_t[ss];
    n = gr_n[ss];

    for (j = 0; j < 3; ++j) {
      ny[j] = vset1(n[j]);
    }

    for (i = 0; i < 3; ++i) {	// x-, y-, z- component
      for (j = 0; j < 3; ++j) {	// vertex A, B, C
	vs[i][j] = vload1(gr_x[tri_s[j]] + i);
      }
    }

    for (i = 0; i < cols; ++i) {

      sum = 0.0;

      sum_r = vsetzero();
      sum_i = vsetzero();

      for (j = 0; j < 3; ++j) {
	z[j] = vset1(Z[i][j]);
	nz[j] = vset1(N[i][j]);
      }

      for (q = 0; q < vnq; q += VREAL) {
	ty = vload(xq + q);
	sy = vload(yq + q);
	w = vload(wq + q);

	cs[0] = vsub(c_one, ty);
	cs[1] = vsub(ty, sy);
	cs[2] = sy;

	for (j = 0; j < 3; ++j) {
	  y[j] = vdot3(vs[j], cs);
	}

	kernel(z, y, nz, ny, (void *) bem, &eval_r, &eval_i);

	if (nq % VREAL && q >= remainder) {
	  cmp = vcmplt(vadd(vcount, vset1((real) q)), vset1((real) nq));
	  eval_r = vand(cmp, eval_r);
	  eval_i = vand(cmp, eval_i);
	}

	sum_r = vfmadd(w, eval_r, sum_r);
	sum_i = vfmadd(w, eval_i, sum_i);
      }

      sum =
	((vreduce(sum_r) + base) * gs_fac) + (vreduce(sum_i) * gs_fac) * I;

      V->a[s + i * ld] = sum;
    }
  }
}
#endif
#endif

void
fill_dnz_row_l_bem3d(const uint * idx, const real(*Z)[3],
		     const real(*N)[3], pcbem3d bem, pamatrix V,
		     kernel_func3d kernel)
{
  pcsurface3d gr = bem->gr;
  const     real(*gr_x)[3] = (const real(*)[3]) gr->x;
  const     uint(*gr_t)[3] = (const uint(*)[3]) gr->t;
  const     real(*gr_n)[3] = (const real(*)[3]) gr->n;
  const preal gr_g = (const preal) gr->g;
  plistnode *v2t = bem->v2t;
  uint      rows = V->rows;
  uint      cols = V->cols;
  field    *aa = V->a;
  longindex ld = V->ld;
  uint      nq = bem->sq->n_single;
  uint      vnq = ROUNDUP(nq, VREAL);
  real     *xx = bem->sq->x_single;
  real     *yy = bem->sq->y_single;
  real     *ww = bem->sq->w_single;
  real      base = bem->sq->base_single;
  field    *quad;

  ptri_list tl, tl1;
  pvert_list vl;
  const real *A, *B, *C, *n;
  const uint *tri_t;
  plistnode v;
  uint      t, i, j, k, q, rj;
  real      gt_fac, x[3], tx, sx, Ax, Bx, Cx;
  longindex ii, tt, vv;
  field     sum;

  quad = allocfield(vnq);

  clear_amatrix(V);

  tl = NULL;

  rj = 0;
  for (i = 0; i < rows; ++i) {
    ii = (idx == NULL ? i : idx[i]);
    for (v = v2t[ii], vv = v->data; v->next != NULL;
	 v = v->next, vv = v->data) {

      tl1 = tl;
      while (tl1 && tl1->t != vv) {
	tl1 = tl1->next;
      }

      if (tl1 == NULL) {
	tl1 = tl = new_tri_list(tl);
	tl->t = vv;
	rj++;
      }

      tl1->vl = new_vert_list(tl1->vl);
      tl1->vl->v = i;
    }
  }

  for (t = 0, tl1 = tl; t < rj; t++, tl1 = tl1->next) {
    tt = tl1->t;
    tri_t = gr_t[tt];
    gt_fac = gr_g[tt] * bem->kernel_const;
    A = gr_x[gr_t[tt][0]];
    B = gr_x[gr_t[tt][1]];
    C = gr_x[gr_t[tt][2]];
    n = gr_n[tt];

    for (j = 0; j < cols; ++j) {

      for (q = 0; q < nq; ++q) {
	tx = xx[q];
	sx = yy[q];
	Ax = 1.0 - tx;
	Bx = tx - sx;
	Cx = sx;

	x[0] = A[0] * Ax + B[0] * Bx + C[0] * Cx;
	x[1] = A[1] * Ax + B[1] * Bx + C[1] * Cx;
	x[2] = A[2] * Ax + B[2] * Bx + C[2] * Cx;

	quad[q] = kernel(x, Z[j], n, N[j], (void *) bem);
      }

      ww = bem->sq->w_single;
      vl = tl1->vl;
      while (vl) {
	k = vl->v;
	if (k < rows) {
	  ii = idx == NULL ? k : idx[k];
	  for (i = 0; i < 3; ++i) {
	    if (ii == tri_t[i]) {
	      sum = base;

	      for (q = 0; q < nq; ++q) {
		sum += ww[q] * quad[q];
	      }

	      aa[k + j * ld] += sum * gt_fac;
	    }
	    ww += vnq;
	  }
	  ww = bem->sq->w_single;
	}
	vl = vl->next;
      }

    }
  }

  del_tri_list(tl);

  freemem(quad);
}

void
fill_dnz_col_l_bem3d(const uint * idx, const real(*Z)[3],
		     const real(*N)[3], pcbem3d bem, pamatrix V,
		     kernel_func3d kernel)
{
  pcsurface3d gr = bem->gr;
  const     real(*gr_x)[3] = (const real(*)[3]) gr->x;
  const     uint(*gr_t)[3] = (const uint(*)[3]) gr->t;
  const     real(*gr_n)[3] = (const real(*)[3]) gr->n;
  const preal gr_g = (const preal) gr->g;
  plistnode *v2t = bem->v2t;
  uint      rows = V->rows;
  uint      cols = V->cols;
  field    *aa = V->a;
  longindex ld = V->ld;
  uint      nq = bem->sq->n_single;
  uint      vnq = ROUNDUP(nq, VREAL);
  real     *xx = bem->sq->x_single;
  real     *yy = bem->sq->y_single;
  real     *ww = bem->sq->w_single;
  real      base = bem->sq->base_single;
  field    *quad;

  ptri_list tl, tl1;
  pvert_list vl;
  const real *A, *B, *C, *n;
  const uint *tri_s;
  plistnode v;
  uint      s, i, j, k, q, rj;
  real      gs_fac, x[3], tx, sx, Ax, Bx, Cx;
  longindex ii, ss, vv;
  field     sum;

  quad = allocfield(vnq);

  clear_amatrix(V);

  tl = NULL;

  rj = 0;
  for (i = 0; i < rows; ++i) {
    ii = (idx == NULL ? i : idx[i]);
    for (v = v2t[ii], vv = v->data; v->next != NULL;
	 v = v->next, vv = v->data) {

      tl1 = tl;
      while (tl1 && tl1->t != vv) {
	tl1 = tl1->next;
      }

      if (tl1 == NULL) {
	tl1 = tl = new_tri_list(tl);
	tl->t = vv;
	rj++;
      }

      tl1->vl = new_vert_list(tl1->vl);
      tl1->vl->v = i;
    }
  }

  for (s = 0, tl1 = tl; s < rj; s++, tl1 = tl1->next) {
    ss = tl1->t;
    tri_s = gr_t[ss];
    gs_fac = gr_g[ss] * bem->kernel_const;
    A = gr_x[gr_t[ss][0]];
    B = gr_x[gr_t[ss][1]];
    C = gr_x[gr_t[ss][2]];
    n = gr_n[ss];

    for (j = 0; j < cols; ++j) {

      for (q = 0; q < nq; ++q) {
	tx = xx[q];
	sx = yy[q];
	Ax = 1.0 - tx;
	Bx = tx - sx;
	Cx = sx;

	x[0] = A[0] * Ax + B[0] * Bx + C[0] * Cx;
	x[1] = A[1] * Ax + B[1] * Bx + C[1] * Cx;
	x[2] = A[2] * Ax + B[2] * Bx + C[2] * Cx;

	quad[q] = kernel(Z[j], x, N[j], n, (void *) bem);
      }

      ww = bem->sq->w_single;
      vl = tl1->vl;
      while (vl) {
	k = vl->v;
	if (k < rows) {
	  ii = idx == NULL ? k : idx[k];
	  for (i = 0; i < 3; ++i) {
	    if (ii == tri_s[i]) {
	      sum = base;

	      for (q = 0; q < nq; ++q) {
		sum += ww[q] * quad[q];
	      }

	      aa[k + j * ld] += sum * gs_fac;
	    }
	    ww += vnq;
	  }
	  ww = bem->sq->w_single;
	}
	vl = vl->next;
      }

    }
  }

  del_tri_list(tl);

  freemem(quad);
}

/****************************************************
 * Setup approximation schemes
 ****************************************************/

static void
setup_interpolation_bem3d(paprxbem3d aprx, uint m)
{
  real      e;
  uint      i;

  uninit_interpolation_bem3d(aprx);

  aprx->x_inter = allocreal(m);
  aprx->m_inter = m;
  aprx->k_inter = m * m * m;

  /* build tschebyscheff-points */
  e = 1.0 / (2.0 * m);

  for (i = 0; i < m; ++i) {
    aprx->x_inter[i] = cos(M_PI * (2.0 * i * e + e));
  }

}

static void
setup_green_bem3d(paprxbem3d aprx, uint m, uint l, real delta,
		  quadpoints3d quadpoints)
{
  uint      i, j;
  real      h, c;
  real     *s, *ht, *hw;

  uninit_green_bem3d(aprx);

  s = allocreal(l + 1);
  ht = allocreal(m);
  hw = allocreal(m);

  aprx->m_green = m;
  aprx->l_green = l;
  aprx->delta_green = delta;
  aprx->ml_green = m * l;
  if (quadpoints == build_bem3d_cube_quadpoints) {
    aprx->k_green = 12 * m * l * m * l;
  }

  aprx->quadpoints = quadpoints;

  /* quadrature points and weights */
  if (m == 1) {
    ht[0] = 0.0;
    hw[0] = 2.0;
  }
  else {
    assemble_gauss(m, ht, hw);
  }

  /* partitioning the intervall */
  c = 2.0 / (real) l;
  for (i = 0; i <= l; i++) {
    s[i] = -1.0 + (c * i);
  }

  aprx->t_green = allocreal(aprx->ml_green);
  aprx->w_green = allocreal(aprx->ml_green);

  /* adjust quadrature points and weights */
  for (j = 1; j <= l; j++) {
    h = 0.5 * (s[j] - s[j - 1]);
    for (i = 0; i < m; i++) {
      aprx->t_green[i + (j - 1) * m] = h * ht[i] + 0.5 * (s[j] + s[j - 1]);
      aprx->w_green[i + (j - 1) * m] = h * hw[i];
    }
  }

  freemem(s);
  freemem(ht);
  freemem(hw);
}

static void
setup_aca_bem3d(paprxbem3d aprx, real accur)
{
  assert(accur > 0.0);

  uninit_aca_bem3d(aprx);

  aprx->accur_aca = accur;
}

void
build_bem3d_cube_quadpoints(pcbem3d bem, const real a[3], const real b[3],
			    const real delta, real(**Z)[3], real(**N)[3])
{
  pcaprxbem3d aprx = bem->aprx;
  const real *t = aprx->t_green;
  const real *w = aprx->w_green;
  const uint ml = aprx->ml_green;
  const uint k2 = aprx->k_green * 0.5;

  uint      gammanumber, mu1, mu2, nu;
  real      velo, factor, u, v;

  *N = (real(*)[3]) allocreal(3 * k2);
  *Z = (real(*)[3]) allocreal(3 * k2);

  nu = 0;

  for (gammanumber = 1; gammanumber <= 6; gammanumber++) {

    switch (gammanumber) {
    case 1:
      velo = (b[0] - a[0] + 2.0 * delta) * (b[1] - a[1] + 2.0 * delta) * 0.25;
      break;
    case 2:
      velo = (b[1] - a[1] + 2.0 * delta) * (b[2] - a[2] + 2.0 * delta) * 0.25;
      break;
    case 3:
      velo = (b[0] - a[0] + 2.0 * delta) * (b[1] - a[1] + 2.0 * delta) * 0.25;
      break;
    case 4:
      velo = (b[1] - a[1] + 2.0 * delta) * (b[2] - a[2] + 2.0 * delta) * 0.25;
      break;
    case 5:
      velo = (b[0] - a[0] + 2.0 * delta) * (b[2] - a[2] + 2.0 * delta) * 0.25;
      break;
    case 6:
      velo = (b[0] - a[0] + 2.0 * delta) * (b[2] - a[2] + 2.0 * delta) * 0.25;
      break;
    default:
      printf("ERROR: unknown gammanumber!\n");
      exit(0);
      break;
    }

    for (mu1 = 0; mu1 < ml; mu1++) {
      u = 0.5 * t[mu1];
      for (mu2 = 0; mu2 < ml; mu2++) {
	v = 0.5 * t[mu2];
	factor = velo * w[mu1] * w[mu2];

	(*N)[nu][0] = 0.0;
	(*N)[nu][1] = 0.0;
	(*N)[nu][2] = 0.0;

	switch (gammanumber) {
	case 1:
	  (*Z)[nu][0] =
	    (0.5 * (b[0] + a[0])) + (b[0] - a[0] + 2.0 * delta) * u;
	  (*Z)[nu][1] =
	    (0.5 * (b[1] + a[1])) + (b[1] - a[1] + 2.0 * delta) * v;
	  (*Z)[nu][2] = a[2] - delta;
	  (*N)[nu][2] = -factor;
	  break;
	case 2:
	  (*Z)[nu][0] = b[0] + delta;
	  (*Z)[nu][1] =
	    (0.5 * (b[1] + a[1])) + (b[1] - a[1] + 2.0 * delta) * u;
	  (*Z)[nu][2] =
	    (0.5 * (b[2] + a[2])) + (b[2] - a[2] + 2.0 * delta) * v;
	  (*N)[nu][0] = factor;
	  break;
	case 3:
	  (*Z)[nu][0] =
	    (0.5 * (b[0] + a[0])) + (b[0] - a[0] + 2.0 * delta) * u;
	  (*Z)[nu][1] =
	    (0.5 * (b[1] + a[1])) + (b[1] - a[1] + 2.0 * delta) * v;
	  (*Z)[nu][2] = b[2] + delta;
	  (*N)[nu][2] = factor;
	  break;
	case 4:
	  (*Z)[nu][0] = a[0] - delta;
	  (*Z)[nu][1] =
	    (0.5 * (b[1] + a[1])) + (b[1] - a[1] + 2.0 * delta) * u;
	  (*Z)[nu][2] =
	    (0.5 * (b[2] + a[2])) + (b[2] - a[2] + 2.0 * delta) * v;
	  (*N)[nu][0] = -factor;
	  break;
	case 5:
	  (*Z)[nu][0] =
	    (0.5 * (b[0] + a[0])) + (b[0] - a[0] + 2.0 * delta) * u;
	  (*Z)[nu][1] = a[1] - delta;
	  (*Z)[nu][2] =
	    (0.5 * (b[2] + a[2])) + (b[2] - a[2] + 2.0 * delta) * v;
	  (*N)[nu][1] = -factor;
	  break;
	case 6:
	  (*Z)[nu][0] =
	    (0.5 * (b[0] + a[0])) + (b[0] - a[0] + 2.0 * delta) * u;
	  (*Z)[nu][1] = b[1] + delta;
	  (*Z)[nu][2] =
	    (0.5 * (b[2] + a[2])) + (b[2] - a[2] + 2.0 * delta) * v;
	  (*N)[nu][1] = factor;
	  break;
	}

	nu++;
      }
    }
  }

  assert(nu == k2);
}

static void
assemble_interpoints3d_array(pcbem3d bem, pcreal bmin, pcreal bmax,
			     real(*X)[3])
{
  pcaprxbem3d aprx = bem->aprx;
  real     *x = aprx->x_inter;
  uint      m = aprx->m_inter;
  uint      k = aprx->k_inter;
  real      ax = bmin[0];
  real      bx = bmax[0];
  real      ay = bmin[1];
  real      by = bmax[1];
  real      az = bmin[2];
  real      bz = bmax[2];

  real      cx, dx, cy, dy, cz, dz;
  uint      i, j, l, index;

  if (bx - ax < INTERPOLATION_EPS_BEM3D) {
    bx += INTERPOLATION_EPS_BEM3D;
    ax -= INTERPOLATION_EPS_BEM3D;
  }
  if (by - ay < INTERPOLATION_EPS_BEM3D) {
    by += INTERPOLATION_EPS_BEM3D;
    ay -= INTERPOLATION_EPS_BEM3D;
  }
  if (bz - az < INTERPOLATION_EPS_BEM3D) {
    bz += INTERPOLATION_EPS_BEM3D;
    az -= INTERPOLATION_EPS_BEM3D;
  }

  /*
   * Tschebyscheff points
   */

  cx = (bx + ax) * 0.5;
  dx = (bx - ax) * 0.5;
  cy = (by + ay) * 0.5;
  dy = (by - ay) * 0.5;
  cz = (bz + az) * 0.5;
  dz = (bz - az) * 0.5;

  index = 0;
  for (i = 0; i < m; ++i) {
    for (j = 0; j < m; ++j) {
      for (l = 0; l < m; ++l) {
	X[index][0] = cx + dx * x[i];
	X[index][1] = cy + dy * x[j];
	X[index][2] = cz + dz * x[l];
	index++;
      }
    }
  }

  assert(index == k);
}

static void
assemble_interpoints3d_realavector(pcbem3d bem, pcreal bmin,
				   pcreal bmax, prealavector px,
				   prealavector py, prealavector pz)
{
  pcaprxbem3d aprx = bem->aprx;
  real     *x = aprx->x_inter;
  uint      m = aprx->m_inter;
  real      ax = bmin[0];
  real      bx = bmax[0];
  real      ay = bmin[1];
  real      by = bmax[1];
  real      az = bmin[2];
  real      bz = bmax[2];

  real      cx, dx, cy, dy, cz, dz;
  uint      i;

  if (bx - ax < INTERPOLATION_EPS_BEM3D) {
    bx += INTERPOLATION_EPS_BEM3D;
    ax -= INTERPOLATION_EPS_BEM3D;
  }
  if (by - ay < INTERPOLATION_EPS_BEM3D) {
    by += INTERPOLATION_EPS_BEM3D;
    ay -= INTERPOLATION_EPS_BEM3D;
  }
  if (bz - az < INTERPOLATION_EPS_BEM3D) {
    bz += INTERPOLATION_EPS_BEM3D;
    az -= INTERPOLATION_EPS_BEM3D;
  }

  /*
   * Tschebyscheff points
   */

  cx = (bx + ax) * 0.5;
  dx = (bx - ax) * 0.5;
  cy = (by + ay) * 0.5;
  dy = (by - ay) * 0.5;
  cz = (bz + az) * 0.5;
  dz = (bz - az) * 0.5;

  for (i = 0; i < m; ++i) {
    px->v[i] = cx + dx * x[i];
    py->v[i] = cy + dy * x[i];
    pz->v[i] = cz + dz * x[i];
  }
}

static void
assemble_bem3d_inter_row_rkmatrix(pccluster rc, uint rname,
				  pccluster cc, uint cname, pcbem3d bem,
				  prkmatrix R)
{
  paprxbem3d aprx = bem->aprx;
  pkernelbem3d kernels = bem->kernels;
  pamatrix  A = &R->A;
  pamatrix  B = &R->B;
  const uint rows = rc->size;
  const uint cols = cc->size;
  const uint m = aprx->m_inter;
  const uint k = aprx->k_inter;
  real(*z)[3];
  prealavector px, py, pz;

  (void) rname;
  (void) cname;

  z = (real(*)[3]) allocreal((size_t) (3 * k));
  px = new_realavector(m);
  py = new_realavector(m);
  pz = new_realavector(m);

  assemble_interpoints3d_array(bem, rc->bmin, rc->bmax, z);
  assemble_interpoints3d_realavector(bem, rc->bmin, rc->bmax, px, py, pz);

  resize_rkmatrix(R, rows, cols, k);

  kernels->kernel_col(cc->idx, (const real(*)[3]) z, bem, B);
  conjugate_amatrix(B);
  kernels->lagrange_row(rc->idx, px, py, pz, bem, A);

  del_realavector(px);
  del_realavector(py);
  del_realavector(pz);
  freemem(z);
}

static void
assemble_bem3d_inter_col_rkmatrix(pccluster rc, uint rname,
				  pccluster cc, uint cname, pcbem3d bem,
				  prkmatrix R)
{
  paprxbem3d aprx = bem->aprx;
  pkernelbem3d kernels = bem->kernels;
  pamatrix  A = &R->A;
  pamatrix  B = &R->B;
  const uint rows = rc->size;
  const uint cols = cc->size;
  const uint m = aprx->m_inter;
  const uint k = aprx->k_inter;
  real(*z)[3];
  prealavector px, py, pz;

  (void) rname;
  (void) cname;

  z = (real(*)[3]) allocreal((size_t) (3 * k));
  px = new_realavector(m);
  py = new_realavector(m);
  pz = new_realavector(m);

  assemble_interpoints3d_array(bem, cc->bmin, cc->bmax, z);
  assemble_interpoints3d_realavector(bem, cc->bmin, cc->bmax, px, py, pz);

  resize_rkmatrix(R, rows, cols, k);

  kernels->fundamental_row(rc->idx, (const real(*)[3]) z, bem, A);
  kernels->lagrange_col(cc->idx, px, py, pz, bem, B);

  del_realavector(px);
  del_realavector(py);
  del_realavector(pz);
  freemem(z);
}

static void
assemble_bem3d_inter_mixed_rkmatrix(pccluster rc, uint rname,
				    pccluster cc, uint cname, pcbem3d bem,
				    prkmatrix R)
{
  if (getdiam_2_cluster(rc) < getdiam_2_cluster(cc)) {
    assemble_bem3d_inter_row_rkmatrix(rc, rname, cc, cname, bem, R);
  }
  else {
    assemble_bem3d_inter_col_rkmatrix(rc, rname, cc, cname, bem, R);
  }
}

static void
assemble_bem3d_green_row_rkmatrix(pccluster rc, uint rname,
				  pccluster cc, uint cname, pcbem3d bem,
				  prkmatrix R)
{
  paprxbem3d aprx = bem->aprx;
  pkernelbem3d kernels = bem->kernels;
  pamatrix  A = &R->A;
  pamatrix  B = &R->B;
  real      delta = aprx->delta_green;
  real     *a = rc->bmin;
  real     *b = rc->bmax;

  pamatrix  T;
  real(*Z)[3], (*N)[3];
  real      diam;
  uint      nu, k2;

  (void) rname;
  (void) cname;

  k2 = 0.5 * aprx->k_green;
  diam =
    (aprx->quadpoints == build_bem3d_cube_quadpoints) ?
    getdiam_max_cluster(rc) : getdiam_2_cluster(rc);
  delta = delta * diam;

  T = new_amatrix(0, 0);

  resize_rkmatrix(R, rc->size, cc->size, 2 * k2);

  aprx->quadpoints(bem, a, b, delta, &Z, &N);

  init_sub_amatrix(T, A, rc->size, 0, k2, 0);
  kernels->fundamental_row(rc->idx, (const real(*)[3]) Z, (pcbem3d) bem, T);
  uninit_amatrix(T);

  init_sub_amatrix(T, B, cc->size, 0, k2, 0);
  kernels->dnz_kernel_col(cc->idx, (const real(*)[3]) Z, (const real(*)[3]) N,
			  (pcbem3d) bem, T);
  uninit_amatrix(T);

  for (nu = 0; nu < k2; ++nu) {
    N[nu][0] *= -1.0;
    N[nu][1] *= -1.0;
    N[nu][2] *= -1.0;
  }

  init_sub_amatrix(T, A, rc->size, 0, k2, k2);
  kernels->dnz_fundamental_row(rc->idx, (const real(*)[3]) Z,
			       (const real(*)[3]) N, (pcbem3d) bem, T);
  uninit_amatrix(T);

  init_sub_amatrix(T, B, cc->size, 0, k2, k2);
  kernels->kernel_col(cc->idx, (const real(*)[3]) Z, (pcbem3d) bem, T);
  uninit_amatrix(T);

  conjugate_amatrix(B);

  del_amatrix(T);
  freemem(Z);
  freemem(N);

}

static void
assemble_bem3d_green_col_rkmatrix(pccluster rc, uint rname,
				  pccluster cc, uint cname, pcbem3d bem,
				  prkmatrix R)
{
  paprxbem3d aprx = bem->aprx;
  pkernelbem3d kernels = bem->kernels;
  pamatrix  A = &R->A;
  pamatrix  B = &R->B;
  real      delta = aprx->delta_green;
  real     *a = cc->bmin;
  real     *b = cc->bmax;

  pamatrix  T;
  real(*Z)[3], (*N)[3];
  real      diam;
  uint      k2, nu;

  (void) rname;
  (void) cname;

  k2 = 0.5 * aprx->k_green;
  diam =
    (aprx->quadpoints == build_bem3d_cube_quadpoints) ?
    getdiam_max_cluster(cc) : getdiam_2_cluster(cc);
  delta = delta * diam;

  T = new_amatrix(0, 0);

  resize_rkmatrix(R, rc->size, cc->size, 2 * k2);

  aprx->quadpoints(bem, a, b, delta, &Z, &N);

  init_sub_amatrix(T, B, cc->size, 0, k2, 0);
  kernels->kernel_col(cc->idx, (const real(*)[3]) Z, (pcbem3d) bem, T);
  uninit_amatrix(T);

  init_sub_amatrix(T, A, rc->size, 0, k2, 0);
  kernels->dnz_fundamental_row(rc->idx, (const real(*)[3]) Z,
			       (const real(*)[3]) N, (pcbem3d) bem, T);
  uninit_amatrix(T);

  for (nu = 0; nu < k2; ++nu) {
    N[nu][0] *= -1.0;
    N[nu][1] *= -1.0;
    N[nu][2] *= -1.0;
  }

  init_sub_amatrix(T, B, cc->size, 0, k2, k2);
  kernels->dnz_kernel_col(cc->idx, (const real(*)[3]) Z, (const real(*)[3]) N,
			  (pcbem3d) bem, T);
  uninit_amatrix(T);

  init_sub_amatrix(T, A, rc->size, 0, k2, k2);
  kernels->fundamental_row(rc->idx, (const real(*)[3]) Z, (pcbem3d) bem, T);
  uninit_amatrix(T);

  conjugate_amatrix(B);

  del_amatrix(T);
  freemem(Z);
  freemem(N);
}

static void
assemble_bem3d_green_mixed_rkmatrix(pccluster rc, uint rname,
				    pccluster cc, uint cname, pcbem3d bem,
				    prkmatrix R)
{
  paprxbem3d aprx = bem->aprx;
  real      diamt, diams;

  diamt =
    (aprx->quadpoints == build_bem3d_cube_quadpoints) ?
    getdiam_max_cluster(rc) : getdiam_2_cluster(rc);
  diams =
    (aprx->quadpoints == build_bem3d_cube_quadpoints) ?
    getdiam_max_cluster(cc) : getdiam_2_cluster(cc);

  if (diamt < diams) {
    assemble_bem3d_green_row_rkmatrix(rc, rname, cc, cname, bem, R);
  }
  else {
    assemble_bem3d_green_col_rkmatrix(rc, rname, cc, cname, bem, R);
  }
}

static void
assemble_row_greencluster3d(pcbem3d bem, pgreencluster3d gc)
{
  paprxbem3d aprx = bem->aprx;
  pkernelbem3d kernels = bem->kernels;
  pccluster c = gc->t;
  uint      rows = c->size;
  real      eps = aprx->accur_aca;
  uint      rank = aprx->k_green;
  real      delta = aprx->delta_green;
  uint      k = aprx->k_green;
  real     *a = c->bmin;
  real     *b = c->bmax;

  prkmatrix R;
  pamatrix  A_t, RC, T;
  real(*Z)[3], (*N)[3];
  real      diam;
  uint     *xi;
  uint      i, k2;

  k2 = 0.5 * rank;
  diam =
    (aprx->quadpoints == build_bem3d_cube_quadpoints) ?
    getdiam_max_cluster(c) : getdiam_2_cluster(c);
  delta = delta * diam;

  A_t = new_amatrix(rows, k);

  T = new_amatrix(0, 0);

  aprx->quadpoints(bem, a, b, delta, &Z, &N);

  init_sub_amatrix(T, A_t, c->size, 0, k2, k2);
  kernels->dnz_fundamental_row(c->idx, (const real(*)[3]) Z,
			       (const real(*)[3]) N, (pcbem3d) bem, T);
  uninit_amatrix(T);

  init_sub_amatrix(T, A_t, c->size, 0, k2, 0);
  kernels->fundamental_row(c->idx, (const real(*)[3]) Z, (pcbem3d) bem, T);
  uninit_amatrix(T);

  R = new_rkmatrix(rows, k, 0);
  decomp_fullaca_rkmatrix(A_t, eps, &xi, NULL, R);
  rank = R->k;

  RC = new_amatrix(rank, rank);
  copy_lower_aca_amatrix(true, &R->A, xi, RC);
  triangularsolve_amatrix(true, true, true, RC, true, &R->A);

  gc->xi = xi;
  gc->xihat = allocuint(rank);
  for (i = 0; i < rank; ++i) {
    gc->xihat[i] = c->idx[xi[i]];
  }

  resize_amatrix(gc->V, rows, rank);
  copy_amatrix(false, &R->A, gc->V);

  del_rkmatrix(R);
  del_amatrix(RC);
  del_amatrix(A_t);
  del_amatrix(T);
  freemem(Z);
  freemem(N);
}

static void
assemble_bem3d_greenhybrid_row_rkmatrix(pccluster rc, uint rname,
					pccluster cc, uint cname, pcbem3d bem,
					prkmatrix R)
{
  pparbem3d par = bem->par;
  pamatrix  A = &R->A;
  pamatrix  B = &R->B;
  uint      rows = rc->size;
  uint      cols = cc->size;

  pamatrix  V;
  pgreencluster3d grc;
  uint     *xihat;
  uint      rank;

  (void) cname;

#ifdef USE_OPENMP
#pragma omp critical(greenhybrid)
#endif
  {
    grc = par->grcn[rname];

    if (grc == NULL) {
      grc = par->grcn[rname] = new_greencluster3d(rc);
      assemble_row_greencluster3d(bem, grc);
    }
  }

  V = grc->V;
  rank = V->cols;
  xihat = grc->xihat;

  resize_rkmatrix(R, rows, cols, rank);

  copy_amatrix(false, V, A);
  bem->nearfield_far(xihat, cc->idx, bem, true, B);
}

static void
assemble_col_greencluster3d(pcbem3d bem, pgreencluster3d gc)
{
  paprxbem3d aprx = bem->aprx;
  pkernelbem3d kernels = bem->kernels;
  pccluster c = gc->t;
  uint      cols = c->size;
  real      eps = aprx->accur_aca;
  uint      rank = aprx->k_green;
  real      delta = aprx->delta_green;
  real     *a = c->bmin;
  real     *b = c->bmax;

  prkmatrix R;
  pamatrix  A_t, RC, T;
  real(*Z)[3], (*N)[3];
  uint     *xi;
  real      diam;
  uint      i, k2;

  k2 = 0.5 * rank;
  diam =
    (aprx->quadpoints == build_bem3d_cube_quadpoints) ?
    getdiam_max_cluster(c) : getdiam_2_cluster(c);
  delta = delta * diam;

  A_t = new_amatrix(cols, 2 * k2);
  T = new_amatrix(0, 0);

  aprx->quadpoints(bem, a, b, delta, &Z, &N);

  init_sub_amatrix(T, A_t, c->size, 0, k2, k2);
  kernels->dnz_kernel_col(c->idx, (const real(*)[3]) Z, (const real(*)[3]) N,
			  (pcbem3d) bem, T);
  uninit_amatrix(T);

  init_sub_amatrix(T, A_t, c->size, 0, k2, 0);
  kernels->kernel_col(c->idx, (const real(*)[3]) Z, (pcbem3d) bem, T);
  uninit_amatrix(T);

  R = new_rkmatrix(cols, rank, 0);
  decomp_fullaca_rkmatrix(A_t, eps, &xi, NULL, R);
  rank = R->k;

  RC = new_amatrix(rank, rank);
  copy_lower_aca_amatrix(true, &R->A, xi, RC);
  triangularsolve_amatrix(true, true, true, RC, true, &R->A);

  gc->xi = xi;
  gc->xihat = allocuint(rank);
  for (i = 0; i < rank; ++i) {
    gc->xihat[i] = c->idx[xi[i]];
  }

  resize_amatrix(gc->V, cols, rank);
  copy_amatrix(false, &R->A, gc->V);

  del_rkmatrix(R);
  del_amatrix(RC);
  del_amatrix(A_t);
  del_amatrix(T);
  freemem(Z);
  freemem(N);
}

static void
assemble_bem3d_greenhybrid_col_rkmatrix(pccluster rc, uint rname,
					pccluster cc, uint cname, pcbem3d bem,
					prkmatrix R)
{
  pparbem3d par = bem->par;
  pamatrix  A = &R->A;
  pamatrix  B = &R->B;
  uint      rows = rc->size;
  uint      cols = cc->size;

  pamatrix  V;
  pgreencluster3d gcc;
  uint     *xihat;
  uint      rank;

  (void) rname;

#ifdef USE_OPENMP
#pragma omp critical(greenhybrid)
#endif
  {
    gcc = par->gccn[cname];

    if (gcc == NULL) {
      gcc = par->gccn[cname] = new_greencluster3d(cc);
      assemble_col_greencluster3d(bem, gcc);
    }
  }

  V = gcc->V;
  rank = V->cols;
  xihat = gcc->xihat;

  resize_amatrix(A, rows, rank);
  resize_amatrix(B, cols, rank);
  R->k = rank;

  bem->nearfield_far(rc->idx, xihat, bem, false, A);
  copy_amatrix(false, V, B);
  conjugate_amatrix(B);
}

static void
assemble_bem3d_greenhybrid_mixed_rkmatrix(pccluster rc, uint rname,
					  pccluster cc, uint cname,
					  pcbem3d bem, prkmatrix R)
{
  pparbem3d par = bem->par;
  pamatrix  A = &R->A;
  pamatrix  B = &R->B;
  uint      rows = rc->size;
  uint      cols = cc->size;

  pamatrix  V, W;
  pgreencluster3d grc, gcc;
  uint     *xihatV, *xihatW;
  uint      rankV, rankW;

#ifdef USE_OPENMP
#pragma omp critical(greenhybrid)
#endif
  {
    grc = par->grcn[rname];

    if (grc == NULL) {
      grc = par->grcn[rname] = new_greencluster3d(rc);
      assemble_row_greencluster3d(bem, grc);
    }
  }

#ifdef USE_OPENMP
#pragma omp critical(greenhybrid)
#endif
  {
    gcc = par->gccn[cname];

    if (gcc == NULL) {
      gcc = par->gccn[cname] = new_greencluster3d(cc);
      assemble_col_greencluster3d(bem, gcc);
    }
  }

  rankV = grc->V->cols;
  rankW = gcc->V->cols;

  if (cols * rankV <= rows * rankW) {
    V = grc->V;
    xihatV = grc->xihat;

    resize_amatrix(A, rows, rankV);
    resize_amatrix(B, cols, rankV);
    R->k = rankV;

    copy_amatrix(false, V, A);
    bem->nearfield_far(xihatV, cc->idx, bem, true, B);
  }
  else {
    W = gcc->V;
    xihatW = gcc->xihat;

    resize_amatrix(A, rows, rankW);
    resize_amatrix(B, cols, rankW);
    R->k = rankW;

    bem->nearfield_far(rc->idx, xihatW, bem, false, A);
    copy_amatrix(false, W, B);
  }

}

static void
assemble_bem3d_ACA_rkmatrix(pccluster rc, uint rname, pccluster cc,
			    uint cname, pcbem3d bem, prkmatrix R)
{
  paprxbem3d aprx = bem->aprx;
  const real accur = aprx->accur_aca;
  const uint *ridx = rc->idx;
  const uint *cidx = cc->idx;
  const uint rows = rc->size;
  const uint cols = cc->size;

  pamatrix  G;

  (void) rname;
  (void) cname;

  G = new_amatrix(rows, cols);
  bem->nearfield_far(ridx, cidx, bem, false, G);

  decomp_fullaca_rkmatrix(G, accur, NULL, NULL, R);

  del_amatrix(G);
}

static void
assemble_bem3d_PACA_rkmatrix(pccluster rc, uint rname, pccluster cc,
			     uint cname, pcbem3d bem, prkmatrix R)
{
  paprxbem3d aprx = bem->aprx;
  const real accur = aprx->accur_aca;
  matrixentry_t entry = (matrixentry_t) bem->nearfield_far;
  const uint *ridx = rc->idx;
  const uint *cidx = cc->idx;
  const uint rows = rc->size;
  const uint cols = cc->size;

  (void) rname;
  (void) cname;

  decomp_partialaca_rkmatrix(entry, (void *) bem, ridx, rows, cidx, cols,
			     accur, NULL, NULL, R);
}

static void
assemble_bem3d_HCA_rkmatrix(pccluster rc, uint rname, pccluster cc,
			    uint cname, pcbem3d bem, prkmatrix R)
{
  paprxbem3d aprx = bem->aprx;
  pkernelbem3d kernels = bem->kernels;
  const real accur = aprx->accur_aca;
  const uint k = aprx->k_inter;
  const uint *ridx = rc->idx;
  const uint *cidx = cc->idx;
  const uint rows = rc->size;
  const uint cols = cc->size;

  prkmatrix R2;
  pamatrix  S, C, D;
  real(*IT)[3], (*IS)[3], (*ITk)[3], (*ISk)[3];
  uint     *I_k, *J_k;
  uint      i, j, rank;

  (void) rname;
  (void) cname;

  IT = (real(*)[3]) allocreal(3 * k);
  IS = (real(*)[3]) allocreal(3 * k);

  assemble_interpoints3d_array(bem, rc->bmin, rc->bmax, IT);
  assemble_interpoints3d_array(bem, cc->bmin, cc->bmax, IS);

  S = new_amatrix(k, k);
  R2 = new_rkmatrix(k, k, 0);

  kernels->fundamental(bem, (const real(*)[3]) IT, (const real(*)[3]) IS, S);
  decomp_fullaca_rkmatrix(S, accur, &I_k, &J_k, R2);
  rank = R2->k;
  C = &R2->A;
  D = &R2->B;

  conjugate_amatrix(D);

  ITk = (real(*)[3]) allocreal(3 * rank);
  ISk = (real(*)[3]) allocreal(3 * rank);

  for (i = 0; i < rank; ++i) {
    for (j = 0; j < 3; ++j) {
      ITk[i][j] = IT[I_k[i]][j];
      ISk[i][j] = IS[J_k[i]][j];
    }
  }

  resize_amatrix(S, rank, rank);

  copy_lower_aca_amatrix(true, C, I_k, S);
  copy_upper_aca_amatrix(false, D, J_k, S);

  resize_rkmatrix(R, rows, cols, rank);

  kernels->kernel_row(ridx, (const real(*)[3]) ISk, bem, &R->A);
  kernels->kernel_col(cidx, (const real(*)[3]) ITk, bem, &R->B);

  conjugate_amatrix(&R->B);

  if (rows < cols) {
    triangularsolve_amatrix(false, false, true, S, true, &R->A);
    triangularsolve_amatrix(true, true, true, S, true, &R->A);
  }
  else {
    triangularsolve_amatrix(true, true, false, S, true, &R->B);
    triangularsolve_amatrix(false, false, false, S, true, &R->B);
  }

  del_rkmatrix(R2);
  del_amatrix(S);
  freemem(IT);
  freemem(IS);
  freemem(ITk);
  freemem(ISk);
  freemem(I_k);
  freemem(J_k);
}

/****************************************************
 * H^2 Interpolation
 ****************************************************/

static void
assemble_bem3d_inter_row_clusterbasis(uint rname, pcbem3d bem)
{
  paprxbem3d aprx = bem->aprx;
  pkernelbem3d kernels = bem->kernels;
  pparbem3d par = bem->par;
  pclusterbasis rb = par->rbn[rname];
  pamatrix  V = &rb->V;
  pccluster t = rb->t;
  const uint m = aprx->m_inter;
  const uint k = aprx->k_inter;

  prealavector px, py, pz;

  px = new_realavector(m);
  py = new_realavector(m);
  pz = new_realavector(m);

  assemble_interpoints3d_realavector(bem, t->bmin, t->bmax, px, py, pz);

  resize_amatrix(V, t->size, k);
  rb->k = k;
  update_clusterbasis(rb);

  kernels->lagrange_row(t->idx, px, py, pz, bem, V);

  del_realavector(px);
  del_realavector(py);
  del_realavector(pz);
}

static void
assemble_bem3d_inter_transfer_row_clusterbasis(uint rname, pcbem3d bem)
{
  paprxbem3d aprx = bem->aprx;
  pparbem3d par = bem->par;
  pclusterbasis rb = par->rbn[rname];
  pccluster t = rb->t;
  uint      sons = t->sons;
  const uint m = aprx->m_inter;
  const uint k = aprx->k_inter;

  pamatrix  E;
  real(*X)[3];
  prealavector px, py, pz;
  uint      s;

  X = (real(*)[3]) allocmem(3 * k * sizeof(real));
  px = new_realavector(m);
  py = new_realavector(m);
  pz = new_realavector(m);

  assemble_interpoints3d_realavector(bem, t->bmin, t->bmax, px, py, pz);

  resize_clusterbasis(rb, k);

  for (s = 0; s < sons; ++s) {
    E = &rb->son[s]->E;

    assemble_interpoints3d_array(bem, rb->son[s]->t->bmin,
				 rb->son[s]->t->bmax, X);

    assemble_bem3d_lagrange_amatrix((const real(*)[3]) X, px, py, pz, bem, E);
  }

  del_realavector(px);
  del_realavector(py);
  del_realavector(pz);
  freemem(X);
}

static void
assemble_bem3d_inter_col_clusterbasis(uint cname, pcbem3d bem)
{
  paprxbem3d aprx = bem->aprx;
  pkernelbem3d kernels = bem->kernels;
  pparbem3d par = bem->par;
  pclusterbasis cb = par->cbn[cname];
  pamatrix  V = &cb->V;
  pccluster t = cb->t;
  const uint m = aprx->m_inter;
  const uint k = aprx->k_inter;

  prealavector px, py, pz;

  px = new_realavector(m);
  py = new_realavector(m);
  pz = new_realavector(m);

  assemble_interpoints3d_realavector(bem, t->bmin, t->bmax, px, py, pz);

  resize_amatrix(V, t->size, k);
  cb->k = k;
  update_clusterbasis(cb);

  kernels->lagrange_col(t->idx, px, py, pz, bem, V);

  del_realavector(px);
  del_realavector(py);
  del_realavector(pz);
}

static void
assemble_bem3d_inter_transfer_col_clusterbasis(uint cname, pcbem3d bem)
{
  paprxbem3d aprx = bem->aprx;
  pparbem3d par = bem->par;
  pclusterbasis cb = par->cbn[cname];
  pccluster t = cb->t;
  uint      sons = t->sons;
  const uint m = aprx->m_inter;
  const uint k = aprx->k_inter;

  pamatrix  E;
  real(*X)[3];
  prealavector px, py, pz;
  uint      s;

  X = (real(*)[3]) allocmem(3 * k * sizeof(real));
  px = new_realavector(m);
  py = new_realavector(m);
  pz = new_realavector(m);

  assemble_interpoints3d_realavector(bem, t->bmin, t->bmax, px, py, pz);

  resize_clusterbasis(cb, k);

  for (s = 0; s < sons; ++s) {
    E = &cb->son[s]->E;

    assemble_interpoints3d_array(bem, cb->son[s]->t->bmin,
				 cb->son[s]->t->bmax, X);

    assemble_bem3d_lagrange_amatrix((const real(*)[3]) X, px, py, pz, bem, E);
  }

  del_realavector(px);
  del_realavector(py);
  del_realavector(pz);
  freemem(X);
}

static void
assemble_bem3d_inter_uniform(uint rname, uint cname, uint bname, pcbem3d bem)
{
  pkernelbem3d kernels = bem->kernels;
  pparbem3d par = bem->par;
  puniform  U = par->h2n[bname]->u;
  pccluster rc = U->rb->t;
  pccluster cc = U->cb->t;
  const uint kr = U->rb->k;
  const uint kc = U->cb->k;
  pamatrix  S = &U->S;

  real(*xi_r)[3], (*xi_c)[3];

  (void) rname;
  (void) cname;

  resize_amatrix(S, kr, kc);

  xi_r = (real(*)[3]) allocreal(3 * kr);
  xi_c = (real(*)[3]) allocreal(3 * kc);

  assemble_interpoints3d_array(bem, rc->bmin, rc->bmax, xi_r);
  assemble_interpoints3d_array(bem, cc->bmin, cc->bmax, xi_c);

  kernels->fundamental(bem, (const real(*)[3]) xi_r, (const real(*)[3]) xi_c,
		       S);

  freemem(xi_r);
  freemem(xi_c);
}

static void
update_pivotelements_greenclusterbasis3d(pgreenclusterbasis3d * grbn,
					 pcclusterbasis cb, uint * I_t,
					 uint k)
{
  uint      sons = cb->sons;
  pgreenclusterbasis3d grb = *grbn;

  pgreenclusterbasis3d grb1;
  uint      i, j, rname1, size, rank;

  uint     *xi;

  if (sons > 0) {
    xi = allocuint(k);

    for (j = 0; j < k; ++j) {
      rname1 = 1;
      grb1 = grbn[rname1];
      size = 0;
      rank = 0;
      for (i = 0; i < sons; ++i) {
	grb1 = grbn[rname1];
	if (I_t[j] < rank + cb->son[i]->k) {
	  break;
	}

	rname1 += cb->son[i]->t->desc;
	rank += cb->son[i]->k;
	size += cb->son[i]->t->size;
      }
      xi[j] = grb1->xi[I_t[j] - rank] + size;
    }
  }
  else {
    xi = I_t;
  }

  grb->xi = xi;
  grb->xihat = allocuint(k);
  for (i = 0; i < k; ++i) {
    grb->xihat[i] = cb->t->idx[xi[i]];
  }
}

static uint *
collect_pivotelements_greenclusterbasis3d(pgreenclusterbasis3d * grbn,
					  uint * rank)
{
  pcclusterbasis cb = (*grbn)->cb;
  uint      sons = cb->sons;

  pgreenclusterbasis3d grb1;
  uint      size, rname1, i, j, r;

  uint     *idx;

  r = 0;
  for (i = 0; i < sons; ++i) {
    assert(cb->son[i]->k > 0);
    r += cb->son[i]->k;
  }

  idx = allocuint(r);

  r = 0;
  size = 0;
  rname1 = 1;
  for (i = 0; i < sons; ++i) {
    grb1 = grbn[rname1];
    for (j = 0; j < cb->son[i]->k; ++j) {
      idx[r + j] = cb->t->idx[size + grb1->xi[j]];
    }
    r += cb->son[i]->k;
    size += cb->son[i]->t->size;
    rname1 += cb->son[i]->t->desc;
  }

  *rank = r;

  return idx;
}

static void
assemble_bem3d_greenhybrid_leaf_row_clusterbasis(uint rname, pcbem3d bem)
{
  paprxbem3d aprx = bem->aprx;
  pkernelbem3d kernels = bem->kernels;
  pparbem3d par = bem->par;
  pclusterbasis rb = par->rbn[rname];
  pccluster c = rb->t;
  const uint rows = c->size;
  const real eps = aprx->accur_aca;
  real      delta = aprx->delta_green;
  uint      rank = aprx->k_green;

  prkmatrix R;
  pamatrix  A_t, T, RC;
  pgreenclusterbasis3d grb;
  real(*Z)[3], (*N)[3];
  uint     *xi;
  uint      k2;

  k2 = rank * 0.5;
  delta = aprx->delta_green
    * ((aprx->quadpoints == build_bem3d_cube_quadpoints) ?
       getdiam_max_cluster(c) : getdiam_2_cluster(c));

  aprx->quadpoints(bem, c->bmin, c->bmax, delta, &Z, &N);

  T = new_amatrix(0, 0);

  grb = par->grbn[rname];
  if (grb == NULL) {
    grb = par->grbn[rname] = new_greenclusterbasis3d(rb);
  }

  A_t = new_amatrix(rows, rank);

  init_sub_amatrix(T, A_t, rows, 0, k2, 0);
  kernels->fundamental_row(c->idx, (const real(*)[3]) Z, (pcbem3d) bem, T);
  uninit_amatrix(T);

  init_sub_amatrix(T, A_t, rows, 0, k2, k2);
  kernels->dnz_fundamental_row(c->idx, (const real(*)[3]) Z,
			       (const real(*)[3]) N, (pcbem3d) bem, T);
  uninit_amatrix(T);

  R = new_rkmatrix(rows, rank, 0);
  decomp_fullaca_rkmatrix(A_t, eps, &xi, NULL, R);
  rank = R->k;
  resize_amatrix(&rb->V, rows, rank);
  copy_amatrix(false, &R->A, &rb->V);

  rb->k = rank;
  update_clusterbasis(rb);

  RC = new_amatrix(rank, rank);
  copy_lower_aca_amatrix(true, &rb->V, xi, RC);
  triangularsolve_amatrix(true, true, true, RC, true, &rb->V);

  update_pivotelements_greenclusterbasis3d(par->grbn + rname, rb, xi, rank);

  del_rkmatrix(R);
  del_amatrix(A_t);
  del_amatrix(RC);
  del_amatrix(T);
  freemem(Z);
  freemem(N);
}

static void
assemble_bem3d_greenhybrid_transfer_row_clusterbasis(uint rname, pcbem3d bem)
{
  paprxbem3d aprx = bem->aprx;
  pkernelbem3d kernels = bem->kernels;
  pparbem3d par = bem->par;
  pclusterbasis rb = par->rbn[rname];
  pccluster c = rb->t;
  const real eps = aprx->accur_aca;
  real      delta = aprx->delta_green;
  uint      rank = aprx->k_green;
  uint      sons = c->sons;

  prkmatrix R;
  pamatrix  A_t, T, RC, E;
  pgreenclusterbasis3d grb;
  real(*Z)[3], (*N)[3];
  uint     *I_t, *idx;
  uint      i, k2, newrank;

  k2 = rank * 0.5;
  delta = aprx->delta_green
    * ((aprx->quadpoints == build_bem3d_cube_quadpoints) ?
       getdiam_max_cluster(c) : getdiam_2_cluster(c));

  aprx->quadpoints(bem, c->bmin, c->bmax, delta, &Z, &N);

  T = new_amatrix(0, 0);

  grb = par->grbn[rname];
  if (grb == NULL) {
    grb = par->grbn[rname] = new_greenclusterbasis3d(rb);
  }

  idx = collect_pivotelements_greenclusterbasis3d(par->grbn + rname, &rank);

  A_t = new_amatrix(rank, 2 * k2);

  init_sub_amatrix(T, A_t, rank, 0, k2, 0);
  kernels->fundamental_row(idx, (const real(*)[3]) Z, (pcbem3d) bem, T);
  uninit_amatrix(T);

  init_sub_amatrix(T, A_t, rank, 0, k2, k2);
  kernels->dnz_fundamental_row(idx, (const real(*)[3]) Z,
			       (const real(*)[3]) N, (pcbem3d) bem, T);
  uninit_amatrix(T);

  R = new_rkmatrix(rank, 2 * k2, 0);
  decomp_fullaca_rkmatrix(A_t, eps, &I_t, NULL, R);
  newrank = R->k;

  RC = new_amatrix(newrank, newrank);
  copy_lower_aca_amatrix(true, &R->A, I_t, RC);
  triangularsolve_amatrix(true, true, true, RC, true, &R->A);

  resize_clusterbasis(rb, newrank);

  rank = 0;
  for (i = 0; i < sons; ++i) {
    E = &rb->son[i]->E;
    init_sub_amatrix(T, &R->A, rb->son[i]->k, rank, newrank, 0);
    copy_amatrix(false, T, E);
    uninit_amatrix(T);
    rank += rb->son[i]->k;
  }

  /* save local and global row pivot indices. */
  update_pivotelements_greenclusterbasis3d(par->grbn + rname, rb, I_t,
					   newrank);

  freemem(idx);
  freemem(I_t);
  del_rkmatrix(R);
  del_amatrix(RC);
  ;
  del_amatrix(A_t);
  del_amatrix(T);
  freemem(Z);
  freemem(N);
}

static void
assemble_bem3d_greenhybrid_leaf_col_clusterbasis(uint cname, pcbem3d bem)
{
  paprxbem3d aprx = bem->aprx;
  pkernelbem3d kernels = bem->kernels;
  pparbem3d par = bem->par;
  pclusterbasis cb = par->cbn[cname];
  pccluster c = cb->t;
  const uint rows = c->size;
  const real eps = aprx->accur_aca;
  real      delta = aprx->delta_green;
  uint      rank = aprx->k_green;

  prkmatrix R;
  pamatrix  A_t, T, RC;
  pgreenclusterbasis3d gcb;
  real(*Z)[3], (*N)[3];
  uint     *xi;
  uint      k2;

  k2 = rank * 0.5;
  delta = aprx->delta_green
    * ((aprx->quadpoints == build_bem3d_cube_quadpoints) ?
       getdiam_max_cluster(c) : getdiam_2_cluster(c));

  aprx->quadpoints(bem, c->bmin, c->bmax, delta, &Z, &N);

  T = new_amatrix(0, 0);

  gcb = par->gcbn[cname];
  if (gcb == NULL) {
    gcb = par->gcbn[cname] = new_greenclusterbasis3d(cb);
  }

  A_t = new_amatrix(rows, rank);

  init_sub_amatrix(T, A_t, rows, 0, k2, 0);
  kernels->kernel_col(c->idx, (const real(*)[3]) Z, (pcbem3d) bem, T);
  uninit_amatrix(T);

  init_sub_amatrix(T, A_t, rows, 0, k2, k2);
  kernels->dnz_kernel_col(c->idx, (const real(*)[3]) Z, (const real(*)[3]) N,
			  (pcbem3d) bem, T);
  uninit_amatrix(T);

  R = new_rkmatrix(rows, rank, 0);
  decomp_fullaca_rkmatrix(A_t, eps, &xi, NULL, R);
  rank = R->k;
  resize_amatrix(&cb->V, rows, rank);
  copy_amatrix(false, &R->A, &cb->V);

  cb->k = rank;
  update_clusterbasis(cb);

  RC = new_amatrix(rank, rank);
  copy_lower_aca_amatrix(true, &cb->V, xi, RC);
  triangularsolve_amatrix(true, true, true, RC, true, &cb->V);
  conjugate_amatrix(&cb->V);

  update_pivotelements_greenclusterbasis3d(par->gcbn + cname, cb, xi, rank);

  del_rkmatrix(R);
  del_amatrix(A_t);
  del_amatrix(RC);
  del_amatrix(T);
  freemem(Z);
  freemem(N);
}

static void
assemble_bem3d_greenhybrid_transfer_col_clusterbasis(uint cname, pcbem3d bem)
{
  paprxbem3d aprx = bem->aprx;
  pkernelbem3d kernels = bem->kernels;
  pparbem3d par = bem->par;
  pclusterbasis cb = par->cbn[cname];
  pccluster c = cb->t;
  const real eps = aprx->accur_aca;
  real      delta = aprx->delta_green;
  uint      rank = aprx->k_green;
  uint      sons = c->sons;

  prkmatrix R;
  pamatrix  A_t, T, RC, E;
  pgreenclusterbasis3d gcb;
  real(*Z)[3], (*N)[3];
  uint     *I_t, *idx;
  uint      i, k2, newrank;

  k2 = rank * 0.5;
  delta = aprx->delta_green
    * ((aprx->quadpoints == build_bem3d_cube_quadpoints) ?
       getdiam_max_cluster(c) : getdiam_2_cluster(c));

  aprx->quadpoints(bem, c->bmin, c->bmax, delta, &Z, &N);

  T = new_amatrix(0, 0);

  gcb = par->gcbn[cname];
  if (gcb == NULL) {
    gcb = par->gcbn[cname] = new_greenclusterbasis3d(cb);
  }

  idx = collect_pivotelements_greenclusterbasis3d(par->gcbn + cname, &rank);

  A_t = new_amatrix(rank, 2 * k2);

  init_sub_amatrix(T, A_t, rank, 0, k2, 0);
  kernels->kernel_col(idx, (const real(*)[3]) Z, (pcbem3d) bem, T);
  uninit_amatrix(T);

  init_sub_amatrix(T, A_t, rank, 0, k2, k2);
  kernels->dnz_kernel_col(idx, (const real(*)[3]) Z, (const real(*)[3]) N,
			  (pcbem3d) bem, T);
  uninit_amatrix(T);

  R = new_rkmatrix(rank, 2 * k2, 0);
  decomp_fullaca_rkmatrix(A_t, eps, &I_t, NULL, R);
  newrank = R->k;

  RC = new_amatrix(newrank, newrank);
  copy_lower_aca_amatrix(true, &R->A, I_t, RC);
  triangularsolve_amatrix(true, true, true, RC, true, &R->A);
  conjugate_amatrix(&R->A);

  resize_clusterbasis(cb, newrank);

  rank = 0;
  for (i = 0; i < sons; ++i) {
    E = &cb->son[i]->E;
    init_sub_amatrix(T, &R->A, cb->son[i]->k, rank, newrank, 0);
    copy_amatrix(false, T, E);
    uninit_amatrix(T);
    rank += cb->son[i]->k;
  }

  /* save local and global row pivot indices. */
  update_pivotelements_greenclusterbasis3d(par->gcbn + cname, cb, I_t,
					   newrank);

  freemem(I_t);
  freemem(idx);
  del_rkmatrix(R);
  del_amatrix(RC);
  del_amatrix(A_t);
  del_amatrix(T);
  freemem(Z);
  freemem(N);
}

static void
assemble_bem3d_greenhybrid_uniform(uint rname, uint cname,
				   uint bname, pcbem3d bem)
{
  pparbem3d par = bem->par;
  puniform  U = par->h2n[bname]->u;
  const uint kr = U->rb->k;
  const uint kc = U->cb->k;
  pamatrix  S = &U->S;

  pgreenclusterbasis3d grb, gcb;
  uint     *xihatV, *xihatW;

  grb = par->grbn[rname];
  gcb = par->gcbn[cname];

  assert(grb != NULL);
  assert(gcb != NULL);

  xihatV = grb->xihat;
  xihatW = gcb->xihat;

  resize_amatrix(S, kr, kc);

  bem->nearfield_far(xihatV, xihatW, bem, false, S);
}

static void
assemble_bem3d_greenhybrid_ortho_leaf_row_clusterbasis(uint rname,
						       pcbem3d bem)
{
  paprxbem3d aprx = bem->aprx;
  pkernelbem3d kernels = bem->kernels;
  pparbem3d par = bem->par;
  pclusterbasis rb = par->rbn[rname];
  pccluster t = rb->t;
  pamatrix  V = &rb->V;
  const uint rows = t->size;
  const real eps = aprx->accur_aca;
  real      delta = aprx->delta_green;
  uint      rank = aprx->k_green;

  prkmatrix R;
  pamatrix  A_t, T, RC;
  pavector  tau;
  pgreenclusterbasis3d grb;
  real(*Z)[3], (*N)[3];
  uint     *xi;
  uint      k2;

  k2 = rank * 0.5;
  delta = aprx->delta_green
    * ((aprx->quadpoints == build_bem3d_cube_quadpoints) ?
       getdiam_max_cluster(t) : getdiam_2_cluster(t));

  aprx->quadpoints(bem, t->bmin, t->bmax, delta, &Z, &N);

  T = new_amatrix(0, 0);

  grb = par->grbn[rname];
  if (grb == NULL) {
    grb = par->grbn[rname] = new_greenclusterbasis3d(rb);
  }

  /* setup up A_t from green's formula. */
  A_t = new_amatrix(rows, rank);

  init_sub_amatrix(T, A_t, rows, 0, k2, 0);
  kernels->fundamental_row(t->idx, (const real(*)[3]) Z, (pcbem3d) bem, T);
  uninit_amatrix(T);

  init_sub_amatrix(T, A_t, rows, 0, k2, k2);
  kernels->dnz_fundamental_row(t->idx, (const real(*)[3]) Z,
			       (const real(*)[3]) N, (pcbem3d) bem, T);
  uninit_amatrix(T);

  /* Do ACA of A_t to obtain row pivots. */
  R = new_rkmatrix(rows, rank, 0);
  decomp_fullaca_rkmatrix(A_t, eps, &xi, NULL, R);
  rank = R->k;
  assert(rows >= rank);

  resize_clusterbasis(rb, rank);

  /* Compute \hat V_t = C_t * (R_t * C_t)^-1. */
  RC = new_amatrix(rank, rank);
  copy_lower_aca_amatrix(true, &R->A, xi, RC);
  triangularsolve_amatrix(true, true, true, RC, true, &R->A);

  /* Compute QR factorization of \hat V_t. */
  tau = new_avector(rank);
  resize_amatrix(V, rows, rank);
  grb->Qinv = new_amatrix(rank, rank);

  qrdecomp_amatrix(&R->A, tau);
  qrexpand_amatrix(&R->A, tau, V);
  copy_upper_amatrix(&R->A, false, grb->Qinv);

  /* save local and global row pivot indices. */
  update_pivotelements_greenclusterbasis3d(par->grbn + rname, rb, xi, rank);

  /* clean up */
  del_rkmatrix(R);
  del_amatrix(A_t);
  del_amatrix(RC);
  del_amatrix(T);
  del_avector(tau);
  freemem(Z);
  freemem(N);
}

static void
assemble_bem3d_greenhybrid_ortho_transfer_row_clusterbasis(uint rname,
							   pcbem3d bem)
{
  paprxbem3d aprx = bem->aprx;
  pkernelbem3d kernels = bem->kernels;
  pparbem3d par = bem->par;
  pclusterbasis rb = par->rbn[rname];
  pccluster t = rb->t;
  const real eps = aprx->accur_aca;
  real      delta = aprx->delta_green;
  uint      rank = aprx->k_green;
  uint      sons = t->sons;

  prkmatrix R;
  pamatrix  A_t, T, RC, E, Q;
  pavector  tau;
  pgreenclusterbasis3d grb, grb1;
  real(*Z)[3], (*N)[3];
  uint     *I_t, *idx;
  uint      i, k2, newrank, rname1;

  k2 = rank * 0.5;
  delta = aprx->delta_green
    * ((aprx->quadpoints == build_bem3d_cube_quadpoints) ?
       getdiam_max_cluster(t) : getdiam_2_cluster(t));

  aprx->quadpoints(bem, t->bmin, t->bmax, delta, &Z, &N);

  T = new_amatrix(0, 0);

  grb = par->grbn[rname];
  if (grb == NULL) {
    grb = par->grbn[rname] = new_greenclusterbasis3d(rb);
  }

  idx = collect_pivotelements_greenclusterbasis3d(par->grbn + rname, &rank);

  /* setup up A_t from green's formula with reduced row indices. */
  A_t = new_amatrix(rank, 2 * k2);

  init_sub_amatrix(T, A_t, rank, 0, k2, 0);
  kernels->fundamental_row(idx, (const real(*)[3]) Z, (pcbem3d) bem, T);
  uninit_amatrix(T);

  init_sub_amatrix(T, A_t, rank, 0, k2, k2);
  kernels->dnz_fundamental_row(idx, (const real(*)[3]) Z,
			       (const real(*)[3]) N, (pcbem3d) bem, T);
  uninit_amatrix(T);

  /* Do ACA of \tilde A_t. */
  R = new_rkmatrix(rank, 2 * k2, 0);
  decomp_fullaca_rkmatrix(A_t, eps, &I_t, NULL, R);
  newrank = R->k;
  assert(rank >= newrank);
  resize_clusterbasis(rb, newrank);

  /* Compute \hat V_t = C_t * (R_t * C_t)^-1. */
  RC = new_amatrix(newrank, newrank);
  copy_lower_aca_amatrix(true, &R->A, I_t, RC);
  triangularsolve_amatrix(true, true, true, RC, true, &R->A);

  /* Multiply with triangular matrices of the sons from the left. */
  rank = 0;
  rname1 = rname + 1;
  for (i = 0; i < sons; ++i) {
    grb1 = par->grbn[rname1];
    init_sub_amatrix(T, &R->A, rb->son[i]->k, rank, newrank, 0);
    triangulareval_amatrix(false, false, false, grb1->Qinv, false, T);
    uninit_amatrix(T);
    rank += rb->son[i]->k;
    rname1 += grb1->cb->t->desc;
  }

  /* Compute QR factorization of \hat V_t. */
  tau = new_avector(newrank);
  Q = new_amatrix(rank, newrank);
  grb->Qinv = new_amatrix(newrank, newrank);

  qrdecomp_amatrix(&R->A, tau);
  qrexpand_amatrix(&R->A, tau, Q);
  copy_upper_amatrix(&R->A, false, grb->Qinv);

  /* copy transfer matrices to their location. */
  rank = 0;
  for (i = 0; i < sons; ++i) {
    E = &rb->son[i]->E;
    init_sub_amatrix(T, Q, rb->son[i]->k, rank, newrank, 0);
    copy_amatrix(false, T, E);
    uninit_amatrix(T);
    rank += rb->son[i]->k;
  }

  /* save local and global row pivot indices. */
  update_pivotelements_greenclusterbasis3d(par->grbn + rname, rb, I_t,
					   newrank);

  /* clean up */
  del_rkmatrix(R);
  del_amatrix(RC);
  del_amatrix(A_t);
  del_amatrix(T);
  del_amatrix(Q);
  del_avector(tau);
  freemem(Z);
  freemem(N);
  freemem(idx);
  freemem(I_t);
}

static void
assemble_bem3d_greenhybrid_ortho_leaf_col_clusterbasis(uint cname,
						       pcbem3d bem)
{
  paprxbem3d aprx = bem->aprx;
  pkernelbem3d kernels = bem->kernels;
  pparbem3d par = bem->par;
  pclusterbasis cb = par->cbn[cname];
  pccluster t = cb->t;
  pamatrix  V = &cb->V;
  const uint rows = t->size;
  const real eps = aprx->accur_aca;
  real      delta = aprx->delta_green;
  uint      rank = aprx->k_green;

  prkmatrix R;
  pamatrix  A_t, T, RC;
  pavector  tau;
  pgreenclusterbasis3d gcb;
  real(*Z)[3], (*N)[3];
  uint     *xi;
  uint      k2;

  k2 = rank * 0.5;
  delta = aprx->delta_green
    * ((aprx->quadpoints == build_bem3d_cube_quadpoints) ?
       getdiam_max_cluster(t) : getdiam_2_cluster(t));

  aprx->quadpoints(bem, t->bmin, t->bmax, delta, &Z, &N);

  T = new_amatrix(0, 0);

  gcb = par->gcbn[cname];
  if (gcb == NULL) {
    gcb = par->gcbn[cname] = new_greenclusterbasis3d(cb);
  }

  /* setup up A_t from green's formula. */
  A_t = new_amatrix(rows, rank);

  init_sub_amatrix(T, A_t, rows, 0, k2, 0);
  kernels->kernel_col(t->idx, (const real(*)[3]) Z, (pcbem3d) bem, T);
  uninit_amatrix(T);

  init_sub_amatrix(T, A_t, rows, 0, k2, k2);
  kernels->dnz_kernel_col(t->idx, (const real(*)[3]) Z, (const real(*)[3]) N,
			  (pcbem3d) bem, T);
  uninit_amatrix(T);

  /* Do ACA of A_t to obtain row pivots. */
  R = new_rkmatrix(rows, rank, 0);
  decomp_fullaca_rkmatrix(A_t, eps, &xi, NULL, R);
  rank = R->k;
  assert(rows >= rank);
  resize_clusterbasis(cb, rank);

  /* Compute \hat V_t = C_t * (R_t * C_t)^-1. */
  RC = new_amatrix(rank, rank);
  copy_lower_aca_amatrix(true, &R->A, xi, RC);
  triangularsolve_amatrix(true, true, true, RC, true, &R->A);

  /* Compute QR factorization of \hat V_t. */
  tau = new_avector(rank);
  resize_amatrix(V, rows, rank);
  gcb->Qinv = new_amatrix(rank, rank);

  qrdecomp_amatrix(&R->A, tau);
  qrexpand_amatrix(&R->A, tau, V);
  copy_upper_amatrix(&R->A, false, gcb->Qinv);

  /* save local and global row pivot indices. */
  update_pivotelements_greenclusterbasis3d(par->gcbn + cname, cb, xi, rank);

  /* clean up */
  del_rkmatrix(R);
  del_amatrix(A_t);
  del_amatrix(RC);
  del_amatrix(T);
  del_avector(tau);
  freemem(Z);
  freemem(N);
}

static void
assemble_bem3d_greenhybrid_ortho_transfer_col_clusterbasis(uint cname,
							   pcbem3d bem)
{
  paprxbem3d aprx = bem->aprx;
  pkernelbem3d kernels = bem->kernels;
  pparbem3d par = bem->par;
  pclusterbasis cb = par->cbn[cname];
  pccluster t = cb->t;
  const real eps = aprx->accur_aca;
  real      delta = aprx->delta_green;
  uint      rank = aprx->k_green;
  uint      sons = t->sons;

  prkmatrix R;
  pamatrix  A_t, T, RC, E, Q;
  pavector  tau;
  pgreenclusterbasis3d gcb, gcb1;
  real(*Z)[3], (*N)[3];
  uint     *I_t, *idx;
  uint      i, k2, newrank, cname1;

  k2 = rank * 0.5;
  delta = aprx->delta_green
    * ((aprx->quadpoints == build_bem3d_cube_quadpoints) ?
       getdiam_max_cluster(t) : getdiam_2_cluster(t));

  aprx->quadpoints(bem, t->bmin, t->bmax, delta, &Z, &N);

  T = new_amatrix(0, 0);

  gcb = par->gcbn[cname];
  if (gcb == NULL) {
    gcb = par->gcbn[cname] = new_greenclusterbasis3d(cb);
  }

  idx = collect_pivotelements_greenclusterbasis3d(par->gcbn + cname, &rank);

  /* setup up A_t from green's formula with reduced row indices. */
  A_t = new_amatrix(rank, 2 * k2);

  init_sub_amatrix(T, A_t, rank, 0, k2, 0);
  kernels->kernel_col(idx, (const real(*)[3]) Z, (pcbem3d) bem, T);
  uninit_amatrix(T);

  init_sub_amatrix(T, A_t, rank, 0, k2, k2);
  kernels->dnz_kernel_col(idx, (const real(*)[3]) Z, (const real(*)[3]) N,
			  (pcbem3d) bem, T);
  uninit_amatrix(T);

  /* Do ACA of \tilde A_t. */
  R = new_rkmatrix(rank, 2 * k2, 0);
  decomp_fullaca_rkmatrix(A_t, eps, &I_t, NULL, R);
  newrank = R->k;
  assert(rank >= newrank);
  resize_clusterbasis(cb, newrank);

  /* Compute \hat V_t = C_t * (R_t * C_t)^-1. */
  RC = new_amatrix(newrank, newrank);
  copy_lower_aca_amatrix(true, &R->A, I_t, RC);
  triangularsolve_amatrix(true, true, true, RC, true, &R->A);

  rank = 0;
  cname1 = cname + 1;
  for (i = 0; i < sons; ++i) {
    gcb1 = par->gcbn[cname1];
    init_sub_amatrix(T, &R->A, cb->son[i]->k, rank, newrank, 0);
    triangulareval_amatrix(false, false, false, gcb1->Qinv, false, T);
    uninit_amatrix(T);
    rank += cb->son[i]->k;
    cname1 += gcb1->cb->t->desc;
  }

  /* Compute QR factorization of \hat V_t. */
  tau = new_avector(newrank);
  Q = new_amatrix(rank, newrank);
  gcb->Qinv = new_amatrix(newrank, newrank);

  qrdecomp_amatrix(&R->A, tau);
  qrexpand_amatrix(&R->A, tau, Q);
  copy_upper_amatrix(&R->A, false, gcb->Qinv);

  /* copy transfer matrices to their location. */
  rank = 0;
  for (i = 0; i < sons; ++i) {
    E = &cb->son[i]->E;
    init_sub_amatrix(T, Q, cb->son[i]->k, rank, newrank, 0);
    copy_amatrix(false, T, E);
    uninit_amatrix(T);
    rank += cb->son[i]->k;
  }

  /* save local and global row pivot indices. */
  update_pivotelements_greenclusterbasis3d(par->gcbn + cname, cb, I_t,
					   newrank);

  /* clean up */
  del_rkmatrix(R);
  del_amatrix(RC);
  del_amatrix(A_t);
  del_amatrix(T);
  del_amatrix(Q);
  del_avector(tau);
  freemem(Z);
  freemem(N);
  freemem(idx);
  freemem(I_t);
}

static void
assemble_bem3d_greenhybrid_ortho_uniform(uint rname, uint cname,
					 uint bname, pcbem3d bem)
{
  pparbem3d par = bem->par;
  puniform  U = par->h2n[bname]->u;
  const uint kr = U->rb->k;
  const uint kc = U->cb->k;
  pamatrix  S = &U->S;

  pgreenclusterbasis3d grb, gcb;
  uint     *xihatV, *xihatW;

  grb = par->grbn[rname];
  gcb = par->gcbn[cname];

  assert(grb != NULL);
  assert(gcb != NULL);

  xihatV = grb->xihat;
  xihatW = gcb->xihat;

  resize_amatrix(S, kr, kc);

  bem->nearfield_far(xihatV, xihatW, bem, false, S);

  triangulareval_amatrix(false, false, false, grb->Qinv, false, S);
  triangulareval_amatrix(false, false, false, gcb->Qinv, true, S);
}

/****************************************************
 * Directional H^2 Interpolation
 ****************************************************/

static void
assemble_bem3d_inter_nowave_leaf_row_dclusterbasis(uint rname, pcbem3d bem)
{
  pkernelbem3d kernels = bem->kernels;
  pparbem3d par = bem->par;
  pdclusterbasis rb = par->drbn[rname];
  pamatrix  V = &rb->V[0];
  pcdcluster t = rb->t;
  uint      m = bem->aprx->m_inter;

  prealavector px, py, pz;

  px = new_realavector(m);
  py = new_realavector(m);
  pz = new_realavector(m);

  assemble_interpoints3d_realavector(bem, t->bmin, t->bmax, px, py, pz);

  kernels->lagrange_row(t->idx, px, py, pz, bem, V);

  del_realavector(px);
  del_realavector(py);
  del_realavector(pz);
}

static void
assemble_bem3d_inter_nowave_leaf_col_dclusterbasis(uint cname, pcbem3d bem)
{
  pkernelbem3d kernels = bem->kernels;
  pparbem3d par = bem->par;
  pdclusterbasis cb = par->dcbn[cname];
  pamatrix  V = &cb->V[0];
  pcdcluster t = cb->t;
  uint      m = bem->aprx->m_inter;

  prealavector px, py, pz;

  px = new_realavector(m);
  py = new_realavector(m);
  pz = new_realavector(m);

  assemble_interpoints3d_realavector(bem, t->bmin, t->bmax, px, py, pz);

  kernels->lagrange_col(t->idx, px, py, pz, bem, V);

  del_realavector(px);
  del_realavector(py);
  del_realavector(pz);
}

/* Approximate \int_\Omega exp(i k <y, c_\iota>) L_j(y) \varphi_i(y) dy */
static void
assemble_bem3d_inter_wave_leaf_row_dclusterbasis(uint rname, pcbem3d bem)
{
  pkernelbem3d kernels = bem->kernels;
  pdclusterbasis rb = bem->par->drbn[rname];
  pcdcluster t = rb->t;
  uint      m = bem->aprx->m_inter;
  uint      k = bem->aprx->k_inter;
  pamatrix  V;
  uint      iota;
  prealavector px, py, pz;

  px = new_realavector(m);
  py = new_realavector(m);
  pz = new_realavector(m);

  assemble_interpoints3d_realavector(bem, t->bmin, t->bmax, px, py, pz);

  for (iota = 0; iota < rb->directions; iota++) {
    V = rb->V + iota;
    if (V->cols > 0) {
      assert(V->rows == t->size);
      assert(V->cols == k);

      kernels->lagrange_wave_row(t->idx, px, py, pz, t->dir[iota], bem, V);
    }
  }

  del_realavector(px);
  del_realavector(py);
  del_realavector(pz);
}

/* Approximate \int_\Omega exp(i k <y, c_\iota>) L_j(y) \varphi_i(y) dy */
static void
assemble_bem3d_inter_wave_leaf_col_dclusterbasis(uint cname, pcbem3d bem)
{
  pkernelbem3d kernels = bem->kernels;
  pdclusterbasis cb = bem->par->dcbn[cname];
  pcdcluster t = cb->t;
  uint      m = bem->aprx->m_inter;
  uint      k = bem->aprx->k_inter;
  pamatrix  V;
  uint      iota;
  prealavector px, py, pz;

  px = new_realavector(m);
  py = new_realavector(m);
  pz = new_realavector(m);

  assemble_interpoints3d_realavector(bem, t->bmin, t->bmax, px, py, pz);

  for (iota = 0; iota < cb->directions; iota++) {
    V = cb->V + iota;
    if (V->cols > 0) {
      assert(V->rows == t->size);
      assert(V->cols == k);

      kernels->lagrange_wave_col(t->idx, px, py, pz, t->dir[iota], bem, V);
    }
  }

  del_realavector(px);
  del_realavector(py);
  del_realavector(pz);
}

static void
assemble_bem3d_inter_nowave_transfer_row_dclusterbasis(uint rname,
						       pcbem3d bem)
{
  paprxbem3d aprx = bem->aprx;
  pparbem3d par = bem->par;
  pdclusterbasis rb = par->drbn[rname];
  pcdcluster t = rb->t;
  uint      sons = t->sons;
  uint      m = aprx->m_inter;
  uint      k = aprx->k_inter;

  pamatrix  E;
  real(*X)[3];
  prealavector px, py, pz;
  uint      s;

  X = (real(*)[3]) allocmem(3 * k * sizeof(real));
  px = new_realavector(m);
  py = new_realavector(m);
  pz = new_realavector(m);

  assemble_interpoints3d_realavector(bem, t->bmin, t->bmax, px, py, pz);

  for (s = 0; s < sons; ++s) {
    assert(rb->son[s]->directions == 1);
    assert(rb->son[s]->t->directions == 0);
    E = rb->E[s];

    if (E->cols > 0) {
      assemble_interpoints3d_array(bem, rb->son[s]->t->bmin,
				   rb->son[s]->t->bmax, X);

      assemble_bem3d_lagrange_amatrix((const real(*)[3]) X, px, py, pz, bem,
				      E);
    }
  }

  del_realavector(px);
  del_realavector(py);
  del_realavector(pz);
  freemem(X);
}

static void
assemble_bem3d_inter_nowave_transfer_col_dclusterbasis(uint cname,
						       pcbem3d bem)
{
  paprxbem3d aprx = bem->aprx;
  pparbem3d par = bem->par;
  pdclusterbasis cb = par->dcbn[cname];
  pcdcluster t = cb->t;
  uint      sons = t->sons;
  uint      m = aprx->m_inter;
  uint      k = aprx->k_inter;

  pamatrix  E;
  real(*X)[3];
  prealavector px, py, pz;
  uint      s;

  X = (real(*)[3]) allocmem(3 * k * sizeof(real));
  px = new_realavector(m);
  py = new_realavector(m);
  pz = new_realavector(m);

  assemble_interpoints3d_realavector(bem, t->bmin, t->bmax, px, py, pz);

  for (s = 0; s < sons; ++s) {
    assert(cb->son[s]->directions == 1);
    assert(cb->son[s]->t->directions == 0);
    E = cb->E[s];

    if (E->cols > 0) {
      assemble_interpoints3d_array(bem, cb->son[s]->t->bmin,
				   cb->son[s]->t->bmax, X);

      assemble_bem3d_lagrange_amatrix((const real(*)[3]) X, px, py, pz, bem,
				      E);
    }
  }

  del_realavector(px);
  del_realavector(py);
  del_realavector(pz);
  freemem(X);
}

static void
assemble_bem3d_inter_wave_transfer_row_dclusterbasis(uint rname, pcbem3d bem)
{
  pdclusterbasis rb = bem->par->drbn[rname];
  pcdcluster fc = rb->t;
  uint      ni = bem->aprx->m_inter;
  pcdcluster sc;
  pamatrix  E;
  uint      j, iota;
  prealavector px, py, pz;
  real(*X)[3];

  X = (real(*)[3]) allocreal(3 * ni * ni * ni);
  px = new_realavector(ni);
  py = new_realavector(ni);
  pz = new_realavector(ni);

  /* Find interpolation points for the father */

  assemble_interpoints3d_realavector(bem, fc->bmin, fc->bmax, px, py, pz);

  for (j = 0; j < rb->sons; j++) {
    sc = rb->son[j]->t;

    /* Find interpolation points for the son */
    assemble_interpoints3d_array(bem, sc->bmin, sc->bmax, X);

    if (rb->son[j]->t->directions == 0) {
      /* Father uses wave-Lagrange basis, son uses standard Lagrange */
      assert(rb->son[j]->directions == 1);

      for (iota = 0; iota < rb->directions; iota++) {
	E = rb->E[j] + iota;
	if (E->cols > 0) {
	  assert(E->rows > 0);
	  assert(E->rows == rb->k[iota]);
//          assert(E->cols == rb->son[j]->k[0]);                        /* Isn't true for direct orthogonalization */

	  assemble_bem3d_lagrange_wave_amatrix((const real(*)[3]) X, px, py,
					       pz, fc->dir[iota], bem, E);
	}
      }
    }
  }

  /* Clean up */
  del_realavector(px);
  del_realavector(py);
  del_realavector(pz);
  freemem(X);
}

static void
assemble_bem3d_inter_wave_transfer_col_dclusterbasis(uint cname, pcbem3d bem)
{
  pdclusterbasis cb = bem->par->dcbn[cname];
  pcdcluster fc = cb->t;
  uint      ni = bem->aprx->m_inter;
  pcdcluster sc;
  pamatrix  E;
  uint      j, iota;
  prealavector px, py, pz;
  real(*X)[3];

  X = (real(*)[3]) allocreal(3 * ni * ni * ni);
  px = new_realavector(ni);
  py = new_realavector(ni);
  pz = new_realavector(ni);

  /* Find interpolation points for the father */

  assemble_interpoints3d_realavector(bem, fc->bmin, fc->bmax, px, py, pz);

  for (j = 0; j < cb->sons; j++) {
    sc = cb->son[j]->t;

    /* Find interpolation points for the son */
    assemble_interpoints3d_array(bem, sc->bmin, sc->bmax, X);

    if (cb->son[j]->t->directions == 0) {
      /* Father uses wave-Lagrange basis, son uses standard Lagrange */
      assert(cb->son[j]->directions == 1);

      for (iota = 0; iota < cb->directions; iota++) {
	E = cb->E[j] + iota;
	if (E->cols > 0) {
	  assert(E->rows > 0);
	  assert(E->rows == cb->k[iota]);
//          assert(E->cols == cb->son[j]->k[0]);                  /* Isn't true for direct orthogonalization */

	  assemble_bem3d_lagrange_wave_amatrix((const real(*)[3]) X, px, py,
					       pz, fc->dir[iota], bem, E);
	}
      }
    }
  }

  /* Clean up */
  del_realavector(px);
  del_realavector(py);
  del_realavector(pz);
  freemem(X);
}

static void
assemble_bem3d_inter_wavewave_transfer_row_dclusterbasis(uint rname,
							 pcbem3d bem)
{
  pdclusterbasis rb = bem->par->drbn[rname];
  pcdcluster fc = rb->t;
  uint      ni = bem->aprx->m_inter;

  pcdcluster sc;
  pamatrix  E;
  uint      j, iota, iota1;
  prealavector px, py, pz;
  real(*X)[3];
  real      dir[3];

  X = (real(*)[3]) allocreal(3 * ni * ni * ni);
  px = new_realavector(ni);
  py = new_realavector(ni);
  pz = new_realavector(ni);

  /* Find interpolation points for the father */
  assemble_interpoints3d_realavector(bem, fc->bmin, fc->bmax, px, py, pz);

  for (j = 0; j < rb->sons; j++) {
    sc = rb->son[j]->t;

    /* Find interpolation points for the son */
    assemble_interpoints3d_array(bem, sc->bmin, sc->bmax, X);

    if (rb->son[j]->t->directions > 0) {
      /* Father and son use wave-Lagrange basis */
      assert(rb->son[j]->directions == rb->son[j]->t->directions);

      for (iota = 0; iota < rb->directions; iota++) {
	/* Obtain corresponding direction in the son */
	iota1 = rb->t->dirson[j][iota];
	E = rb->E[j] + iota;

	/* Compute wave-wave transfer matrix */
	if (rb->E[j][iota].cols > 0) {
	  assert(E->rows > 0);
	  assert(iota1 < rb->t->son[j]->directions);
	  assert(E->rows == rb->k[iota]);
	  //         assert(E->cols == rb->son[j]->k[iota1]);                        /* Isn't true for direct orthogonalization */

	  dir[0] = fc->dir[iota][0] - sc->dir[iota1][0];
	  dir[1] = fc->dir[iota][1] - sc->dir[iota1][1];
	  dir[2] = fc->dir[iota][2] - sc->dir[iota1][2];

	  assemble_bem3d_lagrange_wave_amatrix((const real(*)[3]) X, px, py,
					       pz, dir, bem, E);
	}
      }
    }
  }

  /* Clean up */
  del_realavector(px);
  del_realavector(py);
  del_realavector(pz);
  freemem(X);
}

static void
assemble_bem3d_inter_wavewave_transfer_col_dclusterbasis(uint cname,
							 pcbem3d bem)
{
  pdclusterbasis cb = bem->par->dcbn[cname];
  pcdcluster fc = cb->t;
  uint      ni = bem->aprx->m_inter;

  pcdcluster sc;
  pamatrix  E;
  uint      j, iota, iota1;
  prealavector px, py, pz;
  real(*X)[3];
  real      dir[3];

  X = (real(*)[3]) allocreal(3 * ni * ni * ni);
  px = new_realavector(ni);
  py = new_realavector(ni);
  pz = new_realavector(ni);

  /* Find interpolation points for the father */
  assemble_interpoints3d_realavector(bem, fc->bmin, fc->bmax, px, py, pz);

  for (j = 0; j < cb->sons; j++) {
    sc = cb->son[j]->t;

    /* Find interpolation points for the son */
    assemble_interpoints3d_array(bem, sc->bmin, sc->bmax, X);

    if (cb->son[j]->t->directions > 0) {
      /* Father and son use wave-Lagrange basis */
      assert(cb->son[j]->directions == cb->son[j]->t->directions);

      for (iota = 0; iota < cb->directions; iota++) {
	/* Obtain corresponding direction in the son */
	iota1 = cb->t->dirson[j][iota];
	E = cb->E[j] + iota;

	/* Compute wave-wave transfer matrix */
	if (cb->E[j][iota].cols > 0) {
	  assert(E->rows > 0);
	  assert(iota1 < cb->t->son[j]->directions);
	  assert(E->rows == cb->k[iota]);
	  //         assert(E->cols == cb->son[j]->k[iota1]);                  /* Isn't true for direct orthogonalization */

	  dir[0] = fc->dir[iota][0] - sc->dir[iota1][0];
	  dir[1] = fc->dir[iota][1] - sc->dir[iota1][1];
	  dir[2] = fc->dir[iota][2] - sc->dir[iota1][2];

	  assemble_bem3d_lagrange_wave_amatrix((const real(*)[3]) X, px, py,
					       pz, dir, bem, E);
	}
      }
    }
  }

  /* Clean up */
  del_realavector(px);
  del_realavector(py);
  del_realavector(pz);
  freemem(X);
}

static void
assemble_bem3d_inter_nowave_duniform(uint rname, uint cname,
				     uint bname, pcbem3d bem)
{
  pkernelbem3d kernels = bem->kernels;
  pparbem3d par = bem->par;
  pduniform U = par->dh2n[bname]->u;
  pcdcluster rc = U->rb->t;
  pcdcluster cc = U->cb->t;
  uint      k = bem->aprx->k_inter;
  pamatrix  S = &U->S;

  real(*xi_r)[3], (*xi_c)[3];

  (void) rname;
  (void) cname;

  /* Find interpolation points for the row cluster */
  xi_r = (real(*)[3]) allocreal(3 * k);
  assemble_interpoints3d_array(bem, rc->bmin, rc->bmax, xi_r);

  /* Find interpolation points for the column cluster */
  xi_c = (real(*)[3]) allocreal(3 * k);
  assemble_interpoints3d_array(bem, cc->bmin, cc->bmax, xi_c);

  kernels->fundamental(bem, (const real(*)[3]) xi_r, (const real(*)[3]) xi_c,
		       S);

  freemem(xi_r);
  freemem(xi_c);
}

static void
assemble_bem3d_inter_wave_duniform(uint rname, uint cname,
				   uint bname, pcbem3d bem)
{
  pduniform u = bem->par->dh2n[bname]->u;
  pcdcluster rc = u->rb->t;
  pcdcluster cc = u->cb->t;
  uint      rd = u->rd;
  uint      cd = u->cd;
  pamatrix  S = &u->S;
  uint      k = bem->aprx->k_inter;
  real(*xi_r)[3], (*xi_c)[3];

  (void) rname;
  (void) cname;

  assert(S->rows == k);
  assert(S->cols == k);
  assert(cc->dir[cd][0] == rc->dir[rd][0]);
  assert(cc->dir[cd][1] == rc->dir[rd][1]);
  assert(cc->dir[cd][2] == rc->dir[rd][2]);

  /* Find interpolation points for the row cluster */
  xi_r = (real(*)[3]) allocreal(3 * k);
  assemble_interpoints3d_array(bem, rc->bmin, rc->bmax, xi_r);

  /* Find interpolation points for the column cluster */
  xi_c = (real(*)[3]) allocreal(3 * k);
  assemble_interpoints3d_array(bem, cc->bmin, cc->bmax, xi_c);

  /* Evaluate kernel function multiplied by plane wave */
  bem->kernels->fundamental_wave(bem, (const real(*)[3]) xi_r,
				 (const real(*)[3]) xi_c, rc->dir[rd], S);

  /* Clean up */
  freemem(xi_r);
  freemem(xi_c);
}

static void
assemble_bem3d_inter_ortho_nowave_duniform(uint rname, uint cname,
					   uint bname, pcbem3d bem)
{

  pkernelbem3d kernels = bem->kernels;
  pparbem3d par = bem->par;
  pduniform U = par->dh2n[bname]->u;
  pcdcluster rc = U->rb->t;
  pcdcluster cc = U->cb->t;
  uint      rd = U->rd;
  uint      cd = U->cd;
  uint      k = bem->aprx->k_inter;
  pamatrix  S, T;
  amatrix   tmp1, tmp2;

  real(*xi_r)[3], (*xi_c)[3];

  assert(rc == par->ron[rname]->t);
  assert(cc == par->con[cname]->t);

  S = init_amatrix(&tmp1, k, k);
  clear_amatrix(S);

  /* Find interpolation points for the row cluster */
  xi_r = (real(*)[3]) allocreal(3 * k);
  assemble_interpoints3d_array(bem, rc->bmin, rc->bmax, xi_r);

  /* Find interpolation points for the column cluster */
  xi_c = (real(*)[3]) allocreal(3 * k);
  assemble_interpoints3d_array(bem, cc->bmin, cc->bmax, xi_c);

  kernels->fundamental(bem, (const real(*)[3]) xi_r, (const real(*)[3]) xi_c,
		       S);

  /* Resize */
  T = init_amatrix(&tmp2, k, par->con[cname]->C[cd].rows);
  clear_amatrix(T);
  addmul_amatrix(1.0, false, S, true, &par->con[cname]->C[cd], T);
  uninit_amatrix(S);
  resize_amatrix(&U->S, par->ron[rname]->C[rd].rows,
		 par->con[cname]->C[cd].rows);
  clear_amatrix(&U->S);
  addmul_amatrix(1.0, false, &par->ron[rname]->C[rd], false, T, &U->S);
  uninit_amatrix(T);

  /* Clean up */
  freemem(xi_r);
  freemem(xi_c);
}

static void
assemble_bem3d_inter_ortho_wave_duniform(uint rname, uint cname,
					 uint bname, pcbem3d bem)
{

  pparbem3d par = bem->par;
  pduniform U = par->dh2n[bname]->u;
  pcdcluster rc = U->rb->t;
  pcdcluster cc = U->cb->t;
  uint      rd = U->rd;
  uint      cd = U->cd;
  pamatrix  S, T;
  amatrix   tmp1, tmp2;
  uint      k = bem->aprx->k_inter;
  real(*xi_r)[3], (*xi_c)[3];

  assert(rc == par->ron[rname]->t);
  assert(cc == par->con[cname]->t);

  S = init_amatrix(&tmp1, k, k);
  clear_amatrix(S);
  assert(cc->dir[cd][0] == rc->dir[rd][0]);
  assert(cc->dir[cd][1] == rc->dir[rd][1]);
  assert(cc->dir[cd][2] == rc->dir[rd][2]);

  /* Find interpolation points for the row cluster */
  xi_r = (real(*)[3]) allocreal(3 * k);
  assemble_interpoints3d_array(bem, rc->bmin, rc->bmax, xi_r);

  /* Find interpolation points for the column cluster */
  xi_c = (real(*)[3]) allocreal(3 * k);
  assemble_interpoints3d_array(bem, cc->bmin, cc->bmax, xi_c);

  /* Evaluate kernel function multiplied by plane wave */
  bem->kernels->fundamental_wave(bem, (const real(*)[3]) xi_r,
				 (const real(*)[3]) xi_c, rc->dir[rd], S);

  /* Resize */
  T = init_amatrix(&tmp2, k, par->con[cname]->C[cd].rows);
  clear_amatrix(T);
  addmul_amatrix(1.0, false, S, true, &par->con[cname]->C[cd], T);
  uninit_amatrix(S);
  resize_amatrix(&U->S, par->ron[rname]->C[rd].rows,
		 par->con[cname]->C[cd].rows);
  clear_amatrix(&U->S);
  addmul_amatrix(1.0, false, &par->ron[rname]->C[rd], false, T, &U->S);
  uninit_amatrix(T);

  /* Clean up */
  freemem(xi_r);
  freemem(xi_c);
}

/* ------------------------------------------------------------
 lagrange-polynomials
 ------------------------------------------------------------ */

static inline real
eval_lagrange1d(real x, pcreal xi, uint mu, uint m)
{
  uint      l;

  real      lagr;

  lagr = 1.0;

  for (l = 0; l < mu; ++l) {
    lagr *= (x - xi[l]) / (xi[mu] - xi[l]);
  }

  for (l = mu + 1; l < m; ++l) {
    lagr *= (x - xi[l]) / (xi[mu] - xi[l]);
  }

  return lagr;
}

static inline real
eval_dn_lagrange1d(real x, pcreal xi, uint mu, uint m)
{
  uint      l, l2;
  real      tmp;

  real      lagr;

  lagr = 0.0;

  for (l = 0; l < mu; ++l) {
    tmp = 1.0;
    for (l2 = 0; l2 < l; ++l2) {
      tmp *= (x - xi[l2]) / (xi[mu] - xi[l2]);
    }
    for (l2 = l + 1; l2 < mu; ++l2) {
      tmp *= (x - xi[l2]) / (xi[mu] - xi[l2]);
    }
    for (l2 = mu + 1; l2 < m; ++l2) {
      tmp *= (x - xi[l2]) / (xi[mu] - xi[l2]);
    }
    lagr += tmp / (xi[mu] - xi[l]);
  }

  for (l = mu + 1; l < m; ++l) {
    tmp = 1.0;
    for (l2 = 0; l2 < mu; ++l2) {
      tmp *= (x - xi[l2]) / (xi[mu] - xi[l2]);
    }
    for (l2 = mu + 1; l2 < l; ++l2) {
      tmp *= (x - xi[l2]) / (xi[mu] - xi[l2]);
    }
    for (l2 = l + 1; l2 < m; ++l2) {
      tmp *= (x - xi[l2]) / (xi[mu] - xi[l2]);
    }
    lagr += tmp / (xi[mu] - xi[l]);
  }

  return lagr;
}

void
assemble_bem3d_lagrange_c_amatrix(const uint * idx, pcrealavector px,
				  pcrealavector py, pcrealavector pz,
				  pcbem3d bem, pamatrix V)
{

  pcsurface3d gr = bem->gr;
  const     real(*gr_x)[3] = (const real(*)[3]) gr->x;
  const     uint(*gr_t)[3] = (const uint(*)[3]) gr->t;
  const preal gr_g = (const preal) gr->g;
  uint      rows = V->rows;
  longindex ld = V->ld;
  uint      nq = bem->sq->n_single;
  uint      vnq = ROUNDUP(nq, VREAL);
  real     *xx = bem->sq->x_single;
  real     *yy = bem->sq->y_single;
  real     *ww = bem->sq->w_single + 3 * vnq;
  uint      mx = px->dim;
  uint      my = py->dim;
  uint      mz = pz->dim;

  const real *A, *B, *C;
  uint      t, tt, jx, jy, jz, q, index;
  real      gt, sum, lagr, x[3], tx, sx, Ax, Bx, Cx;

  /*
   * integrate Lagrange polynomials with constant basisfunctions
   */

  index = 0;
  for (jx = 0; jx < mx; ++jx) {
    for (jy = 0; jy < my; ++jy) {
      for (jz = 0; jz < mz; ++jz) {
	for (t = 0; t < rows; ++t) {
	  tt = (idx == NULL ? t : idx[t]);
	  gt = gr_g[tt];
	  A = gr_x[gr_t[tt][0]];
	  B = gr_x[gr_t[tt][1]];
	  C = gr_x[gr_t[tt][2]];

	  sum = 0.0;

	  for (q = 0; q < nq; ++q) {
	    tx = xx[q];
	    sx = yy[q];
	    Ax = 1.0 - tx;
	    Bx = tx - sx;
	    Cx = sx;

	    x[0] = A[0] * Ax + B[0] * Bx + C[0] * Cx;
	    x[1] = A[1] * Ax + B[1] * Bx + C[1] * Cx;
	    x[2] = A[2] * Ax + B[2] * Bx + C[2] * Cx;

	    lagr = eval_lagrange1d(x[0], px->v, jx, mx);
	    lagr *= eval_lagrange1d(x[1], py->v, jy, my);
	    lagr *= eval_lagrange1d(x[2], pz->v, jz, mz);

	    sum += ww[q] * lagr;

	  }
	  V->a[t + index * ld] = gt * sum;
	}
	index++;
      }
    }
  }
}

void
assemble_bem3d_lagrange_wave_c_amatrix(const uint * idx, pcrealavector px,
				       pcrealavector py, pcrealavector pz,
				       pcreal dir, pcbem3d bem, pamatrix V)
{
  pcsurface3d gr = bem->gr;
  const     real(*gr_x)[3] = (const real(*)[3]) gr->x;
  const     uint(*gr_t)[3] = (const uint(*)[3]) gr->t;
  const preal gr_g = (const preal) gr->g;
  uint      rows = V->rows;
  longindex ld = V->ld;
  uint      nq = bem->sq->n_single;
  uint      vnq = ROUNDUP(nq, VREAL);
  real     *xx = bem->sq->x_single;
  real     *yy = bem->sq->y_single;
  real     *ww = bem->sq->w_single + 3 * vnq;
  uint      mx = px->dim;
  uint      my = py->dim;
  uint      mz = pz->dim;
  real      wave_k = bem->k;

  const real *A, *B, *C;
  uint      t, tt, jx, jy, jz, q, index;
  real      x[3], sum_r, sum_i, gt, lagr, angle, tx, sx, Ax, Bx, Cx;

  /*
   * integrate Lagrange polynomials with constant basisfunctions
   */

  index = 0;
  for (jx = 0; jx < mx; ++jx) {
    for (jy = 0; jy < my; ++jy) {
      for (jz = 0; jz < mz; ++jz) {
	for (t = 0; t < rows; ++t) {
	  tt = (idx == NULL ? t : idx[t]);
	  gt = gr_g[tt];
	  A = gr_x[gr_t[tt][0]];
	  B = gr_x[gr_t[tt][1]];
	  C = gr_x[gr_t[tt][2]];

	  sum_r = 0.0;
	  sum_i = 0.0;

	  for (q = 0; q < nq; ++q) {
	    tx = xx[q];
	    sx = yy[q];
	    Ax = 1.0 - tx;
	    Bx = tx - sx;
	    Cx = sx;

	    x[0] = A[0] * Ax + B[0] * Bx + C[0] * Cx;
	    x[1] = A[1] * Ax + B[1] * Bx + C[1] * Cx;
	    x[2] = A[2] * Ax + B[2] * Bx + C[2] * Cx;

	    lagr = eval_lagrange1d(x[0], px->v, jx, mx);
	    lagr *= eval_lagrange1d(x[1], py->v, jy, my);
	    lagr *= eval_lagrange1d(x[2], pz->v, jz, mz);

	    lagr = lagr * ww[q];

	    angle = wave_k * REAL_DOT3(x, dir);

	    sum_r += lagr * REAL_COS(angle);
	    sum_i += lagr * REAL_SIN(angle);

	  }
	  V->a[t + index * ld] = (gt * sum_r) + (gt * sum_i) * I;
	}
	index++;
      }
    }
  }
}

void
assemble_bem3d_lagrange_l_amatrix(const uint * idx, pcrealavector px,
				  pcrealavector py, pcrealavector pz,
				  pcbem3d bem, pamatrix V)
{
  pcsurface3d gr = bem->gr;
  const     real(*gr_x)[3] = (const real(*)[3]) gr->x;
  const     uint(*gr_t)[3] = (const uint(*)[3]) gr->t;
  const preal gr_g = (const preal) gr->g;
  plistnode *v2t = bem->v2t;
  uint      rows = V->rows;
  field    *aa = V->a;
  longindex ld = V->ld;
  uint      nq = bem->sq->n_single;
  uint      vnq = ROUNDUP(nq, VREAL);
  real     *xx = bem->sq->x_single;
  real     *yy = bem->sq->y_single;
  real     *ww = bem->sq->w_single;
  real      base = bem->sq->base_single;
  real     *quad;
  uint      mx = px->dim;
  uint      my = py->dim;
  uint      mz = pz->dim;

  ptri_list tl, tl1;
  pvert_list vl;
  const real *A, *B, *C;
  uint      tri_sp[3];
  plistnode v;
  uint      s, ss, i, jx, jy, jz, k, q, cj, index;
  real      gs, sum, lagr, x[3], tx, sx, Ax, Bx, Cx;
  longindex ii, vv;

  quad = allocreal(vnq);

  clear_amatrix(V);

  tl = NULL;

  cj = 0;
  for (i = 0; i < rows; ++i) {
    ii = (idx == NULL ? i : idx[i]);
    for (v = v2t[ii], vv = v->data; v->next != NULL;
	 v = v->next, vv = v->data) {

      tl1 = tl;
      while (tl1 && tl1->t != vv) {
	tl1 = tl1->next;
      }

      if (tl1 == NULL) {
	tl1 = tl = new_tri_list(tl);
	tl->t = vv;
	cj++;
      }

      tl1->vl = new_vert_list(tl1->vl);
      tl1->vl->v = i;
    }
  }

  for (s = 0, tl1 = tl; s < cj; s++, tl1 = tl1->next) {
    ss = tl1->t;
    gs = gr_g[ss];
    A = gr_x[gr_t[ss][0]];
    B = gr_x[gr_t[ss][1]];
    C = gr_x[gr_t[ss][2]];

    for (i = 0; i < 3; ++i) {
      tri_sp[i] = gr_t[ss][i];
    }

    index = 0;

    for (jx = 0; jx < mx; ++jx) {
      for (jy = 0; jy < my; ++jy) {
	for (jz = 0; jz < mz; ++jz) {

	  for (q = 0; q < nq; ++q) {
	    tx = xx[q];
	    sx = yy[q];
	    Ax = 1.0 - tx;
	    Bx = tx - sx;
	    Cx = sx;

	    x[0] = A[0] * Ax + B[0] * Bx + C[0] * Cx;
	    x[1] = A[1] * Ax + B[1] * Bx + C[1] * Cx;
	    x[2] = A[2] * Ax + B[2] * Bx + C[2] * Cx;

	    lagr = eval_lagrange1d(x[0], px->v, jx, mx);
	    lagr *= eval_lagrange1d(x[1], py->v, jy, my);
	    lagr *= eval_lagrange1d(x[2], pz->v, jz, mz);

	    quad[q] = lagr;
	  }

	  ww = bem->sq->w_single;
	  vl = tl1->vl;
	  while (vl) {
	    k = vl->v;
	    if (k < rows) {
	      ii = idx == NULL ? k : idx[k];
	      for (i = 0; i < 3; ++i) {
		if (ii == tri_sp[i]) {
		  sum = base;

		  for (q = 0; q < nq; ++q) {
		    sum += ww[q] * quad[q];
		  }

		  aa[k + index * ld] += sum * gs;
		}
		ww += vnq;
	      }
	      ww = bem->sq->w_single;
	    }
	    vl = vl->next;
	  }

	  index++;
	}
      }
    }
  }

  del_tri_list(tl);
  freemem(quad);
}

void
assemble_bem3d_lagrange_wave_l_amatrix(const uint * idx, pcrealavector px,
				       pcrealavector py, pcrealavector pz,
				       pcreal dir, pcbem3d bem, pamatrix V)
{
  pcsurface3d gr = bem->gr;
  const     real(*gr_x)[3] = (const real(*)[3]) gr->x;
  const     uint(*gr_t)[3] = (const uint(*)[3]) gr->t;
  const preal gr_g = (const preal) gr->g;
  plistnode *v2t = bem->v2t;
  uint      rows = V->rows;
  field    *aa = V->a;
  longindex ld = V->ld;
  uint      nq = bem->sq->n_single;
  uint      vnq = ROUNDUP(nq, VREAL);
  real     *xx = bem->sq->x_single;
  real     *yy = bem->sq->y_single;
  real     *ww = bem->sq->w_single;
  real      base = bem->sq->base_single;
  real     *quad_r, *quad_i;
  uint      mx = px->dim;
  uint      my = py->dim;
  uint      mz = pz->dim;
  real      wave_k = bem->k;

  ptri_list tl, tl1;
  pvert_list vl;
  const real *A, *B, *C;
  uint      tri_sp[3];
  plistnode v;
  uint      s, ss, i, jx, jy, jz, k, q, cj, index;
  real      gs, sum_r, sum_i, lagr, angle, x[3], tx, sx, Ax, Bx, Cx;
  longindex ii, vv;

  quad_r = allocreal(vnq);
  quad_i = allocreal(vnq);

  clear_amatrix(V);

  tl = NULL;

  cj = 0;
  for (i = 0; i < rows; ++i) {
    ii = (idx == NULL ? i : idx[i]);
    for (v = v2t[ii], vv = v->data; v->next != NULL;
	 v = v->next, vv = v->data) {

      tl1 = tl;
      while (tl1 && tl1->t != vv) {
	tl1 = tl1->next;
      }

      if (tl1 == NULL) {
	tl1 = tl = new_tri_list(tl);
	tl->t = vv;
	cj++;
      }

      tl1->vl = new_vert_list(tl1->vl);
      tl1->vl->v = i;
    }
  }

  for (s = 0, tl1 = tl; s < cj; s++, tl1 = tl1->next) {
    ss = tl1->t;
    gs = gr_g[ss];
    A = gr_x[gr_t[ss][0]];
    B = gr_x[gr_t[ss][1]];
    C = gr_x[gr_t[ss][2]];

    for (i = 0; i < 3; ++i) {
      tri_sp[i] = gr_t[ss][i];
    }

    index = 0;

    for (jx = 0; jx < mx; ++jx) {
      for (jy = 0; jy < my; ++jy) {
	for (jz = 0; jz < mz; ++jz) {

	  for (q = 0; q < nq; ++q) {
	    tx = xx[q];
	    sx = yy[q];
	    Ax = 1.0 - tx;
	    Bx = tx - sx;
	    Cx = sx;

	    x[0] = A[0] * Ax + B[0] * Bx + C[0] * Cx;
	    x[1] = A[1] * Ax + B[1] * Bx + C[1] * Cx;
	    x[2] = A[2] * Ax + B[2] * Bx + C[2] * Cx;

	    lagr = eval_lagrange1d(x[0], px->v, jx, mx);
	    lagr *= eval_lagrange1d(x[1], py->v, jy, my);
	    lagr *= eval_lagrange1d(x[2], pz->v, jz, mz);

	    angle = wave_k * REAL_DOT3(x, dir);

	    quad_r[q] = lagr * REAL_COS(angle);
	    quad_i[q] = lagr * REAL_SIN(angle);
	  }

	  ww = bem->sq->w_single;
	  vl = tl1->vl;
	  while (vl) {
	    k = vl->v;
	    if (k < rows) {
	      ii = idx == NULL ? k : idx[k];
	      for (i = 0; i < 3; ++i) {
		if (ii == tri_sp[i]) {
		  sum_r = base;
		  sum_i = 0.0;

		  for (q = 0; q < nq; ++q) {
		    sum_r += ww[q] * quad_r[q];
		    sum_i += ww[q] * quad_i[q];
		  }

		  aa[k + index * ld] += (sum_r * gs) + (sum_i * gs) * I;
		}
		ww += vnq;
	      }
	      ww = bem->sq->w_single;
	    }
	    vl = vl->next;
	  }

	  index++;
	}
      }
    }
  }

  del_tri_list(tl);
  freemem(quad_r);
  freemem(quad_i);
}

void
assemble_bem3d_dn_lagrange_c_amatrix(const uint * idx, pcrealavector px,
				     pcrealavector py, pcrealavector pz,
				     pcbem3d bem, pamatrix V)
{
  pcsurface3d gr = bem->gr;
  const     real(*gr_x)[3] = (const real(*)[3]) gr->x;
  const     uint(*gr_t)[3] = (const uint(*)[3]) gr->t;
  const     real(*gr_n)[3] = (const real(*)[3]) gr->n;
  const preal gr_g = (const preal) gr->g;
  uint      rows = V->rows;
  longindex ld = V->ld;

  uint      nq = bem->sq->n_single;
  uint      vnq = ROUNDUP(nq, VREAL);
  real     *xx = bem->sq->x_single;
  real     *yy = bem->sq->y_single;
  real     *ww = bem->sq->w_single + 3 * vnq;
  uint      mx = px->dim;
  uint      my = py->dim;
  uint      mz = pz->dim;

  const real *A, *B, *C, *nt;
  real      lagr[3];
  uint      t, tt, jx, jy, jz, q, index;
  real      gt, sum, lagrx, lagry, lagrz, x[3], tx, sx, Ax, Bx, Cx;

  /*
   * integrate Lagrange polynomials with constant basisfunctions
   */

  for (t = 0; t < rows; ++t) {
    tt = (idx == NULL ? t : idx[t]);
    gt = gr_g[tt];
    nt = gr_n[tt];
    A = gr_x[gr_t[tt][0]];
    B = gr_x[gr_t[tt][1]];
    C = gr_x[gr_t[tt][2]];

    index = 0;

    for (jx = 0; jx < mx; ++jx) {
      for (jy = 0; jy < my; ++jy) {
	for (jz = 0; jz < mz; ++jz) {

	  sum = 0.0;

	  for (q = 0; q < nq; ++q) {
	    tx = xx[q];
	    sx = yy[q];
	    Ax = 1.0 - tx;
	    Bx = tx - sx;
	    Cx = sx;

	    x[0] = A[0] * Ax + B[0] * Bx + C[0] * Cx;
	    x[1] = A[1] * Ax + B[1] * Bx + C[1] * Cx;
	    x[2] = A[2] * Ax + B[2] * Bx + C[2] * Cx;

	    lagrx = eval_lagrange1d(x[0], px->v, jx, mx);
	    lagry = eval_lagrange1d(x[1], py->v, jy, my);
	    lagrz = eval_lagrange1d(x[2], pz->v, jz, mz);

	    lagr[0] = eval_dn_lagrange1d(x[0], px->v, jx, mx) * lagry * lagrz;
	    lagr[1] = lagrx * eval_dn_lagrange1d(x[1], py->v, jy, my) * lagrz;
	    lagr[2] = lagrx * lagry * eval_dn_lagrange1d(x[2], pz->v, jz, mz);

	    sum += ww[q] * REAL_DOT3(nt, lagr);
	  }

	  V->a[t + index * ld] = gt * sum;
	  index++;
	}
      }
    }
  }
}

void
assemble_bem3d_dn_lagrange_wave_c_amatrix(const uint * idx,
					  pcrealavector px, pcrealavector py,
					  pcrealavector pz, pcreal dir,
					  pcbem3d bem, pamatrix V)
{
  pcsurface3d gr = bem->gr;
  const     real(*gr_x)[3] = (const real(*)[3]) gr->x;
  const     uint(*gr_t)[3] = (const uint(*)[3]) gr->t;
  const     real(*gr_n)[3] = (const real(*)[3]) gr->n;
  const preal gr_g = (const preal) gr->g;
  uint      rows = V->rows;
  longindex ld = V->ld;

  uint      nq = bem->sq->n_single;
  uint      vnq = ROUNDUP(nq, VREAL);
  real     *xx = bem->sq->x_single;
  real     *yy = bem->sq->y_single;
  real     *ww = bem->sq->w_single + 3 * vnq;
  uint      mx = px->dim;
  uint      my = py->dim;
  uint      mz = pz->dim;
  real      wave_k = bem->k;

  const real *A, *B, *C, *nt;
  real      lagr[3];
  uint      t, tt, jx, jy, jz, q, index;
  real      gt, lagrx, lagry, lagrz, angle, x[3], tx, sx, Ax, Bx, Cx, sum_r,
    sum_i;
  real      cs, sn, re[3], im[3];

  /*
   * integrate Lagrange polynomials with constant basisfunctions
   */

  index = 0;

  for (jx = 0; jx < mx; ++jx) {
    for (jy = 0; jy < my; ++jy) {
      for (jz = 0; jz < mz; ++jz) {
	for (t = 0; t < rows; ++t) {
	  tt = (idx == NULL ? t : idx[t]);
	  gt = gr_g[tt];
	  nt = gr_n[tt];
	  A = gr_x[gr_t[tt][0]];
	  B = gr_x[gr_t[tt][1]];
	  C = gr_x[gr_t[tt][2]];

	  sum_r = 0.0;
	  sum_i = 0.0;

	  for (q = 0; q < nq; ++q) {
	    tx = xx[q];
	    sx = yy[q];
	    Ax = 1.0 - tx;
	    Bx = tx - sx;
	    Cx = sx;

	    x[0] = A[0] * Ax + B[0] * Bx + C[0] * Cx;
	    x[1] = A[1] * Ax + B[1] * Bx + C[1] * Cx;
	    x[2] = A[2] * Ax + B[2] * Bx + C[2] * Cx;

	    lagrx = eval_lagrange1d(x[0], px->v, jx, mx);
	    lagry = eval_lagrange1d(x[1], py->v, jy, my);
	    lagrz = eval_lagrange1d(x[2], pz->v, jz, mz);

	    lagr[0] = eval_dn_lagrange1d(x[0], px->v, jx, mx) * lagry * lagrz;
	    lagr[1] = lagrx * eval_dn_lagrange1d(x[1], py->v, jy, my) * lagrz;
	    lagr[2] = lagrx * lagry * eval_dn_lagrange1d(x[2], pz->v, jz, mz);

	    angle = wave_k * REAL_DOT3(x, dir);
	    cs = REAL_COS(angle);
	    sn = REAL_SIN(angle);

	    angle = wave_k * lagrx * lagry * lagrz;
	    im[0] = angle * dir[0];
	    im[1] = angle * dir[1];
	    im[2] = angle * dir[2];

	    re[0] = cs * lagr[0] - sn * im[0];
	    re[1] = cs * lagr[1] - sn * im[1];
	    re[2] = cs * lagr[2] - sn * im[2];

	    im[0] = cs * im[0] + sn * lagr[0];
	    im[1] = cs * im[1] + sn * lagr[1];
	    im[2] = cs * im[2] + sn * lagr[2];

	    sum_r += ww[q] * REAL_DOT3(nt, re);
	    sum_i += ww[q] * REAL_DOT3(nt, im);
	  }

	  V->a[t + index * ld] = (gt * sum_r) + (gt * sum_i) * I;
	}
	index++;
      }
    }
  }
}

void
assemble_bem3d_dn_lagrange_l_amatrix(const uint * idx, pcrealavector px,
				     pcrealavector py, pcrealavector pz,
				     pcbem3d bem, pamatrix V)
{
  pcsurface3d gr = bem->gr;
  const     real(*gr_x)[3] = (const real(*)[3]) gr->x;
  const     uint(*gr_t)[3] = (const uint(*)[3]) gr->t;
  const     real(*gr_n)[3] = (const real(*)[3]) gr->n;
  const preal gr_g = (const preal) gr->g;
  plistnode *v2t = bem->v2t;
  uint      rows = V->rows;
  field    *aa = V->a;
  longindex ld = V->ld;
  uint      nq = bem->sq->n_single;
  uint      vnq = ROUNDUP(nq, VREAL);
  real     *xx = bem->sq->x_single;
  real     *yy = bem->sq->y_single;
  real     *ww = bem->sq->w_single;
  real      base = bem->sq->base_single;
  real     *quad;
  uint      mx = px->dim;
  uint      my = py->dim;
  uint      mz = pz->dim;

  ptri_list tl, tl1;
  pvert_list vl;
  const real *A, *B, *C, *ns;
  real      lagr[3];
  uint      tri_sp[3];
  plistnode v;
  uint      s, ss, i, jx, jy, jz, k, q, cj, index;
  real      gs, sum, lagrx, lagry, lagrz, x[3], tx, sx, Ax, Bx, Cx;
  longindex ii, vv;

  quad = allocreal(vnq);

  clear_amatrix(V);

  tl = NULL;

  cj = 0;
  for (i = 0; i < rows; ++i) {
    ii = (idx == NULL ? i : idx[i]);
    for (v = v2t[ii], vv = v->data; v->next != NULL;
	 v = v->next, vv = v->data) {

      tl1 = tl;
      while (tl1 && tl1->t != vv) {
	tl1 = tl1->next;
      }

      if (tl1 == NULL) {
	tl1 = tl = new_tri_list(tl);
	tl->t = vv;
	cj++;
      }

      tl1->vl = new_vert_list(tl1->vl);
      tl1->vl->v = i;
    }
  }

  for (s = 0, tl1 = tl; s < cj; s++, tl1 = tl1->next) {
    ss = tl1->t;
    gs = gr_g[ss];
    ns = gr_n[ss];
    A = gr_x[gr_t[ss][0]];
    B = gr_x[gr_t[ss][1]];
    C = gr_x[gr_t[ss][2]];

    for (i = 0; i < 3; ++i) {
      tri_sp[i] = gr_t[ss][i];
    }

    index = 0;

    for (jx = 0; jx < mx; ++jx) {
      for (jy = 0; jy < my; ++jy) {
	for (jz = 0; jz < mz; ++jz) {

	  for (q = 0; q < nq; ++q) {
	    tx = xx[q];
	    sx = yy[q];
	    Ax = 1.0 - tx;
	    Bx = tx - sx;
	    Cx = sx;

	    x[0] = A[0] * Ax + B[0] * Bx + C[0] * Cx;
	    x[1] = A[1] * Ax + B[1] * Bx + C[1] * Cx;
	    x[2] = A[2] * Ax + B[2] * Bx + C[2] * Cx;

	    lagrx = eval_lagrange1d(x[0], px->v, jx, mx);
	    lagry = eval_lagrange1d(x[1], py->v, jy, my);
	    lagrz = eval_lagrange1d(x[2], pz->v, jz, mz);

	    lagr[0] = eval_dn_lagrange1d(x[0], px->v, jx, mx) * lagry * lagrz;
	    lagr[1] = lagrx * eval_dn_lagrange1d(x[1], py->v, jy, my) * lagrz;
	    lagr[2] = lagrx * lagry * eval_dn_lagrange1d(x[2], pz->v, jz, mz);

	    quad[q] = REAL_DOT3(ns, lagr);
	  }

	  ww = bem->sq->w_single;
	  vl = tl1->vl;
	  while (vl) {
	    k = vl->v;
	    if (k < rows) {
	      ii = idx == NULL ? k : idx[k];
	      for (i = 0; i < 3; ++i) {
		if (ii == tri_sp[i]) {
		  sum = base;

		  for (q = 0; q < nq; ++q) {
		    sum += ww[q] * quad[q];
		  }

		  aa[k + index * ld] += sum * gs;
		}
		ww += vnq;
	      }
	      ww = bem->sq->w_single;
	    }
	    vl = vl->next;
	  }

	  index++;
	}
      }
    }
  }

  del_tri_list(tl);
  freemem(quad);
}

void
assemble_bem3d_dn_lagrange_wave_l_amatrix(const uint * idx,
					  pcrealavector px, pcrealavector py,
					  pcrealavector pz, pcreal dir,
					  pcbem3d bem, pamatrix V)
{
  pcsurface3d gr = bem->gr;
  const     real(*gr_x)[3] = (const real(*)[3]) gr->x;
  const     uint(*gr_t)[3] = (const uint(*)[3]) gr->t;
  const     real(*gr_n)[3] = (const real(*)[3]) gr->n;
  const preal gr_g = (const preal) gr->g;
  plistnode *v2t = bem->v2t;
  uint      rows = V->rows;
  field    *aa = V->a;
  longindex ld = V->ld;
  uint      nq = bem->sq->n_single;
  uint      vnq = ROUNDUP(nq, VREAL);
  real     *xx = bem->sq->x_single;
  real     *yy = bem->sq->y_single;
  real     *ww = bem->sq->w_single;
  real      base = bem->sq->base_single;
  real     *quad_r, *quad_i;
  uint      mx = px->dim;
  uint      my = py->dim;
  uint      mz = pz->dim;
  real      wave_k = bem->k;

  ptri_list tl, tl1;
  pvert_list vl;
  const real *A, *B, *C, *ns;
  real      lagr[3];
  uint      tri_sp[3];
  plistnode v;
  uint      s, ss, i, jx, jy, jz, k, q, cj, index;
  real      gs, sum_r, sum_i, lagrx, lagry, lagrz, angle, cs, sn, x[3], re[3],
    im[3], tx, sx, Ax, Bx, Cx;
  longindex ii, vv;

  quad_r = allocreal(vnq);
  quad_i = allocreal(vnq);

  clear_amatrix(V);

  tl = NULL;

  cj = 0;
  for (i = 0; i < rows; ++i) {
    ii = (idx == NULL ? i : idx[i]);
    for (v = v2t[ii], vv = v->data; v->next != NULL;
	 v = v->next, vv = v->data) {

      tl1 = tl;
      while (tl1 && tl1->t != vv) {
	tl1 = tl1->next;
      }

      if (tl1 == NULL) {
	tl1 = tl = new_tri_list(tl);
	tl->t = vv;
	cj++;
      }

      tl1->vl = new_vert_list(tl1->vl);
      tl1->vl->v = i;
    }
  }

  for (s = 0, tl1 = tl; s < cj; s++, tl1 = tl1->next) {
    ss = tl1->t;
    gs = gr_g[ss];
    ns = gr_n[ss];
    A = gr_x[gr_t[ss][0]];
    B = gr_x[gr_t[ss][1]];
    C = gr_x[gr_t[ss][2]];

    for (i = 0; i < 3; ++i) {
      tri_sp[i] = gr_t[ss][i];
    }

    index = 0;

    for (jx = 0; jx < mx; ++jx) {
      for (jy = 0; jy < my; ++jy) {
	for (jz = 0; jz < mz; ++jz) {

	  for (q = 0; q < nq; ++q) {
	    tx = xx[q];
	    sx = yy[q];
	    Ax = 1.0 - tx;
	    Bx = tx - sx;
	    Cx = sx;

	    x[0] = A[0] * Ax + B[0] * Bx + C[0] * Cx;
	    x[1] = A[1] * Ax + B[1] * Bx + C[1] * Cx;
	    x[2] = A[2] * Ax + B[2] * Bx + C[2] * Cx;

	    lagrx = eval_lagrange1d(x[0], px->v, jx, mx);
	    lagry = eval_lagrange1d(x[1], py->v, jy, my);
	    lagrz = eval_lagrange1d(x[2], pz->v, jz, mz);

	    lagr[0] = eval_dn_lagrange1d(x[0], px->v, jx, mx) * lagry * lagrz;
	    lagr[1] = lagrx * eval_dn_lagrange1d(x[1], py->v, jy, my) * lagrz;
	    lagr[2] = lagrx * lagry * eval_dn_lagrange1d(x[2], pz->v, jz, mz);

	    angle = wave_k * REAL_DOT3(x, dir);
	    cs = REAL_COS(angle);
	    sn = REAL_SIN(angle);

	    angle = wave_k * lagrx * lagry * lagrz;
	    im[0] = angle * dir[0];
	    im[1] = angle * dir[1];
	    im[2] = angle * dir[2];

	    re[0] = cs * lagr[0] - sn * im[0];
	    re[1] = cs * lagr[1] - sn * im[1];
	    re[2] = cs * lagr[2] - sn * im[2];

	    im[0] = cs * im[0] + sn * lagr[0];
	    im[1] = cs * im[1] + sn * lagr[1];
	    im[2] = cs * im[2] + sn * lagr[2];

	    quad_r[q] = REAL_DOT3(ns, re);
	    quad_i[q] = REAL_DOT3(ns, im);
	  }

	  ww = bem->sq->w_single;
	  vl = tl1->vl;
	  while (vl) {
	    k = vl->v;
	    if (k < rows) {
	      ii = idx == NULL ? k : idx[k];
	      for (i = 0; i < 3; ++i) {
		if (ii == tri_sp[i]) {
		  sum_r = base;
		  sum_i = 0.0;

		  for (q = 0; q < nq; ++q) {
		    sum_r += ww[q] * quad_r[q];
		    sum_i += ww[q] * quad_i[q];
		  }

		  aa[k + index * ld] += (sum_r * gs) + (sum_i * gs) * I;
		}
		ww += vnq;
	      }
	      ww = bem->sq->w_single;
	    }
	    vl = vl->next;
	  }

	  index++;
	}
      }
    }
  }

  del_tri_list(tl);
  freemem(quad_r);
  freemem(quad_i);
}

void
assemble_bem3d_lagrange_amatrix(const real(*X)[3], pcrealavector px,
				pcrealavector py, pcrealavector pz,
				pcbem3d bem, pamatrix V)
{
  const uint rows = V->rows;
  const uint cols = V->cols;
  const longindex ld = V->ld;
  const uint mx = px->dim;
  const uint my = py->dim;
  const uint mz = pz->dim;

  uint      jx, jy, jz, i, index;
  real      lagr;

  (void) bem;

  /*
   * Eval Lagrange polynomials at points X
   */

  assert(mx * my * mz == cols);

  index = 0;
  for (jx = 0; jx < mx; ++jx) {
    for (jy = 0; jy < my; ++jy) {
      for (jz = 0; jz < mz; ++jz) {
	for (i = 0; i < rows; ++i) {
	  lagr = eval_lagrange1d(X[i][0], px->v, jx, mx);
	  lagr *= eval_lagrange1d(X[i][1], py->v, jy, my);
	  lagr *= eval_lagrange1d(X[i][2], pz->v, jz, mz);

	  V->a[i + index * ld] = lagr;
	}
	index++;
      }
    }
  }
}

void
assemble_bem3d_lagrange_wave_amatrix(const real(*X)[3], pcrealavector px,
				     pcrealavector py, pcrealavector pz,
				     pcreal dir, pcbem3d bem, pamatrix V)
{
  const uint rows = V->rows;
  const uint cols = V->cols;
  const longindex ld = V->ld;
  const uint mx = px->dim;
  const uint my = py->dim;
  const uint mz = pz->dim;
  real      wave_k = bem->k;

  uint      jx, jy, jz, i, index;
  real      lagr, angle;

  /*
   * Eval Lagrange polynomials at points X
   */

  assert(mx * my * mz == cols);

  index = 0;
  for (jx = 0; jx < mx; ++jx) {
    for (jy = 0; jy < my; ++jy) {
      for (jz = 0; jz < mz; ++jz) {
	for (i = 0; i < rows; ++i) {
	  lagr = eval_lagrange1d(X[i][0], px->v, jx, mx);
	  lagr *= eval_lagrange1d(X[i][1], py->v, jy, my);
	  lagr *= eval_lagrange1d(X[i][2], pz->v, jz, mz);

	  angle = wave_k * REAL_DOT3(X[i], dir);

	  V->a[i + index * ld] = (lagr * REAL_COS(angle))
	    + (lagr * REAL_SIN(angle)) * I;
	}
	index++;
      }
    }
  }
}

void
integrate_bem3d_c_avector(pbem3d bem, boundary_func3d f, pavector w,
			  void *data)
{
  pcsurface3d gr = bem->gr;
  const     real(*gr_x)[3] = (const real(*)[3]) gr->x;
  const     uint(*gr_t)[3] = (const uint(*)[3]) gr->t;
  const     real(*gr_n)[3] = (const real(*)[3]) gr->n;
  const real *gr_g = gr->g;
  uint      nq = bem->sq->n_single;
  uint      vnq = ROUNDUP(nq, VREAL);
  real     *xx = bem->sq->x_single;
  real     *yy = bem->sq->y_single;
  real     *ww = bem->sq->w_single + 3 * vnq;
  uint      n = w->dim;

  const real *A, *B, *C, *N;
  real      x[3];
  real      tx, sx, Ax, Bx, Cx;
  uint      t, q;
  field     sum;

  /*
   *  integrate function with constant basisfunctions
   */

  for (t = 0; t < n; ++t) {
    A = gr_x[gr_t[t][0]];
    B = gr_x[gr_t[t][1]];
    C = gr_x[gr_t[t][2]];
    N = gr_n[t];

    sum = 0.0;

    for (q = 0; q < nq; ++q) {
      tx = xx[q];
      sx = yy[q];
      Ax = 1.0 - tx;
      Bx = tx - sx;
      Cx = sx;

      x[0] = A[0] * Ax + B[0] * Bx + C[0] * Cx;
      x[1] = A[1] * Ax + B[1] * Bx + C[1] * Cx;
      x[2] = A[2] * Ax + B[2] * Bx + C[2] * Cx;

      sum += ww[q] * f(x, N, data);
    }

    w->v[t] = sum * gr_g[t];
  }
}

void
integrate_bem3d_l_avector(pbem3d bem, boundary_func3d f, pavector w,
			  void *data)
{
  pcsurface3d gr = bem->gr;
  const     real(*gr_x)[3] = (const real(*)[3]) gr->x;
  const     uint(*gr_t)[3] = (const uint(*)[3]) gr->t;
  const     real(*gr_n)[3] = (const real(*)[3]) gr->n;
  const real *gr_g = (const real(*)) gr->g;
  const uint triangles = gr->triangles;
  const uint vertices = gr->vertices;
  uint      nq = bem->sq->n_single;
  uint      vnq = ROUNDUP(nq, VREAL);
  real     *xx = bem->sq->x_single;
  real     *yy = bem->sq->y_single;
  real     *ww = bem->sq->w_single;
  real      base = bem->sq->base_single;

  const real *A, *B, *C, *N;
  field    *quad, sum;
  real      x[3];
  const uint *tri_t;
  real      tx, sx, Ax, Bx, Cx, gt_fac;
  uint      t, q, i;
  longindex ii;

  assert(vertices == w->dim);

  quad = allocfield(vnq);
  clear_avector(w);

  for (t = 0; t < triangles; t++) {
    tri_t = gr_t[t];
    gt_fac = gr_g[t];
    A = gr_x[tri_t[0]];
    B = gr_x[tri_t[1]];
    C = gr_x[tri_t[2]];
    N = gr_n[t];

    for (q = 0; q < nq; ++q) {
      tx = xx[q];
      sx = yy[q];
      Ax = 1.0 - tx;
      Bx = tx - sx;
      Cx = sx;

      x[0] = A[0] * Ax + B[0] * Bx + C[0] * Cx;
      x[1] = A[1] * Ax + B[1] * Bx + C[1] * Cx;
      x[2] = A[2] * Ax + B[2] * Bx + C[2] * Cx;

      quad[q] = f(x, N, data);
    }

    ww = bem->sq->w_single;

    for (i = 0; i < 3; ++i) {
      ii = tri_t[i];
      sum = base;

      for (q = 0; q < nq; ++q) {
	sum += ww[q] * quad[q];
      }

      assert(ii < vertices);
      w->v[ii] += sum * gt_fac;

      ww += vnq;
    }
  }

  freemem(quad);
}

real
normL2_bem3d(pbem3d bem, boundary_func3d rhs, void *data)
{
  pcsurface3d gr = bem->gr;
  const     real(*gr_x)[3] = (const real(*)[3]) gr->x;
  const     uint(*gr_t)[3] = (const uint(*)[3]) gr->t;
  const     real(*gr_n)[3] = (const real(*)[3]) gr->n;
  const real *gr_g = gr->g;
  uint      n = gr->triangles;
  real     *xq = bem->sq->x_single;
  real     *yq = bem->sq->y_single;
  uint      nq = bem->sq->n_single;
  uint      vnq = ROUNDUP(nq, VREAL);
  real     *wq = bem->sq->w_single + 3 * vnq;

  const real *A, *B, *C, *N;
  real      X[3];
  real      sum, tx, sx, Ax, Bx, Cx;
  uint      t, q;

  real      norm;

  norm = 0.0;
  for (t = 0; t < n; ++t) {
    A = gr_x[gr_t[t][0]];
    B = gr_x[gr_t[t][1]];
    C = gr_x[gr_t[t][2]];
    N = gr_n[t];

    sum = 0.0;
    for (q = 0; q < nq; ++q) {
      tx = xq[q];
      sx = yq[q];
      Ax = 1.0 - tx;
      Bx = tx - sx;
      Cx = sx;

      X[0] = A[0] * Ax + B[0] * Bx + C[0] * Cx;
      X[1] = A[1] * Ax + B[1] * Bx + C[1] * Cx;
      X[2] = A[2] * Ax + B[2] * Bx + C[2] * Cx;

      sum += wq[q] * ABSSQR(rhs(X, N, data));
    }
    norm += gr_g[t] * sum;
  }
  norm = REAL_SQRT(norm);

  return norm;
}

real
normL2_c_bem3d(pbem3d bem, pavector x)
{
  uint      n = x->dim;
  const real *gr_g = bem->gr->g;
  field    *xv = x->v;

  uint      t;

  real      norm;

  assert(n == bem->gr->triangles);

  norm = 0.0;
  for (t = 0; t < n; ++t) {
    norm += gr_g[t] * 0.5 * ABSSQR(xv[t]);
  }
  norm = REAL_SQRT(norm);

  return norm;
}

real
normL2diff_c_bem3d(pbem3d bem, pavector x, boundary_func3d rhs, void *data)
{
  uint      n = x->dim;
  pcsurface3d gr = bem->gr;
  const     real(*gr_x)[3] = (const real(*)[3]) gr->x;
  const     uint(*gr_t)[3] = (const uint(*)[3]) gr->t;
  const     real(*gr_n)[3] = (const real(*)[3]) gr->n;
  const real *gr_g = gr->g;
  real     *xq = bem->sq->x_single;
  real     *yq = bem->sq->y_single;
  uint      nq = bem->sq->n_single;
  uint      vnq = ROUNDUP(nq, VREAL);
  real     *wq = bem->sq->w_single + 3 * vnq;
  field    *xv = x->v;

  const real *A, *B, *C, *N;
  real      X[3];
  real      sum, tx, sx, Ax, Bx, Cx;
  uint      t, q;

  real      norm;

  norm = 0.0;
  for (t = 0; t < n; ++t) {
    A = gr_x[gr_t[t][0]];
    B = gr_x[gr_t[t][1]];
    C = gr_x[gr_t[t][2]];
    N = gr_n[t];

    sum = 0.0;
    for (q = 0; q < nq; ++q) {
      tx = xq[q];
      sx = yq[q];
      Ax = 1.0 - tx;
      Bx = tx - sx;
      Cx = sx;

      X[0] = A[0] * Ax + B[0] * Bx + C[0] * Cx;
      X[1] = A[1] * Ax + B[1] * Bx + C[1] * Cx;
      X[2] = A[2] * Ax + B[2] * Bx + C[2] * Cx;

      sum += wq[q] * ABSSQR(rhs(X, N, data) - xv[t]);
    }
    norm += gr_g[t] * sum;
  }
  norm = REAL_SQRT(norm);

  return norm;
}

void
projectL2_bem3d_c_avector(pbem3d bem, boundary_func3d f, pavector w,
			  void *data)
{
  pcsurface3d gr = bem->gr;
  const     real(*gr_x)[3] = (const real(*)[3]) gr->x;
  const     uint(*gr_t)[3] = (const uint(*)[3]) gr->t;
  const     real(*gr_n)[3] = (const real(*)[3]) gr->n;
  uint      nq = bem->sq->n_single;
  uint      vnq = ROUNDUP(nq, VREAL);
  real     *xx = bem->sq->x_single;
  real     *yy = bem->sq->y_single;
  real     *ww = bem->sq->w_single + 3 * vnq;
  uint      n = w->dim;

  const real *A, *B, *C, *N;
  real      x[3];
  real      tx, sx, Ax, Bx, Cx;
  uint      t, q;
  field     sum;

  /*
   *  integrate function with constant basisfunctions
   */

  for (t = 0; t < n; ++t) {
    A = gr_x[gr_t[t][0]];
    B = gr_x[gr_t[t][1]];
    C = gr_x[gr_t[t][2]];
    N = gr_n[t];

    sum = 0.0;

    for (q = 0; q < nq; ++q) {
      tx = xx[q];
      sx = yy[q];
      Ax = 1.0 - tx;
      Bx = tx - sx;
      Cx = sx;

      x[0] = A[0] * Ax + B[0] * Bx + C[0] * Cx;
      x[1] = A[1] * Ax + B[1] * Bx + C[1] * Cx;
      x[2] = A[2] * Ax + B[2] * Bx + C[2] * Cx;

      sum += ww[q] * f(x, N, data);
    }

    w->v[t] = 2.0 * sum;
  }
}

static void
addeval_mass_linear_bem3d(field alpha, void *A, pcavector x, pavector y)
{
  pcbem3d   bem = (pcbem3d) A;
  pcsurface3d gr = bem->gr;
  plistnode *v2t = bem->v2t;
  preal     gr_g = gr->g;
  uint      n = x->dim;

  uint     *tri_k;
  plistnode v;
  uint      i, j;
  real      c1, c2;
  longindex vv;

  c1 = 1.0 / 24.0;
  c2 = 1.0 / 12.0;

  assert(y->dim == n);

  for (i = 0; i < n; ++i) {
    for (v = v2t[i], vv = v->data; v->next != NULL; v = v->next, vv = v->data) {
      tri_k = gr->t[vv];
      for (j = 0; j < 3; ++j) {
	if (i != tri_k[j]) {
	  y->v[i] += alpha * x->v[tri_k[j]] * gr_g[vv] * c1;
	}
      }
      y->v[i] += alpha * x->v[i] * gr_g[vv] * c2;
    }
  }
}

real
normL2_l_bem3d(pbem3d bem, pavector x)
{
  pcsurface3d gr = bem->gr;
  uint      triangles = gr->triangles;
  const     uint(*gr_t)[3] = (const uint(*)[3]) gr->t;
  const real *gr_g = gr->g;
  real     *xq = bem->sq->x_single;
  real     *yq = bem->sq->y_single;
  uint      nq = bem->sq->n_single;
  uint      vnq = ROUNDUP(nq, VREAL);
  real     *wq = bem->sq->w_single + 3 * vnq;
  field    *xv = x->v;

  real      sum, bf, tx, sx, Ax, Bx, Cx;
  uint      t, q;

  real      norm;

  assert(x->dim == gr->vertices);

  norm = 0.0;
  for (t = 0; t < triangles; ++t) {

    sum = 0.0;
    for (q = 0; q < nq; ++q) {
      tx = xq[q];
      sx = yq[q];
      Ax = 1.0 - tx;
      Bx = tx - sx;
      Cx = sx;

      bf = xv[gr_t[t][0]] * Ax + xv[gr_t[t][1]] * Bx + xv[gr_t[t][2]] * Cx;

      sum += wq[q] * REAL_SQR(bf);
    }
    norm += gr_g[t] * sum;
  }
  norm = REAL_SQRT(norm);

  return norm;
}

real
normL2diff_l_bem3d(pbem3d bem, pavector x, boundary_func3d rhs, void *data)
{
  pcsurface3d gr = bem->gr;
  uint      triangles = gr->triangles;
  const     real(*gr_x)[3] = (const real(*)[3]) gr->x;
  const     uint(*gr_t)[3] = (const uint(*)[3]) gr->t;
  const     real(*gr_n)[3] = (const real(*)[3]) gr->n;
  const real *gr_g = gr->g;
  real     *xq = bem->sq->x_single;
  real     *yq = bem->sq->y_single;
  uint      nq = bem->sq->n_single;
  uint      vnq = ROUNDUP(nq, VREAL);
  real     *wq = bem->sq->w_single + 3 * vnq;
  field    *xv = x->v;

  const real *A, *B, *C, *N;
  real      X[3];
  real      sum, bf, tx, sx, Ax, Bx, Cx;
  uint      t, q;

  real      norm;

  assert(x->dim == gr->vertices);

  norm = 0.0;
  for (t = 0; t < triangles; ++t) {
    A = gr_x[gr_t[t][0]];
    B = gr_x[gr_t[t][1]];
    C = gr_x[gr_t[t][2]];
    N = gr_n[t];

    sum = 0.0;
    for (q = 0; q < nq; ++q) {
      tx = xq[q];
      sx = yq[q];
      Ax = 1.0 - tx;
      Bx = tx - sx;
      Cx = sx;

      X[0] = A[0] * Ax + B[0] * Bx + C[0] * Cx;
      X[1] = A[1] * Ax + B[1] * Bx + C[1] * Cx;
      X[2] = A[2] * Ax + B[2] * Bx + C[2] * Cx;

      bf = xv[gr_t[t][0]] * Ax + xv[gr_t[t][1]] * Bx + xv[gr_t[t][2]] * Cx;

      sum += wq[q] * REAL_SQR(rhs(X, N, data) - bf);
    }
    norm += gr_g[t] * sum;
  }
  norm = REAL_SQRT(norm);

  return norm;
}

void
projectL2_bem3d_l_avector(pbem3d bem, boundary_func3d f, pavector w,
			  void *data)
{
  pcsurface3d gr = bem->gr;
  const     real(*gr_x)[3] = (const real(*)[3]) gr->x;
  const     uint(*gr_t)[3] = (const uint(*)[3]) gr->t;
  const     real(*gr_n)[3] = (const real(*)[3]) gr->n;
  const real *gr_g = (const real(*)) gr->g;
  const uint triangles = gr->triangles;
  const uint vertices = gr->vertices;
  uint      nq = bem->sq->n_single;
  uint      vnq = ROUNDUP(nq, VREAL);
  real     *xx = bem->sq->x_single;
  real     *yy = bem->sq->y_single;
  real     *ww = bem->sq->w_single;
  real      base = bem->sq->base_single;

  pavector  v;
  const real *A, *B, *C, *N;
  field    *quad, sum;
  real      x[3];
  const uint *tri_t;
  real      tx, sx, Ax, Bx, Cx, gt_fac;
  uint      t, q, i;
  longindex ii;

  assert(vertices == w->dim);

  quad = allocfield(vnq);
  v = new_avector(vertices);
  clear_avector(v);

  for (t = 0; t < triangles; t++) {
    tri_t = gr_t[t];
    gt_fac = gr_g[t];
    A = gr_x[tri_t[0]];
    B = gr_x[tri_t[1]];
    C = gr_x[tri_t[2]];
    N = gr_n[t];

    for (q = 0; q < nq; ++q) {
      tx = xx[q];
      sx = yy[q];
      Ax = 1.0 - tx;
      Bx = tx - sx;
      Cx = sx;

      x[0] = A[0] * Ax + B[0] * Bx + C[0] * Cx;
      x[1] = A[1] * Ax + B[1] * Bx + C[1] * Cx;
      x[2] = A[2] * Ax + B[2] * Bx + C[2] * Cx;

      quad[q] = f(x, N, data);
    }

    ww = bem->sq->w_single;

    for (i = 0; i < 3; ++i) {
      ii = tri_t[i];
      sum = base;

      for (q = 0; q < nq; ++q) {
	sum += ww[q] * quad[q];
      }

      assert(ii < vertices);
      v->v[ii] += sum * gt_fac;

      ww += vnq;
    }
  }

  random_real_avector(w);

  solve_cg_avector((void *) bem, (addeval_t) addeval_mass_linear_bem3d, v, w,
		   1.0e-15, 500);

  del_avector(v);
  freemem(quad);
}

prkmatrix
build_bem3d_rkmatrix(pccluster row, pccluster col, void *data)
{
  pcbem3d   bem = (pcbem3d) data;

  prkmatrix R;

  R = new_rkmatrix(row->size, col->size, 0);
  bem->farfield_rk(row, 0, col, 0, bem, R);

  return R;
}

pamatrix
build_bem3d_amatrix(pccluster row, pccluster col, void *data)
{
  pcbem3d   bem = (pcbem3d) data;

  pamatrix  N;

  N = new_amatrix(row->size, col->size);
  bem->nearfield(row->idx, col->idx, bem, false, N);

  return N;
}

/* ------------------------------------------------------------
 Initializerfunctions for h-matrix approximations
 ------------------------------------------------------------ */

void
setup_hmatrix_recomp_bem3d(pbem3d bem, bool recomp, real accur_recomp,
			   bool coarsen, real accur_coarsen)
{
  paprxbem3d aprx = bem->aprx;

  assert(accur_recomp >= 0.0);
  assert(accur_coarsen >= 0.0);

  uninit_recompression_bem3d(aprx);

  aprx->recomp = recomp;
  aprx->accur_recomp = accur_recomp;
  aprx->coarsen = coarsen;
  aprx->accur_coarsen = accur_coarsen;
}

/* ------------------------------------------------------------
 Interpolation
 ------------------------------------------------------------ */

void
setup_hmatrix_aprx_inter_row_bem3d(pbem3d bem, pccluster rc, pccluster cc,
				   pcblock tree, uint m)
{

  (void) rc;
  (void) cc;
  (void) tree;

  assert(bem->kernels->kernel_col != NULL);
  assert(bem->kernels->lagrange_row != NULL);

  setup_interpolation_bem3d(bem->aprx, m);

  bem->farfield_rk = assemble_bem3d_inter_row_rkmatrix;
  bem->farfield_u = NULL;
  bem->farfield_wave_u = NULL;

  bem->leaf_row = NULL;
  bem->leaf_wave_row = NULL;
  bem->leaf_col = NULL;
  bem->leaf_wave_col = NULL;
  bem->transfer_row = NULL;
  bem->transfer_wave_row = NULL;
  bem->transfer_wave_wave_row = NULL;
  bem->transfer_col = NULL;
  bem->transfer_wave_col = NULL;
  bem->transfer_wave_wave_col = NULL;
}

void
setup_hmatrix_aprx_inter_col_bem3d(pbem3d bem, pccluster rc, pccluster cc,
				   pcblock tree, uint m)
{

  (void) rc;
  (void) cc;
  (void) tree;

  assert(bem->kernels->fundamental_row != NULL);
  assert(bem->kernels->lagrange_col != NULL);

  setup_interpolation_bem3d(bem->aprx, m);

  bem->farfield_rk = assemble_bem3d_inter_col_rkmatrix;
  bem->farfield_u = NULL;
  bem->farfield_wave_u = NULL;

  bem->leaf_row = NULL;
  bem->leaf_wave_row = NULL;
  bem->leaf_col = NULL;
  bem->leaf_wave_col = NULL;
  bem->transfer_row = NULL;
  bem->transfer_wave_row = NULL;
  bem->transfer_wave_wave_row = NULL;
  bem->transfer_col = NULL;
  bem->transfer_wave_col = NULL;
  bem->transfer_wave_wave_col = NULL;
}

void
setup_hmatrix_aprx_inter_mixed_bem3d(pbem3d bem, pccluster rc,
				     pccluster cc, pcblock tree, uint m)
{

  (void) rc;
  (void) cc;
  (void) tree;

  assert(bem->kernels->kernel_col != NULL);
  assert(bem->kernels->lagrange_row != NULL);
  assert(bem->kernels->fundamental_row != NULL);
  assert(bem->kernels->lagrange_col != NULL);

  setup_interpolation_bem3d(bem->aprx, m);

  bem->farfield_rk = assemble_bem3d_inter_mixed_rkmatrix;
  bem->farfield_u = NULL;
  bem->farfield_wave_u = NULL;

  bem->leaf_row = NULL;
  bem->leaf_wave_row = NULL;
  bem->leaf_col = NULL;
  bem->leaf_wave_col = NULL;
  bem->transfer_row = NULL;
  bem->transfer_wave_row = NULL;
  bem->transfer_wave_wave_row = NULL;
  bem->transfer_col = NULL;
  bem->transfer_wave_col = NULL;
  bem->transfer_wave_wave_col = NULL;
}

/* ------------------------------------------------------------
 Green
 ------------------------------------------------------------ */

void
setup_hmatrix_aprx_green_row_bem3d(pbem3d bem, pccluster rc, pccluster cc,
				   pcblock tree, uint m, uint l, real delta,
				   quadpoints3d quadpoints)
{

  (void) rc;
  (void) cc;
  (void) tree;

  assert(quadpoints != NULL);
  assert(bem->kernels->fundamental_row != NULL);
  assert(bem->kernels->kernel_col != NULL);
  assert(bem->kernels->dnz_fundamental_row != NULL);
  assert(bem->kernels->dnz_kernel_col != NULL);

  setup_green_bem3d(bem->aprx, m, l, delta, quadpoints);

  bem->farfield_rk = assemble_bem3d_green_row_rkmatrix;
  bem->farfield_u = NULL;
  bem->farfield_wave_u = NULL;

  bem->leaf_row = NULL;
  bem->leaf_wave_row = NULL;
  bem->leaf_col = NULL;
  bem->leaf_wave_col = NULL;
  bem->transfer_row = NULL;
  bem->transfer_wave_row = NULL;
  bem->transfer_wave_wave_row = NULL;
  bem->transfer_col = NULL;
  bem->transfer_wave_col = NULL;
  bem->transfer_wave_wave_col = NULL;
}

void
setup_hmatrix_aprx_green_col_bem3d(pbem3d bem, pccluster rc, pccluster cc,
				   pcblock tree, uint m, uint l, real delta,
				   quadpoints3d quadpoints)
{

  (void) rc;
  (void) cc;
  (void) tree;

  assert(quadpoints != NULL);
  assert(bem->kernels->fundamental_row != NULL);
  assert(bem->kernels->kernel_col != NULL);
  assert(bem->kernels->dnz_fundamental_row != NULL);
  assert(bem->kernels->dnz_kernel_col != NULL);

  setup_green_bem3d(bem->aprx, m, l, delta, quadpoints);

  bem->farfield_rk = assemble_bem3d_green_col_rkmatrix;
  bem->farfield_u = NULL;
  bem->farfield_wave_u = NULL;

  bem->leaf_row = NULL;
  bem->leaf_wave_row = NULL;
  bem->leaf_col = NULL;
  bem->leaf_wave_col = NULL;
  bem->transfer_row = NULL;
  bem->transfer_wave_row = NULL;
  bem->transfer_wave_wave_row = NULL;
  bem->transfer_col = NULL;
  bem->transfer_wave_col = NULL;
  bem->transfer_wave_wave_col = NULL;
}

void
setup_hmatrix_aprx_green_mixed_bem3d(pbem3d bem, pccluster rc,
				     pccluster cc, pcblock tree, uint m,
				     uint l, real delta,
				     quadpoints3d quadpoints)
{

  (void) rc;
  (void) cc;
  (void) tree;

  assert(quadpoints != NULL);
  assert(bem->kernels->fundamental_row != NULL);
  assert(bem->kernels->kernel_col != NULL);
  assert(bem->kernels->dnz_fundamental_row != NULL);
  assert(bem->kernels->dnz_kernel_col != NULL);

  setup_green_bem3d(bem->aprx, m, l, delta, quadpoints);

  bem->farfield_rk = assemble_bem3d_green_mixed_rkmatrix;
  bem->farfield_u = NULL;
  bem->farfield_wave_u = NULL;

  bem->leaf_row = NULL;
  bem->leaf_wave_row = NULL;
  bem->leaf_col = NULL;
  bem->leaf_wave_col = NULL;
  bem->transfer_row = NULL;
  bem->transfer_wave_row = NULL;
  bem->transfer_wave_wave_row = NULL;
  bem->transfer_col = NULL;
  bem->transfer_wave_col = NULL;
  bem->transfer_wave_wave_col = NULL;
}

/* ------------------------------------------------------------
 Greenhybrid
 ------------------------------------------------------------ */

void
setup_hmatrix_aprx_greenhybrid_row_bem3d(pbem3d bem, pccluster rc,
					 pccluster cc, pcblock tree, uint m,
					 uint l, real delta, real accur,
					 quadpoints3d quadpoints)
{
  pparbem3d par = bem->par;

  uint      n, i;

  (void) cc;
  (void) tree;

  assert(quadpoints != NULL);
  assert(bem->kernels->fundamental_row != NULL);
  assert(bem->kernels->dnz_fundamental_row != NULL);
  assert(bem->nearfield_far != NULL);

  setup_green_bem3d(bem->aprx, m, l, delta, quadpoints);
  setup_aca_bem3d(bem->aprx, accur);

  bem->farfield_rk = assemble_bem3d_greenhybrid_row_rkmatrix;
  bem->farfield_u = NULL;
  bem->farfield_wave_u = NULL;

  bem->leaf_row = NULL;
  bem->leaf_wave_row = NULL;
  bem->leaf_col = NULL;
  bem->leaf_wave_col = NULL;
  bem->transfer_row = NULL;
  bem->transfer_wave_row = NULL;
  bem->transfer_wave_wave_row = NULL;
  bem->transfer_col = NULL;
  bem->transfer_wave_col = NULL;
  bem->transfer_wave_wave_col = NULL;

  n = par->grcnn;

  if (par->grcn != NULL && par->grcnn != 0) {
    for (i = 0; i < n; ++i) {
      if (par->grcn[i] != NULL) {
	del_greencluster3d(par->grcn[i]);
      }
    }
    freemem(par->grcn);
  }

  n = rc->desc;

  par->grcn = (pgreencluster3d *) allocmem((size_t) n *
					   sizeof(pgreencluster3d));

  for (i = 0; i < n; ++i) {
    par->grcn[i] = NULL;
  }

  par->grcnn = n;
}

void
setup_hmatrix_aprx_greenhybrid_col_bem3d(pbem3d bem, pccluster rc,
					 pccluster cc, pcblock tree, uint m,
					 uint l, real delta, real accur,
					 quadpoints3d quadpoints)
{
  pparbem3d par = bem->par;

  uint      i, n;

  (void) rc;
  (void) tree;

  assert(quadpoints != NULL);
  assert(bem->kernels->kernel_col != NULL);
  assert(bem->kernels->dnz_kernel_col != NULL);
  assert(bem->nearfield_far != NULL);

  setup_green_bem3d(bem->aprx, m, l, delta, quadpoints);
  setup_aca_bem3d(bem->aprx, accur);

  bem->farfield_rk = assemble_bem3d_greenhybrid_col_rkmatrix;
  bem->farfield_u = NULL;
  bem->farfield_wave_u = NULL;

  bem->leaf_row = NULL;
  bem->leaf_wave_row = NULL;
  bem->leaf_col = NULL;
  bem->leaf_wave_col = NULL;
  bem->transfer_row = NULL;
  bem->transfer_wave_row = NULL;
  bem->transfer_wave_wave_row = NULL;
  bem->transfer_col = NULL;
  bem->transfer_wave_col = NULL;
  bem->transfer_wave_wave_col = NULL;

  n = par->gccnn;

  if (par->gccn != NULL && par->gccnn != 0) {
    for (i = 0; i < n; ++i) {
      if (par->gccn[i] != NULL) {
	del_greencluster3d(par->gccn[i]);
      }
    }
    freemem(par->gccn);
  }

  n = cc->desc;

  par->gccn = (pgreencluster3d *) allocmem((size_t) n *
					   sizeof(pgreencluster3d));

  for (i = 0; i < n; ++i) {
    par->gccn[i] = NULL;
  }

  par->gccnn = n;
}

void
setup_hmatrix_aprx_greenhybrid_mixed_bem3d(pbem3d bem, pccluster rc,
					   pccluster cc, pcblock tree, uint m,
					   uint l, real delta, real accur,
					   quadpoints3d quadpoints)
{
  pparbem3d par = bem->par;

  uint      i, n;

  (void) tree;

  assert(quadpoints != NULL);
  assert(bem->kernels->fundamental_row != NULL);
  assert(bem->kernels->dnz_fundamental_row != NULL);
  assert(bem->kernels->kernel_col != NULL);
  assert(bem->kernels->dnz_kernel_col != NULL);
  assert(bem->nearfield_far != NULL);

  setup_green_bem3d(bem->aprx, m, l, delta, quadpoints);
  setup_aca_bem3d(bem->aprx, accur);

  bem->farfield_rk = assemble_bem3d_greenhybrid_mixed_rkmatrix;
  bem->farfield_u = NULL;
  bem->farfield_wave_u = NULL;

  bem->leaf_row = NULL;
  bem->leaf_wave_row = NULL;
  bem->leaf_col = NULL;
  bem->leaf_wave_col = NULL;
  bem->transfer_row = NULL;
  bem->transfer_wave_row = NULL;
  bem->transfer_wave_wave_row = NULL;
  bem->transfer_col = NULL;
  bem->transfer_wave_col = NULL;
  bem->transfer_wave_wave_col = NULL;

  n = par->grcnn;

  if (par->grcn != NULL && par->grcnn != 0) {
    for (i = 0; i < n; ++i) {
      if (par->grcn[i] != NULL) {
	del_greencluster3d(par->grcn[i]);
      }
    }
    freemem(par->grcn);
  }

  n = rc->desc;

  par->grcn = (pgreencluster3d *) allocmem((size_t) n *
					   sizeof(pgreencluster3d));

  for (i = 0; i < n; ++i) {
    par->grcn[i] = NULL;
  }

  par->grcnn = n;

  n = par->gccnn;

  if (par->gccn != NULL && par->gccnn != 0) {
    for (i = 0; i < n; ++i) {
      if (par->gccn[i] != NULL) {
	del_greencluster3d(par->gccn[i]);
      }
    }
    freemem(par->gccn);
  }

  n = cc->desc;

  par->gccn = (pgreencluster3d *) allocmem((size_t) n *
					   sizeof(pgreencluster3d));

  for (i = 0; i < n; ++i) {
    par->gccn[i] = NULL;
  }

  par->gccnn = n;
}

/* ------------------------------------------------------------
 ACA
 ------------------------------------------------------------ */

void
setup_hmatrix_aprx_aca_bem3d(pbem3d bem, pccluster rc, pccluster cc,
			     pcblock tree, real accur)
{

  (void) rc;
  (void) cc;
  (void) tree;

  assert(bem->nearfield_far != NULL);

  setup_aca_bem3d(bem->aprx, accur);

  bem->farfield_rk = assemble_bem3d_ACA_rkmatrix;
  bem->farfield_u = NULL;
  bem->farfield_wave_u = NULL;

  bem->leaf_row = NULL;
  bem->leaf_wave_row = NULL;
  bem->leaf_col = NULL;
  bem->leaf_wave_col = NULL;
  bem->transfer_row = NULL;
  bem->transfer_wave_row = NULL;
  bem->transfer_wave_wave_row = NULL;
  bem->transfer_col = NULL;
  bem->transfer_wave_col = NULL;
  bem->transfer_wave_wave_col = NULL;
}

void
setup_hmatrix_aprx_paca_bem3d(pbem3d bem, pccluster rc, pccluster cc,
			      pcblock tree, real accur)
{

  (void) rc;
  (void) cc;
  (void) tree;

  assert(bem->nearfield_far != NULL);

  setup_aca_bem3d(bem->aprx, accur);

  bem->farfield_rk = assemble_bem3d_PACA_rkmatrix;
  bem->farfield_u = NULL;
  bem->farfield_wave_u = NULL;

  bem->leaf_row = NULL;
  bem->leaf_wave_row = NULL;
  bem->leaf_col = NULL;
  bem->leaf_wave_col = NULL;
  bem->transfer_row = NULL;
  bem->transfer_wave_row = NULL;
  bem->transfer_wave_wave_row = NULL;
  bem->transfer_col = NULL;
  bem->transfer_wave_col = NULL;
  bem->transfer_wave_wave_col = NULL;
}

/* ------------------------------------------------------------
 HCA
 ------------------------------------------------------------ */

void
setup_hmatrix_aprx_hca_bem3d(pbem3d bem, pccluster rc, pccluster cc,
			     pcblock tree, uint m, real accur)
{

  (void) rc;
  (void) cc;
  (void) tree;

  assert(bem->kernels->fundamental != NULL);
  assert(bem->kernels->kernel_row != NULL);
  assert(bem->kernels->kernel_col != NULL);

  setup_interpolation_bem3d(bem->aprx, m);
  setup_aca_bem3d(bem->aprx, accur);

  bem->farfield_rk = assemble_bem3d_HCA_rkmatrix;
  bem->farfield_u = NULL;
  bem->farfield_wave_u = NULL;

  bem->leaf_row = NULL;
  bem->leaf_wave_row = NULL;
  bem->leaf_col = NULL;
  bem->leaf_wave_col = NULL;
  bem->transfer_row = NULL;
  bem->transfer_wave_row = NULL;
  bem->transfer_wave_wave_row = NULL;
  bem->transfer_col = NULL;
  bem->transfer_wave_col = NULL;
  bem->transfer_wave_wave_col = NULL;
}

/* ------------------------------------------------------------
 Initializerfunctions for h2-matrix approximations
 ------------------------------------------------------------ */

void
setup_h2matrix_recomp_bem3d(pbem3d bem, bool hiercomp, real accur_hiercomp)
{
  paprxbem3d aprx = bem->aprx;

  assert(accur_hiercomp >= 0.0);

  aprx->hiercomp = hiercomp;
  aprx->accur_hiercomp = accur_hiercomp;
  aprx->tm = new_releucl_truncmode();
}

void
setup_h2matrix_aprx_inter_bem3d(pbem3d bem, pcclusterbasis rb,
				pcclusterbasis cb, pcblock tree, uint m)
{

  (void) rb;
  (void) cb;
  (void) tree;

  assert(bem->kernels->lagrange_row != NULL);
  assert(bem->kernels->lagrange_col != NULL);
  assert(bem->kernels->fundamental != NULL);

  setup_interpolation_bem3d(bem->aprx, m);

  bem->farfield_rk = NULL;
  bem->farfield_u = assemble_bem3d_inter_uniform;
  bem->farfield_wave_u = NULL;

  bem->leaf_row = assemble_bem3d_inter_row_clusterbasis;
  bem->leaf_wave_row = NULL;
  bem->leaf_col = assemble_bem3d_inter_col_clusterbasis;
  bem->leaf_wave_col = NULL;
  bem->transfer_row = assemble_bem3d_inter_transfer_row_clusterbasis;
  bem->transfer_wave_row = NULL;
  bem->transfer_wave_wave_row = NULL;
  bem->transfer_col = assemble_bem3d_inter_transfer_col_clusterbasis;
  bem->transfer_wave_col = NULL;
  bem->transfer_wave_wave_col = NULL;
}

void
setup_h2matrix_aprx_greenhybrid_bem3d(pbem3d bem, pcclusterbasis rb,
				      pcclusterbasis cb, pcblock tree, uint m,
				      uint l, real delta, real accur,
				      quadpoints3d quadpoints)
{
  pparbem3d par = bem->par;

  uint      i, n;

  (void) tree;

  assert(quadpoints != NULL);
  assert(bem->kernels->fundamental_row != NULL);
  assert(bem->kernels->dnz_fundamental_row != NULL);
  assert(bem->kernels->kernel_col != NULL);
  assert(bem->kernels->dnz_kernel_col != NULL);
  assert(bem->nearfield_far != NULL);

  setup_green_bem3d(bem->aprx, m, l, delta, quadpoints);
  setup_aca_bem3d(bem->aprx, accur);

  bem->farfield_rk = NULL;
  bem->farfield_u = assemble_bem3d_greenhybrid_uniform;
  bem->farfield_wave_u = NULL;

  bem->leaf_row = assemble_bem3d_greenhybrid_leaf_row_clusterbasis;
  bem->leaf_wave_row = NULL;
  bem->leaf_col = assemble_bem3d_greenhybrid_leaf_col_clusterbasis;
  bem->leaf_wave_col = NULL;
  bem->transfer_row = assemble_bem3d_greenhybrid_transfer_row_clusterbasis;
  bem->transfer_wave_row = NULL;
  bem->transfer_wave_wave_row = NULL;
  bem->transfer_col = assemble_bem3d_greenhybrid_transfer_col_clusterbasis;
  bem->transfer_wave_col = NULL;
  bem->transfer_wave_wave_col = NULL;

  n = par->grbnn;

  if (par->grbn != NULL && par->grbnn != 0) {
    for (i = 0; i < n; ++i) {
      if (par->grbn[i] != NULL) {
	del_greenclusterbasis3d(par->grbn[i]);
      }
    }
    freemem(par->grbn);
  }

  n = rb->t->desc;

  par->grbn = (pgreenclusterbasis3d *) allocmem((size_t) n *
						sizeof(pgreenclusterbasis3d));

  for (i = 0; i < n; ++i) {
    par->grbn[i] = NULL;
  }

  par->grbnn = n;

  n = par->gcbnn;

  if (par->gcbn != NULL && par->gcbnn != 0) {
    for (i = 0; i < n; ++i) {
      if (par->gcbn[i] != NULL) {
	del_greenclusterbasis3d(par->gcbn[i]);
      }
    }
    freemem(par->gcbn);
  }

  n = cb->t->desc;

  par->gcbn = (pgreenclusterbasis3d *) allocmem((size_t) n *
						sizeof(pgreenclusterbasis3d));

  for (i = 0; i < n; ++i) {
    par->gcbn[i] = NULL;
  }

  par->gcbnn = n;
}

void
setup_h2matrix_aprx_greenhybrid_ortho_bem3d(pbem3d bem, pcclusterbasis rb,
					    pcclusterbasis cb, pcblock tree,
					    uint m, uint l, real delta,
					    real accur,
					    quadpoints3d quadpoints)
{
  pparbem3d par = bem->par;

  uint      i, n;

  (void) tree;

  assert(quadpoints != NULL);
  assert(bem->kernels->fundamental_row != NULL);
  assert(bem->kernels->dnz_fundamental_row != NULL);
  assert(bem->kernels->kernel_col != NULL);
  assert(bem->kernels->dnz_kernel_col != NULL);
  assert(bem->nearfield_far != NULL);

  setup_green_bem3d(bem->aprx, m, l, delta, quadpoints);
  setup_aca_bem3d(bem->aprx, accur);

  bem->farfield_rk = NULL;
  bem->farfield_u = assemble_bem3d_greenhybrid_ortho_uniform;
  bem->farfield_wave_u = NULL;

  bem->leaf_row = assemble_bem3d_greenhybrid_ortho_leaf_row_clusterbasis;
  bem->leaf_wave_row = NULL;
  bem->leaf_col = assemble_bem3d_greenhybrid_ortho_leaf_col_clusterbasis;
  bem->leaf_wave_col = NULL;
  bem->transfer_row =
    assemble_bem3d_greenhybrid_ortho_transfer_row_clusterbasis;
  bem->transfer_wave_row = NULL;
  bem->transfer_wave_wave_row = NULL;
  bem->transfer_col =
    assemble_bem3d_greenhybrid_ortho_transfer_col_clusterbasis;
  bem->transfer_wave_col = NULL;
  bem->transfer_wave_wave_col = NULL;

  n = par->grbnn;

  if (par->grbn != NULL && par->grbnn != 0) {
    for (i = 0; i < n; ++i) {
      if (par->grbn[i] != NULL) {
	del_greenclusterbasis3d(par->grbn[i]);
      }
    }
    freemem(par->grbn);
  }

  n = rb->t->desc;

  par->grbn = (pgreenclusterbasis3d *) allocmem((size_t) n *
						sizeof(pgreenclusterbasis3d));

  for (i = 0; i < n; ++i) {
    par->grbn[i] = NULL;
  }

  par->grbnn = n;

  n = par->gcbnn;

  if (par->gcbn != NULL && par->gcbnn != 0) {
    for (i = 0; i < n; ++i) {
      if (par->gcbn[i] != NULL) {
	del_greenclusterbasis3d(par->gcbn[i]);
      }
    }
    freemem(par->gcbn);
  }

  n = cb->t->desc;

  par->gcbn = (pgreenclusterbasis3d *) allocmem((size_t) n *
						sizeof(pgreenclusterbasis3d));

  for (i = 0; i < n; ++i) {
    par->gcbn[i] = NULL;
  }

  par->gcbnn = n;
}

/* ------------------------------------------------------------
 * Initializer functions for DH2-matrix approximations
 * ------------------------------------------------------------ */

HEADER_PREFIX void
setup_dh2matrix_aprx_inter_bem3d(pbem3d bem,
				 pcdclusterbasis rb, pcdclusterbasis cb,
				 pcdblock tree, uint m)
{

  (void) rb;
  (void) cb;
  (void) tree;

  assert(bem->kernels->lagrange_row != NULL);
  assert(bem->kernels->lagrange_wave_row != NULL);
  assert(bem->kernels->lagrange_col != NULL);
  assert(bem->kernels->lagrange_wave_col != NULL);
  assert(bem->kernels->fundamental != NULL);
  assert(bem->kernels->fundamental_wave != NULL);

  setup_interpolation_bem3d(bem->aprx, m);

  bem->farfield_rk = NULL;
  bem->farfield_u = assemble_bem3d_inter_nowave_duniform;
  bem->farfield_wave_u = assemble_bem3d_inter_wave_duniform;

  bem->leaf_row = assemble_bem3d_inter_nowave_leaf_row_dclusterbasis;
  bem->leaf_wave_row = assemble_bem3d_inter_wave_leaf_row_dclusterbasis;
  bem->leaf_col = assemble_bem3d_inter_nowave_leaf_col_dclusterbasis;
  bem->leaf_wave_col = assemble_bem3d_inter_wave_leaf_col_dclusterbasis;
  bem->transfer_row = assemble_bem3d_inter_nowave_transfer_row_dclusterbasis;
  bem->transfer_wave_row =
    assemble_bem3d_inter_wave_transfer_row_dclusterbasis;
  bem->transfer_wave_wave_row =
    assemble_bem3d_inter_wavewave_transfer_row_dclusterbasis;
  bem->transfer_col = assemble_bem3d_inter_nowave_transfer_col_dclusterbasis;
  bem->transfer_wave_col =
    assemble_bem3d_inter_wave_transfer_col_dclusterbasis;
  bem->transfer_wave_wave_col =
    assemble_bem3d_inter_wavewave_transfer_col_dclusterbasis;
}

/* Version for orthogonale cluster basis */

void
setup_dh2matrix_aprx_inter_ortho_bem3d(pbem3d bem, pcdclusterbasis rb,
				       pcdclusterbasis cb, pcdblock tree,
				       uint m)
{

  (void) rb;
  (void) cb;
  (void) tree;

  assert(bem->kernels->lagrange_row != NULL);
  assert(bem->kernels->lagrange_wave_row != NULL);
  assert(bem->kernels->lagrange_col != NULL);
  assert(bem->kernels->lagrange_wave_col != NULL);
  assert(bem->kernels->fundamental != NULL);
  assert(bem->kernels->fundamental_wave != NULL);

  setup_interpolation_bem3d(bem->aprx, m);

  bem->farfield_rk = NULL;
  bem->farfield_u = assemble_bem3d_inter_ortho_nowave_duniform;
  bem->farfield_wave_u = assemble_bem3d_inter_ortho_wave_duniform;

  bem->leaf_row = assemble_bem3d_inter_nowave_leaf_row_dclusterbasis;
  bem->leaf_wave_row = assemble_bem3d_inter_wave_leaf_row_dclusterbasis;
  bem->leaf_col = assemble_bem3d_inter_nowave_leaf_col_dclusterbasis;
  bem->leaf_wave_col = assemble_bem3d_inter_wave_leaf_col_dclusterbasis;
  bem->transfer_row = assemble_bem3d_inter_nowave_transfer_row_dclusterbasis;
  bem->transfer_wave_row =
    assemble_bem3d_inter_wave_transfer_row_dclusterbasis;
  bem->transfer_wave_wave_row =
    assemble_bem3d_inter_wavewave_transfer_row_dclusterbasis;
  bem->transfer_col = assemble_bem3d_inter_nowave_transfer_col_dclusterbasis;
  bem->transfer_wave_col =
    assemble_bem3d_inter_wave_transfer_col_dclusterbasis;
  bem->transfer_wave_wave_col =
    assemble_bem3d_inter_wavewave_transfer_col_dclusterbasis;
}

void
setup_dh2matrix_aprx_inter_recomp_bem3d(pbem3d bem, pcdclusterbasis rb,
					pcdclusterbasis cb, pcdblock tree,
					uint m, ptruncmode tm, real eps)
{

  (void) rb;
  (void) cb;
  (void) tree;

  /* for recompression */
  bem->aprx->recomp = true;
  bem->aprx->accur_recomp = eps;
  bem->aprx->tm = tm;

  assert(bem->kernels->lagrange_row != NULL);
  assert(bem->kernels->lagrange_wave_row != NULL);
  assert(bem->kernels->lagrange_col != NULL);
  assert(bem->kernels->lagrange_wave_col != NULL);
  assert(bem->kernels->fundamental != NULL);
  assert(bem->kernels->fundamental_wave != NULL);

  setup_interpolation_bem3d(bem->aprx, m);

  bem->farfield_rk = NULL;
  bem->farfield_u = assemble_bem3d_inter_ortho_nowave_duniform;
  bem->farfield_wave_u = assemble_bem3d_inter_ortho_wave_duniform;

  bem->leaf_row = assemble_bem3d_inter_nowave_leaf_row_dclusterbasis;
  bem->leaf_wave_row = assemble_bem3d_inter_wave_leaf_row_dclusterbasis;
  bem->leaf_col = assemble_bem3d_inter_nowave_leaf_col_dclusterbasis;
  bem->leaf_wave_col = assemble_bem3d_inter_wave_leaf_col_dclusterbasis;
  bem->transfer_row = assemble_bem3d_inter_nowave_transfer_row_dclusterbasis;
  bem->transfer_wave_row =
    assemble_bem3d_inter_wave_transfer_row_dclusterbasis;
  bem->transfer_wave_wave_row =
    assemble_bem3d_inter_wavewave_transfer_row_dclusterbasis;
  bem->transfer_col = assemble_bem3d_inter_nowave_transfer_col_dclusterbasis;
  bem->transfer_wave_col =
    assemble_bem3d_inter_wave_transfer_col_dclusterbasis;
  bem->transfer_wave_wave_col =
    assemble_bem3d_inter_wavewave_transfer_col_dclusterbasis;
}

/* ------------------------------------------------------------
 Fill hmatrix
 ------------------------------------------------------------ */

void
assemble_bem3d_amatrix(pbem3d bem, pamatrix G)
{
  uint      rows = G->rows;
  uint      cols = G->rows;

  uint     *idx = NULL, *jdx = NULL;
  uint      i;

  idx = allocuint(rows);
  for (i = 0; i < rows; ++i) {
    idx[i] = i;
  }

  if (rows != cols) {
    jdx = allocuint(cols);
    for (i = 0; i < cols; ++i) {
      jdx[i] = i;
    }
  }

  bem->nearfield(idx, jdx, bem, false, G);

  if (rows != cols) {
    freemem(jdx);
  }
  freemem(idx);
}

/* ------------------------------------------------------------
 Fill hmatrix
 ------------------------------------------------------------ */

static void
assemble_bem3d_block_hmatrix(pcblock b, uint bname, uint rname,
			     uint cname, uint pardepth, void *data)
{
  pbem3d    bem = (pbem3d) data;
  paprxbem3d aprx = bem->aprx;
  pparbem3d par = bem->par;
  phmatrix *hn = par->hn;
  phmatrix  G = hn[bname];

  (void) b;
  (void) pardepth;

  if (G->r) {
    bem->farfield_rk(G->rc, rname, G->cc, cname, bem, G->r);
    if (aprx->recomp == true) {
      trunc_rkmatrix(0, aprx->accur_recomp, G->r);
    }
  }
  else if (G->f) {
    bem->nearfield(G->rc->idx, G->cc->idx, bem, false, G->f);
  }
}

static void
assemblecoarsen_bem3d_block_hmatrix(pcblock b, uint bname,
				    uint rname, uint cname, uint pardepth,
				    void *data)
{
  pbem3d    bem = (pbem3d) data;
  paprxbem3d aprx = bem->aprx;
  pparbem3d par = bem->par;
  phmatrix *hn = par->hn;
  phmatrix  G = hn[bname];

  (void) b;
  (void) pardepth;

  if (G->r) {
    bem->farfield_rk(G->rc, rname, G->cc, cname, bem, G->r);
    if (aprx->recomp == true) {
      trunc_rkmatrix(NULL, aprx->accur_recomp, G->r);
    }
  }
  else if (G->f) {
    bem->nearfield(G->rc->idx, G->cc->idx, bem, false, G->f);
  }
  else {
    assert(G->son != NULL);
    coarsen_hmatrix(G, NULL, aprx->accur_coarsen, false);
  }
}

static void
assemble_bem3d_nearfield_block_hmatrix(pcblock b, uint bname,
				       uint rname, uint cname, uint pardepth,
				       void *data)
{
  pbem3d    bem = (pbem3d) data;
  pparbem3d par = bem->par;
  phmatrix *hn = par->hn;
  phmatrix  G = hn[bname];

  (void) b;
  (void) rname;
  (void) cname;
  (void) pardepth;

  if (G->f) {
    bem->nearfield(G->rc->idx, G->cc->idx, bem, false, G->f);
  }
}

static void
assemble_bem3d_farfield_block_hmatrix(pcblock b, uint bname,
				      uint rname, uint cname, uint pardepth,
				      void *data)
{
  pbem3d    bem = (pbem3d) data;
  paprxbem3d aprx = bem->aprx;
  pparbem3d par = bem->par;
  phmatrix *hn = par->hn;
  phmatrix  G = hn[bname];

  (void) b;
  (void) pardepth;

  if (G->r) {
    bem->farfield_rk(G->rc, rname, G->cc, cname, bem, G->r);
    if (aprx->recomp == true) {
      trunc_rkmatrix(0, aprx->accur_recomp, G->r);
    }
  }
}

void
assemble_bem3d_hmatrix(pbem3d bem, pblock b, phmatrix G)
{
  pparbem3d par = bem->par;
  par->hn = enumerate_hmatrix(b, G);

  iterate_byrow_block(b, 0, 0, 0, max_pardepth, NULL,
		      assemble_bem3d_block_hmatrix, bem);

  freemem(par->hn);
  par->hn = NULL;
}

void
assemblecoarsen_bem3d_hmatrix(pbem3d bem, pblock b, phmatrix G)
{
  pparbem3d par = bem->par;
  par->hn = enumerate_hmatrix(b, G);

  iterate_byrow_block(b, 0, 0, 0, max_pardepth, NULL,
		      assemblecoarsen_bem3d_block_hmatrix, bem);

  freemem(par->hn);
  par->hn = NULL;
}

void
assemble_bem3d_nearfield_hmatrix(pbem3d bem, pblock b, phmatrix G)
{
  pparbem3d par = bem->par;
  par->hn = enumerate_hmatrix(b, G);

  iterate_byrow_block(b, 0, 0, 0, max_pardepth, NULL,
		      assemble_bem3d_nearfield_block_hmatrix, bem);

  freemem(par->hn);
  par->hn = NULL;
}

void
assemble_bem3d_farfield_hmatrix(pbem3d bem, pblock b, phmatrix G)
{
  pparbem3d par = bem->par;
  par->hn = enumerate_hmatrix(b, G);

  iterate_byrow_block(b, 0, 0, 0, max_pardepth, NULL,
		      assemble_bem3d_farfield_block_hmatrix, bem);

  freemem(par->hn);
  par->hn = NULL;
}

/* ------------------------------------------------------------
 Fill h2-matrix
 ------------------------------------------------------------ */

static void
assemble_bem3d_block_h2matrix(ph2matrix G, uint bname, uint rname,
			      uint cname, uint pardepth, void *data)
{
  pbem3d    bem = (pbem3d) data;

  (void) pardepth;

  if (G->u) {
    bem->farfield_u(rname, cname, bname, bem);
  }
  else if (G->f) {
    bem->nearfield(G->rb->t->idx, G->cb->t->idx, bem, false, G->f);
  }
}

static void
assemble_bem3d_nearfield_block_h2matrix(ph2matrix G, uint bname,
					uint rname, uint cname, uint pardepth,
					void *data)
{
  pbem3d    bem = (pbem3d) data;

  (void) bname;
  (void) rname;
  (void) cname;
  (void) pardepth;

  if (G->f) {
    bem->nearfield(G->rb->t->idx, G->cb->t->idx, bem, false, G->f);
  }
}

static void
assemble_bem3d_farfield_block_h2matrix(ph2matrix G, uint bname,
				       uint rname, uint cname, uint pardepth,
				       void *data)
{
  pbem3d    bem = (pbem3d) data;

  (void) bname;
  (void) rname;
  (void) cname;
  (void) pardepth;

  if (G->u) {
    bem->farfield_u(rname, cname, bname, bem);
  }
}

static void
assemblehiercomp_bem3d_block_h2matrix(pcblock b, uint bname,
				      uint rname, uint cname, uint pardepth,
				      void *data)
{
  pbem3d    bem = (pbem3d) data;
  pparbem3d par = bem->par;
  paprxbem3d aprx = bem->aprx;
  ph2matrix *h2n = par->h2n;
  ph2matrix G = h2n[bname];
  pclusteroperator *rwn = par->rwn;
  pclusteroperator *cwn = par->cwn;
  uint     *leveln = par->leveln;
  uint      rsons = G->rsons;
  uint      csons = G->csons;
  pccluster rc = G->rb->t;
  pccluster cc = G->cb->t;
  uint      rows = rc->size;
  uint      cols = cc->size;
  ptruncmode tm = aprx->tm;
  real      eps = aprx->accur_hiercomp;

  prkmatrix R;
  pclusteroperator *rw1, *cw1;
  uint      bname1, i, j;

  (void) b;
  (void) pardepth;

  assert(rwn[bname] == 0);
  assert(cwn[bname] == 0);

  if (G->son) {
    rw1 = (pclusteroperator *) allocmem((size_t) sizeof(pclusteroperator) *
					rsons * csons);
    cw1 =
      (pclusteroperator *) allocmem((size_t) sizeof(pclusteroperator) *
				    rsons * csons);

    bname1 = bname + 1;

    for (j = 0; j < csons; j++) {
      for (i = 0; i < rsons; i++) {
	assert(rwn[bname1]);
	assert(cwn[bname1]);

	rw1[i + j * rsons] = rwn[bname1];
	cw1[i + j * rsons] = cwn[bname1];

	bname1 += b->son[i + j * rsons]->desc;
      }
    }

    assert(bname1 == bname + b->desc);

    /* Unify submatrices */
    unify_h2matrix(G, rw1, cw1, tm, eps * pow(tm->zeta_level, leveln[bname]),
		   rwn + bname, cwn + bname);

    /* Clean up */
    for (j = 0; j < csons; j++) {
      for (i = 0; i < rsons; i++) {
	del_clusteroperator(cw1[i + j * rsons]);
	del_clusteroperator(rw1[i + j * rsons]);
      }
    }
    freemem(cw1);
    freemem(rw1);
  }

  if (G->u) {
    R = new_rkmatrix(rows, cols, 0);

    bem->farfield_rk(rc, rname, cc, cname, bem, R);
    convert_rkmatrix_uniform(R, G->u, tm, rwn + bname, cwn + bname);

    ref_clusterbasis(&G->rb, G->u->rb);
    ref_clusterbasis(&G->cb, G->u->cb);

    del_rkmatrix(R);
  }
  else if (G->f) {
    bem->nearfield(rc->idx, cc->idx, bem, false, G->f);

    rwn[bname] = new_leaf_clusteroperator(rc);
    resize_clusteroperator(rwn[bname], 0, G->rb->k);
    cwn[bname] = new_leaf_clusteroperator(cc);
    resize_clusteroperator(cwn[bname], 0, G->cb->k);
  }
}

void
assemble_bem3d_h2matrix(pbem3d bem, ph2matrix G)
{
  bem->par->h2n = enumerate_h2matrix(G);

  iterate_h2matrix(G, 0, 0, 0, max_pardepth, NULL,
		   assemble_bem3d_block_h2matrix, bem);

  freemem(bem->par->h2n);
  bem->par->h2n = NULL;
}

void
assemble_bem3d_nearfield_h2matrix(pbem3d bem, ph2matrix G)
{
  bem->par->h2n = enumerate_h2matrix(G);

  iterate_h2matrix(G, 0, 0, 0, max_pardepth, NULL,
		   assemble_bem3d_nearfield_block_h2matrix, bem);

  freemem(bem->par->h2n);
  bem->par->h2n = NULL;
}

void
assemble_bem3d_farfield_h2matrix(pbem3d bem, ph2matrix G)
{
  bem->par->h2n = enumerate_h2matrix(G);

  iterate_h2matrix(G, 0, 0, 0, max_pardepth, NULL,
		   assemble_bem3d_farfield_block_h2matrix, bem);

  freemem(bem->par->h2n);
  bem->par->h2n = NULL;
}

void
assemblehiercomp_bem3d_h2matrix(pbem3d bem, pblock b, ph2matrix G)
{
  pparbem3d par = bem->par;
  paprxbem3d aprx = bem->aprx;
  uint      blocks = b->desc;

  real      s;
  uint      i;

  par->h2n = enumerate_h2matrix(G);
  par->rwn = (pclusteroperator *) allocmem((size_t) sizeof(pclusteroperator) *
					   blocks);
  par->cwn =
    (pclusteroperator *) allocmem((size_t) sizeof(pclusteroperator) * blocks);
  par->leveln = enumerate_level_block(b);

  for (i = 0; i < blocks; ++i) {
    par->rwn[i] = NULL;
    par->cwn[i] = NULL;
  }

  s = REAL_POW(aprx->tm->zeta_level, getdepth_block(b));
  aprx->accur_hiercomp /= s;

  iterate_byrow_block(b, 0, 0, 0, max_pardepth, NULL,
		      assemblehiercomp_bem3d_block_h2matrix, bem);

  aprx->accur_hiercomp *= s;

  freemem(par->h2n);
  par->h2n = NULL;
  del_clusteroperator(par->rwn[0]);
  del_clusteroperator(par->cwn[0]);
  freemem(par->rwn);
  par->rwn = NULL;
  freemem(par->cwn);
  par->cwn = NULL;
  freemem(par->leveln);
  par->leveln = NULL;
}

static void
assemble_bem3d_h2matrix_cluster_row_clusterbasis(pcclusterbasis rb,
						 uint rname, void *data)
{
  pcbem3d   bem = (pcbem3d) data;

  if (rb->sons > 0) {
    assert(bem->transfer_row);
    bem->transfer_row(rname, bem);
  }
  else {
    assert(bem->leaf_row);
    bem->leaf_row(rname, bem);
  }
}

void
assemble_bem3d_h2matrix_row_clusterbasis(pcbem3d bem, pclusterbasis rb)
{
  bem->par->rbn = enumerate_clusterbasis(rb->t, rb);

  iterate_parallel_clusterbasis((pcclusterbasis) rb, 0, max_pardepth, NULL,
				assemble_bem3d_h2matrix_cluster_row_clusterbasis,
				(void *) bem);

  freemem(bem->par->rbn);
  bem->par->rbn = NULL;
}

static void
assemble_bem3d_h2matrix_cluster_col_clusterbasis(pcclusterbasis cb,
						 uint cname, void *data)
{
  pcbem3d   bem = (pcbem3d) data;

  if (cb->sons > 0) {
    assert(bem->transfer_col);
    bem->transfer_col(cname, bem);
  }
  else {
    assert(bem->leaf_col);
    bem->leaf_col(cname, bem);
  }
}

void
assemble_bem3d_h2matrix_col_clusterbasis(pcbem3d bem, pclusterbasis cb)
{
  bem->par->cbn = enumerate_clusterbasis(cb->t, cb);

  iterate_parallel_clusterbasis((pcclusterbasis) cb, 0, max_pardepth, NULL,
				assemble_bem3d_h2matrix_cluster_col_clusterbasis,
				(void *) bem);

  freemem(bem->par->cbn);
  bem->par->cbn = NULL;
}

/* ------------------------------------------------------------
 Fill DH2-matrix
 ------------------------------------------------------------ */

static void
assemble_bem3d_dh2matrix_cluster_row_dclusterbasis(pdclusterbasis rb,
						   uint rname, uint pardepth,
						   void *data)
{
  pcbem3d   bem = (pcbem3d) data;

  (void) pardepth;

  if (rb->sons > 0) {
    if (rb->t->directions > 0) {
      assert(rb->directions == rb->t->directions);

      bem->transfer_wave_wave_row(rname, bem);

      bem->transfer_wave_row(rname, bem);
    }
    else {
      /* Father and son use standard Lagrange basis */
      assert(rb->directions == 1);

      bem->transfer_row(rname, bem);
    }

  }
  else {
    if (rb->t->directions > 0) {
      /* Cluster uses wave-Lagrange basis */
      assert(rb->directions == rb->t->directions);

      bem->leaf_wave_row(rname, bem);
    }
    else {
      /* Cluster uses standard Lagrange basis */
      assert(rb->directions == 1);

      if (rb->V[0].cols > 0)
	bem->leaf_row(rname, bem);
    }
  }
}

static void
assemble_bem3d_dh2matrix_cluster_col_dclusterbasis(pdclusterbasis cb,
						   uint cname, uint pardepth,
						   void *data)
{
  pcbem3d   bem = (pcbem3d) data;

  (void) pardepth;

  if (cb->sons > 0) {
    if (cb->t->directions > 0) {
      assert(cb->directions == cb->t->directions);

      bem->transfer_wave_wave_col(cname, bem);

      bem->transfer_wave_col(cname, bem);
    }
    else {
      /* Father and son use standard Lagrange basis */
      assert(cb->directions == 1);

      bem->transfer_col(cname, bem);
    }

  }
  else {
    if (cb->t->directions > 0) {
      /* Cluster uses wave-Lagrange basis */
      assert(cb->directions == cb->t->directions);

      bem->leaf_wave_col(cname, bem);
    }
    else {
      /* Cluster uses standard Lagrange basis */
      assert(cb->directions == 1);

      if (cb->V[0].cols > 0)
	bem->leaf_col(cname, bem);
    }
  }
}

/* Orthogonoalize a single cluster basis matrix */

static void
ortho_leaf(pdclusterbasis cb, uint name, bool row, pcbem3d bem)
{

  uint      i;
  pamatrix  Vhat;
  amatrix   tmp1;
  pavector  tau;
  avector   b1;
  uint      rows, dim;
  pdclusteroperator co;

  if (row == true) {		/* Row cluster basis case */
    co = bem->par->ron[name];
  }
  else {			/* Column cluster basis case */
    co = bem->par->con[name];
  }

  rows = cb->t->size;		/* Get size for all V[i] */

  assert(co->t == cb->t);
  for (i = 0; i < cb->directions; i++) {	/* Run through all directions and evaluate QR decomposition */
    Vhat = init_amatrix(&tmp1, rows, cb->k[i]);	/* Set up Vhat rows = size of cluster */
    /* Column = rank of underlying direction */
    copy_amatrix(false, &cb->V[i], Vhat);

    dim = UINT_MIN(rows, cb->k[i]);	/* Set up tau */
    tau = init_avector(&b1, dim);

    qrdecomp_amatrix(Vhat, tau);

    resize_dclusteroperator(co, dim, cb->k[i], i);	/* Change size of clusteroperator (Attention co->V will be destroyed!) */
    setrank_dclusterbasis(cb, i, dim);	/* Change size of clusterbasis */
    copy_upper_amatrix(Vhat, false, &co->C[i]);	/* Reearn R */
    qrexpand_amatrix(Vhat, tau, &cb->V[i]);	/* Reearn Q */

    uninit_avector(tau);
    uninit_amatrix(Vhat);
  }
  update_dclusterbasis(cb);

}

static void
ortho_transfer_row(pdclusterbasis cb, uint rname, pcbem3d bem)
{

  amatrix   tmp1, tmp2;
  avector   b1;
  pamatrix  Vhat, Qhat, Vhat1, Qhat1;
  pavector  tau;
  uint      i, j, roff, rows, dim;
  uint      temp;
  uint      direction = cb->directions;
  uint      dname;
  pdclusteroperator co = bem->par->ron[rname];

  rows = cb->t->size;
  assert(cb->t == co->t);

  /* For every single direction */
  for (j = 0; j < direction; j++) {
    rows = 0;
    dname = rname + 1;
    for (i = 0; i < cb->sons; i++) {
      assert(dname - rname <= co->t->desc);
      rows += bem->par->ron[dname]->krow[cb->dirson[i][j]];
      dname += bem->par->ron[dname]->t->desc;
    }

    Vhat = init_amatrix(&tmp1, rows, cb->k[j]);	/* Initialize matrix with required size */
    clear_amatrix(Vhat);

    roff = 0;
    dname = rname + 1;
    /* Compute new Vhat_t = R_tE_t for every single son and save this matrix as part of the Vhat */
    for (i = 0; i < cb->sons; i++) {
      temp = cb->dirson[i][j];	/* For one term save the needed direction */
      assert(dname - rname <= co->t->desc);

      Vhat1 =
	init_sub_amatrix(&tmp2, Vhat, cb->son[i]->k[temp], roff, cb->k[j], 0);
      clear_amatrix(Vhat1);

      assert(bem->par->ron[dname]->krow[temp] == cb->son[i]->k[temp]);	/* It should be the same size */
      assert(bem->par->ron[dname]->kcol[temp] == cb->E[i][j].rows);
      addmul_amatrix(1.0, false, &bem->par->ron[dname]->C[temp], false,
		     &cb->E[i][j], Vhat1);
      uninit_amatrix(Vhat1);
      roff += cb->son[i]->k[temp];
      dname += bem->par->ron[dname]->t->desc;
    }

    assert(roff == rows);
    dim = UINT_MIN(rows, cb->k[j]);

    tau = init_avector(&b1, dim);
    qrdecomp_amatrix(Vhat, tau);	/* Solve qr decomposition */

    resize_dclusteroperator(co, dim, Vhat->cols, j);	/* Prepare for restructure  Vhat->cols should be cb->k[j] */
    copy_upper_amatrix(Vhat, false, &co->C[j]);	/* Rewrite C[j] with our R_t */

    setrank_dclusterbasis(cb, j, dim);	/* Next step will be Qhat so, that we get the E_t */
    Qhat = init_amatrix(&tmp2, rows, dim);
    qrexpand_amatrix(Vhat, tau, Qhat);
    uninit_amatrix(Vhat);
    uninit_avector(tau);

    roff = 0;			/* Write every single E_t for the direction j of this clusterbasis */
    for (i = 0; i < cb->sons; i++) {
      temp = cb->dirson[i][j];
      Qhat1 =
	init_sub_amatrix(&tmp1, Qhat, cb->son[i]->k[temp], roff, dim, 0);
      assert(cb->E[i][j].rows == cb->son[i]->k[temp]);
      assert(cb->E[i][j].cols == cb->k[j]);
      copy_amatrix(false, Qhat1, &cb->E[i][j]);
      uninit_amatrix(Qhat1);
      roff += cb->son[i]->k[temp];
    }
    assert(roff == rows);

    uninit_amatrix(Qhat);
    update_dclusterbasis(cb);
  }
}

static void
ortho_transfer_col(pdclusterbasis cb, uint cname, pcbem3d bem)
{

  amatrix   tmp1, tmp2;
  avector   b1;
  pamatrix  Vhat, Qhat, Vhat1, Qhat1;
  pavector  tau;
  uint      i, j, roff, rows, dim;
  uint      temp;
  uint      direction = cb->directions;
  uint      dname;
  pdclusteroperator co = bem->par->con[cname];

  rows = cb->t->size;
  assert(cb->t == co->t);

  /* For every single direction */
  for (j = 0; j < direction; j++) {
    rows = 0;
    dname = cname + 1;
    for (i = 0; i < cb->sons; i++) {
      assert(dname - cname <= co->t->desc);
      rows += bem->par->con[dname]->krow[cb->dirson[i][j]];
      dname += bem->par->con[dname]->t->desc;
    }

    Vhat = init_amatrix(&tmp1, rows, cb->k[j]);	/* Initialize matrix with required size */
    clear_amatrix(Vhat);

    roff = 0;
    dname = cname + 1;
    /* Compute new Vhat_t = R_tE_t for every single son and save this matrix as part of the Vhat */
    for (i = 0; i < cb->sons; i++) {
      temp = cb->dirson[i][j];	/* For one term save the needed direction */
      assert(dname - cname <= co->t->desc);

      Vhat1 =
	init_sub_amatrix(&tmp2, Vhat, cb->son[i]->k[temp], roff, cb->k[j], 0);
      clear_amatrix(Vhat1);
      assert(bem->par->con[dname]->krow[temp] == cb->son[i]->k[temp]);	/* It should be the same size */
      assert(bem->par->con[dname]->kcol[temp] == cb->E[i][j].rows);
      addmul_amatrix(1.0, false, &bem->par->con[dname]->C[temp], false,
		     &cb->E[i][j], Vhat1);
      uninit_amatrix(Vhat1);
      roff += cb->son[i]->k[temp];
      dname += bem->par->con[dname]->t->desc;
    }

    assert(roff == rows);
    dim = UINT_MIN(rows, cb->k[j]);

    tau = init_avector(&b1, dim);
    qrdecomp_amatrix(Vhat, tau);	/* Solve qr decomposition */

    resize_dclusteroperator(co, dim, Vhat->cols, j);	/* Prepare for restructure  Vhat->cols should be cb->k[j] */
    copy_upper_amatrix(Vhat, false, &co->C[j]);	/* Rewrite C[j] with our R_t */

    setrank_dclusterbasis(cb, j, dim);	/* Next step will be Qhat so, that we get the E_t */
    Qhat = init_amatrix(&tmp2, rows, dim);
    qrexpand_amatrix(Vhat, tau, Qhat);
    uninit_amatrix(Vhat);
    uninit_avector(tau);

    roff = 0;			/* Write every single E_t for the direction j of this clusterbasis */
    for (i = 0; i < cb->sons; i++) {
      temp = cb->dirson[i][j];
      Qhat1 =
	init_sub_amatrix(&tmp1, Qhat, cb->son[i]->k[temp], roff, dim, 0);
      assert(cb->E[i][j].rows == cb->son[i]->k[temp]);
      assert(cb->E[i][j].cols == cb->k[j]);
      copy_amatrix(false, Qhat1, &cb->E[i][j]);
      uninit_amatrix(Qhat1);
      roff += cb->son[i]->k[temp];
    }
    assert(roff == rows);

    uninit_amatrix(Qhat);
    update_dclusterbasis(cb);
  }
}

static void
assemble_bem3d_dh2matrix_ortho_cluster_row_dclusterbasis(pdclusterbasis rb,
							 uint rname,
							 uint pardepth,
							 void *data)
{

  pcbem3d   bem = (pcbem3d) data;

  (void) pardepth;

  if (rb->sons > 0) {
    if (rb->t->directions > 0) {
      assert(rb->directions == rb->t->directions);

      bem->transfer_wave_wave_row(rname, bem);

      bem->transfer_wave_row(rname, bem);
    }
    else {
      /* Father and son use standard Lagrange basis */
      assert(rb->directions == 1);

      bem->transfer_row(rname, bem);
    }

    ortho_transfer_row(rb, rname, bem);
  }
  else {			/*Leaf case */
    if (rb->t->directions > 0) {
      /* Cluster uses wave-Lagrange basis */
      assert(rb->directions == rb->t->directions);

      bem->leaf_wave_row(rname, bem);
    }
    else {
      /* Cluster uses standard Lagrange basis */
      assert(rb->directions == 1);

      if (rb->V[0].cols > 0) {
	bem->leaf_row(rname, bem);
      }
    }
    ortho_leaf(rb, rname, true, bem);
  }
}

static void
assemble_bem3d_dh2matrix_ortho_cluster_col_dclusterbasis(pdclusterbasis cb,
							 uint cname,
							 uint pardepth,
							 void *data)
{

  pcbem3d   bem = (pcbem3d) data;

  (void) pardepth;

  if (cb->sons > 0) {
    if (cb->t->directions > 0) {
      assert(cb->directions == cb->t->directions);

      bem->transfer_wave_wave_col(cname, bem);

      bem->transfer_wave_col(cname, bem);
    }
    else {
      /* Father and son use standard Lagrange basis */
      assert(cb->directions == 1);

      bem->transfer_col(cname, bem);
    }
    ortho_transfer_col(cb, cname, bem);
  }
  else {			/*Leaf case */
    if (cb->t->directions > 0) {
      /* Cluster uses wave-Lagrange basis */
      assert(cb->directions == cb->t->directions);

      bem->leaf_wave_col(cname, bem);
    }
    else {
      /* Cluster uses standard Lagrange basis */
      assert(cb->directions == 1);

      if (cb->V[0].cols > 0) {
	bem->leaf_col(cname, bem);
      }
    }
    ortho_leaf(cb, cname, false, bem);
  }
}

void
assemble_bem3d_dh2matrix_row_dclusterbasis(pcbem3d bem, pdclusterbasis rb)
{
  pparbem3d par = bem->par;
  par->drbn = enumerate_dclusterbasis(rb->t, rb);

  iterate_dclusterbasis(rb, 0, max_pardepth, NULL,
			assemble_bem3d_dh2matrix_cluster_row_dclusterbasis,
			(void *) bem);

  freemem(par->drbn);
  par->drbn = NULL;
}

void
assemble_bem3d_dh2matrix_col_dclusterbasis(pcbem3d bem, pdclusterbasis cb)
{
  pparbem3d par = bem->par;
  par->dcbn = enumerate_dclusterbasis(cb->t, cb);

  iterate_dclusterbasis(cb, 0, max_pardepth, NULL,
			assemble_bem3d_dh2matrix_cluster_col_dclusterbasis,
			(void *) bem);

  freemem(par->dcbn);
  par->dcbn = NULL;
}

void
assemble_bem3d_dh2matrix_ortho_row_dclusterbasis(pcbem3d bem,
						 pdclusterbasis rb,
						 pdclusteroperator ro)
{

  pparbem3d par = bem->par;
  par->drbn = enumerate_dclusterbasis(rb->t, rb);
  par->ron = enumerate_dclusteroperator(rb->t, ro);

  iterate_dclusterbasis(rb, 0, max_pardepth, NULL,
			assemble_bem3d_dh2matrix_cluster_row_dclusterbasis,
			(void *) bem);

  iterate_dclusterbasis(rb, 0, max_pardepth, NULL,
			assemble_bem3d_dh2matrix_ortho_cluster_row_dclusterbasis,
			(void *) bem);

  freemem(par->drbn);
  par->drbn = NULL;
}

void
assemble_bem3d_dh2matrix_ortho_col_dclusterbasis(pcbem3d bem,
						 pdclusterbasis cb,
						 pdclusteroperator co)
{

  pparbem3d par = bem->par;
  par->dcbn = enumerate_dclusterbasis(cb->t, cb);
  par->con = enumerate_dclusteroperator(cb->t, co);

  iterate_dclusterbasis(cb, 0, max_pardepth, NULL,
			assemble_bem3d_dh2matrix_cluster_col_dclusterbasis,
			(void *) bem);

  iterate_dclusterbasis(cb, 0, max_pardepth, NULL,
			assemble_bem3d_dh2matrix_ortho_cluster_col_dclusterbasis,
			(void *) bem);

  freemem(par->dcbn);
  par->dcbn = NULL;
}

/*---------------------------------------------------------------------
 Stuff for parallel recompression
 ----------------------------------------------------------------------*/

typedef struct _admisblock admisblock;

struct _admisblock {
  uint      name;		/* number of the corresponding block */
  uint      rname;		/* number of row cluster */
  uint      cname;		/* number of col cluster */
  uint      father;		/* number of father */
  uint      son;		/* son number */
  uint      length;		/* length of the list, only updated in first list entry */
  struct _admisblock *next;	/* Next admissible block */
};

typedef struct _compdata compdata;
typedef compdata *pcompdata;
struct _compdata {

  pdclusteroperator *nco;
  pdclusteroperator *nro;

  pdclusteroperator *noro;
  pdclusteroperator *noco;
  pdclusterbasis *ncb;
  pdclusterbasis *nrb;
  pdblock  *nb;

  pcbem3d   bem;
  bool      rows;

  admisblock **cblock;
  admisblock **rblock;

};

static admisblock *
create_newadmisblock(uint bname, uint rname, uint cname, uint father)
{

  admisblock *ab = (admisblock *) allocmem(sizeof(admisblock));

  ab->name = bname;
  ab->rname = rname;
  ab->cname = cname;
  ab->father = father;
  ab->next = NULL;
  return ab;
}

static void
create_weight_lists(uint bname, uint rname, uint cname,
		    uint rfather, uint cfather, uint rson, uint cson,
		    void *data, admisblock ** rlist, admisblock ** clist)
{

  pcompdata cdata = (pcompdata) data;
  uint      i, j;
  uint      bname1, rname1, cname1;
  pdblock   b = cdata->nb[bname];
  admisblock *abr = *rlist;
  admisblock *abc = *clist;

  if (b->csons + b->rsons > 0) {

    abr->father = rfather;
    abc->father = cfather;
    abr->son = rson;
    abc->son = cson;
    bname1 = bname + 1;
    cname1 = (b->csons > 0 ? cname + 1 : cname);

    for (j = 0; j < b->csons; j++) {
      rname1 = (b->rsons > 0 ? rname + 1 : rname);
      for (i = 0; i < b->rsons; i++) {
	create_weight_lists(bname1, rname1, cname1, rname, cname, i, j, data,
			    &(cdata->rblock[rname1]),
			    &(cdata->cblock[cname1]));

	rname1 += cdata->nb[bname1]->rc->desc;
	bname1 += cdata->nb[bname1]->desc;
      }
      assert(rname1 == rname + b->rc->desc);
      cname1 += b->son[j * b->rsons]->cc->desc;
    }
    assert(cname1 == cname + b->cc->desc);
    assert(bname1 == bname + b->desc);
  }
  else {
    abr->father = rfather;
    abc->father = cfather;
    abr->son = rson;
    abc->son = cson;

    if (b->adm == true) {
      admisblock *newr = create_newadmisblock(bname, rname, cname, rfather);
      /* save number bname in list for rows */
      if (abr->length == 0) {
	abr->name = bname;
	abr->rname = rname;
	abr->cname = cname;
	abr->father = rfather;
	abr->son = rson;
	abr->length = 1;
      }
      else {
	while (abr->next != NULL) {
	  abr = abr->next;
	}
	abr->next = newr;
	cdata->rblock[rname]->length += 1;
      }
      /* save number bname in list for cols */
      admisblock *newc = create_newadmisblock(bname, rname, cname, cfather);
      if (abc->length == 0) {
	abc->name = bname;
	abc->rname = rname;
	abc->cname = cname;
	abc->father = cfather;
	abc->son = cson;
	abc->length = 1;
      }
      else {
	while (abc->next != NULL) {
	  abc = abc->next;
	}
	abc->next = newc;
	cdata->cblock[cname]->length += 1;
      }
    }
  }
}

static void
local_weight(pdclusterbasis cb, uint tname, uint pardepth, void *data)
{

  pcompdata cdata = (pcompdata) data;
  pdclusteroperator yco = cdata->nco[tname];
  pdclusteroperator yro = cdata->nro[tname];
  pdclusteroperator oco, oro;
  pcbem3d   bem = cdata->bem;
  pdblock   b;
  bool      rows = cdata->rows;
  admisblock *ab;
  uint    **information;
  uint      i, j;
  uint      n, off;
  uint      all, size, length;
  uint      rd, cd;
  amatrix   tmp1, tmp2, tmp3, tmp4;
  pamatrix  What, What1;
  uint      dim, m;
  avector   b1;
  pavector  tau;
  real      alpha, norm;
  pamatrix  S, Shat;
  pamatrix  Yhat, Yhat1;
  real      zeta_age;
  pkernelbem3d kernels = bem->kernels;
  uint      k = bem->aprx->k_inter;
  real(*xi_r)[3], (*xi_c)[3];
  if (rows == true) {
    ab = cdata->rblock[tname];
  }
  else {
    ab = cdata->cblock[tname];
  }
  if (ab->length > 0) {
    oro = cdata->noro[ab->rname];
    oco = cdata->noco[ab->cname];
    alpha = 1.0;
    if (rows == true) {		/* Row version */

      length = ab->length;
      information = (uint **) malloc(sizeof(uint *) * (length));
      for (j = 0; j < length; j++) {
	information[j] = (uint *) malloc(sizeof(uint) * 2);
	information[j][0] = 0;	/* direction */
	information[j][1] = 0;	/* all clusteroperators or basis together ....for counting */
      }
      all = 1;			/* number of different cluster for computing weights */
      /* First entry */
      information[0][0] = cdata->nb[ab->name]->rd;
      information[0][1] = oco->krow[cdata->nb[ab->name]->cd];

      while (ab->next != NULL) {
	ab = ab->next;
	/* Find direction if already used or next empty array */
	j = 0;
	while ((information[j][1] > 0)
	       && (information[j][0] != cdata->nb[ab->name]->cd)) {
	  j += 1;
	  assert(j < length);
	}
	oco = cdata->noco[ab->cname];
	information[j][0] = cdata->nb[ab->name]->rd;
	information[j][1] += oco->krow[cdata->nb[ab->name]->cd];
	all += 1;
      }
      /* Set up matrices */
      j = 0;
      while (all > 0) {
	/* size of matrix */
	rd = information[j][0];
	size = information[j][1];
	off = 0;
	ab = cdata->rblock[tname];
	/* Compute coupling matrix and local weight */
	while ((off < size) && (ab != NULL)) {
	  if (cdata->nb[ab->name]->rd == rd) {	/* element of list belongs to this direction */
	    /* set up */
	    cd = cdata->nb[ab->name]->cd;
	    b = cdata->nb[ab->name];
	    oco = cdata->noco[ab->cname];
	    assert(b->cd == cd);
	    assert(b->rd == rd);
	    assert(tname == ab->rname);
	    n = yro->krow[rd] + oco->krow[cd];
	    What = init_amatrix(&tmp1, n, cb->k[rd]);	/* In this case cb is the row cluster basis */
	    clear_amatrix(What);

	    /* Compute S */
	    S = init_amatrix(&tmp3, k, k);
	    clear_amatrix(S);

	    /* Find interpolation points for the row cluster */
	    xi_r = (real(*)[3]) allocreal(3 * k);
	    assemble_interpoints3d_array(bem, b->rc->bmin, b->rc->bmax, xi_r);

	    /* Find interpolation points for the column cluster */
	    xi_c = (real(*)[3]) allocreal(3 * k);
	    assemble_interpoints3d_array(bem, b->cc->bmin, b->cc->bmax, xi_c);

	    if (cb->t->directions > 0) {
	      assert(cdata->ncb[ab->cname]->t->directions > 0);
	      kernels->fundamental_wave(bem, (const real(*)[3]) xi_r,
					(const real(*)[3]) xi_c,
					b->rc->dir[rd], S);
	    }
	    else {
	      kernels->fundamental(bem, (const real(*)[3]) xi_r,
				   (const real(*)[3]) xi_c, S);
	    }
	    Shat = init_amatrix(&tmp4, cb->k[rd], k);
	    clear_amatrix(Shat);
	    addmul_amatrix(1.0, false, &oro->C[rd], false, S, Shat);
	    uninit_amatrix(S);
	    S = init_amatrix(&tmp3, cb->k[rd], cdata->ncb[ab->cname]->k[cd]);
	    clear_amatrix(S);
	    addmul_amatrix(1.0, false, Shat, true, &oco->C[cd], S);
	    uninit_amatrix(Shat);

	    alpha = 1.0;
	    /* Factor for Truncation, if necessary */
	    if (bem->aprx->tm && bem->aprx->tm->blocks) {
	      if (bem->aprx->tm->frobenius) {	// case: Frobenius norm
		norm = normfrob_amatrix(S);
	      }
	      else {		// case: Euclidean norm
		norm = norm2_amatrix(S);
	      }
	      alpha = (norm > 0.0 ? 1.0 / norm : 1.0);
	    }
	    /* If we already have something */
	    if (yro->krow[rd] > 0) {
	      What1 = init_sub_amatrix(&tmp2, What, yro->krow[rd], 0,
				       yro->C[rd].cols, 0);
	      clear_amatrix(What1);
	      copy_amatrix(false, &yro->C[rd], What1);
	      uninit_amatrix(What1);
	    }
	    What1 =
	      init_sub_amatrix(&tmp2, What, oco->krow[cd], yro->krow[rd],
			       cb->k[rd], 0);
	    clear_amatrix(What1);
	    add_amatrix(alpha, true, S, What1);
	    uninit_amatrix(What1);

	    uninit_amatrix(S);
	    off += oco->krow[cd];

	    /* QR */
	    dim = UINT_MIN(cb->k[rd], n);
	    tau = init_avector(&b1, dim);
	    qrdecomp_amatrix(What, tau);
	    uninit_avector(tau);

	    resize_dclusteroperator(yro, dim, cb->k[rd], rd);
	    copy_upper_amatrix(What, false, &yro->C[rd]);
	    uninit_amatrix(What);
	  }
	  ab = ab->next;
	}
	/* prepare for next round */
	all -= 1;
	j += 1;			/* next element of array */
      }
    }
    else {
      length = ab->length;
      information = (uint **) malloc(sizeof(uint *) * (length));
      for (j = 0; j < length; j++) {
	information[j] = (uint *) malloc(sizeof(uint) * 2);
	information[j][0] = 0;	/* direction */
	information[j][1] = 0;	/* size clusteroperator or clusterbasis */
      }
      all = 1;			/* number of different cluster for computing weights */

      /* First entry */
      information[0][0] = cdata->nb[ab->name]->cd;
      information[0][1] = oro->krow[cdata->nb[ab->name]->rd];

      while (ab->next != NULL) {
	ab = ab->next;
	/* Find direction if already used or next empty array */
	j = 0;
	while ((information[j][1] > 0)
	       && (information[j][0] != cdata->nb[ab->name]->rd)) {
	  j += 1;
	  assert(j < length);
	}
	oro = cdata->noro[ab->rname];
	information[j][0] = cdata->nb[ab->name]->cd;
	information[j][1] += oro->krow[cdata->nb[ab->name]->rd];
	all += 1;
      }
      /* Set up matrices */
      j = 0;
      while (all > 0) {
	/* size of matrix */
	cd = information[j][0];
	size = information[j][1];
	off = 0;
	ab = cdata->cblock[tname];
	/* Compute */
	while ((off < size) && (ab != NULL)) {
	  if (cdata->nb[ab->name]->cd == cd) {	/* element of list belongs to this direction */
	    b = cdata->nb[ab->name];
	    rd = b->rd;
	    oro = cdata->noro[ab->rname];
	    assert(b->cd == cd);
	    assert(b->rd == rd);
	    assert(tname == ab->cname);
	    n = oro->krow[rd] + yco->krow[cd];
	    What = init_amatrix(&tmp1, n, cb->k[cd]);
	    clear_amatrix(What);
	    /* Compute S */
	    S = init_amatrix(&tmp3, k, k);
	    clear_amatrix(S);

	    /* Find interpolation points for the row cluster */
	    xi_r = (real(*)[3]) allocreal(3 * k);
	    assemble_interpoints3d_array(bem, b->rc->bmin, b->rc->bmax, xi_r);

	    /* Find interpolation points for the column cluster */
	    xi_c = (real(*)[3]) allocreal(3 * k);
	    assemble_interpoints3d_array(bem, b->cc->bmin, b->cc->bmax, xi_c);

	    if (cdata->nrb[ab->rname]->t->directions > 0) {
	      assert(cb->t->directions > 0);
	      kernels->fundamental_wave(bem, (const real(*)[3]) xi_r,
					(const real(*)[3]) xi_c,
					b->rc->dir[rd], S);
	    }
	    else {
	      kernels->fundamental(bem, (const real(*)[3]) xi_r,
				   (const real(*)[3]) xi_c, S);
	    }

	    Shat = init_amatrix(&tmp4, cdata->nrb[ab->rname]->k[rd], k);
	    clear_amatrix(Shat);
	    addmul_amatrix(1.0, false, &oro->C[rd], false, S, Shat);
	    uninit_amatrix(S);
	    S = init_amatrix(&tmp3, cdata->nrb[ab->rname]->k[rd], cb->k[cd]);
	    clear_amatrix(S);
	    addmul_amatrix(1.0, false, Shat, true, &oco->C[cd], S);
	    uninit_amatrix(Shat);
	    alpha = 1.0;
	    /* Factor for Truncation, if necessary */

	    if (bem->aprx->tm && bem->aprx->tm->blocks) {
	      if (bem->aprx->tm->frobenius) {	// case Frobenius norm
		norm = normfrob_amatrix(S);
	      }
	      else {		// case Euclidean norm
		norm = norm2_amatrix(S);
	      }
	      alpha = (norm > 0.0 ? 1.0 / norm : 1.0);
	    }

	    if (yco->krow[cd] > 0) {
	      What1 =
		init_sub_amatrix(&tmp2, What, yco->krow[cd], 0, cb->k[cd], 0);
	      copy_amatrix(false, &yco->C[cd], What1);
	      uninit_amatrix(What1);
	    }

	    What1 =
	      init_sub_amatrix(&tmp2, What, oro->krow[rd], yco->krow[cd],
			       cb->k[cd], 0);
	    clear_amatrix(What1);
	    add_amatrix(alpha, false, S, What1);
	    uninit_amatrix(What1);

	    off += oro->krow[rd];
	    uninit_amatrix(S);

	    /* QR */
	    dim = UINT_MIN(cb->k[cd], n);
	    tau = init_avector(&b1, dim);
	    qrdecomp_amatrix(What, tau);
	    uninit_avector(tau);

	    resize_dclusteroperator(yco, dim, cb->k[cd], cd);
	    copy_upper_amatrix(What, false, &yco->C[cd]);
	    uninit_amatrix(What);
	  }
	  ab = ab->next;
	}
	/* prepare for next round */
	all -= 1;
	j += 1;			/* next element of array */
      }
    }
    /*cleaning up */
    for (j = 0; j < length; j++) {
      freemem(information[j]);
    }
    freemem(information);
  }

  /* Now combine to total weights */
  if ((rows == true) && (tname != 0)) {
    ab = cdata->rblock[tname];
    zeta_age = (bem->aprx->tm ? bem->aprx->tm->zeta_age : 1.0);
    for (i = 0; i < cb->directions; i++) {
      /* Finding all corresponding father directions and save size */
      m = 0;
      for (j = 0; j < cdata->nrb[ab->father]->directions; j++) {
	m = (cdata->nrb[ab->father]->dirson[ab->son][j] == i ?
	     m + cdata->nro[ab->father]->C[j].rows : m);
      }
      /* Set up matrix */
      Yhat = init_amatrix(&tmp1, m + yro->C[i].rows, yro->C[i].cols);
      /* Evaluate ZE* from father and collect in upper part of Yhat */
      off = 0;
      if (m > 0) {
	for (j = 0; j < cdata->nrb[ab->father]->directions; j++) {
	  if (cdata->nrb[ab->father]->dirson[ab->son][j] == i) {
	    Yhat1 = init_sub_amatrix(&tmp2, Yhat,
				     cdata->nro[ab->father]->C[j].rows, off,
				     cdata->nrb[ab->father]->E[ab->son][j].
				     rows, 0);
	    clear_amatrix(Yhat1);
	    assert(cdata->nro[ab->father]->C[j].cols ==
		   cdata->nrb[ab->father]->E[ab->son][j].cols);
	    addmul_amatrix(zeta_age, false, &cdata->nro[ab->father]->C[j],
			   true, &cdata->nrb[ab->father]->E[ab->son][j],
			   Yhat1);
	    off += cdata->nro[ab->father]->C[j].rows;
	    uninit_amatrix(Yhat1);
	  }
	}
	assert(off == m);
      }
      Yhat1 = init_sub_amatrix(&tmp2, Yhat, yro->C[i].rows, m, yro->C[i].cols,
			       0);
      copy_amatrix(false, &yro->C[i], Yhat1);
      uninit_amatrix(Yhat1);

      dim = UINT_MIN(m + yro->C[i].rows, yro->C[i].cols);
      tau = init_avector(&b1, dim);
      qrdecomp_amatrix(Yhat, tau);
      uninit_avector(tau);

      resize_dclusteroperator(yro, dim, yro->C[i].cols, i);
      copy_upper_amatrix(Yhat, false, &yro->C[i]);
      uninit_amatrix(Yhat);
    }
  }
  else {
    if (tname != 0) {
      zeta_age = (bem->aprx->tm ? bem->aprx->tm->zeta_age : 1.0);
      ab = cdata->cblock[tname];
      /* For all his direction take a look at the cluster operator */
      for (i = 0; i < cb->directions; i++) {

	/* Finding all father directions and save size */
	m = 0;
	for (j = 0; j < cdata->ncb[ab->father]->directions; j++) {
	  m = (cdata->ncb[ab->father]->dirson[ab->son][j] == i ?
	       m + cdata->nco[ab->father]->C[j].rows : m);
	}
	/* Set up matrix */
	Yhat = init_amatrix(&tmp1, m + yco->C[i].rows, yco->C[i].cols);
	/* Evaluate ZE* from father and collect in upper part of Yhat */
	off = 0;
	if (m > 0) {
	  for (j = 0; j < cdata->ncb[ab->father]->directions; j++) {
	    if (cdata->ncb[ab->father]->dirson[ab->son][j] == i) {
	      Yhat1 = init_sub_amatrix(&tmp2, Yhat,
				       cdata->nco[ab->father]->C[j].rows, off,
				       cdata->ncb[ab->father]->E[ab->son][j].
				       rows, 0);
	      clear_amatrix(Yhat1);
	      assert(cdata->nco[ab->father]->C[j].cols ==
		     cdata->ncb[ab->father]->E[ab->son][j].cols);
	      addmul_amatrix(zeta_age, false, &cdata->nco[ab->father]->C[j],
			     true, &cdata->ncb[ab->father]->E[ab->son][j],
			     Yhat1);
	      off += cdata->nco[ab->father]->C[j].rows;
	      uninit_amatrix(Yhat1);
	    }
	  }
	  assert(off == m);
	}
	Yhat1 =
	  init_sub_amatrix(&tmp2, Yhat, yco->C[i].rows, m, yco->C[i].cols, 0);
	copy_amatrix(false, &yco->C[i], Yhat1);
	uninit_amatrix(Yhat1);

	dim = UINT_MIN(m + yco->C[i].rows, yco->C[i].cols);
	tau = init_avector(&b1, dim);
	qrdecomp_amatrix(Yhat, tau);
	uninit_avector(tau);

	resize_dclusteroperator(yco, dim, yco->C[i].cols, i);
	copy_upper_amatrix(Yhat, false, &yco->C[i]);
	uninit_amatrix(Yhat);
      }
    }
  }
}

static void
recomp_dweights(pcdblock b, pdclusterbasis rb, pdclusterbasis cb,
		pdclusteroperator oro, pdclusteroperator oco,
		pdclusteroperator yro, pdclusteroperator yco, pcbem3d bem)
{

  compdata  cdata;
  uint      i;
  uint      maxcname = cb->t->desc;
  uint      maxrname = rb->t->desc;
  admisblock **rab;
  admisblock **cab;

  rab = (admisblock **) allocmem(sizeof(admisblock *) * maxrname);
  for (i = 0; i < maxrname; i++) {
    rab[i] = (admisblock *) allocmem(sizeof(admisblock));
    rab[i]->next = NULL;
    rab[i]->length = 0;
  }

  cab = (admisblock **) allocmem(sizeof(admisblock *) * maxcname);
  for (i = 0; i < maxcname; i++) {
    cab[i] = (admisblock *) allocmem(sizeof(admisblock));
    cab[i]->next = NULL;
    cab[i]->length = 0;
  }

  cdata.nb = enumerate_dblock((pdblock) b);
  cdata.bem = bem;
  cdata.rblock = rab;
  cdata.cblock = cab;

  create_weight_lists(0, 0, 0, 0, 0, 0, 0, (void *) &cdata, &rab[0], &cab[0]);
  if (oro) {
    cdata.noro = enumerate_dclusteroperator(rb->t, oro);
  }
  else {
    cdata.noro = NULL;
  }
  if (oco) {
    cdata.noco = enumerate_dclusteroperator(cb->t, oco);
  }
  else {
    cdata.noco = NULL;
  }

  cdata.nco = enumerate_dclusteroperator(cb->t, yco);
  cdata.nro = enumerate_dclusteroperator(rb->t, yro);
  cdata.nrb = enumerate_dclusterbasis(rb->t, rb);
  cdata.ncb = enumerate_dclusterbasis(cb->t, cb);

  /* Compute weights */
  cdata.rows = true;
  iterate_dclusterbasis(rb, 0, max_pardepth, local_weight, NULL,
			(void *) &cdata);
  cdata.rows = false;
  iterate_dclusterbasis(cb, 0, max_pardepth, local_weight, NULL,
			(void *) &cdata);

  freemem(cdata.nco);
  freemem(cdata.noco);
  freemem(cdata.nro);
  freemem(cdata.noro);
  freemem(cdata.nb);

  for (i = 0; i < maxcname; i++) {
    freemem(cab[i]);
  }
  freemem(cab);

  for (i = 0; i < maxrname; i++) {
    freemem(rab[i]);
  }
  freemem(rab);
}

void
assemble_bem3d_dh2matrix_recomp_both_dclusterbasis(pcbem3d bem,
						   pdclusterbasis rb,
						   pdclusteroperator bro,
						   pdclusterbasis cb,
						   pdclusteroperator bco,
						   pcdblock broot)
{

  pparbem3d par = bem->par;
  pdclusteroperator oro, oco;
  pdclusteroperator yco, yro;
  pdclusterbasis trb, tcb;

  /* row cluster basis */
  oro = build_from_dclusterbasis_dclusteroperator(rb);
  par->drbn = enumerate_dclusterbasis(rb->t, rb);
  par->ron = enumerate_dclusteroperator(rb->t, oro);

  iterate_dclusterbasis(rb, 0, max_pardepth, NULL,
			assemble_bem3d_dh2matrix_cluster_row_dclusterbasis,
			(void *) bem);

  iterate_dclusterbasis(rb, 0, max_pardepth, NULL,
			assemble_bem3d_dh2matrix_ortho_cluster_row_dclusterbasis,
			(void *) bem);

  freemem(par->drbn);
  par->drbn = NULL;
  freemem(par->ron);
  par->ron = NULL;

  /* column cluster basis */
  oco = build_from_dclusterbasis_dclusteroperator(cb);
  par->dcbn = enumerate_dclusterbasis(cb->t, cb);
  par->con = enumerate_dclusteroperator(cb->t, oco);

  iterate_dclusterbasis(cb, 0, max_pardepth, NULL,
			assemble_bem3d_dh2matrix_cluster_col_dclusterbasis,
			(void *) bem);

  iterate_dclusterbasis(cb, 0, max_pardepth, NULL,
			assemble_bem3d_dh2matrix_ortho_cluster_col_dclusterbasis,
			(void *) bem);

  freemem(par->dcbn);
  par->dcbn = NULL;
  freemem(par->con);
  par->con = NULL;

  yco = build_from_dclusterbasis_dclusteroperator(cb);
  yro = build_from_dclusterbasis_dclusteroperator(rb);

  /* weights */
  recomp_dweights(broot, rb, cb, oro, oco, yro, yco, bem);

  tcb = duplicate_dclusterbasis(cb);
  trb = duplicate_dclusterbasis(rb);

  /* truncate */
  truncate_dclusterbasis(trb, rb, yro, bro, bem->aprx->tm,
			 bem->aprx->accur_recomp);
  del_dclusterbasis(trb);
  del_dclusteroperator(yro);

  merge_dclusteropertator(oro, bro);
  del_dclusteroperator(oro);
  par->ron = enumerate_dclusteroperator(rb->t, bro);

  truncate_dclusterbasis(tcb, cb, yco, bco, bem->aprx->tm,
			 bem->aprx->accur_recomp);
  del_dclusterbasis(tcb);
  del_dclusteroperator(yco);

  merge_dclusteropertator(oco, bco);
  del_dclusteroperator(oco);
  par->con = enumerate_dclusteroperator(cb->t, bco);

}

static void
assemble_bem3d_block_dh2matrix(pdh2matrix G, uint bname, uint rname,
			       uint cname, uint pardepth, void *data)
{
  pbem3d    bem = (pbem3d) data;

  pduniform u;

  (void) pardepth;

  if (G->f) {
    bem->nearfield(G->rb->t->idx, G->cb->t->idx, bem, false, G->f);
  }
  else if (G->u) {
    u = G->u;

    if (u->rb->t->directions > 0) {
      assert(u->cb->t->directions > 0);
      bem->farfield_wave_u(rname, cname, bname, bem);
    }
    else {
      bem->farfield_u(rname, cname, bname, bem);
    }
  }
}

static void
assemble_bem3d_farfield_block_dh2matrix(pdh2matrix G, uint bname,
					uint rname, uint cname, uint pardepth,
					void *data)
{
  pbem3d    bem = (pbem3d) data;

  pduniform u;

  (void) pardepth;

  if (G->u) {
    u = G->u;

    if (u->rb->t->directions > 0) {
      assert(u->cb->t->directions > 0);
      bem->farfield_wave_u(rname, cname, bname, bem);
    }
    else {
      bem->farfield_u(rname, cname, bname, bem);
    }
  }
}

static void
assemble_bem3d_nearfield_block_dh2matrix(pdh2matrix G, uint bname,
					 uint rname, uint cname,
					 uint pardepth, void *data)
{
  pbem3d    bem = (pbem3d) data;

  (void) pardepth;
  (void) bname;
  (void) rname;
  (void) cname;

  if (G->f) {
    bem->nearfield(G->rb->t->idx, G->cb->t->idx, bem, false, G->f);
  }
}

void
assemble_bem3d_dh2matrix(pbem3d bem, pdh2matrix G)
{
  pparbem3d par = bem->par;
  par->dh2n = enumerate_dh2matrix(G);

  iterate_dh2matrix(G, 0, 0, 0, max_pardepth, assemble_bem3d_block_dh2matrix,
		    NULL, (void *) bem);

  freemem(par->dh2n);
  par->dh2n = NULL;
}

void
assemble_bem3d_farfield_dh2matrix(pbem3d bem, pdh2matrix G)
{
  pparbem3d par = bem->par;
  par->dh2n = enumerate_dh2matrix(G);

  iterate_dh2matrix(G, 0, 0, 0, max_pardepth,
		    assemble_bem3d_farfield_block_dh2matrix,
		    NULL, (void *) bem);

  freemem(par->dh2n);
  par->dh2n = NULL;
}

void
assemble_bem3d_nearfield_dh2matrix(pbem3d bem, pdh2matrix G)
{
  pparbem3d par = bem->par;
  par->dh2n = enumerate_dh2matrix(G);

  iterate_dh2matrix(G, 0, 0, 0, max_pardepth,
		    assemble_bem3d_nearfield_block_dh2matrix,
		    NULL, (void *) bem);

  freemem(par->dh2n);
  par->dh2n = NULL;
}

/* ------------------------------------------------------------
 * Surface curl
 * ------------------------------------------------------------ */

void
build_bem3d_curl_sparsematrix(pcbem3d bm, psparsematrix * C0,
			      psparsematrix * C1, psparsematrix * C2)
{
  pcsurface3d gr = bm->gr;
  const     real(*x)[3] = (const real(*)[3]) gr->x;
  const     uint(*t)[3] = (const uint(*)[3]) gr->t;
  const real *g = (const real *) gr->g;
  uint      vertices = gr->vertices;
  uint      triangles = gr->triangles;
  psparsepattern sp;
  psparsematrix lC0, lC1, lC2;
  uint      i, j;

  sp = new_sparsepattern(triangles, vertices);

  for (i = 0; i < triangles; i++)
    for (j = 0; j < 3; j++)
      addnz_sparsepattern(sp, i, t[i][j]);

  *C0 = lC0 = new_zero_sparsematrix(sp);
  *C1 = lC1 = new_zero_sparsematrix(sp);
  *C2 = lC2 = new_zero_sparsematrix(sp);

  del_sparsepattern(sp);

  for (i = 0; i < triangles; i++) {
    for (j = 0; j < 3; j++) {
      addentry_sparsematrix(lC0, i, t[i][j],
			    (x[t[i][(j + 2) % 3]][0] -
			     x[t[i][(j + 1) % 3]][0]) / g[i]);
      addentry_sparsematrix(lC1, i, t[i][j],
			    (x[t[i][(j + 2) % 3]][1] -
			     x[t[i][(j + 1) % 3]][1]) / g[i]);
      addentry_sparsematrix(lC2, i, t[i][j],
			    (x[t[i][(j + 2) % 3]][2] -
			     x[t[i][(j + 1) % 3]][2]) / g[i]);
    }
  }
}
