/* ------------------------------------------------------------
 This is the file "bem2d.c" of the H2Lib package.
 All rights reserved, Sven Christophersen 2011
 ------------------------------------------------------------ */

/* C STD LIBRARY */
/* CORE 0 */
#include "basic.h"
/* CORE 1 */
/* CORE 2 */
/* CORE 3 */
/* SIMPLE */
/* PARTICLES */
/* BEM */
#include "bem2d.h"

/* ------------------------------------------------------------
 Structs and typedefs
 ------------------------------------------------------------ */

/*
 * This constant defines the minimal width of an interval, that is used for
 * interpolation. If an interval @f$[a,b]@f$ is smaller than \ref INTERPOLATION_EPS_BEM2D
 * it is expanded to @f$[a - 0.5 \cdot \texttt{INTERPOLATION\_EPS\_BEM2D},
 * b + 0.5 \cdot \texttt{INTERPOLATION\_EPS\_BEM2D} ]@f$ .
 */
#define INTERPOLATION_EPS_BEM2D 5.0e-3

/*
 * Just an abbreviation for the struct _greencluster2d .
 */
typedef struct _greencluster2d greencluster2d;
/*
 * Pointer to a @ref greencluster2d object.
 */
typedef greencluster2d *pgreencluster2d;
/*
 * Pointer to a constant @ref greencluster2d object.
 */
typedef const greencluster2d *pcgreencluster2d;

/*
 * Just an abbreviation for the struct _greencluster2d .
 */
typedef struct _greenclusterbasis2d greenclusterbasis2d;
/*
 * Pointer to a @ref greenclusterbasis2d object.
 */
typedef greenclusterbasis2d *pgreenclusterbasis2d;
/*
 * Pointer to a constant @ref greenclusterbasis2d object.
 */
typedef const greenclusterbasis2d *pcgreenclusterbasis2d;

/*
 * @brief Substructure used for approximating @ref _hmatrix "h-", @ref
 * _uniformhmatrix "uniformh-" and @ref _h2matrix "h2matrices".
 *
 * This struct contains many needed parameters for the different approximation
 * techniques such as interpolation, green based methods or multipole.
 */
struct _aprxbem2d {

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
   * the rank can vary. In case of rectangle parameterization @ref
   * build_bem2d_rect_quadpoints the rank applies to
   * @f$ k = 8 \, m \cdot l @f$ . When using circle parameterization
   * @ref build_bem2d_circle_quadpoints @f$ k @f$ can
   * be reduced to @f$ k = 2 \, m \cdot l @f$ .
   */
  uint      k_green;

  /*
   * @brief @f$ \delta @f$ controls the distance between the bounding box of the current
   * cluster and the parameterization.
   *
   * This parameter is used as a relative
   * value. In fact the distance will be computed as @f$ \delta = \widetilde
   * \delta \cdot \operatorname{diam\_max\_cluser(t)} @f$ if rectangle parameterization
   * @ref build_bem2d_rect_quadpoints ist used. Otherwise while using circle parameterization
   * @ref build_bem2d_circle_quadpoints the distance will be computed as @f$ \delta =
   * \widetilde \delta \cdot \operatorname{diam\_2\_cluster(t)} @f$ .
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
   * @brief This is a callback function of type @ref quadpoints2d.
   *
   * It defines which
   * parameterization should be used for green based approximation techniques.
   */
  quadpoints2d quadpoints;

  /*
   * @brief Additional Information for greenhybrid methods.
   *
   * When using the greenhybrid methods, one needs to save the pivot elements
   * and the matrix @f$ V_t @f$ for each cluster. There we clone the structure of a
   * clustertree in a new struct _greencluster2d "greencluster2d" and add
   * these information to the clusters. <tt>grc_green</tt> is responsible for
   * the row clustertree.
   */
  pgreencluster2d grc_green;

  /*
   * @brief Additional Information for greenhybrid methods.
   *
   * When using the greenhybrid methods, one needs to save the pivot elements
   * and the matrix @f$ W_t @f$ for each cluster. There we clone the structure of a
   * clustertree in a new struct _greencluster2d "greencluster2d" and add
   * these information to the clusters. <tt>gcc_green</tt> is responsible for
   * the column clustertree.
   */
  pgreencluster2d gcc_green;

  /*
   * @brief Additional Information for greenhybrid methods.
   *
   * When using the greenhybrid methods, one needs to save the pivot elements
   * and the matrix @f$ V_t @f$ for each cluster. There we clone the structure of a
   * clusterbasis in a new struct _greenclusterbasis2d "greencluster2d" and add
   * these information to the clusters. <tt>grb_green</tt> is responsible for
   * the row clusterbasis. The matrix @f$ V_t @f$ is stored in the usual way inside
   * the clusterbasis itself.
   */
  pgreenclusterbasis2d grb_green;

  /*
   * @brief Additional Information for greenhybrid methods.
   *
   * When using the greenhybrid methods, one needs to save the pivot elements
   * and the matrix @f$ W_t @f$ for each cluster. There we clone the structure of a
   * clusterbasis in a new struct _greenclusterbasis2d "greencluster2d" and add
   * these information to the clusters. <tt>gcb_green</tt> is responsible for
   * the column clusterbasis. The matrix @f$ W_t @f$ is stored in the usual way inside
   * the clusterbasis itself.
   */
  pgreenclusterbasis2d gcb_green;

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

struct _parbem2d {
  /*
   * special members for different H- and H2-matrix approximations
   */

  phmatrix *hn;			/* temporary enumerated list of hmatrices */
  ph2matrix *h2n;		/* temporary enumerated list of h2matrices */
  pclusteroperator *rwn;	/* temporary enumerated list of clusteroperators for row-cluster */
  pclusteroperator *cwn;	/* temporary enumerated list of clusteroperators for col-cluster */
  uint     *leveln;		/* temporary list of levelnumber for each block in blocktree. */
  pgreencluster2d *grcn;
  uint      grcnn;
  pgreencluster2d *gccn;
  uint      gccnn;
  pgreenclusterbasis2d *grbn;
  uint      grbnn;
  pgreenclusterbasis2d *gcbn;
  uint      gcbnn;
};

struct _greencluster2d {
  uint     *xi;
  /** local indices of pivot elements */
  uint     *xihat;
  /** global indices of pivot elements */
  pamatrix  V;
  /** */
  pccluster t; /** corresponding cluster */
  uint      sons;
/** number of sons for current cluster */
};

struct _greenclusterbasis2d {
  uint     *xi;
  /** local indices of pivot elements */
  uint     *xihat;
  /** global indices of pivot elements */
  pamatrix  Qinv;
  /** Triangular factor of QR-decomposition */
  pcclusterbasis cb; /** corresponding clusterbasis */
  uint      sons;
/** number of sons for current clusterbasis */
};

/* ------------------------------------------------------------
 Constructors and destructors
 ------------------------------------------------------------ */

static void
uninit_interpolation_bem2d(paprxbem2d aprx)
{
  if (aprx->x_inter != NULL) {
    freemem(aprx->x_inter);
    aprx->x_inter = NULL;
  }
  aprx->m_inter = 0;
  aprx->k_inter = 0;
}

static void
uninit_green_bem2d(paprxbem2d aprx)
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
uninit_aca_bem2d(paprxbem2d aprx)
{
  aprx->accur_aca = 0.0;
}

static void
uninit_recompression_bem2d(paprxbem2d aprx)
{
  aprx->recomp = false;
  aprx->accur_recomp = 0.0;
  aprx->coarsen = false;
  aprx->accur_coarsen = 0.0;
  aprx->hiercomp = false;
  aprx->accur_hiercomp = 0.0;
  if (aprx->tm != NULL) {
    del_truncmode(aprx->tm);
    aprx->tm = NULL;
  }
}

static    pgreencluster2d
new_greencluster2d(pccluster c)
{
  pgreencluster2d gc;
  uint      sons = c->sons;

  gc = (pgreencluster2d) allocmem(sizeof(greencluster2d));

  gc->V = new_amatrix(0, 0);
  gc->xi = NULL;
  gc->xihat = NULL;

  gc->t = c;
  gc->sons = sons;

  return gc;
}

static void
del_greencluster2d(pgreencluster2d gc)
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

static    pgreenclusterbasis2d
new_greenclusterbasis2d(pcclusterbasis cb)
{
  pgreenclusterbasis2d gcb;
  uint      sons = cb->sons;

  gcb = (pgreenclusterbasis2d) allocmem(sizeof(greenclusterbasis2d));

  gcb->xi = NULL;
  gcb->xihat = NULL;
  gcb->Qinv = NULL;

  gcb->cb = cb;
  gcb->sons = sons;

  return gcb;
}

static void
del_greenclusterbasis2d(pgreenclusterbasis2d gcb)
{
  if (gcb->xi != NULL) {
    freemem(gcb->xi);
  }
  if (gcb->xihat != NULL) {
    freemem(gcb->xihat);
  }
  if (gcb->Qinv != NULL) {
    del_amatrix(gcb->Qinv);
  }

  freemem(gcb);
}

static    paprxbem2d
new_aprxbem2d()
{
  paprxbem2d aprx;

  aprx = (paprxbem2d) allocmem(sizeof(aprxbem2d));

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
del_aprxbem2d(paprxbem2d aprx)
{
  uninit_interpolation_bem2d(aprx);
  uninit_green_bem2d(aprx);
  uninit_aca_bem2d(aprx);
  uninit_recompression_bem2d(aprx);

  freemem(aprx);
}

static    pkernelbem2d
new_kernelbem2d()
{
  pkernelbem2d kernels;

  kernels = allocmem(sizeof(kernelbem2d));

  kernels->kernel_row = NULL;
  kernels->kernel_col = NULL;
  kernels->dnz_kernel_row = NULL;
  kernels->dnz_kernel_col = NULL;
  kernels->fundamental_row = NULL;
  kernels->fundamental_col = NULL;
  kernels->dnz_fundamental_row = NULL;
  kernels->dnz_fundamental_col = NULL;
  kernels->lagrange_row = NULL;
  kernels->lagrange_col = NULL;

  return kernels;
}

static void
del_kernelbem2d(pkernelbem2d kernels)
{
  freemem(kernels);
}

static    pparbem2d
new_parbem2d()
{
  pparbem2d par;

  par = (pparbem2d) allocmem(sizeof(parbem2d));

  par->hn = NULL;
  par->h2n = NULL;
  par->rwn = NULL;
  par->cwn = NULL;
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
del_parbem2d(pparbem2d par)
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
	del_greencluster2d(par->grcn[i]);
      }
    }
    freemem(par->grcn);
  }

  n = par->gccnn;
  if (par->gccn != NULL && n != 0) {
    for (i = 0; i < n; ++i) {
      if (par->gccn[i] != NULL) {
	del_greencluster2d(par->gccn[i]);
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
	del_greenclusterbasis2d(par->grbn[i]);
      }
    }
    freemem(par->grbn);
  }

  n = par->gcbnn;
  if (par->gcbn != NULL && n != 0) {
    for (i = 0; i < n; ++i) {
      if (par->gcbn[i] != NULL) {
	del_greenclusterbasis2d(par->gcbn[i]);
      }
    }
    freemem(par->gcbn);
  }

  freemem(par);
}

pbem2d
new_bem2d(pccurve2d gr)
{
  pbem2d    bem;

  bem = (pbem2d) allocmem(sizeof(bem2d));

  bem->gr = gr;

  bem->mass = NULL;
  bem->alpha = 0.0;

  bem->N_neumann = 0;
  bem->basis_neumann = BASIS_NONE_BEM2D;
  bem->N_dirichlet = 0;
  bem->basis_dirichlet = BASIS_NONE_BEM2D;

  bem->nearfield = NULL;
  bem->farfield_rk = NULL;
  bem->farfield_u = NULL;
  bem->leaf_row = NULL;
  bem->leaf_col = NULL;
  bem->transfer_row = NULL;
  bem->transfer_col = NULL;

  bem->aprx = new_aprxbem2d();
  bem->kernels = new_kernelbem2d();
  bem->par = new_parbem2d();

  return bem;
}

void
del_bem2d(pbem2d bem)
{

  if (bem->mass != NULL) {
    freemem(bem->mass);
  }

  if (bem->sq != NULL) {
    del_singquad1d(bem->sq);
  }

  del_aprxbem2d(bem->aprx);
  del_kernelbem2d(bem->kernels);
  del_parbem2d(bem->par);

  freemem(bem);

}

/* ------------------------------------------------------------
 Methods to build clustertrees
 ------------------------------------------------------------ */

pclustergeometry
build_bem2d_const_clustergeometry(pcbem2d bem, uint ** idx)
{
  pccurve2d gr = bem->gr;
  const     real(*x)[2] = (const real(*)[2]) gr->x;
  const     uint(*e)[2] = (const uint(*)[2]) gr->e;
  uint      edges = gr->edges;

  pclustergeometry cg;
  uint      i;

  cg = new_clustergeometry(2, edges);
  *idx = (uint *) allocmem(edges * sizeof(uint));

  for (i = 0; i < edges; i++) {
    (*idx)[i] = i;

    /* Center of gravity as characteristic point */
    cg->x[i][0] = (x[e[i][0]][0] + x[e[i][1]][0]) * 0.5;
    cg->x[i][1] = (x[e[i][0]][1] + x[e[i][1]][1]) * 0.5;

    /* Lower front left corner of bounding box */
    cg->smin[i][0] = REAL_MIN(x[e[i][0]][0], x[e[i][1]][0]);
    cg->smin[i][1] = REAL_MIN(x[e[i][0]][1], x[e[i][1]][1]);

    /* Upper back right corner of bounding box */
    cg->smax[i][0] = REAL_MAX(x[e[i][0]][0], x[e[i][1]][0]);
    cg->smax[i][1] = REAL_MAX(x[e[i][0]][1], x[e[i][1]][1]);
  }
  return cg;
}

pclustergeometry
build_bem2d_linear_clustergeometry(pcbem2d bem, uint ** idx)
{
  pccurve2d gr = bem->gr;
  const     real(*x)[2] = (const real(*)[2]) gr->x;
  const     uint(*e)[2] = (const uint(*)[2]) gr->e;
  uint      edges = gr->edges;
  uint      vertices = gr->vertices;

  pclustergeometry cg;
  uint      tv0, tv1, i, j;
  real      tmin[2], tmax[2];

  cg = new_clustergeometry(2, vertices);
  *idx = (uint *) allocmem(vertices * sizeof(uint));

  for (i = 0; i < vertices; ++i) {
    (*idx)[i] = i;

    /* Vertices as characteristic points */
    cg->x[i][0] = x[i][0];
    cg->x[i][1] = x[i][1];
    cg->smin[i][0] = x[i][0];
    cg->smin[i][1] = x[i][1];
    cg->smax[i][0] = x[i][0];
    cg->smax[i][1] = x[i][1];
  }

  for (i = 0; i < edges; i++) {
    tv0 = e[i][0];
    tv1 = e[i][1];

    /* Lower front left corner of bounding box for current triangle */
    tmin[0] = REAL_MIN(x[tv0][0], x[tv1][0]);
    tmin[1] = REAL_MIN(x[tv0][1], x[tv1][1]);

    /* Upper back right corner of bounding box for current triangle */
    tmax[0] = REAL_MAX(x[tv0][0], x[tv1][0]);
    tmax[1] = REAL_MAX(x[tv0][1], x[tv1][1]);

    /* update the bounding box for every vertex of the current triangle */
    for (j = 0; j < 2; ++j) {
      tv0 = e[i][j];

      /* Lower front left corner of bounding box */
      cg->smin[tv0][0] = REAL_MIN(cg->smin[tv0][0], tmin[0]);
      cg->smin[tv0][1] = REAL_MIN(cg->smin[tv0][1], tmin[1]);

      /* Upper back right corner of bounding box */
      cg->smax[tv0][0] = REAL_MAX(cg->smax[tv0][0], tmax[0]);
      cg->smax[tv0][1] = REAL_MAX(cg->smax[tv0][1], tmax[1]);
    }
  }

  return cg;
}

pclustergeometry
build_bem2d_clustergeometry(pcbem2d bem, uint ** idx,
			    basisfunctionbem2d basis)
{
  pclustergeometry cg;

  if (basis == BASIS_CONSTANT_BEM2D) {
    cg = build_bem2d_const_clustergeometry(bem, idx);
  }
  else {
    assert(basis == BASIS_LINEAR_BEM2D);
    cg = build_bem2d_linear_clustergeometry(bem, idx);
  }

  return cg;
}

pcluster
build_bem2d_cluster(pcbem2d bem, uint clf, basisfunctionbem2d basis)
{
  pclustergeometry cg;
  pcluster  c;
  uint     *idx;
  uint      n;

  if (basis == BASIS_CONSTANT_BEM2D) {
    cg = build_bem2d_const_clustergeometry(bem, &idx);
    n = bem->gr->edges;
  }
  else {
    assert(basis == BASIS_LINEAR_BEM2D);
    cg = build_bem2d_linear_clustergeometry(bem, &idx);
    n = bem->gr->vertices;
  }

  c = build_adaptive_cluster(cg, n, idx, clf);

  del_clustergeometry(cg);

  return c;
}

/* ------------------------------------------------------------
 Init basis approximation schemes
 ------------------------------------------------------------ */

static void
setup_interpolation_bem2d(paprxbem2d aprx, uint m)
{
  real      e;
  uint      i;

  uninit_interpolation_bem2d(aprx);

  aprx->x_inter = allocreal(m);
  aprx->m_inter = m;
  aprx->k_inter = m * m;

  /* build tschebyscheff-points */
  e = 1.0 / (2.0 * m);

  for (i = 0; i < m; ++i) {
    aprx->x_inter[i] = cos(M_PI * (2.0 * i * e + e));
  }

}

static void
setup_green_bem2d(paprxbem2d aprx, uint m, uint l, real delta,
		  quadpoints2d quadpoints)
{
  uint      i, j;
  real      h, c;
  real     *s, *ht, *hw;

  uninit_green_bem2d(aprx);

  s = allocreal(l + 1);
  ht = allocreal(m);
  hw = allocreal(m);

  aprx->m_green = m;
  aprx->l_green = l;
  aprx->delta_green = delta;
  aprx->ml_green = m * l;
  aprx->quadpoints = quadpoints;

  /* quadrature points and weights */
  if (quadpoints == build_bem2d_rect_quadpoints) {
    aprx->k_green = 8 * m * l;

    if (m > 1) {
      assemble_gauss(m, ht, hw);
    }
    else {
      assert(m == 1);
      ht[0] = 0.0;
      hw[0] = 2.0;
    }
  }
  else {
    assert(quadpoints == build_bem2d_circle_quadpoints);

    aprx->k_green = 2 * m * l;
    if (m > 1) {
      for (i = 0; i < m; ++i) {
	ht[i] = -1.0 + (2.0 * i / m);
	hw[i] = 2.0 / m;
      }
    }
    else {
      assert(m == 1);
      ht[0] = 0.0;
      hw[0] = 2.0;
    }
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

  aprx->grc_green = NULL;
  aprx->gcc_green = NULL;
  aprx->grb_green = NULL;
  aprx->gcb_green = NULL;

  freemem(s);
  freemem(ht);
  freemem(hw);
}

static void
setup_aca_bem2d(paprxbem2d aprx, real accur)
{
  assert(accur > 0.0);

  uninit_aca_bem2d(aprx);

  aprx->accur_aca = accur;
}

/* ------------------------------------------------------------
 Parameterizations for green's method
 ------------------------------------------------------------ */

void
build_bem2d_rect_quadpoints(pcbem2d bem, const real a[2], const real b[2],
			    const real delta, real(**Z)[2], real(**N)[2])
{
  paprxbem2d aprx = bem->aprx;
  const real *t = aprx->t_green;
  const real *w = aprx->w_green;
  const uint ml = aprx->ml_green;
  const uint k2 = 4 * ml;

  uint      mu, nu;
  real      velo, factor, u, c, d[2];

  *N = (real(*)[2]) allocreal(2 * k2);
  *Z = (real(*)[2]) allocreal(2 * k2);

  nu = 0;

  velo = (b[0] - a[0] + 2.0 * delta) * 0.5;
  c = 0.5 * (b[0] + a[0]);
  d[0] = a[1] - delta;
  d[1] = b[1] + delta;

  for (mu = 0; mu < ml; mu++) {
    u = velo * t[mu];
    factor = velo * w[mu];

    (*Z)[nu][0] = c + u;
    (*Z)[nu][1] = d[0];
    (*N)[nu][0] = 0.0;
    (*N)[nu][1] = -factor;
    nu++;

    (*Z)[nu][0] = c + u;
    (*Z)[nu][1] = d[1];
    (*N)[nu][0] = 0.0;
    (*N)[nu][1] = factor;
    nu++;

  }

  velo = (b[1] - a[1] + 2.0 * delta) * 0.5;
  c = 0.5 * (b[1] + a[1]);
  d[0] = b[0] + delta;
  d[1] = a[0] - delta;

  for (mu = 0; mu < ml; mu++) {
    u = velo * t[mu];
    factor = velo * w[mu];

    (*Z)[nu][0] = d[0];
    (*Z)[nu][1] = c + u;
    (*N)[nu][0] = factor;
    (*N)[nu][1] = 0.0;
    nu++;

    (*N)[nu][0] = 0.0;
    (*N)[nu][1] = 0.0;

    (*Z)[nu][0] = d[1];
    (*Z)[nu][1] = c + u;
    (*N)[nu][0] = -factor;
    (*N)[nu][1] = 0.0;
    nu++;

  }

  assert(nu == k2);
}

void
build_bem2d_circle_quadpoints(pcbem2d bem, const real a[2],
			      const real b[2], const real delta, real(**Z)[2],
			      real(**N)[2])
{
  paprxbem2d aprx = bem->aprx;
  const real *t = aprx->t_green;
  const real *w = aprx->w_green;
  const uint ml = aprx->ml_green;

  real      center[2];
  real      R, factor, c, s;
  uint      mu;

  *N = (real(*)[2]) allocreal(2 * ml);
  *Z = (real(*)[2]) allocreal(2 * ml);

  center[0] = 0.5 * (b[0] + a[0]);
  center[1] = 0.5 * (b[1] + a[1]);
  R = 0.5 * REAL_SQRT(REAL_SQR(b[0] - a[0]) + REAL_SQR(b[1] - a[1])) + delta;

  for (mu = 0; mu < ml; mu++) {
    factor = R * M_PI * w[mu];
    c = cos(M_PI * t[mu]);
    s = sin(M_PI * t[mu]);

    (*Z)[mu][0] = center[0] + R * c;
    (*Z)[mu][1] = center[1] + R * s;
    (*N)[mu][0] = factor * c;
    (*N)[mu][1] = factor * s;

  }
}

/* ------------------------------------------------------------
 Transfer Interpolation points to bounding box
 ------------------------------------------------------------ */

static void
assemble_interpoints2d_array(pcbem2d bem, pccluster t, real(*X)[2])
{
  pcaprxbem2d aprx = bem->aprx;
  real     *x = aprx->x_inter;
  uint      m = aprx->m_inter;
  uint      k = aprx->k_inter;
  real      ax = t->bmin[0];
  real      bx = t->bmax[0];
  real      ay = t->bmin[1];
  real      by = t->bmax[1];

  real      cx, dx, cy, dy;
  uint      i, j, index;

  if (bx - ax < INTERPOLATION_EPS_BEM2D) {
    bx += INTERPOLATION_EPS_BEM2D;
    ax -= INTERPOLATION_EPS_BEM2D;
  }
  if (by - ay < INTERPOLATION_EPS_BEM2D) {
    by += INTERPOLATION_EPS_BEM2D;
    ay -= INTERPOLATION_EPS_BEM2D;
  }

  /*
   * Tschebyscheff points
   */

  cx = (bx + ax) * 0.5;
  dx = (bx - ax) * 0.5;
  cy = (by + ay) * 0.5;
  dy = (by - ay) * 0.5;

  index = 0;
  for (i = 0; i < m; ++i) {
    for (j = 0; j < m; ++j) {
      X[index][0] = cx + dx * x[i];
      X[index][1] = cy + dy * x[j];
      index++;
    }
  }

  assert(index == k);
}

static void
assemble_interpoints2d_avector(pcbem2d bem, pccluster t,
			       pavector px, pavector py)
{
  pcaprxbem2d aprx = bem->aprx;
  real     *x = aprx->x_inter;
  uint      m = aprx->m_inter;
  real      ax = t->bmin[0];
  real      bx = t->bmax[0];
  real      ay = t->bmin[1];
  real      by = t->bmax[1];

  real      cx, dx, cy, dy;
  uint      i;

  if (bx - ax < INTERPOLATION_EPS_BEM2D) {
    bx += INTERPOLATION_EPS_BEM2D;
    ax -= INTERPOLATION_EPS_BEM2D;
  }
  if (by - ay < INTERPOLATION_EPS_BEM2D) {
    by += INTERPOLATION_EPS_BEM2D;
    ay -= INTERPOLATION_EPS_BEM2D;
  }

  /*
   * Tschebyscheff points
   */

  cx = (bx + ax) * 0.5;
  dx = (bx - ax) * 0.5;
  cy = (by + ay) * 0.5;
  dy = (by - ay) * 0.5;

  for (i = 0; i < m; ++i) {
    px->v[i] = cx + dx * x[i];
    py->v[i] = cy + dy * x[i];

  }
}

/* ------------------------------------------------------------
 Methods to fill rank-k-matrices with several techniques
 ------------------------------------------------------------ */

static void
assemble_bem2d_inter_row_rkmatrix(pccluster rc, uint rname,
				  pccluster cc, uint cname, pcbem2d bem,
				  prkmatrix R)
{
  paprxbem2d aprx = bem->aprx;
  pkernelbem2d kernels = bem->kernels;
  const uint rows = rc->size;
  const uint cols = cc->size;
  const uint m = aprx->m_inter;
  const uint k = aprx->k_inter;
  real(*z)[2];
  pavector  px, py;

  (void) rname;
  (void) cname;

  z = (real(*)[2]) allocreal((size_t) (2 * k));
  px = new_avector(m);
  py = new_avector(m);

  assemble_interpoints2d_array(bem, rc, z);
  assemble_interpoints2d_avector(bem, rc, px, py);

  resize_rkmatrix(R, rows, cols, k);

  kernels->lagrange_row(rc->idx, px, py, bem, &R->A);
  kernels->kernel_col(cc->idx, (const real(*)[2]) z, bem, &R->B);

  del_avector(px);
  del_avector(py);
  freemem(z);
}

static void
assemble_bem2d_inter_col_rkmatrix(pccluster rc, uint rname,
				  pccluster cc, uint cname, pcbem2d bem,
				  prkmatrix R)
{
  paprxbem2d aprx = bem->aprx;
  pkernelbem2d kernels = bem->kernels;
  const uint rows = rc->size;
  const uint cols = cc->size;
  const uint m = aprx->m_inter;
  const uint k = aprx->k_inter;
  real(*z)[2];
  pavector  px, py;

  (void) rname;
  (void) cname;

  z = (real(*)[2]) allocreal((size_t) (2 * k));
  px = new_avector(m);
  py = new_avector(m);

  assemble_interpoints2d_array(bem, cc, z);
  assemble_interpoints2d_avector(bem, cc, px, py);

  resize_rkmatrix(R, rows, cols, k);

  kernels->fundamental_row(rc->idx, (const real(*)[2]) z, bem, &R->A);
  kernels->lagrange_col(cc->idx, px, py, bem, &R->B);

  del_avector(px);
  del_avector(py);
  freemem(z);
}

static void
assemble_bem2d_inter_mixed_rkmatrix(pccluster rc, uint rname,
				    pccluster cc, uint cname, pcbem2d bem,
				    prkmatrix R)
{
  if (getdiam_2_cluster(rc) < getdiam_2_cluster(cc)) {
    assemble_bem2d_inter_row_rkmatrix(rc, rname, cc, cname, bem, R);
  }
  else {
    assemble_bem2d_inter_col_rkmatrix(rc, rname, cc, cname, bem, R);
  }
}

static void
assemble_bem2d_green_row_rkmatrix(pccluster rc, uint rname,
				  pccluster cc, uint cname, pcbem2d bem,
				  prkmatrix R)
{
  paprxbem2d aprx = bem->aprx;
  pkernelbem2d kernels = bem->kernels;
  pamatrix  A = &R->A;
  pamatrix  B = &R->B;
  real      delta = aprx->delta_green;
  real     *a = rc->bmin;
  real     *b = rc->bmax;

  pamatrix  T;
  real(*Z)[2], (*N)[2];
  real      diam;
  uint      nu, k2;

  (void) rname;
  (void) cname;

  k2 = 0.5 * aprx->k_green;
  diam =
    (aprx->quadpoints == build_bem2d_rect_quadpoints) ?
    getdiam_max_cluster(rc) : getdiam_2_cluster(rc);
  delta = delta * diam;

  T = new_amatrix(0, 0);

  resize_rkmatrix(R, rc->size, cc->size, 2 * k2);

  aprx->quadpoints(bem, a, b, delta, &Z, &N);

  init_sub_amatrix(T, A, rc->size, 0, k2, 0);
  kernels->fundamental_row(rc->idx, (const real(*)[2]) Z, (pcbem2d) bem, T);
  uninit_amatrix(T);

  init_sub_amatrix(T, B, cc->size, 0, k2, 0);
  kernels->dnz_kernel_col(cc->idx, (const real(*)[2]) Z, (const real(*)[2]) N,
			  (pcbem2d) bem, T);
  uninit_amatrix(T);

  for (nu = 0; nu < k2; ++nu) {
    N[nu][0] *= -1.0;
    N[nu][1] *= -1.0;
  }

  init_sub_amatrix(T, A, rc->size, 0, k2, k2);
  kernels->dnz_fundamental_row(rc->idx, (const real(*)[2]) Z,
			       (const real(*)[2]) N, (pcbem2d) bem, T);
  uninit_amatrix(T);

  init_sub_amatrix(T, B, cc->size, 0, k2, k2);
  kernels->kernel_col(cc->idx, (const real(*)[2]) Z, (pcbem2d) bem, T);
  uninit_amatrix(T);

  del_amatrix(T);
  freemem(Z);
  freemem(N);

}

static void
assemble_bem2d_green_col_rkmatrix(pccluster rc, uint rname,
				  pccluster cc, uint cname, pcbem2d bem,
				  prkmatrix R)
{
  paprxbem2d aprx = bem->aprx;
  pkernelbem2d kernels = bem->kernels;
  pamatrix  A = &R->A;
  pamatrix  B = &R->B;
  real      delta = aprx->delta_green;
  uint      k = aprx->k_green;
  real     *a = cc->bmin;
  real     *b = cc->bmax;

  pamatrix  T;
  real(*Z)[2], (*N)[2];
  real      diam;
  uint      k2, nu;

  (void) rname;
  (void) cname;

  k2 = 0.5 * k;
  diam =
    (aprx->quadpoints == build_bem2d_rect_quadpoints) ?
    getdiam_max_cluster(cc) : getdiam_2_cluster(cc);
  delta = delta * diam;

  T = new_amatrix(0, 0);

  resize_rkmatrix(R, rc->size, cc->size, 2 * k2);

  aprx->quadpoints(bem, a, b, delta, &Z, &N);

  init_sub_amatrix(T, B, cc->size, 0, k2, 0);
  kernels->kernel_col(cc->idx, (const real(*)[2]) Z, (pcbem2d) bem, T);
  uninit_amatrix(T);

  init_sub_amatrix(T, A, rc->size, 0, k2, 0);
  kernels->dnz_fundamental_row(rc->idx, (const real(*)[2]) Z,
			       (const real(*)[2]) N, (pcbem2d) bem, T);
  uninit_amatrix(T);

  for (nu = 0; nu < k2; ++nu) {
    N[nu][0] *= -1.0;
    N[nu][1] *= -1.0;
  }

  init_sub_amatrix(T, B, cc->size, 0, k2, k2);
  kernels->dnz_kernel_col(cc->idx, (const real(*)[2]) Z, (const real(*)[2]) N,
			  (pcbem2d) bem, T);
  uninit_amatrix(T);

  init_sub_amatrix(T, A, rc->size, 0, k2, k2);
  kernels->fundamental_row(rc->idx, (const real(*)[2]) Z, (pcbem2d) bem, T);
  uninit_amatrix(T);

  del_amatrix(T);
  freemem(Z);
  freemem(N);

}

static void
assemble_bem2d_green_mixed_rkmatrix(pccluster rc, uint rname,
				    pccluster cc, uint cname, pcbem2d bem,
				    prkmatrix R)
{
  paprxbem2d aprx = bem->aprx;

  real      diamt, diams;

  diamt =
    (aprx->quadpoints == build_bem2d_rect_quadpoints) ?
    getdiam_max_cluster(rc) : getdiam_2_cluster(rc);
  diams =
    (aprx->quadpoints == build_bem2d_rect_quadpoints) ?
    getdiam_max_cluster(cc) : getdiam_2_cluster(cc);

  if (diamt < diams) {
    assemble_bem2d_green_row_rkmatrix(rc, rname, cc, cname, bem, R);
  }
  else {
    assemble_bem2d_green_col_rkmatrix(rc, rname, cc, cname, bem, R);
  }
}

static void
assemble_row_greencluster2d(pcbem2d bem, pgreencluster2d gc)
{
  paprxbem2d aprx = bem->aprx;
  pkernelbem2d kernels = bem->kernels;
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
  real(*Z)[2], (*N)[2];
  real      diam;
  uint     *xi;
  uint      i, k2;

  k2 = 0.5 * k;
  diam =
    (aprx->quadpoints == build_bem2d_rect_quadpoints) ?
    getdiam_max_cluster(c) : getdiam_2_cluster(c);
  delta = delta * diam;

  A_t = new_amatrix(rows, k);
  T = new_amatrix(0, 0);

  aprx->quadpoints(bem, a, b, delta, &Z, &N);

  init_sub_amatrix(T, A_t, c->size, 0, k2, 0);
  kernels->fundamental_row(c->idx, (const real(*)[2]) Z, (pcbem2d) bem, T);
  uninit_amatrix(T);

  init_sub_amatrix(T, A_t, c->size, 0, k2, k2);
  kernels->dnz_fundamental_row(c->idx, (const real(*)[2]) Z,
			       (const real(*)[2]) N, (pcbem2d) bem, T);
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
assemble_bem2d_greenhybrid_row_rkmatrix(pccluster rc, uint rname,
					pccluster cc, uint cname, pcbem2d bem,
					prkmatrix R)
{
  pparbem2d par = bem->par;
  pamatrix  A = &R->A;
  pamatrix  B = &R->B;
  uint      rows = rc->size;
  uint      cols = cc->size;

  pamatrix  V;
  pgreencluster2d grc;
  uint     *xihat;
  uint      rank;

  (void) cname;

#ifdef USE_OPENMP
#pragma omp critical(greenhybrid)
#endif
  {
    grc = par->grcn[rname];

    if (grc == NULL) {
      grc = par->grcn[rname] = new_greencluster2d(rc);
      assemble_row_greencluster2d(bem, grc);
    }
  }

  V = grc->V;
  rank = V->cols;
  xihat = grc->xihat;

  resize_rkmatrix(R, rows, cols, rank);

  copy_amatrix(false, V, A);
  bem->nearfield(xihat, cc->idx, bem, true, B);
}

static void
assemble_col_greencluster2d(pcbem2d bem, pgreencluster2d gc)
{
  paprxbem2d aprx = bem->aprx;
  pkernelbem2d kernels = bem->kernels;
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
  real(*Z)[2], (*N)[2];
  uint     *xi;
  real      diam;
  uint      i, k2;

  k2 = 0.5 * k;
  diam =
    (aprx->quadpoints == build_bem2d_rect_quadpoints) ?
    getdiam_max_cluster(c) : getdiam_2_cluster(c);
  delta = delta * diam;

  A_t = new_amatrix(rows, k);

  T = new_amatrix(0, 0);

  aprx->quadpoints(bem, a, b, delta, &Z, &N);

  init_sub_amatrix(T, A_t, c->size, 0, k2, 0);
  kernels->kernel_col(c->idx, (const real(*)[2]) Z, (pcbem2d) bem, T);
  uninit_amatrix(T);

  init_sub_amatrix(T, A_t, c->size, 0, k2, k2);
  kernels->dnz_kernel_col(c->idx, (const real(*)[2]) Z, (const real(*)[2]) N,
			  (pcbem2d) bem, T);
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
assemble_bem2d_greenhybrid_col_rkmatrix(pccluster rc, uint rname,
					pccluster cc, uint cname, pcbem2d bem,
					prkmatrix R)
{
  pparbem2d par = bem->par;
  pamatrix  A = &R->A;
  pamatrix  B = &R->B;
  uint      rows = rc->size;
  uint      cols = cc->size;

  pamatrix  V;
  pgreencluster2d gcc;
  uint     *xihat;
  uint      rank;

  (void) rname;

#ifdef USE_OPENMP
#pragma omp critical(greenhybrid)
#endif
  {
    gcc = par->gccn[cname];

    if (gcc == NULL) {
      gcc = par->gccn[cname] = new_greencluster2d(cc);
      assemble_col_greencluster2d(bem, gcc);
    }
  }

  V = gcc->V;
  rank = V->cols;
  xihat = gcc->xihat;

  resize_rkmatrix(R, rows, cols, rank);

  bem->nearfield(rc->idx, xihat, bem, false, A);
  copy_amatrix(false, V, B);
}

static void
assemble_bem2d_greenhybrid_mixed_rkmatrix(pccluster rc, uint rname,
					  pccluster cc, uint cname,
					  pcbem2d bem, prkmatrix R)
{
  pparbem2d par = bem->par;
  pamatrix  A = &R->A;
  pamatrix  B = &R->B;
  uint      rows = rc->size;
  uint      cols = cc->size;

  pamatrix  V, W;
  pgreencluster2d grc, gcc;
  uint     *xihatV, *xihatW;
  uint      rankV, rankW;

#ifdef USE_OPENMP
#pragma omp critical(greenhybrid)
#endif
  {
    grc = par->grcn[rname];

    if (grc == NULL) {
      grc = par->grcn[rname] = new_greencluster2d(rc);
      assemble_row_greencluster2d(bem, grc);
    }
  }

#ifdef USE_OPENMP
#pragma omp critical(greenhybrid)
#endif
  {
    gcc = par->gccn[cname];

    if (gcc == NULL) {
      gcc = par->gccn[cname] = new_greencluster2d(cc);
      assemble_col_greencluster2d(bem, gcc);
    }
  }

  rankV = grc->V->cols;
  rankW = gcc->V->cols;

  if (cols * rankV <= rows * rankW) {
    V = grc->V;
    xihatV = grc->xihat;

    resize_rkmatrix(R, rows, cols, rankV);

    copy_amatrix(false, V, A);
    bem->nearfield(xihatV, cc->idx, bem, true, B);
  }
  else {
    W = gcc->V;
    xihatW = gcc->xihat;

    resize_rkmatrix(R, rows, cols, rankW);

    bem->nearfield(rc->idx, xihatW, bem, false, A);
    copy_amatrix(false, W, B);
  }

}

static void
assemble_bem2d_ACA_rkmatrix(pccluster rc, uint rname, pccluster cc,
			    uint cname, pcbem2d bem, prkmatrix R)
{
  paprxbem2d aprx = bem->aprx;
  const real accur = aprx->accur_aca;
  const uint *ridx = rc->idx;
  const uint *cidx = cc->idx;
  const uint rows = rc->size;
  const uint cols = cc->size;

  pamatrix  G;

  (void) rname;
  (void) cname;

  G = new_amatrix(rows, cols);
  bem->nearfield(ridx, cidx, bem, false, G);

  decomp_fullaca_rkmatrix(G, accur, NULL, NULL, R);

  del_amatrix(G);
}

static void
assemble_bem2d_PACA_rkmatrix(pccluster rc, uint rname, pccluster cc,
			     uint cname, pcbem2d bem, prkmatrix R)
{
  paprxbem2d aprx = bem->aprx;
  const real accur = aprx->accur_aca;
  matrixentry_t entry = (matrixentry_t) bem->nearfield;
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
assemble_bem2d_HCA_rkmatrix(pccluster rc, uint rname, pccluster cc,
			    uint cname, pcbem2d bem, prkmatrix R)
{
  paprxbem2d aprx = bem->aprx;
  pkernelbem2d kernels = bem->kernels;
  const real accur = aprx->accur_aca;
  const uint k = aprx->k_inter;
  const uint *ridx = rc->idx;
  const uint *cidx = cc->idx;
  const uint rows = rc->size;
  const uint cols = cc->size;

  prkmatrix R2;
  pamatrix  S, C, D;
  real(*IT)[2], (*IS)[2], (*ITk)[2], (*ISk)[2];
  uint     *I_k, *J_k;
  uint      i, j, rank;

  (void) rname;
  (void) cname;

  IT = (real(*)[2]) allocreal(2 * k);
  IS = (real(*)[2]) allocreal(2 * k);

  assemble_interpoints2d_array(bem, rc, IT);
  assemble_interpoints2d_array(bem, cc, IS);

  S = new_amatrix(k, k);
  R2 = new_rkmatrix(k, k, 0);

  kernels->fundamental((const real(*)[2]) IT, (const real(*)[2]) IS, S);
  decomp_fullaca_rkmatrix(S, accur, &I_k, &J_k, R2);
  rank = R2->k;
  C = &R2->A;
  D = &R2->B;

  ITk = (real(*)[2]) allocreal(2 * rank);
  ISk = (real(*)[2]) allocreal(2 * rank);

  for (i = 0; i < rank; ++i) {
    for (j = 0; j < 2; ++j) {
      ITk[i][j] = IT[I_k[i]][j];
      ISk[i][j] = IS[J_k[i]][j];
    }
  }

  resize_amatrix(S, rank, rank);

  copy_lower_aca_amatrix(true, C, I_k, S);
  copy_upper_aca_amatrix(false, D, J_k, S);

  resize_rkmatrix(R, rows, cols, rank);

  kernels->kernel_row(ridx, (const real(*)[2]) ISk, bem, &R->A);
  kernels->kernel_col(cidx, (const real(*)[2]) ITk, bem, &R->B);

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

/* ------------------------------------------------------------
 Methods to fill clusterbasis and h2matrix
 ------------------------------------------------------------ */

static void
assemble_bem2d_inter_leaf_row_clusterbasis(pcbem2d bem,
					   pclusterbasis rb, uint rname)
{
  paprxbem2d aprx = bem->aprx;
  pkernelbem2d kernels = bem->kernels;
  pamatrix  V = &rb->V;
  pccluster t = rb->t;
  const uint m = aprx->m_inter;
  const uint k = aprx->k_inter;

  pavector  px, py;

  (void) rname;

  px = new_avector(m);
  py = new_avector(m);

  assemble_interpoints2d_avector(bem, t, px, py);

  resize_amatrix(V, t->size, k);
  rb->k = k;
  update_clusterbasis(rb);

  kernels->lagrange_row(t->idx, px, py, bem, V);

  del_avector(px);
  del_avector(py);
}

static void
assemble_bem2d_inter_tranfer_clusterbasis(pcbem2d bem,
					  pclusterbasis cb, uint rname)
{
  paprxbem2d aprx = bem->aprx;
  pccluster t = cb->t;
  uint      sons = t->sons;
  const uint m = aprx->m_inter;
  const uint k = aprx->k_inter;

  pamatrix  E;
  real(*X)[2];
  pavector  px, py;
  uint      s;

  (void) rname;

  X = (real(*)[2]) allocmem(2 * k * sizeof(real));
  px = new_avector(m);
  py = new_avector(m);

  assemble_interpoints2d_avector(bem, t, px, py);

  resize_clusterbasis(cb, k);

  for (s = 0; s < sons; ++s) {
    E = &cb->son[s]->E;

    assemble_interpoints2d_array(bem, cb->son[s]->t, X);

    assemble_bem2d_lagrange_amatrix((const real(*)[2]) X, px, py, E);
  }

  del_avector(px);
  del_avector(py);
  freemem(X);
}

static void
assemble_bem2d_inter_leaf_col_clusterbasis(pcbem2d bem,
					   pclusterbasis cb, uint cname)
{
  paprxbem2d aprx = bem->aprx;
  pkernelbem2d kernels = bem->kernels;
  pamatrix  V = &cb->V;
  pccluster t = cb->t;
  const uint m = aprx->m_inter;
  const uint k = aprx->k_inter;

  pavector  px, py;

  (void) cname;

  px = new_avector(m);
  py = new_avector(m);

  assemble_interpoints2d_avector(bem, t, px, py);

  resize_amatrix(V, t->size, k);
  cb->k = k;
  update_clusterbasis(cb);

  kernels->lagrange_col(t->idx, px, py, bem, V);

  del_avector(px);
  del_avector(py);
}

static void
assemble_bem2d_inter_uniform(uint rname, uint cname, pcbem2d bem, puniform U)
{
  pkernelbem2d kernels = bem->kernels;
  pccluster rc = U->rb->t;
  pccluster cc = U->cb->t;
  const uint kr = U->rb->k;
  const uint kc = U->cb->k;
  pamatrix  S = &U->S;

  real(*xi_r)[2], (*xi_c)[2];

  (void) rname;
  (void) cname;

  resize_amatrix(S, kr, kc);

  xi_r = (real(*)[2]) allocreal(2 * kr);
  xi_c = (real(*)[2]) allocreal(2 * kc);

  assemble_interpoints2d_array(bem, rc, xi_r);
  assemble_interpoints2d_array(bem, cc, xi_c);

  kernels->fundamental((const real(*)[2]) xi_r, (const real(*)[2]) xi_c, S);

  freemem(xi_r);
  freemem(xi_c);
}

static void
update_pivotelements_greenclusterbasis2d(pgreenclusterbasis2d * grbn,
					 pcclusterbasis cb, uint * I_t,
					 uint k)
{
  uint      sons = cb->sons;
  pgreenclusterbasis2d grb = *grbn;

  pgreenclusterbasis2d grb1;
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
collect_pivotelements_greenclusterbasis2d(pgreenclusterbasis2d * grbn,
					  uint * rank)
{
  pcclusterbasis cb = (*grbn)->cb;
  uint      sons = cb->sons;

  pgreenclusterbasis2d grb1;
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
assemble_bem2d_greenhybrid_leaf_row_clusterbasis(pcbem2d bem,
						 pclusterbasis rb, uint rname)
{
  paprxbem2d aprx = bem->aprx;
  pkernelbem2d kernels = bem->kernels;
  pparbem2d par = bem->par;
  pccluster c = rb->t;
  const uint rows = c->size;
  const real eps = aprx->accur_aca;
  real      delta = aprx->delta_green;
  uint      rank = aprx->k_green;

  prkmatrix R;
  pamatrix  A_t, T, RC;
  pgreenclusterbasis2d grb;
  real(*Z)[2], (*N)[2];
  uint     *xi;
  uint      k2;

  k2 = rank * 0.5;
  delta = delta * getdiam_max_cluster(c);

  aprx->quadpoints(bem, c->bmin, c->bmax, delta, &Z, &N);

  T = new_amatrix(0, 0);

  grb = par->grbn[rname];
  if (grb == NULL) {
    grb = par->grbn[rname] = new_greenclusterbasis2d(rb);
  }

  A_t = new_amatrix(rows, rank);

  init_sub_amatrix(T, A_t, rows, 0, k2, 0);
  kernels->fundamental_row(c->idx, (const real(*)[2]) Z, (pcbem2d) bem, T);
  uninit_amatrix(T);

  init_sub_amatrix(T, A_t, rows, 0, k2, k2);
  kernels->dnz_fundamental_row(c->idx, (const real(*)[2]) Z,
			       (const real(*)[2]) N, (pcbem2d) bem, T);
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

  update_pivotelements_greenclusterbasis2d(par->grbn + rname, rb, xi, rank);

  del_rkmatrix(R);
  del_amatrix(A_t);
  del_amatrix(RC);
  del_amatrix(T);
  freemem(Z);
  freemem(N);
}

static void
assemble_bem2d_greenhybrid_transfer_row_clusterbasis(pcbem2d bem,
						     pclusterbasis rb,
						     uint rname)
{
  paprxbem2d aprx = bem->aprx;
  pkernelbem2d kernels = bem->kernels;
  pparbem2d par = bem->par;
  pccluster c = rb->t;
  const real eps = aprx->accur_aca;
  real      delta = aprx->delta_green;
  uint      rank = aprx->k_green;
  uint      sons = c->sons;

  prkmatrix R;
  pamatrix  A_t, T, RC, E;
  pgreenclusterbasis2d grb;
  real(*Z)[2], (*N)[2];
  uint     *I_t, *idx;
  uint      i, k2, newrank;

  k2 = rank * 0.5;
  delta = delta * getdiam_max_cluster(c);

  aprx->quadpoints(bem, c->bmin, c->bmax, delta, &Z, &N);

  T = new_amatrix(0, 0);

  grb = par->grbn[rname];
  if (grb == NULL) {
    grb = par->grbn[rname] = new_greenclusterbasis2d(rb);
  }

  idx = collect_pivotelements_greenclusterbasis2d(par->grbn + rname, &rank);
  A_t = new_amatrix(rank, 2 * k2);

  init_sub_amatrix(T, A_t, rank, 0, k2, 0);
  kernels->fundamental_row(idx, (const real(*)[2]) Z, (pcbem2d) bem, T);
  uninit_amatrix(T);

  init_sub_amatrix(T, A_t, rank, 0, k2, k2);
  kernels->dnz_fundamental_row(idx, (const real(*)[2]) Z,
			       (const real(*)[2]) N, (pcbem2d) bem, T);
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

  update_pivotelements_greenclusterbasis2d(par->grbn + rname, rb, I_t,
					   newrank);

  freemem(idx);
  freemem(I_t);
  del_rkmatrix(R);
  del_amatrix(RC);
  del_amatrix(A_t);
  del_amatrix(T);
  freemem(Z);
  freemem(N);
}

static void
assemble_bem2d_greenhybrid_leaf_col_clusterbasis(pcbem2d bem,
						 pclusterbasis cb, uint cname)
{
  paprxbem2d aprx = bem->aprx;
  pkernelbem2d kernels = bem->kernels;
  pparbem2d par = bem->par;
  pccluster c = cb->t;
  const uint rows = c->size;
  const real eps = aprx->accur_aca;
  real      delta = aprx->delta_green;
  uint      rank = aprx->k_green;

  prkmatrix R;
  pamatrix  A_t, T, RC;
  pgreenclusterbasis2d gcb;
  real(*Z)[2], (*N)[2];
  uint     *xi;
  uint      k2;

  k2 = rank * 0.5;
  delta = delta * getdiam_max_cluster(c);

  aprx->quadpoints(bem, c->bmin, c->bmax, delta, &Z, &N);

  T = new_amatrix(0, 0);

  gcb = par->gcbn[cname];
  if (gcb == NULL) {
    gcb = par->gcbn[cname] = new_greenclusterbasis2d(cb);
  }

  A_t = new_amatrix(rows, rank);

  init_sub_amatrix(T, A_t, rows, 0, k2, 0);
  kernels->kernel_col(c->idx, (const real(*)[2]) Z, (pcbem2d) bem, T);
  uninit_amatrix(T);

  init_sub_amatrix(T, A_t, rows, 0, k2, k2);
  kernels->dnz_kernel_col(c->idx, (const real(*)[2]) Z, (const real(*)[2]) N,
			  (pcbem2d) bem, T);
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

  update_pivotelements_greenclusterbasis2d(par->gcbn + cname, cb, xi, rank);

  del_rkmatrix(R);
  del_amatrix(A_t);
  del_amatrix(RC);
  del_amatrix(T);
  freemem(Z);
  freemem(N);
}

static void
assemble_bem2d_greenhybrid_transfer_col_clusterbasis(pcbem2d bem,
						     pclusterbasis cb,
						     uint cname)
{
  paprxbem2d aprx = bem->aprx;
  pkernelbem2d kernels = bem->kernels;
  pparbem2d par = bem->par;
  pccluster c = cb->t;
  const real eps = aprx->accur_aca;
  real      delta = aprx->delta_green;
  uint      rank = aprx->k_green;
  uint      sons = c->sons;

  prkmatrix R;
  pamatrix  A_t, T, RC, E;
  pgreenclusterbasis2d gcb;
  real(*Z)[2], (*N)[2];
  uint     *I_t, *idx;
  uint      i, k2, newrank;

  k2 = rank * 0.5;
  delta = delta * getdiam_max_cluster(c);

  aprx->quadpoints(bem, c->bmin, c->bmax, delta, &Z, &N);

  T = new_amatrix(0, 0);

  gcb = par->gcbn[cname];
  if (gcb == NULL) {
    gcb = par->gcbn[cname] = new_greenclusterbasis2d(cb);
  }

  idx = collect_pivotelements_greenclusterbasis2d(par->gcbn + cname, &rank);
  A_t = new_amatrix(rank, 2 * k2);

  init_sub_amatrix(T, A_t, rank, 0, k2, 0);
  kernels->kernel_col(idx, (const real(*)[2]) Z, (pcbem2d) bem, T);
  uninit_amatrix(T);

  init_sub_amatrix(T, A_t, rank, 0, k2, k2);
  kernels->dnz_kernel_col(idx, (const real(*)[2]) Z, (const real(*)[2]) N,
			  (pcbem2d) bem, T);
  uninit_amatrix(T);

  R = new_rkmatrix(rank, 2 * k2, 0);
  decomp_fullaca_rkmatrix(A_t, eps, &I_t, NULL, R);
  newrank = R->k;

  RC = new_amatrix(newrank, newrank);
  copy_lower_aca_amatrix(true, &R->A, I_t, RC);
  triangularsolve_amatrix(true, true, true, RC, true, &R->A);

  resize_clusterbasis(cb, newrank);

  rank = 0;
  for (i = 0; i < sons; ++i) {
    E = &cb->son[i]->E;
    init_sub_amatrix(T, &R->A, cb->son[i]->k, rank, newrank, 0);
    copy_amatrix(false, T, E);
    uninit_amatrix(T);
    rank += cb->son[i]->k;
  }

  update_pivotelements_greenclusterbasis2d(par->gcbn + cname, cb, I_t,
					   newrank);

  freemem(idx);
  freemem(I_t);
  del_rkmatrix(R);
  del_amatrix(RC);
  del_amatrix(A_t);
  del_amatrix(T);
  freemem(Z);
  freemem(N);
}

static void
assemble_bem2d_greenhybrid_uniform(uint rname, uint cname,
				   pcbem2d bem, puniform U)
{
  pparbem2d par = bem->par;
  const uint kr = U->rb->k;
  const uint kc = U->cb->k;
  pamatrix  S = &U->S;

  pgreenclusterbasis2d grb, gcb;
  uint     *xihatV, *xihatW;

  grb = par->grbn[rname];
  gcb = par->gcbn[cname];

  assert(grb != NULL);
  assert(gcb != NULL);

  xihatV = grb->xihat;
  xihatW = gcb->xihat;

  resize_amatrix(S, kr, kc);

  bem->nearfield(xihatV, xihatW, bem, false, S);
}

static void
assemble_bem2d_greenhybridortho_leaf_row_clusterbasis(pcbem2d bem,
						      pclusterbasis rb,
						      uint rname)
{
  paprxbem2d aprx = bem->aprx;
  pkernelbem2d kernels = bem->kernels;
  pparbem2d par = bem->par;
  pccluster t = rb->t;
  pamatrix  V = &rb->V;
  const uint rows = t->size;
  const real eps = aprx->accur_aca;
  real      delta = aprx->delta_green;
  uint      rank = aprx->k_green;

  prkmatrix R;
  pamatrix  A_t, T, RC;
  pavector  tau;
  pgreenclusterbasis2d grb;
  real(*Z)[2], (*N)[2];
  uint     *xi;
  uint      k2;

  k2 = rank * 0.5;
  delta = aprx->delta_green
    * ((aprx->quadpoints == build_bem2d_rect_quadpoints) ?
       getdiam_max_cluster(t) : getdiam_2_cluster(t));

  aprx->quadpoints(bem, t->bmin, t->bmax, delta, &Z, &N);

  T = new_amatrix(0, 0);

  grb = par->grbn[rname];
  if (grb == NULL) {
    grb = par->grbn[rname] = new_greenclusterbasis2d(rb);
  }

  /* setup up A_t from green's formula. */
  A_t = new_amatrix(rows, rank);

  init_sub_amatrix(T, A_t, rows, 0, k2, 0);
  kernels->fundamental_row(t->idx, (const real(*)[2]) Z, (pcbem2d) bem, T);
  uninit_amatrix(T);

  init_sub_amatrix(T, A_t, rows, 0, k2, k2);
  kernels->dnz_fundamental_row(t->idx, (const real(*)[2]) Z,
			       (const real(*)[2]) N, (pcbem2d) bem, T);
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
  update_pivotelements_greenclusterbasis2d(par->grbn + rname, rb, xi, rank);

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
assemble_bem2d_greenhybridortho_transfer_row_clusterbasis(pcbem2d bem,
							  pclusterbasis rb,
							  uint rname)
{
  paprxbem2d aprx = bem->aprx;
  pkernelbem2d kernels = bem->kernels;
  pparbem2d par = bem->par;
  pccluster t = rb->t;
  const real eps = aprx->accur_aca;
  real      delta = aprx->delta_green;
  uint      rank = aprx->k_green;
  uint      sons = t->sons;

  prkmatrix R;
  pamatrix  A_t, T, RC, E, Q;
  pavector  tau;
  pgreenclusterbasis2d grb, grb1;
  real(*Z)[2], (*N)[2];
  uint     *I_t, *idx;
  uint      i, k2, newrank, rname1;

  k2 = rank * 0.5;
  delta = aprx->delta_green
    * ((aprx->quadpoints == build_bem2d_rect_quadpoints) ?
       getdiam_max_cluster(t) : getdiam_2_cluster(t));

  aprx->quadpoints(bem, t->bmin, t->bmax, delta, &Z, &N);

  T = new_amatrix(0, 0);

  grb = par->grbn[rname];
  if (grb == NULL) {
    grb = par->grbn[rname] = new_greenclusterbasis2d(rb);
  }

  idx = collect_pivotelements_greenclusterbasis2d(par->grbn + rname, &rank);

  /* setup up A_t from green's formula with reduced row indices. */
  A_t = new_amatrix(rank, 2 * k2);

  init_sub_amatrix(T, A_t, rank, 0, k2, 0);
  kernels->fundamental_row(idx, (const real(*)[2]) Z, (pcbem2d) bem, T);
  uninit_amatrix(T);

  init_sub_amatrix(T, A_t, rank, 0, k2, k2);
  kernels->dnz_fundamental_row(idx, (const real(*)[2]) Z,
			       (const real(*)[2]) N, (pcbem2d) bem, T);
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
  update_pivotelements_greenclusterbasis2d(par->grbn + rname, rb, I_t,
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
assemble_bem2d_greenhybridortho_leaf_col_clusterbasis(pcbem2d bem,
						      pclusterbasis cb,
						      uint cname)
{
  paprxbem2d aprx = bem->aprx;
  pkernelbem2d kernels = bem->kernels;
  pparbem2d par = bem->par;
  pccluster t = cb->t;
  pamatrix  V = &cb->V;
  const uint rows = t->size;
  const real eps = aprx->accur_aca;
  real      delta = aprx->delta_green;
  uint      rank = aprx->k_green;

  prkmatrix R;
  pamatrix  A_t, T, RC;
  pavector  tau;
  pgreenclusterbasis2d gcb;
  real(*Z)[2], (*N)[2];
  uint     *xi;
  uint      k2;

  k2 = rank * 0.5;
  delta = aprx->delta_green
    * ((aprx->quadpoints == build_bem2d_rect_quadpoints) ?
       getdiam_max_cluster(t) : getdiam_2_cluster(t));

  aprx->quadpoints(bem, t->bmin, t->bmax, delta, &Z, &N);

  T = new_amatrix(0, 0);

  gcb = par->gcbn[cname];
  if (gcb == NULL) {
    gcb = par->gcbn[cname] = new_greenclusterbasis2d(cb);
  }

  /* setup up A_t from green's formula. */
  A_t = new_amatrix(rows, rank);

  init_sub_amatrix(T, A_t, rows, 0, k2, 0);
  kernels->kernel_col(t->idx, (const real(*)[2]) Z, (pcbem2d) bem, T);
  uninit_amatrix(T);

  init_sub_amatrix(T, A_t, rows, 0, k2, k2);
  kernels->dnz_kernel_col(t->idx, (const real(*)[2]) Z, (const real(*)[2]) N,
			  (pcbem2d) bem, T);
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
  update_pivotelements_greenclusterbasis2d(par->gcbn + cname, cb, xi, rank);

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
assemble_bem2d_greenhybridortho_transfer_col_clusterbasis(pcbem2d bem,
							  pclusterbasis cb,
							  uint cname)
{
  paprxbem2d aprx = bem->aprx;
  pkernelbem2d kernels = bem->kernels;
  pparbem2d par = bem->par;
  pccluster t = cb->t;
  const real eps = aprx->accur_aca;
  real      delta = aprx->delta_green;
  uint      rank = aprx->k_green;
  uint      sons = t->sons;

  prkmatrix R;
  pamatrix  A_t, T, RC, E, Q;
  pavector  tau;
  pgreenclusterbasis2d gcb, gcb1;
  real(*Z)[2], (*N)[2];
  uint     *I_t, *idx;
  uint      i, k2, newrank, cname1;

  k2 = rank * 0.5;
  delta = aprx->delta_green
    * ((aprx->quadpoints == build_bem2d_rect_quadpoints) ?
       getdiam_max_cluster(t) : getdiam_2_cluster(t));

  aprx->quadpoints(bem, t->bmin, t->bmax, delta, &Z, &N);

  T = new_amatrix(0, 0);

  gcb = par->gcbn[cname];
  if (gcb == NULL) {
    gcb = par->gcbn[cname] = new_greenclusterbasis2d(cb);
  }

  idx = collect_pivotelements_greenclusterbasis2d(par->gcbn + cname, &rank);

  /* setup up A_t from green's formula with reduced row indices. */
  A_t = new_amatrix(rank, 2 * k2);

  init_sub_amatrix(T, A_t, rank, 0, k2, 0);
  kernels->kernel_col(idx, (const real(*)[2]) Z, (pcbem2d) bem, T);
  uninit_amatrix(T);

  init_sub_amatrix(T, A_t, rank, 0, k2, k2);
  kernels->dnz_kernel_col(idx, (const real(*)[2]) Z, (const real(*)[2]) N,
			  (pcbem2d) bem, T);
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
  update_pivotelements_greenclusterbasis2d(par->gcbn + cname, cb, I_t,
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
assemble_bem2d_greenhybridortho_uniform(uint rname, uint cname,
					pcbem2d bem, puniform U)
{
  pparbem2d par = bem->par;
  const uint kr = U->rb->k;
  const uint kc = U->cb->k;
  pamatrix  S = &U->S;

  pgreenclusterbasis2d grb, gcb;
  uint     *xihatV, *xihatW;

  grb = par->grbn[rname];
  gcb = par->gcbn[cname];

  assert(grb != NULL);
  assert(gcb != NULL);

  xihatV = grb->xihat;
  xihatW = gcb->xihat;

  resize_amatrix(S, kr, kc);

  bem->nearfield(xihatV, xihatW, bem, false, S);

  triangulareval_amatrix(false, false, false, grb->Qinv, false, S);
  triangulareval_amatrix(false, false, false, gcb->Qinv, true, S);
}

/* ------------------------------------------------------------
 Methods to evaluate or integrate lagrange polynomials
 ------------------------------------------------------------ */

void
assemble_bem2d_lagrange_amatrix(const real(*X)[2], pcavector px,
				pcavector py, pamatrix V)
{
  const uint rows = V->rows;
  const uint cols = V->cols;
  const uint ld = V->ld;
  const uint mx = px->dim;
  const uint my = py->dim;

  real     *denomx, *denomy;
  uint      jx, jy, i, l, index;
  real      lagr, denom;

  denomx = allocreal(mx);
  denomy = allocreal(my);

  /*
   * Eval Lagrange polynomials at points X
   */

  assert(mx * my == cols);

  for (jx = 0; jx < mx; ++jx) {
    denom = 1.0;
    for (l = 0; l < jx; ++l) {
      denom *= (px->v[jx] - px->v[l]);
    }
    for (l = jx + 1; l < mx; ++l) {
      denom *= (px->v[jx] - px->v[l]);
    }
    denomx[jx] = 1.0 / denom;
  }

  for (jy = 0; jy < my; ++jy) {
    denom = 1.0;
    for (l = 0; l < jy; ++l) {
      denom *= (py->v[jy] - py->v[l]);
    }
    for (l = jy + 1; l < my; ++l) {
      denom *= (py->v[jy] - py->v[l]);
    }
    denomy[jy] = 1.0 / denom;
  }

  index = 0;
  for (jx = 0; jx < mx; ++jx) {
    for (jy = 0; jy < my; ++jy) {
      denom = denomx[jx] * denomy[jy];
      for (i = 0; i < rows; ++i) {

	lagr = 1.0;

	for (l = 0; l < jx; ++l) {
	  lagr *= (X[i][0] - px->v[l]);
	}
	for (l = jx + 1; l < mx; ++l) {
	  lagr *= (X[i][0] - px->v[l]);
	}

	for (l = 0; l < jy; ++l) {
	  lagr *= (X[i][1] - py->v[l]);
	}
	for (l = jy + 1; l < my; ++l) {
	  lagr *= (X[i][1] - py->v[l]);
	}

	V->a[i + index * ld] = lagr * denom;
      }
      index++;
    }
  }

  freemem(denomx);
  freemem(denomy);
}

void
assemble_bem2d_lagrange_const_amatrix(const uint * idx, pcavector px,
				      pcavector py, pcbem2d bem, pamatrix V)
{
  pccurve2d gr = bem->gr;
  const     real(*gr_x)[2] = (const real(*)[2]) gr->x;
  const     uint(*gr_e)[2] = (const uint(*)[2]) gr->e;
  const preal gr_g = (const preal) gr->g;
  uint      rows = V->rows;
  uint      ld = V->ld;
  uint      nq = bem->sq->n_single;
  real     *xx = bem->sq->x_single;
  real     *ww = bem->sq->w_single;
  uint      mx = px->dim;
  uint      my = py->dim;

  const real *A, *B;
  real     *denomx, *denomy;
  uint      t, tt, jx, jy, q, l, index;
  real      gt, sum, lagr, denom, x, y, tx, Ax, Bx;

  denomx = allocreal(mx);
  denomy = allocreal(my);

  /*
   * integrate Lagrange polynomials with constant basisfunctions
   */

  for (jx = 0; jx < mx; ++jx) {
    denom = 1.0;
    for (l = 0; l < jx; ++l) {
      denom *= (px->v[jx] - px->v[l]);
    }
    for (l = jx + 1; l < mx; ++l) {
      denom *= (px->v[jx] - px->v[l]);
    }
    denomx[jx] = 1.0 / denom;
  }

  for (jy = 0; jy < my; ++jy) {
    denom = 1.0;
    for (l = 0; l < jy; ++l) {
      denom *= (py->v[jy] - py->v[l]);
    }
    for (l = jy + 1; l < my; ++l) {
      denom *= (py->v[jy] - py->v[l]);
    }
    denomy[jy] = 1.0 / denom;
  }

  index = 0;
  for (jx = 0; jx < mx; ++jx) {
    for (jy = 0; jy < my; ++jy) {
      denom = denomx[jx] * denomy[jy];
      for (t = 0; t < rows; ++t) {
	tt = (idx == NULL ? t : idx[t]);
	gt = gr_g[tt];
	A = gr_x[gr_e[tt][0]];
	B = gr_x[gr_e[tt][1]];

	sum = 0.0;

	for (q = 0; q < nq; ++q) {
	  tx = xx[q];
	  Ax = 1.0 - tx;
	  Bx = tx;

	  x = A[0] * Ax + B[0] * Bx;
	  y = A[1] * Ax + B[1] * Bx;

	  lagr = 1.0;

	  for (l = 0; l < jx; ++l) {
	    lagr *= (x - px->v[l]);
	  }
	  for (l = jx + 1; l < mx; ++l) {
	    lagr *= (x - px->v[l]);
	  }

	  for (l = 0; l < jy; ++l) {
	    lagr *= (y - py->v[l]);
	  }
	  for (l = jy + 1; l < my; ++l) {
	    lagr *= (y - py->v[l]);
	  }

	  sum += ww[q] * lagr;

	}
	V->a[t + index * ld] = gt * sum * denom;
      }
      index++;
    }
  }

  freemem(denomx);
  freemem(denomy);
}

void
assemble_bem2d_dn_lagrange_const_amatrix(const uint * idx, pcavector px,
					 pcavector py, pcbem2d bem,
					 pamatrix V)
{
  pccurve2d gr = bem->gr;
  const     real(*gr_x)[2] = (const real(*)[2]) gr->x;
  const     uint(*gr_e)[2] = (const uint(*)[2]) gr->e;
  const     real(*gr_n)[2] = (const real(*)[2]) gr->n;
  const preal gr_g = (const preal) gr->g;
  uint      rows = V->rows;
  uint      ld = V->ld;
  uint      nq = bem->sq->n_single;
  real     *xx = bem->sq->x_single;
  real     *ww = bem->sq->w_single;
  uint      mx = px->dim;
  uint      my = py->dim;
  uint      cols = V->cols;

  const real *A, *B, *nt;
  real      lagr[2];
  uint      t, tt, jx, jy, q, l, index;
  real      gt, sum, lagrx, lagry, x, y, tx, Ax, Bx;

  assert(cols == mx * my);

  /*
   * integrate Lagrange polynomials with constant basisfunctions
   */

  for (t = 0; t < rows; ++t) {
    tt = (idx == NULL ? t : idx[t]);
    gt = gr_g[tt];
    nt = gr_n[tt];
    A = gr_x[gr_e[tt][0]];
    B = gr_x[gr_e[tt][1]];

    index = 0;

    for (jx = 0; jx < mx; ++jx) {
      for (jy = 0; jy < my; ++jy) {

	sum = 0.0;

	for (q = 0; q < nq; ++q) {
	  tx = xx[q];
	  Ax = 1.0 - tx;
	  Bx = tx;

	  x = A[0] * Ax + B[0] * Bx;
	  y = A[1] * Ax + B[1] * Bx;

	  lagrx = 1.0;
	  lagry = 1.0;

	  /* TODO optimize lagrange polynomial eval */
	  for (l = 0; l < jx; ++l) {
	    lagrx *= (x - px->v[l]) / (px->v[jx] - px->v[l]);
	  }
	  for (l = jx + 1; l < mx; ++l) {
	    lagrx *= (x - px->v[l]) / (px->v[jx] - px->v[l]);
	  }

	  for (l = 0; l < jy; ++l) {
	    lagry *= (y - py->v[l]) / (py->v[jy] - py->v[l]);
	  }
	  for (l = jy + 1; l < my; ++l) {
	    lagry *= (y - py->v[l]) / (py->v[jy] - py->v[l]);
	  }

	  lagr[0] = 0.0;
	  lagr[1] = 0.0;

	  for (l = 0; l < jx; ++l) {
	    lagr[0] += lagrx / (x - px->v[l]);
	  }
	  for (l = jx + 1; l < mx; ++l) {
	    lagr[0] += lagrx / (x - px->v[l]);
	  }

	  for (l = 0; l < jy; ++l) {
	    lagr[1] += lagry / (y - py->v[l]);
	  }
	  for (l = jy + 1; l < my; ++l) {
	    lagr[1] += lagry / (y - py->v[l]);
	  }

	  lagr[0] *= lagry;
	  lagr[1] *= lagrx;

	  sum += ww[q] * (nt[0] * lagr[0] + nt[1] * lagr[1]);

	}

	V->a[t + index * ld] = -1.0 * gt * sum;
	index++;

      }
    }
  }
}

void
projectl2_bem2d_const_avector(pbem2d bem,
			      field(*rhs) (const real * x, const real * n),
			      pavector f)
{
  pccurve2d gr = bem->gr;
  const     real(*gr_x)[2] = (const real(*)[2]) gr->x;
  const     uint(*gr_e)[2] = (const uint(*)[2]) gr->e;
  const     real(*gr_n)[2] = (const real(*)[2]) gr->n;
  uint      nq = bem->sq->n_single;
  real     *xx = bem->sq->x_single;
  real     *ww = bem->sq->w_single;
  uint      cols = f->dim;

  const real *A, *B, *n;
  real      x[2];
  real      tx, Ax, Bx;
  field     sum;
  uint      t, q;

  /*
   *  integrate function with constant basisfunctions
   */

  for (t = 0; t < cols; ++t) {
    A = gr_x[gr_e[t][0]];
    B = gr_x[gr_e[t][1]];
    n = gr_n[t];

    sum = 0.0;

    for (q = 0; q < nq; ++q) {
      tx = xx[q];
      Ax = 1.0 - tx;
      Bx = tx;

      x[0] = A[0] * Ax + B[0] * Bx;
      x[1] = A[1] * Ax + B[1] * Bx;

      sum += ww[q] * rhs(x, n);
    }

    f->v[t] = sum;
  }
}

prkmatrix
build_bem2d_rkmatrix(pccluster row, pccluster col, void *data)
{
  pbem2d    bem = (pbem2d) data;

  prkmatrix R;

  R = new_rkmatrix(row->size, col->size, 0);
  bem->farfield_rk(row, 0, col, 0, bem, R);

  return R;
}

pamatrix
build_bem2d_amatrix(pccluster row, pccluster col, void *data)
{
  pbem2d    bem = (pbem2d) data;

  pamatrix  N;

  N = new_amatrix(row->size, col->size);
  bem->nearfield(row->idx, col->idx, bem, false, N);

  return N;
}

/* ------------------------------------------------------------
 Methods to initialize an approximation technique
 ------------------------------------------------------------ */

void
setup_hmatrix_recomp_bem2d(pbem2d bem, bool recomp, real accur_recomp,
			   bool coarsen, real accur_coarsen)
{
  paprxbem2d aprx = bem->aprx;
  assert(accur_recomp >= 0.0);
  assert(accur_coarsen >= 0.0);

  uninit_recompression_bem2d(aprx);

  aprx->recomp = recomp;
  aprx->accur_recomp = accur_recomp;
  aprx->coarsen = coarsen;
  aprx->accur_coarsen = accur_coarsen;
}

void
setup_hmatrix_aprx_inter_row_bem2d(pbem2d bem, pccluster rc, pccluster cc,
				   pcblock tree, uint m)
{
  (void) rc;
  (void) cc;
  (void) tree;

  assert(bem->kernels->kernel_col != NULL);
  assert(bem->kernels->lagrange_row != NULL);

  setup_interpolation_bem2d(bem->aprx, m);

  bem->farfield_rk = assemble_bem2d_inter_row_rkmatrix;
  bem->farfield_u = NULL;

  bem->leaf_row = NULL;
  bem->leaf_col = NULL;
  bem->transfer_row = NULL;
  bem->transfer_col = NULL;
}

void
setup_hmatrix_aprx_inter_col_bem2d(pbem2d bem, pccluster rc, pccluster cc,
				   pcblock tree, uint m)
{

  (void) rc;
  (void) cc;
  (void) tree;

  assert(bem->kernels->fundamental_row != NULL);
  assert(bem->kernels->lagrange_col != NULL);

  setup_interpolation_bem2d(bem->aprx, m);

  bem->farfield_rk = assemble_bem2d_inter_col_rkmatrix;
  bem->farfield_u = NULL;

  bem->leaf_row = NULL;
  bem->leaf_col = NULL;
  bem->transfer_row = NULL;
  bem->transfer_col = NULL;
}

void
setup_hmatrix_aprx_inter_mixed_bem2d(pbem2d bem, pccluster rc,
				     pccluster cc, pcblock tree, uint m)
{

  (void) rc;
  (void) cc;
  (void) tree;

  assert(bem->kernels->kernel_col != NULL);
  assert(bem->kernels->lagrange_row != NULL);
  assert(bem->kernels->fundamental_row != NULL);
  assert(bem->kernels->lagrange_col != NULL);

  setup_interpolation_bem2d(bem->aprx, m);

  bem->farfield_rk = assemble_bem2d_inter_mixed_rkmatrix;
  bem->farfield_u = NULL;

  bem->leaf_row = NULL;
  bem->leaf_col = NULL;
  bem->transfer_row = NULL;
  bem->transfer_col = NULL;
}

void
setup_hmatrix_aprx_green_row_bem2d(pbem2d bem, pccluster rc, pccluster cc,
				   pcblock tree, uint m, uint l, real delta,
				   quadpoints2d quadpoints)
{

  (void) rc;
  (void) cc;
  (void) tree;

  assert(quadpoints != NULL);
  assert(bem->kernels->fundamental_row != NULL);
  assert(bem->kernels->kernel_col != NULL);
  assert(bem->kernels->dnz_fundamental_row != NULL);
  assert(bem->kernels->dnz_kernel_col != NULL);

  setup_green_bem2d(bem->aprx, m, l, delta, quadpoints);

  bem->farfield_rk = assemble_bem2d_green_row_rkmatrix;
  bem->farfield_u = NULL;

  bem->leaf_row = NULL;
  bem->leaf_col = NULL;
  bem->transfer_row = NULL;
  bem->transfer_col = NULL;
}

void
setup_hmatrix_aprx_green_col_bem2d(pbem2d bem, pccluster rc, pccluster cc,
				   pcblock tree, uint m, uint l, real delta,
				   quadpoints2d quadpoints)
{

  (void) rc;
  (void) cc;
  (void) tree;

  assert(quadpoints != NULL);
  assert(bem->kernels->fundamental_row != NULL);
  assert(bem->kernels->kernel_col != NULL);
  assert(bem->kernels->dnz_fundamental_row != NULL);
  assert(bem->kernels->dnz_kernel_col != NULL);

  setup_green_bem2d(bem->aprx, m, l, delta, quadpoints);

  bem->farfield_rk = assemble_bem2d_green_col_rkmatrix;
  bem->farfield_u = NULL;

  bem->leaf_row = NULL;
  bem->leaf_col = NULL;
  bem->transfer_row = NULL;
  bem->transfer_col = NULL;
}

void
setup_hmatrix_aprx_green_mixed_bem2d(pbem2d bem, pccluster rc,
				     pccluster cc, pcblock tree, uint m,
				     uint l, real delta,
				     quadpoints2d quadpoints)
{

  (void) rc;
  (void) cc;
  (void) tree;

  assert(quadpoints != NULL);
  assert(bem->kernels->fundamental_row != NULL);
  assert(bem->kernels->kernel_col != NULL);
  assert(bem->kernels->dnz_fundamental_row != NULL);
  assert(bem->kernels->dnz_kernel_col != NULL);

  setup_green_bem2d(bem->aprx, m, l, delta, quadpoints);

  bem->farfield_rk = assemble_bem2d_green_mixed_rkmatrix;
  bem->farfield_u = NULL;

  bem->leaf_row = NULL;
  bem->leaf_col = NULL;
  bem->transfer_row = NULL;
  bem->transfer_col = NULL;
}

void
setup_hmatrix_aprx_greenhybrid_row_bem2d(pbem2d bem, pccluster rc,
					 pccluster cc, pcblock tree, uint m,
					 uint l, real delta, real accur,
					 quadpoints2d quadpoints)
{
  pparbem2d par = bem->par;

  uint      n, i;

  (void) cc;
  (void) tree;

  assert(quadpoints != NULL);
  assert(bem->kernels->fundamental_row != NULL);
  assert(bem->kernels->dnz_fundamental_row != NULL);
  assert(bem->nearfield != NULL);

  setup_green_bem2d(bem->aprx, m, l, delta, quadpoints);
  setup_aca_bem2d(bem->aprx, accur);

  bem->farfield_rk = assemble_bem2d_greenhybrid_row_rkmatrix;
  bem->farfield_u = NULL;

  bem->leaf_row = NULL;
  bem->leaf_col = NULL;
  bem->transfer_row = NULL;
  bem->transfer_col = NULL;

  n = par->grcnn;

  if (par->grcn != NULL && par->grcnn != 0) {
    for (i = 0; i < n; ++i) {
      if (par->grcn[i] != NULL) {
	del_greencluster2d(par->grcn[i]);
      }
    }
    freemem(par->grcn);
  }

  n = rc->desc;

  par->grcn = (pgreencluster2d *) allocmem((size_t) n *
					   sizeof(pgreencluster2d));

  for (i = 0; i < n; ++i) {
    par->grcn[i] = NULL;
  }

  par->grcnn = n;
}

void
setup_hmatrix_aprx_greenhybrid_col_bem2d(pbem2d bem, pccluster rc,
					 pccluster cc, pcblock tree, uint m,
					 uint l, real delta, real accur,
					 quadpoints2d quadpoints)
{
  pparbem2d par = bem->par;

  uint      i, n;

  (void) rc;
  (void) tree;

  assert(quadpoints != NULL);
  assert(bem->kernels->kernel_col != NULL);
  assert(bem->kernels->dnz_kernel_col != NULL);
  assert(bem->nearfield != NULL);

  setup_green_bem2d(bem->aprx, m, l, delta, quadpoints);
  setup_aca_bem2d(bem->aprx, accur);

  bem->farfield_rk = assemble_bem2d_greenhybrid_col_rkmatrix;
  bem->farfield_u = NULL;

  bem->leaf_row = NULL;
  bem->leaf_col = NULL;
  bem->transfer_row = NULL;
  bem->transfer_col = NULL;

  n = par->gccnn;

  if (par->gccn != NULL && par->gccnn != 0) {
    for (i = 0; i < n; ++i) {
      if (par->gccn[i] != NULL) {
	del_greencluster2d(par->gccn[i]);
      }
    }
    freemem(par->gccn);
  }

  n = cc->desc;

  par->gccn = (pgreencluster2d *) allocmem((size_t) n *
					   sizeof(pgreencluster2d));

  for (i = 0; i < n; ++i) {
    par->gccn[i] = NULL;
  }

  par->gccnn = n;
}

void
setup_hmatrix_aprx_greenhybrid_mixed_bem2d(pbem2d bem, pccluster rc,
					   pccluster cc, pcblock tree, uint m,
					   uint l, real delta, real accur,
					   quadpoints2d quadpoints)
{
  pparbem2d par = bem->par;

  uint      i, n;

  (void) tree;

  assert(quadpoints != NULL);
  assert(bem->kernels->fundamental_row != NULL);
  assert(bem->kernels->dnz_fundamental_row != NULL);
  assert(bem->kernels->kernel_col != NULL);
  assert(bem->kernels->dnz_kernel_col != NULL);
  assert(bem->nearfield != NULL);

  setup_green_bem2d(bem->aprx, m, l, delta, quadpoints);
  setup_aca_bem2d(bem->aprx, accur);

  bem->farfield_rk = assemble_bem2d_greenhybrid_mixed_rkmatrix;
  bem->farfield_u = NULL;

  bem->leaf_row = NULL;
  bem->leaf_col = NULL;
  bem->transfer_row = NULL;
  bem->transfer_col = NULL;

  n = par->grcnn;

  if (par->grcn != NULL && par->grcnn > 0) {
    for (i = 0; i < n; ++i) {
      if (par->grcn[i] != NULL) {
	del_greencluster2d(par->grcn[i]);
      }
    }
    freemem(par->grcn);
  }

  n = rc->desc;

  par->grcn = (pgreencluster2d *) allocmem((size_t) n *
					   sizeof(pgreencluster2d));

  for (i = 0; i < n; ++i) {
    par->grcn[i] = NULL;
  }

  par->grcnn = n;

  n = par->gccnn;

  if (par->gccn != NULL && par->gccnn > 0) {
    for (i = 0; i < n; ++i) {
      if (par->gccn[i] != NULL) {
	del_greencluster2d(par->gccn[i]);
      }
    }
    freemem(par->gccn);
  }

  n = cc->desc;

  par->gccn = (pgreencluster2d *) allocmem((size_t) n *
					   sizeof(pgreencluster2d));

  for (i = 0; i < n; ++i) {
    par->gccn[i] = NULL;
  }

  par->gccnn = n;
}

void
setup_hmatrix_aprx_aca_bem2d(pbem2d bem, pccluster rc, pccluster cc,
			     pcblock tree, real accur)
{

  (void) rc;
  (void) cc;
  (void) tree;

  assert(bem->nearfield != NULL);

  setup_aca_bem2d(bem->aprx, accur);

  bem->farfield_rk = assemble_bem2d_ACA_rkmatrix;
  bem->farfield_u = NULL;

  bem->leaf_row = NULL;
  bem->leaf_col = NULL;
  bem->transfer_row = NULL;
  bem->transfer_col = NULL;
}

void
setup_hmatrix_aprx_paca_bem2d(pbem2d bem, pccluster rc, pccluster cc,
			      pcblock tree, real accur)
{

  (void) rc;
  (void) cc;
  (void) tree;

  assert(bem->nearfield != NULL);

  setup_aca_bem2d(bem->aprx, accur);

  bem->farfield_rk = assemble_bem2d_PACA_rkmatrix;
  bem->farfield_u = NULL;

  bem->leaf_row = NULL;
  bem->leaf_col = NULL;
  bem->transfer_row = NULL;
  bem->transfer_col = NULL;
}

void
setup_hmatrix_aprx_hca_bem2d(pbem2d bem, pccluster rc, pccluster cc,
			     pcblock tree, uint m, real accur)
{

  (void) rc;
  (void) cc;
  (void) tree;

  assert(bem->kernels->fundamental != NULL);
  assert(bem->kernels->kernel_row != NULL);
  assert(bem->kernels->kernel_col != NULL);

  setup_interpolation_bem2d(bem->aprx, m);
  setup_aca_bem2d(bem->aprx, accur);

  bem->farfield_rk = assemble_bem2d_HCA_rkmatrix;
  bem->farfield_u = NULL;

  bem->leaf_row = NULL;
  bem->leaf_col = NULL;
  bem->transfer_row = NULL;
  bem->transfer_col = NULL;
}

void
setup_h2matrix_recomp_bem2d(pbem2d bem, bool hiercomp, real accur_hiercomp)
{
  paprxbem2d aprx = bem->aprx;

  assert(accur_hiercomp >= 0.0);

  aprx->hiercomp = hiercomp;
  aprx->accur_hiercomp = accur_hiercomp;
  aprx->tm = new_releucl_truncmode();
}

void
setup_h2matrix_aprx_inter_bem2d(pbem2d bem, pcclusterbasis rb,
				pcclusterbasis cb, pcblock tree, uint m)
{

  (void) rb;
  (void) cb;
  (void) tree;

  assert(bem->kernels->lagrange_row != NULL);
  assert(bem->kernels->lagrange_col != NULL);
  assert(bem->kernels->fundamental != NULL);

  setup_interpolation_bem2d(bem->aprx, m);

  bem->farfield_rk = NULL;
  bem->farfield_u = assemble_bem2d_inter_uniform;

  bem->leaf_row = assemble_bem2d_inter_leaf_row_clusterbasis;
  bem->leaf_col = assemble_bem2d_inter_leaf_col_clusterbasis;
  bem->transfer_row = assemble_bem2d_inter_tranfer_clusterbasis;
  bem->transfer_col = assemble_bem2d_inter_tranfer_clusterbasis;
}

void
setup_h2matrix_aprx_greenhybrid_bem2d(pbem2d bem, pcclusterbasis rb,
				      pcclusterbasis cb, pcblock tree, uint m,
				      uint l, real delta, real accur,
				      quadpoints2d quadpoints)
{
  pparbem2d par = bem->par;

  uint      i, n;

  (void) tree;

  assert(quadpoints != NULL);
  assert(bem->kernels->fundamental_row != NULL);
  assert(bem->kernels->dnz_fundamental_row != NULL);
  assert(bem->kernels->kernel_col != NULL);
  assert(bem->kernels->dnz_kernel_col != NULL);
  assert(bem->nearfield != NULL);

  setup_green_bem2d(bem->aprx, m, l, delta, quadpoints);
  setup_aca_bem2d(bem->aprx, accur);

  bem->farfield_rk = NULL;
  bem->farfield_u = assemble_bem2d_greenhybrid_uniform;

  bem->leaf_row = assemble_bem2d_greenhybrid_leaf_row_clusterbasis;
  bem->leaf_col = assemble_bem2d_greenhybrid_leaf_col_clusterbasis;
  bem->transfer_row = assemble_bem2d_greenhybrid_transfer_row_clusterbasis;
  bem->transfer_col = assemble_bem2d_greenhybrid_transfer_col_clusterbasis;

  n = par->grbnn;

  if (par->grbn != NULL && par->grbnn != 0) {
    for (i = 0; i < n; ++i) {
      if (par->grbn[i] != NULL) {
	del_greenclusterbasis2d(par->grbn[i]);
      }
    }
    freemem(par->grbn);
  }

  n = rb->t->desc;

  par->grbn = (pgreenclusterbasis2d *) allocmem((size_t) n *
						sizeof(pgreenclusterbasis2d));

  for (i = 0; i < n; ++i) {
    par->grbn[i] = NULL;
  }

  par->grbnn = n;

  n = par->gcbnn;

  if (par->gcbn != NULL && par->gcbnn != 0) {
    for (i = 0; i < n; ++i) {
      if (par->gcbn[i] != NULL) {
	del_greenclusterbasis2d(par->gcbn[i]);
      }
    }
    freemem(par->gcbn);
  }

  n = cb->t->desc;

  par->gcbn = (pgreenclusterbasis2d *) allocmem((size_t) n *
						sizeof(pgreenclusterbasis2d));

  for (i = 0; i < n; ++i) {
    par->gcbn[i] = NULL;
  }

  par->gcbnn = n;
}

void
setup_h2matrix_aprx_greenhybrid_ortho_bem2d(pbem2d bem, pcclusterbasis rb,
					    pcclusterbasis cb, pcblock tree,
					    uint m, uint l, real delta,
					    real accur,
					    quadpoints2d quadpoints)
{
  pparbem2d par = bem->par;

  uint      i, n;

  (void) tree;

  assert(quadpoints != NULL);
  assert(bem->kernels->fundamental_row != NULL);
  assert(bem->kernels->dnz_fundamental_row != NULL);
  assert(bem->kernels->kernel_col != NULL);
  assert(bem->kernels->dnz_kernel_col != NULL);
  assert(bem->nearfield != NULL);

  setup_green_bem2d(bem->aprx, m, l, delta, quadpoints);
  setup_aca_bem2d(bem->aprx, accur);

  bem->farfield_rk = NULL;
  bem->farfield_u = assemble_bem2d_greenhybridortho_uniform;

  bem->leaf_row = assemble_bem2d_greenhybridortho_leaf_row_clusterbasis;
  bem->leaf_col = assemble_bem2d_greenhybridortho_leaf_col_clusterbasis;
  bem->transfer_row =
    assemble_bem2d_greenhybridortho_transfer_row_clusterbasis;
  bem->transfer_col =
    assemble_bem2d_greenhybridortho_transfer_col_clusterbasis;
  n = par->grbnn;

  if (par->grbn != NULL && par->grbnn != 0) {
    for (i = 0; i < n; ++i) {
      if (par->grbn[i] != NULL) {
	del_greenclusterbasis2d(par->grbn[i]);
      }
    }
    freemem(par->grbn);
  }

  n = rb->t->desc;

  par->grbn = (pgreenclusterbasis2d *) allocmem((size_t) n *
						sizeof(pgreenclusterbasis2d));

  for (i = 0; i < n; ++i) {
    par->grbn[i] = NULL;
  }

  par->grbnn = n;

  n = par->gcbnn;

  if (par->gcbn != NULL && par->gcbnn != 0) {
    for (i = 0; i < n; ++i) {
      if (par->gcbn[i] != NULL) {
	del_greenclusterbasis2d(par->gcbn[i]);
      }
    }
    freemem(par->gcbn);
  }

  n = cb->t->desc;

  par->gcbn = (pgreenclusterbasis2d *) allocmem((size_t) n *
						sizeof(pgreenclusterbasis2d));

  for (i = 0; i < n; ++i) {
    par->gcbn[i] = NULL;
  }

  par->gcbnn = n;
}

/* ------------------------------------------------------------
 Methods to fill hmatrices
 ------------------------------------------------------------ */

static void
assemble_bem2d_block_hmatrix(pcblock b, uint bname, uint rname,
			     uint cname, uint pardepth, void *data)
{
  pbem2d    bem = (pbem2d) data;
  pparbem2d par = bem->par;
  paprxbem2d aprx = bem->aprx;
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
assemblecoarsen_bem2d_block_hmatrix(pcblock b, uint bname,
				    uint rname, uint cname, uint pardepth,
				    void *data)
{
  pbem2d    bem = (pbem2d) data;
  pparbem2d par = bem->par;
  paprxbem2d aprx = bem->aprx;
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

void
assemble_bem2d_hmatrix(pbem2d bem, pblock b, phmatrix G)
{
  pparbem2d par = bem->par;

  par->hn = enumerate_hmatrix(b, G);

  iterate_byrow_block(b, 0, 0, 0, max_pardepth, NULL,
		      assemble_bem2d_block_hmatrix, bem);

  freemem(par->hn);
  par->hn = NULL;
}

void
assemblecoarsen_bem2d_hmatrix(pbem2d bem, pblock b, phmatrix G)
{
  pparbem2d par = bem->par;

  par->hn = enumerate_hmatrix(b, G);

  iterate_byrow_block(b, 0, 0, 0, max_pardepth, NULL,
		      assemblecoarsen_bem2d_block_hmatrix, bem);

  freemem(par->hn);
  par->hn = NULL;
}

/* ------------------------------------------------------------
 Methods to fill clusterbasis and h2matrices.
 ------------------------------------------------------------ */

static void
assemble_bem2d_block_h2matrix(pcblock b, uint bname, uint rname,
			      uint cname, uint pardepth, void *data)
{
  pbem2d    bem = (pbem2d) data;
  pparbem2d par = bem->par;
  ph2matrix *h2n = par->h2n;
  ph2matrix G = h2n[bname];

  (void) b;
  (void) pardepth;

  if (G->u) {
    bem->farfield_u(rname, cname, bem, G->u);
  }
  else if (G->f) {
    bem->nearfield(G->rb->t->idx, G->cb->t->idx, bem, false, G->f);
  }
}

static void
assemblehiercomp_bem2d_block_h2matrix(pcblock b, uint bname,
				      uint rname, uint cname, uint pardepth,
				      void *data)
{
  pbem2d    bem = (pbem2d) data;
  paprxbem2d aprx = bem->aprx;
  pparbem2d par = bem->par;
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
assemble_bem2d_h2matrix(pbem2d bem, pblock b, ph2matrix G)
{
  pparbem2d par = bem->par;
  par->h2n = enumerate_h2matrix(G);

  iterate_byrow_block(b, 0, 0, 0, max_pardepth, NULL,
		      assemble_bem2d_block_h2matrix, bem);

  freemem(par->h2n);
  par->h2n = NULL;
}

void
assemblehiercomp_bem2d_h2matrix(pbem2d bem, pblock b, ph2matrix G)
{
  paprxbem2d aprx = bem->aprx;
  pparbem2d par = bem->par;
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
		      assemblehiercomp_bem2d_block_h2matrix, bem);

  aprx->accur_hiercomp *= s;

  freemem(par->h2n);
  par->h2n = NULL;
  del_clusteroperator(par->rwn[0]);
  freemem(par->rwn);
  par->rwn = NULL;
  del_clusteroperator(par->cwn[0]);
  freemem(par->cwn);
  par->cwn = NULL;
  freemem(par->leveln);
  par->leveln = NULL;
}

static void
assemble_h2matrix_row_clusterbasis(pcclusterbasis rb, uint rname, void *data)
{
  pcbem2d   bem = (pcbem2d) data;

  if (rb->sons > 0) {
    assert(bem->transfer_row);
    bem->transfer_row(bem, (pclusterbasis) rb, rname);
  }
  else {
    assert(bem->leaf_row);
    bem->leaf_row(bem, (pclusterbasis) rb, rname);
  }
}

void
assemble_bem2d_h2matrix_row_clusterbasis(pcbem2d bem, pclusterbasis rb)
{
  iterate_parallel_clusterbasis((pcclusterbasis) rb, 0, max_pardepth, NULL,
				assemble_h2matrix_row_clusterbasis,
				(void *) bem);
}

static void
assemble_h2matrix_col_clusterbasis(pcclusterbasis cb, uint cname, void *data)
{
  pcbem2d   bem = (pcbem2d) data;

  if (cb->sons > 0) {
    assert(bem->transfer_col);
    bem->transfer_col(bem, (pclusterbasis) cb, cname);
  }
  else {
    assert(bem->leaf_col);
    bem->leaf_col(bem, (pclusterbasis) cb, cname);
  }
}

void
assemble_bem2d_h2matrix_col_clusterbasis(pcbem2d bem, pclusterbasis cb)
{
  iterate_parallel_clusterbasis((pcclusterbasis) cb, 0, max_pardepth, NULL,
				assemble_h2matrix_col_clusterbasis,
				(void *) bem);
}
