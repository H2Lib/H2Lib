
/* ------------------------------------------------------------
 * This is the file "tri2dp1.h" of the H2Lib package.
 * All rights reserved, Nadine Albrecht 2015
 * ------------------------------------------------------------ */

/** @file tri2dp1.h
 *  @author Nadine Albrecht */

#ifndef TRI2DP1_H
#define TRI2DP1_H

/** @defgroup tri2dp1 tri2p1
 *  @brief Piecewise linear space with nodal basis functions on a
 *  two-dimensional triangular mesh.
 *  @{ */

/** @brief Representation of a trial space with piecewise linear
 *  nodal basis functions on a two-dimensional triangular mesh. */
typedef struct _tri2dp1 tri2dp1;

/** @brief Pointer to a @ref tri2dp1 object */
typedef tri2dp1 *ptri2dp1;

/** @brief Pointer to a constant @ref tri2dp1 object */
typedef const tri2dp1 *pctri2dp1;

#include "sparsematrix.h"
#include "tri2d.h"
#include "clustergeometry.h"

/** @brief Representation of a trial space with piecewise linear
 *  nodal basis functions on a two-dimensional triangular mesh.
 *
 *  Since the mesh itself is represented by a @ref tri2d object,
 *  we only have to keep track of which vertices are degrees of
 *  freedom and provide them with consecutive indices for the
 *  algebraic functions.
 *
 *  The remaining vertices are also provided with consecutive
 *  indices, allowing us to construct <em>interaction matrices</em>
 *  that describe how non-zero values in the fixed vertices influence
 *  the right-hand side of a variational equation.
 *  This mechanism is useful, e.g., for handling inhomogeneous Dirichlet
 *  boundary conditions. */
struct _tri2dp1
{
  /** @brief Triangular Mesh, represented by @ref tri2d object.*/
  pctri2d t2;     
  
  /** @brief Number of degrees of freedom.*/
  uint ndof;      
  
  /** @brief Number of fixed vertices.*/
  uint nfix;     
  
  /** @brief Determines whether a vertex is a degree of freedom.*/
  bool *is_dof;  
  
  /** @brief Consecutive indices for all degrees of freedom and all 
   *fixed vertices.*/
  uint *idx2dof; 
};

/** @brief Create a @ref tri2dp1 object using a @ref tri2d mesh.
 *
 *  The boundary flags <tt>gr->xb</tt> are used to decide whether
 *  a vertex corresponds to a degree of freedom or is fixed.
 *
 *  @param gr Mesh
 *  @returns @ref tri2dp1 object matching <tt>gr</tt>. */
HEADER_PREFIX ptri2dp1
new_tri2dp1(pctri2d gr);

/** @brief Delete a @ref tri2dp1 object.
 *
 *  @param dc Object to be deleted. */
HEADER_PREFIX void
del_tri2dp1(ptri2dp1 dc);

/** @brief Create a @ref sparsematrix with a sparsity pattern matching
 *  the nodal basis functions representing the system matrix.
 *
 *  Creates a sparse matrix matching the nodal basis functions:
 *  if the supports of basis functions @f$\varphi_i@f$ and
 *  @f$\varphi_j@f$ share at least one edge, the entry
 *  @f$a_{ij}@f$ appears in the graph of the matrix.
 *
 *  Here, both @f$i@f$ and @f$j@f$ refer to degrees of freedom,
 *  i.e., there indices are between <tt>0</tt> and <tt>ndof-1</tt>.
 *
 *  @param dc @ref tri2dp1 object describing the space.
 *  @returns Sparse matrix with <tt>dc->ndof</tt> rows and columns. */
HEADER_PREFIX psparsematrix
build_tri2dp1_sparsematrix(pctri2dp1 dc);

/** @brief Create a @ref sparsematrix with a sparsity pattern matching
 *  the nodal basis functions representing the interaction matrix.
 *
 *  Creates a sparse matrix matching the nodal basis functions:
 *  if the supports of basis functions @f$\varphi_i@f$ and
 *  @f$\varphi_j@f$ share at least one edge, the entry
 *  @f$a_{ij}@f$ appears in the graph of the matrix.
 *
 *  Here, only @f$i@f$ refers to a degree of freedom, while
 *  @f$j@f$ is the index of a fixed vertex, i.e., has to be
 *  between <tt>0</tt> and <tt>nfix-1</tt>.
 *
 *  @param dc @ref tri2dp1 object describing the space.
 *  @returns Sparse matrix with <tt>dc->ndof</tt> rows and
 *    <tt>dc->nfix</tt> columns. */
HEADER_PREFIX psparsematrix
build_tri2dp1_interaction_sparsematrix(pctri2dp1 dc);

/** @brief Create a prolongation matrix.
 *
 *  <tt>dcoarse</tt> and <tt>dfine</tt> have to correspond to a
 *  coarse mesh and its refinement created by @ref refine_tri2d.
 *  The refinement relationship is described by <tt>rf</tt>.
 *
 *  This function not only creates the matrix, but also fills
 *  it with correct coefficients.
 *
 *  @param dfine @ref tri2dp1 object corresponding to the refined mesh.
 *  @param dcoarse @ref tri2dp1 object corresponding to the coarse mesh.
 *  @param rf Refinement relationship.
 *  @returns @ref sparsematrix containing the prolongation matrix. */
HEADER_PREFIX psparsematrix
build_tri2dp1_prolongation_sparsematrix(pctri2dp1 dfine, pctri2dp1 dcoarse,
                               pctri2dref rf);
   
/** @brief Assemble stiffness matrix.
 *
 *  Element matrices for all triangles are computed and added
 *  to the corresponding entries in <tt>A</tt> and <tt>Af</tt>.
 *
 *  @param dc @ref tri2dp1 object describing the trial space.
 *  @param A Target matrix for degrees of freedom.
 *  @param Af Target matrix for interactions between fixed vertices
 *    and degrees of freedom. */
HEADER_PREFIX void
assemble_tri2dp1_laplace_sparsematrix(pctri2dp1 dc, psparsematrix A, psparsematrix Af);

/** @brief Assemble mass matrix.
 *
 *  Element matrices for all triangles are computed and added
 *  to the corresponding entries in <tt>M</tt> and <tt>Mf</tt>.
 *
 *  @param dc @ref tri2dp1 object describing the trial space.
 *  @param M Target matrix for degrees of freedom.
 *  @param Mf Target matrix for interactions between fixed vertices
 *    and degrees of freedom. */
HEADER_PREFIX void
assemble_tri2dp1_mass_sparsematrix(pctri2dp1 dc, psparsematrix M, psparsematrix Mf);

/** @brief Discretize Dirichlet boundary values.
 *
 *  For each fixed vertex, the function @f$f@f$ is evaluated and
 *  the result is written to the corresponding entry in the result.
 *  In the standard case, this provides us with a coefficient vector
 *  corresponding to the piecewise linear interpolation of the
 *  Dirichlet boundary values.
 *  Subtracting @f$A_f d@f$ from the right-hand side of our equation,
 *  where @f$A_f@f$ is the interaction matrix, yields the right-hand
 *  side of the problem with homogeneous Dirichlet boundary conditions.
 *
 *  @param dc @ref tri2dp1 object describing the trial space.
 *  @param d Dirichlet boundary values.
 *  @param fdata Additional data for <tt>f</tt>.
 *  @param dv Target vector of dimension <tt>dc->nfix</tt>. */
HEADER_PREFIX void
assemble_tri2dp1_dirichlet_avector(pctri2dp1 dc,
		  field (*d)(const real *x, void *fdata), void *fdata,
		  pavector dv);

/** @brief Discretize an @f$L^2@f$-functional.
 *
 *  For each triangle, the element vector is computed by applying
 *  the edge-midpoint quadrature. Its components are then added to
 *  the corresponding entries in @f$b@f$.
 *  Since @f$b@f$ is cleared before this step, the result is given by
 *  @f$b_i = \int f(x) \varphi_i(x) \,dx@f$.
 *
 *  @param dc @ref tri2dp1 object describing the trial space.
 *  @param f Right-hand side function.
 *  @param fdata Additional data for <tt>f</tt>.
 *  @param fv Target vector of dimension <tt>dc->ndof</tt>. */
HEADER_PREFIX void
assemble_tri2dp1_functional_avector(ptri2dp1 dc,
                    field (*f)(const real *x, void *fdata), void *fdata,
                    pavector fv);

/** @brief Compute the vertex-wise maximum norm of the discretization error.
 *
 *  Evaluate <tt>f</tt> in all non-fixed vertices and compute the
 *  difference with the corresponding component of <tt>xs</tt>.
 *  The maximum of the absolute values gives the vertex-wise
 *  @f$\ell^\infty@f$-norm of the discretization error.
 *
 *  @param dc @ref tri2dp1 object describing the trial space.
 *  @param u Function.
 *  @param fdata Additional data for <tt>f</tt>.
 *  @param xs Coefficients of discrete function.
 *  @returns @f$\ell^\infty@f$-norm of discretization error. */
HEADER_PREFIX real
normmax_tri2dp1(pctri2dp1 dc,
		field (*u)(const real *x, void *fdata), void *fdata,
		pavector xs);


/** @brief Compute the @f$L^2@f$-norm of the discretization error.
 *
 *  Compare the function <tt>f</tt> to the discrete function
 *  described by the vector <tt>xs</tt> (for the degrees of freedom)
 *  and <tt>xf</tt> (for the fixed vertices).
 *  Measure the difference in the @f$L^2@f$-norm.
 *
 *  @param dc @ref tri2dp1 object describing the trial space.
 *  @param f Function.
 *  @param fdata Additional data for <tt>f</tt>.
 *  @param xs Coefficients of discrete function.
 *  @param xf Coefficients of fixed part of discrete function
 *  @returns @f$L^2@f$-norm of discretization error. */
HEADER_PREFIX real
norml2_tri2dp1(pctri2dp1 dc,
	       field (*f)(const real *x, void *fdata), void *fdata,
	       pcavector xs, pcavector xf);


/* ------------------------------------------------------------
   Clustergeometry
   ------------------------------------------------------------ */

/** @brief Create a @ref clustergeometry object for a FEM-discretisation with p1-functions
 * 
 * Create a @ref clustergeometry object basing on a triangular \ref tri2d mesh 
 * and linear nodal basis function.
 * 
 * @param p1 tri2dp1 object.
 * @param idx Array of indices.
 * @returns Clustergeometry object.
 */
HEADER_PREFIX pclustergeometry
build_tri2dp1_clustergeometry(pctri2dp1 p1, uint *idx);


/** @} */

#endif
