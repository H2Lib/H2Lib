
/* ------------------------------------------------------------
   This is the file "krylov.h" of the H2Lib package.
   All rights reserved, Steffen Boerm 2010
   ------------------------------------------------------------ */

/** @file krylov.h
 *  @author Steffen B&ouml;rm
 */
 
#ifndef KRYLOV_H
#define KRYLOV_H

#include "settings.h"
#include "avector.h"
#include "factorizations.h"

/** @defgroup krylov krylov
 *  @brief Iterative solvers of Krylov type
 *
 *  @{ */
 
/** @brief Matrix callback.
 *
 *  Used to evaluate the system matrix @f$A@f$ or its adjoint,
 *  i.e., to perform @f$y \gets y + \alpha A x@f$.
 *
 *  Functions like @ref addeval_amatrix_avector,
 *  @ref addevaltrans_amatrix_avector, @ref addeval_hmatrix_avector,
 *  @ref addevaltrans_hmatrix_avector, @ref addeval_h2matrix_avector,
 *  @ref addevaltrans_h2matrix_avector or @ref addeval_sparsematrix_avector
 *  can be cast to <tt>addeval_t</tt>.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param matrix Matrix data describing @f$A@f$.
 *  @param x Source vector @f$x@f$.
 *  @param y Target vector @f$y@f$. */
typedef void (*addeval_t)(field alpha, void *matrix,
			  pcavector x, pavector y);

/** @brief Matrix callback.
 *
 *  Used to evaluate the system matrix @f$A@f$ or its adjoint,
 *  i.e., to perform @f$y \gets y + \alpha A x@f$.
 *
 *  Compared to @ref addeval_t, callbacks of this type allow us
 *  to evaluate both the matrix and its adjoint.
 *
 *  Functions like @ref mvm_amatrix_avector,
 *  @ref mvm_hmatrix_avector, @ref mvm_h2matrix_avector, or
 *  @ref addeval_sparsematrix_avector can be cast to <tt>mvm_t</tt>.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param trans Set if @f$A^*@f$ is to be used instead of @f$A@f$.
 *  @param matrix Matrix data describing @f$A@f$.
 *  @param x Source vector @f$x@f$.
 *  @param y Target vector @f$y@f$. */
typedef void (*mvm_t)(field alpha, bool trans, void *matrix,
		      pcavector x, pavector y);

/** @brief Preconditioner callback.
 *
 *  Used to apply a precondtioner to a vector, i.e., to perform
 *  @f$r \gets N r@f$.
 *
 *  Functions like @ref lrsolve_amatrix_avector, @ref cholsolve_amatrix_avector,
 *  @ref ldltsolve_amatrix_avector, @ref qrsolve_amatrix_avector,
 *  @ref lrsolve_hmatrix_avector or @ref cholsolve_hmatrix_avector can be
 *  cast to <tt>prcd_t</tt>.
 *
 *  @param pdata Preconditioner data describing @f$N@f$.
 *  @param r Source vector, will be overwritten by result. */
typedef void (*prcd_t)(void *pdata, pavector r);

/* ------------------------------------------------------------
 * Conjugate gradients method (cg)
 * ------------------------------------------------------------ */

/** @brief Initialize a standard conjugate gradient method to
 *  solve @f$A x = b@f$.
 *
 *  The matrix @f$A@f$ has to be self-adjoint and positive definite
 *  (or at least positive semidefinite).
 *
 *  @param addeval Callback function name
 *  @param matrix untyped pointer to matrix data
 *  @param b Right-hand side.
 *  @param x Approximate solution.
 *  @param r Residual @f$b-Ax@f$.
 *  @param p Search direction.
 *  @param a Auxiliary vector.
 */
HEADER_PREFIX void
init_cg(addeval_t addeval,
	void *matrix,
	pcavector b,		/* Right-hand side */
	pavector x,		/* Approximate solution */
	pavector r,		/* Residual b-Ax */
	pavector p,		/* Search direction */
	pavector a);

/** @brief One step of a standard conjugate gradient method
 *
 *  @param addeval Callback function name
 *  @param matrix untyped pointer to matrix data (A)
 *  @param b Right-hand side.
 *  @param x Approximate solution.
 *  @param r Residual b-Ax
 *  @param p Search direction
 *  @param a auxiliary vector
 */
HEADER_PREFIX void
step_cg(addeval_t addeval,
	void *matrix,
	pcavector b,		/* Right-hand side */
	pavector x,		/* Approximate solution */
	pavector r,		/* Residual b-Ax */
	pavector p,		/* Search direction */
	pavector a);

/** @brief error information of a standard conjugate gradient method
 *
 *  @param addeval Callback function name
 *  @param matrix untyped pointer to matrix data (A)
 *  @param b Right-hand side.
 *  @param x Approximate solution.
 *  @param r Residual @f$b-Ax@f$
 *  @returns @f$-<Ax,x>_2/2@f$
 */
HEADER_PREFIX real
evalfunctional_cg(addeval_t addeval,
		  void *matrix,
		  pcavector b,
		  pcavector x,
		  pcavector r);

/* ------------------------------------------------------------
 * Preconditioned conjugate gradients method (pcg)
 * ------------------------------------------------------------ */

/** @brief Initialize a preconditioned conjugate gradient method
 *  to solve @f$N^{1/2} A N^{1/2} \widehat{x} = N^{1/2} b@f$
 *  with @f$x = N^{1/2} \widehat{x}@f$.
 *
 *  The matrix @f$A@f$ has to be self-adjoint and positive definite
 *  (or at least positive semidefinite).
 *  The matrix @f$N@f$ has to be self-adjoint and positive definite.
 *
 *  @param addeval Callback function name
 *  @param matrix untyped pointer to matrix data
 *  @param prcd Callback function name
 *  @param pdata untyped pointer to matrix data
 *  @param b Right-hand side.
 *  @param x Approximate solution.
 *  @param r Residual b-Ax
 *  @param q Preconditioned residual
 *  @param p Search direction
 *  @param a auxiliary vector
 */
HEADER_PREFIX void
init_pcg(addeval_t addeval,
	 void *matrix,
	 prcd_t prcd,
	 void *pdata,
	 pcavector b,		/* Right-hand side */
	 pavector x,		/* Approximate solution */
	 pavector r,		/* Residual b-Ax */
	 pavector q,		/* Preconditioned residual */
	 pavector p,		/* Search direction */
	 pavector a);

/** @brief One step of a preconditioned conjugate gradient method
 *
 *  @param addeval Callback function name
 *  @param matrix untyped pointer to matrix data
 *  @param prcd Callback function name
 *  @param pdata untyped pointer to matrix data
 *  @param b Right-hand side.
 *  @param x Approximate solution.
 *  @param r Residual b-Ax
 *  @param q Preconditioned residual
 *  @param p Search direction
 *  @param a auxiliary vector
 */
HEADER_PREFIX void
step_pcg(addeval_t addeval,
	 void *matrix,
	 prcd_t prcd,
	 void *pdata,
	 pcavector b,		/* Right-hand side */
	 pavector x,		/* Approximate solution */
	 pavector r,		/* Residual b-Ax */
	 pavector q,		/* Preconditioned residual */
	 pavector p,		/* Search direction */
	 pavector a);

/* ------------------------------------------------------------
 * Generalized minimal residual method (GMRES)
 * ------------------------------------------------------------ */

/** @brief Initialize GMRES.
 *
 *  The parameters are prepared for solving @f$A x = b@f$ with
 *  the GMRES method.
 *
 *  The maximal dimension of the Krylov subspace is determined
 *  by the number of columns of <tt>qr</tt>:
 *  for a <tt>k</tt>-dimensional subspace, <tt>qr->cols==k+1</tt>
 *  is required.
 *
 *  @param addeval Callback function representing the matrix @f$A@f$.
 *  @param matrix Data for the <tt>addeval</tt> callback.
 *  @param b Right-hand side vector @f$b@f$.
 *  @param x Initial guess for the solution @f$x@f$, will eventually
 *         be replaced by an improved approximation.
 *  @param rhat Transformed residual. The absolute value of
 *         <tt>rhat[*kk]</tt> is the Euclidean norm of the residual.
 *  @param q Next vector of the Krylov basis, constructed by
 *         Householder's elementary reflectors.
 *  @param kk Pointer to current dimension of the Krylov space.
 *  @param qr Representation of the Arnoldi basis @f$Q_{k+1}@f$ and the
 *         transformed matrix @f$Q_{k+1}^* A Q_k@f$.
 *  @param tau Scaling factors of elementary reflectors,
 *         provided by @ref qrdecomp_amatrix. */
HEADER_PREFIX void
init_gmres(addeval_t addeval,
	   void *matrix,
	   pcavector b, pavector x,
	   pavector rhat, pavector q,
	   uint *kk, pamatrix qr, pavector tau);

/** @brief One step of the GMRES method.
 *
 *  If <tt>*kk+1 >= qr->cols</tt>, there is no room for the next
 *  Arnoldi basis vector and the function returns immediately.
 *  It can be restarted using @ref finish_gmres.
 *
 *  Otherwise a new vector is added to the Arnoldi basis and the
 *  matrix @f$Q_{k+1}^* A Q_k@f$ is updated.
 *
 *  @remark This function currently makes no use of <tt>b</tt>
 *  and does not update <tt>x</tt>. The current residual can be
 *  tracked via <tt>rhat[*kk]</tt>. Once it is sufficiently small,
 *  the improved solution <tt>x</tt> can be obtained by using
 *  @ref finish_gmres.
 *
 *  @param addeval Callback function representing the matrix @f$A@f$.
 *  @param matrix Data for the <tt>addeval</tt> callback.
 *  @param b Right-hand side vector @f$b@f$.
 *  @param x Initial guess for the solution @f$x@f$, will eventually
 *         be replaced by an improved approximation.
 *  @param rhat Transformed residual. The absolute value of
 *         <tt>rhat[*kk]</tt> is the Euclidean norm of the residual.
 *  @param q Next vector of the Krylov basis, constructed by
 *         Householder's elementary reflectors.
 *  @param kk Pointer to current dimension of the Krylov space.
 *  @param qr Representation of the Arnoldi basis @f$Q_{k+1}@f$ and the
 *         transformed matrix @f$Q_{k+1}^* A Q_k@f$.
 *  @param tau Scaling factors of elementary reflectors,
 *         provided by @ref qrdecomp_amatrix. */
HEADER_PREFIX void
step_gmres(addeval_t addeval,
	   void *matrix,
	   pcavector b, pavector x,
	   pavector rhat, pavector q,
	   uint *kk, pamatrix qr, pavector tau);

/** @brief Completes or restarts the GMRES method.
 *
 *  Solves the least-squares problem
 *  @f$Q_{k+1}^* A Q_k \widehat{x} = Q_{k+1}^* r@f$ and performs the
 *  update @f$x \gets x + Q_k \widehat{x}@f$.
 *
 *  The function calls @ref init_gmres to reset the iteration and
 *  prepare for a restart.
 *  
 *  @param addeval Callback function representing the matrix @f$A@f$.
 *  @param matrix Data for the <tt>addeval</tt> callback.
 *  @param b Right-hand side vector @f$b@f$.
 *  @param x Initial guess for the solution @f$x@f$, will eventually
 *         be replaced by an improved approximation.
 *  @param rhat Transformed residual. The absolute value of
 *         <tt>rhat[*kk]</tt> is the Euclidean norm of the residual.
 *  @param q Next vector of the Krylov basis, constructed by
 *         Householder's elementary reflectors.
 *  @param kk Pointer to current dimension of the Krylov space.
 *  @param qr Representation of the Arnoldi basis @f$Q_{k+1}@f$ and the
 *         transformed matrix @f$Q_{k+1}^* A Q_k@f$.
 *  @param tau Scaling factors of elementary reflectors,
 *         provided by @ref qrdecomp_amatrix. */
HEADER_PREFIX void
finish_gmres(addeval_t addeval,
	     void *matrix,
	     pcavector b, pavector x,
	     pavector rhat, pavector q,
	     uint *kk, pamatrix qr, pavector tau);

/** @brief Returns norm of current residual vector.
 *
 *  @param rhat Transformed residual. The absolute value of
 *         <tt>rhat[k]</tt> is the Euclidean norm of the residual.
 *  @param k Current dimension of the Krylov space.
 *  @returns Norm of the residual. */
HEADER_PREFIX real
residualnorm_gmres(pcavector rhat, uint k);

/* ------------------------------------------------------------
 * Preconditioned generalized minimal residual method (GMRES)
 * ------------------------------------------------------------ */

/** @brief Initialize preconditioned GMRES.
 *
 *  The parameters are prepared for solving @f$N A x = N b@f$ with
 *  the GMRES method.
 *
 *  The maximal dimension of the Krylov subspace is determined
 *  by the number of columns of <tt>qr</tt>:
 *  for a <tt>k</tt>-dimensional subspace, <tt>qr->cols==k+1</tt>
 *  is required.
 *
 *  @param addeval Callback function representing the matrix @f$A@f$.
 *  @param matrix Data for the <tt>addeval</tt> callback.
 *  @param prcd Callback function representing the preconditioner @f$N@f$.
 *  @param pdata Data for the <tt>prcd</tt> callback.
 *  @param b Right-hand side vector @f$b@f$.
 *  @param x Initial guess for the solution @f$x@f$, will eventually
 *         be replaced by an improved approximation.
 *  @param rhat Transformed residual. The absolute value of
 *         <tt>rhat[*kk]</tt> is the Euclidean norm of the
 *         preconditioned residual.
 *  @param q Next vector of the Krylov basis, constructed by
 *         Householder's elementary reflectors.
 *  @param kk Pointer to current dimension of the Krylov space.
 *  @param qr Representation of the Arnoldi basis @f$Q_{k+1}@f$ and the
 *         transformed matrix @f$Q_{k+1}^* A Q_k@f$.
 *  @param tau Scaling factors of elementary reflectors,
 *         provided by @ref qrdecomp_amatrix. */
HEADER_PREFIX void
init_pgmres(addeval_t addeval,
	    void *matrix,
	    prcd_t prcd,
	    void *pdata,
	    pcavector b, pavector x,
	    pavector rhat, pavector q,
	    uint *kk, pamatrix qr, pavector tau);

/** @brief One step of the preconditioned GMRES method.
 *
 *  If <tt>*kk+1 >= qr->cols</tt>, there is no room for the next
 *  Arnoldi basis vector and the function returns immediately.
 *  It can be restarted using @ref finish_gmres.
 *
 *  Otherwise a new vector is added to the Arnoldi basis and the
 *  matrix @f$Q_{k+1}^* N A Q_k@f$ is updated.
 *
 *  @remark This function currently makes no use of <tt>b</tt>
 *  and does not update <tt>x</tt>. The current preconditioned residual
 *  can be tracked via <tt>rhat[*kk]</tt>. Once it is sufficiently small,
 *  the improved solution <tt>x</tt> can be obtained by using
 *  @ref finish_gmres.
 *
 *  @param addeval Callback function representing the matrix @f$A@f$.
 *  @param matrix Data for the <tt>addeval</tt> callback.
 *  @param prcd Callback function representing the preconditioner @f$N@f$.
 *  @param pdata Data for the <tt>prcd</tt> callback.
 *  @param b Right-hand side vector @f$b@f$.
 *  @param x Initial guess for the solution @f$x@f$, will eventually
 *         be replaced by an improved approximation.
 *  @param rhat Transformed residual. The absolute value of
 *         <tt>rhat[*kk]</tt> is the Euclidean norm of the
 *         preconditioned residual.
 *  @param q Next vector of the Krylov basis, constructed by
 *         Householder's elementary reflectors.
 *  @param kk Pointer to current dimension of the Krylov space.
 *  @param qr Representation of the Arnoldi basis @f$Q_{k+1}@f$ and the
 *         transformed matrix @f$Q_{k+1}^* A Q_k@f$.
 *  @param tau Scaling factors of elementary reflectors,
 *         provided by @ref qrdecomp_amatrix. */
HEADER_PREFIX void
step_pgmres(addeval_t addeval,
	    void *matrix,
	    prcd_t prcd,
	    void *pdata,
	    pcavector b, pavector x,
	    pavector rhat, pavector q,
	    uint *kk, pamatrix qr, pavector tau);

/** @brief Completes or restarts the preconditioned GMRES method.
 *
 *  Solves the least-squares problem
 *  @f$Q_{k+1}^* N A Q_k \widehat{x} = Q_{k+1}^* N r@f$ and performs the
 *  update @f$x \gets x + Q_k \widehat{x}@f$.
 *
 *  The function calls @ref init_pgmres to reset the iteration and
 *  prepare for a restart.
 *  
 *  @param addeval Callback function representing the matrix @f$A@f$.
 *  @param matrix Data for the <tt>addeval</tt> callback.
 *  @param prcd Callback function representing the preconditioner @f$N@f$.
 *  @param pdata Data for the <tt>prcd</tt> callback.
 *  @param b Right-hand side vector @f$b@f$.
 *  @param x Initial guess for the solution @f$x@f$, will eventually
 *         be replaced by an improved approximation.
 *  @param rhat Transformed residual. The absolute value of
 *         <tt>rhat[*kk]</tt> is the Euclidean norm of the
 *         preconditioned residual.
 *  @param q Next vector of the Krylov basis, constructed by
 *         Householder's elementary reflectors.
 *  @param kk Pointer to current dimension of the Krylov space.
 *  @param qr Representation of the Arnoldi basis @f$Q_{k+1}@f$ and the
 *         transformed matrix @f$Q_{k+1}^* A Q_k@f$.
 *  @param tau Scaling factors of elementary reflectors,
 *         provided by @ref qrdecomp_amatrix. */
HEADER_PREFIX void
finish_pgmres(addeval_t addeval,
	      void *matrix,
	      prcd_t prcd,
	      void *pdata,
	      pcavector b, pavector x,
	      pavector rhat, pavector q,
	      uint *kk, pamatrix qr, pavector tau);

/** @brief Returns norm of current preconditioned residual vector.
 *
 *  @param rhat Transformed residual. The absolute value of
 *         <tt>rhat[k]</tt> is the Euclidean norm of the
 *         preconditioned residual.
 *  @param k Current dimension of the Krylov space.
 *  @returns Norm of the residual. */
HEADER_PREFIX real
residualnorm_pgmres(pcavector rhat, uint k);

/** @} */

#endif
