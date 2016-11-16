/* ------------------------------------------------------------
 * This is the file "krylovsolvers.h" of the H2Lib package.
 * All rights reserved, Steffen Boerm 2016
 * ------------------------------------------------------------ */

/** @file krylovsolvers.h
 *  @author Steffen B&ouml;rm */

#ifndef KRYLOVSOLVERS_H
#define KRYLOVSOLVERS_H

#include "amatrix.h"
#include "sparsematrix.h"
#include "hmatrix.h"
#include "h2matrix.h"
#include "dh2matrix.h"
#include "krylov.h"

/** @defgroup krylovsolvers krylovsolvers
 *  @brief Convenience functions for solving linear systems by
 *  Krylov methods.
 *  @{ */

/** @brief Solve a self-adjoint positive definite system @f$Ax=b@f$
 *  with the conjugate gradient method and a general matrix type <tt>A</tt>.
 *
 *  @param A System matrix, has to be self-adjoint and positive definite.
 *  @param addeval_A General callback function for evaluation of a matrix
 *         <tt>A</tt>.
 *  @param b Right-hand side vector.
 *  @param x Initial guess, will be overwritten by approximate solution.
 *  @param eps Relative accuracy @f$\epsilon@f$, the method stops if
 *         @f$\|Ax-b\|_2 \leq \epsilon \|b\|_2@f$.
 *  @param maxiter Maximal number of iterations. <tt>maxiter=0</tt>
 *         means that the number of iterations is not bounded.
 *  @returns Number of iterations. */
HEADER_PREFIX uint
solve_cg_avector(void* A, addeval_t addeval_A, pcavector b, pavector x,
    real eps, uint maxiter);

/** @brief Solve a self-adjoint positive definite system @f$Ax=b@f$
 *  with the conjugate gradient method.
 *
 *  @param A System matrix, has to be self-adjoint and positive definite.
 *  @param b Right-hand side vector.
 *  @param x Initial guess, will be overwritten by approximate solution.
 *  @param eps Relative accuracy @f$\epsilon@f$, the method stops if
 *         @f$\|Ax-b\|_2 \leq \epsilon \|b\|_2@f$.
 *  @param maxiter Maximal number of iterations. <tt>maxiter=0</tt>
 *         means that the number of iterations is not bounded.
 *  @returns Number of iterations. */
HEADER_PREFIX uint
solve_cg_amatrix_avector(pcamatrix A, pcavector b, pavector x, real eps,
    uint maxiter);

/** @brief Solve a self-adjoint positive definite system @f$Ax=b@f$
 *  with the conjugate gradient method.
 *
 *  @param A System matrix, has to be self-adjoint and positive definite.
 *  @param b Right-hand side vector.
 *  @param x Initial guess, will be overwritten by approximate solution.
 *  @param eps Relative accuracy @f$\epsilon@f$, the method stops if
 *         @f$\|Ax-b\|_2 \leq \epsilon \|b\|_2@f$.
 *  @param maxiter Maximal number of iterations. <tt>maxiter=0</tt>
 *         means that the number of iterations is not bounded.
 *  @returns Number of iterations. */
HEADER_PREFIX uint
solve_cg_sparsematrix_avector(pcsparsematrix A, pcavector b, pavector x,
    real eps, uint maxiter);

/** @brief Solve a self-adjoint positive definite system @f$Ax=b@f$
 *  with the conjugate gradient method.
 *
 *  @param A System matrix, has to be self-adjoint and positive definite.
 *  @param b Right-hand side vector.
 *  @param x Initial guess, will be overwritten by approximate solution.
 *  @param eps Relative accuracy @f$\epsilon@f$, the method stops if
 *         @f$\|Ax-b\|_2 \leq \epsilon \|b\|_2@f$.
 *  @param maxiter Maximal number of iterations. <tt>maxiter=0</tt>
 *         means that the number of iterations is not bounded.
 *  @returns Number of iterations. */
HEADER_PREFIX uint
solve_cg_hmatrix_avector(pchmatrix A, pcavector b, pavector x, real eps,
    uint maxiter);

/** @brief Solve a self-adjoint positive definite system @f$Ax=b@f$
 *  with the conjugate gradient method.
 *
 *  @param A System matrix, has to be self-adjoint and positive definite.
 *  @param b Right-hand side vector.
 *  @param x Initial guess, will be overwritten by approximate solution.
 *  @param eps Relative accuracy @f$\epsilon@f$, the method stops if
 *         @f$\|Ax-b\|_2 \leq \epsilon \|b\|_2@f$.
 *  @param maxiter Maximal number of iterations. <tt>maxiter=0</tt>
 *         means that the number of iterations is not bounded.
 *  @returns Number of iterations. */
HEADER_PREFIX uint
solve_cg_h2matrix_avector(pch2matrix A, pcavector b, pavector x, real eps,
    uint maxiter);

/** @brief Solve a self-adjoint positive definite system @f$Ax=b@f$
 *  with the conjugate gradient method.
 *
 *  @param A System matrix, has to be self-adjoint and positive definite.
 *  @param b Right-hand side vector.
 *  @param x Initial guess, will be overwritten by approximate solution.
 *  @param eps Relative accuracy @f$\epsilon@f$, the method stops if
 *         @f$\|Ax-b\|_2 \leq \epsilon \|b\|_2@f$.
 *  @param maxiter Maximal number of iterations. <tt>maxiter=0</tt>
 *         means that the number of iterations is not bounded.
 *  @returns Number of iterations. */
HEADER_PREFIX uint
solve_cg_dh2matrix_avector(pcdh2matrix A, pcavector b, pavector x, real eps,
    uint maxiter);

/** @brief Solve a self-adjoint positive definite system @f$Ax=b@f$
 *  with the preconditioned conjugate gradient method.
 *
 *  @param A System matrix, has to be self-adjoint and positive definite.
 *  @param addeval_A General callback function for evaluation of a matrix
 *         <tt>A</tt>.
 *  @param prcd Callback function for preconditioner.
 *  @param pdata Data for <tt>prcd</tt> callback function.
 *  @param b Right-hand side vector.
 *  @param x Initial guess, will be overwritten by approximate solution.
 *  @param eps Relative accuracy @f$\epsilon@f$, the method stops if
 *         @f$\|Ax-b\|_2 \leq \epsilon \|b\|_2@f$.
 *  @param maxiter Maximal number of iterations. <tt>maxiter=0</tt>
 *         means that the number of iterations is not bounded.
 *  @returns Number of iterations. */
HEADER_PREFIX uint
solve_pcg_avector(void *A, addeval_t addeval_A, prcd_t prcd, void *pdata,
    pcavector b, pavector x, real eps, uint maxiter);

/** @brief Solve a self-adjoint positive definite system @f$Ax=b@f$
 *  with the preconditioned conjugate gradient method.
 *
 *  @param A System matrix, has to be self-adjoint and positive definite.
 *  @param prcd Callback function for preconditioner.
 *  @param pdata Data for <tt>prcd</tt> callback function.
 *  @param b Right-hand side vector.
 *  @param x Initial guess, will be overwritten by approximate solution.
 *  @param eps Relative accuracy @f$\epsilon@f$, the method stops if
 *         @f$\|Ax-b\|_2 \leq \epsilon \|b\|_2@f$.
 *  @param maxiter Maximal number of iterations. <tt>maxiter=0</tt>
 *         means that the number of iterations is not bounded.
 *  @returns Number of iterations. */
HEADER_PREFIX uint
solve_pcg_amatrix_avector(pcamatrix A, prcd_t prcd, void *pdata, pcavector b,
    pavector x, real eps, uint maxiter);

/** @brief Solve a self-adjoint positive definite system @f$Ax=b@f$
 *  with the preconditioned conjugate gradient method.
 *
 *  @param A System matrix, has to be self-adjoint and positive definite.
 *  @param prcd Callback function for preconditioner.
 *  @param pdata Data for <tt>prcd</tt> callback function.
 *  @param b Right-hand side vector.
 *  @param x Initial guess, will be overwritten by approximate solution.
 *  @param eps Relative accuracy @f$\epsilon@f$, the method stops if
 *         @f$\|Ax-b\|_2 \leq \epsilon \|b\|_2@f$.
 *  @param maxiter Maximal number of iterations. <tt>maxiter=0</tt>
 *         means that the number of iterations is not bounded.
 *  @returns Number of iterations. */
HEADER_PREFIX uint
solve_pcg_sparsematrix_avector(pcsparsematrix A, prcd_t prcd, void *pdata,
    pcavector b, pavector x, real eps, uint maxiter);

/** @brief Solve a self-adjoint positive definite system @f$Ax=b@f$
 *  with the preconditioned conjugate gradient method.
 *
 *  @param A System matrix, has to be self-adjoint and positive definite.
 *  @param prcd Callback function for preconditioner.
 *  @param pdata Data for <tt>prcd</tt> callback function.
 *  @param b Right-hand side vector.
 *  @param x Initial guess, will be overwritten by approximate solution.
 *  @param eps Relative accuracy @f$\epsilon@f$, the method stops if
 *         @f$\|Ax-b\|_2 \leq \epsilon \|b\|_2@f$.
 *  @param maxiter Maximal number of iterations. <tt>maxiter=0</tt>
 *         means that the number of iterations is not bounded.
 *  @returns Number of iterations. */
HEADER_PREFIX uint
solve_pcg_hmatrix_avector(pchmatrix A, prcd_t prcd, void *pdata, pcavector b,
    pavector x, real eps, uint maxiter);

/** @brief Solve a self-adjoint positive definite system @f$Ax=b@f$
 *  with the preconditioned conjugate gradient method.
 *
 *  @param A System matrix, has to be self-adjoint and positive definite.
 *  @param prcd Callback function for preconditioner.
 *  @param pdata Data for <tt>prcd</tt> callback function.
 *  @param b Right-hand side vector.
 *  @param x Initial guess, will be overwritten by approximate solution.
 *  @param eps Relative accuracy @f$\epsilon@f$, the method stops if
 *         @f$\|Ax-b\|_2 \leq \epsilon \|b\|_2@f$.
 *  @param maxiter Maximal number of iterations. <tt>maxiter=0</tt>
 *         means that the number of iterations is not bounded.
 *  @returns Number of iterations. */
HEADER_PREFIX uint
solve_pcg_h2matrix_avector(pch2matrix A, prcd_t prcd, void *pdata, pcavector b,
    pavector x, real eps, uint maxiter);

/** @brief Solve a self-adjoint positive definite system @f$Ax=b@f$
 *  with the preconditioned conjugate gradient method.
 *
 *  @param A System matrix, has to be self-adjoint and positive definite.
 *  @param prcd Callback function for preconditioner.
 *  @param pdata Data for <tt>prcd</tt> callback function.
 *  @param b Right-hand side vector.
 *  @param x Initial guess, will be overwritten by approximate solution.
 *  @param eps Relative accuracy @f$\epsilon@f$, the method stops if
 *         @f$\|Ax-b\|_2 \leq \epsilon \|b\|_2@f$.
 *  @param maxiter Maximal number of iterations. <tt>maxiter=0</tt>
 *         means that the number of iterations is not bounded.
 *  @returns Number of iterations. */
HEADER_PREFIX uint
solve_pcg_dh2matrix_avector(pcdh2matrix A, prcd_t prcd, void *pdata,
    pcavector b, pavector x, real eps, uint maxiter);

/** @brief Solve a linear system @f$Ax=b@f$ with the
 *  generalized minimal residual method.
 *
 *  @param A System matrix, should be invertible.
 *  @param addeval_A General callback function for evaluation of a matrix
 *         <tt>A</tt>.
 *  @param b Right-hand side vector.
 *  @param x Initial guess, will be overwritten by approximate solution.
 *  @param eps Relative accuracy @f$\epsilon@f$, the method stops if
 *         @f$\|Ax-b\|_2 \leq \epsilon \|b\|_2@f$.
 *  @param maxiter Maximal number of iterations. <tt>maxiter=0</tt>
 *         means that the number of iterations is not bounded.
 *  @param kmax Maximal dimension of Krylov subspace.
 *  @returns Number of iterations. */
HEADER_PREFIX uint
solve_gmres_avector(void* A, addeval_t addeval_A, pcavector b, pavector x,
    real eps, uint maxiter, uint kmax);

/** @brief Solve a linear system @f$Ax=b@f$ with the
 *  generalized minimal residual method.
 *
 *  @param A System matrix, should be invertible.
 *  @param b Right-hand side vector.
 *  @param x Initial guess, will be overwritten by approximate solution.
 *  @param eps Relative accuracy @f$\epsilon@f$, the method stops if
 *         @f$\|Ax-b\|_2 \leq \epsilon \|b\|_2@f$.
 *  @param maxiter Maximal number of iterations. <tt>maxiter=0</tt>
 *         means that the number of iterations is not bounded.
 *  @param kmax Maximal dimension of Krylov subspace.
 *  @returns Number of iterations. */
HEADER_PREFIX uint
solve_gmres_amatrix_avector(pcamatrix A, pcavector b, pavector x, real eps,
    uint maxiter, uint kmax);

/** @brief Solve a linear system @f$Ax=b@f$ with the
 *  generalized minimal residual method.
 *
 *  @param A System matrix, should be invertible.
 *  @param b Right-hand side vector.
 *  @param x Initial guess, will be overwritten by approximate solution.
 *  @param eps Relative accuracy @f$\epsilon@f$, the method stops if
 *         @f$\|Ax-b\|_2 \leq \epsilon \|b\|_2@f$.
 *  @param maxiter Maximal number of iterations. <tt>maxiter=0</tt>
 *         means that the number of iterations is not bounded.
 *  @param kmax Maximal dimension of Krylov subspace.
 *  @returns Number of iterations. */
HEADER_PREFIX uint
solve_gmres_sparsematrix_avector(pcsparsematrix A, pcavector b, pavector x,
    real eps, uint maxiter, uint kmax);

/** @brief Solve a linear system @f$Ax=b@f$ with the
 *  generalized minimal residual method.
 *
 *  @param A System matrix, should be invertible.
 *  @param b Right-hand side vector.
 *  @param x Initial guess, will be overwritten by approximate solution.
 *  @param eps Relative accuracy @f$\epsilon@f$, the method stops if
 *         @f$\|Ax-b\|_2 \leq \epsilon \|b\|_2@f$.
 *  @param maxiter Maximal number of iterations. <tt>maxiter=0</tt>
 *         means that the number of iterations is not bounded.
 *  @param kmax Maximal dimension of Krylov subspace.
 *  @returns Number of iterations. */
HEADER_PREFIX uint
solve_gmres_hmatrix_avector(pchmatrix A, pcavector b, pavector x, real eps,
    uint maxiter, uint kmax);

/** @brief Solve a linear system @f$Ax=b@f$ with the
 *  generalized minimal residual method.
 *
 *  @param A System matrix, should be invertible.
 *  @param b Right-hand side vector.
 *  @param x Initial guess, will be overwritten by approximate solution.
 *  @param eps Relative accuracy @f$\epsilon@f$, the method stops if
 *         @f$\|Ax-b\|_2 \leq \epsilon \|b\|_2@f$.
 *  @param maxiter Maximal number of iterations. <tt>maxiter=0</tt>
 *         means that the number of iterations is not bounded.
 *  @param kmax Maximal dimension of Krylov subspace.
 *  @returns Number of iterations. */
HEADER_PREFIX uint
solve_gmres_h2matrix_avector(pch2matrix A, pcavector b, pavector x, real eps,
    uint maxiter, uint kmax);

/** @brief Solve a linear system @f$Ax=b@f$ with the
 *  generalized minimal residual method.
 *
 *  @param A System matrix, should be invertible.
 *  @param b Right-hand side vector.
 *  @param x Initial guess, will be overwritten by approximate solution.
 *  @param eps Relative accuracy @f$\epsilon@f$, the method stops if
 *         @f$\|Ax-b\|_2 \leq \epsilon \|b\|_2@f$.
 *  @param maxiter Maximal number of iterations. <tt>maxiter=0</tt>
 *         means that the number of iterations is not bounded.
 *  @param kmax Maximal dimension of Krylov subspace.
 *  @returns Number of iterations. */
HEADER_PREFIX uint
solve_gmres_dh2matrix_avector(pcdh2matrix A, pcavector b, pavector x, real eps,
    uint maxiter, uint kmax);

/** @brief Solve a linear system @f$Ax=b@f$ with the
 *  preconditioned generalized minimal residual method.
 *
 *  @param A System matrix, should be invertible.
 *  @param addeval_A General callback function for evaluation of a matrix
 *         <tt>A</tt>.
 *  @param prcd Callback function for preconditioner @f$N@f$.
 *  @param pdata Data for <tt>prcd</tt> callback function.
 *  @param b Right-hand side vector.
 *  @param x Initial guess, will be overwritten by approximate solution.
 *  @param eps Relative accuracy @f$\epsilon@f$, the method stops if
 *         @f$\|N(Ax-b)\|_2 \leq \epsilon \|N b\|_2@f$.
 *  @param maxiter Maximal number of iterations. <tt>maxiter=0</tt>
 *         means that the number of iterations is not bounded.
 *  @param kmax Maximal dimension of Krylov subspace.
 *  @returns Number of iterations. */
uint
solve_pgmres_avector(void* A, addeval_t addeval_A, prcd_t prcd, void *pdata,
    pcavector b, pavector x, real eps, uint maxiter, uint kmax);

/** @brief Solve a linear system @f$Ax=b@f$ with the
 *  preconditioned generalized minimal residual method.
 *
 *  @param A System matrix, should be invertible.
 *  @param prcd Callback function for preconditioner @f$N@f$.
 *  @param pdata Data for <tt>prcd</tt> callback function.
 *  @param b Right-hand side vector.
 *  @param x Initial guess, will be overwritten by approximate solution.
 *  @param eps Relative accuracy @f$\epsilon@f$, the method stops if
 *         @f$\|N(Ax-b)\|_2 \leq \epsilon \|N b\|_2@f$.
 *  @param maxiter Maximal number of iterations. <tt>maxiter=0</tt>
 *         means that the number of iterations is not bounded.
 *  @param kmax Maximal dimension of Krylov subspace.
 *  @returns Number of iterations. */
uint
solve_pgmres_amatrix_avector(pcamatrix A, prcd_t prcd, void *pdata, pcavector b,
    pavector x, real eps, uint maxiter, uint kmax);

/** @brief Solve a linear system @f$Ax=b@f$ with the
 *  preconditioned generalized minimal residual method.
 *
 *  @param A System matrix, should be invertible.
 *  @param prcd Callback function for preconditioner @f$N@f$.
 *  @param pdata Data for <tt>prcd</tt> callback function.
 *  @param b Right-hand side vector.
 *  @param x Initial guess, will be overwritten by approximate solution.
 *  @param eps Relative accuracy @f$\epsilon@f$, the method stops if
 *         @f$\|N(Ax-b)\|_2 \leq \epsilon \|N b\|_2@f$.
 *  @param maxiter Maximal number of iterations. <tt>maxiter=0</tt>
 *         means that the number of iterations is not bounded.
 *  @param kmax Maximal dimension of Krylov subspace.
 *  @returns Number of iterations. */
uint
solve_pgmres_sparsematrix_avector(pcsparsematrix A, prcd_t prcd, void *pdata,
    pcavector b, pavector x, real eps, uint maxiter, uint kmax);

/** @brief Solve a linear system @f$Ax=b@f$ with the
 *  preconditioned generalized minimal residual method.
 *
 *  @param A System matrix, should be invertible.
 *  @param prcd Callback function for preconditioner @f$N@f$.
 *  @param pdata Data for <tt>prcd</tt> callback function.
 *  @param b Right-hand side vector.
 *  @param x Initial guess, will be overwritten by approximate solution.
 *  @param eps Relative accuracy @f$\epsilon@f$, the method stops if
 *         @f$\|N(Ax-b)\|_2 \leq \epsilon \|N b\|_2@f$.
 *  @param maxiter Maximal number of iterations. <tt>maxiter=0</tt>
 *         means that the number of iterations is not bounded.
 *  @param kmax Maximal dimension of Krylov subspace.
 *  @returns Number of iterations. */
uint
solve_pgmres_hmatrix_avector(pchmatrix A, prcd_t prcd, void *pdata, pcavector b,
    pavector x, real eps, uint maxiter, uint kmax);

/** @brief Solve a linear system @f$Ax=b@f$ with the
 *  preconditioned generalized minimal residual method.
 *
 *  @param A System matrix, should be invertible.
 *  @param prcd Callback function for preconditioner @f$N@f$.
 *  @param pdata Data for <tt>prcd</tt> callback function.
 *  @param b Right-hand side vector.
 *  @param x Initial guess, will be overwritten by approximate solution.
 *  @param eps Relative accuracy @f$\epsilon@f$, the method stops if
 *         @f$\|N(Ax-b)\|_2 \leq \epsilon \|N b\|_2@f$.
 *  @param maxiter Maximal number of iterations. <tt>maxiter=0</tt>
 *         means that the number of iterations is not bounded.
 *  @param kmax Maximal dimension of Krylov subspace.
 *  @returns Number of iterations. */
uint
solve_pgmres_h2matrix_avector(pch2matrix A, prcd_t prcd, void *pdata,
    pcavector b, pavector x, real eps, uint maxiter, uint kmax);

/** @brief Solve a linear system @f$Ax=b@f$ with the
 *  preconditioned generalized minimal residual method.
 *
 *  @param A System matrix, should be invertible.
 *  @param prcd Callback function for preconditioner @f$N@f$.
 *  @param pdata Data for <tt>prcd</tt> callback function.
 *  @param b Right-hand side vector.
 *  @param x Initial guess, will be overwritten by approximate solution.
 *  @param eps Relative accuracy @f$\epsilon@f$, the method stops if
 *         @f$\|N(Ax-b)\|_2 \leq \epsilon \|N b\|_2@f$.
 *  @param maxiter Maximal number of iterations. <tt>maxiter=0</tt>
 *         means that the number of iterations is not bounded.
 *  @param kmax Maximal dimension of Krylov subspace.
 *  @returns Number of iterations. */
uint
solve_pgmres_dh2matrix_avector(pcdh2matrix A, prcd_t prcd, void *pdata,
    pcavector b, pavector x, real eps, uint maxiter, uint kmax);

/** @} */

#endif
