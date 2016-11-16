/* ------------------------------------------------------------
 This is the file "laplacebem3d.h" of the H2Lib package.
 All rights reserved, Sven Christophersen 2011
 ------------------------------------------------------------ */

/**
 * @file laplacebem3d.h
 * @author Sven Christophersen
 * @date 2011
 */

#ifndef LAPLACEBEM3D_H_
#define LAPLACEBEM3D_H_

/* C STD LIBRARY */
/* CORE 0 */
#include "settings.h"
#include "parameters.h"
/* CORE 1 */
/* CORE 2 */
/* CORE 3 */
/* SIMPLE */
/* PARTICLES */
/* BEM */
#include "bem3d.h"

/** \defgroup laplacebem3d laplacebem3d
 *  @brief This module contains functions to setup and solve boundary integral
 *  equations for the Laplace operator in 3D.
 *  @{ */

/** @brief Constant that originates from the fundamental solution of the
 * Laplace equation. The value is @f$ 1 / 4 \pi @f$
 */
#define KERNEL_CONST_LAPLACEBEM3D 0.0795774715459476679

/* ------------------------------------------------------------
 Constructors and destructors
 ------------------------------------------------------------ */

/**
 * @brief Creates a new @ref _bem3d "bem3d"-object for computation of
 * single layer potential matrix of the Laplace equation.
 *
 * After calling this function the resulting @ref _bem3d "bem"-object will
 * provide all functionality that is necessary to build up fully populated
 * single layer potential matrix @f$ V \in \mathbb R^{\mathcal I \times
 * \mathcal I}@f$ and also @ref _hmatrix "hmatrix"
 * or @ref _h2matrix "h2matrix" approximation of this matrix.
 *
 * @param gr Surface mesh
 * @param q_regular Order of gaussian quadrature used within computation of matrix
 *        entries for single integrals and regular double integrals.
 * @param q_singular Order of gaussian quadrature used within computation of matrix
 *        entries singular double integrals.
 * @param row_basis Type of basis functions that are used for the test space.
 *        Can be one of the values defined in @ref basisfunctionbem3d.
 * @param col_basis Type of basis functions that are used for the trial space.
 *        Can be one of the values defined in @ref basisfunctionbem3d.
 *
 * @return Returns a @ref _bem3d "bem"-object that can compute fully populated
 * slp matrices @f$ V @f$ for the Laplace equation.
 */
HEADER_PREFIX pbem3d new_slp_laplace_bem3d(pcsurface3d gr, uint q_regular,
    uint q_singular, basisfunctionbem3d row_basis, basisfunctionbem3d col_basis);

/**
 * @brief Creates a new @ref _bem3d "bem3d"-object for computation of
 * double layer potential matrix plus a scalar times the mass matrix
 * of the Laplace equation.
 *
 * After calling this function the resulting @ref _bem3d "bem"-object will
 * provide all functionality that is necessary to build up fully populated
 * double layer potential matrix @f$ K + \frac{1}{2} M \in \mathbb
 * R^{\mathcal I \times \mathcal J}@f$ and also @ref _hmatrix "hmatrix"
 * or @ref _h2matrix "h2matrix"
 * approximation of this matrix.
 *
 * @param gr Surface mesh.
 * @param q_regular Order of gaussian quadrature used within computation of matrix
 *        entries for single integrals and regular double integrals.
 * @param q_singular Order of gaussian quadrature used within computation of matrix
 *        entries singular double integrals.
 * @param row_basis Type of basis functions that are used for the test space.
 *        Can be one of the values defined in @ref basisfunctionbem3d.
 * @param col_basis Type of basis functions that are used for the trial space.
 *        Can be one of the values defined in @ref basisfunctionbem3d.
 * @param alpha Double layer operator + @f$\alpha@f$ mass matrix.
 *
 * @return Returns a @ref _bem3d "bem"-object that can compute fully populated
 * dlp matrices @f$ K + \frac{1}{2} M @f$ for the Laplace equation.
 */
HEADER_PREFIX pbem3d new_dlp_laplace_bem3d(pcsurface3d gr, uint q_regular,
    uint q_singular, basisfunctionbem3d row_basis, basisfunctionbem3d col_basis,
    field alpha);

/**
 * @brief Creates a new @ref _bem3d "bem3d"-object for computation of
 * adjoint double layer potential matrix plus a scalar times the mass matrix
 * of the Helmholtz equation.
 *
 * After calling this function the resulting @ref _bem3d "bem"-object will
 * provide all functionality that is necessary to build up fully populated
 * adjoint double layer potential matrix @f$ K' + \alpha M \in \mathbb
 * R^{\mathcal I \times \mathcal J}@f$ and also @ref _hmatrix "hmatrix"
 * or @ref _h2matrix "h2matrix"
 * approximation of this matrix.
 *
 * @param gr Surface mesh.
 * @param q_regular Order of gaussian quadrature used within computation of matrix
 *        entries for single integrals and regular double integrals.
 * @param q_singular Order of gaussian quadrature used within computation of matrix
 *        entries singular double integrals.
 * @param row_basis Type of basis functions that are used for the test space.
 *        Can be one of the values defined in @ref basisfunctionbem3d.
 * @param col_basis Type of basis functions that are used for the trial space.
 *        Can be one of the values defined in @ref basisfunctionbem3d.
 * @param alpha Adjoint double layer operator + @f$\alpha@f$ mass matrix.
 *
 * @return Returns a @ref _bem3d "bem"-object that can compute fully populated
 * adlp matrices @f$ K' + \alpha M @f$ for the Laplace equation.
 */
HEADER_PREFIX pbem3d
new_adlp_laplace_bem3d(pcsurface3d gr, uint q_regular, uint q_singular,
    basisfunctionbem3d row_basis, basisfunctionbem3d col_basis, field alpha);

/**
 * @brief Delete a @ref _bem3d "bem3d" object for the Laplace equation
 *
 * @param bem Object to be deleted.
 */
HEADER_PREFIX void del_laplace_bem3d(pbem3d bem);

/* ------------------------------------------------------------
 Examples for Dirichlet- / Neumann-data to test linear system
 ------------------------------------------------------------ */

/**
 * @brief A simple linear harmonic function that will serve as dirichlet values.
 *
 * When computing the neumann data out of the dirichlet data one can use this
 * function as test data which will generate dirichlet values of with the following
 * values:
 * @f[
 * f(\vec x, \, \vec n) = x_1 + x_2 + x_3
 * @f]
 * Corresponding neumann data can be generated by using
 * @ref eval_neumann_linear_laplacebem3d.
 * <br>
 * To build up an appropriate dirichlet data coefficient vector one needs the
 * @f$ L_2@f$-projection. This can be done by passing this function to
 * @ref projectL2_bem3d_c_avector for piecewise constant basis functions or
 * to @ref projectL2_bem3d_l_avector for piecewise linear basis functions.
 *
 * @param x Evaluation point.
 * @param n Normal vector to current evaluation point.
 * @param data Additional data for evaluating the functional
 * @return returns the function value of @f$ f(\vec x, \, \vec n) @f$.
 */
HEADER_PREFIX field eval_dirichlet_linear_laplacebem3d(const real *x,
    const real *n, void *data);

/**
 * @brief A simple linear harmonic function that will serve as neumann values.
 *
 * When computing the neumann data out of the dirichlet data one can use this
 * function as test data which will generate neumann values of with the following
 * values:
 * @f[
 * f(\vec x, \, \vec n) = n_1 + n_2 + n_3
 * @f]
 * Corresponding dirichlet data can be generated by using
 * @ref eval_dirichlet_linear_laplacebem3d.
 * <br>
 * To build up an appropriate neumann data coefficient vector one needs the
 * @f$ L_2@f$-projection. This can be done by passing this function to
 * @ref projectL2_bem3d_c_avector for piecewise constant basis functions or
 * to @ref projectL2_bem3d_l_avector for piecewise linear basis functions.
 *
 * @param x Evaluation point.
 * @param n Normal vector to current evaluation point.
 * @param data Additional data for evaluating the functional
 * @return returns the function value of @f$ f(\vec x, \, \vec n) @f$.
 */
HEADER_PREFIX field eval_neumann_linear_laplacebem3d(const real *x,
    const real *n, void *data);

/**
 * @brief A simple quadratic harmonic function that will serve as dirichlet values.
 *
 * When computing the neumann data out of the dirichlet data one can use this
 * function as test data which will generate dirichlet values of with the following
 * values:
 * @f[
 * f(\vec x, \, \vec n) = x_1^2 - 0.5 \, x_2^2 - 0.5 \, x_3^2
 * @f]
 * Corresponding neumann data can be generated by using
 * @ref eval_neumann_quadratic_laplacebem3d.
 * <br>
 * To build up an appropriate dirichlet data coefficient vector one needs the
 * @f$ L_2@f$-projection. This can be done by passing this function to
 * @ref projectL2_bem3d_c_avector for piecewise constant basis functions or
 * to @ref projectL2_bem3d_l_avector for piecewise linear basis functions.
 *
 * @param x Evaluation point.
 * @param n Normal vector to current evaluation point.
 * @param data Additional data for evaluating the functional
 * @return returns the function value of @f$ f(\vec x, \, \vec n) @f$.
 */
HEADER_PREFIX field eval_dirichlet_quadratic_laplacebem3d(const real *x,
    const real *n, void *data);

/**
 * @brief A simple quadratic harmonic function that will serve as neumann values.
 *
 * When computing the neumann data out of the dirichlet data one can use this
 * function as test data which will generate neumann values of with the following
 * values:
 * @f[
 * f(\vec x, \, \vec n) = 2 \, \left( n_1 \, x_1 - 0.5 \, n_2  \, x_2 -
 * 0.5 \, n_3 \, x_3 \right)
 * @f]
 * Corresponding dirichlet data can be generated by using
 * @ref eval_dirichlet_quadratic_laplacebem3d.
 * <br>
 * To build up an appropriate neumann data coefficient vector one needs the
 * @f$ L_2@f$-projection. This can be done by passing this function to
 * @ref projectL2_bem3d_c_avector for piecewise constant basis functions or
 * to @ref projectL2_bem3d_l_avector for piecewise linear basis functions.
 *
 * @param x Evaluation point.
 * @param n Normal vector to current evaluation point.
 * @param data Additional data for evaluating the functional
 * @return returns the function value of @f$ f(\vec x, \, \vec n) @f$.
 */
HEADER_PREFIX field eval_neumann_quadratic_laplacebem3d(const real *x,
    const real *n, void *data);

/**
 * @brief A harmonic function based upon the fundamental solution,
 * that will serve as dirichlet values.
 *
 * When computing the neumann data out of the dirichlet data one can use this
 * function as test data which will generate dirichlet values of with the following
 * values:
 * @f[
 * f(\vec x, \, \vec n) = \frac{1}{\lVert \vec x - \vec z \rVert}
 * @f]
 * with @f$ \vec z = \left( 1.2, 1.2, 1.2 \right)^* @f$.
 * Corresponding neumann data can be generated by using
 * @ref eval_neumann_fundamental_laplacebem3d.
 * <br>
 * To build up an appropriate dirichlet data coefficient vector one needs the
 * @f$ L_2@f$-projection. This can be done by passing this function to
 * @ref projectL2_bem3d_c_avector for piecewise constant basis functions or
 * to @ref projectL2_bem3d_l_avector for piecewise linear basis functions.
 *
 * @param x Evaluation point.
 * @param n Normal vector to current evaluation point.
 * @param data Additional data for evaluating the functional
 * @return returns the function value of @f$ f(\vec x, \, \vec n) @f$.
 */
HEADER_PREFIX field eval_dirichlet_fundamental_laplacebem3d(const real *x,
    const real *n, void *data);
/**
 * @brief A harmonic function based upon the fundamental solution,
 * that will serve as neumann values.
 *
 * When computing the neumann data out of the dirichlet data one can use this
 * function as test data which will generate neumann values of with the following
 * values:
 * @f[
 * f(\vec x, \, \vec n) = -\frac{<\vec n, \vec x - \vec z>}{\lVert \vec x - \vec z \rVert^3}
 * @f]
 * with @f$ \vec z = \left( 1.2, 1.2, 1.2 \right)^* @f$.
 * Corresponding dirichlet data can be generated by using
 * @ref eval_dirichlet_fundamental_laplacebem3d.
 * <br>
 * To build up an appropriate neumann data coefficient vector one needs the
 * @f$ L_2@f$-projection. This can be done by passing this function to
 * @ref projectL2_bem3d_c_avector for piecewise constant basis functions or
 * to @ref projectL2_bem3d_l_avector for piecewise linear basis functions.
 *
 * @param x Evaluation point.
 * @param n Normal vector to current evaluation point.
 * @param data Additional data for evaluating the functional
 * @return returns the function value of @f$ f(\vec x, \, \vec n) @f$.
 */
HEADER_PREFIX field eval_neumann_fundamental_laplacebem3d(const real *x,
    const real *n, void *data);

/**
 * @brief A harmonic function based upon the fundamental solution,
 * that will serve as dirichlet values.
 *
 * When computing the neumann data out of the dirichlet data one can use this
 * function as test data which will generate dirichlet values of with the following
 * values:
 * @f[
 * f(\vec x, \, \vec n) = \frac{1}{\lVert \vec x - \vec z \rVert}
 * @f]
 * with @f$ \vec z = \left( 0.25, 0.25, 0.25 \right)^* @f$.
 * Corresponding neumann data can be generated by using
 * @ref eval_neumann_fundamental2_laplacebem3d.
 * <br>
 * To build up an appropriate dirichlet data coefficient vector one needs the
 * @f$ L_2@f$-projection. This can be done by passing this function to
 * @ref projectL2_bem3d_c_avector for piecewise constant basis functions or
 * to @ref projectL2_bem3d_l_avector for piecewise linear basis functions.
 *
 * @param x Evaluation point.
 * @param n Normal vector to current evaluation point.
 * @param data Additional data for evaluating the functional
 * @return returns the function value of @f$ f(\vec x, \, \vec n) @f$.
 */
HEADER_PREFIX field eval_dirichlet_fundamental2_laplacebem3d(const real *x,
    const real *n, void *data);
/**
 * @brief A harmonic function based upon the fundamental solution,
 * that will serve as neumann values.
 *
 * When computing the neumann data out of the dirichlet data one can use this
 * function as test data which will generate neumann values of with the following
 * values:
 * @f[
 * f(\vec x, \, \vec n) = -\frac{<\vec n, \vec x - \vec z>}{\lVert \vec x - \vec z \rVert^3}
 * @f]
 * with @f$ \vec z = \left( 0.25, 0.25, 0.25 \right)^* @f$.
 * Corresponding dirichlet data can be generated by using
 * @ref eval_dirichlet_fundamental_laplacebem3d.
 * <br>
 * To build up an appropriate neumann data coefficient vector one needs the
 * @f$ L_2@f$-projection. This can be done by passing this function to
 * @ref projectL2_bem3d_c_avector for piecewise constant basis functions or
 * to @ref projectL2_bem3d_l_avector for piecewise linear basis functions.
 *
 * @param x Evaluation point.
 * @param n Normal vector to current evaluation point.
 * @param data Additional data for evaluating the functional
 * @return returns the function value of @f$ f(\vec x, \, \vec n) @f$.
 */
HEADER_PREFIX field eval_neumann_fundamental2_laplacebem3d(const real *x,
    const real *n, void *data);

/** @} */

#endif /* LAPLACEBEM3D_H_ */
