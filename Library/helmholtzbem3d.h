/* ------------------------------------------------------------
 * This is the file "helmholtzbem3d.h" of the H2Lib package.
 * All rights reserved, Sven Christophersen 2015
 * ------------------------------------------------------------ */

/**
 * @file helmholtzbem3d.h
 * @author Sven Christophersen
 * @date 2015
 */

#ifndef HELMHOLTZBEM3D_H_
#define HELMHOLTZBEM3D_H_

/* C STD LIBRARY */
/* CORE 0 */
/* CORE 1 */
/* CORE 2 */
/* CORE 3 */
/* SIMPLE */
/* PARTICLES */
/* BEM */
#include "bem3d.h"

#ifdef USE_COMPLEX

/** \defgroup helmholtzbem3d helmholtzbem3d
 *  @brief This module contains functions to setup and solve boundary integral
 *  equations for the Helmholtz problem in 3D.
 *  @{ */

/** @brief Constant that originates from the fundamental solution of the
 * Helmholtz equation. The value is @f$ 1 / 4 \pi @f$
 */
#define KERNEL_CONST_HELMHOLTZBEM3D 0.0795774715459476679

/**
 * @brief Simple struct that containts the wavevector @f$\vec \kappa@f$ and some
 * source point.
 */
struct _helmholtz_data {
  /** @brief The wavevector @f$\vec \kappa@f$. */
  real *kvec;

  /** @brief coupling parameter eta. */
  real eta;

  /** @brief Some source point. */
  real *source;
};

/**
 * @brief Typedef of the struct @ref _helmholtz_data.
 */
typedef struct _helmholtz_data helmholtz_data;

/**
 * @brief Creates a new @ref _bem3d "bem3d"-object for computation of
 * single layer potential matrix of the Helmholtz equation.
 *
 * After calling this function the resulting @ref _bem3d "bem"-object will
 * provide all functionality that is necessary to build up fully populated
 * single layer potential matrix @f$ V \in \mathbb R^{\mathcal I \times
 * \mathcal I}@f$ and also @ref _hmatrix "hmatrix"
 * or @ref _h2matrix "h2matrix" approximation of this matrix.
 *
 * @param k Wavenumber @f$\kappa@f$.
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
 * slp matrices @f$ V @f$ for the helmholtz equation.
 */
HEADER_PREFIX pbem3d
new_slp_helmholtz_bem3d(field k, pcsurface3d gr, uint q_regular,
    uint q_singular, basisfunctionbem3d row_basis, basisfunctionbem3d col_basis);

/**
 * @brief Creates a new @ref _bem3d "bem3d"-object for computation of
 * double layer potential matrix plus a scalar times the mass matrix
 * of the Helmholtz equation.
 *
 * After calling this function the resulting @ref _bem3d "bem"-object will
 * provide all functionality that is necessary to build up fully populated
 * double layer potential matrix @f$ K + \alpha M \in \mathbb
 * R^{\mathcal I \times \mathcal J}@f$ and also @ref _hmatrix "hmatrix"
 * or @ref _h2matrix "h2matrix"
 * approximation of this matrix.
 *
 * @param k Wavenumber @f$\kappa@f$.
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
 * dlp matrices @f$ K + \alpha M @f$ for the Helmholtz equation.
 */
HEADER_PREFIX pbem3d
new_dlp_helmholtz_bem3d(field k, pcsurface3d gr, uint q_regular,
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
 * @param k Wavenumber @f$\kappa@f$.
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
 * adlp matrices @f$ K' + \alpha M @f$ for the Helmholtz equation.
 */
HEADER_PREFIX pbem3d
new_adlp_helmholtz_bem3d(field k, pcsurface3d gr, uint q_regular,
    uint q_singular, basisfunctionbem3d row_basis, basisfunctionbem3d col_basis,
    field alpha);

/**
 * @brief Delete a @ref _bem3d "bem3d" object for the Helmholtz equation
 *
 * @param bem Object to be deleted.
 */
HEADER_PREFIX void
del_helmholtz_bem3d(pbem3d bem);

/**
 * @brief A function based upon the fundamental solution,
 * that will serve as Dirichlet values.
 *
 * When computing the Neumann data out of the Dirichlet data one can use this
 * function as test data which will generate Dirichlet values of with the following
 * values:
 * @f[
 * f(\vec x, \, \vec n) = g(\vec x, \vec z)
 * @f]
 * where @f$\vec z@f$ is some arbitrarily chosen source point and can be set
 * via <tt>data->source</tt>.
 * Corresponding Neumann data can be generated by using
 * @ref rhs_neumann_point_helmholtzbem3d.
 * <br>
 * To build up an appropriate Dirichlet data coefficient vector one needs the
 * @f$ L_2@f$-projection. This can be done by passing this function to
 * @ref projectL2_bem3d_c_avector for piecewise constant basis functions or
 * to @ref projectL2_bem3d_l_avector for piecewise linear basis functions.
 *
 * @param x Evaluation point.
 * @param n Normal vector to current evaluation point.
 * @param data Additional data for evaluating the functional. Here a pointer to
 *   a @ref helmholtz_data object is expected.
 * @return returns the function value of @f$ f(\vec x, \, \vec n) @f$.
 */
HEADER_PREFIX field
rhs_dirichlet_point_helmholtzbem3d(const real *x, const real *n,
    const void *data);

/**
 * @brief A function based upon the fundamental solution,
 * that will serve as Neumann values.
 *
 * When computing the Neumann data out of the Dirichlet data one can use this
 * function as test data which will generate Dirichlet values of with the following
 * values:
 * @f[
 * f(\vec x, \, \vec n) = \frac{\partial g}{\partial \vec n}(\vec x, \vec z)
 * @f]
 * where @f$\vec z@f$ is some arbitrarily chosen source point and can be set
 * via <tt>data->source</tt>.
 * Corresponding Dirichlet data can be generated by using
 * @ref rhs_dirichlet_point_helmholtzbem3d.
 * <br>
 * To build up an appropriate Neumann data coefficient vector one needs the
 * @f$ L_2@f$-projection. This can be done by passing this function to
 * @ref projectL2_bem3d_c_avector for piecewise constant basis functions or
 * to @ref projectL2_bem3d_l_avector for piecewise linear basis functions.
 *
 * @param x Evaluation point.
 * @param n Normal vector to current evaluation point.
 * @param data Additional data for evaluating the functional. Here a pointer to
 *   a @ref helmholtz_data object is expected.
 * @return returns the function value of @f$ f(\vec x, \, \vec n) @f$.
 */
HEADER_PREFIX field
rhs_neumann_point_helmholtzbem3d(const real *x, const real *n, const void *data);

/**
 * @brief A function based upon the fundamental solution,
 * that will serve for Robin boundary conditions.
 *
 * This function is a linear combination of @ref rhs_dirichlet_point_helmholtzbem3d
 * and @ref rhs_neumann_point_helmholtzbem3d and serves as Robin boundary
 * condition.
 * The function is evaluated as:
 * @f[
 * f(\vec x, \, \vec n) = \frac{\partial g}{\partial \vec n}(\vec x, \vec z)
 *   - i \eta g(\vec x, \vec z)
 * @f]
 * where @f$\vec z@f$ is some arbitrarily chosen source point and can be set
 * via <tt>data->source</tt>.
 * The coefficient @f$\eta@f$ can be set via <tt>data->eta</tt>.
 *
 * To build up an appropriate coefficient vector one needs the some integration
 * on the boundary. This can be done by passing this function to
 * @ref integrate_bem3d_c_avector for piecewise constant basis functions or
 * to @ref integrate_bem3d_l_avector for piecewise linear basis functions.
 *
 * @param x Evaluation point.
 * @param n Normal vector to current evaluation point.
 * @param data Additional data for evaluating the functional. Here a pointer to
 *   a @ref helmholtz_data object is expected.
 * @return returns the function value of @f$ f(\vec x, \, \vec n) @f$.
 */
HEADER_PREFIX field
rhs_robin_point_helmholtzbem3d(const real *x, const real *n, const void *data);

/**
 * @brief A function based upon a plane wave,
 * that will serve as Dirichlet values.
 *
 * When computing the Neumann data out of the Dirichlet data one can use this
 * function as test data which will generate Dirichlet values of with the following
 * values:
 * @f[
 * f(\vec x, \, \vec n) = e^{i \langle \vec c, \vec x - \vec s \rangle}
 * @f]
 * where @f$\vec s@f$ is some arbitrarily chosen source point and can be set
 * via <tt>data->source</tt>.
 * Corresponding Neumann data can be generated by using
 * @ref rhs_neumann_plane_helmholtzbem3d.
 * <br>
 * To build up an appropriate Dirichlet data coefficient vector one needs the
 * @f$ L_2@f$-projection. This can be done by passing this function to
 * @ref projectL2_bem3d_c_avector for piecewise constant basis functions or
 * to @ref projectL2_bem3d_l_avector for piecewise linear basis functions.
 *
 * @param x Evaluation point.
 * @param n Normal vector to current evaluation point.
 * @param data Additional data for evaluating the functional. Here a pointer to
 *   a @ref helmholtz_data object is expected.
 * @return returns the function value of @f$ f(\vec x, \, \vec n) @f$.
 */
HEADER_PREFIX field
rhs_dirichlet_plane_helmholtzbem3d(const real *x, const real *n,
    const void *data);

/**
 * @brief A function based upon the plane wave,
 * that will serve as Neumann values.
 *
 * When computing the Neumann data out of the Dirichlet data one can use this
 * function as test data which will generate Dirichlet values of with the following
 * values:
 * @f[
 * f(\vec x, \, \vec n) = \frac{\partial}{\partial \vec n}
 e^{i \langle \vec c, \vec x - \vec s \rangle}
 * @f]
 * where @f$\vec s@f$ is some arbitrarily chosen source point and can be set
 * via <tt>data->source</tt>.
 * Corresponding Dirichlet data can be generated by using
 * @ref rhs_dirichlet_plane_helmholtzbem3d.
 * <br>
 * To build up an appropriate Neumann data coefficient vector one needs the
 * @f$ L_2@f$-projection. This can be done by passing this function to
 * @ref projectL2_bem3d_c_avector for piecewise constant basis functions or
 * to @ref projectL2_bem3d_l_avector for piecewise linear basis functions.
 *
 * @param x Evaluation point.
 * @param n Normal vector to current evaluation point.
 * @param data Additional data for evaluating the functional. Here a pointer to
 *   a @ref helmholtz_data object is expected.
 * @return returns the function value of @f$ f(\vec x, \, \vec n) @f$.
 */
HEADER_PREFIX field
rhs_neumann_plane_helmholtzbem3d(const real *x, const real *n, const void *data);

/**
 * @brief A function based upon the fundamental solution,
 * that will serve for Robin boundary conditions.
 *
 * This function is a linear combination of @ref rhs_dirichlet_plane_helmholtzbem3d
 * and @ref rhs_neumann_plane_helmholtzbem3d and serves as Robin boundary
 * condition.
 * The function is evaluated as:
 * @f[
 * f(\vec x, \, \vec n) = \frac{\partial}{\partial \vec n}
 *    e^{i \langle \vec c, \vec x - \vec s \rangle}
 *   - i \eta e^{i \langle \vec c, \vec x - \vec s \rangle}
 * @f]
 * where @f$\vec s@f$ is some arbitrarily chosen source point and can be set
 * via <tt>data->source</tt>.
 * The coefficient @f$\eta@f$ can be set via <tt>data->eta</tt>.
 *
 * To build up an appropriate coefficient vector one needs the some integration
 * on the boundary. This can be done by passing this function to
 * @ref integrate_bem3d_c_avector for piecewise constant basis functions or
 * to @ref integrate_bem3d_l_avector for piecewise linear basis functions.
 *
 * @param x Evaluation point.
 * @param n Normal vector to current evaluation point.
 * @param data Additional data for evaluating the functional. Here a pointer to
 *   a @ref helmholtz_data object is expected.
 * @return returns the function value of @f$ f(\vec x, \, \vec n) @f$.
 */
HEADER_PREFIX field
rhs_robin_plane_helmholtzbem3d(const real *x, const real *n, const void *data);

/**
 * @}
 */

#endif

#endif /* HELMHOLTZBEM3D_H_ */
