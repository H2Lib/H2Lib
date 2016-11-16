/* ------------------------------------------------------------
 This is the file "laplaceoclbem3d.h" of the H2Lib package.
 All rights reserved, Sven Christophersen 2015
 ------------------------------------------------------------ */

/**
 * @file laplaceoclbem3d.h
 * @author Sven Christophersen
 * @date 2015
 */

#ifndef LIBRARY_LAPLACEOCLBEM3D_H_
#define LIBRARY_LAPLACEOCLBEM3D_H_

/**
 * @addtogroup laplacebem3d
 * @{
 */

#ifdef USE_OPENMP
#ifdef USE_OPENCL

/* C STD LIBRARY */
/* CORE 0 */
/* CORE 1 */
/* CORE 2 */
/* CORE 3 */
/* SIMPLE */
/* PARTICLES */
/* BEM */
#include "oclbem3d.h"
#include "laplacebem3d.h"
#include "laplaceoclbem3d.c"

/**
 * @brief Creates a new @ref _bem3d "bem3d"-object for computation of
 * single layer potential matrix supported by OpenCL.
 *
 * After calling this function the resulting @ref _bem3d "bem"-object will
 * provide all functionality that is necessary to build up fully populated
 * single layer potential matrix @f$ V \in \mathbb R^{\mathcal I \times
 * \mathcal I}@f$ and also @ref _hmatrix "hmatrix"
 * or @ref _h2matrix "h2matrix" approximation of this matrix.
 *
 * The computation intense parts like @ref _bem3d.nearfield are carried out
 * to the GPU, if one uses this constructor and enabled the use of OpenCL.
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
HEADER_PREFIX pbem3d new_slp_laplace_ocl_bem3d(pcsurface3d gr, uint q_regular,
    uint q_singular, basisfunctionbem3d row_basis, basisfunctionbem3d col_basis);

/**
 * @brief Creates a new @ref _bem3d "bem3d"-object for computation of
 * double layer potential matrix supported by OpenCL.
 *
 * After calling this function the resulting @ref _bem3d "bem"-object will
 * provide all functionality that is necessary to build up fully populated
 * double layer potential matrix @f$ K + \alpha M \in \mathbb
 * R^{\mathcal I \times \mathcal J}@f$ and also @ref _hmatrix "hmatrix"
 * or @ref _h2matrix "h2matrix"
 * approximation of this matrix.
 *
 * The computation intense parts like @ref _bem3d.nearfield are carried out
 * to the GPU, if one uses this constructor and enabled the use of OpenCL.
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
 * dlp matrices @f$ K + \alpha M @f$ for the Laplace equation.
 */
HEADER_PREFIX pbem3d new_dlp_laplace_ocl_bem3d(pcsurface3d gr, uint q_regular,
    uint q_singular, basisfunctionbem3d row_basis, basisfunctionbem3d col_basis,
    field alpha);

/**
 * @brief Delete a @ref _bem3d "bem3d" object for the Laplace equation with
 * OpenCL support.
 *
 * @param bem Object to be deleted.
 */
HEADER_PREFIX void del_laplace_ocl_bem3d(pbem3d bem);

#endif
#endif

/**
 * @}
 */

#endif /* LIBRARY_LAPLACEOCLBEM3D_H_ */
