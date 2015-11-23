/* ------------------------------------------------------------
 This is the file "oclbem3d.h" of the H2Lib package.
 All rights reserved, Sven Christophersen 2015
 ------------------------------------------------------------ */

/**
 * @file oclbem3d.h
 * @author Sven Christophersen
 * @date 2015
 */

#ifndef LIBRARY_OCLBEM3D_H_
#define LIBRARY_OCLBEM3D_H_

#ifdef USE_OPENCL
#ifdef USE_OPENMP

/** @defgroup oclbem3d oclbem3d
 *  @brief BEM-Calculations on OpenCL hardware
 *
 *  This modules defines various functions that are needed to compute integrals
 *  that arise in the BEM computations on OpenCL hardware, mainly for GPGPUs.
 *
 *  @{*/

/* C STD LIBRARY */
#include <string.h>
/* CORE 0 */
#include "opencl.h"
/* CORE 1 */
/* CORE 2 */
/* CORE 3 */
/* SIMPLE */
/* PARTICLES */
/* BEM */
#include "bem3d.h"

/**
 * @brief This structure contains geometry information for the OpenCL computation
 *        aswell as all OpenCL kernels and buffers for transferring work packages
 *        between main memory and GPU memory.
 */
struct _bem3d_ocl {
  /** Total number of OpenCL kernels stored in @ref _bem3d_ocl.num_kernels. */
  cl_uint num_kernels;
  /** Array of all OpenCL kernels */
  cl_kernel *kernels;
  /** Number of gaussian quadrature points in 1D for regular cases. */
  uint nq;
  /** Gaussian quadrature points and weights in 1D for regular cases. */
  cl_mem *mem_q_xw;
  /** Number of gaussian quadrature points in 1D for singular cases. */
  uint nq2;
  /** Gaussian quadrature points and weights in 1D for singular cases. */
  cl_mem *mem_q2_xw;
  /** Coordinates for all vertices. */
  cl_mem *mem_gr_x;
  /** Number of triangles. */
  uint triangles;
  /** Array with number of vertices contained in the each triangle. */
  cl_mem *mem_gr_t;
  /** Buffer that contains resulting matrix entries. */
  cl_mem *mem_N;
  /** Buffer for row indices. */
  cl_mem *mem_ridx;
  /** Buffer for column indices. */
  cl_mem *mem_cidx;
};

/**
 * @brief Abbreviation for struct @ref _bem3d_ocl.
 */
typedef struct _bem3d_ocl bem3d_ocl;

/**
 * @brief Abbreviation for a pointer to a @ref bem3d_ocl object.
 */
typedef bem3d_ocl *pbem3d_ocl;

/**
 * @brief Simple Enum to distinguish between different boundary integral operators.
 */
enum _op_type {
  SLP, //!< SLP Single layer potential
  DLP //!< DLP Double layer potential
};

/**
 * @brief Abbreviation for enum @ref _op_type.
 */
typedef enum _op_type op_type;

/**
 * @brief Struct of callback functions for different quadrature cases that arise
 * in the nearfield computation.
 */
struct _op_callbacks {
  void (*dist)(void*); //!< Callback for distant triangles.
  void (*vert)(void*); //!< Callback for triangles with common vertex.
  void (*edge)(void*); //!< Callback for triangles with common edge.
  void (*iden)(void*); //!< Callback for identical triangles.
  void (*sing_cpu)(void*); //!< Callback CPU computation.
};

/**
 * @brief Abbreviation for struct @ref _op_callbacks.
 */
typedef struct _op_callbacks op_callbacks;

/**
 * @brief Struct to store the parameters of the nearfield integration
 *   routine.
 */
struct _nearfield_args {
  const uint *ridx; //!< Array of row indices.
  const uint *cidx; //!< Array of column indices.
  bool dist; //!< Flag that indices if triangle are disjoint or not.
  pcbem3d bem; //!< @ref bem3d object
  bool ntrans; //!< Flag that indicates the entries are stored in normal or in transposed manner.
  pamatrix N; //!< The matrix where the value should be stored in.
};

/**
 * @brief Abbreviation for struct _nearfield_args.
 */
typedef struct _nearfield_args nearfield_args;

/**
 * @brief Struct to store input data for multiple calls to nearfield integration
 * routine.
 */
struct _merge_data_nf {
  pcbem3d bem; //!< @ref bem3d object.
  field *N; //!< Matrix to store all matrix coefficients.
  uint *ridx; //!< Array of all row indices.
  uint *cidx; //!< Array of all column indices.
  bool *trans; //!< Array of all 'trans' flags
  /**
   * @brief Array of memory addresses, that point to the locations where every
   *  single entry has to be stored.
   */
  field **addr;
  size_t pos; //!< current fill position in all above mentioned arrays.
  op_callbacks *op_cb; //!< Set of callback functions for the current operator.
};

/**
 * @brief Abbreviation for struct @ref _merge_data_nf.
 */
typedef struct _merge_data_nf merge_data_nf;

/**
 * @brief Static variable that contains OpenCL buffers and kernels for
 * bem computations.
 */
extern bem3d_ocl ocl_bem3d;

/**
 * @brief Global variable for general nearfield tasks.
 */
extern ptaskgroup nf;
/**
 * @brief OpenMP lock for @ref nf.
 */
extern omp_lock_t nf_lock;

/**
 * @brief Global variable for distant triangle pairs.
 */
extern ptaskgroup nf_dist;
/**
 * @brief OpenMP lock for @ref nf_dist.
 */
extern omp_lock_t nf_dist_lock;

/**
 * @brief Global variable for common vertex triangle pairs.
 */
extern ptaskgroup nf_vert;
/**
 * @brief OpenMP lock for @ref nf_vert.
 */
extern omp_lock_t nf_vert_lock;

/**
 * @brief Global variable for common edge triangle pairs.
 */
extern ptaskgroup nf_edge;
/**
 * @brief OpenMP lock for @ref nf_edge.
 */
extern omp_lock_t nf_edge_lock;

/**
 * @brief Global variable for identical triangle pairs.
 */
extern ptaskgroup nf_iden;
/**
 * @brief OpenMP lock for @ref nf_iden.
 */
extern omp_lock_t nf_iden_lock;

/**
 * @brief Global variable for integral operator callback functions.
 */
extern op_callbacks op_cb;

/**
 * @brief returns the current size of @ref taskgroup for nearfield computations.
 *
 * @param data Pointer to a struct of @ref nearfield_args.
 * @return The current size of a @ref taskgroup.
 */
HEADER_PREFIX size_t
getsize_nf(void *data);

/**
 * @brief Cleanup intermediate storage for nearfield computation task.
 *
 * @param data Pointer to a struct of @ref nearfield_args.
 */
HEADER_PREFIX void
cleanup_nf_task(void *data);

/**
 * @brief Cleanup intermediate storage from a nearfield computation merge task.
 *
 * @param data Pointer to a struct of @ref nearfield_args.
 */
HEADER_PREFIX void
cleanup_nf_merge(void *data);

/**
 * @brief Splits a neafield computation task into smaller tasks.
 *
 * @param data Pointer to a struct of @ref nearfield_args.
 * @param split Array of smaller subtasks is returned via @p splits.
 * @param n Number of subtasks is returned by @p n.
 */
HEADER_PREFIX void
split_nf_task(void *data, void ***split, uint * n);

/**
 * @brief returns the current size of @ref taskgroup for singular nearfield computations.
 *
 * @param data Pointer to a struct of @ref merge_data_nf.
 * @return The current size of a @ref taskgroup.
 */
HEADER_PREFIX size_t
getsize_nf_sing(void *data);

/**
 * @brief Cleanup intermediate storage for singular nearfield computation task.
 *
 * @param data Pointer to a struct of @ref merge_data_nf.
 */
HEADER_PREFIX void
cleanup_nf_sing_task(void *data);

/**
 * @brief Cleanup intermediate storage from a singular nearfield computation merge task.
 *
 * @param data Pointer to a struct of @ref merge_data_nf.
 */
HEADER_PREFIX void
cleanup_nf_sing_merge(void *data);

/**
 * @brief Splits a singular neafield computation task into smaller tasks.
 *
 * @param data Pointer to a struct of @ref merge_data_nf.
 * @param split Array of smaller subtasks is returned via @p splits.
 * @param n Number of subtasks is returned by @p n.
 */
HEADER_PREFIX void
split_nf_sing_task(void *data, void ***split, uint * n);

/**
 * @brief Close collection of task for disjoint triangle pairs.
 * @param tg @ref taskgroup object to be closed.
 */
HEADER_PREFIX void
close_nf_dist(ptaskgroup tg);

/**
 * @brief Close collection of task for triangle pairs with common vertex.
 * @param tg @ref taskgroup object to be closed.
 */
HEADER_PREFIX void
close_nf_vert(ptaskgroup tg);

/**
 * @brief Close collection of task for triangle pairs with common edge.
 * @param tg @ref taskgroup object to be closed.
 */
HEADER_PREFIX void
close_nf_edge(ptaskgroup tg);

/**
 * @brief Close collection of task for identical triangle pairs.
 * @param tg @ref taskgroup object to be closed.
 */
HEADER_PREFIX void
close_nf_iden(ptaskgroup tg);

/**
 * @brief Merge a list of singular nearfield computation tasks into a single object.
 *
 * @param tg @ref taskgroup to be merged.
 * @param data Object of type @ref merge_data_nf that is being returned.
 */
HEADER_PREFIX void
merge_nf_sing(ptaskgroup tg, void **data);

/**
 * @brief Distribute computational results to their correct positions in memory.
 *
 * @param tg @ref taskgroup to be distributed.
 * @param data Results, that have to be distributed to the tasks in @p tg.
 */
HEADER_PREFIX void
distribute_nf_sing(ptaskgroup tg, void *data);

/**
 * @brief Adds a single singular entry to the corresponding list.
 *
 * @param bem @ref bem3d object.
 * @param op_cb Struct of callback functions for all quadrature cases for the
 *   current operator.
 * @param c Number of common vertices for the triangle @p tt and @p ss.
 * @param tt First triangle.
 * @param ss Second triangle.
 * @param trans Transposed flag
 * @param entry Memory position where this entry should be stored lateron.
 */
HEADER_PREFIX void
add_nf_entry(pcbem3d bem, op_callbacks *op_cb, uint c, uint tt, uint ss,
    bool trans, field * entry);

/**
 * @brief Close collection of task for nearfield calculations.
 * @param tg @ref taskgroup object to be closed.
 */
HEADER_PREFIX void
close_nf(ptaskgroup tg);

/**
 * @brief Merge a list of nearfield computation tasks into a single object.
 *
 * @param tg @ref taskgroup to be merged.
 * @param data Object of type @ref merge_data_nf that is being returned.
 */
HEADER_PREFIX void
merge_nf(ptaskgroup tg, void **data);

/**
 * @brief Distribute computational results to their correct positions in memory.
 *
 * @param tg @ref taskgroup to be distributed.
 * @param data Results, that have to be distributed to the tasks in @p tg.
 */
HEADER_PREFIX void
distribute_nf(ptaskgroup tg, void *data);

/**
 * @}
 */

#endif
#endif

#endif /* LIBRARY_OCLBEM3D_H_ */
