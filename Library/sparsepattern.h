
/* ------------------------------------------------------------
   This is the file "sparsepattern.h" of the H2Lib package.
   All rights reserved, Steffen Boerm 2012
   ------------------------------------------------------------ */

/** @file sparsepattern.h
 *  @author Steffen B&ouml;rm */

#ifndef SPARSEPATTERN_H
#define SPARSEPATTERN_H

/** @defgroup sparsepattern sparsepattern
 *  @brief Representation of the sparsity pattern of a matrix.
 *
 *  @{ */

/** @brief Representation of the sparsity pattern of a matrix. */
typedef struct _sparsepattern sparsepattern;

/** @brief Pointer to @ref sparsepattern object. */
typedef sparsepattern *psparsepattern;

/** @brief Pointer to constant @ref sparsepattern object. */
typedef const sparsepattern *pcsparsepattern;

/** @brief Representation of one non-zero entry in a @ref sparsepattern. */
typedef struct _patentry patentry;

/** @brief Pointer to @ref patentry. */
typedef patentry *ppatentry;

/** @brief Pointer to constant @ref patentry. */
typedef const patentry *pcpatentry;

#include "settings.h"

/** @brief Representation of the sparsity pattern of a matrix.
 *
 *  The pattern is represented by an array of lists:
 *  for the @f$i@f$-th row, @c row[i] is the pointer to the head
 *  of a list of non-zero entries. */
struct _sparsepattern {
  /** @brief Number of rows. */
  uint rows;
  /** @brief Number of columns. */
  uint cols;

  /** @brief Pointers to non-zero lists for all rows. */
  ppatentry *row;
};

/** @brief Representation of one non-zero entry in a @ref sparsepattern.
 *
 *  The non-zero entries are represented by a linked list for each
 *  matrix row. */
struct _patentry {
  /** @brief Row index. */
  uint row;
  /** @brief Column index. */
  uint col;

  /** @brief Next entry in this row. */
  struct _patentry *next;
};

/* ------------------------------------------------------------ *
 * Constructor and destructor                                   *
 * ------------------------------------------------------------ */

/** @brief Create an empty @ref sparsepattern object
 *
 *  @param rows Number of rows of the sparse matrix.
 *  @param cols Number of columns of the sparse matrix.
 *  @returns Empty @ref sparsepattern object. */
HEADER_PREFIX psparsepattern
new_sparsepattern(uint rows, uint cols);

/** @brief Delete a @ref sparsepattern object.
 *
 *  Releases the storage corresponding to this object.
 *
 *  @param sp Object to be deleted. */
HEADER_PREFIX void
del_sparsepattern(psparsepattern sp);

/* ------------------------------------------------------------ *
 * Setting up the sparsity pattern                              *
 * ------------------------------------------------------------ */

/** @brief Remove all non-zero entries from a @ref sparsepattern object
 *
 *  @param sp Object to be cleared. */
HEADER_PREFIX void
clear_sparsepattern(psparsepattern sp);

/** @brief Add a non-zero entry to a @ref sparsepattern object,
 *  unless it already exists.
 *
 *  @param sp Target @ref sparsepattern object.
 *  @param row Row of the non-zero entry.
 *  @param col Column of the non-zero entry. */
HEADER_PREFIX void
addnz_sparsepattern(psparsepattern sp, uint row, uint col);

/* ------------------------------------------------------------ *
 * Simple utility functions                                     *
 * ------------------------------------------------------------ */

/** @brief Print a sparsity pattern.
 *
 *  @param sp Source sparsity pattern. */
HEADER_PREFIX void
print_sparsepattern(pcsparsepattern sp);

/** @} */

#endif
