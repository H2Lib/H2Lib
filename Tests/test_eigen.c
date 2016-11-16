
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "eigensolvers.h"
#include "factorizations.h"

static uint problems = 0;

#ifdef USE_FLOAT
static const real tolerance = 5.0e-5;
#else
static const real tolerance = 1.0e-12;
#endif

int
main()
{
  ptridiag  T, Tcopy;
  pamatrix  A, Acopy, Q, U, Vt;
  pavector  work;
  prealavector sigma, lambda;
  real      error;
  uint      rows, cols, mid;
  uint      i, n, iter;
  int       info;

  /* ------------------------------------------------------------
   * Testing symmetric tridiagonal eigenvalue solver
   * ------------------------------------------------------------ */

  n = 6;

  /* Testing symmetric tridiagonal eigenvalue solver */

  (void) printf("==================================================\n"
		"Testing symmetric tridiagonal eigenvalue solver\n"
		"==================================================\n"
		"Setting up %u x %u tridiagonal matrix\n", n, n);
  T = new_tridiag(n);
  for (i = 0; i < n; i++)
    T->d[i] = 2.0;
  for (i = 0; i < n - 1; i++)
    T->l[i] = T->u[i] = -1.0;
  Tcopy = new_tridiag(n);
  copy_tridiag(T, Tcopy);

  A = new_amatrix(n, n);
  clear_amatrix(A);
  for (i = 0; i < n; i++)
    A->a[i + i * A->ld] = T->d[i];
  for (i = 0; i < n - 1; i++) {
    A->a[(i + 1) + i * A->ld] = T->l[i];
    A->a[i + (i + 1) * A->ld] = T->u[i];
  }
  Acopy = clone_amatrix(A);

  Q = new_identity_amatrix(n, n);
  U = new_amatrix(n, n);

  (void) printf("Performing self-made implicit QR iteration\n");
  iter = sb_muleig_tridiag(T, Q, 8 * n);

  (void) printf("  %u iterations\n", iter);

  (void) printf("Checking accuracy\n");
  error = check_ortho_amatrix(false, Q);
  (void) printf("  Orthogonality Q %g, %sokay\n", error,
		(error < tolerance ? "" : "    NOT "));
  if (error >= tolerance)
    problems++;

  copy_amatrix(false, Q, U);
  diageval_tridiag_amatrix(1.0, true, T, true, U);
  addmul_amatrix(-1.0, false, U, true, Q, A);
  error = normfrob_amatrix(A);
  (void) printf("  Accuracy %g. %sokay\n", error,
		(error < tolerance ? "" : "    NOT "));
  if (error >= tolerance)
    problems++;

  (void) printf("Performing default implicit QR iteration\n");
  identity_amatrix(Q);
  i = muleig_tridiag(Tcopy, Q);
  if (i == 0)
    (void) printf("  Success\n");
  else {
    (void) printf("  Failure\n");
    problems++;
  }

  (void) printf("Checking accuracy\n");
  error = check_ortho_amatrix(false, Q);
  (void) printf("  Orthogonality Q %g, %sokay\n", error,
		(error < tolerance ? "" : "    NOT "));
  if (error >= tolerance)
    problems++;

  copy_amatrix(false, Q, U);
  diageval_tridiag_amatrix(1.0, true, Tcopy, true, U);
  addmul_amatrix(-1.0, false, U, true, Q, Acopy);
  error = normfrob_amatrix(Acopy);
  (void) printf("  Accuracy %g. %sokay\n", error,
		(error < tolerance ? "" : "    NOT "));
  if (error >= tolerance)
    problems++;

  del_amatrix(U);
  del_amatrix(Q);
  del_amatrix(Acopy);
  del_amatrix(A);
  del_tridiag(Tcopy);
  del_tridiag(T);

  (void) printf("--------------------------------------------------\n"
		"Setting up random %u x %u tridiagonal matrix\n", n, n);
  T = new_tridiag(n);
  for (i = 0; i < n; i++)
    T->d[i] = 2.0 * rand() / RAND_MAX - 1.0;
  for (i = 0; i < n - 1; i++) {
    T->l[i] = 2.0 * rand() / RAND_MAX - 1.0;
    T->u[i] = CONJ(T->l[i]);
  }
  A = new_amatrix(n, n);
  clear_amatrix(A);
  for (i = 0; i < n; i++)
    A->a[i + i * A->ld] = T->d[i];
  for (i = 0; i < n - 1; i++) {
    A->a[(i + 1) + i * A->ld] = T->l[i];
    A->a[i + (i + 1) * A->ld] = T->u[i];
  }
  Tcopy = new_tridiag(n);
  copy_tridiag(T, Tcopy);
  Acopy = clone_amatrix(A);

  Q = new_identity_amatrix(n, n);
  U = new_amatrix(n, n);

  (void) printf("Performing implicit QR iteration\n");
  iter = sb_muleig_tridiag(T, Q, 8 * n);

  (void) printf("  %u iterations\n", iter);

  (void) printf("Checking accuracy\n");
  error = check_ortho_amatrix(false, Q);
  (void) printf("  Orthogonality Q %g, %sokay\n", error,
		(error < tolerance ? "" : "NOT "));
  if (error >= tolerance)
    problems++;

  copy_amatrix(false, Q, U);
  diageval_tridiag_amatrix(1.0, true, T, true, U);
  addmul_amatrix(-1.0, false, U, true, Q, A);
  error = normfrob_amatrix(A);
  (void) printf("  Accuracy %g, %sokay\n", error,
		(error < tolerance ? "" : "    NOT "));
  if (error >= tolerance)
    problems++;

  (void) printf("Using default eigenvalue solver\n");
  identity_amatrix(Q);
  muleig_tridiag(Tcopy, Q);

  (void) printf("Checking accuracy\n");
  error = check_ortho_amatrix(false, Q);
  (void) printf("  Orthogonality Q %g, %sokay\n", error,
		(error < tolerance ? "" : "NOT "));
  if (error >= tolerance)
    problems++;

  copy_amatrix(false, Q, U);
  diageval_tridiag_amatrix(1.0, true, Tcopy, true, U);
  addmul_amatrix(-1.0, false, U, true, Q, Acopy);
  error = normfrob_amatrix(Acopy);
  (void) printf("  Accuracy %g, %sokay\n", error,
		(error < tolerance ? "" : "    NOT "));
  if (error >= tolerance)
    problems++;

  del_amatrix(U);
  del_amatrix(Q);
  del_amatrix(Acopy);
  del_amatrix(A);
  del_tridiag(Tcopy);
  del_tridiag(T);

  /* ------------------------------------------------------------
   * Testing self-adjoint matrix eigenvalue solver
   * ------------------------------------------------------------ */

  (void) printf("==================================================\n"
		"Testing self-adjoint matrix eigenvalue solver\n"
		"==================================================\n"
		"Setting up random %u x %u self-adjoint matrix\n", n, n);
  A = new_amatrix(n, n);
  random_selfadjoint_amatrix(A);

  Acopy = new_amatrix(n, n);
  copy_amatrix(false, A, Acopy);

  lambda = new_realavector(n);
  Q = new_identity_amatrix(n, n);

  (void) printf("Performing implicit QR iteration\n");
  iter = sb_eig_amatrix(A, lambda, Q, 8 * n);

  (void) printf("  %u iterations\n", iter);

  (void) printf("Checking accuracy\n");
  error = check_ortho_amatrix(false, Q);
  (void) printf("  Orthogonality Q %g, %sokay\n", error,
		(error < tolerance ? "" : "NOT "));
  if (error >= tolerance)
    problems++;

  U = new_amatrix(n, n);
  copy_amatrix(false, Q, U);
  copy_amatrix(false, Acopy, A);
  diageval_realavector_amatrix(1.0, true, lambda, true, U);
  addmul_amatrix(-1.0, false, U, true, Q, A);
  error = normfrob_amatrix(A);
  (void) printf("  Accuracy %g, %sokay\n", error,
		(error < tolerance ? "" : "NOT "));
  if (error >= tolerance)
    problems++;

  (void) printf("Using default eigenvalue solver\n");
  copy_amatrix(false, Acopy, A);
  info = eig_amatrix(A, lambda, Q);
  assert(info == 0);

  (void) printf("Checking accuracy\n");
  error = check_ortho_amatrix(false, Q);
  (void) printf("  Orthogonality Q %g, %sokay\n", error,
		(error < tolerance ? "" : "NOT "));
  if (error >= tolerance)
    problems++;

  copy_amatrix(false, Q, U);
  diageval_realavector_amatrix(1.0, true, lambda, true, U);
  addmul_amatrix(-1.0, false, U, true, Q, Acopy);
  error = normfrob_amatrix(Acopy);
  (void) printf("  Accuracy %g, %sokay\n", error,
		(error < tolerance ? "" : "NOT "));
  if (error >= tolerance)
    problems++;

  del_amatrix(U);
  del_amatrix(Q);
  del_realavector(lambda);
  del_amatrix(Acopy);
  del_amatrix(A);

  /* ------------------------------------------------------------
   * Testing bidiagonal SVD solver
   * ------------------------------------------------------------ */

  (void) printf("==================================================\n"
		"Testing bidiagonal SVD solver\n"
		"==================================================\n"
		"Setting up bidiagonal %u x %u matrix\n", n, n);
  T = new_tridiag(n);
  for (i = 0; i < n; i++)
    T->d[i] = i + 1.0;
  for (i = 0; i < n - 1; i++) {
    T->l[i] = 1.0;
    T->u[i] = 0.0;
  }
  A = new_amatrix(n, n);
  clear_amatrix(A);
  for (i = 0; i < n; i++)
    A->a[i + i * A->ld] = T->d[i];
  for (i = 0; i < n - 1; i++)
    A->a[(i + 1) + i * A->ld] = T->l[i];

  Tcopy = new_tridiag(n);
  copy_tridiag(T, Tcopy);
  Acopy = new_amatrix(n, n);
  copy_amatrix(false, A, Acopy);

  U = new_identity_amatrix(n, n);
  Vt = new_identity_amatrix(n, n);

  (void) printf("Performing self-made implicit SVD iteration\n");
  iter = sb_mulsvd_tridiag(T, U, Vt, 8 * n);
  (void) printf("  %u iterations\n", iter);

  (void) printf("Checking accuracy\n");
  error = check_ortho_amatrix(false, U);
  (void) printf("  Orthogonality U %g, %sokay\n", error,
		(error < tolerance ? "" : "NOT "));
  if (error >= tolerance)
    problems++;

  error = check_ortho_amatrix(true, Vt);
  (void) printf("  Orthogonality Vt %g, %sokay\n", error,
		(error < tolerance ? "" : "NOT "));
  if (error >= tolerance)
    problems++;

  diageval_tridiag_amatrix(1.0, true, T, true, U);
  addmul_amatrix(-1.0, false, U, false, Vt, A);
  error = normfrob_amatrix(A);
  (void) printf("  Accuracy %g, %sokay\n", error,
		(error < tolerance ? "" : "NOT "));
  if (error >= tolerance)
    problems++;

  (void) printf("Using default SVD solver\n");
  copy_tridiag(Tcopy, T);
  identity_amatrix(U);
  identity_amatrix(Vt);
  mulsvd_tridiag(T, U, Vt);

  (void) printf("Checking accuracy\n");
  error = check_ortho_amatrix(false, U);
  (void) printf("  Orthogonality U %g, %sokay\n", error,
		(error < tolerance ? "" : "NOT "));
  if (error >= tolerance)
    problems++;

  error = check_ortho_amatrix(true, Vt);
  (void) printf("  Orthogonality Vt %g, %sokay\n", error,
		(error < tolerance ? "" : "NOT "));
  if (error >= tolerance)
    problems++;

  copy_amatrix(false, Acopy, A);
  diageval_tridiag_amatrix(1.0, true, T, true, U);
  addmul_amatrix(-1.0, false, U, false, Vt, A);
  error = normfrob_amatrix(A);
  (void) printf("  Accuracy %g, %sokay\n", error,
		(error < tolerance ? "" : "NOT "));
  if (error >= tolerance)
    problems++;

  del_amatrix(Vt);
  del_amatrix(U);
  del_amatrix(Acopy);
  del_tridiag(Tcopy);
  del_amatrix(A);
  del_tridiag(T);

  (void) printf("--------------------------------------------------\n"
		"Setting up random bidiagonal %u x %u matrix\n", n, n);
  T = new_tridiag(n);
  for (i = 0; i < n; i++) {
    T->d[i] = 2.0 * rand() / RAND_MAX - 1.0;
  }
  for (i = 0; i < n - 1; i++) {
    T->l[i] = 2.0 * rand() / RAND_MAX - 1.0;
    T->u[i] = 0.0;
  }
  A = new_amatrix(n, n);
  clear_amatrix(A);
  for (i = 0; i < n; i++)
    A->a[i + i * A->ld] = T->d[i];
  for (i = 0; i < n - 1; i++)
    A->a[(i + 1) + i * A->ld] = T->l[i];

  Tcopy = new_tridiag(n);
  copy_tridiag(T, Tcopy);
  Acopy = clone_amatrix(A);

  U = new_identity_amatrix(n, n);
  Vt = new_identity_amatrix(n, n);

  (void) printf("Performing implicit SVD iteration\n");
  iter = sb_mulsvd_tridiag(T, U, Vt, 8 * n);
  (void) printf("  %u iterations\n", iter);

  (void) printf("Checking accuracy\n");
  error = check_ortho_amatrix(false, U);
  (void) printf("  Orthogonality U %g, %sokay\n", error,
		(error < tolerance ? "" : "NOT "));
  if (error >= tolerance)
    problems++;

  error = check_ortho_amatrix(true, Vt);
  (void) printf("  Orthogonality Vt %g, %sokay\n", error,
		(error < tolerance ? "" : "NOT "));
  if (error >= tolerance)
    problems++;

  diageval_tridiag_amatrix(1.0, true, T, true, U);
  addmul_amatrix(-1.0, false, U, false, Vt, A);
  error = normfrob_amatrix(A);
  (void) printf("  Accuracy %g, %sokay\n", error,
		(error < tolerance ? "" : "NOT "));
  if (error >= tolerance)
    problems++;

  (void) printf("Using default SVD solver\n");
  copy_tridiag(Tcopy, T);
  copy_amatrix(false, Acopy, A);
  identity_amatrix(U);
  identity_amatrix(Vt);
  mulsvd_tridiag(T, U, Vt);

  (void) printf("Checking accuracy\n");
  error = check_ortho_amatrix(false, U);
  (void) printf("  Orthogonality U %g, %sokay\n", error,
		(error < tolerance ? "" : "NOT "));
  if (error >= tolerance)
    problems++;

  error = check_ortho_amatrix(true, Vt);
  (void) printf("  Orthogonality Vt %g, %sokay\n", error,
		(error < tolerance ? "" : "NOT "));
  if (error >= tolerance)
    problems++;

  diageval_tridiag_amatrix(1.0, true, T, true, U);
  addmul_amatrix(-1.0, false, U, false, Vt, A);
  error = normfrob_amatrix(A);
  (void) printf("  Accuracy %g, %sokay\n", error,
		(error < tolerance ? "" : "NOT "));
  if (error >= tolerance)
    problems++;

  del_amatrix(Vt);
  del_amatrix(U);
  del_amatrix(Acopy);
  del_tridiag(Tcopy);
  del_amatrix(A);
  del_tridiag(T);

  /* ------------------------------------------------------------
   * Testing Golub-Kahan bidiagonalization
   * ------------------------------------------------------------ */

  rows = 10;
  cols = 7;
  mid = UINT_MIN(rows, cols);
  (void) printf("==================================================\n"
		"Testing Golub-Kahan bidiagonalization\n"
		"==================================================\n"
		"Setting up random %u x %u matrix\n", rows, cols);
  A = new_amatrix(rows, cols);
  random_amatrix(A);
  Acopy = new_amatrix(rows, cols);
  copy_amatrix(false, A, Acopy);
  U = new_amatrix(rows, mid);
  Vt = new_amatrix(mid, cols);
  T = new_tridiag(mid);

  (void) printf("Bidiagonalizing\n");
  bidiagonalize_amatrix(A, T, U, Vt);

  (void) printf("Checking accuracy\n");
  error = check_ortho_amatrix(false, U);
  (void) printf("  Orthogonality U %g, %sokay\n", error,
		(error < tolerance ? "" : "NOT "));
  if (error >= tolerance)
    problems++;

  error = check_ortho_amatrix(true, Vt);
  (void) printf("  Orthogonality Vt %g, %sokay\n", error,
		(error < tolerance ? "" : "NOT "));
  if (error >= tolerance)
    problems++;

  lowereval_tridiag_amatrix(1.0, true, T, true, U);
  addmul_amatrix(-1.0, false, U, false, Vt, Acopy);
  error = normfrob_amatrix(Acopy);
  (void) printf("  Accuracy %g, %sokay\n", error,
		(error < tolerance ? "" : "NOT "));
  if (error >= tolerance)
    problems++;

  del_tridiag(T);
  del_amatrix(Vt);
  del_amatrix(U);
  del_amatrix(Acopy);
  del_amatrix(A);

  rows = 8;
  cols = 15;
  mid = UINT_MIN(rows, cols);
  (void) printf("--------------------------------------------------\n"
		"Setting up %u x %u matrix\n", rows, cols);
  A = new_amatrix(rows, cols);
  random_amatrix(A);
  Acopy = new_amatrix(rows, cols);
  copy_amatrix(false, A, Acopy);
  U = new_amatrix(rows, mid);
  Vt = new_amatrix(mid, cols);
  T = new_tridiag(mid);

  (void) printf("Bidiagonalizing\n");
  bidiagonalize_amatrix(A, T, U, Vt);

  (void) printf("Checking accuracy\n");
  error = check_ortho_amatrix(false, U);
  (void) printf("  Orthogonality U %g, %sokay\n", error,
		(error < tolerance ? "" : "NOT "));
  if (error >= tolerance)
    problems++;

  error = check_ortho_amatrix(true, Vt);
  (void) printf("  Orthogonality Vt %g, %sokay\n", error,
		(error < tolerance ? "" : "NOT "));
  if (error >= tolerance)
    problems++;

  lowereval_tridiag_amatrix(1.0, true, T, true, U);
  addmul_amatrix(-1.0, false, U, false, Vt, Acopy);
  error = normfrob_amatrix(Acopy);
  (void) printf("  Accuracy %g, %sokay\n", error,
		(error < tolerance ? "" : "NOT "));
  if (error >= tolerance)
    problems++;

  del_tridiag(T);
  del_amatrix(Vt);
  del_amatrix(U);
  del_amatrix(Acopy);
  del_amatrix(A);

  /* ------------------------------------------------------------
   * Testing general SVD solver
   * ------------------------------------------------------------ */

  (void) printf("==================================================\n"
		"Testing general SVD solver\n"
		"==================================================\n"
		"Setting up 3 x 4 matrix\n");
  A = new_amatrix(3, 4);
  setentry_amatrix(A, 0, 0, 1.0);
  setentry_amatrix(A, 1, 0, 2.0);
  setentry_amatrix(A, 2, 0, 3.0);
  setentry_amatrix(A, 0, 1, 2.0);
  setentry_amatrix(A, 1, 1, 4.0);
  setentry_amatrix(A, 2, 1, 6.0);
  setentry_amatrix(A, 0, 2, 2.0);
  setentry_amatrix(A, 1, 2, 5.0);
  setentry_amatrix(A, 2, 2, 8.0);
  setentry_amatrix(A, 0, 3, 1.0);
  setentry_amatrix(A, 1, 3, 4.0);
  setentry_amatrix(A, 2, 3, 7.0);

  Acopy = new_amatrix(A->rows, A->cols);
  copy_amatrix(false, A, Acopy);

  U = new_identity_amatrix(3, 3);
  Vt = new_identity_amatrix(3, 4);
  sigma = new_realavector(3);
  work = new_avector(3 * 3);

  (void) printf("Running self-made SVD solver\n");
  iter = sb_svd_amatrix(A, sigma, U, Vt, 24);

  (void) printf("  %u iterations\n", iter);

  (void) printf("Checking accuracy\n");
  error = check_ortho_amatrix(false, U);
  (void) printf("  Orthogonality U %g, %sokay\n", error,
		(error < tolerance ? "" : "NOT "));
  if (error >= tolerance)
    problems++;

  error = check_ortho_amatrix(true, Vt);
  (void) printf("  Orthogonality Vt %g, %sokay\n", error,
		(error < tolerance ? "" : "NOT "));
  if (error >= tolerance)
    problems++;

  diageval_realavector_amatrix(1.0, true, sigma, true, U);
  addmul_amatrix(-1.0, false, U, false, Vt, Acopy);
  error = normfrob_amatrix(Acopy);
  (void) printf("  Accuracy %g, %sokay\n", error,
		(error < tolerance ? "" : "NOT "));
  if (error >= tolerance)
    problems++;

  del_avector(work);
  del_realavector(sigma);
  del_amatrix(Vt);
  del_amatrix(U);
  del_amatrix(Acopy);
  del_amatrix(A);

  (void) printf("--------------------------------------------------\n"
		"Setting up 4 x 3 matrix\n");
  A = new_amatrix(4, 3);
  setentry_amatrix(A, 0, 0, 1.0);
  setentry_amatrix(A, 0, 1, 2.0);
  setentry_amatrix(A, 0, 2, 3.0);
  setentry_amatrix(A, 1, 0, 2.0);
  setentry_amatrix(A, 1, 1, 4.0);
  setentry_amatrix(A, 1, 2, 6.0);
  setentry_amatrix(A, 2, 0, 2.0);
  setentry_amatrix(A, 2, 1, 5.0);
  setentry_amatrix(A, 2, 2, 8.0);
  setentry_amatrix(A, 3, 0, 1.0);
  setentry_amatrix(A, 3, 1, 4.0);
  setentry_amatrix(A, 3, 2, 7.0);

  Acopy = new_amatrix(A->rows, A->cols);
  copy_amatrix(false, A, Acopy);

  U = new_amatrix(4, 3);
  Vt = new_amatrix(3, 3);
  sigma = new_realavector(3);
  work = new_avector(3 * 3);

  (void) printf("Running self-made SVD solver\n");
  iter = sb_svd_amatrix(A, sigma, U, Vt, 24);
  (void) printf("  %u iterations\n", iter);

  (void) printf("Checking accuracy\n");
  error = check_ortho_amatrix(false, U);
  (void) printf("  Orthogonality U %g, %sokay\n", error,
		(error < tolerance ? "" : "NOT "));
  if (error >= tolerance)
    problems++;

  error = check_ortho_amatrix(true, Vt);
  (void) printf("  Orthogonality V %g, %sokay\n", error,
		(error < tolerance ? "" : "NOT "));
  if (error >= tolerance)
    problems++;

  diageval_realavector_amatrix(1.0, true, sigma, true, U);
  addmul_amatrix(-1.0, false, U, false, Vt, Acopy);
  error = normfrob_amatrix(Acopy);
  (void) printf("  Accuracy %g, %sokay\n", error,
		(error < tolerance ? "" : "NOT "));
  if (error >= tolerance)
    problems++;

  del_avector(work);
  del_realavector(sigma);
  del_amatrix(Vt);
  del_amatrix(U);
  del_amatrix(Acopy);
  del_amatrix(A);

  (void) printf("--------------------------------------------------\n"
		"Setting up 4 x 3 matrix\n");
  A = new_amatrix(4, 3);
  setentry_amatrix(A, 0, 0, 1.0);
  setentry_amatrix(A, 0, 1, 2.0);
  setentry_amatrix(A, 0, 2, 3.0);
  setentry_amatrix(A, 1, 0, 2.0);
  setentry_amatrix(A, 1, 1, 4.0);
  setentry_amatrix(A, 1, 2, 6.0);
  setentry_amatrix(A, 2, 0, 2.0);
  setentry_amatrix(A, 2, 1, 5.0);
  setentry_amatrix(A, 2, 2, 8.0);
  setentry_amatrix(A, 3, 0, 1.0);
  setentry_amatrix(A, 3, 1, 4.0);
  setentry_amatrix(A, 3, 2, 7.0);

  Acopy = clone_amatrix(A);

  U = new_amatrix(4, 3);
  Vt = new_amatrix(3, 3);
  sigma = new_realavector(3);
  work = new_avector(3 * 3);

  (void) printf("Running self-made SVD solver\n");
  iter = sb_svd_amatrix(A, sigma, U, Vt, 24);
  (void) printf("  %u iterations\n", iter);

  (void) printf("Checking accuracy\n");
  error = check_ortho_amatrix(false, U);
  (void) printf("  Orthogonality U %g, %sokay\n", error,
		(error < tolerance ? "" : "NOT "));
  if (error >= tolerance)
    problems++;

  error = check_ortho_amatrix(true, Vt);
  (void) printf("  Orthogonality Vt %g, %sokay\n", error,
		(error < tolerance ? "" : "NOT "));
  if (error >= tolerance)
    problems++;

  diageval_realavector_amatrix(1.0, true, sigma, true, U);
  addmul_amatrix(-1.0, false, U, false, Vt, Acopy);
  error = normfrob_amatrix(Acopy);
  (void) printf("  Accuracy %g, %sokay\n", error,
		(error < tolerance ? "" : "NOT "));
  if (error >= tolerance)
    problems++;

  del_avector(work);
  del_realavector(sigma);
  del_amatrix(Vt);
  del_amatrix(U);
  del_amatrix(Acopy);
  del_amatrix(A);

  rows = 9;
  cols = 7;
  mid = UINT_MIN(rows, cols);
  (void) printf("--------------------------------------------------\n"
		"Setting up random %u x %u matrix\n", rows, cols);
  A = new_amatrix(rows, cols);
  random_amatrix(A);

  Acopy = new_amatrix(A->rows, A->cols);
  copy_amatrix(false, A, Acopy);

  U = new_amatrix(rows, mid);
  Vt = new_amatrix(mid, cols);
  sigma = new_realavector(mid);
  work = new_avector(3 * mid);

  (void) printf("Running self-made SVD solver\n");
  iter = sb_svd_amatrix(A, sigma, U, Vt, 10 * mid);
  (void) printf("  %u iterations\n", iter);

  (void) printf("Checking accuracy\n");
  error = check_ortho_amatrix(false, U);
  (void) printf("  Orthogonality U %g, %sokay\n", error,
		(error < tolerance ? "" : "NOT "));
  if (error >= tolerance)
    problems++;

  error = check_ortho_amatrix(true, Vt);
  (void) printf("  Orthogonality Vt %g, %sokay\n", error,
		(error < tolerance ? "" : "NOT "));
  if (error >= tolerance)
    problems++;

  diageval_realavector_amatrix(1.0, true, sigma, true, U);
  addmul_amatrix(-1.0, false, U, false, Vt, Acopy);
  error = normfrob_amatrix(Acopy);
  (void) printf("  Accuracy %g, %sokay\n", error,
		(error < tolerance ? "" : "NOT "));
  if (error >= tolerance)
    problems++;

  del_avector(work);
  del_realavector(sigma);
  del_amatrix(Vt);
  del_amatrix(U);
  del_amatrix(Acopy);
  del_amatrix(A);

  rows = 10;
  cols = 6;
  mid = UINT_MIN(rows, cols);
  (void) printf("--------------------------------------------------\n"
		"Setting up random %u x %u matrix\n", rows, cols);
  A = new_amatrix(rows, cols);
  random_amatrix(A);

  Acopy = new_amatrix(A->rows, A->cols);
  copy_amatrix(false, A, Acopy);

  U = new_amatrix(rows, mid);
  Vt = new_amatrix(mid, cols);
  sigma = new_realavector(mid);
  work = new_avector(3 * mid);

  (void) printf("Running self-made SVD solver\n");
  iter = sb_svd_amatrix(A, sigma, U, Vt, 10 * mid);
  (void) printf("  %u iterations\n", iter);

  (void) printf("Checking accuracy\n");
  error = check_ortho_amatrix(false, U);
  (void) printf("  Orthogonality U %g, %sokay\n", error,
		(error < tolerance ? "" : "NOT "));
  if (error >= tolerance)
    problems++;

  error = check_ortho_amatrix(true, Vt);
  (void) printf("  Orthogonality Vt %g, %sokay\n", error,
		(error < tolerance ? "" : "NOT "));
  if (error >= tolerance)
    problems++;

  copy_amatrix(false, Acopy, A);
  diageval_realavector_amatrix(1.0, true, sigma, true, U);
  addmul_amatrix(-1.0, false, U, false, Vt, A);
  error = normfrob_amatrix(A);
  (void) printf("  Accuracy %g, %sokay\n", error,
		(error < tolerance ? "" : "NOT "));
  if (error >= tolerance)
    problems++;

  (void) printf("Running default SVD solver\n");
  copy_amatrix(false, Acopy, A);
  svd_amatrix(A, sigma, U, Vt);

  (void) printf("Checking accuracy\n");
  error = check_ortho_amatrix(false, U);
  (void) printf("  Orthogonality U %g, %sokay\n", error,
		(error < tolerance ? "" : "NOT "));
  if (error >= tolerance)
    problems++;

  error = check_ortho_amatrix(true, Vt);
  (void) printf("  Orthogonality Vt %g, %sokay\n", error,
		(error < tolerance ? "" : "NOT "));
  if (error >= tolerance)
    problems++;

  copy_amatrix(false, Acopy, A);
  diageval_realavector_amatrix(1.0, true, sigma, true, U);
  addmul_amatrix(-1.0, false, U, false, Vt, A);
  error = normfrob_amatrix(A);
  (void) printf("  Accuracy %g, %sokay\n", error,
		(error < tolerance ? "" : "NOT "));
  if (error >= tolerance)
    problems++;

  del_avector(work);
  del_realavector(sigma);
  del_amatrix(Vt);
  del_amatrix(U);
  del_amatrix(Acopy);
  del_amatrix(A);

  rows = 7;
  cols = 13;
  mid = UINT_MIN(rows, cols);
  (void) printf("--------------------------------------------------\n"
		"Setting up random %u x %u matrix\n", rows, cols);
  A = new_amatrix(rows, cols);
  random_amatrix(A);

  Acopy = new_amatrix(A->rows, A->cols);
  copy_amatrix(false, A, Acopy);

  U = new_amatrix(rows, mid);
  Vt = new_amatrix(mid, cols);
  sigma = new_realavector(mid);
  work = new_avector(3 * mid);

  (void) printf("Running self-made SVD solver\n");
  iter = sb_svd_amatrix(A, sigma, U, Vt, 10 * mid);
  (void) printf("  %u iterations\n", iter);

  (void) printf("Checking accuracy\n");
  error = check_ortho_amatrix(false, U);
  (void) printf("  Orthogonality U %g, %sokay\n", error,
		(error < tolerance ? "" : "NOT "));
  if (error >= tolerance)
    problems++;

  error = check_ortho_amatrix(true, Vt);
  (void) printf("  Orthogonality Vt %g, %sokay\n", error,
		(error < tolerance ? "" : "NOT "));
  if (error >= tolerance)
    problems++;

  copy_amatrix(false, Acopy, A);
  diageval_realavector_amatrix(1.0, true, sigma, true, U);
  addmul_amatrix(-1.0, false, U, false, Vt, A);
  error = normfrob_amatrix(A);
  (void) printf("  Accuracy %g, %sokay\n", error,
		(error < tolerance ? "" : "NOT "));
  if (error >= tolerance)
    problems++;

  (void) printf("Running default SVD solver\n");
  copy_amatrix(false, Acopy, A);
  svd_amatrix(A, sigma, U, Vt);

  (void) printf("Checking accuracy\n");
  error = check_ortho_amatrix(false, U);
  (void) printf("  Orthogonality U %g, %sokay\n", error,
		(error < tolerance ? "" : "NOT "));
  if (error >= tolerance)
    problems++;

  error = check_ortho_amatrix(true, Vt);
  (void) printf("  Orthogonality Vt %g, %sokay\n", error,
		(error < tolerance ? "" : "NOT "));
  if (error >= tolerance)
    problems++;

  copy_amatrix(false, Acopy, A);
  diageval_realavector_amatrix(1.0, true, sigma, true, U);
  addmul_amatrix(-1.0, false, U, false, Vt, A);
  error = normfrob_amatrix(A);
  (void) printf("  Accuracy %g, %sokay\n", error,
		(error < tolerance ? "" : "NOT "));
  if (error >= tolerance)
    problems++;

  del_avector(work);
  del_realavector(sigma);
  del_amatrix(Vt);
  del_amatrix(U);
  del_amatrix(Acopy);
  del_amatrix(A);

  printf("----------------------------------------\n"
	 "  %u matrices and\n"
	 "  %u vectors still active\n"
	 "  %u errors found\n", getactives_amatrix(), getactives_avector(),
	 problems);

  return problems;
}
