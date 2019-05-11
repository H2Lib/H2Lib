#include "matrixnorms.h"

/****************************************************
 * Norm2diff for amatrix
 ****************************************************/

real
norm2diff_amatrix_rkmatrix(pcrkmatrix a, pcamatrix b)
{
  return norm2diff_matrix((mvm_t) mvm_rkmatrix_avector, (void *) a,
			  (mvm_t) mvm_amatrix_avector, (void *) b, a->A.rows,
			  a->B.rows);
}

real
norm2diff_amatrix_sparsematrix(pcsparsematrix a, pcamatrix b)
{
  return norm2diff_matrix((mvm_t) mvm_sparsematrix_avector, (void *) a,
			  (mvm_t) mvm_amatrix_avector, (void *) b, a->rows,
			  a->cols);
}

real
norm2diff_amatrix_hmatrix(pchmatrix a, pcamatrix b)
{
  return norm2diff_matrix((mvm_t) mvm_hmatrix_avector, (void *) a,
			  (mvm_t) mvm_amatrix_avector, (void *) b,
			  a->rc->size, a->cc->size);
}

real
norm2diff_amatrix_h2matrix(pch2matrix a, pcamatrix b)
{
  return norm2diff_matrix((mvm_t) mvm_h2matrix_avector, (void *) a,
			  (mvm_t) mvm_amatrix_avector, (void *) b,
			  a->rb->t->size, a->cb->t->size);
}

real
norm2diff_amatrix_dh2matrix(pcdh2matrix a, pcamatrix b)
{
  return norm2diff_matrix((mvm_t) mvm_dh2matrix_avector, (void *) a,
			  (mvm_t) mvm_amatrix_avector, (void *) b,
			  a->rb->t->size, a->cb->t->size);
}

/****************************************************
 * Norm2diff for rkmatrix
 ****************************************************/

real
norm2diff_rkmatrix_sparsematrix(pcsparsematrix a, pcrkmatrix b)
{
  return norm2diff_matrix((mvm_t) mvm_sparsematrix_avector, (void *) a,
			  (mvm_t) mvm_rkmatrix_avector, (void *) b, a->rows,
			  a->rows);
}

real
norm2diff_rkmatrix_hmatrix(pchmatrix a, pcrkmatrix b)
{
  return norm2diff_matrix((mvm_t) mvm_hmatrix_avector, (void *) a,
			  (mvm_t) mvm_rkmatrix_avector, (void *) b,
			  a->rc->size, a->cc->size);
}

real
norm2diff_rkmatrix_h2matrix(pch2matrix a, pcrkmatrix b)
{
  return norm2diff_matrix((mvm_t) mvm_h2matrix_avector, (void *) a,
			  (mvm_t) mvm_rkmatrix_avector, (void *) b,
			  a->rb->t->size, a->cb->t->size);
}

real
norm2diff_rkmatrix_dh2matrix(pcdh2matrix a, pcrkmatrix b)
{
  return norm2diff_matrix((mvm_t) mvm_dh2matrix_avector, (void *) a,
			  (mvm_t) mvm_rkmatrix_avector, (void *) b,
			  a->rb->t->size, a->cb->t->size);
}

/****************************************************
 * Norm2diff for sparsematrix
 ****************************************************/

real
norm2diff_sparsematrix_hmatrix(pchmatrix a, pcsparsematrix b)
{
  return norm2diff_matrix((mvm_t) mvm_hmatrix_avector, (void *) a,
			  (mvm_t) mvm_sparsematrix_avector, (void *) b,
			  a->rc->size, a->cc->size);
}

real
norm2diff_sparsematrix_h2matrix(pch2matrix a, pcsparsematrix b)
{
  return norm2diff_matrix((mvm_t) mvm_h2matrix_avector, (void *) a,
			  (mvm_t) mvm_sparsematrix_avector, (void *) b,
			  a->rb->t->size, a->cb->t->size);
}

real
norm2diff_sparsematrix_dh2matrix(pcdh2matrix a, pcsparsematrix b)
{
  return norm2diff_matrix((mvm_t) mvm_dh2matrix_avector, (void *) a,
			  (mvm_t) mvm_sparsematrix_avector, (void *) b,
			  a->rb->t->size, a->cb->t->size);
}

/****************************************************
 * Norm2diff for hmatrix
 ****************************************************/

real
norm2diff_hmatrix_h2matrix(pch2matrix a, pchmatrix b)
{
  return norm2diff_matrix((mvm_t) mvm_h2matrix_avector, (void *) a,
			  (mvm_t) mvm_hmatrix_avector, (void *) b,
			  a->rb->t->size, a->cb->t->size);
}

real
norm2diff_hmatrix_dh2matrix(pcdh2matrix a, pchmatrix b)
{
  return norm2diff_matrix((mvm_t) mvm_dh2matrix_avector, (void *) a,
			  (mvm_t) mvm_hmatrix_avector, (void *) b,
			  a->rb->t->size, a->cb->t->size);
}

/****************************************************
 * Norm2diff for h2matrix
 ****************************************************/

real
norm2diff_h2matrix_dh2matrix(pcdh2matrix a, pch2matrix b)
{
  return norm2diff_matrix((mvm_t) mvm_dh2matrix_avector, (void *) a,
			  (mvm_t) mvm_h2matrix_avector, (void *) b,
			  a->rb->t->size, a->cb->t->size);
}

/****************************************************
 * Norm2diff_pre_matrix for amatrix
 ****************************************************/

real
norm2diff_lr_amatrix(pcamatrix A, pcamatrix LR)
{
  return norm2diff_pre_matrix((mvm_t) mvm_amatrix_avector, (void *) A,
			      (prcd_t) lreval_n_amatrix_avector,
			      (prcd_t) lreval_t_amatrix_avector, (void *) LR,
			      A->rows, A->cols);
}

real
norm2diff_lr_amatrix_hmatrix(pcamatrix A, pchmatrix LR)
{
  return norm2diff_pre_matrix((mvm_t) mvm_amatrix_avector, (void *) A,
			      (prcd_t) lreval_n_hmatrix_avector,
			      (prcd_t) lreval_t_hmatrix_avector, (void *) LR,
			      A->rows, A->cols);
}

real
norm2diff_chol_amatrix(pcamatrix A, pcamatrix chol)
{
  return norm2diff_pre_matrix((mvm_t) mvm_amatrix_avector, (void *) A,
			      (prcd_t) choleval_amatrix_avector,
			      (prcd_t) choleval_amatrix_avector,
			      (void *) chol, A->rows, A->cols);
}

real
norm2diff_chol_amatrix_hmatrix(pcamatrix A, pchmatrix chol)
{
  return norm2diff_pre_matrix((mvm_t) mvm_amatrix_avector, (void *) A,
			      (prcd_t) choleval_hmatrix_avector,
			      (prcd_t) choleval_hmatrix_avector,
			      (void *) chol, A->rows, A->cols);
}

/****************************************************
 * Norm2diff_pre_matrix for hmatrix
 ****************************************************/

real
norm2diff_lr_hmatrix_amatrix(pchmatrix A, pcamatrix LR)
{
  return norm2diff_pre_matrix((mvm_t) mvm_hmatrix_avector, (void *) A,
			      (prcd_t) lreval_n_amatrix_avector,
			      (prcd_t) lreval_t_amatrix_avector, (void *) LR,
			      A->rc->size, A->cc->size);
}

real
norm2diff_lr_hmatrix(pchmatrix A, pchmatrix LR)
{
  return norm2diff_pre_matrix((mvm_t) mvm_hmatrix_avector, (void *) A,
			      (prcd_t) lreval_n_hmatrix_avector,
			      (prcd_t) lreval_t_hmatrix_avector, (void *) LR,
			      A->rc->size, A->cc->size);
}

real
norm2diff_chol_hmatrix_amatrix(pchmatrix A, pcamatrix chol)
{
  return norm2diff_pre_matrix((mvm_t) mvm_hmatrix_avector, (void *) A,
			      (prcd_t) choleval_amatrix_avector,
			      (prcd_t) choleval_amatrix_avector,
			      (void *) chol, A->rc->size, A->cc->size);
}

real
norm2diff_chol_hmatrix(pchmatrix A, pchmatrix chol)
{
  return norm2diff_pre_matrix((mvm_t) mvm_hmatrix_avector, (void *) A,
			      (prcd_t) choleval_hmatrix_avector,
			      (prcd_t) choleval_hmatrix_avector,
			      (void *) chol, A->rc->size, A->cc->size);
}

/****************************************************
 * Norm2diff_pre_matrix for sparsematrix
 ****************************************************/

real
norm2diff_lr_sparsematrix_amatrix(pcsparsematrix A, pcamatrix LR)
{
  return norm2diff_pre_matrix((mvm_t) mvm_sparsematrix_avector, (void *) A,
			      (prcd_t) lreval_n_amatrix_avector,
			      (prcd_t) lreval_t_amatrix_avector, (void *) LR,
			      A->rows, A->cols);
}

real
norm2diff_lr_sparsematrix_hmatrix(pcsparsematrix A, pchmatrix LR)
{
  return norm2diff_pre_matrix((mvm_t) mvm_sparsematrix_avector, (void *) A,
			      (prcd_t) lreval_n_hmatrix_avector,
			      (prcd_t) lreval_t_hmatrix_avector, (void *) LR,
			      A->rows, A->cols);
}

real
norm2diff_chol_sparsematrix_amatrix(pcsparsematrix A, pcamatrix chol)
{
  return norm2diff_pre_matrix((mvm_t) mvm_sparsematrix_avector, (void *) A,
			      (prcd_t) choleval_amatrix_avector,
			      (prcd_t) choleval_amatrix_avector,
			      (void *) chol, A->rows, A->cols);
}

real
norm2diff_chol_sparsematrix_hmatrix(pcsparsematrix A, pchmatrix chol)
{
  return norm2diff_pre_matrix((mvm_t) mvm_sparsematrix_avector, (void *) A,
			      (prcd_t) choleval_hmatrix_avector,
			      (prcd_t) choleval_hmatrix_avector,
			      (void *) chol, A->rows, A->cols);
}

/****************************************************
 * Norm2diff_pre_matrix for h2matrix
 ****************************************************/

real
norm2diff_lr_h2matrix_amatrix(pch2matrix A, pcamatrix LR)
{
  return norm2diff_pre_matrix((mvm_t) mvm_h2matrix_avector, (void *) A,
			      (prcd_t) lreval_n_amatrix_avector,
			      (prcd_t) lreval_t_amatrix_avector, (void *) LR,
			      A->rb->t->size, A->cb->t->size);
}

real
norm2diff_lr_h2matrix_hmatrix(pch2matrix A, pchmatrix LR)
{
  return norm2diff_pre_matrix((mvm_t) mvm_h2matrix_avector, (void *) A,
			      (prcd_t) lreval_n_hmatrix_avector,
			      (prcd_t) lreval_t_hmatrix_avector, (void *) LR,
			      A->rb->t->size, A->cb->t->size);
}

real
norm2diff_chol_h2matrix_amatrix(pch2matrix A, pcamatrix chol)
{
  return norm2diff_pre_matrix((mvm_t) mvm_h2matrix_avector, (void *) A,
			      (prcd_t) choleval_amatrix_avector,
			      (prcd_t) choleval_amatrix_avector,
			      (void *) chol, A->rb->t->size, A->cb->t->size);
}

real
norm2diff_chol_h2matrix_hmatrix(pch2matrix A, pchmatrix chol)
{
  return norm2diff_pre_matrix((mvm_t) mvm_h2matrix_avector, (void *) A,
			      (prcd_t) choleval_hmatrix_avector,
			      (prcd_t) choleval_hmatrix_avector,
			      (void *) chol, A->rb->t->size, A->cb->t->size);
}

/****************************************************
 * Norm2diff_pre_matrix for dh2matrix
 ****************************************************/

real
norm2diff_lr_dh2matrix_amatrix(pcdh2matrix A, pcamatrix LR)
{
  return norm2diff_pre_matrix((mvm_t) mvm_dh2matrix_avector, (void *) A,
			      (prcd_t) lreval_n_amatrix_avector,
			      (prcd_t) lreval_t_amatrix_avector, (void *) LR,
			      A->rb->t->size, A->cb->t->size);
}

real
norm2diff_lr_dh2matrix_hmatrix(pcdh2matrix A, pchmatrix LR)
{
  return norm2diff_pre_matrix((mvm_t) mvm_dh2matrix_avector, (void *) A,
			      (prcd_t) lreval_n_hmatrix_avector,
			      (prcd_t) lreval_t_hmatrix_avector, (void *) LR,
			      A->rb->t->size, A->cb->t->size);
}

real
norm2diff_chol_dh2matrix_amatrix(pcdh2matrix A, pcamatrix chol)
{
  return norm2diff_pre_matrix((mvm_t) mvm_dh2matrix_avector, (void *) A,
			      (prcd_t) choleval_amatrix_avector,
			      (prcd_t) choleval_amatrix_avector,
			      (void *) chol, A->rb->t->size, A->cb->t->size);
}

real
norm2diff_chol_dh2matrix_hmatrix(pcdh2matrix A, pchmatrix chol)
{
  return norm2diff_pre_matrix((mvm_t) mvm_dh2matrix_avector, (void *) A,
			      (prcd_t) choleval_hmatrix_avector,
			      (prcd_t) choleval_hmatrix_avector,
			      (void *) chol, A->rb->t->size, A->cb->t->size);
}

/****************************************************
 * Norm2diff_id_pre_matrix for amatrix
 ****************************************************/

real
norm2diff_id_lr_amatrix(pcamatrix A, pcamatrix LR)
{
  return norm2diff_id_pre_matrix((mvm_t) mvm_amatrix_avector, (void *) A,
				 (prcd_t) lrsolve_n_amatrix_avector,
				 (prcd_t) lrsolve_t_amatrix_avector,
				 (void *) LR, A->rows, A->cols);
}

real
norm2diff_id_lr_amatrix_hmatrix(pcamatrix A, pchmatrix LR)
{
  return norm2diff_id_pre_matrix((mvm_t) mvm_amatrix_avector, (void *) A,
				 (prcd_t) lrsolve_n_hmatrix_avector,
				 (prcd_t) lrsolve_t_hmatrix_avector,
				 (void *) LR, A->rows, A->cols);
}

real
norm2diff_id_chol_amatrix(pcamatrix A, pcamatrix chol)
{
  return norm2diff_id_pre_matrix((mvm_t) mvm_amatrix_avector, (void *) A,
				 (prcd_t) cholsolve_amatrix_avector,
				 (prcd_t) cholsolve_amatrix_avector,
				 (void *) chol, A->rows, A->cols);
}

real
norm2diff_id_chol_amatrix_hmatrix(pcamatrix A, pchmatrix chol)
{
  return norm2diff_id_pre_matrix((mvm_t) mvm_amatrix_avector, (void *) A,
				 (prcd_t) cholsolve_hmatrix_avector,
				 (prcd_t) cholsolve_hmatrix_avector,
				 (void *) chol, A->rows, A->cols);
}

/****************************************************
 * Norm2diff_id_pre_matrix for hmatrix
 ****************************************************/

real
norm2diff_id_lr_hmatrix_amatrix(pchmatrix A, pcamatrix LR)
{
  return norm2diff_id_pre_matrix((mvm_t) mvm_hmatrix_avector, (void *) A,
				 (prcd_t) lrsolve_n_amatrix_avector,
				 (prcd_t) lrsolve_t_amatrix_avector,
				 (void *) LR, A->rc->size, A->cc->size);
}

real
norm2diff_id_lr_hmatrix(pchmatrix A, pchmatrix LR)
{
  return norm2diff_id_pre_matrix((mvm_t) mvm_hmatrix_avector, (void *) A,
				 (prcd_t) lrsolve_n_hmatrix_avector,
				 (prcd_t) lrsolve_t_hmatrix_avector,
				 (void *) LR, A->rc->size, A->cc->size);
}

real
norm2diff_id_chol_hmatrix_amatrix(pchmatrix A, pcamatrix chol)
{
  return norm2diff_id_pre_matrix((mvm_t) mvm_hmatrix_avector, (void *) A,
				 (prcd_t) cholsolve_amatrix_avector,
				 (prcd_t) cholsolve_amatrix_avector,
				 (void *) chol, A->rc->size, A->cc->size);
}

real
norm2diff_id_chol_hmatrix(pchmatrix A, pchmatrix chol)
{
  return norm2diff_id_pre_matrix((mvm_t) mvm_hmatrix_avector, (void *) A,
				 (prcd_t) cholsolve_hmatrix_avector,
				 (prcd_t) cholsolve_hmatrix_avector,
				 (void *) chol, A->rc->size, A->cc->size);
}

/****************************************************
 * Norm2diff_id_pre_matrix for sparsematrix
 ****************************************************/

real
norm2diff_id_lr_sparsematrix_amatrix(pcsparsematrix A, pcamatrix LR)
{
  return norm2diff_id_pre_matrix((mvm_t) mvm_sparsematrix_avector, (void *) A,
				 (prcd_t) lrsolve_n_amatrix_avector,
				 (prcd_t) lrsolve_t_amatrix_avector,
				 (void *) LR, A->rows, A->cols);
}

real
norm2diff_id_lr_sparsematrix_hmatrix(pcsparsematrix A, pchmatrix LR)
{
  return norm2diff_id_pre_matrix((mvm_t) mvm_sparsematrix_avector, (void *) A,
				 (prcd_t) lrsolve_n_hmatrix_avector,
				 (prcd_t) lrsolve_t_hmatrix_avector,
				 (void *) LR, A->rows, A->cols);
}

real
norm2diff_id_chol_sparsematrix_amatrix(pcsparsematrix A, pcamatrix chol)
{
  return norm2diff_id_pre_matrix((mvm_t) mvm_sparsematrix_avector, (void *) A,
				 (prcd_t) cholsolve_amatrix_avector,
				 (prcd_t) cholsolve_amatrix_avector,
				 (void *) chol, A->rows, A->cols);
}

real
norm2diff_id_chol_sparsematrix_hmatrix(pcsparsematrix A, pchmatrix chol)
{
  return norm2diff_id_pre_matrix((mvm_t) mvm_sparsematrix_avector, (void *) A,
				 (prcd_t) cholsolve_hmatrix_avector,
				 (prcd_t) cholsolve_hmatrix_avector,
				 (void *) chol, A->rows, A->cols);
}

/****************************************************
 * Norm2diff_id_pre_matrix for h2matrix
 ****************************************************/

real
norm2diff_id_lr_h2matrix_amatrix(pch2matrix A, pcamatrix LR)
{
  return norm2diff_id_pre_matrix((mvm_t) mvm_h2matrix_avector, (void *) A,
				 (prcd_t) lrsolve_n_amatrix_avector,
				 (prcd_t) lrsolve_t_amatrix_avector,
				 (void *) LR, A->rb->t->size, A->cb->t->size);
}

real
norm2diff_id_lr_h2matrix_hmatrix(pch2matrix A, pchmatrix LR)
{
  return norm2diff_id_pre_matrix((mvm_t) mvm_h2matrix_avector, (void *) A,
				 (prcd_t) lrsolve_n_hmatrix_avector,
				 (prcd_t) lrsolve_t_hmatrix_avector,
				 (void *) LR, A->rb->t->size, A->cb->t->size);
}

real
norm2diff_id_chol_h2matrix_amatrix(pch2matrix A, pcamatrix chol)
{
  return norm2diff_id_pre_matrix((mvm_t) mvm_h2matrix_avector, (void *) A,
				 (prcd_t) cholsolve_amatrix_avector,
				 (prcd_t) cholsolve_amatrix_avector,
				 (void *) chol, A->rb->t->size,
				 A->cb->t->size);
}

real
norm2diff_id_chol_h2matrix_hmatrix(pch2matrix A, pchmatrix chol)
{
  return norm2diff_id_pre_matrix((mvm_t) mvm_h2matrix_avector, (void *) A,
				 (prcd_t) cholsolve_hmatrix_avector,
				 (prcd_t) cholsolve_hmatrix_avector,
				 (void *) chol, A->rb->t->size,
				 A->cb->t->size);
}

/****************************************************
 * Norm2diff_id_pre_matrix for dh2matrix
 ****************************************************/

real
norm2diff_id_lr_dh2matrix_amatrix(pcdh2matrix A, pcamatrix LR)
{
  return norm2diff_id_pre_matrix((mvm_t) mvm_dh2matrix_avector, (void *) A,
				 (prcd_t) lrsolve_n_amatrix_avector,
				 (prcd_t) lrsolve_t_amatrix_avector,
				 (void *) LR, A->rb->t->size, A->cb->t->size);
}

real
norm2diff_id_lr_dh2matrix_hmatrix(pcdh2matrix A, pchmatrix LR)
{
  return norm2diff_id_pre_matrix((mvm_t) mvm_dh2matrix_avector, (void *) A,
				 (prcd_t) lrsolve_n_hmatrix_avector,
				 (prcd_t) lrsolve_t_hmatrix_avector,
				 (void *) LR, A->rb->t->size, A->cb->t->size);
}

real
norm2diff_id_chol_dh2matrix_amatrix(pcdh2matrix A, pcamatrix chol)
{
  return norm2diff_id_pre_matrix((mvm_t) mvm_dh2matrix_avector, (void *) A,
				 (prcd_t) cholsolve_amatrix_avector,
				 (prcd_t) cholsolve_amatrix_avector,
				 (void *) chol, A->rb->t->size,
				 A->cb->t->size);
}

real
norm2diff_id_chol_dh2matrix_hmatrix(pcdh2matrix A, pchmatrix chol)
{
  return norm2diff_id_pre_matrix((mvm_t) mvm_dh2matrix_avector, (void *) A,
				 (prcd_t) cholsolve_hmatrix_avector,
				 (prcd_t) cholsolve_hmatrix_avector,
				 (void *) chol, A->rb->t->size,
				 A->cb->t->size);
}
