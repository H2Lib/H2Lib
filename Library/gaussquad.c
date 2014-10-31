
/* ------------------------------------------------------------
   This is the file "gaussquad.c" of the H2Lib package.
   All rights reserved, Steffen Boerm 2010
   ------------------------------------------------------------ */

#include "gaussquad.h"

#include "eigensolvers.h"
#include "basic.h"

void
assemble_gauss(uint m, preal x, preal w)
{
  ptridiag  T;
  pamatrix  Q;
  uint      i, info;

  T = new_tridiag(m);
  for (i = 0; i < m; i++)
    T->d[i] = 0.0;
  for (i = 0; i < m - 1; i++)
    T->l[i] = (i + 1.0) / REAL_SQRT((2.0 * i + 1.0) * (2.0 * i + 3.0));

  Q = new_amatrix(m, m);

  info = eig_tridiag(T, Q);
  assert(info == 0);

  for (i = 0; i < m; i++) {
    x[i] = T->d[i];
    w[i] = 2.0 * ABSSQR(getentry_amatrix(Q, 0, i));
  }

  del_amatrix(Q);
  del_tridiag(T);
}
