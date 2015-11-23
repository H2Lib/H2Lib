/* ------------------------------------------------------------
 This is the file "truncation.c" of the H2Lib package.
 All rights reserved, Steffen Boerm 2014
 ------------------------------------------------------------ */

#include "truncation.h"
#include "basic.h"

/* ------------------------------------------------------------
 Constructors for standard truncation strategies
 ------------------------------------------------------------ */

ptruncmode
new_truncmode()
{
  ptruncmode tm;

  tm = (ptruncmode) allocmem(sizeof(truncmode));

  tm->frobenius = false;
  tm->absolute = false;
  tm->blocks = false;
  tm->zeta_level = 1.0;
  tm->zeta_age = 1.0;

  return tm;
}

void
del_truncmode(ptruncmode tm)
{
  freemem(tm);
}

ptruncmode
new_releucl_truncmode()
{
  ptruncmode tm;

  tm = new_truncmode();

  return tm;
}

ptruncmode
new_relfrob_truncmode()
{
  ptruncmode tm;

  tm = new_truncmode();

  tm->frobenius = true;

  return tm;
}

ptruncmode
new_blockreleucl_truncmode()
{
  ptruncmode tm;

  tm = new_truncmode();

  tm->blocks = true;
  tm->absolute = true;
  tm->zeta_age = 1.5;

  return tm;
}

ptruncmode
new_blockrelfrob_truncmode()
{
  ptruncmode tm;

  tm = new_truncmode();

  tm->frobenius = true;
  tm->absolute = true;
  tm->blocks = true;
  tm->zeta_age = 1.5;

  return tm;
}

ptruncmode
new_abseucl_truncmode()
{
  ptruncmode tm;

  tm = new_truncmode();

  tm->absolute = true;
  tm->zeta_level = 1.2;
  tm->zeta_age = 1.5;

  return tm;
}

/* ------------------------------------------------------------
 Find minimal acceptable rank
 ------------------------------------------------------------ */

uint
findrank_truncmode(pctruncmode tm, real eps, pcrealavector sigma)
{
  real      norm, sum;
  uint      k;

  if (tm && tm->frobenius) {
    if (tm->absolute) {
      sum = 0.0;
      k = sigma->dim;
      while (k > 0 && (sum += ABSSQR(sigma->v[k - 1])) <= eps * eps)
	k--;
    }
    else {
      norm = 0.0;
      for (k = 0; k < sigma->dim; k++)
	norm += ABSSQR(sigma->v[k]);

      sum = 0.0;
      k = sigma->dim;
      while (k > 0 && (sum += ABSSQR(sigma->v[k - 1])) <= eps * eps * norm)
	k--;
    }
  }
  else {
    if (tm && tm->absolute) {
      k = 0;
      while (k < sigma->dim && sigma->v[k] > eps)
	k++;
    }
    else {
      k = 0;
      while (k < sigma->dim && sigma->v[k] > eps * sigma->v[0])
	k++;
    }
  }

  return k;
}
