
/* ---------------------------------------------------------------
 * This is the file "dclusterbasisoperator.c" of the H2Lib package.
 * All rights reserved, Christina Boerst 2015/16
 * --------------------------------------------------------------- */

/** @file dclusteroperator.c
 *  @author Christina Boerst
 */


#include <stdio.h>

#include "basic.h"
#include "dclusteroperator.h"


static uint active_dclusteroperator = 0;

 /*****************************************************************/

pdclusteroperator
init_dclusteroperator(pdclusteroperator co, pcdcluster t)
{

  uint      i, sons;
  sons = t->sons;

  (void) init_leaf_dclusteroperator(co, t);

  if (sons > 0) {
    co->sons = sons;
    co->son =
      (pdclusteroperator *) allocmem((size_t) sizeof(pdclusteroperator) *
				     sons);
    for (i = 0; i < sons; i++) {
      co->son[i] = NULL;
    }
  }

  return co;
}


pdclusteroperator
init_leaf_dclusteroperator(pdclusteroperator co, pcdcluster t)
{

  uint      i;
  /* Something really stupid happened, if a directional cluster has no direction t-> directions is null, but the associated directional cluster basis cb has cb-> directions =1!!  */
  uint      direction;
  if (t->directions == 0) {
    direction = 1;
  }
  else {
    direction = t->directions;
  }

  assert(co != NULL);

  co->t = t;
  co->dir = direction;
  co->krow = (uint *) allocmem(sizeof(uint) * direction);
  co->kcol = (uint *) allocmem(sizeof(uint) * direction);
  co->C = (pamatrix) allocmem(sizeof(amatrix) * direction);
  for (i = 0; i < direction; i++) {
    co->krow[i] = 0;
    co->kcol[i] = 0;
    init_amatrix(co->C + i, 0, 0);
  }
  co->sons = 0;
  co->son = NULL;
  co->refs = 0;

  active_dclusteroperator++;

  return co;
}


void
uninit_dclusteroperator(pdclusteroperator co)
{

  uint      i;

  assert(co->refs == 0);

  if (co->sons > 0) {
    for (i = 0; i < co->sons; i++) {
      unref_dclusteroperator(co->son[i]);
    }
    freemem(co->son);
  }
  for (i = 0; i < co->dir; i++) {
    uninit_amatrix(&co->C[i]);
  }
  freemem(co->krow);
  freemem(co->kcol);
  freemem(co->C);
  assert(active_dclusteroperator > 0);

  active_dclusteroperator--;
}


pdclusteroperator
new_dclusteroperator(pcdcluster t)
{

  pdclusteroperator co;


  co = (pdclusteroperator) allocmem(sizeof(dclusteroperator));
  init_dclusteroperator(co, t);

  return co;
}


pdclusteroperator
new_leaf_dclusteroperator(pcdcluster t)
{

  pdclusteroperator co;


  co = (pdclusteroperator) allocmem(sizeof(dclusteroperator));
  init_leaf_dclusteroperator(co, t);

  return co;
}


void
del_dclusteroperator(pdclusteroperator co)
{

  uninit_dclusteroperator(co);

  freemem(co);
}


void
ref_dclusteroperator(pdclusteroperator * ptr, pdclusteroperator co)
{

  if (*ptr) {
    unref_dclusteroperator(*ptr);
  }

  *ptr = co;

  if (co) {
    co->refs++;
  }
}


void
unref_dclusteroperator(pdclusteroperator co)
{


  assert(co->refs > 0);

  co->refs--;

  if (co->refs == 0) {
    del_dclusteroperator(co);
  }
}

uint
getactives_dclusteroperator()
{

  return active_dclusteroperator;
}


void
resize_dclusteroperator(pdclusteroperator co, uint krow, uint kcol, uint i)
{

  assert(i <= co->t->directions);

  if (krow != co->krow[i] || kcol != co->kcol[i]) {

    resize_amatrix(&co->C[i], krow, kcol);
    co->krow[i] = krow;
    co->kcol[i] = kcol;
  }

}


pdclusteroperator
build_from_dcluster_dclusteroperator(pcdcluster t)
{

  pdclusteroperator co, co1;
  uint      i;

  co = new_dclusteroperator(t);

  for (i = 0; i < co->sons; i++) {
    co1 = build_from_dcluster_dclusteroperator(t->son[i]);
    ref_dclusteroperator(co->son + i, co1);
  }

  return co;
}


pdclusteroperator
build_from_dclusterbasis_dclusteroperator(pcdclusterbasis cb)
{

  pdclusteroperator co, co1;
  uint      i;

  co =
    (cb->son ? new_dclusteroperator(cb->t) :
     new_leaf_dclusteroperator(cb->t));

//        printf("cb->dir %u, co->dir %u\n", cb->directions, co->dir);
  for (i = 0; i < cb->directions; i++) {
    resize_dclusteroperator(co, 0, cb->k[i], i);
  }
  for (i = 0; i < co->sons; i++) {
    co1 = build_from_dclusterbasis_dclusteroperator(cb->son[i]);
    ref_dclusteroperator(co->son + i, co1);
  }

  return co;
}


void
merge_dclusteropertator(pcdclusteroperator co, pdclusteroperator co2)
{

  uint      i, j;
  amatrix   tmp;
  pamatrix  A;

  assert(co->t == co2->t);

  for (i = 0; i < co->sons; i++) {
    merge_dclusteropertator(co->son[i], co2->son[i]);
  }

  for (j = 0; j < co->dir; j++) {
    A = init_amatrix(&tmp, co2->C[j].rows, co2->C[j].cols);
    copy_amatrix(false, &co2->C[j], A);

    assert(co2->C[j].cols == co->C[j].rows);
    resize_dclusteroperator(co2, co2->C[j].rows, co->C[j].cols, j);
    clear_amatrix(&co2->C[j]);
    addmul_amatrix(1.0, false, A, false, &co->C[j], &co2->C[j]);

    uninit_amatrix(A);
  }
}



static void
print_tree(pcdclusteroperator co, uint level)
{

  uint      i;

  printf("level %u\n", level);
  printf("----------------------------------------------\n");

  for (i = 0; i < co->dir; i++) {
    printf("direction: %u Size: %u %u %.4g\n", i, co->krow[i], co->kcol[i],
	   norm2_amatrix(&(co->C[i])));
  }

  for (i = 0; i < co->sons; i++)
    print_tree(co->son[i], level + 1);
}


void
print_tree_dclusteroperator(pcdclusteroperator co)
{

  print_tree(co, 0);
}


/* -------------------------------------------------- */
/*    Enumeration                                     */
/* -------------------------------------------------- */


static void
enumerate(pcdcluster t, uint t_name, pdclusteroperator co,
	  pdclusteroperator * num)
{

  uint      i;
  uint      t_name1;


  assert(co == 0 || co->t == t);
  /* If it is the right dclusteroperator we will name it */

  num[t_name] = co;
  t_name1 = t_name + 1;

  if (co == 0 || co->son == 0) {
    for (i = 0; i < t->sons; i++) {
      enumerate(t->son[i], t_name1, 0, num);
      t_name1 += t->son[i]->desc;
    }
  }
  else {
    assert(t->sons == co->sons);

    for (i = 0; i < t->sons; i++) {
      enumerate(t->son[i], t_name1, co->son[i], num);
      t_name1 += t->son[i]->desc;
    }
  }
  assert(t_name1 == t_name + t->desc);
}


pdclusteroperator *
enumerate_dclusteroperator(pcdcluster t, pdclusteroperator co)
{

  pdclusteroperator *numbered;

  /* struct in form of the given dclusteroperator which contains the numbers */
  numbered = allocmem((size_t) sizeof(pdclusteroperator) * co->t->desc);

  enumerate(t, 0, co, numbered);

  return numbered;
}
