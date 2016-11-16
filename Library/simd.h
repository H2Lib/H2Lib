#ifndef SIMD_H_
#define SIMD_H_

#ifdef USE_SIMD

#include <immintrin.h>

#ifdef __AVX__
#include "simd_avx.h"
#else
#ifdef __SSE2__
#include "simd_sse2.h"
#else
#ifdef __SSE__
#include "simd_sse.h"
#endif
#endif
#endif

/****************************************************
 * Define vector sizes
 ****************************************************/

#ifdef USE_FLOAT
#define VREAL VFLOAT
#else
#define VREAL VDOUBLE
#endif

#ifdef USE_COMPLEX
#define VFIELD (VREAL/2)
#else
#define VFIELD VREAL
#endif

/****************************************************
 * Define vector types
 ****************************************************/

#ifdef USE_FLOAT
typedef vecf vreal;
#else
typedef vecd vreal;
#endif

typedef vreal vfield;

/****************************************************
 * Define arithmetic operations for real data-type
 ****************************************************/

#ifdef USE_FLOAT
#define vadd vadd_ps
#define vsub vsub_ps
#define vmul vmul_ps
#define vdiv vdiv_ps
#define vsqrt vsqrt_ps
#define vrsqrt vrsqrt_ps
#define vfmadd vfmadd_ps
#define vfmsub vfmsub_ps
#define vfnmadd vfnmadd_ps
#define vfnmsub vfnmsub_ps
#else
#define vadd vadd_pd
#define vsub vsub_pd
#define vmul vmul_pd
#define vdiv vdiv_pd
#define vsqrt vsqrt_pd
#define vrsqrt vrsqrt_pd
#define vfmadd vfmadd_pd
#define vfmsub vfmsub_pd
#define vfnmadd vfnmadd_pd
#define vfnmsub vfnmsub_pd
#endif

/****************************************************
 * Define advanced arithmetic operations for real data-type
 ****************************************************/

#ifdef USE_FLOAT
#define vsin vsin_ps
#define vcos vcos_ps
#define vsincos vsincos_ps
#define vexp vexp_ps
#else
#define vsin vsin_pd
#define vcos vcos_pd
#define vsincos vsincos_pd
#define vexp vexp_pd
#endif

/****************************************************
 * Define load/store operations for real data-type
 ****************************************************/

#ifdef USE_FLOAT
#define vload vload_ps
#define vload1 vload1_ps
#define vloadu vloadu_ps
#define vset1 vset1_ps
#define vsetzero vsetzero_ps
#define vstore vstore_ps
#define vstoreu vstoreu_ps
#else
#define vload vload_pd
#define vload1 vload1_pd
#define vloadu vloadu_pd
#define vset1 vset1_pd
#define vsetzero vsetzero_pd
#define vstore vstore_pd
#define vstoreu vstoreu_pd
#endif

/****************************************************
 * Define compare operations for real data-type
 ****************************************************/

#ifdef USE_FLOAT
#define vcmpeq(a,b) vcmpeq_ps(a,b)
#define vcmpneq(a,b) vcmpneq_ps(a,b)
#define vcmpge(a,b) vcmpge_ps(a,b)
#define vcmpgt(a,b) vcmpgt_ps(a,b)
#define vcmpnge(a,b) vcmpnge_ps(a,b)
#define vcmpngt(a,b) vcmpngt_ps(a,b)
#define vcmple(a,b) vcmple_ps(a,b)
#define vcmplt(a,b) vcmplt_ps(a,b)
#define vcmpnle(a,b) vcmpnle_ps(a,b)
#define vcmpnlt(a,b) vcmpnlt_ps(a,b)
#else
#define vcmpeq(a,b) vcmpeq_pd(a,b)
#define vcmpneq(a,b) vcmpneq_pd(a,b)
#define vcmpge(a,b) vcmpge_pd(a,b)
#define vcmpgt(a,b) vcmpgt_pd(a,b)
#define vcmpnge(a,b) vcmpnge_pd(a,b)
#define vcmpngt(a,b) vcmpngt_pd(a,b)
#define vcmple(a,b) vcmple_pd(a,b)
#define vcmplt(a,b) vcmplt_pd(a,b)
#define vcmpnle(a,b) vcmpnle_pd(a,b)
#define vcmpnlt(a,b) vcmpnlt_pd(a,b)
#endif

/****************************************************
 * Definitions of bit operations for real data-type
 ****************************************************/

#ifdef USE_FLOAT
#define vand vand_ps
#define vandnot vandnot_ps
#define vor vor_ps
#define vxor vxor_ps
#else
#define vand vand_pd
#define vandnot vandnot_pd
#define vor vor_pd
#define vxor vxor_pd
#endif

/****************************************************
 * Define reductions of vector registers for real data-type
 ****************************************************/

#ifdef USE_FLOAT
#define vreduce vreduce_ps
#else
#define vreduce vreduce_pd
#endif

/****************************************************
 * Definition of little helper functions for real-type
 ****************************************************/

#ifdef USE_FLOAT
#define vdot3 vdot3_ps
#else
#define vdot3 vdot3_pd
#endif


#endif

#endif /* SIMD_H_ */
