/**
 * Header to define the real type. Prevent the need to include the entirety of the LP solver includes.
 * Should be compatible, hence splitting in both solver types
 */

#ifdef __SOPLEX__

//#include "soplex.h"
//using soplex::infinity;

typedef double Real;
// defined in lp.cpp
extern const Real infinity;

#else

//#include "lp_lib.h"
//typedef REAL Real;

typedef double Real;
// defined in lp.cpp
extern const Real infinity;

#endif
