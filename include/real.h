/**
 * Header to define the real type. Prevent the need to include the entirety of the LP solver includes.
 * Should be compatible, hence splitting in both solver types
 */

#define _SOPLEX_ 1
#define _LPSOLVE_ 2

#ifdef __LPSOLVER__
#if __LPSOLVER__==_SOPLEX_

//#include "soplex.h"
//using soplex::infinity;

typedef double Real;
// defined in lp.cpp
extern const Real infinity;

#elif __LPSOLVER__==_LPSOLVE_

//#include "lp_lib.h"
//typedef REAL Real;

typedef double Real;
// defined in lp.cpp
extern const Real infinity;

#endif
#else

typedef double Real;
// defined in lp.cpp
extern const Real infinity;

#endif //__LPSOLVER__
