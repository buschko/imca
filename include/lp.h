/*
 * lp.h
 *
 * Interface for LP solvers
 *
 *  Created on: Jul 17, 2013
 *      Author: harold
 */

#ifndef LP_H_
#define LP_H_

// choose between SoPlex and LP_solve
#ifdef __SOPLEX__

#include "soplex.h"
using namespace soplex;

#else

// have to include gmp for REAL arithmetic
#include "gmp.h"
#include <math.h>
// define location in Makefile
#include "lp_lib.h"

#endif

#endif /* LP_H_ */
