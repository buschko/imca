/* 
 * File:   sccs2.h
 * Author: mneox
 *
 * Created on January 5, 2013, 4:37 PM
 */

#ifndef SCCS2_H
#define	SCCS2_H

#include "sparse.h"
#include "soplex.h"
#include <vector>

using namespace soplex;


/**
 * Compute MECs
 * 
 * @param ma file to read MA from
 */
extern SparseMatrixMEC* mEC_decomposition_previous_algorithm(SparseMatrix*);

extern SparseMatrixMEC* mEC_decomposition_previous_algorithm_without_attractor(SparseMatrix*, vector<unsigned long>&);

#endif	/* SCCS2_H */

