/* 
 * File:   sccs2.h
 * Author: mneox
 *
 * Created on January 5, 2013, 4:37 PM
 */

#ifndef SCCS2_H
#define	SCCS2_H

#include "sparse.h"
#include <vector>

#ifdef __SOPLEX__
#include "soplex.h"
#endif

using namespace soplex;


/**
 * Compute MECs
 * 
 * @param ma file to read MA from
 */
extern SparseMatrixMEC* mEC_decomposition_previous_algorithm(SparseMatrix*);

extern SparseMatrixMEC* mEC_decomposition_previous_algorithm_without_attractor(SparseMatrix*, vector<unsigned long>&);

extern void compute_SCC_decomposition_tarjan(SparseMatrix *ma, vector<unsigned long>& scc_states, bool* bad_states, bool* bad_transitions, unsigned long& scc_nr);

#endif	/* SCCS2_H */

