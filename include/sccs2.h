/**
* IMCA is a analyzing tool for unbounded reachability probabilities, expected-
* time, and long-run averages for Interactive Markov Chains and Markov Automata.
* Copyright (C) RWTH Aachen, 2012
*				UTwente, 2013
* 	Author: Dennis Guck
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program.  If not, see <http://www.gnu.org/licenses/>.
* 
* @file scc2.h
* @author Daan van Beek (s0167789)
* @version 1.0
*
* Created on January 5, 2013, 4:37 PM
*/

#ifndef SCCS2_H
#define	SCCS2_H

#include "sparse.h"
#include <vector>

using std::vector;

/**
 * Compute MECs
 * 
 * @param ma file to read MA from
 */
extern SparseMatrixMEC* mEC_decomposition_previous_algorithm(SparseMatrix*);

extern SparseMatrixMEC* mEC_decomposition_previous_algorithm_without_attractor(SparseMatrix*, vector<unsigned long>&);

extern void compute_SCC_decomposition_tarjan(SparseMatrix *ma, vector<unsigned long>& scc_states, bool* bad_states, bool* bad_transitions, unsigned long& scc_nr);

#endif	/* SCCS2_H */

