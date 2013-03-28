/**
* IMCA is a analyzing tool for unbounded reachability probabilities, expected-
* time, and long-run averages for Interactive Markov Chains and Markov Automata.
* Copyright (C) RWTH Aachen, 2012
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
* 
* @file sccs.cpp
* @brief Strongly connected components functions to process an MA
* @author Dennis Guck
* @version 1.0
*
*/

#ifndef SCCS_H
#define SCCS_H

#include "sparse.h"

#ifdef __SOPLEX__
#include "soplex.h"
#endif

using namespace soplex;

/**
 * Compute BSCCs with respect to a set of bad states
 * 
  * @param ma file to read MA from
 */
extern SparseMatrixMEC* compute_bottom_strongly_connected_components(SparseMatrix*);

/**
 * Compute MECs
 * 
 * @param ma file to read MA from
 */
extern SparseMatrixMEC* compute_maximal_end_components(SparseMatrix*);

/**
* Computes locks for unbounded reachability for MA.
*
* @param locks identifier if state can't reach a goal state
* @param ma file to read MA from
*/
extern bool* compute_locks_strong(SparseMatrix*);

/**
* Computes locks expected time for MA.
*
* @param locks identifier if state can't reach a goal state
* @param ma file to read MA from
*/
extern bool* compute_locks_weak(SparseMatrix*);

/**
* Computes locks for unbounded reachability for MA.
*
* @param locks identifier if state can't reach a goal state
* @param ma file to read MA from
*/
extern bool* compute_locks_strong(SparseMatrix*, bool *bad);

/**
* Computes locks expected time for MA.
*
* @param locks identifier if state can't reach a goal state
* @param ma file to read MA from
*/
extern bool* compute_locks_weak(SparseMatrix*, bool *bad);


#endif