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
* 
* @file unbounded.cpp
* @brief Compute the unbounded reachability for an MA
* @author Dennis Guck
* @version 1.0
*
*/

#ifndef UNBOUNDED_H
#define UNBOUNDED_H

#include "sparse.h"

#ifdef __SOPLEX__
#include "soplex.h"
#endif

using namespace soplex;

/**
* Computes unbounded reachability for MA.
*
* @param ma file to read MA from
* @param max identifier for min or max prb
* @return probability for unbounded reachability
*/
extern Real compute_unbounded_reachability(SparseMatrix*, bool);

extern Real unbounded_value_iteration(SparseMatrix*, bool);


#endif
