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
* @file expected_reward.cpp
* @brief Compute the expected time for an MA
* @author Dennis Guck
* @version 1.0
*
*/

#ifndef EXPECTED_REWARD_H
#define EXPECTED_REWARD_H

#include "sparse.h"
#include "real.h"

/**
* Computes expected reward for MA with rewards.
*
* @param ma file to read MA from
* @param max identifier for min or max
* @return expected reward 
*/
extern Real expected_reward_value_iteration(SparseMatrix*, bool);


#endif
