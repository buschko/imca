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
* @file read_file.cpp
* @brief Reads an compatible MA File
* @author Dennis Guck
* @version 1.0
*
*/

#ifndef READ_FILE_H
#define READ_FILE_H

#include "sparse.h"

#ifdef __SOPLEX__
#include "soplex.h"
using namespace soplex;
#endif

#define MAX_LINE_LENGTH 1024
#define MARKOV_ACTION "!"
#define INITIALS "#INITIALS"
#define GOALS "#GOALS"
#define TRANSITIONS "#TRANSITIONS"

extern void print_model(SparseMatrix*,bool);

extern void print_model_info(SparseMatrix*);


#ifdef __SOPLEX__
/**
* print out information about SoPlex LP solve.
*
* @param lp_model the LP
*/
extern void print_lp_info(SoPlex);
#endif

/**
* Reads MA file @a filename.
*
* @param filename file to read MA from
* @return MA read from file
*/
extern SparseMatrix *read_MA_SparseMatrix_file(const char*,bool mrm);
extern void witeToDot(SparseMatrix* ma, std::ostream& outStream);


#endif
