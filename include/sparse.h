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
* @file sparse.cpp
* @brief definition of sparse matrix
* @author Dennis Guck
* @version 1.0
* Note: modified version of MRMCs NDSparseMatrix
*	
*/

#ifndef SPARSE_H
#define SPARSE_H

#include <map>
#include <string>

#include "lp.h"

using std::map;
using std::string;

typedef struct SparseMatrix SparseMatrix;

/**
* The SparseMatrix structure for MAs.
*/
struct SparseMatrix
{
	unsigned long n;			/* # of states */
	unsigned long ms_n;			/* # of Markovian states */
	unsigned long choices_n;	/* # of choices */
	unsigned long non_zero_n;	/* # of transitions */
	map<string,unsigned long> states;	/* state names mapped to # */
	map<unsigned long,string> states_nr;	/* state # mapped to names */
	bool *initials;				/* initial states = true, otherwise false */
	bool *goals;				/* goal states = true, otherwise false */
	bool *isPS;				/* probabilistic state = true, otherwise false */
	
	Real *non_zeros;			/* transition probabilities */
	Real *exit_rates;			/* exit rates for Markovian states */
	Real *rewards;			/* rewards for Markovian states and probabilistic transitions */
	Real max_exit_rate;			/* max exit rate */
	unsigned long *cols;			/* successors */
	unsigned char *row_counts;		/* pointer to state */
	unsigned char *rate_counts;		/* pointer to exit rates */
	unsigned char *choice_counts;		/* pointer to distribution */
};

/**
* The SparseMatrix structure for MECs.
*/
struct SparseMatrixMEC
{
	unsigned long n;				/* # of MECs */
	
	unsigned long *cols;				/* MEC states */
	unsigned char *row_counts;		/* pointer to MEC states */
};

extern SparseMatrix* SparseMatrix_new(unsigned long, map<string,unsigned long>, map<unsigned long,string>);
extern SparseMatrixMEC* SparseMatrixMEC_new(unsigned long, unsigned long);
extern SparseMatrix* SparseMatrixDiscrete_new(SparseMatrix* ma);

extern void SparseMatrix_free(SparseMatrix *);
extern void SparseMatrixMEC_free(SparseMatrixMEC *);

#endif
