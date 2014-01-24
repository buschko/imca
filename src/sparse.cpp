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
* Source description: 
*	Definition of sparse matrix
*/

#include "sparse.h"

#include <stdlib.h>

#include "debug.h"

/**
* Creates a new MA with given number of states.
*
* @param num_states number of states of new MA
* @return new MA
*/
SparseMatrix *SparseMatrix_new(unsigned long num_states, map<string,unsigned long> states, map<unsigned long,string> states_nr)
{
	unsigned long *row_starts = (unsigned long *) calloc((size_t) (num_states + 1),sizeof(unsigned long));
	unsigned long *rate_starts = (unsigned long *) calloc((size_t) (num_states + 1),sizeof(unsigned long));
	bool * initials = (bool *) malloc(num_states * sizeof(bool));
	bool * goals = (bool *) malloc(num_states * sizeof(bool));
	bool * isPS = (bool *) malloc(num_states * sizeof(bool));
	//SparseMatrix *model = (SparseMatrix*)malloc(sizeof(SparseMatrix));
    SparseMatrix *model=new SparseMatrix;
	model->n = num_states;
	//cout << "" ; // This line currently prevents a segmentation fault under OSX. Do not delete it! TODO: fix this.
    model->states = states;
	model->states_nr = states_nr;
    //model->states.insert(states.begin(),states.end());
    //model->states_nr.insert(states_nr.begin(),states_nr.end());
	model->initials = initials;
	model->goals = goals;
	model->isPS = isPS;
	model->rewards = NULL;
	model->non_zeros = NULL;
	model->exit_rates = NULL;
	model->row_counts = (unsigned char *) row_starts;
	model->rate_counts = (unsigned char *) rate_starts;
	model->cols = NULL;
	model->choice_counts = NULL;
	return model;
}

/**
* Creates a new MEC with given number of states.
*
* @param num_states number of states of new MA
* @return new MEC
*/
SparseMatrixMEC *SparseMatrixMEC_new(unsigned long num_states, unsigned long num_mecs)
{
	unsigned long *row_starts = (unsigned long *) calloc((size_t) (num_mecs + 1),sizeof(unsigned long));
	unsigned long * cols = (unsigned long *) malloc(num_states * sizeof(unsigned long));
	//SparseMatrixMEC *mec = (SparseMatrixMEC*)malloc(sizeof(SparseMatrixMEC));
    SparseMatrixMEC *mec =new SparseMatrixMEC;
	mec->n = num_mecs;
	mec->row_counts = (unsigned char *) row_starts;
	mec->cols = cols;
	return mec;
}

/**
* Creates a new MA with given number of states.
*
* @param num_states number of states of new MA
* @return new MA
*/
SparseMatrix *SparseMatrixDiscrete_new(SparseMatrix *ma)
{
	// copy SparseMatrix, except for non_zeros & cols
	unsigned long num_ms=ma->ms_n;
	unsigned long num_states = ma->n;
	unsigned long num_choices = ma->choices_n;
	unsigned long *row_starts = (unsigned long *) calloc((size_t) (num_states + 1),sizeof(unsigned long));
	unsigned long *rate_starts = (unsigned long *) calloc((size_t) (num_states + 1),sizeof(unsigned long));
	unsigned long * choice_starts = (unsigned long *) calloc((size_t) (num_choices + 1), sizeof(unsigned long));
	bool * initials = (bool *) malloc(num_states * sizeof(bool));
	bool * goals = (bool *) malloc(num_states * sizeof(bool));
	bool * isPS = (bool *) malloc(num_states * sizeof(bool));
	Real * exit_rates = (Real *) malloc(num_ms * sizeof(Real));
	//SparseMatrix *model = (SparseMatrix*)malloc(sizeof(SparseMatrix));
    SparseMatrix *model=new SparseMatrix;
	model->n = num_states;
	model->ms_n = num_ms;
	model->choices_n = num_choices;
	//model->states = ma->states;
	//model->states_nr = ma->states_nr;
	model->initials = initials;
	model->goals = goals;
	model->isPS = isPS;
	model->exit_rates = exit_rates;
	model->max_exit_rate=ma->max_exit_rate;
	model->row_counts = (unsigned char *) row_starts;
	model->rate_counts = (unsigned char *) rate_starts;
	model->choice_counts = (unsigned char *) choice_starts;
	model->rewards = NULL;
	
	// values to assign
	initials = (bool *) model->initials;
	goals = (bool *)model->goals;
	isPS = (bool *)model->isPS;
	exit_rates = model->exit_rates;
	rate_starts = (unsigned long *) model->rate_counts;
	row_starts = (unsigned long *) model->row_counts;
	// help values from MA
	unsigned long *ma_row_starts = (unsigned long *) ma->row_counts;
	//unsigned long *ma_rate_starts = (unsigned long *) ma->rate_counts;
	// copy values
	for(unsigned long state_nr=0; state_nr < num_states; state_nr++){
		// copy bool values for initials, goals and isPS
		initials[state_nr] = ma->initials[state_nr];
		goals[state_nr] = ma->goals[state_nr];
		isPS[state_nr] = ma->isPS[state_nr];
		//isPS[state_nr]=true;
		// copy row and rate values
		unsigned long state_start = ma_row_starts[state_nr];
		row_starts[state_nr] = state_start;
		//unsigned long rate_start = ma_rate_starts[state_nr];
		//rate_starts[state_nr] = rate_start;
		rate_starts[state_nr] = 0;
	}
	// copy row and rate values for last state
	unsigned long state_start = ma_row_starts[num_states];
	row_starts[num_states] = state_start;
	//unsigned long rate_start = ma_rate_starts[num_states];
	//rate_starts[num_states] = rate_start;
	rate_starts[num_states] = 0;
	// copy exit rates for Markovian states
	for(unsigned long state_nr=0; state_nr < num_ms; state_nr++){
		exit_rates[state_nr]=ma->exit_rates[state_nr];
	}
	
	// detect new number of non_zeros and choice_counts entries
	unsigned long num_non_zeros = ma->non_zero_n + num_ms;
	Real * non_zeros = (Real *) malloc(num_non_zeros * sizeof(Real));
	unsigned long * cols = (unsigned long *) malloc(num_non_zeros * sizeof(unsigned long));
	Real * rewards = (Real *) malloc(num_choices * sizeof(Real));
	model->non_zeros = non_zeros;
	model->rewards = rewards;
	model->cols = cols;
	model->non_zero_n = num_non_zeros;
	choice_starts = (unsigned long *) model->choice_counts;
	int offset=0;
	unsigned long *ma_choice_starts = (unsigned long *) ma->choice_counts;
	for (unsigned long state_nr = 0; state_nr < ma->n; state_nr++) {
		unsigned long state_start = ma_row_starts[state_nr];
		unsigned long state_end = ma_row_starts[state_nr + 1];
		// Look at Markovian states
		if(!ma->isPS[state_nr])
		{
			for (unsigned long choice_nr = state_start; choice_nr < state_end; choice_nr++) {
				unsigned long i_start = ma_choice_starts[choice_nr];
				unsigned long i_end = ma_choice_starts[choice_nr + 1];
				bool loop=false;
				for (unsigned long i = i_start; i < i_end; i++) {
					if(state_nr == ma -> cols[i])
						loop=true;
				}
				if(!loop) {
					offset++;
					choice_starts[choice_nr] = i_start+offset-1;
				}
				else {
					choice_starts[choice_nr] = i_start+offset;
				}
				choice_starts[choice_nr+1] = i_start+offset;
			}
		}else {
			for (unsigned long choice_nr = state_start; choice_nr < state_end; choice_nr++) {
				unsigned long i_start = ma_choice_starts[choice_nr];
				//unsigned long i_end = ma_choice_starts[choice_nr + 1];
				choice_starts[choice_nr] = i_start+offset;
				choice_starts[choice_nr+1] = i_start+offset;
			}
		}
	}
	choice_starts[num_choices] = ma_choice_starts[num_choices] + offset;
	
	return model;
}


/**
* This function frees the sparse matrix.
* @param sparse the matrix to be freed
*/
void SparseMatrix_free(SparseMatrix *sparse) {
	if (sparse == NULL) {
		dbg_printf("Sparse is free\n");
		return;
	}
	if (sparse->initials != NULL) {
		dbg_printf("free initials\n");
		free(sparse->initials);
	}
	if (sparse->goals != NULL) {
		dbg_printf("free goals\n");
		free(sparse->goals);
	}
	if (sparse->isPS != NULL) {
		dbg_printf("free isPS\n");
		free(sparse->isPS);
	}
	if (sparse->non_zeros != NULL) {
		dbg_printf("free non zeros\n");
		free(sparse->non_zeros);
	}
	if (sparse->rewards != NULL) {
		dbg_printf("free rewards\n");
		free(sparse->rewards);
	}
	if (sparse->exit_rates != NULL) {
		dbg_printf("free exit rates\n");
		free(sparse->exit_rates);
	}
	if (sparse->cols != NULL) {
		dbg_printf("free cols\n");
		free(sparse->cols);
	}
	if (NULL != sparse->row_counts) {
		dbg_printf("free row counts\n");
		free(sparse->row_counts);
	}
	if (NULL != sparse->rate_counts) {
		dbg_printf("free rate counts\n");
		free(sparse->rate_counts);
	}
	if (NULL != sparse->choice_counts) {
		dbg_printf("free choice counts\n");
		free(sparse->choice_counts);
	}
    
    sparse->states.clear();
    sparse->states_nr.clear();
	
	//dbg_printf("free sparse matrix\n");
	//free(sparse);
}

/**
* This function frees the sparse matrix.
* @param sparse the matrix to be freed
*/
void SparseMatrixMEC_free(SparseMatrixMEC *sparse) {
	if (sparse == NULL) {
		return;
	}
	if (sparse->cols != NULL) {
		free(sparse->cols);
	}
	if (NULL != sparse->row_counts) {
		free(sparse->row_counts);
	}
	
	free(sparse);
}
