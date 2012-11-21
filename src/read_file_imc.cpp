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
* Source description: 
*	Reads an compatible IMC File
*/

#include "read_file_imc.h"
#include "read_file.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <map>
#include <string>

using namespace std;

/**
* Read number of states and create a hash table for state names.
*
* @param line_no line number in file (should be 1 before)
* @param error error flag
* @param p MA file
* @param filename filename of MA file
* @param num_states will store number of states here
* @param states will map the state names to state numbers 
*/
static void read_states(unsigned long *line_no, bool *error, FILE *p, const char *filename, unsigned long *num_states, map<string,unsigned long> *states, map<unsigned long,string> *states_nr) 
{
	char s[MAX_LINE_LENGTH];
	unsigned long state_nr = 0;
	pair<map<string,unsigned long>::iterator,bool> ret;
	char state[MAX_LINE_LENGTH];
	char act[MAX_LINE_LENGTH];
	map<string ,unsigned long> tmp;
	map<unsigned long, string> tmp_nr;

	/* go to Transitions */
	if (!*error) {
		if(fgets(s, MAX_LINE_LENGTH, p) == 0)
		{
			fprintf(stderr, "Reading line %ld of file \"%s\" failed.\n",*line_no, filename);
			*error = true;
		}
		++(*line_no);
		sscanf(s, "%s", state);
		while(strcmp(state,TRANSITIONS) != 0 && !*error) {
			if(fgets(s, MAX_LINE_LENGTH, p) == 0)
			{
				fprintf(stderr, "Reading line %ld of file \"%s\" failed.\n",*line_no, filename);
				*error = true;
			}
			sscanf(s, "%s", state);
			++(*line_no);
		}
	}
	
	
	/* TODO: also check lines */
	while (fgets(s, MAX_LINE_LENGTH, p) != NULL && !*error) {
		sscanf(s, "%s%s", state, act);
		ret=tmp.insert(pair<string,unsigned long>(state,state_nr));
		if(ret.second == true){
			tmp_nr.insert(pair<unsigned long,string>(state_nr,state));
			state_nr++;
		} 
		++(*line_no);
	}
	
	(*num_states) = state_nr;	
	(*states)=tmp;
	(*states_nr)=tmp_nr;
}

/**
* Check for deadlocks and add selfloop.
*
* @param line_no line number in file (should be 1 before)
* @param error error flag
* @param p MA file
* @param filename filename of MA file
* @param num_states will store number of states here
* @param states will map the state names to state numbers 
*/
static void check_dedlocks(unsigned long *line_no, bool *error, FILE *p, const char *filename, unsigned long *num_states, map<string,unsigned long> *states, map<unsigned long,string> *states_nr) 
{
	char s[MAX_LINE_LENGTH];
	unsigned long state_nr = (*num_states);
	pair<map<string,unsigned long>::iterator,bool> ret;
	char state[MAX_LINE_LENGTH];
	char star[MAX_LINE_LENGTH];
	map<string ,unsigned long> tmp=*states;
	map<unsigned long, string> tmp_nr=*states_nr;

	/* go to Transitions */
	if (!*error) {
		if(fgets(s, MAX_LINE_LENGTH, p) == 0)
		{
			fprintf(stderr, "Reading line %ld of file \"%s\" failed.\n",*line_no, filename);
			*error = true;
		}
		++(*line_no);
		sscanf(s, "%s", state);
		while(strcmp(state,TRANSITIONS) != 0 && !*error) {
			if(fgets(s, MAX_LINE_LENGTH, p) == 0)
			{
				fprintf(stderr, "Reading line %ld of file \"%s\" failed.\n",*line_no, filename);
				*error = true;
			}
			sscanf(s, "%s", state);
			++(*line_no);
		}
	}
	
	
	/* TODO: also check lines */
	while (fgets(s, MAX_LINE_LENGTH, p) != NULL && !*error) {
		if (s[0] == '*') {
			sscanf(s, "%s%s", star, state);
			ret=tmp.insert(pair<string,unsigned long>(state,state_nr));
			if(ret.second == true){
				tmp_nr.insert(pair<unsigned long,string>(state_nr,state));
				state_nr++;
			} 
			++(*line_no);
		} else {
			++(*line_no);
		}
	}
	
	(*num_states) = state_nr;	
	(*states)=tmp;
	(*states_nr)=tmp_nr;
}

static void init_states(unsigned long *line_no, bool *error, FILE *p, const char *filename, unsigned long *num_ms_states, map<string,unsigned long> states, SparseMatrix *ma) {
	char src[MAX_LINE_LENGTH];
	char dst[MAX_LINE_LENGTH];
	char act[MAX_LINE_LENGTH];
	char s[MAX_LINE_LENGTH];
	unsigned long from;
	unsigned long last_from = 0;
	unsigned long ms_states=0;
	bool *isPS;
	bool isMs=false;
	Real rate;
	
	if (!*error) {
		bool *initials = (bool *) ma->initials;
		bool *goals = (bool *)ma->goals;
		isPS = (bool *)ma->isPS;

		for(from=0; from < ma->n; from++){
			/* iinitialize init and goal and isPS with false */
			initials[from]=false;
			goals[from]=false;
			isPS[from]=false;
		}
		
		/* store initials and goals */
		if(fgets(s, MAX_LINE_LENGTH, p) == 0)
		{
			fprintf(stderr, "Reading line %ld of file \"%s\" failed.\n",*line_no, filename);
			*error = true;
		}
		++(*line_no);
		sscanf(s, "%s", src);
		if(strcmp(src,INITIALS) != 0 && !*error) {
			fprintf(stderr,"No declaration of initial states.\n");
			*error = true;
		} else {
			++(*line_no);
			if(fgets(s, MAX_LINE_LENGTH, p) == 0)
			{
				fprintf(stderr, "Reading line %ld of file \"%s\" failed.\n",*line_no, filename);
				*error = true;
			}
			sscanf(s, "%s", src);
			while(strcmp(src,GOALS) != 0 && !*error)
			{
				from=states.find(src)->second;
				initials[from]=true;
				if(fgets(s, MAX_LINE_LENGTH, p) == 0)
				{
					fprintf(stderr, "Reading line %ld of file \"%s\" failed.\n",*line_no, filename);
					*error = true;
				}
				sscanf(s, "%s", src);
				++(*line_no);
			}
			if(fgets(s, MAX_LINE_LENGTH, p) == 0)
			{
				fprintf(stderr, "Reading line %ld of file \"%s\" failed.\n",*line_no, filename);
				*error = true;
			}
			sscanf(s, "%s", src);
			while(strcmp(src,TRANSITIONS) != 0 && !*error) {
				from=states.find(src)->second;
				goals[from]=true;
				if(fgets(s, MAX_LINE_LENGTH, p) == 0)
				{
					fprintf(stderr, "Reading line %ld of file \"%s\" failed.\n",*line_no, filename);
					*error = true;
				}
				sscanf(s, "%s", src);
				++(*line_no);
			}
		}
	}
	
	while (!*error && (fgets(s, MAX_LINE_LENGTH, p) != NULL)) {
		isMs=false;
		if (sscanf(s, "%s%s%lf", src, dst, &rate) != 3) { /* check if all values are returned */
			if(sscanf(s, "%s%s%s", src, dst, act) != 3)
				*error = true;
		} else {
			isMs=true;
		}
		if(!*error) {
			from=states.find(src)->second;
			/* test consistency of ma file */
			if (from < last_from) {
				fprintf(stderr, "Line %ld: State transitions must be given in continous order for one state.\n", *line_no);
				*error = true;		
			} else {
				if(!isMs)
					isPS[from]=true;
				
				last_from = from;
			}
		}
		++(*line_no);
	}
	
	if(!*error) {
		for(from=0; from<ma->n; from++) {
			if(!isPS[from])
				ms_states++;
		}
	}
	
	if(!*error) {
		(*num_ms_states) = ms_states;
		ma->ms_n = ms_states;
		Real * exit_rates = (Real *) malloc(ms_states * sizeof(Real));
		ma->exit_rates = exit_rates;
	}
}


/**
* Allocate memory needed for MA data structures and stores rates for Markovian states.
*
* @param line_no line number in file
* @param error error flag
* @param p MA file
* @param model MA for which allocations will be done
*/
static void reserve_transition_memory(unsigned long *line_no, bool *error, FILE *p, const char *filename, SparseMatrix *model) {
	char s[MAX_LINE_LENGTH];
	char src[MAX_LINE_LENGTH];
	char dst[MAX_LINE_LENGTH];
	char act[MAX_LINE_LENGTH];
	Real exit_rate=0;
	unsigned long exit_index = 0;
	unsigned long num_choice = 0;
	unsigned long num_non_zeros = 0;
	bool isMs = false;
	bool isMs_last = false;
	Real *exit_rates;
	unsigned long *rate_starts;
	unsigned long from;
	unsigned long last_from = 0;
	map<string ,unsigned long> states;
	bool *isPS;
	Real rate;
	bool first=true;

	/* go to Transitions */
	if (!*error) {
		if(fgets(s, MAX_LINE_LENGTH, p) == 0)
		{
			fprintf(stderr, "Reading line %ld of file \"%s\" failed.\n",*line_no, filename);
			*error = true;
		}
		++(*line_no);
		sscanf(s, "%s", src);
		while(strcmp(src,TRANSITIONS) != 0 && !*error) {
			if(fgets(s, MAX_LINE_LENGTH, p) == 0)
			{
				fprintf(stderr, "Reading line %ld of file \"%s\" failed.\n",*line_no, filename);
				*error = true;
			}
			sscanf(s, "%s", src);
			++(*line_no);
		}
	}

	if (!*error) {
		exit_rates = model->exit_rates;
		rate_starts = (unsigned long *) model->rate_counts;
		states = model->states;
		isPS = model->isPS;
	}
	

	while (!*error && (fgets(s, MAX_LINE_LENGTH, p) != NULL)) {
		isMs=false;
		if (sscanf(s, "%s%s%lf", src, dst, &rate) != 3) { /* check if all values are returned */
			if (sscanf(s, "%s%s%s", src, dst, act) != 3)
				*error = true;
		} else {
			isMs = true;
		}
		if(!*error)
		{
			from=states.find(src)->second;
			/* test consistency of ma file */
			if (from < last_from) {
				fprintf(stderr, "Line %ld: State transitions must be given in continous order for one state.\n", *line_no);
				*error = true;
			} else if(from == last_from){	
				if(isMs && !isPS[from]) {
					exit_rate += rate;
				}
				
				if(isPS[from]) {
					num_choice++;
					first=false;
				} else if(first) {
					num_choice++;
					first=false;
				}
				
				isMs_last = isMs;
				last_from = from;
				/* probabilistic transitions are choosen before markovian transitions */
				if((isPS[from] && !isMs) || (!isPS[from] && isMs)) {
					num_non_zeros++;
				}
			} else if(from > last_from) {
				/* if Markovian state, store exit rate */
				if(isMs_last && !isPS[last_from]) {
					exit_rates[exit_index] = exit_rate;
					exit_index++;
					for (; last_from < from; last_from++) {
						rate_starts[last_from + 1] = rate_starts[last_from + 0];
					}
					rate_starts[from + 0]++;
				} else {
					for (; last_from < from; last_from++) {
						rate_starts[last_from + 1] = rate_starts[last_from + 0];
					}
				}
				
				exit_rate=0;
				isMs_last = isMs;
				
				if(isMs && !isPS[from]) {
					exit_rate += rate;
				}
				
				/* probabilistic transitions are choosen before markovian transitions */
				if((isPS[from] && !isMs) || (!isPS[from] && isMs)) {
					num_non_zeros++;
					num_choice++;
				}
				
				last_from = from;
			}
		}
		//++(*line_no);
		++(*line_no);
	}
	/* probabilistic transitions are choosen before markovian transitions */
	if((isMs)) {
		exit_rates[exit_index] = exit_rate;
		for (; from+0 < model->n; from++) {
			rate_starts[from+1] = rate_starts[from+0];
		}
		rate_starts[from + 0]++;
	} else {
		for (; last_from < from; last_from++) {
			rate_starts[last_from + 1] = rate_starts[last_from + 0];
		}
	}
	
	/* now allocate the memory needed */
	if (!*error) {
		unsigned long * choice_starts = (unsigned long *) calloc((size_t) (num_choice + 1), sizeof(unsigned long));
		Real * non_zeros = (Real *) malloc(num_non_zeros * sizeof(Real));
		unsigned long * cols = (unsigned long *) malloc(num_non_zeros * sizeof(unsigned long));
		model->choice_counts = (unsigned char *) choice_starts;
		model->non_zeros = non_zeros;
		model->cols = cols;
	}
}

/**
* Reads the transitions of a MA.
*
* @param line_no line number in file (should be 1 before)
* @param error error flag
* @param p MA file
* @param filename of MA file
* @param ma the transitions shall be added
*/
static void read_transitions(unsigned long *line_no, bool *error, FILE *p, const char *filename, SparseMatrix *ma) {
	unsigned long choice_index = 0;
	unsigned long choice_size = 0;
	//unsigned long choice_size_old = 0;
	unsigned long nz_index = 0;
	char src[MAX_LINE_LENGTH];
	unsigned long from;
	unsigned long last_from = 0;
	char s[MAX_LINE_LENGTH];
	char dst[MAX_LINE_LENGTH];
	unsigned long to;
	char act[MAX_LINE_LENGTH];
	bool bad=false;
	bool first=true;
	Real rate;

	if (!*error) {
		Real *non_zeros = ma->non_zeros;
		unsigned long *cols = ma->cols;
		unsigned long *choice_starts = (unsigned long *) ma->choice_counts;
		unsigned long *row_starts = (unsigned long *) ma->row_counts;
		map<string ,unsigned long> states = ma->states;
		bool *isPS = ma->isPS;
		
		/* go to Transitions */
		if (!*error) {
			if(fgets(s, MAX_LINE_LENGTH, p) == 0)
			{
				fprintf(stderr, "Reading line %ld of file \"%s\" failed.\n",*line_no, filename);
				*error = true;
			}
			++(*line_no);
			sscanf(s, "%s", src);
			while(strcmp(src,TRANSITIONS) != 0 && !*error) {
				if(fgets(s, MAX_LINE_LENGTH, p) == 0)
				{
					fprintf(stderr, "Reading line %ld of file \"%s\" failed.\n",*line_no, filename);
					*error = true;
				}
				sscanf(s, "%s", src);
				++(*line_no);
			}
		}
		

		while (fgets(s, MAX_LINE_LENGTH, p) != NULL) {
			bad=false;
			if (sscanf(s, "%s%s%lf", src, dst, &rate) != 3) { /* check if all values are returned */
				if (sscanf(s, "%s%s%s", src, dst, act) != 3)
					*error = true;
			} else {
				int tmp=states.find(src)->second;
				/* probabilistic transitions are choosen before markovian transitions */
				if(isPS[tmp])
					bad=true;
			}
			if(!bad && !*error)
			{
				from=states.find(src)->second;
				if(from==last_from) {
					to=states.find(dst)->second;
					if(!isPS[from]) {
						non_zeros[nz_index] = rate;
						choice_size++;
						if(first) {
							row_starts[from + 1]++;
							choice_index++;
							first=false;
						}
					} else {
						first=false;
						non_zeros[nz_index] = 1;
						choice_size=1;
						if (choice_index > 0) {
							choice_starts[choice_index] =
							choice_starts[choice_index - 1] + choice_size;
						}
						choice_size = 0;
						choice_index++;
						row_starts[from + 1]++;
					}
					cols[nz_index] = to;
					nz_index++;
				} else if(from > last_from) {
					if(!isPS[last_from]) {
						if (choice_index > 0) {
							choice_starts[choice_index] =
							choice_starts[choice_index - 1] + choice_size;
						}
						choice_size = 0;
						choice_index++;
					} else {
						choice_size=1;
						if (choice_index > 0) {
							choice_starts[choice_index] =
							choice_starts[choice_index - 1] + choice_size;
						}
						choice_size = 0;
						choice_index++;
					}
					
					/* new state rows start after last one */
					for (; last_from < from; last_from++) {
						row_starts[last_from + 2] = row_starts[last_from + 1];
					}
					last_from = from;
					row_starts[from + 1]++;
					
					to=states.find(dst)->second;
					if(!isPS[from]) {
						non_zeros[nz_index] = rate;
						choice_size++;
					} else {
						non_zeros[nz_index] = 1;
						choice_size=1;
					}
					cols[nz_index] = to;
					nz_index++;
				}
			}
		}
		choice_starts[choice_index] = choice_starts[choice_index - 1] + choice_size;
		for (; from+1 < ma->n; from++) {
			row_starts[from+2] = row_starts[from+1];
		}
	}
	
}

/**
* Reads MA file @a filename.
*
* @param filename file to read MA from
* @return MA read from file
*/
SparseMatrix *read_IMC_SparseMatrix_file(const char *filename)
{
	bool error = false;
	unsigned long line_no;
	SparseMatrix *model = NULL;
	FILE *p = NULL;
	unsigned long num_states = 0;
	unsigned long num_ms_states = 0;
	map<string,unsigned long> states;
	map<unsigned long,string> states_nr;
	
	if (filename == NULL) {
		fprintf(stderr, "Called with filename == NULL\n");
		error = true;
	}
	
	/* first pass on file: create a hash table with state names and get number of states. */
	if (!error) {
		p = fopen(filename, "r");
		if (p == NULL) {
			fprintf(stderr, "Could not open file \"%s\"\n", filename);
			error = true;
		}
	}
	
	line_no = 1;
	
	if(!error)
		read_states(&line_no, &error, p, filename, &num_states, &states, &states_nr);
	
	if (p != NULL) {
		rewind(p);
	}
	
	line_no = 1;
	
	// check for deadlock states and add a selfloop
	if(!error)
		check_dedlocks(&line_no, &error, p, filename, &num_states, &states, &states_nr);
	
	if (p != NULL) {
		rewind(p);
	}
	
	if (!error) {
		model = SparseMatrix_new(num_states, states, states_nr); /* create MA model and reserve state memory */
	}
	
	if (!error) {
		model = SparseMatrix_new(num_states, states, states_nr); /* create MA model and reserve state memory */
	}

	
	/* second pass: count probabilistic states and store initial and goal states */
	line_no = 1;
	if(!error)
		init_states(&line_no, &error, p, filename, &num_ms_states, states, model); 
	
	if (!error) {
		rewind(p);
	}

	
	/* third pass on file: than reserve transition memory and store exit rates */
	line_no = 1;
	
	if (!error) {
		reserve_transition_memory(&line_no, &error, p, filename, model);
	}
	
	if (!error) {
		rewind(p);
	}
	
	
	/* fourth pass on file: save transitions */
	read_transitions(&line_no, &error, p, filename, model);
	
	if (p != NULL) {
		fclose(p);
	}
	
	if (error) {
		/* free the halfly-complete MDP structure if an error has occured */
		SparseMatrix_free(model);
		model = NULL;
	}else{
		//print_model(model);
		print_model_info(model);
	}
	
	//print_model(model);
	
	return model;
}
