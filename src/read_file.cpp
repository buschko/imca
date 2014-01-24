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
*	Reads an compatible MA File
*/

#include "read_file.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <map>
#include <string>
#include <vector>

#include "debug.h"

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
		if (s[0] == '*') {
			++(*line_no);
		} else {
			sscanf(s, "%s%s", state, act);
			ret=tmp.insert(pair<string,unsigned long>(state,state_nr));
			if(ret.second == true){
				tmp_nr.insert(pair<unsigned long,string>(state_nr,state));
				state_nr++;
			} 
			++(*line_no);
		}
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
static void check_dedlocks(unsigned long *line_no, bool *error, FILE *p, const char *filename, unsigned long *num_states, map<string,unsigned long> *states, map<unsigned long,string> *states_nr, vector<unsigned long> *deadlocks) 
{
	char s[MAX_LINE_LENGTH];
	unsigned long state_nr = (*num_states);
	pair<map<string,unsigned long>::iterator,bool> ret;
	char state[MAX_LINE_LENGTH];
	char star[MAX_LINE_LENGTH];
	map<string ,unsigned long> tmp=*states;
	map<unsigned long, string> tmp_nr=*states_nr;
	bool deadlock=false;
	
	/* go to Transitions */
	if(fgets(s, MAX_LINE_LENGTH, p) == 0)
		{
			fprintf(stderr, "Reading line %ld of file \"%s\" failed.\n",*line_no, filename);
			*error = true;
		}
		++(*line_no);
		sscanf(s, "%s", state);
		if(strcmp(state,INITIALS) != 0 && !*error) {
			fprintf(stderr,"No declaration of initial states.\n");
			*error = true;
		} else {
			++(*line_no);
			if(fgets(s, MAX_LINE_LENGTH, p) == 0)
			{
				fprintf(stderr, "Reading line %ld of file \"%s\" failed.\n",*line_no, filename);
				*error = true;
			}
			sscanf(s, "%s", state);
			while(strcmp(state,GOALS) != 0 && !*error)
			{
				ret=tmp.insert(pair<string,unsigned long>(state,state_nr));
				if(ret.second == true){
					//cout << "Deadlock: " << state << "    State nr = " << state_nr << endl;
					deadlock=true;
					tmp_nr.insert(pair<unsigned long,string>(state_nr,state));
					(*deadlocks).push_back(state_nr);
					state_nr++;
				} 
				if(fgets(s, MAX_LINE_LENGTH, p) == 0)
				{
					fprintf(stderr, "Reading line %ld of file \"%s\" failed.\n",*line_no, filename);
					*error = true;
				}
				sscanf(s, "%s", state);
				++(*line_no);
			}
			if(fgets(s, MAX_LINE_LENGTH, p) == 0)
			{
				fprintf(stderr, "Reading line %ld of file \"%s\" failed.\n",*line_no, filename);
				*error = true;
			}
			sscanf(s, "%s", state);
			while(strcmp(state,TRANSITIONS) != 0 && !*error) {
				ret=tmp.insert(pair<string,unsigned long>(state,state_nr));
				if(ret.second == true){
					//cout << "Deadlock: " << state << "    State nr = " << state_nr << endl;
					deadlock=true;
					tmp_nr.insert(pair<unsigned long,string>(state_nr,state));
					(*deadlocks).push_back(state_nr);
					state_nr++;
				} 
				if(fgets(s, MAX_LINE_LENGTH, p) == 0)
				{
					fprintf(stderr, "Reading line %ld of file \"%s\" failed.\n",*line_no, filename);
					*error = true;
				}
				sscanf(s, "%s", state);
				++(*line_no);
			}
		}
	
	// cout << "hey" << endl;
	
	/* TODO: also check lines */
	while (fgets(s, MAX_LINE_LENGTH, p) != NULL && !*error) {
		if (s[0] == '*') {
			sscanf(s, "%s%s", star, state);
			ret=tmp.insert(pair<string,unsigned long>(state,state_nr));
			if(ret.second == true){
				cout << "Deadlock: " << state << "    State nr = " << state_nr << endl;
				//cout << "";
				deadlock=true;
				tmp_nr.insert(pair<unsigned long,string>(state_nr,state));
				(*deadlocks).push_back(state_nr);
				state_nr++;
			}
			//cout << state << endl;
			++(*line_no);
		} else {
			++(*line_no);
		}
	}
	
	if(deadlock) {
		(*num_states) = state_nr;
		(*states)=tmp;
		(*states_nr)=tmp_nr;
	}
}

static void init_states(unsigned long *line_no, bool *error, FILE *p, const char *filename, unsigned long *num_ms_states, map<string,unsigned long> states, SparseMatrix *ma, vector<unsigned long> deadlocks) {
	char src[MAX_LINE_LENGTH];
	char act[MAX_LINE_LENGTH];
	char s[MAX_LINE_LENGTH];
	unsigned long from;
	unsigned long last_from = 0;
	unsigned long ms_states=0;
	bool *isPS;
	bool first=true;
	
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
		if (s[0] != '*') {
			if (sscanf(s, "%s%s", src, act) != 2) { /* check if both values are returned */
				*error = true;
			} else {
				from=states.find(src)->second;
				/* test consistency of ma file */
				if (from < last_from) {
					fprintf(stderr, "Line %ld: State transitions must be given in continuous order for one state.\n", *line_no);
					*error = true;
				} else if(from==last_from && !first){
					if(strcmp(act,MARKOV_ACTION) == 0 && !isPS[from]) {
						*error = true;
						fprintf(stderr, "Line %ld: Markovian state transitions of one state should not be divided.\n", (*line_no));
					} else {
						isPS[from]=true;
					}				
				} else {
					if(strcmp(act,MARKOV_ACTION) != 0)
						isPS[from]=true;
					first=false;
					last_from = from;
				}
			}
		}
		++(*line_no);
	}
	
	vector<unsigned long>::iterator it;
	for(it=deadlocks.begin(); it<deadlocks.end(); it++) {
		isPS[(*it)]=true;
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
static void reserve_transition_memory(unsigned long *line_no, bool *error, FILE *p, const char *filename, SparseMatrix *model, vector<unsigned long> deadlocks, bool mrm) {
	char s[MAX_LINE_LENGTH];
	char src[MAX_LINE_LENGTH];
	char dst[MAX_LINE_LENGTH];
	char act[MAX_LINE_LENGTH];
	char star[MAX_LINE_LENGTH];
	Real exit_rate=0;
	Real max_exit_rate=0;
	Real prob;
    Real reward;
	unsigned long exit_index = 0;
	unsigned long num_choice = 0;
	unsigned long num_non_zeros = 0;
	bool is_ms = false;
	Real *exit_rates;
	unsigned long *rate_starts;
	unsigned long from;
	unsigned long last_from = 0;
	unsigned long to;
	unsigned long last_to = -1;
	map<string ,unsigned long> states;
	bool *isPS;

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
		if (s[0] != '*') {
			if (sscanf(s, "%s%s", src, act) != 2) { /* check if both values are returned */
				*error = true;
			} else {
				from=states.find(src)->second;
				last_to = -1;
				/* test consistency of ma file */
				if (from < last_from) {
					fprintf(stderr, "Line %ld: State transitions must be given in continous order for one state.\n", *line_no);
					*error = true;
				} else {	
					/* if Markovian state, store exit rate */
					if(is_ms) {
						if(max_exit_rate<exit_rate)
							max_exit_rate=exit_rate;
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
					
						
					is_ms=false;
					exit_rate=0;
					if(strcmp(act,MARKOV_ACTION) == 0 && !isPS[from])
					{
						is_ms=true;
					}	
					
					/* probabilistic transitions are choosen before markovian transitions */
					if((isPS[from] && !is_ms) || (!isPS[from] && is_ms)){
						num_choice++;
					}
					
					last_from = from;
				}	
			}
			//++(*line_no);
		} else {
			/* probabilistic transitions are choosen before markovian transitions */
			Real rate;
			sscanf(s, "%s%s%lf", star, dst, &rate);
			to = states.find(dst)->second;
			if(((isPS[from] && !is_ms) || (!isPS[from] && is_ms)) && to != last_to) {
				num_non_zeros++;
				if(is_ms){
					//Real rate;
					//char star[MAX_LINE_LENGTH];
					//sscanf(s, "%s%s%lf", star, dst, &rate);
					exit_rate += rate;
				}
			}
			last_to = to;
		}
		++(*line_no);
	}
	
	//printf("nz %ld  nc %ld\n",num_non_zeros,num_choice);
	
	if(deadlocks.size() == 0) {
		/* probabilistic transitions are choosen before markovian transitions */
		if((is_ms)) {
			exit_rates[exit_index] = exit_rate;
			for (; from+0 < model->n; from++) {
				rate_starts[from+1] = rate_starts[from+0];
			}
			rate_starts[from + 0]++;
		} else {
			//cout << from << " " << last_from << endl;
			for (; last_from < from; last_from++) {
				rate_starts[last_from + 1] = rate_starts[last_from + 0];
			}
		}
	}
	else{
		vector<unsigned long>::iterator it;
		for(it=deadlocks.begin(); it<deadlocks.end(); it++) {
			from=(*it);
			if(is_ms) {
				exit_rates[exit_index] = exit_rate;
				exit_index++;
				if(max_exit_rate<exit_rate)
					max_exit_rate=exit_rate;
				for (; last_from < from; last_from++) {
					rate_starts[last_from + 1] = rate_starts[last_from + 0];
				}
				rate_starts[from + 0]++;
				is_ms=false;
			} else {
				for (; last_from < from; last_from++) {
					rate_starts[last_from + 1] = rate_starts[last_from + 0];
				}
			}
			num_choice++;
			num_non_zeros++;
			last_from = from;
		}
		for (; last_from < from; last_from++) {
			rate_starts[last_from + 1] = rate_starts[last_from + 0];
		}
	}
	
	//printf("nz %ld  nc %ld\n",num_non_zeros,num_choice);
	
	/* now allocate the memory needed */
	if (!*error) {
		unsigned long * choice_starts = (unsigned long *) calloc((size_t) (num_choice + 1), sizeof(unsigned long));
		Real * non_zeros = (Real *) malloc(num_non_zeros * sizeof(Real));
		unsigned long * cols = (unsigned long *) malloc(num_non_zeros * sizeof(unsigned long));
		model->choice_counts = (unsigned char *) choice_starts;
		model->non_zeros = non_zeros;
		model->cols = cols;
		model->choices_n = num_choice;
		model->non_zero_n = num_non_zeros;
		model->max_exit_rate=max_exit_rate;
		if(mrm){
            Real * rewards = (Real *) malloc(num_choice * sizeof(Real));
			model->rewards=rewards;
        }
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
static void read_transitions(unsigned long *line_no, bool *error, FILE *p, const char *filename, SparseMatrix *ma, vector<unsigned long> deadlocks,bool mrm) {
	unsigned long choice_index = 0;
	unsigned long choice_size = 0;
	unsigned long nz_index = 0;
	unsigned long reward_index = 0;
	char src[MAX_LINE_LENGTH];
	unsigned long from;
	unsigned long last_from = 0;
	char s[MAX_LINE_LENGTH];
	char dst[MAX_LINE_LENGTH];
	unsigned long to;
	char act[MAX_LINE_LENGTH];
	bool bad=false;
	unsigned long last_to = -1;
	//Real r=0;
	
	Real tmp = 0;
	unsigned int old_index;

	if (!*error) {
		Real *non_zeros = ma->non_zeros;
		Real *rewards = ma->rewards;
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
			if (s[0] == '*') {
				if(!bad)
				{
					Real rate;
					char star[MAX_LINE_LENGTH];
					sscanf(s, "%s%s%lf", star, dst, &rate);
					to=states.find(dst)->second;
					if(to != last_to) {
						non_zeros[nz_index] = rate;
						cols[nz_index] = to;
						old_index = nz_index;
						nz_index++;
						choice_size++;
					}else if(to == last_to) {
						non_zeros[old_index] += rate;
					}
					last_to = to;
				}
			} else {
				last_to = -1;
				/* just read a line "state act" */
				// if mrm we also read in the reward
				if(mrm){
					Real reward=0;
					sscanf(s, "%s%s%lf", src, act, &reward);
					rewards[reward_index] = reward;
                    //r += reward;
					reward_index++;
				}else {
					sscanf(s, "%s%s", src, act);
				}
				from=states.find(src)->second;
				/* probabilistic transitions are choosen before markovian transitions */
				if(strcmp(act,MARKOV_ACTION) == 0 && isPS[from])
					bad=true;
				else
					bad=false;
				if(!bad) {
					if (from == last_from) {
						row_starts[from + 1]++;
					} else {
						/* new state rows start after last one */
						for (; last_from < from; last_from++) {
							row_starts[last_from + 2] = row_starts[last_from + 1];
						}
						row_starts[from + 1]++;
						last_from = from;
					}
					if (choice_index > 0) {
						choice_starts[choice_index] =
						choice_starts[choice_index - 1] + choice_size;
					}
					choice_size = 0;
					choice_index++;
				}
			}
		}
		if(deadlocks.size() == 0) {
			choice_starts[choice_index] = choice_starts[choice_index - 1] + choice_size;
			for (; from+1 < ma->n; from++) {
				row_starts[from+2] = row_starts[from+1];
			}
		}else {
			if (choice_index > 0) {
				choice_starts[choice_index] =
				choice_starts[choice_index - 1] + choice_size;
			}
			choice_index++;
			vector<unsigned long>::iterator it;
			for(it=deadlocks.begin(); it<deadlocks.end(); it++) {
				from=(*it);
				to=(*it);
				non_zeros[nz_index] = 1;
				cols[nz_index] = to;
				nz_index++;
				choice_size=1;
				for (; last_from < from; last_from++) {
					row_starts[last_from + 2] = row_starts[last_from + 1];
				}
				row_starts[from + 1]++;
				last_from = from;
				if (choice_index > 0) {
					choice_starts[choice_index] =
					choice_starts[choice_index - 1] + choice_size;
				}
				choice_index++;
			}
			
			//choice_starts[choice_index] = choice_starts[choice_index - 1] + choice_size;
			for (; from+1 < ma->n; from++) {
				row_starts[from+2] = row_starts[from+1];
			}
		}
	}
    if(mrm){
        //ma->reward=r;
    }
	
}

void print_model(SparseMatrix *ma, bool mrm)
{
	unsigned long i;
	unsigned long tau=1;
	unsigned long state_nr;
	unsigned long choice_nr;
	map<unsigned long,string> states_nr = ma->states_nr;
	unsigned long *row_starts = (unsigned long *) ma->row_counts;
	unsigned long *rate_starts = (unsigned long *) ma->rate_counts;
	unsigned long *choice_starts = (unsigned long *) ma->choice_counts;
	Real *non_zeros = ma->non_zeros;
	Real *rewards = ma->rewards;
	Real *exit_rates = ma->exit_rates;
	unsigned long *cols = ma->cols;
	Real prob;
	bool *initials = ma->initials;
	bool *goals = ma->goals;
	
	
	for (state_nr = 0; state_nr < ma->n; state_nr++) {
		unsigned long state_start = row_starts[state_nr];
		unsigned long state_end = row_starts[state_nr + 1];
		dbg_printf("row_starts: %li row_ends: %li\n",state_start,state_end);
		for (choice_nr = state_start; choice_nr < state_end; choice_nr++) {
			/* Add up all outgoing rates of the distribution */
			unsigned long i_start = choice_starts[choice_nr];
			unsigned long i_end = choice_starts[choice_nr + 1];
			dbg_printf("choice_starts: %li choice_ends: %li\n",i_start,i_end);
			printf("choice %d\n",tau);
			if(mrm){
				printf("reward: %lg\n",rewards[choice_nr]);
			}
			tau++;
			for (i = i_start; i < i_end; i++) {
				prob=non_zeros[i];
				unsigned long r_start = rate_starts[state_nr];
				unsigned long r_end = rate_starts[state_nr + 1];
				dbg_printf("rate_starts: %li rate_ends: %li\n",r_start,r_end);
				for (unsigned long j = r_start; j < r_end; j++) {
					prob /= exit_rates[j];
				}
				printf("%s - %lg -> %s\n",(states_nr.find(state_nr)->second).c_str(),prob,(states_nr.find(cols[i])->second).c_str());
			}
		}
		tau=1;
		unsigned long i_start = rate_starts[state_nr];
		unsigned long i_end = rate_starts[state_nr + 1];
		for (i = i_start; i < i_end; i++) {
			printf("ExitRate(%s): %lg\n",(states_nr.find(state_nr)->second).c_str(),exit_rates[i]);
		}
		if(initials[state_nr])
			printf("initial state %s\n",(states_nr.find(state_nr)->second).c_str());
		if(goals[state_nr])
			printf("goal state %s\n",(states_nr.find(state_nr)->second).c_str());
	}
}

void print_model_info(SparseMatrix *ma)
{	
	unsigned long state_nr;
	bool *initials = ma->initials;
	bool *goals = ma->goals;
	unsigned long n_init=0;
	unsigned long n_goal=0;
	unsigned long n_trans=0;
	
	unsigned long *row_starts = (unsigned long *) ma->row_counts;
	//unsigned long *rate_starts = (unsigned long *) ma->rate_counts;
	unsigned long *choice_starts = (unsigned long *) ma->choice_counts;
	
	unsigned long nd_states=0;
	for (state_nr = 0; state_nr < ma->n; state_nr++) {
		if(initials[state_nr])
			n_init++;
		if(goals[state_nr])
			n_goal++;
		unsigned long state_start = row_starts[state_nr];
		unsigned long state_end = row_starts[state_nr + 1];
		for (unsigned long choice_nr = state_start; choice_nr < state_end; choice_nr++) {
			unsigned long i_start = choice_starts[choice_nr];
			unsigned long i_end = choice_starts[choice_nr + 1];
			for (unsigned long i = i_start; i < i_end; i++) {
				n_trans++;
			}
		}
		if(((state_end)-state_start)>1)
			nd_states++;
	}
	
	printf("#States: %ld    #MS: %ld    #PS: %ld    #ND: %ld\n",ma->n, ma->ms_n, (ma->n - ma->ms_n),nd_states);
	printf("#Initials: %ld    #Goals: %ld    #Transitions: %ld\n", n_init, n_goal, n_trans);
}

#ifdef __SOPLEX__
/**
* print out onformation about SoPlex LP solve.
*
* @param lp_model the LP
*/
void print_lp_info(SoPlex lp_model) {
	printf("\n");
	printf("SoPlex parameters:\n");
	printf("Delta          = %g\n",lp_model.delta());
	printf("Epsilon Zero   = %g\n", Param::epsilon());  
	printf("Epsilon Factor = %g\n", Param::epsilonFactorization());
	printf("Epsilon Update = %g\n", Param::epsilonUpdate());
	printf("\n");
	printf("algorithm      = %s\n", (lp_model.type() == SPxSolver::ENTER ? "Entering" : "Leaving"));
	printf("representation = %s\n", (lp_model.rep() == SPxSolver::ROW ? "Row" : "Column"));
	printf("piercing       = %s\n", (lp_model.pricing() == SPxSolver::FULL ? "Full" : "Partial"));
	printf("\n");
	printf("SoPlex LP solve information:\n");
	printf("Factorizations : %d\n",lp_model.getFactorCount());
	printf("  Time spent   : %f\n",lp_model.getFactorTime());
	printf("Solves         : %d\n",lp_model.getSolveCount());
	printf("  Time spent   : %f\n",lp_model.getSolveTime());
	printf("Iterations     : %d\n",lp_model.iteration());
}
#endif

/**
* Reads MA file @a filename.
*
* @param filename file to read MA from
* @return MA read from file
*/
SparseMatrix *read_MA_SparseMatrix_file(const char *filename, bool mrm)
{
	bool error = false;
	unsigned long line_no;
	SparseMatrix *model = NULL;
	FILE *p = NULL;
	unsigned long num_states = 0;
	unsigned long num_ms_states = 0;
	map<string,unsigned long> states;
	map<unsigned long,string> states_nr;
    states.clear();
    states_nr.clear();
	
	if (filename == NULL) {
		fprintf(stderr, "Called with filename == NULL\n");
		error = true;
	}
	//cout << "first pass" << endl;
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
	//cout << "deadlock pass" << endl;
	vector<unsigned long> deadlocks;
	deadlocks.clear();
	// check for deadlock states and add a selfloop
	if(!error)
		check_dedlocks(&line_no, &error, p, filename, &num_states, &states, &states_nr, &deadlocks);
	
	if (p != NULL) {
		rewind(p);
	}
	
	line_no = 1;
	
    //cout << "reserve memory" << endl;
	if (!error) {
		model = SparseMatrix_new(num_states, states, states_nr); /* create MA model and reserve state memory */
	}
	//cout << "second pass" << endl;
	/* second pass: count probabilistic states and store initial and goal states */
	line_no = 1;
	if(!error)
		init_states(&line_no, &error, p, filename, &num_ms_states, states, model, deadlocks); 
	
	if (!error) {
		rewind(p);
	}
	//cout << "third pass" << endl;
	/* third pass on file: than reserve transition memory and store exit rates */
	line_no = 1;
	
	if (!error) {
		reserve_transition_memory(&line_no, &error, p, filename, model, deadlocks, mrm);
	}
	
	if (!error) {
		rewind(p);
	}
	//cout << "fourth pass" << endl;
	/* fourth pass on file: save transitions */
	read_transitions(&line_no, &error, p, filename, model, deadlocks,mrm);
    
	if (p != NULL) {
		fclose(p);
	}
    
	if (error) {
		/* free the halfly-complete MA structure if an error has occured */
		SparseMatrix_free(model);
		model = NULL;
	}else{
		//print_model(model,mrm);
		print_model_info(model);
	}
	
	//print_model(model,mrm);
	
	return model;
}