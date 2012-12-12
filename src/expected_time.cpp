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
* 	Compute the expected time for an MA
*/

#include "expected_time.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <map>
#include <string>

#include "soplex.h"
#include "sccs.h"
#include "read_file.h"

using namespace std;
using namespace soplex;

/**
* sets the objective function and bounds for goal states
*
* @param lp_model linear program
* @param ma the MA
* @param max identifier for maximum/minimum
*/
static void set_obj_function_ext(SoPlex& lp_model, SparseMatrix *ma, bool max, bool *locks) {
	unsigned long state_nr;
	bool *goals = ma->goals;
	DSVector dummycol(0);
	double inf = soplex::infinity;
	
	/* set objective function to max, resp. min */
	if(max)
		lp_model.changeSense(SPxLP::MINIMIZE);
	else
		lp_model.changeSense(SPxLP::MAXIMIZE);
	
	/* set objective and bounds resp. to goal states*/
	for (state_nr = 0; state_nr < ma->n; state_nr++) {
		if(!locks[state_nr]) {
			if(goals[state_nr])
				lp_model.addCol(LPCol(1.0, dummycol, 0.0, 0.0));
			else
				lp_model.addCol(LPCol(1.0, dummycol, inf, 0));
		} else {
			lp_model.addCol(LPCol(0.0, dummycol, inf, inf));
		}
	}
}

/**
* sets the constraints for the linear program
*
* @param lp_model linear program
* @param ma the MA
* @param max identifier for maximum/minimum
*/
static void set_constraints_ext(SoPlex& lp_model, SparseMatrix *ma, bool max, bool *locks) {
	unsigned long i;
	unsigned long state_nr;
	unsigned long choice_nr;
	unsigned long states = ma->n;
	bool *goals = ma->goals;
	//map<unsigned long,string> states_nr = ma->states_nr;
	unsigned long *row_starts = (unsigned long *) ma->row_counts;
	unsigned long *rate_starts = (unsigned long *) ma->rate_counts;
	unsigned long *choice_starts = (unsigned long *) ma->choice_counts;
	Real *non_zeros = ma->non_zeros;
	Real *exit_rates = ma->exit_rates;
	unsigned long *cols = ma->cols;
	Real prob;
	Real rate;
	int m=0; // greater equal 0
	if(!max)
		m=2; // less equal 0
	
	
	DSVector row(states);
	bool loop;
	
	for (state_nr = 0; state_nr < ma->n; state_nr++) {
		if(!goals[state_nr] && !locks[state_nr]) {
			unsigned long state_start = row_starts[state_nr];
			unsigned long state_end = row_starts[state_nr + 1];
			for (choice_nr = state_start; choice_nr < state_end; choice_nr++) {
				DSVector row(states);
				loop=false;
				/* Add up all outgoing rates of the distribution */
				unsigned long i_start = choice_starts[choice_nr];
				unsigned long i_end = choice_starts[choice_nr + 1];
				for (i = i_start; i < i_end; i++) {
					prob=non_zeros[i];
					rate=0;
					unsigned long r_start = rate_starts[state_nr];
					unsigned long r_end = rate_starts[state_nr + 1];
					for (unsigned long j = r_start; j < r_end; j++) {
						prob /= exit_rates[j];
						rate = -1/exit_rates[j];
					}
					//printf("%s - %lf -> %s\n",(states_nr.find(state_nr)->second).c_str(),prob,(states_nr.find(cols[i])->second).c_str());
					if(state_nr==cols[i]) {
						loop=true;
						row.add(state_nr,-1.0+prob);
					} else {
						row.add(cols[i],prob);
					}
				}
				if(!loop)
					row.add(state_nr,-1.0);
				lp_model.addRow(LPRow(row,LPRow::Type(m), rate));
				row.~DSVector();
			}
		}
	}
	
	
}


/**
* Computes expected time for MA.
*
* @param ma file to read MA from
* @param max identifier for min or max
* @return expected time
*/
Real compute_expected_time(SparseMatrix* ma, bool max) {
	SoPlex lp_model;
	bool *locks;
	if(max) {
		locks=compute_locks_weak(ma);
	} else {
		locks=compute_locks_strong(ma);
	}
	// printf("LP computation start.\n");
	
	for (unsigned int x = 0; x < ma->n; x++) {
		if(locks[x]){
			cout << x << endl;
		}
	}
	
	/* first step: build the lp model */
	set_obj_function_ext(lp_model,ma,max,locks);
	set_constraints_ext(lp_model,ma,max,locks);
	
	lp_model.writeFile("file.lp", NULL, NULL, NULL);
	// TODO: find out why LP model causes segfault in some cases (temorary BUGFIX: load model from file)
	//lp_model.readFile("file.lp");
	
	/* solve the LP */
	SPxSolver::Status stat;
	stat = lp_model.solve();
	
	//print_lp_info(lp_model);	
	
	/* find the max or min prob. for an initial state */
	Real obj;
	double inf = soplex::infinity;
	if(max)
		obj=0;
	else
		obj=inf;
		
	/* show if optimal solution */
	if( stat == SPxSolver::OPTIMAL ) {
		printf("LP solved to optimality.\n\n");
		//printf("Objective value is %lf.\n",lp_model.objValue());
		
		DVector probs(lp_model.nCols());
		lp_model.getPrimal(probs);
		
		bool *initials = ma->initials;
		unsigned long state_nr;
		for (state_nr = 0; state_nr < ma->n; state_nr++) {
			if(initials[state_nr]){
				Real tmp = probs[state_nr];
				if(max && tmp > obj)
					obj=tmp;
				else if(!max && tmp < obj)
					obj = tmp;
			}
		}
		
	} else if ( stat == SPxSolver::INFEASIBLE) {
		fprintf(stderr, "LP is infeasible.\n\n");
	}
	
	free(locks);

	return obj;
}