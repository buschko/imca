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
*	Compute the unbounded reachability for an MA
*/

#include "unbounded.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <map>
#include <string>

#include "soplex.h"
#include "sccs.h"
#include "read_file.h"
#include "debug.h"

using namespace std;
using namespace soplex;

/**
* sets the objective function and bounds for goal states
*
* @param lp_model linear program
* @param ma the MA
* @param max identifier for maximum/minimum
*/
static void set_obj_function_unb(SoPlex& lp_model, SparseMatrix *ma, bool max, bool *locks) {
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
				lp_model.addCol(LPCol(1.0, dummycol, 1.0, 1.0));
			else
				lp_model.addCol(LPCol(1.0, dummycol, inf, 0));
		} else {
			lp_model.addCol(LPCol(0.0, dummycol, 0, 0));
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
static void set_constraints_unb(SoPlex& lp_model, SparseMatrix *ma, bool max, bool *locks) {
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
	
	int m=0; // greater equal 0
	if(!max)
		m=2; // less equal 0
	
	
	DSVector row(states);
	bool loop;
	
	for (state_nr = 0; state_nr < ma->n; state_nr++) {
		if(!goals[state_nr]  && !locks[state_nr]) {
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
					unsigned long r_start = rate_starts[state_nr];
					unsigned long r_end = rate_starts[state_nr + 1];
					for (unsigned long j = r_start; j < r_end; j++) {
						prob /= exit_rates[j];
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
					row.add(state_nr,-1);
				lp_model.addRow(LPRow(row,LPRow::Type(m), 0));
				row.~DSVector();
			}
		}
	}
	
	
}


/**
* Computes unbounded reachability for MA.
*
* @param ma file to read MA from
* @param max identifier for min or max prb
* @return probability for unbounded reachability
*/
Real compute_unbounded_reachability(SparseMatrix* ma, bool max) {
	dbg_printf("before soplex\n");
	SoPlex lp_model;
	dbg_printf("after soplex\n");
	
	bool *locks;
	
	if(max) {
		locks=compute_locks_strong(ma);
	} else {
		locks=compute_locks_weak(ma);
	}
	
	dbg_printf("LP computation start.\n");
	
	/* first step: build the lp model */
	set_obj_function_unb(lp_model,ma,max,locks);
	set_constraints_unb(lp_model,ma,max,locks);
	
	//lp_model.writeFile("file.lp", NULL, NULL, NULL);
	// TODO: (temorary BUGFIX: load model from file)
	//lp_model.readFile("file.lp");
	
	//cout << lp_model.rowVector(0) << endl;
	
	/* solve the LP */
	SPxSolver::Status stat;
	dbg_printf("solve\n");
	stat = lp_model.solve();
	
	//print_lp_info(lp_model);
	
	/* find the max or min prob. for an initial state */
	Real obj;
	if(max)
		obj=0;
	else
		obj=1;
			
	/* show if optimal solution */
	if( stat == SPxSolver::OPTIMAL ) {
		printf("LP solved to optimality.\n\n");
		//printf("Objective value is %lf.\n",lp_model.objValue());
		
		DVector probs(lp_model.nCols());
		lp_model.getPrimal(probs);
		
		bool *initials = ma->initials;
		unsigned long state_nr;
		for (state_nr = 0; state_nr < ma->n; state_nr++) {
			//cout << "s" << state_nr << " " << probs[state_nr] << endl;
			if(initials[state_nr]){
				Real tmp = probs[state_nr];
				dbg_printf("%li: prob %lf\n",state_nr, tmp);
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