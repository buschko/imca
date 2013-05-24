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
#include <vector>

#ifdef __SOPLEX__
#include "soplex.h"
#endif

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
* computes one step for Markovian states
*
* @param ma the MA
* @param v Markovian vector
* @param u result vector
*/
void compute_ub_markovian_vector(SparseMatrix* ma, vector<Real>& v, const vector<Real> &u, bool* locks){
	unsigned long *row_starts = (unsigned long *) ma->row_counts;
	unsigned long *choice_starts = (unsigned long *) ma->choice_counts;
	unsigned long *rate_starts = (unsigned long *) ma->rate_counts;
	unsigned long *cols = ma->cols;
	Real *non_zeros = ma->non_zeros;
	Real *exit_rates = ma->exit_rates;
	bool *goals = ma->goals;
	for (unsigned long state_nr = 0; state_nr < ma->n; state_nr++) {
		unsigned long state_start = row_starts[state_nr];
		unsigned long state_end = row_starts[state_nr + 1];
		// Look at Markovian states
		if(goals[state_nr]){
			v[state_nr]=1.0;
		}else if(locks[state_nr]){
			v[state_nr]=0.0;
		}else {
			if(!ma->isPS[state_nr])
			{
				for (unsigned long choice_nr = state_start; choice_nr < state_end; choice_nr++) {
					// Add up all outgoing rates of the distribution
					unsigned long i_start = choice_starts[choice_nr];
					unsigned long i_end = choice_starts[choice_nr + 1];
					Real exit_rate;
					unsigned long r_start = rate_starts[state_nr];
					unsigned long r_end = rate_starts[state_nr + 1];
					for (unsigned long j = r_start; j < r_end; j++) {
						 exit_rate = ma->exit_rates[j];
					}
					v[state_nr] = 0;
					for (unsigned long i = i_start; i < i_end; i++) {
						v[state_nr] += (non_zeros[i]/exit_rate) * u[cols[i]];
					}
				}
			}else {
					v[state_nr] = u[state_nr];
			}
		}
	}
}

/**
* computes one step for probabilistic states
*
* @param ma the MA
* @param v Markovian vector (we use it as temporary variable for u vector elements as well, so its value after calling of the function is not valid)
* @param u result vector
* @param max maximum/minimum
*/
void compute_ub_probabilistic_vector(SparseMatrix* ma, vector<Real>& v, vector<Real>& u, bool max, bool* locks){

	const Real precision = 1e-7;
	unsigned long statecount = ma-> n;
	unsigned long *row_starts = (unsigned long *) ma->row_counts;
	unsigned long *choice_starts = (unsigned long *) ma->choice_counts;
	unsigned long *cols = ma->cols;
	Real* non_zeros = ma->non_zeros;

	// Initialization
	for (unsigned long s_idx = 0; s_idx < statecount; s_idx++) {
		if ( ma -> isPS[s_idx] && !ma -> goals[s_idx] )
			v[s_idx] = 0;
		else
			u[s_idx] = v[s_idx];

	}
	// main loop
	bool done = false;
	// done will set to true inside the while loop when we reach fixed point
	while (! done) {
		done = true;


		// do one step of fixed point, it is applied on interactive states
		for (unsigned long s_idx = 0; s_idx < statecount; s_idx++) {
			unsigned long state_start = row_starts[s_idx];
			unsigned long state_end = row_starts[s_idx + 1];
			if ( ma -> isPS[s_idx] && !ma -> goals[s_idx] && !locks[s_idx]){  // do processing only if the state is interactive
				if (max)
					u[s_idx] = 0.0;
				else
					u[s_idx] = 1.0;

				// find the max/min prob. to reach Markovians and store it to tmp
				for (unsigned long choice_nr = state_start; choice_nr < state_end; choice_nr++) {
					Real tmp = 0;
					// Add up all outgoing rates of the distribution
					unsigned long i_start = choice_starts[choice_nr];
					unsigned long i_end = choice_starts[choice_nr + 1];
					for (unsigned long i = i_start; i < i_end; i++) {
						tmp += non_zeros[i] * v[cols[i]];
					}
					if( max ) {
						if(tmp > u[s_idx] )
							u[s_idx] = tmp;
					}
					else {
						if(tmp < u[s_idx] )
							u[s_idx] = tmp;
					}
				}
				if( fabs(u[s_idx] - v[s_idx]) >= precision )
					done = false;
					

			}
		}

 		bool *initials = ma->initials;
		// write back u to v, only for interactive elements
		for (unsigned long s_idx = 0; s_idx < statecount; s_idx++) {
			if ( ma -> isPS[s_idx] && !ma -> goals[s_idx] ){
				v[s_idx] = u[s_idx];
			}
		}

	}
	
}

/**
 * Computes the expected time reachability using a value iteration algorithm
 * Fixpoint is reached if afer 100 iteration steps, the distance is lower 
 * than epsilon=10^-6
 *
 * @return expected time vector
 */
Real unbounded_value_iteration(SparseMatrix* ma, bool max) {
	/* 
	 * Define value iteration likewise to time-bounded reachability
	 * with expected time property
	 */
	
	unsigned long num_states = ma->n;
	vector<Real> v(num_states,0); // Markovian vector
	vector<Real> u(num_states,0); // Probabilistic vector
	vector<Real> tmp(num_states,0); // tmp vector
	// initialize goal states
	bool *goals = ma->goals;
	for (unsigned long state_nr = 0; state_nr < num_states; state_nr++) {
		if(goals[state_nr]){
			v[state_nr]=1;
		}
	}
	
	bool *locks=(bool *)malloc(ma->n * sizeof(bool));
	/*
	if(max) {
		locks=compute_locks_weak(ma);
	} else {
		locks=compute_locks_strong(ma);
	}*/
	
	for(unsigned long i=0; i<ma->n; i++) {
		locks[i]=false;
	}
	
	cout << "start value iteration" << endl;
	
	bool done=false;
	
	while(!done){
		tmp=u;
		//done=true;
		// compute v for Markovian states: from b dwon to a, we make discrete model absorbing
		compute_ub_markovian_vector(ma,v,u,locks);
		// compute u for Probabilistic states
		compute_ub_probabilistic_vector(ma,v,u,max,locks);
		if(tmp==u)
			done=true;
	}

	// find prob. for initial state and return
	Real obj;
	if(max)
		obj=0;
	else
		obj=1;
	bool *initials = ma->initials;
	for (unsigned long state_nr = 0; state_nr < num_states; state_nr++) {
		//cout << (ma->states_nr.find(state_nr)->second).c_str() << ": " << u[state_nr] << endl;
		if(initials[state_nr]){
			if(max){
				if(obj<u[state_nr])
					obj=u[state_nr];
			}else{
				if(obj>u[state_nr])
					obj=u[state_nr];
			}
		}
	}

	 
	return obj;
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
	
	bool *locks=(bool *)malloc(ma->n * sizeof(bool));
	
	if(max) {
		locks=compute_locks_strong(ma);
	} else {
		locks=compute_locks_strong(ma);
	}
	/*
	for(unsigned long i=0; i<ma->n; i++) {
		locks[i]=false;
	}*/

	
	dbg_printf("LP computation start.\n");
	
	printf("LP computation start.\n");
	
	/* first step: build the lp model */
	set_obj_function_unb(lp_model,ma,max,locks);
	set_constraints_unb(lp_model,ma,max,locks);
	
	//lp_model.writeFile("file.lp", NULL, NULL, NULL);
	// TODO: (temorary BUGFIX: load model from file)
	//lp_model.readFile("file.lp");
	
	//cout << lp_model.rowVector(0) << endl;
	
	lp_model.setDelta(1e-4);
	
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