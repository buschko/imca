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
*	Compute the unbounded reachability for an MA
*/


#include "bounded.h"
#include "debug.h"
#include "sccs.h"
#include <math.h>
#include <vector>
#include <iostream>

using std::vector;

/**
* sets the error bound for given epsilon and tb
*
* @param ma the MA
* @param epsilon epsilon bound
* @param tb time bound
*/
Real compute_error_bound(SparseMatrix* ma, Real epsilon, Real tb) {
	Real tau;
	// TODO: add lambda to read in function
	Real lambda=ma->max_exit_rate;
	tau = (2*epsilon)/(tb*lambda*lambda);
	return tau;
}

/**
* discretizes the MA
*
* @param ma the MA
* @param tau discretization factor
* @return discretized MA
*/
SparseMatrix* discretize_model(SparseMatrix* ma, Real tau) {
	SparseMatrix* discrete_ma;
	dbg_printf("memory alloc.\n");
	discrete_ma=SparseMatrixDiscrete_new(ma);
	
	// transitions for MA
	unsigned long *row_starts = (unsigned long *) ma->row_counts;
	unsigned long *rate_starts = (unsigned long *) ma->rate_counts;
	unsigned long *choice_starts = (unsigned long *) ma->choice_counts;
	
	// transition variables for discrete_ma
	//unsigned long *d_choice_starts = (unsigned long *) discrete_ma->choice_counts;
	//unsigned long *d_row_starts = (unsigned long *) discrete_ma->row_counts;
	Real *non_zeros = discrete_ma->non_zeros;
	unsigned long *cols = discrete_ma->cols;
	unsigned long nz_index = 0;
	//unsigned long choice_index = 0;
	//unsigned long choice_size = 0;
	
	dbg_printf("discretization\n");
	// add for each Markovian state a loop transition ==> new memory allocation for non_zero, choices and successors
	for (unsigned long state_nr = 0; state_nr < ma->n; state_nr++) {
		unsigned long state_start = row_starts[state_nr];
		unsigned long state_end = row_starts[state_nr + 1];
		unsigned long r_start = rate_starts[state_nr];
		//unsigned long r_end = rate_starts[state_nr];
		dbg_printf("check if Markovian\n");
		// Look at Markovian states
		if(!ma->isPS[state_nr])
		{
			for (unsigned long choice_nr = state_start; choice_nr < state_end; choice_nr++) {
				// Add up all outgoing rates of the distribution
				unsigned long i_start = choice_starts[choice_nr];
				unsigned long i_end = choice_starts[choice_nr + 1];
				Real exit_rate = ma->exit_rates[r_start];
				//Real tmp=Real(1)-exp(-(exit_rate*tau));
				// precision of Real could be to small (exp_estau)
				Real estau = -(exit_rate * tau);	
				Real exp_estau = exp(estau);
				Real exp_estau_com = Real(1)-exp_estau;
				bool loop=false;
				for (unsigned long i = i_start; i < i_end; i++) {
					Real prob = ma->non_zeros[i]/exit_rate;
					non_zeros[nz_index] = exp_estau_com*prob;
					cols[nz_index] = ma->cols[i];
					if(state_nr==ma->cols[i]){
						loop=true;
						non_zeros[nz_index] += exp_estau;
					}
					nz_index++;
				}
				dbg_printf("add selfloop\n");
				// add selfloop 
				if(!loop) {
					non_zeros[nz_index] = exp_estau;
					cols[nz_index] = state_nr;
					nz_index++;
				}
			}
		}else {
			for (unsigned long choice_nr = state_start; choice_nr < state_end; choice_nr++) {
				// Add up all outgoing rates of the distribution
				unsigned long i_start = choice_starts[choice_nr];
				unsigned long i_end = choice_starts[choice_nr + 1];
				for (unsigned long i = i_start; i < i_end; i++) {
					non_zeros[nz_index] = ma->non_zeros[i];
					cols[nz_index] = ma->cols[i];
					nz_index++;
				}
			}
		}
	}
	return discrete_ma;
}

/**
* computes one step for Markovian states
*
* @param ma the MA
* @param v Markovian vector
* @param u result vector
* @param is_MA_made_absorbing: if true then all goal states are made absorbing, otherwise not
*/
void compute_markovian_vector(SparseMatrix* ma, vector<Real>& v, const vector<Real> &u, bool is_MA_made_absorbing){
	unsigned long *row_starts = (unsigned long *) ma->row_counts;
	unsigned long *choice_starts = (unsigned long *) ma->choice_counts;
	unsigned long *cols = ma->cols;
	Real* non_zeros = ma->non_zeros;
	bool *goals = ma->goals;
	for (unsigned long state_nr = 0; state_nr < ma->n; state_nr++) {
		unsigned long state_start = row_starts[state_nr];
		unsigned long state_end = row_starts[state_nr + 1];
		// Look at Markovian states
		if(goals[state_nr] && is_MA_made_absorbing){
			v[state_nr]=1;
		}else{
			if(!ma->isPS[state_nr])
			{
				for (unsigned long choice_nr = state_start; choice_nr < state_end; choice_nr++) {
					// Add up all outgoing rates of the distribution
					unsigned long i_start = choice_starts[choice_nr];
					unsigned long i_end = choice_starts[choice_nr + 1];
					v[state_nr]=0;
					for (unsigned long i = i_start; i < i_end; i++) {
						v[state_nr] += non_zeros[i] * u[cols[i]];
					}
				}
			}else {
					v[state_nr] = u[state_nr];
			}
		}
	}
}

/**
* computes one step for Interactive states
*
* @param ma the MA
* @param v Markovian vector
* @param u result vector
* @param max maximum/minimum
* @param locks Lock set for states never reach a goal state
* @param reach reachability for interactive states
*/
void compute_interactive_vector(SparseMatrix* ma, vector<Real> v, vector<Real>& u, bool max, bool* locks, vector< vector<unsigned long> > reach) {
	bool *goals = ma->goals;
	vector<unsigned long>::const_iterator r;
	for (unsigned long state_nr = 0; state_nr < ma->n; state_nr++) {
		if(locks[state_nr]){
			u[state_nr] = 0;
		}else if (goals[state_nr]){
			u[state_nr] = 1;
		}else{
			if(max) {
				u[state_nr]=0;
				for(r=reach[state_nr].begin(); r<reach[state_nr].end(); r++) {
					if(u[state_nr] < v[(*r)]){
						u[state_nr] = v[(*r)];
					}
				}
			}else {
				u[state_nr]=1;
				for(r=reach[state_nr].begin(); r<reach[state_nr].end(); r++) {
					if(u[state_nr] > v[(*r)]){
						u[state_nr] = v[(*r)];
					}
				}
			}
		}
	}
}

/**
* computes reachability for interactive states
*
* @param ma the MA
* @param state_nr the actual state
* @param visited states visited before
*/
void interactiveReachability(SparseMatrix* ma, const unsigned long state_nr, vector<unsigned long>& visited) {
	unsigned long dst;
	visited[state_nr]=1;
	
	unsigned long *row_starts = (unsigned long *) ma->row_counts;
	unsigned long *choice_starts = (unsigned long *) ma->choice_counts;
	unsigned long *cols = ma->cols;
	bool *goals = ma->goals;
	unsigned long state_start = row_starts[state_nr];
	unsigned long state_end = row_starts[state_nr + 1];
	
	if(!goals[state_nr]){
		if(ma->isPS[state_nr])
		{
			for (unsigned long choice_nr = state_start; choice_nr < state_end; choice_nr++) {
				// Add up all outgoing rates of the distribution
				unsigned long i_start = choice_starts[choice_nr];
				unsigned long i_end = choice_starts[choice_nr + 1];
				for (unsigned long i = i_start; i < i_end; i++) {
					dst=cols[i];
					if(visited[dst]==0) {
						interactiveReachability(ma,dst,visited);
					}
				}
			}
		}
	}
}

/**
* computes reachability for interactive states
*
* @param ma the MA
* @param state_nr the actual state
*/
vector<unsigned long> interactiveReachability(SparseMatrix* ma, const unsigned long state) {
	vector<unsigned long> visited(ma->n);
	vector<unsigned long> goals;
	unsigned long x;
	for(x=0; x < ma->n; x++) {
		visited[x]=0;
	}
	interactiveReachability(ma,state, visited);
	for(x=0; x < ma->n; x++) {
		if(visited[x]==1) {
		goals.push_back(x);
		}
	}
	return goals;
}

/**
* computes reachability for interactive states
*
* @param ma the MA
*/
vector< vector<unsigned long> > interactiveReachability(SparseMatrix* ma) {
	vector< vector<unsigned long> > r;
	vector<unsigned long> reachSucc;
	unsigned long s_idx;
	unsigned long statecount = ma->n;
	
	for (s_idx = 0; s_idx < statecount; s_idx++) {
		reachSucc.clear();
		reachSucc=interactiveReachability(ma,s_idx);
		r.push_back(reachSucc);
	}
	
	return r;
}

/**
* computes one step for probabilistic states
*
* @param ma the MA
* @param v Markovian vector (we use it as temporary variable for u vector elements as well, so its value after calling of the function is not valid)
* @param u result vector
* @param max maximum/minimum
* @param is_MA_made_absorbing: if true then all goal states are made absorbing, otherwise not
*/
void compute_probabilistic_vector(SparseMatrix* ma, vector<Real>& v, vector<Real>& u, bool max, bool is_MA_made_absorbing){

	const Real precision = 1e-7;
	unsigned long statecount = ma-> n;
	unsigned long *row_starts = (unsigned long *) ma->row_counts;
	unsigned long *choice_starts = (unsigned long *) ma->choice_counts;
	unsigned long *cols = ma->cols;
	Real* non_zeros = ma->non_zeros;

	// Initialization
	for (unsigned long s_idx = 0; s_idx < statecount; s_idx++) {
		if ( ma -> isPS[s_idx] && (!is_MA_made_absorbing || (is_MA_made_absorbing && !ma -> goals[s_idx])) )
			v[s_idx] = 0.0;
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
			if ( ma -> isPS[s_idx] && (!is_MA_made_absorbing || (is_MA_made_absorbing && !ma -> goals[s_idx])) ){  // do processing only if the state is interactive
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

		// write back u to v, only for interactive elements
		for (unsigned long s_idx = 0; s_idx < statecount; s_idx++) {
			if ( ma -> isPS[s_idx] && (!is_MA_made_absorbing || (is_MA_made_absorbing && !ma -> goals[s_idx])) ){
				v[s_idx] = u[s_idx];
			}
		}


	}



}

/**
* computes time-bounded reachability
*
* @param ma the MA
* @param max maximum/minimum
* @param epsilon the given error
* @param tb the given time bound
* @param is_imc indicates if MA is an IMC
*/
Real compute_time_bounded_reachability(SparseMatrix* ma, bool max, Real epsilon, Real ta, Real tb, bool is_imc, Real interval,Real interval_start) {

	// Check whether the given time interval is zero
	if( ta > tb ) {
		printf("WARNING: The given interval is empty (upper bound < lower bound.) The reachability probability is 0.\n");
		return 0.0;
	}

	unsigned long num_states = ma->n;
	vector<Real> v(num_states,0); // Markovian vector
	vector<Real> u(num_states,0); // Probabilistic vector
	// initialize goal states
	bool *goals = ma->goals;
	for (unsigned long state_nr = 0; state_nr < num_states; state_nr++) {
		if(goals[state_nr]){
			v[state_nr]=1;
		}
	}
	// in case MA is an IMC: precomputation of paths for interactive states
	vector< vector<unsigned long> > reach;
	bool *locks;
	if(is_imc) {
		std::cout << "precomputation" << std::endl;
		if(max){
			locks=compute_locks_strong(ma);
			//{s' in S | s ~>i* s'}
			reach = interactiveReachability(ma);
		}else{
			locks=compute_locks_weak(ma);
			//{s' in PS U G | s ~>min* s'}
			reach = interactiveReachability(ma); 
		}
	}
	
	std::cout << "start value iteration" << std::endl;
	if( ta > 0 ) {
		//TODO one unique function that tell you what should be tau1 and what tau2
		// compute the upper bound of discretization step
		Real tau = compute_error_bound(ma, epsilon,tb);

		unsigned long steps; // stores the number of steps should be taken
		Real current_tau;    // stores the current value of tau which is calculated to make sure the number of steps is a integral value
		SparseMatrix* discrete_ma; // discretized model induced from the original ma

		steps = (unsigned long) ceil((tb - ta) / tau); // calculating the number of steps for interval tb down to ta
		current_tau = (tb - ta) / (Real)steps;               // recalculate tau based on the number of steps


		// discretize model with respect to the current value of tau 'current_tau'
		dbg_printf("discretize model for interval [%g,%g] ... \n", ta, tb);
		discrete_ma = discretize_model(ma,current_tau);
		dbg_printf("model discretized\n");
		//print_model(discrete_ma);
		
		//unsigned long interval_step = round(interval/current_tau);
		//unsigned long counter=0;
		//std::cout << "interval step: " << interval_step <<std::endl;


		// value iteration
		std::cout << "iterations: " << (unsigned long) ceil((tb - ta) / tau) + (unsigned long) ceil(ta / tau) << std::endl;
		printf("step duration for interval [%g,%g]: %g\n", ta, tb, current_tau);
		std::cout << "iterations: " << (unsigned long) ceil((tb - ta) / tau) << std::endl;
		
		for(unsigned long i=0; i < steps; i++){
			// compute v for Markovian states: from b dwon to a, we make discrete model absorbing
			compute_markovian_vector(discrete_ma,v,u, true);
			// compute u for Probabilistic states
			if(is_imc){
				// if MA is in fact an IMC we can simplify the computation (slower than the PA computation TODO: find bug)
				compute_interactive_vector(discrete_ma,v,u,max,locks,reach);
			}else {
				compute_probabilistic_vector(discrete_ma,v,u,max, true);
			}
			
			/*
			if(counter==interval_step) {
			
				Real prob;
				if(max)
					prob=0;
				else
					prob=1;
				bool *initials = ma->initials;
				for (unsigned long state_nr = 0; state_nr < num_states; state_nr++) {
					//std::cout << (ma->states_nr.find(state_nr)->second).c_str() << ": " << u[state_nr] << std::endl;
					if(initials[state_nr]){
					if(max){
						if(prob<u[state_nr])
							prob=u[state_nr];
						}else{
						if(prob>u[state_nr])
							prob=u[state_nr];
						}
					}
				}
				
				printf("tb=%.5g Maximal time-bounded reachability probability: %.10g\n", ta+i*tau,prob);
			
				counter=0;
			} else {
				counter++;
			}
			*/
			
		}
		// free memory
		SparseMatrix_free(discrete_ma);

		steps = (unsigned long) ceil( (ta + current_tau) / tau); // calculating the number of steps for interval [0,a]
		current_tau = (ta + current_tau) / (Real)steps;                // recalculate tau based on the number of steps

		// discretize model with respect to the current value of tau 'current_tau'
		dbg_printf("discretize model for interval [0,%g] ... \n", ta);
		discrete_ma = discretize_model(ma,current_tau);
		dbg_printf("model discretized\n");
		//print_model(discrete_ma);

		printf("step duration for interval [0,%g]: %g\n", ta, current_tau);
		std::cout << "iterations: " << (unsigned long) ceil( (ta + current_tau) / tau) << std::endl;

		for(unsigned long i=0; i < steps; ++i){
			// compute v for Markovian states: shift up to a, we don't make discrete model absorbing
			compute_markovian_vector(discrete_ma,v,u, false);
			// compute u for Probabilistic states
			if(is_imc){
				// if MA is in fact an IMC we can simplify the computation
				compute_interactive_vector(discrete_ma,v,u,max,locks,reach);
			}else {
				compute_probabilistic_vector(discrete_ma,v,u,max, false);
			}
			
			/*
			if(counter==interval_step) {
			
				Real prob;
				if(max)
					prob=0;
				else
					prob=1;
				bool *initials = ma->initials;
				for (unsigned long state_nr = 0; state_nr < num_states; state_nr++) {
					//std::cout << (ma->states_nr.find(state_nr)->second).c_str() << ": " << u[state_nr] << std::endl;
					if(initials[state_nr]){
					if(max){
						if(prob<u[state_nr])
							prob=u[state_nr];
						}else{
						if(prob>u[state_nr])
							prob=u[state_nr];
						}
					}
				}
				
				
				printf("tb=%.5g Maximal time-bounded reachability probability: %.10g\n", ta+i*tau,prob);
			
				counter=0;
			} else {
				counter++;
			}
			*/
			
		}
		SparseMatrix_free(discrete_ma);


	} else { // if a == 0
		// compute discretisation step
		Real tau = compute_error_bound(ma, epsilon,tb);
		// discretize model for the fi
		dbg_printf("discretize model\n");
		SparseMatrix* discrete_ma = discretize_model(ma,tau);
		dbg_printf("model discretized\n");
		//print_model(discrete_ma);

		// value iteration
		unsigned long steps_for_interval = (unsigned long)lround(tb/tau);
		std::cout << "iterations: " << steps_for_interval<< std::endl;
		std::cout << "step duration: " << tau <<std::endl;
		unsigned long interval_step = (unsigned long)lround(interval/tau);
		unsigned long interval_start_point = (unsigned long)lround(interval_start/tau);
		unsigned long counter=0;
		Real tmp_step = interval;
		Real tmp_interval=interval_start;
		std::cout << "interval step: " << interval_step <<std::endl;
		std::cout << "interval start: " << interval_start_point << std::endl;
		
		
		for(unsigned long i=0; i <= steps_for_interval; i++){
			// compute v for Markovian states: from b dwon to a, we make discrete model absorbing
			compute_markovian_vector(discrete_ma,v,u, true);
			// compute u for Probabilistic states
			if(is_imc){
				// if MA is in fact an IMC we can simplify the computation
				compute_interactive_vector(discrete_ma,v,u,max,locks,reach);
			}else {
				compute_probabilistic_vector(discrete_ma,v,u,max, true);
			}
			
			if(i >= interval_start_point){
			if((counter==interval_step || i==interval_start_point) && interval != tb) {
			
				Real prob;
				if(max)
					prob=0;
				else
					prob=1;
				bool *initials = ma->initials;
				for (unsigned long state_nr = 0; state_nr < num_states; state_nr++) {
					//std::cout << (ma->states_nr.find(state_nr)->second).c_str() << ": " << u[state_nr] << std::endl;
					if(initials[state_nr]){
					if(max){
						if(prob<u[state_nr])
							prob=u[state_nr];
						}else{
						if(prob>u[state_nr])
							prob=u[state_nr];
						}
					}
				}
				
				Real tmp = (Real)i*tau - tmp_interval;
				tmp = (Real)i*tau - tmp;
				
				printf("tb=%.5g Maximal time-bounded reachability probability: %.10g  (Real tb=%.5g)\n", tmp,prob,(Real)i*tau);
				
				tmp_interval += tmp_step;
				counter=0;
			} else {
				counter++;
			}
			}
		}
		SparseMatrix_free(discrete_ma);
        delete(discrete_ma);

	}
	// find prob. for initial state and return
	Real prob;
	if(max)
		prob=0;
	else
		prob=1;
	bool *initials = ma->initials;
	printf("vector(%lu) [",num_states);
	for (unsigned long state_nr = 0; state_nr < num_states; state_nr++) {
		//std::cout << (ma->states_nr.find(state_nr)->second).c_str() << ": " << u[state_nr] << std::endl;
		printf(" %lf,",u[state_nr]);
		if(initials[state_nr]){
			if(max){
				if(prob<u[state_nr])
					prob=u[state_nr];
			}else{
				if(prob>u[state_nr])
					prob=u[state_nr];
			}
		}
	}
	printf(" ]\n");
	
	
	return prob;
}
