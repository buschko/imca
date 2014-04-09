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


#include "bounded_reward.h"
#include "read_file.h"
#include "debug.h"
#include "sccs.h"
#include <math.h>
#include <vector>

/**
* sets the error bound for given epsilon and tb
*
* @param ma the MA
* @param epsilon epsilon bound
* @param tb time bound
*/
Real compute_error_bound_reward(SparseMatrix* ma, Real epsilon, Real tb) {
	Real tau;
	// TODO: add lambda to read in function
	Real lambda=ma->max_exit_rate;
	Real rho = ma -> max_markovian_reward;
	tau = (2*epsilon)/(tb*lambda*(rho+1));
	return tau;
}

/**
* discretizes the MA
*
* @param ma the MA
* @param tau discretization factor
* @return discretized MA
*/
SparseMatrix* discretize_model_reward(SparseMatrix* ma, Real tau) {
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
	Real *rewards = discrete_ma->rewards;
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
				rewards[choice_nr] = ma->rewards[choice_nr] / exit_rate * exp_estau_com; // Reward of each step for Markovian state s is: Rew(s)/E(s)*(1-exp(-E(s)*tau))
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
				// Write the reward into discrete model
				rewards[choice_nr] = ma->rewards[choice_nr];
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
void compute_markovian_vector_with_reward(SparseMatrix* ma, vector<Real>& v, const vector<Real> &u, bool is_MA_made_absorbing){
	unsigned long *row_starts = (unsigned long *) ma->row_counts;
	unsigned long *choice_starts = (unsigned long *) ma->choice_counts;
	unsigned long *cols = ma->cols;
	Real* non_zeros = ma->non_zeros;
	Real* rewards = ma->rewards;
	bool *goals = ma->goals;
	for (unsigned long state_nr = 0; state_nr < ma->n; state_nr++) {
		unsigned long state_start = row_starts[state_nr];
		unsigned long state_end = row_starts[state_nr + 1];
		// Look at Markovian states
//		if(goals[state_nr] && is_MA_made_absorbing){
//			v[state_nr]=1;
//		}else{
			if(!ma->isPS[state_nr])
			{
				for (unsigned long choice_nr = state_start; choice_nr < state_end; choice_nr++) {
					// Add up all outgoing rates of the distribution
					unsigned long i_start = choice_starts[choice_nr];
					unsigned long i_end = choice_starts[choice_nr + 1];
					v[state_nr] = rewards[choice_nr];
					for (unsigned long i = i_start; i < i_end; i++) {
						v[state_nr] +=  non_zeros[i] * u[cols[i]];
					}
				}
			}else {
					v[state_nr] = u[state_nr];
			}
//		}
	}
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
void compute_probabilistic_vector_with_reward(SparseMatrix* ma, vector<Real>& v, vector<Real>& u, bool max, bool is_MA_made_absorbing){

	const Real precision = 1e-7;
	unsigned long statecount = ma-> n;
	unsigned long *row_starts = (unsigned long *) ma->row_counts;
	unsigned long *choice_starts = (unsigned long *) ma->choice_counts;
	unsigned long *cols = ma->cols;
	Real* non_zeros = ma->non_zeros;
	Real* rewards = ma->rewards;

	// Initialization
	for (unsigned long s_idx = 0; s_idx < statecount; s_idx++) {
		if ( ma -> isPS[s_idx] && (!is_MA_made_absorbing || (is_MA_made_absorbing && !ma -> goals[s_idx])) )
			v[s_idx] = 0.0;
		else
			u[s_idx] = v[s_idx];

	}
	//long int bestChoice=0;

	// main loop
	bool done = false;
	// done will set to true inside the while loop when we are close enough to the fixed point
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
					// add reward
					//tmp += rewards[choice_nr];
					// Add up all outgoing rates of the distribution
					unsigned long i_start = choice_starts[choice_nr];
					unsigned long i_end = choice_starts[choice_nr + 1];
					for (unsigned long i = i_start; i < i_end; i++) {
						tmp += non_zeros[i] * v[cols[i]];
					}
					if( max ) {
						if(tmp > u[s_idx] ){
//							bestChoice=choice_nr;
							u[s_idx] = tmp;
						}
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

		//std::cout<<bestChoice<<std::endl;

	}



}

/**
* computes accumulated gained between time @param ta to @param tb
*
* @param ma the MA
* @param max maximum/minimum
* @param epsilon the given error
* @param tb the given time bound
* @param is_imc indicates if MA is an IMC
*/
Real compute_time_bounded_accumulated_reward(SparseMatrix* ma, bool max, Real epsilon, Real ta, Real tb, bool is_imc, Real interval,Real interval_start) {

	// Check whether the given time interval is zero
	if( ta > tb ) {
		printf("WARNING: The given interval is empty (upper bound < lower bound.) The accumulated reward within an empty interval is ZERO.\n");
		return 0.0;
	}

	unsigned long num_states = ma->n;
	vector<Real> v(num_states,0); // Markovian vector
	vector<Real> u(num_states,0); // Probabilistic vector
	// initialize goal states
	/*
	bool *goals = ma->goals;
	for (unsigned long state_nr = 0; state_nr < num_states; state_nr++) {
		if(goals[state_nr]){
			v[state_nr]=1;
		}
	}
	*/
	// in case MA is an IMC: precomputation of paths for interactive states
	vector< vector<unsigned long> > reach;
	
	cout << "start value iteration" << endl;
	if( ta > 0 ) {

		std::cout<<"***** WARNING *****\nThis computation mode is under development and might give you wrong answer.\n\n";
		//TODO one unique function that tell you what should be tau1 and what tau2
		// compute the upper bound of discretization step
		Real tau = compute_error_bound_reward(ma, epsilon,tb);

		unsigned long steps; // stores the number of steps should be taken
		Real current_tau;    // stores the current value of tau which is calculated to make sure the number of steps is a integral value
		SparseMatrix* discrete_ma; // discretized model induced from the original ma

		steps = (unsigned long) ceil((tb - ta) / tau); // calculating the number of steps for interval tb down to ta
		current_tau = (tb - ta) / steps;               // recalculate tau based on the number of steps


		// discretize model with respect to the current value of tau 'current_tau'
		dbg_printf("discretize model for interval [%g,%g] ... \n", ta, tb);
		discrete_ma = discretize_model_reward(ma,current_tau);
		dbg_printf("model discretized\n");
		//print_model(discrete_ma);
		
		//unsigned long interval_step = round(interval/current_tau);
		//unsigned long counter=0;
		//cout << "interval step: " << interval_step <<endl;


		// value iteration
		cout << "iterations: " << (unsigned long) ceil((tb - ta) / tau) + (unsigned long) ceil(ta / tau) << endl;
		printf("step duration for interval [%g,%g]: %g\n", ta, tb, current_tau);
		
		for(unsigned long i=0; i < steps; i++){
			// compute v for Markovian states: from b dwon to a, we make discrete model absorbing
			compute_markovian_vector_with_reward(discrete_ma,v,u, true);
			// compute u for Probabilistic states
			compute_probabilistic_vector_with_reward(discrete_ma,v,u,max, true);
			
			/*
			if(counter==interval_step) {
			
				Real prob;
				if(max)
					prob=0;
				else
					prob=1;
				bool *initials = ma->initials;
				for (unsigned long state_nr = 0; state_nr < num_states; state_nr++) {
					//cout << (ma->states_nr.find(state_nr)->second).c_str() << ": " << u[state_nr] << endl;
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
		current_tau = (ta + current_tau) / steps;                                // recalculate tau based on the number of steps

		// discretize model with respect to the current value of tau 'current_tau'
		dbg_printf("discretize model for interval [0,%g] ... \n", ta);
		discrete_ma = discretize_model_reward(ma,current_tau);
		dbg_printf("model discretized\n");
		//print_model(discrete_ma);

		printf("step duration for interval [0,%g]: %g\n", ta, current_tau);

		for(unsigned long i=0; i < steps; i++){
			// compute v for Markovian states: shift up to a, we don't make discrete model absorbing
			compute_markovian_vector_with_reward(discrete_ma,v,u, false);
			// compute u for Probabilistic states
			compute_probabilistic_vector_with_reward(discrete_ma,v,u,max, false);
			
			/*
			if(counter==interval_step) {
			
				Real prob;
				if(max)
					prob=0;
				else
					prob=1;
				bool *initials = ma->initials;
				for (unsigned long state_nr = 0; state_nr < num_states; state_nr++) {
					//cout << (ma->states_nr.find(state_nr)->second).c_str() << ": " << u[state_nr] << endl;
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
		Real tau = compute_error_bound_reward(ma, epsilon,tb);
		unsigned long steps_for_interval = round(tb/tau);
		// update tau
		tau = tb / steps_for_interval;
		// discretize model with respect to the given epsilon
		dbg_printf("discretize model\n");
		SparseMatrix* discrete_ma = discretize_model_reward(ma,tau);
		print_model(discrete_ma,true);
		std::cout<<"Exit Rate: "<<ma -> max_exit_rate <<"Reward: "<<ma -> max_markovian_reward;
		dbg_printf("model discretized\n");
		//print_model(discrete_ma);

		// value iteration
		cout << "iterations: " << steps_for_interval<< endl;
		cout << "step duration: " << tau <<endl;
		unsigned long interval_step = round(interval/tau);
		unsigned long interval_start_point = round(interval_start/tau);
		unsigned long counter=0;
		Real tmp_step = interval;
		Real tmp_interval=interval_start;
		cout << "interval step: " << interval_step <<endl;
		cout << "interval start: " << interval_start_point << endl;
		
		// Note: we start the computation from step one, since in contrast to time bounded reachability, the initial vector here is zero.
		for(unsigned long i=1; i <= steps_for_interval; i++){
			// compute v for Markovian states for the current step; we make discrete model absorbing. After this step v contains the updated reward value
			compute_markovian_vector_with_reward(discrete_ma,v,u, true);
			// compute u for Probabilistic states
			compute_probabilistic_vector_with_reward(discrete_ma,v,u,max, true);
			
			if(i >= interval_start_point){
			if((counter==interval_step || i==interval_start_point) && interval != tb) {
			
				Real prob;
				if(max)
					prob=0;
				else
					prob=1;
				bool *initials = ma->initials;
				for (unsigned long state_nr = 0; state_nr < num_states; state_nr++) {
					//cout << (ma->states_nr.find(state_nr)->second).c_str() << ": " << u[state_nr] << endl;
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
				
				Real tmp = i*tau - tmp_interval;
				tmp = i*tau - tmp;
				
				printf("tb=%.5g Maximal time-bounded accumulated reward: %.10g  (Real tb=%.5g)\n", tmp,prob,i*tau);
				
				tmp_interval += tmp_step;
				counter=0;
			} else {
				counter++;
			}
			}
		}
		SparseMatrix_free(discrete_ma);

	}
	// find prob. for initial state and return
	Real prob;
	if(max)
		prob=0;
	else
		prob=1;
	bool *initials = ma->initials;
	for (unsigned long state_nr = 0; state_nr < num_states; state_nr++) {
		//cout << (ma->states_nr.find(state_nr)->second).c_str() << ": " << u[state_nr] << endl;
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
	
	
	return prob;
}
