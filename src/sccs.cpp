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
*	Strongly connected components functions to process an MA
*/

#include "sccs.h"
#include "sccs2.h"

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
 * test if element is in vector
 * 
 */
bool isIn(vector<unsigned long> set, unsigned long element)
{
    vector<unsigned long>::const_iterator e;
    for(e=set.begin(); e != set.end(); e++)
    {
        if((*e) == element)
            return true;
    }
    return false;
}

bool check_if_bad(SparseMatrix *ma,vector<unsigned long> scc_states,unsigned long scc_nr, bool *bad_dist){
	bool new_bad=false;
	
	unsigned long *row_starts = (unsigned long *) ma->row_counts;
	unsigned long *choice_starts = (unsigned long *) ma->choice_counts;
	unsigned long *cols = ma->cols;
	
	for (unsigned long state_nr = 0; state_nr < ma->n; state_nr++) {
		unsigned long state_start = row_starts[state_nr];
		unsigned long state_end = row_starts[state_nr + 1];
		for (unsigned long choice_nr = state_start; choice_nr < state_end; choice_nr++) {
			if(!bad_dist[choice_nr]) {
				/* Add up all outgoing rates of the distribution */
				unsigned long i_start = choice_starts[choice_nr];
				unsigned long i_end = choice_starts[choice_nr + 1];
				unsigned long cmp=scc_states[cols[i_start]];
				for (unsigned long j = i_start+1; j < i_end; j++) {
					if(cmp!=scc_states[cols[j]])
					{
						//printf("%s -> %s\n",(ma->states_nr.find(state_nr)->second).c_str(),(ma->states_nr.find(cols[j])->second).c_str());
						new_bad=true;
						bad_dist[choice_nr]=true;
					}
				}
			}
		}
	}
	
	return new_bad;
}

/**
 * Compute one nontrivial SCC with respect to a set of bad states
 * 
 */
static void strongconnect_weak(SparseMatrix *ma, unsigned long v, vector<long>& index, vector<unsigned long>& lowlink, unsigned long& i, unsigned long& scc_nr,
			       vector<unsigned long>& stack, bool *bad, vector<unsigned long>& scc_states, unsigned long& nr_states, bool& new_bad, bool *bad_dist) {
	// set depth index
	index[v] = i;
	lowlink[v] = i;
	i++;
	stack.push_back(v);
	/*
	printf ("%li ",i);
	if(i==34106)
		printf("test ");
	
	if(i==42459)
		printf("test2 ");
	*/
	
	unsigned long *row_starts = (unsigned long *) ma->row_counts;
	unsigned long *choice_starts = (unsigned long *) ma->choice_counts;
	unsigned long *cols = ma->cols;
	unsigned long dst;
	unsigned long choice_nr,j;
	unsigned long bad_choice = 0;
	unsigned long state_start = row_starts[v];
	unsigned long state_end = row_starts[v + 1];
	unsigned long i_start;
	unsigned long i_end;
	for (choice_nr = state_start; choice_nr < state_end; choice_nr++) {
		/* Add up all outgoing rates of the distribution */
		if(!bad_dist[choice_nr]) {
			i_start = choice_starts[choice_nr];
			i_end = choice_starts[choice_nr + 1];
			for (j = i_start; j < i_end; j++) {
				dst=cols[j];
				//cout << i << endl;
				//printf("%s -> %s\n",(ma->states_nr.find(v)->second).c_str(),(ma->states_nr.find(dst)->second).c_str());
				if(index[(long)dst]< 0 && !bad[dst])
				{
					strongconnect_weak(ma, dst, index, lowlink, i, scc_nr, stack, bad, scc_states, nr_states, new_bad, bad_dist);
					if(lowlink[dst] < lowlink[v])
						lowlink[v]=lowlink[dst];
				}
				else if(isIn(stack,dst) && !bad[dst])
				{
					if(index[dst] < (long)lowlink[v])
						lowlink[v]=index[dst];
				}
				if(bad[dst]){
					j=i_end;
					bad_choice++;
					bad_dist[choice_nr]=true;
				}
			}
		}
	}
	
	
	if(bad_choice == (state_end-state_start)) {
		bad[v]=true;
		new_bad=true;
	}
	
	vector<unsigned long> scc;
	scc.clear();
	
	// consider v is a root node
	if(index[v] == (long)lowlink[v])
	{

		unsigned long w = stack.back();
		stack.pop_back();
		while(w != v)
		{
			scc.push_back(w);
			w = stack.back();
			stack.pop_back();
		}
		scc.push_back(w);
		
		bool selfloop=false;
		
		if(scc.size()==1) {
			state_start = row_starts[scc[0]];
			state_end = row_starts[scc[0] + 1];
			for (choice_nr = state_start; choice_nr < state_end; choice_nr++) {
				/* Add up all outgoing rates of the distribution */
				i_start = choice_starts[choice_nr];
				i_end = choice_starts[choice_nr + 1];
				selfloop=true;
				for (j = i_start; j < i_end; j++) {
					dst=cols[j];
					if(dst!=scc[0])
					{
						selfloop=false;
						break;
					}
				}
			}
			if(!selfloop)
				scc.clear();
		}
	}
	if(!bad[v]){
		for(unsigned long j=0; j<scc.size(); j++) {
			scc_states[scc[j]]=scc_nr;
			nr_states++;
			//cout << ma->states_nr.find(scc[j])->second << " ";
			
		}
		if(scc.size()>0)
			scc_nr++;
	}

}

/**
 * Compute one nontrivial SCC with respect to a set of bad states
 * 
 */
static void strongconnect_strong(SparseMatrix *ma, unsigned long v, vector<long>& index, vector<unsigned long>& lowlink, unsigned long& i, unsigned long& scc_nr,
				 vector<unsigned long>& stack, bool *bad, vector<unsigned long>& scc_states, unsigned long& nr_states, bool& new_bad, bool *bad_dist) {
	// set depth index
	index[v] = i;
	lowlink[v] = i;
	i++;
	stack.push_back(v);
	
	unsigned long *row_starts = (unsigned long *) ma->row_counts;
	unsigned long *choice_starts = (unsigned long *) ma->choice_counts;
	unsigned long *cols = ma->cols;
	unsigned long dst;
	unsigned long choice_nr,j;
	bool bad_choice = false;
		
	unsigned long state_start = row_starts[v];
	unsigned long state_end = row_starts[v + 1];
	unsigned long i_start;
	unsigned long i_end;
	for (choice_nr = state_start; choice_nr < state_end; choice_nr++) {
		/* Add up all outgoing rates of the distribution */
		if(!bad_dist[choice_nr]) {
			i_start = choice_starts[choice_nr];
			i_end = choice_starts[choice_nr + 1];
			for (j = i_start; j < i_end; j++) {
				dst=cols[j];
				//printf("%s -> %s\n",(ma->states_nr.find(v)->second).c_str(),(ma->states_nr.find(dst)->second).c_str());
				if(index[dst]< 0 && !bad[dst])
				{
					strongconnect_strong(ma, dst, index, lowlink, i, scc_nr, stack, bad, scc_states, nr_states,new_bad, bad_dist);
					if(lowlink[dst] < lowlink[v])
						lowlink[v]=lowlink[dst];
				}
				else if(isIn(stack,dst) && !bad[dst])
				{
					if(index[dst] < (long)lowlink[v])
						lowlink[v]=index[dst];
				}
				if(bad[dst]){
					j=i_end;
					bad_choice=true;
					bad_dist[choice_nr]=true;
				}
			}
		}
	}
	
	
	if(bad_choice) {
		bad[v]=true;
		new_bad=true;
	}
	
	vector<unsigned long> scc;
	scc.clear();
	
	// consider v is a root node
	if(index[v] == (long)lowlink[v])
	{

		unsigned long w = stack.back();
		stack.pop_back();
		while(w != v)
		{
			scc.push_back(w);
			w = stack.back();
			stack.pop_back();
		}
		scc.push_back(w);
		
		bool selfloop=false;
		
		if(scc.size()==1) {
			state_start = row_starts[scc[0]];
			state_end = row_starts[scc[0] + 1];
			for (choice_nr = state_start; choice_nr < state_end; choice_nr++) {
				/* Add up all outgoing rates of the distribution */
				i_start = choice_starts[choice_nr];
				i_end = choice_starts[choice_nr + 1];
				for (j = i_start; j < i_end; j++) {
					dst=cols[j];
					if(dst==scc[0])
					{
						selfloop=true;
						break;
					}
				}
			}
			if(!selfloop)
				scc.clear();
		}
	}
	if(!bad[v]) {
		for(unsigned long j=0; j<scc.size(); j++) {
			scc_states[scc[j]]=scc_nr;
			nr_states++;
		}
		if(scc.size()>0)
			scc_nr++;
	}

}

/**
 * Compute BSCCs with respect to a set of bad states
 * 
  * @param ma file to read MA from
 */
SparseMatrixMEC* compute_bottom_strongly_connected_components(SparseMatrix *ma) {
	SparseMatrixMEC *bscc;
	
	unsigned long i;
	unsigned long nr_states=0;
	bool *bad=(bool *) malloc(ma->n * sizeof(bool));
	
	for(i=0; i<ma->n; i++) {
		bad[i]=false;
	}
	
	vector<long> index(ma->n,-1);
	vector<unsigned long> lowlink(ma->n,0);
	vector<unsigned long> stack;
	vector<unsigned long> bscc_states(ma->n,0);
	
	vector<long> indextmp(ma->n,-1);
	vector<unsigned long> lowlinktmp(ma->n,0);
	vector<unsigned long> stacktmp;
	vector<unsigned long> bscc_statestmp(ma->n,0);
	
	unsigned long *dist_starts = (unsigned long *) ma->row_counts;
	unsigned long dist = dist_starts[ma->n];
	
	bool *bad_dist=(bool *) malloc(dist * sizeof(bool));
	for(i=0; i<dist; i++) {
		bad_dist[i]=false;
	}
	
	i = 0;
	bool new_bad=false;
	
	printf("BSCC computation start.\n");
	
	unsigned long idx=0;
	unsigned long scc_nr=1;
	while(idx < ma->n) {
		if(index[idx]<0 && !bad[idx]) {
			strongconnect_strong(ma, idx, index, lowlink, i, scc_nr,stack, bad, bscc_states, nr_states, new_bad, bad_dist);
			if(!new_bad)
				new_bad=check_if_bad(ma,bscc_states,scc_nr,bad_dist);
		}
		idx++;
		if(new_bad)
		{
			index=indextmp;
			lowlink=lowlinktmp;
			stack=stacktmp;
			bscc_states=bscc_statestmp;
			new_bad=false;
			i=0;
			idx=0;
			scc_nr=1;
			nr_states=0;
		}
	}
	
	/* allocate memory for BSCCs and store them */
	bscc=SparseMatrixMEC_new(nr_states,scc_nr-1);
	unsigned long * cols = bscc->cols;
	unsigned long *row_starts = (unsigned long *) bscc->row_counts;
	
	vector<unsigned long>::const_iterator it;
	unsigned long pos=0;
	unsigned long scc_index=0;
	unsigned long scc_size=0;
	unsigned long tmp=0;
	for(idx=0; idx<bscc->n; idx++) {
		i=0;
		for(it=bscc_states.begin(); it != bscc_states.end(); it++)
		{
			if((*it)==idx+1){
				cols[pos]=i;
				pos++;
				tmp++;
			}
			i++;
		}
		if (scc_index > 0) {
			row_starts[scc_index] =
			row_starts[scc_index - 1] + scc_size;
		}
		scc_size = tmp;
		tmp=0;
		scc_index++;
	}
	row_starts[scc_index] = row_starts[scc_index - 1] + scc_size;
	
	
	unsigned long *row = (unsigned long *) bscc->row_counts;
	unsigned long *col = bscc->cols;
	
	for(unsigned long mec_nr=0; mec_nr < bscc->n; mec_nr++) {
		unsigned long mec_start = row[mec_nr];
		unsigned long mec_end = row[mec_nr + 1];
		dbg_printf("MEC %ld: ",mec_nr+1);
		for(unsigned long state_nr=mec_start; state_nr < mec_end; state_nr++) {
			dbg_printf("%s ",ma->states_nr.find(col[state_nr])->second.c_str());
		}
		dbg_printf("\n");
	}
	
	
	free(bad);
	
	return bscc;
}

/**
 * Compute MECs
 * 
 * @param ma file to read MA from
 */
SparseMatrixMEC* compute_maximal_end_components(SparseMatrix *ma) {
	SparseMatrixMEC *mec;
	
	unsigned long i;
	unsigned long nr_states=0;
	bool *bad=(bool *) malloc(ma->n * sizeof(bool));
	
	for(i=0; i<ma->n; i++) {
		bad[i]=false;
	}
	
	vector<long> index(ma->n,-1);
	vector<unsigned long> lowlink(ma->n,0);
	vector<unsigned long> stack;
	vector<unsigned long> mec_states(ma->n,0);
	
	vector<long> indextmp(ma->n,-1);
	vector<unsigned long> lowlinktmp(ma->n,0);
	vector<unsigned long> stacktmp;
	vector<unsigned long> mec_statestmp(ma->n,0);
	
	unsigned long *dist_starts = (unsigned long *) ma->row_counts;
	unsigned long dist = dist_starts[ma->n];
	
	bool *bad_dist=(bool *) malloc(dist * sizeof(bool));
	for(i=0; i<dist; i++) {
		bad_dist[i]=false;
	}
	
	i = 0;
	bool new_bad=false;
	
	printf("MEC computation start.\n");
	
	unsigned long idx=0;
	unsigned long scc_nr=1;
	while(idx < ma->n) {
		if(index[idx]<0 && !bad[idx]) {
			// cout << "start with " << idx << endl;
			strongconnect_weak(ma, idx, index, lowlink, i, scc_nr, stack, bad, mec_states, nr_states, new_bad, bad_dist);
			if(!new_bad)
				new_bad=check_if_bad(ma,mec_states,scc_nr,bad_dist);
		}
		idx++;
		if(new_bad)
		{
			index=indextmp;
			lowlink=lowlinktmp;
			stack=stacktmp;
			mec_states=mec_statestmp;
			new_bad=false;
			i=0;
			idx=0;
			scc_nr=1;
			nr_states=0;
		}
	}
	
	/*
	if(!new_bad)
		new_bad=check_if_bad(ma,mec_states,scc_nr,bad_dist);
	for(unsigned long x=0; x<ma->n; x++) {
		if(bad_dist[x])
			printf("bad choice %ld\n",x);
	}
	*/
	
	/* allocate memory for BSCCs and store them */
	mec=SparseMatrixMEC_new(nr_states,scc_nr-1);
	unsigned long * cols = mec->cols;
	unsigned long *row_starts = (unsigned long *) mec->row_counts;
	
	vector<unsigned long>::const_iterator it;
	unsigned long pos=0;
	unsigned long scc_index=0;
	unsigned long scc_size=0;
	unsigned long tmp=0;
	for(idx=0; idx<mec->n; idx++) {
		i=0;
		for(it=mec_states.begin(); it != mec_states.end(); it++)
		{
			if((*it)==idx+1){
				cols[pos]=i;
				pos++;
				tmp++;
			}
			i++;
		}
		if (scc_index > 0) {
			row_starts[scc_index] =
			row_starts[scc_index - 1] + scc_size;
		}
		scc_size = tmp;
		tmp=0;
		scc_index++;
	}
	row_starts[scc_index] = row_starts[scc_index - 1] + scc_size;
	
	
	unsigned long *row = (unsigned long *) mec->row_counts;
	unsigned long *col = mec->cols;
	
	for(unsigned long mec_nr=0; mec_nr < mec->n; mec_nr++) {
		unsigned long mec_start = row[mec_nr];
		unsigned long mec_end = row[mec_nr + 1];
		dbg_printf("MEC %ld: ",mec_nr+1);
		for(unsigned long state_nr=mec_start; state_nr < mec_end; state_nr++) {
			//if(ma->goals[col[state_nr]])
				dbg_printf("%s ",ma->states_nr.find(col[state_nr])->second.c_str());
		}
		dbg_printf("\n");
	}
	
	
	free(bad);
	
	return mec;
}

/**
* Computes unbounded reachability for MA.
*
* @param locks identifier if state can't reach a goal state
* @param ma file to read MA from
*/
bool* compute_locks_strong(SparseMatrix *ma) {
	bool *goals=ma->goals;
	bool *bad=(bool *) malloc(ma->n * sizeof(bool));
	
	for(unsigned long i=0; i<ma->n; i++) {
		bad[i]=false;
		if(goals[i])
			bad[i]=true;
	}

	bool *tmp=compute_locks_strong(ma, bad);
	free(bad);

	return tmp;
}

/**
* Computes unbounded reachability for MA.
*
* @param locks identifier if state can't reach a goal state
* @param ma file to read MA from
*/
bool* compute_locks_strong(SparseMatrix *ma, bool* bad) {
	bool *locks = (bool *) malloc(ma->n * sizeof(bool));
	
	unsigned long i;
	unsigned long nr_states=0;
	
	for(i=0; i<ma->n; i++) {
		locks[i]=false;
	}
	
	vector<long> index(ma->n,-1);
	vector<unsigned long> lowlink(ma->n,0);
	vector<unsigned long> stack;
	vector<unsigned long> lock_states(ma->n,0);
	
	vector<long> indextmp(ma->n,-1);
	vector<unsigned long> lowlinktmp(ma->n,0);
	vector<unsigned long> stacktmp;
	vector<unsigned long> lock_statestmp(ma->n,0);
	
	unsigned long *row_starts = (unsigned long *) ma->row_counts;
	unsigned long dist = row_starts[ma->n];
	
	bool *bad_dist=(bool *) malloc(dist * sizeof(bool));
	for(i=0; i<dist; i++) {
		bad_dist[i]=false;
	}
	
	
	//unsigned long scc_nr=1;
	//compute_SCC_decomposition_tarjan(ma, lock_states, bad, bad_dist, scc_nr);
	
	
	i = 0;
	bool new_bad=false;
	dbg_printf("SCC strong computation start.\n");
	unsigned long idx=0;
	unsigned long scc_nr=1;
	while(idx < ma->n) {
		if(index[idx]<0 && !bad[idx]) {
			strongconnect_strong(ma, idx, index, lowlink, i, scc_nr, stack, bad, lock_states, nr_states, new_bad, bad_dist);
			//function_strongconnect(ma,v,lock_states,stack,);
			//dbg_printf("check.\n");
			if(!new_bad){
				new_bad=check_if_bad(ma,lock_states,scc_nr,bad_dist);
			}
		}
		idx++;
		//dbg_printf("idx: %d.\n",idx);
		if(new_bad)
		{
			dbg_printf("new.\n");
			index=indextmp;
			lowlink=lowlinktmp;
			stack=stacktmp;
			lock_states=lock_statestmp;
			new_bad=false;
			i=0;
			idx=0;
			scc_nr=1;
			nr_states=0;
		}
	}
	
	unsigned long *choice_starts = (unsigned long *) ma->choice_counts;
	bool check = true;
	bool isLock = false;
	bool strong_lock = false;
	while(check) {
		//cout << "enter" << endl;
		check = false;
		for (unsigned long state_nr = 0; state_nr < ma->n; state_nr++) {
			if(lock_states[state_nr] == 0  and !ma->goals[state_nr]) {
				unsigned long state_start = row_starts[state_nr];
				unsigned long state_end = row_starts[state_nr + 1];
				dbg_printf("row_starts: %li row_ends: %li\n",state_start,state_end);
				strong_lock = true;
				for (unsigned long choice_nr = state_start; choice_nr < state_end; choice_nr++) {
					// Add up all outgoing rates of the distribution 
					unsigned long i_start = choice_starts[choice_nr];
					unsigned long i_end = choice_starts[choice_nr + 1];
					dbg_printf("choice_starts: %li choice_ends: %li\n",i_start,i_end);
					isLock = true;
					for (unsigned long i = i_start; i < i_end; i++) {
						if(!lock_states[ma->cols[i]]>0) {
							isLock = false;
						}
					}
					if(!isLock)
						strong_lock = false;
				}
				if(strong_lock) {
					//cout << state_nr << endl;
					lock_states[state_nr] = 1;
					check = true;
				}
			}
		}
	}
	
	vector<unsigned long>::const_iterator it;
	i=0;
	for(it=lock_states.begin(); it != lock_states.end(); it++)
	{
		if((*it)>0){
			locks[i]=true;
			dbg_printf("%s ",(ma->states_nr.find(i)->second).c_str());
		}
		i++;
	}
	
	dbg_printf("\n");
	
	//free(bad);
	free(bad_dist);
	
	return locks;
}

/**
* Computes unbounded reachability for MA.
*
* @param locks identifier if state can't reach a goal state
* @param ma file to read MA from
*/
bool* compute_locks_weak(SparseMatrix *ma) {
	bool *goals=ma->goals;
	bool *bad=(bool *) malloc(ma->n * sizeof(bool));
	
	for(unsigned long i=0; i<ma->n; i++) {
		bad[i]=false;
		if(goals[i])
			bad[i]=true;
	}
	
	bool* tmp=compute_locks_weak(ma, bad);
	free(bad);

	return tmp;
}

/**
* Computes unbounded reachability for MA.
*
* @param locks identifier if state can't reach a goal state
* @param ma file to read MA from
*/
bool* compute_locks_weak(SparseMatrix *ma, bool* bad) {
	bool *locks = (bool *) malloc(ma->n * sizeof(bool));
	
	unsigned long i;
	unsigned long nr_states=0;
	
	for(i=0; i<ma->n; i++) {
		locks[i]=false;
	}
	
	vector<long> index(ma->n,-1);
	vector<unsigned long> lowlink(ma->n,0);
	vector<unsigned long> stack;
	vector<unsigned long> lock_states(ma->n,0);
	
	vector<long> indextmp(ma->n,-1);
	vector<unsigned long> lowlinktmp(ma->n,0);
	vector<unsigned long> stacktmp;
	vector<unsigned long> lock_statestmp(ma->n,0);
	
	unsigned long *row_starts = (unsigned long *) ma->row_counts;
	unsigned long dist = row_starts[ma->n];
	
	bool *bad_dist=(bool *) malloc(dist * sizeof(bool));
	for(i=0; i<dist; i++) {
		bad_dist[i]=false;
	}
	
	i = 0;
	bool new_bad=false;
	dbg_printf("SCC weak computation start.\n");
	unsigned long idx=0;
	unsigned long scc_nr=1;
	while(idx < ma->n) {
		if(index[idx]<0 && !bad[idx]) {
			strongconnect_weak(ma, idx, index, lowlink, i, scc_nr, stack, bad, lock_states, nr_states, new_bad, bad_dist);
			if(!new_bad)
				new_bad=check_if_bad(ma,lock_states,scc_nr,bad_dist);
		}
		idx++;
		if(new_bad)
		{
			index=indextmp;
			lowlink=lowlinktmp;
			stack=stacktmp;
			lock_states=lock_statestmp;
			new_bad=false;
			i=0;
			idx=0;
			scc_nr=1;
			nr_states=0;
		}
	}
	
	unsigned long *choice_starts = (unsigned long *) ma->choice_counts;
	bool check = true;
	bool isLock = false;
	bool strong_lock = false;
	while(check) {
		//cout << "enter" << endl;
		check = false;
		for (unsigned long state_nr = 0; state_nr < ma->n; state_nr++) {
			if(lock_states[state_nr] == 0 and !ma->goals[state_nr]) {
				unsigned long state_start = row_starts[state_nr];
				unsigned long state_end = row_starts[state_nr + 1];
				strong_lock = false;
				for (unsigned long choice_nr = state_start; choice_nr < state_end; choice_nr++) {
					// Add up all outgoing rates of the distribution 
					unsigned long i_start = choice_starts[choice_nr];
					unsigned long i_end = choice_starts[choice_nr + 1];
					dbg_printf("choice_starts: %li choice_ends: %li\n",i_start,i_end);
					isLock = true;
					for (unsigned long i = i_start; i < i_end; i++) {
						if(!lock_states[ma->cols[i]]>0) {
							isLock = false;
						}
					}
					if(isLock)
						strong_lock = true;
				}
				if(strong_lock) {
					//cout << state_nr << endl;
					lock_states[state_nr] = 1;
					check = true;
				}
			}
		}
	}
	
	vector<unsigned long>::const_iterator it;
	i=0;
	for(it=lock_states.begin(); it != lock_states.end(); it++)
	{
		if((*it)>0){
			locks[i]=true;
			dbg_printf("%s ",(ma->states_nr.find(i)->second).c_str());
		}
		i++;
	}
	
	dbg_printf("\n");
	
    free(bad_dist);
	
	return locks;
	
}
