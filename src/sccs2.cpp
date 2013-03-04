/**
 * New version of the sccs.cpp for an QMA project: Improving MEC Decomposition
 * in the IMCA tool. This class contains a the old, and the two new MEC Decomposition
 * algorithms as proposed in a paper by Chatterjee and Henzinger: Faster and Dynamic
 * Algorithms For Maximal End-Component Decomposition And Related Graph Problems In
 * Probabilistic Verfication.
 * 
 * Author: Daan van Beek (s0167789)
 */
#include "sccs2.h"
#include "sccs.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <map>
#include <string>
#include <vector>
#include <ctime>
#include <time.h>
#include <stack>
#include <pthread.h>

#include "debug.h"

using namespace std;

/**
 * test if element is in vector
 * 
 */
bool isIn_copy(vector<unsigned long> set, unsigned long element)
{
    vector<unsigned long>::const_iterator e;
    for(e=set.begin(); e != set.end(); e++)
    {
        if((*e) == element)
            return true;
    }
    return false;
}

void function_strongconnect(SparseMatrix *ma,unsigned long v, vector<unsigned long>& scc_states, vector<unsigned long>& stack, bool *stack_precense, bool* bad_states, bool* bad_transitions, unsigned long& i, vector<long>& index, vector<unsigned long>& lowlink,unsigned long& scc_nr){
	// Set the depth index for v to the smallest unused index
    index[v] = i; //v.index := index
    lowlink[v] = i; //v.lowlink := index
    i++; //index := index + 1
    stack.push_back(v); //S.push(v) //stack.push_back(v);
    stack_precense[v] = true;
	
	//printf("size: %ld\n",stack.size());
	//printf("capacity: %ld\n",stack.capacity());
	//printf("max size: %ld\n",stack.max_size());
	dbg_printf("function_strongconnect for v: %lu.\n",v);
    dbg_printf("size: %ld\n",stack.size());
	dbg_printf("capacity: %ld\n",stack.capacity());
	dbg_printf("max size: %ld\n",stack.max_size());
    
    //printf("function_strongconnect for v: %lu.\n",v);
	// Consider successors of v
    //for each (v,w) in E do:
    unsigned long *row_starts = (unsigned long *) ma->row_counts;
    unsigned long *choice_starts = (unsigned long *) ma->choice_counts;
    unsigned long *cols = ma->cols;
    unsigned long dst;
    unsigned long row_start = row_starts[v]; //row_start = row_counts[i]
    unsigned long row_end = row_starts[v + 1]; //row_end = row_counts[i+1]
    unsigned long choice_start;
    unsigned long choice_end;
    while(row_start < row_end) {
        //printf("a");
        if(!bad_transitions[row_start]){
            choice_start = choice_starts[row_start];
            choice_end = choice_starts[row_start + 1];
            //printf("b");
            while(choice_start < choice_end){ 
                //printf("c");
                dst=cols[choice_start];
                if(index[(long)dst]< 0 && !bad_states[dst]){
					if(v==47360)
						printf("here\n");
                    //printf("d");
                    // Successor w has not yet been visited; recurse on it
                    function_strongconnect(ma,dst,scc_states,stack,stack_precense,bad_states,bad_transitions,i,index,lowlink,scc_nr); //stronglyconnect(w)
                    if(lowlink[dst] < lowlink[v]){
                        lowlink[v] = lowlink[dst]; //v.lowlink := min(v.lowlink, w.lowlink)
                    }
                }else if(stack_precense[dst] && !bad_states[dst]){ //isIn_copy(stack,dst)  is way too expensive! made the algorithm 200 times slower
                    //printf("e");
                    // Successor W is in stack S and hence in the current SCC
                    if(index[dst] < (long)lowlink[v]){
                        lowlink[v]=index[dst]; //v.lowlink := min(v.lowlink, v.index)
                    }
                }
                choice_start++;
            } 
        }
        row_start++;
    }
    
    vector<unsigned long> scc;
    scc.clear();
    
    // If v is a root node, pop the stack and generate an SCC
    if(index[v] == (long)lowlink[v]){ // if (v.lowlink = v.index) then
        //start a new strongly connected component
        unsigned long w = stack.back();
        stack.pop_back();
        stack_precense[w] = false;
        
        while(w != v)
        {
                scc.push_back(w);
                w = stack.back();
                stack.pop_back();
                stack_precense[w] = false;
        }
        scc.push_back(w);
        

        bool selfloop=false;

        if(scc.size()==1) {
            row_start = row_starts[scc[0]];
            row_end = row_starts[scc[0] + 1];
            while (row_start < row_end) {
                /* Add up all outgoing rates of the distribution */
                choice_start = choice_starts[row_start];
                choice_end = choice_starts[row_start + 1];
                selfloop=true;
                while (choice_start < choice_end) {
                    dst=cols[choice_start];
                    if(dst!=scc[0])
                    {
                            selfloop=false;
                            break;
                    }
                    choice_start++;
                }
                row_start++;
            }
            if(!selfloop)
                scc.clear();
        }            
            
        //output the current strongly connected component
        for(unsigned long j=0; j<scc.size(); j++) {
            scc_states[scc[j]]=scc_nr;		
        }
            
	if(scc.size() > 0){
            scc_nr++;
            dbg_printf("SCC found by decomposition!");
        }
    }
}


/**
 * Computes Tarjan SCC decomposition with respect to a set of bad states
 * @param ma
 * @param answer, A return vector of length ma->n, with per state the SCC it belongs to.
 */
void compute_SCC_decomposition_tarjan(SparseMatrix *ma, vector<unsigned long>& scc_states, bool* bad_states, bool* bad_transitions, unsigned long& scc_nr){
    /* Initialize vectors needed for Tarjan SCC Decomposition algorithm */
    vector<long> index(ma->n,-1); //init v.index
    vector<unsigned long> lowlink(ma->n,0); //init v.lowlink
	
    unsigned long i=0; //index := 0
    vector<unsigned long> stack; //S := empty
    bool *stack_precense=(bool *) malloc(ma->n * sizeof(bool));
    for(unsigned long i = 0; i < ma->n; i++){
        stack_precense[i]=false;
    }
    unsigned long v = 0;
	
    while(v < ma->n){ //for each v in V do
        
        if(index[v]<0 && !bad_states[v]){
            function_strongconnect(ma,v,scc_states, stack, stack_precense, bad_states, bad_transitions,i,index,lowlink,scc_nr); 
        }
        v++;
    }
	
	free(stack_precense);
    
    /*
    for(int j = 0; j < ma->n;j++){
        printf("state %i has number %lu.\n",j,scc_states[j]);
    }*/
    
}

/**
 * Checks if the given SCC is an MEC, answer is the violating set of vertices
 * (which have an random edge leaving the scc).
 * The violating boolean array signifies the result.
 */
void check_mec(SparseMatrix *ma, bool* bad_states, bool* bad_transitions, bool* violating, vector<unsigned long>& scc_states,unsigned long& scc_nr, bool& bad, unsigned long& nr_violating){
    ////printf("check MEC begin for scc_nr %lu\n",scc_nr);
    for(unsigned long i = 0; i < ma->n; i++){
        violating[i] = false;
    }
    /*printf("scc array for scc: %lu, ma->n: %lu\n",scc_nr,ma->n);
    for(unsigned long j = 0; j < ma->n; j++){
        printf("%lu",scc_states[j]);
    }*/
     
    //for each vertex in SCC:
    for(unsigned long v = 0; ((v < ma->n)); v++){ //(v < ma->n) && (scc_states[v] == scc_nr) does not work??
        if((scc_states[v] == scc_nr)){
            //lookup the edges going from this vertex:
            unsigned long *row_starts = (unsigned long *) ma->row_counts;
            unsigned long *choice_starts = (unsigned long *) ma->choice_counts;
            unsigned long *cols = ma->cols;
            unsigned long dst;
            unsigned long row_start = row_starts[v]; //row_start = row_counts[i]
            unsigned long row_end = row_starts[v + 1]; //row_end = row_counts[i+1]
            unsigned long choice_start;
            unsigned long choice_end;


            bool good_transition_found = false;
            bool bad_transition_found = false; //a transition with all choices outside SCC is NOT a bad transition, as its removal does not break up the SCC as otherwise the SCC C would not have been an SCC to begin with.
            bool all_choices_in_SCC = true;
            bool all_choices_outside_SCC = true;
            
            while(row_start < row_end) {
                if(!bad_transitions[row_start]){
                    all_choices_in_SCC = true;
                    all_choices_outside_SCC = true;

                    choice_start = choice_starts[row_start];
                    choice_end = choice_starts[row_start + 1];
                    while(choice_start < choice_end){ 
                        dst=cols[choice_start];
                        ////printf("Transition from state %lu to state %lu. Bad dst state: %s. scc_nr: %lu\n",v,dst,(bad_states[dst])?"true":"false",scc_nr);
                        //if(!bad_states[dst]){
                            ////printf("Good dst state!\n");
                            if(scc_states[dst] == scc_nr){ //this comparison gives a problem if put in the same condition as !bad_states[dst]
                                all_choices_outside_SCC = false;
                                ////printf("all_choices_outside_SCC = false;\n");
                            }else if (scc_states[dst] != scc_nr){
                                all_choices_in_SCC = false;
                                ////printf("all_choices_in_SCC = false;\n");
                            }
                        //}else{
                            //may be possible to mark transition here as bad for efficiency
                            ////printf("Bad dst state! Bad state %lu is: %s.\n",dst,(bad_states[dst])?"true":"false");
                        //}
                        choice_start++;
                    }

                    if(all_choices_in_SCC){
                        good_transition_found = true;
                    }else if(!all_choices_outside_SCC){ //((!all_choices_in_SCC) && (!all_choices_outside_SCC))
                        bad_transition_found = true;
                        bad_transitions[row_start] = true;
                    }
                }
                row_start++;
            }

            if(bad_transition_found && !good_transition_found){
                //only bad transitions found, vertex is violating
                violating[v] = true;
                nr_violating++;
                ////printf("scc %lu contained a violating vertex\n",scc_nr);
            }

            if(bad_transition_found && good_transition_found){
                //bad and good transitions found, bad transitions where already marked
                //SCC decomposition must be redone without the questionable transitions
                bad = true;
                ////printf("A bad transition was found in scc %lu\n",scc_nr);
            }
        }
    }
    ////printf("check_scc done!\n");
}

void attractor_remove(SparseMatrix *ma, bool* u0_states, bool* bad_states, bool* new_bad_states,bool* bad_transitions,vector<unsigned long> scc_states,unsigned long scc_nr,unsigned long& nr_states_removed,bool remove_scc){
    //This function has two modes, separated by boolean scc_remove
    //remove_SCC == true means that the SCC defined by
    //scc_states and scc_nr and its attractor set must be removed,
    //remove_SCC == true then u0_states is the starting set and that only
    // u0_states + its attractor set UNION SCC defined by scc_states and scc_nr must be removed.
    
    //The states that should be removed from bad_states are instead removed from
    //new_bad_states, because if already removed from bad_States (on the end of the function)
    //then the next mec_check would skip those states even though they where not bad in the
    //previous iteration of the whole algorithm. So the bad states are only applied once the
    //iteration is complete.
    
    //initializing the attractor set:
    /*printf("bad_states array:\n");
    for(int j = 0; j < ma->n; j++){
        printf("|%s",bad_states[j]?"true":"false");
    }
    printf("|\n");

    printf("violating array:\n");
    for(int j = 0; j < ma->n; j++){
        printf("|%s",u0_states[j]?"true":"false");
    }
    printf("|\n");
    
    printf("scc_states array:\n");
    for(int j = 0; j < ma->n; j++){
        printf("|%lu",scc_states[j]);
    }
    printf("| and scc_nr: %lu.\n",scc_nr);*/
    
    bool *attractor_set=(bool *) malloc(ma->n * sizeof(bool));
    
    if(remove_scc){
        //initial set is the SCC
        for(int i = 0; i < ma->n; i++){
            if(scc_states[i] == scc_nr){
                attractor_set[i] = true;
            }else{
                attractor_set[i] = false;
            }
        }
    }else{
        for(int i = 0; i < ma->n; i++){
                attractor_set[i] = u0_states[i];           
        }
    }
    
    bool attracted_new_state = true;
    
    while(attracted_new_state){ //continue until no new states have been found
        attracted_new_state = false;
        
        for(unsigned long v = 0; v < ma->n; v++){
            if(!attractor_set[v] && !bad_states[v]){
                //lookup the edges going from this vertex:
                unsigned long *row_starts = (unsigned long *) ma->row_counts;
                unsigned long *choice_starts = (unsigned long *) ma->choice_counts;
                unsigned long *cols = ma->cols;
                unsigned long dst;
                unsigned long row_start = row_starts[v]; //row_start = row_counts[i]
                unsigned long row_end = row_starts[v + 1]; //row_end = row_counts[i+1]
                unsigned long choice_start;
                unsigned long choice_end;


                bool each_transition_choice_to_ui = true;

                while(row_start < row_end) {
                    if(!bad_transitions[row_start]){
                        bool choice_to_ui_found = false;

                        choice_start = choice_starts[row_start];
                        choice_end = choice_starts[row_start + 1];
                        while(choice_start < choice_end){
                            dst=cols[choice_start];
                            //printf("Investigating transition from %lu to %lu\n",v,dst);

                            if(!bad_states[dst] && attractor_set[dst]){
                                choice_to_ui_found = true;
                                //printf("choice_to_ui_found\n");
                            }

                            choice_start++;
                        }

                        if(!choice_to_ui_found){
                            each_transition_choice_to_ui = false;
                            //printf("each_transition_choice_to_ui = false\n");
                        }

                    }
                    row_start++;
                }

                if(each_transition_choice_to_ui){
                    attractor_set[v] = true;
                    attracted_new_state = true;
                }
            }
        }
                
    }
    
    
    
    //now throw away all vertices attracted
    if(remove_scc){
        for(unsigned long i = 0; i < ma->n; i++){
            if(attractor_set[i] && !new_bad_states[i]){
                new_bad_states[i] = true;
                nr_states_removed++;
            }
        }
    }else{
        for(unsigned long i = 0; i < ma->n; i++){
            if(attractor_set[i] && (scc_states[i] == scc_nr) && !new_bad_states[i]){
                new_bad_states[i] = true;
                nr_states_removed++;
            }
        }
    }
    
    /*printf("attractor_set:\n");
    for(int j = 0; j < ma->n; j++){
        printf("|%s",attractor_set[j]?"true":"false");
    }
    printf("|\n");
    
    printf("new_bad_states:\n");
    for(int j = 0; j < ma->n; j++){
        printf("|%s",new_bad_states[j]?"true":"false");
    }
    printf("|\n");*/
}

/**
 * MEC Decomposition as described by Chatterjee and Henzinger as the ¨Old algorithm¨:
 * 
 * while(vertices left){
 * compute SCC decomposition
 *      foreach(SCC C){
 *              if(There are vertices (X) in C with random edges leaving SCC){
 *                      Remove Attr_R(X) U C from the graph
 *              }else{
 *                      Store MEC C, remove Attr_R(C) from the graph
 *              }
 *      }
 * }
 * @param ma
 * @return 
 */
SparseMatrixMEC* mEC_decomposition_previous_algorithm(SparseMatrix *ma){
    //printf("sccs2.ccp: mEC_decomposition_old_algorithm Start!\n");
    
    /*
     double dtstart, dtstop, ttime;
     dtstart = (double)clock()/CLOCKS_PER_SEC;
     */
	 
	 vector<unsigned long> mec_states(ma->n,0);
    
    /* Initialize boolean array bad, signifying which vertices are left to look at */
    
    bool *bad_states=(bool *) malloc(ma->n * sizeof(bool));
    bool *new_bad_states=(bool *) malloc(ma->n * sizeof(bool));
    bool *violating=(bool *) malloc(ma->n * sizeof(bool));
    
    for(unsigned long i=0; i<ma->n; i++) {
        bad_states[i]=false;
        new_bad_states[i]=false;
    }
    /* vector to store bad transition information */
    
    unsigned long *dist_starts = (unsigned long *) ma->row_counts;
    unsigned long nr_transitions = dist_starts[ma->n];
	
    bool *bad_transitions=(bool *) malloc(nr_transitions * sizeof(bool));
    
    for(unsigned long i=0; i<nr_transitions; i++) {
            bad_transitions[i]=false;
    }
    /* vector to store the SCC Decomposition answer */
    vector<unsigned long> scc_states(ma->n,0); //vector to store found SCC
    unsigned long scc_nr = 1;
    unsigned long nr_states_removed = 0;
    
    /* vector to store found mecs */
    //vector<unsigned long> mec_states(ma->n,0);
    
    unsigned long mec_nr = 1; //amount of mecs stored
    
	//printf("lets start\n");
    
    //Start of actual algorithm:
    while(nr_states_removed < ma->n){
        //printf("start main while: nr_states_removed: %lu/%lu\n",nr_states_removed,ma->n);
        
        //reset SCC decomposition answer variables.
        for(unsigned long i = 0; i < ma->n; i++){
            scc_states[i] = 0;
        }
        scc_nr = 1;
        //printf("start SCCD\n");
        compute_SCC_decomposition_tarjan(ma,scc_states,bad_states,bad_transitions,scc_nr);
        //printf("stop SCCD\n");
        
        for(unsigned long z = 1; z < scc_nr; z++){
            bool bad = false; //in order to know if there was a vertex that had a questionable bad transition
            unsigned long nr_violating = 0; //in order to know if violating set is empty
            
            //printf("start MECcheck\n");
            check_mec(ma,bad_states,bad_transitions,violating,scc_states,z,bad,nr_violating);
            //printf("stop MECcheck, nr_violating = %lu, bad = %s \n",nr_violating,(bad)?"true":"false");
            if((nr_violating == 0) && !bad){ //X is empty and bad == false
                //Each vertex in SCC only has good transitions; SCC is an MEC:
                //Store MEC C, 
                /*printf("scc array for scc: %lu\n",i);
                for(int j = 0; j < ma->n; j++){
                    printf("%lu",scc_states[j]);
                }*/
                
                for(int j = 0; j < ma->n; j++){
                    if(scc_states[j]==z){
                        mec_states[j] = mec_nr;
                    }
                }
                mec_nr++;
                
                //remove Attr_R(C) from the graph
                //printf("start attractor remove SCC\n");
                attractor_remove(ma,violating,bad_states,new_bad_states,bad_transitions,scc_states,z,nr_states_removed,true);
                //printf("stop attractor remove SCC\n");
                
                /*/temp
                for(unsigned long i = 0; i < ma->n; i++){
                        bad_states[i] = new_bad_states[i];
                }*/
            }else if(nr_violating > 0){ //X != empty
                //the SCC found called C is not an MEC due to violating vertex set X
                //Remove Attr_R(X) U (set union) C from the graph (add to bad states and bad vector).
                //printf("start attractor remove set union scc\n");
                attractor_remove(ma,violating,bad_states,new_bad_states,bad_transitions,scc_states,z,nr_states_removed,false);
                //printf("stop attractor remove set union scc\n");
                /*/temp
                for(unsigned long i = 0; i < ma->n; i++){
                        bad_states[i] = new_bad_states[i];
                }*/
            }                  
        }
        //iteration complete, apply newly found bad states, if done in the attractor function then check_mec function will not work correct
        //as it then has problems detecting violating vertexes (those that only have bad transitions), because the destination of such a bad transition
        //possibly is already marked as bad, which should only be the case in the next iteration.
        for(unsigned long i = 0; i < ma->n; i++){
            bad_states[i] = new_bad_states[i];
        }
        
        //Also remove all vertices that do not belong to any SCC, as the attractor now works correctly this should not remove any nodes.
        unsigned long counter = 0;
        for(unsigned long i = 0; i < ma->n; i++){
            if(scc_states[i] == 0 && !bad_states[i]){
                bad_states[i] = true;
                nr_states_removed++;
                counter++;
            }
        }
        //printf("Number of vertices removed that do not belong to any SCC %lu\n",counter);
        
        //printf("stop main while: nr_states_removed: %lu/%lu\n",nr_states_removed,ma->n);
    }
    
    unsigned long mec_states_nr = 0;
    for(unsigned long i = 0; i < ma->n; i++){
        if(mec_states[i]!=0){
            mec_states_nr++;
        }
    }
    //printf("Storing MECs, contents of mec_states:\n");
    /*
    for(unsigned long i = 0; i < ma->n; i++){
        printf("%lu|",mec_states[i]);
    }*/
    //printf(" -- contents of mec_nr: %lu, mec_states_nr: %lu\n",mec_nr,mec_states_nr);
    
    /* Store MECs in a new MA,  this is copied from sccs.cpp*/
    /* allocate memory for BSCCs and store them */
    SparseMatrixMEC *mec;
	
	//printf("creating MEC");
    
    mec=SparseMatrixMEC_new(mec_states_nr,mec_nr-1);
    unsigned long * cols = mec->cols;
    unsigned long *row_starts = (unsigned long *) mec->row_counts;

    vector<unsigned long>::const_iterator it;
    unsigned long pos=0;
    unsigned long scc_index=0;
    unsigned long scc_size=0;
    unsigned long tmp=0;
    for(int idx=0; idx<mec->n; idx++) {
            unsigned long i=0;
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
    
    //free allocated memory:
    free(bad_states);
    free(new_bad_states);
    free(violating);
    free(bad_transitions);
    /*
    dtstop = (double)clock()/CLOCKS_PER_SEC;
    ttime = dtstop-dtstart;
    */
    //printf("sccs2.ccp: mEC_decomposition_old_algorithm Done! Time: detailed: %f.\n",ttime);
    return mec;
}

SparseMatrixMEC* mEC_decomposition_previous_algorithm_without_attractor(SparseMatrix *ma,vector<unsigned long>& mec_states){
    //printf("sccs2.ccp: mEC_decomposition_old_algorithm Start!\n");
    
    /*
     double dtstart, dtstop, ttime;
     dtstart = (double)clock()/CLOCKS_PER_SEC;
     */
    
    /* Initialize boolean array bad, signifying which vertices are left to look at */
    
    bool *bad_states=(bool *) malloc(ma->n * sizeof(bool));
    bool *new_bad_states=(bool *) malloc(ma->n * sizeof(bool));
    bool *violating=(bool *) malloc(ma->n * sizeof(bool));
    
    
    for(unsigned long i=0; i<ma->n; i++) {
        bad_states[i]=false;
        new_bad_states[i]=false;
    }
    /* vector to store bad transition information */
    
    unsigned long *dist_starts = (unsigned long *) ma->row_counts;
    unsigned long nr_transitions = dist_starts[ma->n];
	
    bool *bad_transitions=(bool *) malloc(nr_transitions * sizeof(bool));
    
    for(unsigned long i=0; i<nr_transitions; i++) {
            bad_transitions[i]=false;
    }
    /* vector to store the SCC Decomposition answer */
    vector<unsigned long> scc_states(ma->n,0); //vector to store found SCC
    unsigned long scc_nr = 1;
    unsigned long nr_states_removed = 0;
    
    /* vector to store found mecs */
    //vector<unsigned long> mec_states(ma->n,0);
    
    unsigned long mec_nr = 1; //amount of mecs stored
    
    //Start of actual algorithm:
    while(nr_states_removed < ma->n){
        //printf("start main while: nr_states_removed: %lu/%lu\n",nr_states_removed,ma->n);
        
        //reset SCC decomposition answer variables.
        for(unsigned long i = 0; i < ma->n; i++){
            scc_states[i] = 0;
        }
        scc_nr = 1;
        //printf("start SCCD\n");
        compute_SCC_decomposition_tarjan(ma,scc_states,bad_states,bad_transitions,scc_nr);
        //printf("stop SCCD\n");
        
        for(unsigned long z = 1; z < scc_nr; z++){
            bool bad = false; //in order to know if there was a vertex that had a questionable bad transition
            unsigned long nr_violating = 0; //in order to know if violating set is empty
            
            //printf("start MECcheck\n");
            check_mec(ma,bad_states,bad_transitions,violating,scc_states,z,bad,nr_violating);
            /*
            printf("scc array for scc: %lu\n",z);
                for(int j = 0; j < ma->n; j++){
                    printf("%lu",scc_states[j]);
                }
            */
            //printf("stop MECcheck, nr_violating = %lu, bad = %s \n",nr_violating,(bad)?"true":"false");
            if((nr_violating == 0) && !bad){ //X is empty and bad == false
                //Each vertex in SCC only has good transitions; SCC is an MEC:
                //Store MEC C, 
                
                
                for(int j = 0; j < ma->n; j++){
                    if(scc_states[j]==z){
                        mec_states[j] = mec_nr;
                    }
                }
                mec_nr++;
                
                //remove Attr_R(C) from the graph
                
                for(unsigned long i = 0; i < ma->n; i++){
                    if((!new_bad_states[i]) &&(scc_states[i]==z)){
                        new_bad_states[i] = true;
                        nr_states_removed++;
                    }
                }
            }else if(nr_violating > 0){ //X != empty
                //the SCC found called C is not an MEC due to violating vertex set X
                //Remove Attr_R(X) U (set union) C from the graph (add to bad states and bad vector).
             
                for(unsigned long i = 0; i < ma->n; i++){
                    if(!new_bad_states[i]){
                        new_bad_states[i] = violating[i];
                    }
                }
            }                  
        }
        //iteration complete, apply newly found bad states, if done in the attractor function then check_mec function will not work correct
        //as it then has problems detecting violating vertexes (those that only have bad transitions), because the destination of such a bad transition
        //possibly is already marked as bad, which should only be the case in the next iteration.
        for(unsigned long i = 0; i < ma->n; i++){
            bad_states[i] = new_bad_states[i];
        }
        
        //Also remove all vertices that do not belong to any SCC:
        unsigned long counter = 0;
        for(unsigned long i = 0; i < ma->n; i++){
            if(scc_states[i] == 0 && !bad_states[i]){
                bad_states[i] = true;
                nr_states_removed++;
                counter++;
            }
        }
        //printf("Number of vertices removed that do not belong to any SCC %lu\n",counter);
        
        //printf("stop main while: nr_states_removed: %lu/%lu\n",nr_states_removed,ma->n);
    }
    
    unsigned long mec_states_nr = 0;
    for(unsigned long i = 0; i < ma->n; i++){
        if(mec_states[i]!=0){
            mec_states_nr++;
        }
    }
    //printf("Storing MECs, contents of mec_states:\n");
    /*
    for(unsigned long i = 0; i < ma->n; i++){
        printf("%lu|",mec_states[i]);
    }*/
    //printf(" -- contents of mec_nr: %lu, mec_states_nr: %lu\n",mec_nr,mec_states_nr);
    
    /* Store MECs in a new MA,  this is copied from sccs.cpp*/
    /* allocate memory for BSCCs and store them */
    SparseMatrixMEC *mec;
    
    mec=SparseMatrixMEC_new(mec_states_nr,mec_nr-1);
    unsigned long * cols = mec->cols;
    unsigned long *row_starts = (unsigned long *) mec->row_counts;

    vector<unsigned long>::const_iterator it;
    unsigned long pos=0;
    unsigned long scc_index=0;
    unsigned long scc_size=0;
    unsigned long tmp=0;
    for(int idx=0; idx<mec->n; idx++) {
            unsigned long i=0;
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
    
    //free allocated memory:
    free(bad_states);
    free(new_bad_states);
    free(violating);
    free(bad_transitions);
    /*
    dtstop = (double)clock()/CLOCKS_PER_SEC;
    ttime = dtstop-dtstart;
    */
    //printf("sccs2.ccp: mEC_decomposition_old_algorithm Done! Time: detailed: %f.\n",ttime);
    return mec;
}
