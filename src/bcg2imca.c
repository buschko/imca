/*
 * Original Example from http://www.inrialpes.fr/vasy/cadp/man/bcg_read.html, from Hubert Garavel. 
 * Changed version of Silvio De Carolis for his Masterthesis "Zuverlaessigkeitsanalyse auf dynamischen Fehlerbaeumen" at the chair of computer science 2 RWTH-Aachen
 * Original comments in english, Silvio De Carolis comments in german.
 * Extended by Dennis Guck for use in dftcalc -> CADP -> IMCA Tool-chain.auto
 */

#include "bcg_user.h"
#include "bcg_transition.h"
#include <stdio.h>
#include <string.h>
#include "bcg2imca.h"

/* The following function prints information about a BCG graph */
void bcg_print_info (bcg_graph) 
BCG_TYPE_OBJECT_TRANSITION bcg_graph;
{
	printf ("initial state = %lu\n", BCG_OT_INITIAL_STATE (bcg_graph));
	printf ("nb states = %lu\n", BCG_OT_NB_STATES (bcg_graph));
	printf ("nb edges = %lu\n", BCG_OT_NB_EDGES (bcg_graph));
	printf ("nb labels = %u\n", BCG_OT_NB_LABELS (bcg_graph));
}

/* The following function displays an edge 
 * and safes it in the .ma file
 */
void bcg_print_edge (bcg_graph, bcg_state_1, bcg_label_number, bcg_state_2, isMS, ma)
BCG_TYPE_OBJECT_TRANSITION bcg_graph;
BCG_TYPE_STATE_NUMBER bcg_state_1;
BCG_TYPE_LABEL_NUMBER bcg_label_number;
BCG_TYPE_STATE_NUMBER bcg_state_2;
int *isMS;
FILE* ma;
{
	BCG_TYPE_C_STRING bcg_label_string;
	BCG_TYPE_BOOLEAN bcg_visible;
	BCG_TYPE_C_STRING bcg_gate;
	
	bcg_label_string = BCG_OT_LABEL_STRING (bcg_graph, bcg_label_number);
	bcg_visible = BCG_OT_LABEL_VISIBLE (bcg_graph, bcg_label_number);
	
	// print the edge to output
	printf ("transition from state %lu to state %lu\n", bcg_state_1, bcg_state_2);
	printf ("label unique number = %u\n", bcg_label_number);
	printf ("label string = %s\n", bcg_label_string);
	if (bcg_visible) {
		bcg_gate = BCG_OT_LABEL_GATE (bcg_graph, bcg_label_number);
		printf ("visible label (gate = %s)\n", bcg_gate);
	} else {
		bcg_gate = BCG_OT_LABEL_HIDDEN_GATE (bcg_graph, bcg_label_number);
		printf ("hidden label (hidden gate = %s)\n", bcg_gate);
	}
	
	//ignore selfloops (due to DFT condition)
	if(bcg_state_1!=bcg_state_2) {
		// add edge to .ma file
		if(*isMS==1 && strncmp(bcg_label_string,"rate",4)==0){
			strcpy(bcg_label_string+0, bcg_label_string+5); // obtain the rate from string
			fprintf(ma,"* s%lu %s\n",bcg_state_2,bcg_label_string);
		} else if(strncmp(bcg_label_string,"rate",4)==0) {
			*isMS = 1;
			fprintf(ma,"s%lu !\n",bcg_state_1);
			strcpy(bcg_label_string+0, bcg_label_string+5); // obtain the rate from string
			fprintf(ma,"* s%lu %s\n",bcg_state_2,bcg_label_string);
		} else {
			*isMS = 0;
			fprintf(ma,"s%lu %s\n",bcg_state_1,bcg_label_string);
			fprintf(ma,"* s%lu 1\n",bcg_state_2);
		}
	}
}

int main(int argc, char* argv[]) {
	//bcg2imca(argv[1]);
	argc=argc;
	
	BCG_TYPE_OBJECT_TRANSITION bcg_graph;
	BCG_TYPE_STATE_NUMBER bcg_s1;
	BCG_TYPE_LABEL_NUMBER bcg_label_number;
	BCG_TYPE_STATE_NUMBER bcg_s2;
	BCG_TYPE_STATE_NUMBER bcg_nb_states;

	BCG_INIT ();
	
	/* create new File for .ma */
	FILE* ma; // Zeiger auf Datenstrom der Datei
	ma = fopen(argv[2],"w+"); // Datei neu erzeugen bzw. ueberschreiben, wenn es sie schon gibt

   /* The following fragment of code reads and prints all the edges 
      of a BCG graph, in an undefined order */
	BCG_OT_READ_BCG_BEGIN (argv[1], &bcg_graph, 0);
	bcg_print_info (bcg_graph);
	
	/* add initial and goal states in .ma Header */
	fprintf(ma, "#INITIALS\ns%lu\n",BCG_OT_INITIAL_STATE (bcg_graph)); 
	fprintf(ma, "#GOALS\n");
	
	BCG_TYPE_C_STRING bcg_gate;
	
	BCG_OT_ITERATE_PLN (bcg_graph, bcg_s1, bcg_label_number, bcg_s2) {
		bcg_gate = BCG_OT_LABEL_GATE (bcg_graph, bcg_label_number);
		if(strcmp(bcg_gate,argv[3]) == 0)
			fprintf(ma, "s%lu\n",bcg_s2);
	} BCG_OT_END_ITERATE;
	BCG_OT_READ_BCG_END (&bcg_graph);
	
	fprintf(ma,"#TRANSITIONS\n");
	
	/* The following fragment of code reads and prints all the edges 
	 * of a BCG graph, sorted by origin states in increasing order */
	
	int isMS=0; // helper to identify Markovian transitions
	
	BCG_OT_READ_BCG_BEGIN (argv[1], &bcg_graph, 1);
	bcg_nb_states = BCG_OT_NB_STATES (bcg_graph);
	for (bcg_s1 = 0; bcg_s1 < bcg_nb_states; bcg_s1++) {
		printf ("successors of state %lu:\n", bcg_s1);
		BCG_OT_ITERATE_P_LN (bcg_graph, bcg_s1, bcg_label_number, bcg_s2)
		{
			bcg_print_edge (bcg_graph, bcg_s1, bcg_label_number, bcg_s2, &isMS, ma);
		} BCG_OT_END_ITERATE;
	}
	BCG_OT_READ_BCG_END (&bcg_graph);
	
	fclose(ma);
	exit (0);
}
