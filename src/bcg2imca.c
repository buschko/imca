/*
 * Original Example from http://www.inrialpes.fr/vasy/cadp/man/bcg_read.html, from Hubert Garavel. 
 * Changed version of Silvio De Carolis for his Masterthesis "Zuverlaessigkeitsanalyse auf dynamischen Fehlerbaeumen" at the chair of computer science 2 RWTH-Aachen
 * Original comments in english, Silvio De Carolis comments in german.
 * Extended by Dennis Guck for use in dftcalc -> CADP -> IMCA Tool-chain
 * Extended by Axel Belinfante to collect markovian and non-markovian edges seprately,
 *     using dstring, based on make_message example from Linux man page printf(3) 
 */

#include "bcg_user.h"
#include "bcg_transition.h"
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include "bcg2imca.h"

/* ugly, but only enabled in string.h if we set some special macro */
extern char* strdup(const char*);


typedef struct dstring {
	char *s;         /* string buffer */
	char *p;         /* first unused position in string buffer */
	char* beyond;    /* first position beyond allocated string */
} dstring;

void
fatal(char *s) {
	fprintf(stderr, "fatal: %s\n", s);
	exit(1);
}

void
dstring_status(dstring *ds, char * reason) {
	int size = ds->beyond - ds->s;
	int avail = ds->beyond - ds->p;
	int used = ds->p - ds->s;
	fprintf(stderr, "%s dstring %p: size: %d, used:%d, avail:%d,  s: %p, p: %p, beyond: %p\n", reason, ds, size, used, avail, ds->s, ds->p, ds->beyond);
}

dstring*
dstring_new() {
	int size = 3; /* small initial size, to force resizes, for testing */
	dstring *ds = malloc(sizeof(struct dstring));
	if (ds==0)
		fatal("out of memory");
	ds->s = malloc(size);
	if (ds->s==0)
		fatal("out of memory");
        memset(ds->s, 0, size);
	ds->p = ds->s;
	ds->beyond = ds->s + size;
	dstring_status(ds, "new");
	return ds;
}

void
dstring_reset(dstring *ds) {
	int size = ds->beyond - ds->s;
	ds->p = ds->s;
        memset(ds->s, 0, size);
	dstring_status(ds, "reset");
}

char *
dstring_printf(dstring *ds, const char *fmt, ...)
{
	int n, size, used, avail, inc, newsize;
	va_list ap;

	size = ds->beyond - ds->s;
	avail = ds->beyond - ds->p;
	used = ds->p - ds->s;
	dstring_status(ds, "printf");
	//fprintf(stderr, "dstring_printf dstring %p: size: %d, used:%d, avail:%d,  s: %p, p: %p, beyond: %p\n", ds, size, used, avail, ds->s, ds->p, ds->beyond);

	while (1) {
		/* Try to print in the allocated space. */
		va_start(ap, fmt);
		n = vsnprintf(ds->p, avail, fmt, ap);
		va_end(ap);
		/* If that worked, return the string. */
		if (n > -1 && n < avail) {
			ds->p += n;
			size = ds->beyond - ds->s;
			avail = ds->beyond - ds->p;
			used = ds->p - ds->s;
			fprintf(stderr, "printed dstring %p: size:%d, used:%d, avail:%d, s: %p, p: %p, beyond: %p\n", ds, size, used, avail, ds->s, ds->p, ds->beyond);
			return  ds->p - n;
		}
		/* Else try again with more space. */
		if (n > -1)    /* glibc 2.1 */
			inc = n+1; /* precisely what is needed */
		else           /* glibc 2.0 */
			inc = size * 2;  /* twice the old size */
		newsize = size + inc;
		fprintf(stderr, "resize dstring %p: size: %d, newsize: %d\n", ds, size, newsize);
		if ((ds->s = realloc (ds->s, newsize)) == NULL) {
			fatal("out of memory");
		} else {
			ds->p = ds->s + used;
			ds->beyond = ds->s + newsize;
			size = ds->beyond - ds->s;
			avail = ds->beyond - ds->p;
			used = ds->p - ds->s;
			fprintf(stderr, "resized dstring %p: size:%d, used:%d, avail:%d, s: %p, p: %p, beyond: %p\n", ds, size, used, avail, ds->s, ds->p, ds->beyond);
		}
	}
	return NULL;
}


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
 * and safes it to either dstring markovian, or to dstring nonmarkovian
 */
void bcg_print_edge (bcg_graph, bcg_state_1, bcg_label_number, bcg_state_2, markovian, nonmarkovian)
BCG_TYPE_OBJECT_TRANSITION bcg_graph;
BCG_TYPE_STATE_NUMBER bcg_state_1;
BCG_TYPE_LABEL_NUMBER bcg_label_number;
BCG_TYPE_STATE_NUMBER bcg_state_2;
dstring *markovian;
dstring *nonmarkovian;
{
	BCG_TYPE_C_STRING bcg_label_string;
	BCG_TYPE_BOOLEAN bcg_visible;
	BCG_TYPE_C_STRING bcg_gate;

	bcg_label_string = strdup(BCG_OT_LABEL_STRING (bcg_graph, bcg_label_number));
	bcg_label_gate = strdup(BCG_OT_LABEL_GATE (bcg_graph, bcg_label_number));
	bcg_visible = BCG_OT_LABEL_VISIBLE (bcg_graph, bcg_label_number);
	
	// print the edge to output
	printf ("\ttransition from state %lu to state %lu\n", bcg_state_1, bcg_state_2);
	printf ("\t\tlabel unique number = %u\n", bcg_label_number);
	printf ("\t\tlabel string = %s\n", bcg_label_string);
	if (bcg_visible) {
		bcg_gate = strdup(BCG_OT_LABEL_GATE (bcg_graph, bcg_label_number));
		printf ("\t\tvisible label (gate = %s)\n", bcg_gate);
	} else {
		bcg_gate = strdup(BCG_OT_LABEL_HIDDEN_GATE (bcg_graph, bcg_label_number));
		printf ("\t\thidden label (hidden gate = %s)\n", bcg_gate);
	}

	//if(strcmp(bcg_label_string, "rate 0.2")==0)
	//	bcg_label_string = "act3";
	//printf ("\t\tlabel string now = %s\n", bcg_label_string);
	
	//ignore selfloops (due to DFT condition)
	if(bcg_state_1!=bcg_state_2) {
		// add edge
		if(strncmp(bcg_label_string,"rate",4)==0){
			strcpy(bcg_label_string+0, bcg_label_string+5); // obtain the rate from string
			dstring_printf(markovian,"* s%lu %s\n",bcg_state_2,bcg_label_string);
		} else {
			dstring_printf(nonmarkovian,"s%lu %s\n",bcg_state_1,bcg_label_string);
			dstring_printf(nonmarkovian,"* s%lu 1\n",bcg_state_2);
		}
	}
	free(bcg_label_string);
	free(bcg_gate);
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
	
	dstring *markovian = dstring_new();
	dstring *nonmarkovian = dstring_new();
	
	BCG_OT_READ_BCG_BEGIN (argv[1], &bcg_graph, 1);
	bcg_nb_states = BCG_OT_NB_STATES (bcg_graph);
	for (bcg_s1 = 0; bcg_s1 < bcg_nb_states; bcg_s1++) {
		dstring_reset(markovian);
		dstring_reset(nonmarkovian);
		printf ("successors of state %lu:\n", bcg_s1);
		BCG_OT_ITERATE_P_LN (bcg_graph, bcg_s1, bcg_label_number, bcg_s2)
		{
			bcg_print_edge (bcg_graph, bcg_s1, bcg_label_number, bcg_s2, markovian, nonmarkovian);
		} BCG_OT_END_ITERATE;
		if (strcmp(markovian->s, "") != 0) {
			fprintf(ma,"s%lu !\n",bcg_s1);
			fprintf(ma,"%s", markovian->s);
		}
		fprintf(ma,"%s", nonmarkovian->s);
	}
	BCG_OT_READ_BCG_END (&bcg_graph);
	
	fclose(ma);
	exit (0);
}
