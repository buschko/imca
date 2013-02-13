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
* 	Interactive Markov Chain Analyzer
*
* Created by Dennis Guck
* dennis.guck@rwth.aachen.de
*
*/

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <map>
#include <string>
#include "soplex.h"
#include <popt.h>

#ifndef __APPLE__
#include <malloc.h>
#endif

#include "read_file.h"
#include "read_file_imc.h"
#include "unbounded.h"
#include "expected_time.h"
#include "sccs.h"
#include "sccs2.h"
#include "long_run_average.h"
#include "bounded.h"

#ifndef __APPLE__
#include <time.h>
struct timespec tp;
double begin, end;
#endif

/* List of possible models */
#define MA_MODE_STR "ma"
#define IMC_MODE_STR "imc"

/* List of possible file extensions */
#define MA_FILE_EXT ".ma"
#define IMC_FILE_EXT ".imc"

/* This is the list of possible options */
#define UNBOUND_STR "-ub"
#define EXP_TIME_STR "-et"
#define LONG_RUN_AVERAGE_STR "-lra"
#define TIME_BOUNDED_STR "-tb"
#define MIN_MODE_STR "-min"
#define MAX_MODE_STR "-max"
//#define FROM_STR     "--from"
#define TO_STR       "--to"
#define IMC_STR "-imc"

// Boolean macros for option checking
#define IS_LOWER_BOUND(str) (strncmp(str,"--from", 7) == 0 || strncmp(str, "-F", 3 ) == 0)
#define IS_UPPER_BOUND(str) (strncmp(str,"--to", 5) == 0 || strncmp(str, "-T", 3 ) == 0)
#define IS_ERROR_BOUND(str) (strncmp(str,"--error-bound", 14) == 0 || strncmp(str, "-e", 3 ) == 0)
/*
 * The min file name length including extensions:
 *	.ma
 * 	.imc
 */
#define MIN_FILE_NAME_LENGTH 5

/* The min and max length of the extensions:
 *	.ma
 * 	.imc
 */
#define MAX_FILE_EXT_LENGTH 5
#define MIN_FILE_EXT_LENGTH 3

/**
* Boolean flags whether certain files were loaded or not
*/
static bool is_ma_present = false;
static bool is_imc_present = false;
static bool is_max_present = false;
static bool is_min_present = false;
static bool is_unbound_present = false;
static bool is_expected_time_present = false;
static bool is_time_bounded_present = false;
static bool is_lra_present = false;
static bool is_imc = false;

static bool is_lower_bound_present = false;
static bool is_upper_bound_present = false;
static bool is_error_bound_present = false;

/**
* Global variables
*/
static SparseMatrix *ma = NULL;		/* The statespace of an MA and its transitions */
static const char * ma_file  = NULL;	/* pointer to the input file */
static Real epsilon = 1e-6;				/* default error bound for time-bounded reachability */
static Real ta = 0;						/* default lower bound for time interval ( time-bounded reachability ) */
static Real tb = 0;						/* default upper bound for time interval ( time-bounded reachability ) */

using namespace std;

/**
* Programm Intro message
*/
static void print_intro(void) {
	printf(" ------------------------------------------------------------------- \n");
	printf("|                  Interactive Markov Chain Analyzer                |\n");
	printf("|                     IMCA Version 1.4 no segfault                  |\n");
	printf("|                                                                   |\n");
	printf("|                   Copyright (C) RWTH Aachen, 2012                 |\n");
	printf("|                         Author: Dennis Guck                       |\n");
	printf("|           IMCA is distributed under the GPL conditions            |\n");
	printf("|            (GPL stands for GNU General Public License)            |\n");    
	printf("|          This program comes with ABSOLUTELY NO WARRANTY.          |\n");
	printf("|   This is free software, and you are welcome to redistribute it   |\n");
	printf(" ------------------------------------------------------------------- \n");
	printf("\n");
}

/**
* Program usage
*/
static void print_usage(void) {
	printf("Usage: imca <model file> <min/max> <computation>\n");
	printf("	<model file>	- could be one of {.ma, .imc}\n");
	printf("	<min/max>	- could be '-min' or '-max' or both\n");
	printf("	<computation>	- could be one or more of {-ub, -et, -lra, -tb}\n");
	printf("	<options>	- time- and error-bound for tb (default epsilon=1e-6) \n");
	printf("	<model type>	- define if .ma input is an IMC {-imc} \n");
}

/**
* This function should find the extension and then return it and its length.
* @param filename the input parameter which we suspect to be a file name
* @param extension the pointer to an extension string (a return parameter)
* @param ext_length the pointer to the extension string length (a return parameter)
* @return true if we can subtract a file extension which we think is one of required
*/
static bool isValidExtension(const char * filename, char * extension, int * ext_length) {
	bool result = false;
	const char *p;
	int length;

	length = strlen(filename);
	p = strrchr(filename,'.');     /* The last occurance of '.' in the file name */
	*ext_length = length - (p - filename); /* including '.', excluding '\0' */
	/*TODO: change function */
	if( length >= MIN_FILE_NAME_LENGTH && ( *ext_length >= MIN_FILE_EXT_LENGTH) && ( *ext_length <= MAX_FILE_EXT_LENGTH)) {
		/* Get the extension */
		strncpy(extension, p,  max(MIN_FILE_EXT_LENGTH, min(MAX_FILE_EXT_LENGTH, length)));
		result = true;
	}
	

	return result;
}

/**
* This function checks if all required parameters are set
*/
static void checkComputation(){
	if(!is_ma_present && !is_imc_present) {
		printf("ERROR: No model type was set.\n");
		print_usage();
		exit(EXIT_FAILURE);
      
	}
	if(!is_expected_time_present && !is_unbound_present && !is_lra_present && !is_time_bounded_present){
		printf("ERROR: No computation type was set.\n");
		print_usage();
		exit(EXIT_FAILURE);
	}
	if( !is_max_present && !is_min_present ){
		printf("ERROR: No settings for maximum\\minimum computation.\n");
		print_usage();
		exit(EXIT_FAILURE);
	}
	if( is_time_bounded_present && !is_lower_bound_present && !is_upper_bound_present) {
		printf("ERROR: To compute time bounded reachability, for time interval an upper bound (with '--to' or '-T') and optionally a lower bound (with '--from' or '-F' and default value of 0) are required.\n");
		print_usage();
		exit(EXIT_FAILURE);
	}
	if( is_time_bounded_present && !is_lower_bound_present && is_upper_bound_present) {
		printf("WARNING: No lower bound for time interval specified. The default value is %f.\n",ta);
	}
	if( is_time_bounded_present && !is_error_bound_present ) {
		printf("WARNING: No error bound specified. The default value is %f.\n", epsilon);
	}
}

/**
* This function will check the input parameters and validate them
*/
static void parseParams(int argc, char *argv[]) {
	int i = 0, ext_length;
	char extension[MAX_FILE_EXT_LENGTH];

	/* First parse the input parameters */
	for(i = 1; i < argc; i++)
	{
		if( isValidExtension( argv[i], extension, &ext_length ) ){
			if( strcmp(extension, MA_FILE_EXT) == 0 ){
				if( !is_ma_present ){
					is_ma_present = true;
					ma_file = argv[i];
				}else{
					printf("WARNING: A model has already noticed before, skipping the '%s' file.\n", argv[i]);
				}
			} else if( strcmp(extension, IMC_FILE_EXT) == 0 ){
				if( !is_ma_present ){
					is_imc_present = true;
					is_imc = true;
					ma_file = argv[i];
				}else{
					printf("WARNING: A model has already noticed before, skipping the '%s' file.\n", argv[i]);
				}
			} 
		}  else if( strcmp(argv[i], UNBOUND_STR) == 0 ){
			if( !is_unbound_present ){
				is_unbound_present = true;
			}else{
				printf("WARNING: The option has been noticed before, skipping '%s'.\n", argv[i]);
			}
		} else if( strcmp(argv[i], EXP_TIME_STR) == 0 ){
			if( !is_expected_time_present ){
				is_expected_time_present = true;
			}else{
				printf("WARNING: The option has been noticed before, skipping '%s'.\n", argv[i]);
			}
		} else if( strcmp(argv[i], LONG_RUN_AVERAGE_STR) == 0 ){
			if( !is_lra_present ){
				is_lra_present = true;
			}else{
				printf("WARNING: The option has been noticed before, skipping '%s'.\n", argv[i]);
			}
		} else if( strcmp(argv[i], TIME_BOUNDED_STR) == 0 ){
			if( !is_time_bounded_present ){
				is_time_bounded_present = true;
			}else{
				printf("WARNING: The option has been noticed before, skipping '%s'.\n", argv[i]);
			}
		}else if( IS_LOWER_BOUND(argv[i]) ){
			if( is_lower_bound_present ) {
				printf("ERROR: '%s' is repeated.\n", argv[i]);
				exit(EXIT_FAILURE);
			} else if(i+1 >= argc ) {
				printf("ERROR: No lower bound for time interval specified.\n");
				exit(EXIT_FAILURE);
			} else {
				char *toEnd;
				ta=strtod(argv[i+1], &toEnd);
				if( *toEnd == '\0' && ta >= 0 ) {
					is_lower_bound_present = true;
					i++;
				}
				else {
					printf("ERROR: The specified lower bound ('%s') is invalid. After '%s' must be a real value greater than or equal to zero.\n", argv[i+1], argv[i]);
					exit(EXIT_FAILURE);
				}
			}
		}else if( IS_UPPER_BOUND(argv[i]) ){
			if( is_upper_bound_present ) {
				printf("ERROR: '%s' is repeated.\n", argv[i]);
				exit(EXIT_FAILURE);
			} else if(i+1 >= argc ) {
				printf("ERROR: No upper bound for time interval specified.\n");
				exit(EXIT_FAILURE);
			} else {
				char *toEnd;
				tb=strtod(argv[i+1], &toEnd);
				if( *toEnd == '\0' && tb > 0) { //TODO adding support of reachability computation for ta = tb = 0
					is_upper_bound_present = true;
					i++;
				}
				else {
					printf("ERROR: The specified upper bound ('%s') is invalid. After '%s' must be a real value greater than zero.\n", argv[i+1], argv[i]);
					exit(EXIT_FAILURE);
				}
			}
		}else if( IS_ERROR_BOUND(argv[i]) ){
			if( is_error_bound_present ) {
				printf("ERROR: '%s' is repeated.\n", argv[i]);
				exit(EXIT_FAILURE);
			} else if(i+1 >= argc ) {
				printf("ERROR: No error bound specified.\n");
				exit(EXIT_FAILURE);
			} else {
				char *toEnd;
				epsilon = strtod(argv[i+1], &toEnd);
				if( *toEnd == '\0' && epsilon > 0 ) {
					is_error_bound_present = true;
					i++;
				}
				else {
					printf("ERROR: The specified upper bound ('%s') is invalid. After '%s' must be a real value greater than zero.\n", argv[i+1], argv[i]);
					exit(EXIT_FAILURE);
				}
			}
		} else if( strcmp(argv[i], MIN_MODE_STR) == 0 ){
			if( !is_min_present ){
				is_min_present = true;
			}else{
				printf("WARNING: The option has been noticed before, skipping '%s'.\n", argv[i]);
			}
		} else if( strcmp(argv[i], MAX_MODE_STR) == 0 ){
			if( !is_max_present ){
				is_max_present = true;
			}else{
				printf("WARNING: The option has been noticed before, skipping '%s'.\n", argv[i]);
			}
		}else if( strcmp(argv[i], IMC_STR) == 0 ){
			if( !is_imc ){
				is_imc = true;
			}else{
				printf("WARNING: The option has been noticed before, skipping '%s'.\n", argv[i]);
			}
		} 
	}
	
	checkComputation();
}

/**
* Load the .imc file
*/
static void loadIMC(const char *filename) {
	if(is_imc_present) {
		printf("Loading the '%s' file, please wait.\n", filename);
		ma = read_IMC_SparseMatrix_file(filename);
		if(ma == NULL ){
			printf("ERROR: The '%s' file '%s' was not found or is incorrect!\n",IMC_FILE_EXT, filename);
			exit(EXIT_FAILURE);
		}
	}
}

/**
* Load the .ma file
*/
static void loadMA(const char *filename) {
	if(is_ma_present) {
		printf("Loading the '%s' file, please wait.\n", filename);
		ma = read_MA_SparseMatrix_file(filename);
		if(ma == NULL ){
			printf("ERROR: The '%s' file '%s' was not found or is incorrect!\n",MA_FILE_EXT, filename);
			exit(EXIT_FAILURE);
		}
	}
}

int main(int argc, char* argv[]) {

	#ifndef __APPLE__
	// helpers to calculate allocated memory
	unsigned long int sp1=0,sp2=0;
	#endif
	
	#ifndef __APPLE__
	// get memory info before model is loaded
	sp1 = mallinfo().uordblks;
	#endif

	/// print the MAA intro
	print_intro();
	
	/// Parse and validate the input parameters
	parseParams(argc, argv);
	
	if(is_ma_present)
	{
		/// load the MA from file
		loadMA(ma_file);
	} else if(is_imc_present) {
		loadIMC(ma_file);
	}
	
	#ifndef __APPLE__
	// get memory info after model is loaded
	sp2 = mallinfo().uordblks;
	#endif

	#ifndef __APPLE__
	double size = (sp2-sp1);
	printf("The occupied space is %d Bytes", (int)size);
	if((size/1024) < 1024) {
		printf(" (~%f KB).\n\n", (size/1024));
	} else if((size/1024)/1024 < 1024){
		printf(" (~%f MB).\n\n", (size/1024)/1024);
	} else {
		printf(".\n\n");
	}
	#else
	printf("The occupied space is ??? Bytes.\n\n");
	#endif
	
	Real tmp;

	if(is_unbound_present){
		if(is_max_present){
			#ifndef __APPLE__
			clock_gettime(CLOCK_REALTIME, &tp);
			begin = 1e9*tp.tv_sec + tp.tv_nsec;
			#endif
			printf("\nCompute maximal unbounded reachability, please wait.\n");	
			tmp = compute_unbounded_reachability(ma,true);	
			printf("Maximal unbounded reachability: %.10g\n", tmp);
			#ifndef __APPLE__
			clock_gettime(CLOCK_REALTIME, &tp);
			end = 1e9*tp.tv_sec + tp.tv_nsec;
			printf("Computation Time: %f seconds\n", (end-begin)*1e-9);
			#else
			printf("Computation Time: ??? seconds\n");
			#endif
		}
		if(is_min_present){
			#ifndef __APPLE__
			clock_gettime(CLOCK_REALTIME, &tp);
			begin = 1e9*tp.tv_sec + tp.tv_nsec;
			#endif
			printf("\nCompute minimal unbounded reachability, please wait.\n");	
			tmp = compute_unbounded_reachability(ma,false);	
			printf("Minimal unbounded reachability: %.10g\n", tmp);
			#ifndef __APPLE__
			clock_gettime(CLOCK_REALTIME, &tp);
			end = 1e9*tp.tv_sec + tp.tv_nsec;
			printf("Computation Time: %f seconds\n", (end-begin)*1e-9);
			#else
			printf("Computation Time: ??? seconds\n");
			#endif
		}
	}
	if(is_expected_time_present){
		if(is_max_present){
			#ifndef __APPLE__
			clock_gettime(CLOCK_REALTIME, &tp);
			begin = 1e9*tp.tv_sec + tp.tv_nsec;
			#endif
			printf("\nCompute maximal expected time, please wait.\n");
			tmp = compute_expected_time(ma,true);	
			printf("Maximal expected time: %.10g\n", tmp);
			#ifndef __APPLE__
			clock_gettime(CLOCK_REALTIME, &tp);
			end = 1e9*tp.tv_sec + tp.tv_nsec;
			printf("Computation Time: %f seconds\n", (end-begin)*1e-9);
			#else
			printf("Computation Time: ??? seconds\n");
			#endif
		}
		if(is_min_present){
			#ifndef __APPLE__
			clock_gettime(CLOCK_REALTIME, &tp);
			begin = 1e9*tp.tv_sec + tp.tv_nsec;
			#endif
			printf("\nCompute minimal expected time, please wait.\n");	
			tmp = compute_expected_time(ma,false);	
			printf("Minimal expected time: %.10g\n\n", tmp);
			#ifndef __APPLE__
			clock_gettime(CLOCK_REALTIME, &tp);
			end = 1e9*tp.tv_sec + tp.tv_nsec;
			printf("Computation Time: %f seconds\n", (end-begin)*1e-9);
			#else
			printf("Computation Time: ??? seconds\n");
			#endif
		}
	}
	if(is_lra_present){
		if(is_max_present){
			#ifndef __APPLE__
			clock_gettime(CLOCK_REALTIME, &tp);
			begin = 1e9*tp.tv_sec + tp.tv_nsec;
			#endif
			printf("\nCompute maximal LRA, please wait.\n");
			tmp=compute_long_run_average(ma,true);
			printf("Maximal LRA: %.10g\n", tmp);
			#ifndef __APPLE__
			clock_gettime(CLOCK_REALTIME, &tp);
			end = 1e9*tp.tv_sec + tp.tv_nsec;
			printf("Computation Time: %f seconds\n", (end-begin)*1e-9);
			#else
			printf("Computation Time: ??? seconds\n");
			#endif
		}
		if(is_min_present){
			#ifndef __APPLE__
			clock_gettime(CLOCK_REALTIME, &tp);
			begin = 1e9*tp.tv_sec + tp.tv_nsec;
			#endif
			printf("\nCompute minimal LRA, please wait.\n");
			tmp=compute_long_run_average(ma,false);
			printf("Minimal LRA: %.10g\n", tmp);
			#ifndef __APPLE__
			clock_gettime(CLOCK_REALTIME, &tp);
			end = 1e9*tp.tv_sec + tp.tv_nsec;
			printf("Computation Time: %f seconds\n", (end-begin)*1e-9);
			#else
			printf("Computation Time: ??? seconds\n");
			#endif
		}
	}
	if(is_time_bounded_present){
		if(is_max_present){
			#ifndef __APPLE__
			clock_gettime(CLOCK_REALTIME, &tp);
			begin = 1e9*tp.tv_sec + tp.tv_nsec;
			#endif
			printf("\nCompute maximal time-bounded reachability inside interval [%g,%g] with precision %g, please wait.\n", ta, tb, epsilon);
			tmp=compute_time_bounded_reachability(ma,true,epsilon,ta,tb,is_imc);
			printf("Maximal time-bounded reachability probability: %.10g\n", tmp);
			#ifndef __APPLE__
			clock_gettime(CLOCK_REALTIME, &tp);
			end = 1e9*tp.tv_sec + tp.tv_nsec;
			printf("Computation Time: %f seconds\n", (end-begin)*1e-9);
			#else
			printf("Computation Time: ??? seconds\n");
			#endif
		}
		if(is_min_present){
			#ifndef __APPLE__
			clock_gettime(CLOCK_REALTIME, &tp);
			begin = 1e9*tp.tv_sec + tp.tv_nsec;
			#endif
			printf("\nCompute minimal time-bounded reachability inside interval [%g,%g] with precision %g, please wait.\n", ta, tb, epsilon);
			tmp=compute_time_bounded_reachability(ma,false,epsilon,ta,tb,is_imc);
			printf("Minimal time-bounded reachability probability: %.10g\n", tmp);
			#ifndef __APPLE__
			clock_gettime(CLOCK_REALTIME, &tp);
			end = 1e9*tp.tv_sec + tp.tv_nsec;
			printf("Computation Time: %f seconds\n", (end-begin)*1e-9);
			#else
			printf("Computation Time: ??? seconds\n");
			#endif
		}
	}
	
	SparseMatrix_free(ma);
		
	return 0;
}
