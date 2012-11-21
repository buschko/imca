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
* @file bcg2imca.cpp
* @brief Transform CADP .bcg file to IMCA .ma format
* @author Dennis Guck
* @version 1.0
*
* Original code from http://www.inrialpes.fr/vasy/cadp/man/bcg_read.html, from Hubert Garavel. 
* Changed version of Silvio De Carolis for his Masterthesis "Zuverlaessigkeitsanalyse auf dynamischen Fehlerbaeumen" at the chair of computer science 2 at RWTH-Aachen.
* Original comments in english, Silvio De Carolis comments in german.
* Extended by Dennis Guck for use in dftcalc -> CADP -> IMCA tool-chain.
* 
*/

#ifndef __BCG_TO_IMCA
#define __BCG_TO_IMCA
#ifdef __cplusplus
	extern "C" {
#endif 
#include "stdio.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bcg_user.h"



void bcg2imca (char* fname);


#ifdef __cplusplus
	}
#endif 

#endif // __BCG_TO_IMCA
