#-----------------------------------------------------------------------------
#
# IMCA is a analyzing tool for unbounded reachability probabilities, expected-
# time, and long-run averages for Interactive Markov Chains and Markov Automata.
# Copyright (C) RWTH Aachen, 2012
# 	Author: Dennis Guck
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#
# Source description:
# 	Makefile for the Interactive Markov Chain Analyzer (IMCA)
#
# Created by Dennis Guck
# dennis.guck@rwth.aachen.de
#
#-----------------------------------------------------------------------------

VERSION		:=	1.5.beta

#-----------------------------------------------------------------------------
# Architecture
#-----------------------------------------------------------------------------
ARCH		:=	$(shell uname -m | \
			sed \
			-e 's/sun../sparc/' \
			-e 's/i.86/x86/' \
			-e 's/i86pc/x86/' \
			-e 's/[0-9]86/x86/' \
			-e 's/amd64/x86_64/' \
			-e 's/IP../mips/' \
			-e 's/9000..../hppa/' \
			-e 's/Power\ Macintosh/ppc/' \
			-e 's/00........../pwr4/' )
OSTYPE		:=	$(shell uname -s | tr '[:upper:]' '[:lower:]' | \
			sed \
			-e 's/cygwin.*/cygwin/' \
			-e 's/irix../irix/' \
			-e 's/windows.*/windows/' \
			-e 's/mingw.*/mingw/')
HOSTNAME	:=	$(shell uname -n | tr '[:upper:]' '[:lower:]')

UNAME := $(shell uname)

#-----------------------------------------------------------------------------
# Base path for Libaries -- change this according to local install
#-----------------------------------------------------------------------------
LIBPATH		= /libs

#-----------------------------------------------------------------------------
# Soplex Libaries -- path to Soplex library has to be adjusted
#-----------------------------------------------------------------------------
SOPLEXPATH	= $(LIBPATH)/soplex-1.7.1
SOPLEXSRCDIR	= $(SOPLEXPATH)/src/
SOPLEXLIBDIR	= $(SOPLEXPATH)/lib/
SOPLEXLIB	= soplex z

#-----------------------------------------------------------------------------
# lp_solve Libaries -- path to lp_solve library has to be adjusted
#-----------------------------------------------------------------------------
LPSOLVEPATH	= $(LIBPATH)/lp_solve_5.5
LPSOLVESRCDIR	= $(LPSOLVEPATH)
LPSOLVELIBDIR	= $(LPSOLVEPATH)/Release
LPSOLVELIB	= lpsolve55 dl

#-----------------------------------------------------------------------------
# GLPK Libaries -- path to GLPK library has to be adjusted
#-----------------------------------------------------------------------------
GLPKPATH	= $(LIBPATH)/glpk-4.53
GLPKSRCDIR	= $(GLPKPATH)/src/
GLPKLIBDIR	= $(GLPKPATH)/src/.libs
GLPKLIB		= glpk

# set to true to enable lp solver support, otherwise false
USELP		= true
# 1 SOPLEX
# 2 LPSOLVE
# 3 GLPK
LPSOLVER	= 1

ifeq ($(LPSOLVER), 1)
LPPATH		= $(SOPLEXPATH)
LPINCLUDE	= $(SOPLEXSRCDIR)
LPLIBDIR	= $(SOPLEXLIBDIR)
LPLIB		= $(SOPLEXLIB)
else
ifeq ($(LPSOLVER), 2)
LPPATH		= $(LPSOLVEPATH)
LPINCLUDE	= $(LPSOLVESRCDIR)
LPLIBDIR	= $(LPSOLVELIBDIR)
LPLIB		= $(LPSOLVELIB)
else
ifeq ($(LPSOLVER), 3)
LPPATH		= $(SOPLEXPATH)
LPINCLUDE	= $(GLPKSRCDIR)
LPLIBDIR	= $(GLPKLIBDIR)
LPLIB		= $(GLPKLIB)
endif
endif
endif

#-----------------------------------------------------------------------------
# Options
#-----------------------------------------------------------------------------
DEBUG		=	false
VERBOSE		=	false
SHARED		=	false
OPT		=	opt
STATICLIBEXT	=	a
SHAREDLIBEXT	=	so
LIBEXT		=	$(STATICLIBEXT)
EXEEXTENSION	=	#
REPOSIT		=	#

COMP		=	gnu
CXX		=	g++
CXX_c		=	-c #
CXX_o		=	-o #
LINKCXX		=	g++
LINKCXX_L	=	-L
LINKCXX_l	=	-l
LINKCXX_o	=	-o #
AR		=	ar
AR_o		= #
ZLIB		=	-lz
GMPLIB		=	-lgmpxx -lgmp -lstdc++
ifeq ($(UNAME), Linux)
TIMELIB		=	-lrt
else
TIMELIB		=	
endif


SRCDIR		=	src
BINDIR		=	bin
LIBDIR		=	lib
INCLUDEDIR	=	include
LIBOBJ		=	read_file.o read_file_imc.o  sparse.o unbounded.o expected_time.o expected_reward.o bounded_reward.o sccs.o sccs2.o long_run_average.o debug.o bounded.o long_run_reward.o
ifeq ($(USELP),true)
LIBOBJ		+=	lp.o
endif
BINOBJ		=	main.o

NAME		=	imca
BASE		=	$(OSTYPE).$(ARCH).$(COMP).$(OPT)
BINNAME		=	$(NAME)-$(VERSION).$(BASE)
BINSHORTLINK	=	$(BINDIR)/$(NAME)$(EXEEXTENSION)
BINFILE		=	$(BINDIR)/$(BINNAME)
LIBNAME		=	$(NAME)-$(VERSION).$(BASE)
LIBFILE		=	$(LIBDIR)/lib$(LIBNAME).$(LIBEXT)
LIBSHORTLINK	=	$(LIBDIR)/lib$(NAME).$(LIBEXT)
LIBLINK		=	$(LIBDIR)/lib$(NAME).$(BASE).$(LIBEXT)

OBJDIR		=	obj/O.$(BASE)
BINOBJDIR	=	$(OBJDIR)/bin
LIBOBJDIR	=	$(OBJDIR)/lib
BINOBJFILES	=	$(addprefix $(BINOBJDIR)/,$(BINOBJ))
LIBOBJFILES	=	$(addprefix $(LIBOBJDIR)/,$(LIBOBJ))

TMPDIR		=	tmp

LIBBUILD	=	$(AR)
LIBBUILD_o	=	$(AR_o)
LIBBUILDFLAGS	=	$(ARFLAGS)

ifeq ($(DEBUG),true)
FLAGS		=	-g -Wall -D NDEBUG
else
FLAGS		=	-O2
endif
ifeq ($(USELP),true)
CPPFLAGS	=	-I $(INCLUDEDIR) -I $(LPINCLUDE) 
else
CPPFLAGS	=	-I $(INCLUDEDIR)
endif
CXXFLAGS	=	
BINOFLAGS	=	
LIBOFLAGS	=	
ifeq ($(USELP),true)
LDFLAGS		+=	-L$(LPLIBDIR) $(addprefix -l,$(LPLIB)) 
endif
ARFLAGS		=	cr
DFLAGS		=	-MM
VFLAGS		=	--tool=memcheck --leak-check=yes --show-reachable=yes #--gen-suppressions=yes

LN_s		=	ln -s

ifeq ($(USELP),true)
FLAGS		+=	-D__LPSOLVER__=$(LPSOLVER)
endif


ifeq ($(VERBOSE),false)
.SILENT:	$(LIBLINK) $(LIBSHORTLINK) $(BINLINK) $(BINSHORTLINK) $(BINFILE) $(LIBFILE) $(BINOBJFILES) $(LIBOBJFILES)
endif

all: $(TMPDIR) $(LIBFILE) $(BINFILE) $(LIBLINK) $(LIBSHORTLINK) $(BINLINK) $(BINSHORTLINK) $(SOPLEXLINK) 

#-----------------------------------------------------------------------------
# SHARED Libaries
#-----------------------------------------------------------------------------
ifeq ($(SHARED),true)
CPPFLAGS	+= -fPIC
LIBEXT		=	$(SHAREDLIBEXT)
LIBBUILD	=	$(LINKCXX)
LIBBUILDFLAGS	+=	-shared
LIBBUILD_o	= 	-o # the trailing space is important
ARFLAGS		=
RANLIB		=
endif

CXXFLAGS	+=	$(USRCXXFLAGS)
LDFLAGS		+=	$(USRLDFLAGS)
ARFLAGS		+=	$(USRARFLAGS)
DFLAGS		+=	$(USRDFLAGS)

#-----------------------------------------------------------------------------
# Create directories
#-----------------------------------------------------------------------------
$(OBJDIR):
		@-mkdir -p $(OBJDIR)

$(BINOBJDIR):	$(OBJDIR)
		@-mkdir -p $(BINOBJDIR)

$(LIBOBJDIR):	$(OBJDIR)
		@-mkdir -p $(LIBOBJDIR)

$(BINDIR):
		@-mkdir -p $(BINDIR)

$(LIBDIR):
		@-mkdir -p $(LIBDIR)
		
$(TMPDIR):	
		@-mkdir -p $(TMPDIR)
		
#-----------------------------------------------------------------------------
# making library
#-----------------------------------------------------------------------------		
$(LIBLINK) $(LIBSHORTLINK):	$(LIBFILE)
		@rm -f $@
		cd $(dir $@) && $(LN_s) $(notdir $(LIBFILE)) $(notdir $@)

$(BINLINK) $(BINSHORTLINK):	$(BINFILE)
		@rm -f $@
		cd $(dir $@) && $(LN_s) $(notdir $(BINFILE)) $(notdir $@)

$(BINFILE):	$(BINDIR) $(BINOBJDIR) $(LIBOBJFILES) $(BINOBJFILES) 
		@echo "-> linking $@"
		$(LINKCXX) $(BINOBJFILES) $(LIBOBJFILES) \
		$(LDFLAGS) $(LINKCXX_o)$@

$(LIBFILE):	$(LIBDIR) $(LIBOBJDIR) $(LIBOBJFILES)
		@echo "-> generating library $@"
		-rm -f $(LIBFILE)
		$(LIBBUILD) $(LIBBUILDFLAGS) $(LIBBUILD_o)$@ $(LIBOBJFILES) $(REPOSIT)
ifneq ($(RANLIB),)
		$(RANLIB) $@
endif

#-----------------------------------------------------------------------------
# object compiling and linking
#-----------------------------------------------------------------------------	
$(BINOBJDIR)/%.o:	$(SRCDIR)/%.cpp 
		@-mkdir -p $(BINOBJDIR)
		@echo "-> compiling $@"
		$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(BINOFLAGS) $(FLAGS) $(CXX_c)$< $(CXX_o)$@

$(LIBOBJDIR)/%.o:	$(SRCDIR)/%.cpp $(INCLUDEDIR)/%.h
		@-mkdir -p $(LIBOBJDIR)
		@echo "-> compiling $@"
		$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LIBOFLAGS) $(FLAGS) $(CXX_c)$< $(CXX_o)$@
	
#-----------------------------------------------------------------------------
# cleaning
#-----------------------------------------------------------------------------		
.PHONY: cleanbin
cleanbin:       $(BINDIR)
		@echo "remove binary $(BINFILE)"
		@-rm -f $(BINFILE) $(BINLINK) $(BINSHORTLINK)

.PHONY: cleanlib
cleanlib:       $(LIBDIR)
		@echo "remove library $(LIBFILE)" 
		@-rm -f $(LIBFILE) $(LIBLINK) $(LIBSHORTLINK)

.PHONY: clean
clean:          cleanlib cleanbin $(LIBOBJDIR) $(BINOBJDIR) $(OBJDIR)
		@echo "remove objective files" 
ifneq ($(LIBOBJDIR),)
		@-rm -f $(LIBOBJDIR)/*.o && rmdir $(LIBOBJDIR)
endif
ifneq ($(BINOBJDIR),)
		@-rm -f $(BINOBJDIR)/*.o && rmdir $(BINOBJDIR)
endif
ifneq ($(OBJDIR),)
		@-rm -f $(LASTSETTINGS)
		@-rmdir $(OBJDIR)
endif

