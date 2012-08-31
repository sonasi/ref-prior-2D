#-----------------------------------------------------------------------
# File:	Makefile refpriors
#
# Description: This file contains the instructions to compile 
#              and link the refpriors programs
#
# Created: 24-Jul-2008	Harrison B. Prosper
#-----------------------------------------------------------------------
ifndef ROOTSYS
$(error ROOTSYS not defined - please set up Root)
endif

# Sub-directories
srcdir	:= src
tmpdir	:= tmp
bindir	:= bin
libdir	:= lib
testdir	:= test

# Set this equal to the @ symbol to suppress display of instructions
# while make executes
ifdef verbose
AT 	:=
else
AT	:= @
endif

# Get list of sources to be compiled into applications
appsrcs	:= $(wildcard $(testdir)/*.cc)
# Construct names of object models from list of sources
appobjs	:= $(subst $(testdir)/,$(tmpdir)/,$(appsrcs:.cc=.o))
# Construct list of applications
applications := $(subst $(testdir)/,$(bindir)/,$(appsrcs:.cc=))

# Get list of sources to be made into shared libraries
srcs	:= $(wildcard $(srcdir)/*.cc)
objs	:= $(subst $(srcdir)/,$(tmpdir)/,$(srcs:.cc=.o))
sharedlibs := $(libdir)/libBayesianLikelihood.so

# Display list of applications to be built
say	:= $(shell echo -e "Apps: $(applications)" >& 2)
say	:= $(shell echo -e "Libs: $(sharedlibs)" >& 2)
say     := $(shell echo -e "Srcs: $(srcs)" >& 2)
say     := $(shell echo -e "Objs: $(objs)" >& 2)
 #$(error bye!) 

# Libraries
libs	:=  \
$(shell root-config --libs) -lMathMore -lMathCore -lRooFit \
-L$(libdir) \
-lMinuit
#-----------------------------------------------------------------------

# 	Define which compilers and linkers to use

# 	C++ Compiler
ifdef GCC_DIR
CXX	:= $(GCC_DIR)/bin/g++
LD	:= $(GCC_DIR)/bin/g++
LDSHARED:= $(GCC_DIR)/bin/g++ -shared
else
CXX	:= g++
LD	:= g++
LDSHARED:= g++ -shared
endif
#	C++ Linker

# 	Define paths to be searched for C++ header files (#include ....)

CPPFLAGS:= -I. -I include $(shell root-config --cflags) 

# 	Define compiler flags to be used
#	-c		perform compilation step only 
#	-g		include debug information in the executable file
#	-O2		use 2nd level of optimization
#	-ansi		require strict adherance to C++ standard
#	-Wall		warn if source uses any non-standard C++
#	-pipe		communicate via different stages of compilation
#			using pipes rather than temporary files

CXXFLAGS:= -c -g -O -ansi -Wall -pipe -fPIC

#	Linker flags

LDFLAGS := -g

# 	Libraries
#	System libraries reside in /usr/lib and /lib. They are usually 
#	shared libraries, identified by the file extension ".so".

LIBS	:= -lm $(libs)

#	Rules
#	The structure of a rule is
#	target : source
#		command
#	The command makes a target from the source. 
#	$@ refers to the target
#	$< refers to the source

all:	$(sharedlibs) $(applications)

bin:	$(applications)

lib:	$(sharedlibs)

# Syntax:
# list of targets : target pattern : source pattern

$(sharedlibs)	: $(objs)
	@echo "---> Linking $@"
	@echo  $(AT)$(LDSHARED) $(LDFLAGS) -fPIC $(objs) -o $@
	$(AT)$(LDSHARED) $(LDFLAGS) -fPIC $(objs) -o $@

# Make applications depend on shared libraries to force the latter
# to be built first

$(applications)	: $(bindir)/%	: $(tmpdir)/%.o  $(sharedlibs)
	@echo "---> Linking $@"
	echo $(AT)$(LD) $(objs) $(LDFLAGS) $< $(LIBS) -o $@
	$(AT)$(LD) $(objs) $(LDFLAGS) $< $(LIBS) -o $@

$(appobjs)	: $(tmpdir)/%.o	: $(testdir)/%.cc
	@echo "---> Compiling `basename $<`" 
	@echo $(AT)$(CXX) $(CXXFLAGS) $(CPPFLAGS)  $< -o $@ #>& $*.FAILED
	$(AT)$(CXX) $(CXXFLAGS) $(CPPFLAGS)  $< -o $@ #>& $*.FAILED
	$(AT)rm -rf $*.FAILED

$(objs)	: $(tmpdir)/%.o	: $(srcdir)/%.cc
	@echo "---> Compiling `basename $<`" 
	@echo $(AT)$(CXX) $(CXXFLAGS) $(CPPFLAGS)  $< -o $@ #>& $*.FAILED
	$(AT)$(CXX) $(CXXFLAGS) $(CPPFLAGS)  $< -o $@ #>& $*.FAILED
	$(AT)rm -rf $*.FAILED

# 	Define clean up rules
clean   :
	rm -rf $(tmpdir)/*.o

veryclean   :
	rm -rf $(tmpdir)/*.o $(bindir)/* $(libdir)/*.so


