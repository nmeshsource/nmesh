# Makefile
# Wolfgang Tichy 1/2019
# Builds the nmesh executable based on the file MyConfig
# See http://www.gnu.org/software/make/manual for the manual of GNU make

# where am I
UNAME := $(shell uname)
TOP   := $(shell pwd)

# name the fruit of our labor
EXEC = nmesh
EXECDIR = $(TOP)/exe


# variables common to all setups
CC = gcc	# gcc or icc
CXX =		# g++ or icc
CLINKER =	# will be used only in src/main/main/Makefile for linking
AR = ar		# ar command we use to build library from object files

INCS = -I$(TOP)/src/main/main
LIBS = -L$(TOP)/lib
SPECIALINCS =
SPECIALLIBS =
libsys = -lm

DFLAGS =
OFLAGS =
WARN = # -Wall


# --------------------------------------------------------------------------
# different machines and environments

# Linux
OFLAGS = -g


# --------------------------------------------------------------------------
# some libraries are currently required
libpaths = src/main/amr
libpaths += src/main/basis src/main/nMPI


# --------------------------------------------------------------------------
# the user chooses the libraries and some options in the file MyConfig
include MyConfig

# --------------------------------------------------------------------------
# some more libraries are currently required, those need to be last
#libpaths += src/utility/NumericUtils

# --------------------------------------------------------------------------
# set CXX and CLINKER to CC if they are not set in MyConfig
ifeq ($(CXX),)
CXX = $(CC)
endif

ifeq ($(CLINKER),)
CLINKER = $(CC)
endif

# --------------------------------------------------------------------------
# manage how the nmesh sources are compiled

# note that the order matters, e.g.
# main has to go last since it has to be compiled last
libpaths += src/main/main

# extract the list of directory names
libdirs = $(dir $(libpaths))

# extract list of names
libnames = $(notdir $(libpaths))

# make the list of libraries
liblist := $(foreach libname,$(libnames),-l$(libname))

# remove -lmain from that list
liblist := $(subst -lmain,,$(liblist))

# make final list of libraries that is passed to the linker
# uses standard hack to resolve interdependencies by repeating the libraries
# system libraries go in the end
LIBS += $(MPIDIRL) $(liblist) $(liblist)
LIBS += $(SPECIALLIBS) $(libsys) $(MPILIBS)


# make the list of include files that will be automatically included for each
# module
libincludes := $(foreach libpath,$(libpaths),\
	$(libpath)/nmesh_$(notdir $(libpath)).h)

# define the automatic configuration files
autoinclude = src/main/main/nmesh_automatic_include.h
autoinitial = src/main/main/nmesh_automatic_initialize.c
autotext    = \/\* automatically generated from MyConfig \*\/



# --------------------------------------------------------------------------
# some of the above variables are meant to be global, so pass them on
# to the shell 
CFLAGS = $(DFLAGS) $(OFLAGS) $(INCS) $(MPIDIRI) $(HDF5DIRI) $(SPECIALINCS) $(WARN)
export


# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# default target
nmesh: $(autoinclude) $(autoinitial)
	@echo CC=$(CC)
	@echo CXX=$(CXX)
	@echo CLINKER=$(CLINKER)
	for X in $(libnames); do mkdir -p lib/obj/$$X; done
	for X in $(libpaths); do $(MAKE) -C $$X; done


# --------------------------------------------------------------------------
# other targets 

# if there is no MyConfig file, use the example provided in doc
MyConfig:
	-if test ! -f MyConfig; then cp doc/MyConfig.example MyConfig; fi 

# automatic configuration files
$(autoinclude): MyConfig
	echo $(autotext) > $(autoinclude) 
	for X in $(libincludes); do \
	  echo \#include \"$(TOP)/$$X\" >> $(autoinclude); \
	done

$(autoinitial): MyConfig
	echo $(autotext) > $(autoinitial) 
	for X in $(libnames); do \
	  echo int nmesh\_$$X\(struct tMESH *\)\; >> $(autoinitial); \
	done
	echo "/* call nmesh initialization functions: */" >> $(autoinitial); \
	for X in $(libnames); do \
	  echo nmesh\_$$X\(mesh\)\; >> $(autoinitial); \
	done


# create tar file
tar:
	cd ..; tar czf nmesh.tgz --exclude lib --exclude exe ./nmesh

# take a fresh look at things
clean:
	-rm -r lib
	-rm $(autoinclude)
	-rm $(autoinitial)

# remove emacs backup files
cleantilde:
	find . -name "*~" -exec rm {} \;

# remove code that is not needed once the corresponding libs have been built
rm_MemoryMan_code:
	find src/main/MemoryMan/ -name "*.c*" -print -exec rm -rf '{}' \;
	find src/main/MemoryMan/ -name "*.m*" -print -exec rm -rf '{}' \;
	echo -e "ls:\n\tls *.h" > donothing_Makefile
	find src/main/MemoryMan/ -name "Makefile" -print -exec cp -f donothing_Makefile '{}' \;
	rm -f donothing_Makefile

rm_Math_code:
	rm -rf src/Math

rm_utility_code:
	find src/utility/ -name "*.c*" -print -exec rm -rf '{}' \;
	find src/utility/ -name "*.m*" -print -exec rm -rf '{}' \;
	echo -e "ls:\n\tls *.h" > donothing_Makefile
	find src/utility/ -name "Makefile" -print -exec cp -f donothing_Makefile '{}' \;
	rm -f donothing_Makefile

rm_projects_code:
	find src/projects/ -name "*.c*" -print -exec rm -rf '{}' \;
	find src/projects/ -name "*.m*" -print -exec rm -rf '{}' \;
	echo -e "ls:\n\tls *.h" > donothing_Makefile
	find src/projects/ -name "Makefile" -print -exec cp -f donothing_Makefile '{}' \;
	rm -f donothing_Makefile

rm_code:
	make rm_MemoryMan_code
	make rm_Math_code
	make rm_utility_code
	make rm_projects_code
