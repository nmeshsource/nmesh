# Makefile
# Wolfgang Tichy 1/2019, Bernd Bruegmann 5/02
# Builds the nmesh executable by including the file MyConfig
# See http://www.gnu.org/software/make/manual for the manual of GNU make

# top level dir
TOP := $(shell pwd)

# name of program, location of executable and extra projects
EXEC = nmesh
EXECDIR = $(TOP)/exe
RELAPROJECTDIR = src/projects
PROJECTDIR = $(TOP)/$(RELAPROJECTDIR)

# default for variables used in all cases
CC = gcc	# gcc or icc
CXX =		# g++ or icpc
CLINKER =	# will be used only in src/main/main/Makefile for linking
AR = ar		# ar command we use to build library from object files

# defaults for option-flags and defined-flags for compiler
OFLAGS = -g
WARN = -Wall
DFLAGS =

# some include and library defaults
INCS = -I$(TOP)/src/main/main
LIBS = -L$(TOP)/lib
SPECIALINCS =
SPECIALLIBS =
libsys = -lm

# dirs for MPI given with -I and -L to compiler, and MPI lib given with -l
MPIDIRI =	# -I/usr/lib/x86_64-linux-gnu/openmpi/include	#for openmpi
MPIDIRL =	# -L/usr/lib/x86_64-linux-gnu/openmpi/lib	#for openmpi
MPILIBS =	# -lmpi						#for openmpi

# --------------------------------------------------------------------------
# some nmesh libraries are required
libpaths = src/main/amr src/main/nMPI
libpaths += src/main/basis src/main/coordinates src/main/evolve
libpaths += src/main/dg src/main/limiter src/main/checkpoint
libpaths += src/utility/output src/utility/numerics src/utility/units

# --------------------------------------------------------------------------
# we can choose more libraries and options in the file MyConfig
projects =#
include MyConfig

# --------------------------------------------------------------------------
# set variable projectnames from git targets in projects. NOTE: there must
# be a / just before the actual projectname (e.g. mars.fau.edu:/DNSdata)
projectbasenames = $(basename $(projects))
projectnames = $(notdir $(projectbasenames))

# set projectpaths and add them to libpaths
projectpaths = $(addprefix $(RELAPROJECTDIR)/,$(projectnames))
libpaths += $(projectpaths)

# --------------------------------------------------------------------------
# some more required libraries that need to be last
#libpaths += src/utility/NumericUtils

# --------------------------------------------------------------------------
# remove all paths in libpaths that do not contain a corresponding
# nmesh_LIB.c file:

# get all nmesh_LIB.c C-files we can find in libpaths
libpathCfiles := $(foreach X, $(libpaths),\
                   $(wildcard $(X)/nmesh_$(notdir $(X)).c))

# use list of C files to make paths, afterwards remove trailing slashes
libpaths := $(dir $(libpathCfiles))
libpaths := $(patsubst %/, %, $(libpaths))

# make list of all non-main modules to find possible physics parameters
compiledphysics := $(patsubst src/main/%,,$(libpaths))
compiledphysics := $(patsubst src/utility/%,,$(compiledphysics))
compiledphysics := $(notdir $(compiledphysics))

# --------------------------------------------------------------------------
# set CXX and CLINKER to CC if they are not set in MyConfig
ifeq ($(CXX),)
CXX = $(CC)
endif

ifeq ($(CLINKER),)
CLINKER = $(CC)
endif

# --------------------------------------------------------------------------
# manage how the nmesh sources are compiled and linked

# note that the order matters, e.g.
# main has to be last since it has to be linked last
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
# we resolve library interdependencies by repeating the libraries
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
autotext    = \/\* automatically generated from Makefile and MyConfig \*\/

# important header files (when they change we want a recompile)
mainheaders = src/main/amr/nmesh_amr.h src/main/amr/thread_defs.h
mainheaders += src/main/main/skeleton.h src/main/main/variables.h
mainheaders += src/main/nMPI/nMPI_defs.h src/main/evolve/evosys.h
mainheaders += src/main/coordinates/CI.h
mainheaders += src/main/amr/nmesh_amr_defs.h src/main/amr/nmesh_amr_loops.h

# --------------------------------------------------------------------------
# some of the above variables are meant to be global, so we pass them on
# to the shell with export
CFLAGS = $(DFLAGS) $(OFLAGS) $(WARN) $(INCS) $(MPIDIRI) $(HDF5DIRI) $(SPECIALINCS)
export

# --------------------------------------------------------------------------
# Make Targets
# --------------------------------------------------------------------------
# default target
nmesh: $(autoinclude) $(autoinitial) .git/hooks/pre-commit
	@echo
	@echo ======================= Compiling nmesh ========================
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
$(autoinclude): MyConfig $(mainheaders)
	@echo ==================== Auto generating files =====================
	@echo $(autotext) > $(autoinclude)
	for X in $(libincludes); do \
	  echo \#include \"$(TOP)/$$X\" >> $(autoinclude); \
	done
	@for X in $(libnames); do \
	  echo int nmesh\_$$X\(struct tMESH *\)\; >> $(autoinclude); \
	done
	touch $(TOP)/src/main/main/nmesh.h

$(autoinitial): MyConfig
	@echo $(autotext) > $(autoinitial)
	@echo "/* set compiledphysics string: */" >> $(autoinitial);
	@echo compiledphysics=\"$(compiledphysics)\"\; >> $(autoinitial);
	@echo "/* call nmesh initialization functions: */" >> $(autoinitial);
	for X in $(libnames); do \
	  echo nmesh\_$$X\(mesh\)\; >> $(autoinitial); \
	done


# create tar file
tar:
	cd ..; tar czf nmesh.tgz --exclude lib --exclude exe ./nmesh

# take a fresh look at things
cleanonly:
	@echo ====================== Removing files ==========================
	-rm -rf lib
	-rm -f $(autoinclude)
	-rm -f $(autoinitial)

clean: install_git_hooks cleanonly

# remove joe/emacs backup files
cleantilde:
	find . -name "*~" -exec rm {} \;

# print some vars
printvars:
	@echo projects=$(projects)
	@echo PROJECTDIR=$(PROJECTDIR)
	@echo projectpaths=$(projectpaths)
	@echo projectnames=$(projectnames)
#	@echo libpathCfiles=$(libpathCfiles)
	@echo libpaths=$(libpaths)
	@echo compiledphysics=$(compiledphysics)


# targets to get git projects
git_clone:
	@echo ==================== Cloning nmesh projects ====================
	-for X in $(projects); do N=$$(basename $$X .git); printf "==== %s ====\n" $$N; git clone $$X $(PROJECTDIR)/$$N; done
	@$(MAKE) install_git_hooks

git_pull: install_git_hooks
	@echo ====================== main part of nmesh ======================
	git pull
	@echo ======================== nmesh projects ========================
	@for X in $(projects); do N=$$(basename $$X .git); if [ -d "$(PROJECTDIR)/$$N" ]; then printf "==== %s ====\n" $$N; cd $(PROJECTDIR)/$$N; git pull; fi done

git_status:
	@echo ====================== main part of nmesh ======================
	git status -uno
	@echo ======================== nmesh projects ========================
	@for X in $(projects); do N=$$(basename $$X .git); if [ -d "$(PROJECTDIR)/$$N" ]; then printf "==== %s ====\n" $$N; cd $(PROJECTDIR)/$$N; git status -uno; fi done

# targets for git hooks
.git/hooks/pre-commit: git_hooks/pre-commit
	@$(MAKE) install_git_hooks

install_git_hooks:
	@echo ==================== Installing git hooks ======================
	cp git_hooks/pre-commit .git/hooks
	for X in $(projects); do N=$$(basename $$X .git); if [ -d "$(PROJECTDIR)/$$N/.git/hooks" ]; then cp git_hooks/pre-commit $(PROJECTDIR)/$$N/.git/hooks; fi done


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
