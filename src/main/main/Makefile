# main/Makefile

NAME := main

OBJS := nmesh_main.o main.o parameters.o skeleton.o tensors.o utilities.o \
	variables.o endianIO.o lists.o timer.o nan_checker.o norms.o

############################################################################
# NOTE: main is special, so we do not use Makefile.subdirs here
############################################################################

# directory of libraries
LIBDIR  = $(TOP)/lib

# directory for the object files of this library
OBJDIR  = $(LIBDIR)/obj/$(NAME)
OBJSWITHPATH = $(addprefix $(OBJDIR)/,$(OBJS))


# main target: build the executable with the linker
# note that $@ is the matched target
$(EXECDIR)/$(EXEC): $(OBJSWITHPATH) $(wildcard $(LIBDIR)/*)
	$(CLINKER) $(CFLAGS) -o $@ $(OBJSWITHPATH) $(LIBS)

# secondary target: compile C files into object files
# note that $@ is the matched target, 
# $* is the stem of the target defined by %
# MAINREV can be used to print the git revision in main
$(OBJDIR)/%.o: %.c nmesh.h nmesh_$(NAME).h $(NAME).h Makefile
	$(CC) $(CFLAGS) -o $@ -c $*.c -DMAINREV="$$(git rev-parse --short HEAD)"


# dependencies on auto-generated files
nmesh.h: nmesh_automatic_include.h
	touch nmesh.h

main.o: nmesh_automatic_initialize.c
