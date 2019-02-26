# main/Makefile

NAME := main

OBJS := nmesh_main.o main.o parameters.o skeleton.o tensors.o utilities.o \
	variables.o endianIO.o


# main is special, so we do not use Makefile.subdirs here

# directory of the libraries
LIBDIR  = $(TOP)/lib

# directory for the object files of this library
OBJDIR  = $(LIBDIR)/obj/$(NAME)
OBJSWITHPATH = $(addprefix $(OBJDIR)/,$(OBJS))


# main target is the executable
# note that $@ is the matched target
$(EXECDIR)/$(EXEC): $(OBJSWITHPATH) $(wildcard $(LIBDIR)/*)
	$(CLINKER) $(CFLAGS) -o $@ $(OBJSWITHPATH) $(LIBS)

# secondary target: build object files
# note that $@ is the matched target, 
# $* is the stem of the target defined by %
$(OBJDIR)/%.o: %.c nmesh.h nmesh_$(NAME).h $(NAME).h Makefile
	$(CC) $(CFLAGS) -o $@ -c $*.c -DMAINREV="$$(git rev-parse --short HEAD)"


# dependencies for automatic files
nmesh.h: nmesh_automatic_include.h
main.o: nmesh_automatic_initialize.c
