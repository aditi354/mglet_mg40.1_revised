#
#  Makefile fuer Library QPack
#  dated:  Mi Jul 20 13:41:03 CEST 2005
#
.SUFFIXES: .o .f90
#
####################  Welche Shell fuer Kommandos:
SHELL = /bin/sh
#######################################################
#include ../mglet.altix
include ../mglet.pc
#F77=mpif90 -model=F1 -lp64 -64
#FFLAGS= -c -oplist -msg=w -conti199 -nolimit -I../../sources \
#           -precexp=4  -64 -noprecexpio  -halt=e 
#OPTIONS= -opt=4 -noscope -optlog 
#F77=/opt/mpi/bin/mpif90
FC= $(F90)
#FFLAGS= -c -I../../sources +DD64

#OPTIONS= -I../../sources +DD64

# FFLAGS= -c -oplist -msg=w -conti199 -nolimit -I../MGLET -precexp=4 -procnum=8 -64 -noprecexpio 
#FLAGS= -c -oplist -msg=w -conti199 -nolimit -I../MGLET -precexp=4 -64 -noprecexpio -noparallel
#FFLAGS= -c -oplist -nolimit -msg=w -conti199 -I../MGLET -precexp=4 -procnum=8 -64 
#NOPTIONS= -opt=0
# OPTIONS= -opt=4 -noscope -optlog 
#OPTIONS= -opt=ss -noscope -optlog 
# LOPTIONS=+BTLB -rdma -64 
#LOPTIONS=+BTLB -rdma -lpl -64 -noparallel
#LIBRARY= -L./ -lQPack
#
OBJECTS= \
blockgrid.o \
blockquad.o \
blockquad_searcharea.o \
blockquad_geomlim.o \
getquad.o \
xpolquad.o \
blocking_xtrpol.o \
blocking.o \
blockingadd.o \
bubvbw.o \
calcsdiv.o \
compflux.o \
outflux.o \
rlxold.o \
setold.o \
setzero.o \
stencil.o \
stencil_output.o \
stencil_determine.o \
stencil_least_square.o \
stencil_least_square_adiabat.o \
stencil_lagrange.o \
stencil_spline3.o \
xtrpol.o \
xtrpolp.o \
solint.o
#
###############################  Compilier-Aufruf:
#
lib:	objs
	ar ruv libQPack.a $(OBJECTS)
	ranlib libQPack.a

objs:	$(OBJECTS)

.f90.o:
	$(FC) $(FFLAGS) $(OPTIONS) -I../../sources $<

.DEFAULT:
	@echo $@": not known"
clean:
	rm *.o libQPack.a
