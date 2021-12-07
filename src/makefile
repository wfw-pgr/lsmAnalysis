target  = ../main
object  = \
  variablesMod.o \
  utilitiesMod.o \
  lstStructMod.o \
  allocatorMod.o \
  bTriPrismMod.o \
  cglsdenseMod.o \
  invSolverMod.o \
  updateEqsMod.o \
  initiatorMod.o \
  ioUtilityMod.o \
  main.o

FC = gfortran

FFLAGS  = -O3 -llapack -lblas
LDFLAGS =


.SUFFIXES: .o .f90

%.o: %.f90
	$(FC) -c $< $(FFLAGS) $(LDFLAGS)

%.mod: %.f90 %.o
	@:

$(target): $(object)
	$(FC) -o $@ $(object) $(FFLAGS) $(LDFLAGS)


.PHONY:clean all

clean:
	$(RM) $(target) $(object) *.mod

all: clean $(target)
