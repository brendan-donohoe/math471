# An advanced makefile for the differentiation demo example
# Here we use the macros FC, FFLAGS, F90FLAGS and general rules for compiling 
# .f90 and .f files
# We also use the Macro OBJECTS to not have to keep track of the 
# object files in more than one place.
# The PHONY target clean has also been added

FC = gfortran
LD = gfortran
LDFLAGS = -fopenmp
FFLAGS =
F90FLAGS = -fopenmp
EX = ./homework4.x


OBJECTS = homework4.o differentiate.o printdble.o xycoord.o trap.o
MODULES = xycoord.mod

.PHONY: clean

# Compile, run, process and open.
$(EX): $(MODULES) $(OBJECTS) 
	$(LD) $(LDFLAGS) $(OBJECTS) -o $(EX)

graph_it: $(EX)
	$(EX) > out.txt
	nohup matlab -nosplash -nodisplay < plotgrid.m > output.txt
	open -a preview grid.eps

run_it: $(EX)
	$(EX)

%.o : %.f90
	$(FC) $(F90FLAGS) -c $<

%.o : %.f
	$(FC) $(FFLAGS) -c $<

%.mod : %.f90
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f $(OBJECTS) $(MODULES) $(EX) *~ erroranalysis.o erroranalysis.x nohup.out errorvals.txt errorgraph.png grid.eps output.txt out.txt x.txt y.txt