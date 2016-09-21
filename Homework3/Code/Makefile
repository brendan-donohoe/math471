trap.o: trap.f90
	gfortran -c trap.f90

trap.x: trap.o
	gfortran -o trap.x trap.o

trap.txt: trap.x
	  ./trap.x > trap.txt

graphtrap: trap.txt
	   nohup matlab -nodisplay -nosplash -nodesktop < trapgraph.m

clean:
	-rm *.o
	-rm *.x
	-rm *.out
