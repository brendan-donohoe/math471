trap.o: trap.f90
	gfortran -c trap.f90

trap.x: trap.o
	gfortran -o trap.x trap.o

trap.txt: trap.x
	  -rm trap.txt
	  ./trap.x > trap.txt

graphtrap: trap.txt
	   nohup matlab -nodisplay -nosplash -nodesktop < trapgraph.m

gaussQuad.o: gaussQuad.f90
	     gfortran -c gaussQuad.f90

lglnodes.o: lglnodes.f90
	    gfortran -c lglnodes.f90

gaussQuad.x: gaussQuad.o lglnodes.o
	     gfortran -o gaussQuad.x gaussQuad.o lglnodes.o

gaussQuad.txt: gaussQuad.x
	       -rm gaussQuad.txt
	       ./gaussQuad.x > gaussQuad.txt

graphgauss: gaussQuad.txt
	    nohup matlab -nodisplay -nosplash -nodesktop < gaussgraph.m

clean:
	-rm *.o
	-rm *.x
	-rm *.out
