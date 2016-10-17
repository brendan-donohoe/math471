++++++++++++++++++++++++++++++++++++++++++++++++++++
README file for Tim Chavez, Brendan Donohoe Math 471
++++++++++++++++++++++++++++++++++++++++++++++++++++

1. Homework 1:

  - The report is in the directory Homework1/Report and is called homework1.pdf.

  - The code for this homework is located in Homework1/Code and is documented in the appendix of the report.

  - See the report for details on how to compile and run the code.

2. Homework 2:

  - The report for Homework 2 is in the directory Homework2/Report and is called homework2.pdf.

  - The code for the 2nd homework is in Homework2/Code and is included in the appendix of the report.

  - To run our modified Newton's method from the directory Homework2/Code, type:
  $ perl newtonS.p

3. Homework 3:

  - As usual, the report for Homework 3 is in the directory Homework3/Report and is called
  homework3.pdf.

  - The code for Homework 3 is in Homework3/Code.

  - To run our Trapezoidal rule code and graph the results, use the commands in the makefile:
  $ make clean
  $ make graphtrap

  - To run our Gauss Quadrature code and graph the results, use the commands in the makefile:
  $ make clean
  $ make graphgauss

4. Homework 4

  - The report for Homework 4 is in the directory Homework4/Report and is named homework4.pdf.

  - All of the code for Homework 4 is in Homework4/Code.

  - To run the code and graph a function on curvilinear coordinates, use the commands:
  $ make clean
  $ make graph_it

  - To run the code and plot the error, use the commands:
  $ make clean
  $ make graph_error

  - To make changes to the functions being graphed and analyzed for error, edit:
  The FUNCTIONS area in homework4.f90 (lines 60-70)
  The FUNCTIONS area in erroranalysis.f90 (lines 82-100)
  The FUNCTIONS areas in xycoord.f90 (lines 10-20 for x_coord and lines 27-37 for y_coord)

5. Homework 5

  - The report for Homework 5 is in the directory Homework5/Report and is named homework5.pdf. The movies discussed in the report are contained in the AngryBirdsMovies.zip folder.

  - All of the code for Homework 5 is in Homework5/Code and is split into 4 seperate directories for the different parts of the homework.

  - To run one part of the code, run setup.m in the desired folder in the Code folder. Then, play the AVI movie file that is generated called setup1.avi.
