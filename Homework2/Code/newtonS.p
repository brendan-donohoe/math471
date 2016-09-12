#!/usr/apps/bin/perl
#
# perl program to try the convergence for different functions 
# when using Newtoon's method for f(x) = 0
# run with : perl newton.p
#
#
# Here is the generic file
$cmdFile="./newtonS.f90.Template";
$outFile="./newtonS.f90";

# Stuff to converge over

@array_f = ("x", "x*x", "sin(x)+cos(x*x)", "cos(x)-1.d0", "exp(x)-4.d0*x+8.d0*log(2.d0)-4.d0");
@array_fp = ("1.d0", "2.d0*x", "cos(x)-2.d0*x*sin(x*x)", "-sin(x)", "exp(x)-4.d0");
@array_guess = ("1.d0", "1.d0", "1.d0", "1.d0", "1.d0");

for( $m=0; $m < 4; $m = $m+1){ #TODO: Change m to 3 when tests are finished.
    # Open the Template file and the output file. 
    open(FILE,"$cmdFile") || die "cannot open file $cmdFile!" ;
    open(OUTFILE,"> $outFile") || die "cannot open file!" ;
    # Setup the outfile based on the template
    # read one line at a time.
    while( $line = <FILE> )
    {
	# Replace the the stings by using substitution
	# s
	$line =~ s/\bFFFF\b/$array_f[$m]/g;
	$line =~ s/\bFPFP\b/$array_fp[$m]/g;
	print OUTFILE $line;
        # You can always print to secreen to see what is happening.
        # print $line;
    }
    # Close the files
    close( OUTFILE );
    close( FILE );
    
    # Run the shell commands to compile and run the program
    system("gfortran $outFile");
    system("./a.out > apa.txt");

    open(FILE,"apa.txt") || die "cannot open file" ;
    # Setup the outfile based on the template
    # read one line at a time.
    while( $line = <FILE> )
    {
        $line =~ s/\s+/ , /g;
	print substr($line, 2, -2) . "\n";
    }
    close( FILE );
}

exit

