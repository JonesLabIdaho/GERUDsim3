Using GerudSim3 on Linux:

1.  Create a GerudSim3 folder.
2.  Place the GerudSim3 ".cpp" and ".h" files within this folder.
3.  Make sure you have the g++ compiler installed (most easily accomplished by using the Synaptic Package Manager -- look for g++ on the list -- Version 4.6 works on Ubuntu). Compile GerudSim3 within this folder with the following command:

	>g++ -o GS3 GerudSim3.cpp

	You should see a file called GS3 appear.  This file is the executable program.

4.  Prepare the allele frequency file.  The name and format of the file must be the same as the included sample file.  The file name should be "GS3allelefrequencies.txt", and it is case sensitive.  

The first line of the file should indicate the number of loci in the file, as in "4 loci".

Each locus name should be preceded by an exclamation point "!" and should occupy its own line, followed by one line for each allele at that locus.  The line corresponding to each allele should start with the allele name, represented by a three-digit integer, followed by a tab character, followed by the allele frequency.  THERE SHOULD BE NO SPACES ON THESE LINES, AND ONLY ONE TAB PER LINE!!!!!  The allele frequency should be immediately followed by a carriage return.

5.  Prepare the parameters file.  The name of the parameter file must be GS3parameters.txt.  

Line 1: Should say either "Known Mother" or "Unknown Mother", depending on whether or not the mother's genotype is known.

Line 2: Indicate the number of loci to use.  If this line says "2 loci", for instance, then the simulation will use the first two loci from the allele frequency file.  GerudSim3 can use up to 50 loci.

Line 3: Indicate the number of iterations for each combination of paramters.  "1000 iterations" will run 1000 simulations for each combination of parameters.

Line 4: Indicate the number of runs.  Each run will be performed under whatever parameters are specified in lines 6 onward.  A given run will be repeated as many times as you specified on line 3.  The runs must be repeated to some degree, preferably hundreds or thousands of time, because summary statistics are calculated from the iterations.  However, for many loci or many parents, runs could take a very long time.

Line 5: This line is ignored, but it's useful for organizing the parameter values for each run.

Lines 6 and onward:  Each line from line 6 onward provides the information about the progeny array for each run.  Columns must be separated by a single tab character.  The first column merely numbers the runs and is ignored.  The second column indicates the number of progeny genotyped from the progeny array.  The remaining columns indicate the number of offspring sired by each male in the progeny array.  In the sample file, for instance, 60 offspring were genotypes and two males each contributed 100 offspring to the brood in run 1.  If any male's contribution is set to zero, then subsequent males will be ignored.  A maximum of 10 males can be simulated using GerudSim3.


6.  NOW YOU ARE READY TO RUN THE PROGRAM!  Make sure your current directory is the one containing the following files: "GS3allelefrequencies.txt" and "GS3parameters.txt".  You can either put the compiled GS3 program in your path (your "bin" folder for instance) or you can have it in the same folder as the allele frequency and parameters files.

You have several options regarding how to run GerudSim3, but it takes no parameters, so it's pretty easy.  I'll assume that you don't want to know too much about Linux, so you've compiled the program in the same folder in which you want to run it.

Gerudsim3 will write one new file, called "GS3output.txt".  If there is already a file with that name, GerudSim3 will overwrite it, so if you want to run the program twice, be sure to rename the file between runs or run it in different directories each time.

Gerudsim3 will also output some information to the screen, including the parameter files and allele frequency data.  Be sure to look at this information to make sure there were no errors with your files.  If Gerudsim3 encounters any problems, it will also output errors here.

If you want to just run it interactively, which will make your console non-responsive while it's running, type:

./GS3

If you want the output from the screen to instead be written to a log file, type:

./GS3 > GS3log.txt

This command will redirect the output to a file called GS3log.txt

If you want the output from the screen to be written to a log file, and you want the program to run in the background, type:

./GS3 > GS3log.txt &

If you have any trouble, email me at agjones@tamu.edu.

 


