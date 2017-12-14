
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <ctime>
#include "AlleleFrequencyData.h"
#include "GerudSimParameters.h"
#include "GS3functionsclasses.h"
#include "MTwisterFunctions.h"
using namespace std;

// This file is part of GerudSim3.  This version is version 3.1, updated to compile
// on Visual Studio 2012.  Any changes should be indicated here.  This version is ANSI compliant and
// should compile on any standard C++ compiler.  This file was last altered on October 13, 2015.


// Large Global Arrays Need to be declared globally here
int WinnerDad[1000000][10]; // This is an array of dad indices compatible with the progeny array [number of lists of winning dads][max # dads in a solution]
int NumWinDadsMax = 999999;
int EmbryoPat1[50][1000];   // [numberloci][numberembryos]
int EmbryoPat2[50][1000];	// allows up to 50 loci and 1000 embryos in the progeny array (assayed)
int EmbryoAllele1[50][1000];
int EmbryoAllele2[50][1000];
int EmbryoPaternal1[50][1000];
int EmbryoPaternal2[50][1000];

int main( int argc, char* argv[] )
{
	// Set default filenames 

	string tempallelefile, tempparamfile, tempoutfile, temponearrayfile;
	string tempstring1, tempstring2;
	bool verbose;
	bool save_one_array;

	char* AlleleFrequencyFilename = new char[256];
	char* ParameterFilename = new char[256];
	char* OutputFilename = new char[256];
	char* OneArrayFilename = new char[256];

	tempallelefile = "Gerud3allelefrequencies.txt";
	tempoutfile = "Gerudsim3output.txt";
	tempparamfile = "Gerudsim3parameters.txt";
	temponearrayfile = "GS3onearray.txt";
	verbose = true;
	save_one_array = false;

	int iii, jjj;
	for (iii = 0; iii < argc-1; iii++)
	{
		tempstring1 = argv[iii];
		tempstring2 = argv[iii+1];
		if (tempstring1 == "-a")
			tempallelefile = tempstring2;
		if (tempstring1 == "-o")
			tempoutfile = tempstring2;
		if (tempstring1 == "-p")
			tempparamfile = tempstring2;
		if (tempstring1 == "-v")
		{
			if (tempstring2[0] == 'n' || tempstring2[0] == 'N')
				verbose = false;
		}
		if (tempstring1 == "-one_array")
		{
			temponearrayfile = tempstring2;
			save_one_array = true;
		}
	}

	// convert filenames from strings to character arrays
	for (iii = 0; iii < static_cast<int>(tempallelefile.length()); iii++)
		AlleleFrequencyFilename[iii] = tempallelefile[iii];
	AlleleFrequencyFilename[iii] = '\0';

	for (iii = 0; iii < static_cast<int>(tempoutfile.length()); iii++)
		OutputFilename[iii] = tempoutfile[iii];
	OutputFilename[iii] = '\0';

	for (iii = 0; iii < static_cast<int>(tempparamfile.length()); iii++)
		ParameterFilename[iii] = tempparamfile[iii];
	ParameterFilename[iii] = '\0';

	for (iii = 0; iii < static_cast<int>(temponearrayfile.length()); iii++)
		OneArrayFilename[iii] = temponearrayfile[iii];
	OneArrayFilename[iii] = '\0';

	if (verbose)
	{
		cout << "\n" << "Filenames:\n";
		cout << "Allele Frequencies:\t" << AlleleFrequencyFilename << "\n";
		cout << "Parameter File:    \t" << ParameterFilename << "\n";
		cout << "Output File:       \t" << OutputFilename << "\n";
	}

	int AFflag;
	int Paramflag;
	GerudSimParameters GSparam;
	Paramflag = GSparam.LoadParameters(ParameterFilename);

	AlleleFrequencies AlleleFreq;
	AFflag = AlleleFreq.LoadAlleleFreqs(AlleleFrequencyFilename);

	bool AllelesLoaded;
	AllelesLoaded = true;

	// Open output file for writing

	ofstream o_file;
	o_file.open(OutputFilename);

	o_file << "GERUDsim3 RESULTS:\n\n";

	if (AFflag == 0)
	{
		o_file << "Error loading allele frequencies.";
		o_file.close();
		return 0;
	}

	if (Paramflag = 0)
	{
		o_file << "Error loading parameters.";
		o_file.close();
		return 0;
	}

	if (save_one_array)
	{
		o_file << "No output -- one progeny array saved in a separate file.";
		o_file.close();
	}
	

	// The rest of the program is pretty linear, so most of it 
	// can go in main().  Things are a little weird, because this is
	// a modification of an older program.  Thus, there's a little bit
	// of redefining and moving of things in ways that don't make
	// perfect sense.

	// Output allele frequencies

	if (verbose)
	{
		cout << "\nAllele Frequency File:\n";
		cout << AlleleFreq.NumberLoci << " Loci\n";
		for (iii = 0; iii < AlleleFreq.NumberLoci; iii++)
		{
			cout << AlleleFreq.LocusName[iii] << "\n";
			for (jjj = 0; jjj < AlleleFreq.NumberAlleles[iii]; jjj++)
			{
				cout << AlleleFreq.AlleleName[jjj][iii] << "\t" << AlleleFreq.AlleleFreq[jjj][iii] << "\n";
			}
		}
	}

	// Output Parameter File

	if (verbose)
	{
		cout << "\nParameter File Contents:\n";
		if (GSparam.KnownMother)
			cout << "Known Mother\n";
		else
			cout << "Unknown Mother\n";

		cout << GSparam.NumberLoci << " Loci\n";
		cout << GSparam.NumberIterations << " Iterations\n";
		cout << GSparam.NumberRuns << " Runs\n";

		for (iii = 0; iii < GSparam.NumberRuns; iii++)
		{
			cout << "Run " << iii+1 << "\t" << GSparam.EmbryosAssayed[iii] << " Embryos Assayed\t";
			for (jjj = 0; jjj < 10; jjj++)
			{
				cout << GSparam.EmbryosFromMale[iii][jjj] << " from Male" << jjj+1 << "\t";
			}
			cout << "\n";
		}
	}

	// Variables:

	int MaxNumberLoci;
	int MaxNumberProgenyAssayed;

	MaxNumberLoci = 50;
	MaxNumberProgenyAssayed = 1000;

	Dad RealDad[10];
	Dad SolvedDad[10];
	ParentSet RealSolution;
	Locus *Loci;
	Dad *GoodDads;
	Dad *PossDads;
	Dad *SomeDads;
	Dad *DadBuff;
	Dad *Mothers;
	ParentSet *Solutions;
	ParentSet *BestSolutions;
	ParentSet *BestSolutionsOld;

	 // variables for simulations
     int RealMaleNumber[100];	// [Number of Runs]
     int NumTarB[100];			// [Number of Runs]
     int NumTar[100];			// [Number of Runs]
     int NumIterations;
     int NumFromMale[10][100];   // [Max No Fathers in Progeny Array][Number of Runs]
     bool DoRun[100];
     bool ErrorGP, ErrorPP;
     bool BigDadsSolved;

     int MAtotal[50];			//[Number of Loci]
     int MAalleles[50][200];	//[Number of Loci][Number of Maternal Alleles in Progeny Array]
     int TotalNumberSolutions;
     int NumCompMoms;
     bool TotalNumberSolutionsError;
     bool MAnumberError;
     bool MomKnown;
     int RealMomAllele1[50];  // [Number of Loci]
     int RealMomAllele2[50];  // [Number of Loci]
     int MaxNumDadsEx;

     int NumLociLoaded;
     int NumberEms[100];   // [Number of Runs]
     int NumFreqLoci;
     int NumberLoci;
     int EmbryoRealFather[1000];  // [Number Embryos Assayed]
     int NumberEmbryos;
     int run, iteration;
     int MomAllele1[50];		// [Number of Loci]
     int MomAllele2[50];		// same
     int EmsSampFromMale[10];   // [Max Number of Fathers]

     int MinDadsOneLocus, MinDadsOL;
     int MinDadsOLLess, MinDadsOLMore;
     int MinPatAllNum[50];   // [Number of Loci]
     int PatAllList[50][200]; // [Number of Loci][Number of Paternal Alleles in the Progeny Array] 
     int PatAllNum[50];			// [Number of Loci]
     int NumGPGens;
     bool DadsSolved;
     int DadTestNum;
     int NumWinDads;
     int BestDads[51][10];        // [Maximum Number of Ranked Solutions][Max Number of Fathers]
     long double BestDadsPriority[51]; // [Max Number of Ranked Solutions]
     int NumBestDads;
     long double TopPriority;
     int BestByPriority;

     int NumLoci, ii;

     // checking results variables
     bool CorrectDadNumber;
     int CumCorrDadNum;
     bool OverDadNumber;
     int CumOverDadNum;
     bool UnderDadNumber;
     int CumUnderDadNum;
     int NumberDadsCorrect, CumNumDadsCorr;
     int NumberDadsWrong, CumNumDadsWrong;
     int NumberFiveCorrect, CumNumFiveCorr;
     int NumberFiveWrong, CumNumFiveWrong;
     int NumberTenCorrect, CumNumTenCorr;
     int NumberTenWrong, CumNumTenWrong;
     int NumberTwentCorrect, CumNumTwentCorr;
     int NumberTwentWrong, CumNumTwentWrong;
     int NumberMoreCorrect, CumNumMoreCorr;
     int NumberMoreWrong, CumNumMoreWrong;
     int NumberNotSampledWrong, CumNumNSampWrong;
     int NumberNotSampledCorrect, CumNumNSampCorr;
     bool ExactlyRight;
     int CumExactlyRight, CumGenosRight, CumGenosWrong;
     bool GenosRight;
     int NumOneSolnCorr, NumOneSolnWrong;
     int NumTenSolnsCorr, NumTenSolnsWrong;
     int NumMoreSolnsCorr, NumMoreSolnsWrong;
     int NumNotSolved;
     int CumMotherCorrect;
     bool MotherCorrect;
	 bool NumWinDadError;

	 // The loop will run through the different parameter values
	 // For each set of parameter values, the model will simulate a
	 // progeny array, run the GERUD algorithm on it, and collect the results


	 // At this point, the simulation parameters have been loaded
	 // as have the allele frequencies.

	 // First allocate memory for the various pointers that will be needed.
	 // Realistically, the program doesn't need that much memory.

	 Loci = new Locus[50];  // allow up to 50 loci, although this will almost never be feasible
	 GoodDads = new Dad[1000000];  // allow up to 1,000,000 dads
	 PossDads = new Dad[2000000];
	 long int GPMax =  999999;
	 long int PPMax = 1999999;
	 Mothers = new Dad[200000];
	 long int MotherNumberMax = 199999;
	 Solutions = new ParentSet[100000];
	 long int TotalNumberSolutionsMax = 99999;
	 BestSolutions = new ParentSet[100];
	 BestSolutionsOld = new ParentSet[100];
	 Dad * PossMoms = new Dad[200000];
	 long double *WDChiSqr = new long double[TotalNumberSolutionsMax+1];
	 long double *WDFreq = new long double[TotalNumberSolutionsMax+1];
	 long double *WDPriority = new long double[TotalNumberSolutionsMax+1];


	 int i, j, k, l, m;
	 Dad tempDads[10];
	 double db1, db2, rndbl, dbsum;
	 string tempLine, tempLine1;

	 // Also need to pass the allele frequency data to the locus class

	 NumFreqLoci = AlleleFreq.NumberLoci;
	 for (i = 0; i < NumFreqLoci; i++)
	 {
		 Loci[i].NumberofAlleles = AlleleFreq.NumberAlleles[i];
		 for (j = 0; j < Loci[i].NumberofAlleles; j++)
		 {
			 Loci[i].AlleleName[j] = AlleleFreq.AlleleName[j][i];
			 Loci[i].AlleleFreq[j] = AlleleFreq.AlleleFreq[j][i];
		 }
		 
		 for (j = 0; j < static_cast<int>(AlleleFreq.LocusName[i].length()); j++)
			 Loci[i].Name[j] = AlleleFreq.LocusName[i][j];
		 Loci[i].Name[j] = '\0';
	 }
	 NumLociLoaded = NumFreqLoci;


	 // Copy the parameters in the parameter class into the old variables used by GerudSim


	 NumIterations = GSparam.NumberIterations;

	 for (i = 0; i < GSparam.NumberRuns; i++)
	 {
		 NumberEms[i] = GSparam.EmbryosAssayed[i];
		 for (j = 0; j < 10; j++)
			 NumFromMale[j][i] = GSparam.EmbryosFromMale[i][j];
		 
		 RealMaleNumber[i] = 10;
		 if (NumFromMale[9][i] == 0)
			 RealMaleNumber[i] = 9;
		 if (NumFromMale[8][i] == 0)
			 RealMaleNumber[i] = 8;
		 if (NumFromMale[7][i] == 0)
			 RealMaleNumber[i] = 7;
		 if (NumFromMale[6][i] == 0)
			 RealMaleNumber[i] = 6;
		 if (NumFromMale[5][i] == 0)
			 RealMaleNumber[i] = 5;
		 if (NumFromMale[4][i] == 0)
			 RealMaleNumber[i] = 4;
		 if (NumFromMale[3][i] == 0)
			 RealMaleNumber[i] = 3;
		 if (NumFromMale[2][i] == 0)
			 RealMaleNumber[i] = 2;
		 if (NumFromMale[1][i] == 0)
			 RealMaleNumber[i] = 1;
		 
	 }  // end of i loop

	 for (j = 0; j < 100; j++)
	 {
		 k = 0;
		 for (i = 0; i < RealMaleNumber[j]; i++)
			k = k + NumFromMale[i][j];

		 if (k < NumberEms[j])
			 NumberEms[j] = k;
	 }


	 if (GSparam.NumberLoci > NumLociLoaded)
		 GSparam.NumberLoci = NumLociLoaded;
	 NumberLoci = GSparam.NumberLoci;
	 NumLoci = GSparam.NumberLoci;
	 NumFreqLoci = GSparam.NumberLoci;
	 MomKnown = GSparam.KnownMother;

	 // This is the main loop, which goes through all the proposed runs and the reps of these runs

	 ErrorGP = false;
	 ErrorPP = false;
	 TotalNumberSolutionsError = false;
	 MAnumberError = false;
	 // seed the random number generator
	 sgenrand(static_cast<unsigned long int>(time(0)));

	 cout << "\nRunning...\n";

	 for (run = 0; run < GSparam.NumberRuns; run++)  // up to 100 runs with 100 different parameter combinations
	 {
		 
			CumCorrDadNum = 0;
			CumOverDadNum = 0;
			CumUnderDadNum = 0;
			CumNumDadsCorr = 0;
			CumNumDadsWrong = 0;
			CumNumFiveCorr = 0;
			CumNumFiveWrong = 0;
			CumNumTenCorr = 0;
			CumNumTenWrong = 0;
			CumNumTwentCorr = 0;
			CumNumTwentWrong = 0;
			CumNumMoreCorr = 0;
			CumNumMoreWrong = 0;
			CumNumNSampWrong = 0;
			CumNumNSampCorr = 0;
			CumExactlyRight = 0;
			CumGenosRight = 0;
			CumGenosWrong = 0;
			MinDadsOL = 0;
			MinDadsOLMore = 0;
			MinDadsOLLess = 0;
			NumOneSolnCorr = 0;
			NumTenSolnsCorr = 0;
			NumMoreSolnsCorr = 0;
			NumOneSolnWrong = 0;
			NumTenSolnsWrong = 0;
			NumMoreSolnsWrong = 0;
			NumNotSolved = 0;
			CumMotherCorrect = 0;		

			for (iteration = 0; iteration < NumIterations; iteration++)
			{
				MinDadsOneLocus = 10;
				BigDadsSolved = false;

				if (verbose && (iteration+1)%100 == 0)
					cout << "Run: " << (run + 1) << "/" << GSparam.NumberRuns << ", Iteration: " << (iteration + 1) << "/" << NumIterations << "...\n";

				//  Create Array Here -------------------------------------------------------

				Dad tempDads[10];
				for (i = 0; i < RealMaleNumber[run]; i++)
				{
					for (j = 0; j < NumberLoci; j++)
					{
						RealDad[i].Allele1[j] = 0;
						RealDad[i].Allele2[j] = 0;
						rndbl = genrand();
						dbsum = 0;
						for (k = 0; k < Loci[j].NumberofAlleles; k++)
						{
							dbsum = dbsum + Loci[j].AlleleFreq[k];
							if (rndbl < dbsum && RealDad[i].Allele1[j] == 0)
							{
								RealDad[i].Allele1[j] = Loci[j].AlleleName[k];
							}
						} // end of k loop

						rndbl = genrand();
						dbsum = 0;
						for (k = 0; k < Loci[j].NumberofAlleles; k++)
						{
							dbsum = dbsum + Loci[j].AlleleFreq[k];
							if (rndbl < dbsum && RealDad[i].Allele2[j] == 0)
							{
								RealDad[i].Allele2[j] = Loci[j].AlleleName[k];
							}
						} // end of k loop
					} // end of j loop
				} // end of i loop

				for (j = 0; j < NumberLoci; j++)
				{
					MomAllele1[j] = 0;
					MomAllele2[j] = 0;
					rndbl = genrand();
					dbsum = 0;
					for (k = 0; k < Loci[j].NumberofAlleles; k++)
					{
						dbsum = dbsum + Loci[j].AlleleFreq[k];
						if (rndbl < dbsum && MomAllele1[j] == 0)
						{
							MomAllele1[j] = Loci[j].AlleleName[k];
							RealMomAllele1[j] = MomAllele1[j];
						}
					} // end of k loop

					rndbl = genrand();
					dbsum = 0;
					for (k = 0; k < Loci[j].NumberofAlleles; k++)
					{
						dbsum = dbsum + Loci[j].AlleleFreq[k];
						if (rndbl < dbsum && MomAllele2[j] == 0)
						{
							MomAllele2[j] = Loci[j].AlleleName[k];
							RealMomAllele2[j] = MomAllele2[j];
						}
					} // end of k loop

				} // end of j loop

				NumberEmbryos = NumberEms[run];
				int tNumFromMale[10];
				int rndnum, sumnum, totalnum, coinflip;

				for (i = 0; i < 10; i++)
				{
					tNumFromMale[i] = NumFromMale[i][run];
					EmsSampFromMale[i] = 0;
				}

				for (i = 0; i < NumberEmbryos; i++)
				{
					EmbryoRealFather[i] = 99;
					totalnum = 0;
					sumnum = 0;
					for (j = 0; j < RealMaleNumber[run]; j++)
						totalnum = totalnum + tNumFromMale[j];

					rndnum = randnum(totalnum);

					for (j = 0; j < RealMaleNumber[run]; j++)
					{
						sumnum = sumnum + tNumFromMale[j];
						if (rndnum < sumnum && EmbryoRealFather[i] == 99)
						{
							EmbryoRealFather[i] = j;
							EmsSampFromMale[j]++;
							tNumFromMale[j] = tNumFromMale[j] - 1;
							for (k = 0; k < NumberLoci; k++)
							{
								coinflip = randnum(2);
								if (coinflip == 0)
									EmbryoAllele1[k][i] = MomAllele1[k];
								else
									EmbryoAllele1[k][i] = MomAllele2[k];
								coinflip = randnum(2);
								if (coinflip == 0)
									EmbryoAllele2[k][i] = RealDad[j].Allele1[k];
								else
									EmbryoAllele2[k][i] = RealDad[j].Allele2[k];
							} // end of k loop
						}   //end of if rndnum

					} // end of j loop

				} // end of i loop

				if (save_one_array)
				{
					ofstream oa_file;
					oa_file.open(OneArrayFilename);

					if (MomKnown)
						oa_file << "Known Mother\n";
					else
						oa_file << "Unknown Mother\n";

					oa_file << NumberEmbryos << " Embryos\n";
					oa_file << NumberLoci << " Loci\n";

					oa_file << "\t";
					for (j = 0; j < NumberLoci; j++)
						oa_file << Loci[j].Name << "\t";
					oa_file << "\n";

					if (MomKnown)
					{
						oa_file << "Mom\t";
						for (j = 0; j < NumberLoci; j++)
						{
							oa_file << MomAllele1[j] << "/" << MomAllele2[j] << "\t";
						}// end of j loop
						oa_file << "\n";
					}
							
					for (i = 0; i < NumberEmbryos; i++)
					{
						oa_file << "Emb_" << i+1 << "\t";
						for (j = 0; j < NumberLoci; j++)
						{
							oa_file << EmbryoAllele1[j][i] << "/" << EmbryoAllele2[j][i] << "\t";
						} // end of j loop
						oa_file << "\n";
					}  // end of i loop */

					oa_file << "\nReal Solution:\n";
					oa_file << "Mom\t";
					for (j = 0; j < NumberLoci; j++)
					{
						oa_file << MomAllele1[j] << "/" << MomAllele2[j] << "\t";
					}// end of j loop
					oa_file << "\n";

					for (i = 0; i < RealMaleNumber[run]; i++)
					{
						oa_file << "Dad " << i+1 << "\t";
						for (j = 0; j < NumberLoci; j++)
						{
							oa_file << RealDad[i].Allele1[j] << "/" << RealDad[i].Allele2[j] << "\t";
						}// end of j loop
						oa_file << "No. Emb:\t" << EmsSampFromMale[i];
						oa_file << "\n";
						}  // end of i loop
					oa_file.close();
					return 0; // Saves one array and then exits the program

				}

				//  End of Create Array ------------------------------------------------------


							


				// Do Maternal Alleles Routine --------------------------------------------------------------------

					// This loop will establish all possible maternal genotypes that are consistent
					// with the progeny array.

					bool MAflag;

					int n, o, p, q, r, s, t;
					int MAnumbermoms, MAnewnumber;

					int MAnumsols[50], m;
					int PMAllele1[50][1000];
					int PMAllele2[50][1000];


					for (j = 0; j < NumberLoci; j++)
					{
						MAalleles[j][0] = EmbryoAllele1[j][0];
						MAtotal[j] = 1;
					}

					for (i = 0; i < NumberEmbryos; i++)    // first make a list of all alleles
					{
						for (j = 0; j < NumberLoci; j++)
						{
							MAflag = true;
							for (k = 0; k < MAtotal[j]; k++)
							{
								if (EmbryoAllele1[j][i] == MAalleles[j][k])
								MAflag = false;
							} // end of k loop
							if (MAflag)
							{
								MAalleles[j][MAtotal[j]] = EmbryoAllele1[j][i];
								MAtotal[j]++;
							}
							MAflag = true;
							for (k = 0; k < MAtotal[j]; k++)
							{
								if (EmbryoAllele2[j][i] == MAalleles[j][k])
								MAflag = false;
							} // end of k loop
							if (MAflag)
							{
								MAalleles[j][MAtotal[j]] = EmbryoAllele2[j][i];
								MAtotal[j]++;
							}
						} // end of j loop
					} // end of i loop

					if (!MomKnown)
					{
						for (j = 0; j < 50; j++)
							MAnumsols[j] = 0;

						for (j = 0; j < NumLoci; j++)
						{
							for (i = 0; i < MAtotal[j]; i++)
							{
								for (k = i; k < MAtotal[j]; k++)
								{
									PMAllele1[j][MAnumsols[j]] = MAalleles[j][i];
									PMAllele2[j][MAnumsols[j]] = MAalleles[j][k];
									MAflag = true;
									for (m = 0; m < NumberEmbryos; m++)
									{
										if (EmbryoAllele1[j][m] != PMAllele1[j][MAnumsols[j]] &&
											EmbryoAllele1[j][m] != PMAllele2[j][MAnumsols[j]] &&
											EmbryoAllele2[j][m] != PMAllele1[j][MAnumsols[j]] &&
											EmbryoAllele2[j][m] != PMAllele2[j][MAnumsols[j]])
												MAflag = false;
									}  // end of m loop
									if (MAflag)
										MAnumsols[j]++;
								} // end of k loop
							} // end of i loop
						} // end of j loop

						MAnumbermoms = 0;
						MAnewnumber = 0;

						for (i = 0; i < MAnumsols[0]; i++)
						{
							Mothers[i].Allele1[0] = PMAllele1[0][i];
							Mothers[i].Allele2[0] = PMAllele2[0][i];
							MAnumbermoms++;
						} // end of i loop

						for (i = 0; i < MAnumbermoms; i++)
							PossMoms[i] = Mothers[i];


						for (j = 1; j < NumberLoci; j++)      // This loop creates all possible compatible genotypes
						{
							MAnewnumber = 0;
							for (i = 0; i < MAnumsols[j]; i++)
							{
								for (k = 0; k < MAnumbermoms; k++)
								{
									for (m = 0; m < j; m++)
									{
										Mothers[MAnewnumber].Allele1[m] = PossMoms[k].Allele1[m];
										Mothers[MAnewnumber].Allele2[m] = PossMoms[k].Allele2[m];
									} // end of m loop

									Mothers[MAnewnumber].Allele1[j] = PMAllele1[j][i];
									Mothers[MAnewnumber].Allele2[j] = PMAllele2[j][i];
									MAnewnumber++;
									if (MAnewnumber > MotherNumberMax)
									{
										MAnewnumber = MotherNumberMax;
										MAnumberError = true;
									}
								} // end of k loop
							} // end of i loop

							MAnumbermoms = MAnewnumber;
							for (i = 0; i < MAnewnumber; i++)
								PossMoms[i] = Mothers[i];

						} // end of j loop

					} // end of if (!MomKnown)

					if (MomKnown)
					{
						MAnewnumber = 1;
						MAnumbermoms = 1;
						for (j = 0; j < NumberLoci; j++)
						{
							Mothers[0].Allele1[j] = MomAllele1[j];
							Mothers[0].Allele2[j] = MomAllele2[j];
						}
					} // end of if MomKnown


					NumCompMoms = MAnumbermoms;
					


					/*for (i = 0; i < MAnumbermoms; i++)
					{
						cout << "Mom " << i+1 << '\t';
						for (j = 0; j < NumberLoci; j++)
							cout << Mothers[i].Allele1[j] << "/" << Mothers[i].Allele2[j] << '\t';
						cout << "\n";
					} */

				// End of Maternal Alleles Routine ------------------------------------------




				// Do Big Search ------------------------------------------------------------

					// Big Search has several parts, which are conducted in order.
					// First, it lists possible paternal alleles.  Second, it produces all
					// paternal genotypes that can explain at least one progeny.  Finally, it runs
					// the exhaustive Gerud algorithm.

					int BSi, BSj, BSk, BSl, BSm;
					int tempTNS;

					MaxNumDadsEx = 10;
					BestDadsPriority[0] = 0;
					TotalNumberSolutions = 0;
					tempTNS = 0;

					for (BSi = 0; BSi < NumCompMoms; BSi++)
					{

						for (BSj = 0; BSj < NumberLoci; BSj++)
						{
							MomAllele1[BSj] = Mothers[BSi].Allele1[BSj];
							MomAllele2[BSj] = Mothers[BSi].Allele2[BSj];
						}

					// Paternal Alleles ---------------------------------------------------------------------
					// --------------------------------------------------------------------------------------
						int ma1, ma2, ea1, ea2, i, j;
						for (i = 0; i < NumberEmbryos; i++)
						{
							for (j = 0; j < NumberLoci; j++)
							{
								ma1 = MomAllele1[j];
								ma2 = MomAllele2[j];
								ea1 = EmbryoAllele1[j][i];
								ea2 = EmbryoAllele2[j][i];

								if (compatible(ma1, ma2, ea1, ea2))
								{
									if (ea1 == ma1 || ea1 == ma2) EmbryoPaternal1[j][i] = ea2;
									if (ea2 == ma1 || ea2 == ma2) EmbryoPaternal1[j][i] = ea1;
									EmbryoPaternal2[j][i] = 0;
								}
								else
								{
									EmbryoPaternal1[j][i] = 9999;
									EmbryoPaternal2[j][i] = 9999;
								}

								if (identical(ma1, ma2, ea1, ea2) && ma1 != ma2)
								{
									EmbryoPaternal1[j][i] = ea1;
									EmbryoPaternal2[j][i] = ea2;
								}
							} // end of j loop
						} // end of i loop


						// This next section will tally up all of the paternal alleles
						int ep1, ep2, flag, k;

						for (j = 0; j < NumberLoci; j++)
						{
							PatAllList[j][0] = 9999;
							flag = 0;

							while (PatAllList[j][0] == 9999)
							{
								PatAllList[j][0] = EmbryoPaternal1[j][flag];
								flag++;
							}

							PatAllNum[j] = 1;

							for (i = 0; i < NumberEmbryos; i++)
							{
								ep1 = EmbryoPaternal1[j][i];
								ep2 = EmbryoPaternal2[j][i];

							if (ep1 != 0 && ep1 != 9999)
							{
								flag = 0;
								for (k = 0; k < PatAllNum[j]; k++)
								{
									if (ep1 == PatAllList[j][k])
										flag = 1;
								} // end of k loop
								if (flag == 0)
								{
									PatAllList[j][PatAllNum[j]] = ep1;
									PatAllNum[j]++;
								}
							}// end of if

							if (ep2 != 0 && ep2 != 9999)
							{
								flag = 0;
								for (k = 0; k < PatAllNum[j]; k++)
								{
									if (ep2 == PatAllList[j][k])
									flag = 1;
								} // end of k loop
								if (flag == 0)
								{
									PatAllList[j][PatAllNum[j]] = ep2;
									PatAllNum[j]++;
								}
							}// end of if
							} // end of i loop
						} // end of j loop

  
						for (i = 0; i < 1000; i++)
						{
							for (j = 0; j < 50; j++)
							{
								EmbryoPat1[j][i] = EmbryoPaternal1[j][i];
								EmbryoPat2[j][i] = EmbryoPaternal2[j][i];
							}
						}

					// Paternal Alleles Done ----------------------------------------------------------------
					// --------------------------------------------------------------------------------------



					// PaternGenos --------------------------------------------------------------------------
					// --------------------------------------------------------------------------------------

						bool ArrayErrGP;
						bool ArrayErrPP;
						bool tempGoodGen;

						ArrayErrGP = false;
						ArrayErrPP = false;

						long int NumPPGens;

						NumPPGens = 0;
						NumGPGens = 0;

						// The following loop comes up with all possible paternal genotypes at
						// the first locus, and the next loop keeps the genotypes that are
						// consistent with the progeny array (of course they all are for 1 locus)

						for (i = 0; i < PatAllNum[0]; i++)
						{
							for (j = i; j < PatAllNum[0]; j++)
							{
								PossDads[NumPPGens].Allele1[0] = PatAllList[0][i];
								PossDads[NumPPGens].Allele2[0] = PatAllList[0][j];
								NumPPGens++;
								if (NumPPGens > PPMax)
								{
									NumPPGens = PPMax;
									ArrayErrPP = true;
								}
							} // end of j loop
						} // end of i loop

						for (i = 0; i < NumPPGens; i++)
						{
							GoodDads[NumGPGens].Allele1[0] = PossDads[i].Allele1[0];
							GoodDads[NumGPGens].Allele2[0] = PossDads[i].Allele2[0];
							NumGPGens++;
							if (NumGPGens > GPMax)
							{
								NumGPGens = GPMax;
								ArrayErrGP = true;
							}
						}

						bool GoodGenotype;

						for (j = 1; j < NumberLoci; j++)
						{

							NumPPGens = 0;
							// The following series of loops generates all of the possible paternal
							// genotypes for the next locus given the good genotypes that were saved
							// for all previous loci

							for (k = 0; k < NumGPGens; k++)
							{
								for (l = 0; l < PatAllNum[j]; l++)
								{
									for (m = l; m < PatAllNum[j]; m++)
									{
										for (n = 0; n < j; n++)
										{
											PossDads[NumPPGens].Allele1[n] = GoodDads[k].Allele1[n];
											PossDads[NumPPGens].Allele2[n] = GoodDads[k].Allele2[n];
										} // end of n loop

										PossDads[NumPPGens].Allele1[j] = PatAllList[j][l];
										PossDads[NumPPGens].Allele2[j] = PatAllList[j][m];

										NumPPGens++;
										if (NumPPGens > PPMax)
										{
											NumPPGens = PPMax;
											ArrayErrPP = true;
										}
									}// end of m loop
								} // end of l loop
							} // end of k loop

							//  These loops will keep only those possible paternal genotypes that
							// are consistent with the progeny array

							NumGPGens = 0;

							for (k = 0; k < NumPPGens; k++)
							{
								GoodGenotype = false;
								PossDads[k].NumberOffspring = 0;
								for (l = 0; l < NumberEmbryos; l++)
								{
									tempGoodGen = true;
									for (m = 0; m < j+1; m++)
									{
										if (!compatible(PossDads[k].Allele1[m], PossDads[k].Allele2[m],
											EmbryoPat1[m][l], EmbryoPat2[m][l]))
												tempGoodGen = false;
									}

									if (tempGoodGen == true)
									{
										GoodGenotype = true;
										PossDads[k].NumberOffspring++; 
									}
								} // end of l loop

								if (GoodGenotype == true)
								{
									for (n = 0; n < j + 1; n++)
									{
										GoodDads[NumGPGens].Allele1[n] = PossDads[k].Allele1[n];
										GoodDads[NumGPGens].Allele2[n] = PossDads[k].Allele2[n];
									}// end of n loop
									GoodDads[NumGPGens].NumberOffspring = PossDads[k].NumberOffspring;
									NumGPGens++;
									if (NumGPGens > GPMax)
									{
										NumGPGens = GPMax;
										ArrayErrGP = true;
									}
								}
							} // end of k loop

							// Now get rid of redundancy caused by homozygotes that explain as many offspring as almost identical
							// heterozygotes.  For example, if 120/120  143/145 could be the father of 20 offspring,
							// then so could 120/122 143/145 or 120/124 143/145.  So we should only keep the homozygotes and
							// consider them 120/unknown allele (as long as more offspring are not explained by adding the second allele).

							for (i = 0; i < NumGPGens; i++)
							{
								for (m = 0; m < j+1; m++)
								{
									if (GoodDads[i].Allele1[m] == GoodDads[i].Allele2[m])
									{
										for (k = 0; k < NumGPGens; k++)
										{
											tempGoodGen = true;
											if (k != i && GoodDads[i].NumberOffspring == GoodDads[k].NumberOffspring)
											{
												if (GoodDads[i].Allele1[m] == GoodDads[k].Allele1[m] || GoodDads[i].Allele1[m] == GoodDads[k].Allele2[m])
												{
													tempGoodGen = false;
													for (l = 0; l < j + 1; l++)
													{
														if (l != m)
														{
															if (!identical(GoodDads[i].Allele1[l], GoodDads[i].Allele2[l], GoodDads[k].Allele1[l], GoodDads[k].Allele2[l])
																&& GoodDads[i].Allele1[l] != GoodDads[i].Allele2[l])
																	tempGoodGen = true;
															if (GoodDads[i].Allele1[l] == GoodDads[i].Allele2[l] &&
																!compatible(GoodDads[i].Allele1[l], GoodDads[i].Allele2[l], GoodDads[k].Allele1[l], GoodDads[k].Allele2[l]))
																	tempGoodGen = true;
														}
													} // end of l loop
												} // end of if statement
											} // end of if statement

											if (!tempGoodGen)
												GoodDads[k].NumberOffspring = 0;

										}// end of k loop
									} // end of first if
								}// end of m loop
							}// end of i loop

						} // end of the j loop

						NumPPGens = 0;
						for (i = 0; i < NumGPGens; i++)
						{
							if (GoodDads[i].NumberOffspring != 0)
							{
								PossDads[NumPPGens] = GoodDads[i];
								NumPPGens++;
								if (NumPPGens > PPMax)
								{
									NumPPGens = PPMax;
									ArrayErrPP = true;
								}
							}
						} // end of i loop

						NumGPGens = NumPPGens;

						for (i = 0; i < NumGPGens; i++)
							GoodDads[i] = PossDads[i];

						if (ArrayErrPP)
							ErrorPP = true;

						if (ArrayErrGP)
							ErrorGP = true;




					// PaternGenos Done ---------------------------------------------------------------------
					// --------------------------------------------------------------------------------------



					// Exhaustive ---------------------------------------------------------------------------
					// --------------------------------------------------------------------------------------

					// The next section will test the paternal genotypes to see which combinations of fathers
					// can explain the progeny array.

					
					bool CompatFather[10];
					NumWinDadError = false;

					DadsSolved = false;
					DadTestNum = 0;
					NumWinDads = 0;

					// Can 1 dad explain the data?
					for (i = 0; i < NumGPGens; i++)
					{
						if (GoodDads[i].NumberOffspring == NumberEmbryos)
						{
							DadsSolved = true;
							DadTestNum = 1;
							WinnerDad[NumWinDads][0] = i;
							NumWinDads++;
							if (NumWinDads > NumWinDadsMax)
							{
								NumWinDads = NumWinDadsMax;
								NumWinDadError = true;
							}
						}
					} // end of i

					// end of 1 dad section

					// Can 2 dad's explain the data?
					if (!DadsSolved && MaxNumDadsEx > 1)
					{
						for (i = 0; i < NumGPGens; i++)
						{
							for (j = i; j < NumGPGens; j++)
							{
								if (GoodDads[i].NumberOffspring + GoodDads[j].NumberOffspring > NumberEmbryos - 1)
								{
									tempGoodGen = true;
									for (m = 0; m < NumberEmbryos; m++)
									{
										CompatFather[0] = true;
										for (n = 0; n < NumberLoci; n++) 
										{
											if (!compatible(GoodDads[i].Allele1[n], GoodDads[i].Allele2[n],
												EmbryoPat1[n][m], EmbryoPat2[n][m]))
													CompatFather[0] = false;
										}

										CompatFather[1] = true;
										for (n = 0; n < NumberLoci; n++) 
										{
											if (!compatible(GoodDads[j].Allele1[n], GoodDads[j].Allele2[n],
												EmbryoPat1[n][m], EmbryoPat2[n][m]))
													CompatFather[1] = false;
										}

										if (!CompatFather[0] && !CompatFather[1])
										{
											tempGoodGen = false;
											m = NumberEmbryos;
										}

									}// end of m loop
									if (tempGoodGen == true)
									{
										DadsSolved = true;
										DadTestNum = 2;
										WinnerDad[NumWinDads][0] = i;
										WinnerDad[NumWinDads][1] = j;
										NumWinDads++;
										if (NumWinDads > NumWinDadsMax)
											{
												NumWinDads = NumWinDadsMax;
												NumWinDadError = true;
											}
									} // end of if
								} // end of if
							}// end of j loop
						}// end of i loop
					} // end of 2 dads (if statement) section

					// Can three dads explain the data?
					int o;
					if (!DadsSolved && MaxNumDadsEx > 2)
					{
						for (i = 0; i < NumGPGens; i++)
						{
							for (j = i; j < NumGPGens; j++)
							{
								for (o = j; o < NumGPGens; o++)
								{
									if (GoodDads[i].NumberOffspring + GoodDads[j].NumberOffspring
										+ GoodDads[o].NumberOffspring > NumberEmbryos - 1)
									{
										tempGoodGen = true;
										for (m = 0; m < NumberEmbryos; m++)
										{
											CompatFather[0] = true;
											for (n = 0; n < NumberLoci; n++) 
											{
												if (!compatible(GoodDads[i].Allele1[n], GoodDads[i].Allele2[n],
													EmbryoPat1[n][m], EmbryoPat2[n][m]))
														CompatFather[0] = false;
											}

											CompatFather[1] = true;
											for (n = 0; n < NumberLoci; n++) 
											{
												if (!compatible(GoodDads[j].Allele1[n], GoodDads[j].Allele2[n],
													EmbryoPat1[n][m], EmbryoPat2[n][m]))
														CompatFather[1] = false;
											}

											CompatFather[2] = true;
											for (n = 0; n < NumberLoci; n++) 
											{
												if (!compatible(GoodDads[o].Allele1[n], GoodDads[o].Allele2[n],
													EmbryoPat1[n][m], EmbryoPat2[n][m]))
														CompatFather[2] = false;
											}

											if (!CompatFather[0] && !CompatFather[1] && !CompatFather[2])
											{
												tempGoodGen = false;
												m = NumberEmbryos;
											}

										}// end of m loop
										if (tempGoodGen == true)
										{
											DadsSolved = true;
											DadTestNum = 3;
											WinnerDad[NumWinDads][0] = i;
											WinnerDad[NumWinDads][1] = j;
											WinnerDad[NumWinDads][2] = o;
											NumWinDads++;
											if (NumWinDads > NumWinDadsMax)
											{
												NumWinDads = NumWinDadsMax;
												NumWinDadError = true;
											}
										} // end of if
									} // end of if
								} // end of o loop
							}// end of j loop
						}// end of i loop
					} // end of 3 dads (if statement) section

					// Can four dads explain the data?
					int p;
					if (!DadsSolved && MaxNumDadsEx > 3)
					{
						for (i = 0; i < NumGPGens; i++)
						{
							for (j = i; j < NumGPGens; j++)
							{
								for (o = j; o < NumGPGens; o++)
								{
									for (p = o; p < NumGPGens; p++)
									{
										if (GoodDads[i].NumberOffspring + GoodDads[j].NumberOffspring
											+ GoodDads[o].NumberOffspring + GoodDads[p].NumberOffspring > NumberEmbryos - 1)
										{
											tempGoodGen = true;
											for (m = 0; m < NumberEmbryos; m++)
											{
												CompatFather[0] = true;
												for (n = 0; n < NumberLoci; n++) 
												{
													if (!compatible(GoodDads[i].Allele1[n], GoodDads[i].Allele2[n],
														EmbryoPat1[n][m], EmbryoPat2[n][m]))
															CompatFather[0] = false;
												}

												CompatFather[1] = true;
												for (n = 0; n < NumberLoci; n++) 
												{
													if (!compatible(GoodDads[j].Allele1[n], GoodDads[j].Allele2[n],
														EmbryoPat1[n][m], EmbryoPat2[n][m]))
															CompatFather[1] = false;
												}

												CompatFather[2] = true;
												for (n = 0; n < NumberLoci; n++) 
												{
													if (!compatible(GoodDads[o].Allele1[n], GoodDads[o].Allele2[n],
														EmbryoPat1[n][m], EmbryoPat2[n][m]))
															CompatFather[2] = false;
												}

												CompatFather[3] = true;
												for (n = 0; n < NumberLoci; n++) 
												{
													if (!compatible(GoodDads[p].Allele1[n], GoodDads[p].Allele2[n],
														EmbryoPat1[n][m], EmbryoPat2[n][m]))
															CompatFather[3] = false;
												}

												if (!CompatFather[0] && !CompatFather[1] && !CompatFather[2] && !CompatFather[3])
												{
													tempGoodGen = false;
													m = NumberEmbryos;
												}

											}// end of m loop
											if (tempGoodGen == true)
											{
												DadsSolved = true;
												DadTestNum = 4;
												WinnerDad[NumWinDads][0] = i;
												WinnerDad[NumWinDads][1] = j;
												WinnerDad[NumWinDads][2] = o;
												WinnerDad[NumWinDads][3] = p;
												NumWinDads++;
												if (NumWinDads > NumWinDadsMax)
												{
													NumWinDads = NumWinDadsMax;
													NumWinDadError = true;
												}
											} // end of if
										} // end of if
									}  // end of p loop
								} // end of o loop
							}// end of j loop
						}// end of i loop
					} // end of 4 dads (if statement) section

					// Can five dads explain the data?
					int q;
					if (!DadsSolved && MaxNumDadsEx > 4)
					{
						for (i = 0; i < NumGPGens; i++)
						{
							for (j = i; j < NumGPGens; j++)
							{
								for (o = j; o < NumGPGens; o++)
								{
									for (p = o; p < NumGPGens; p++)
									{
										for (q = p; q < NumGPGens; q++)
										{
											if (GoodDads[i].NumberOffspring + GoodDads[j].NumberOffspring
												+ GoodDads[o].NumberOffspring + GoodDads[p].NumberOffspring
												+ GoodDads[q].NumberOffspring > NumberEmbryos - 1)
											{
												tempGoodGen = true;
												for (m = 0; m < NumberEmbryos; m++)
												{
													CompatFather[0] = true;
													for (n = 0; n < NumberLoci; n++) 
													{
														if (!compatible(GoodDads[i].Allele1[n], GoodDads[i].Allele2[n],
															EmbryoPat1[n][m], EmbryoPat2[n][m]))
																CompatFather[0] = false;
													}

													CompatFather[1] = true;
													for (n = 0; n < NumberLoci; n++) 
													{
														if (!compatible(GoodDads[j].Allele1[n], GoodDads[j].Allele2[n],
															EmbryoPat1[n][m], EmbryoPat2[n][m]))
																CompatFather[1] = false;
													}

													CompatFather[2] = true;
													for (n = 0; n < NumberLoci; n++) 
													{
														if (!compatible(GoodDads[o].Allele1[n], GoodDads[o].Allele2[n],
															EmbryoPat1[n][m], EmbryoPat2[n][m]))
																CompatFather[2] = false;
													}

													CompatFather[3] = true;
													for (n = 0; n < NumberLoci; n++) 
													{
														if (!compatible(GoodDads[p].Allele1[n], GoodDads[p].Allele2[n],
															EmbryoPat1[n][m], EmbryoPat2[n][m]))
																CompatFather[3] = false;
													}

													CompatFather[4] = true;
													for (n = 0; n < NumberLoci; n++) 
													{
														if (!compatible(GoodDads[q].Allele1[n], GoodDads[q].Allele2[n],
															EmbryoPat1[n][m], EmbryoPat2[n][m]))
																CompatFather[4] = false;
													}

													if (!CompatFather[0] && !CompatFather[1] && !CompatFather[2] && !CompatFather[3] && !CompatFather[4])
													{
														tempGoodGen = false;
														m = NumberEmbryos;
													}

												}// end of m loop
												if (tempGoodGen == true)
												{
													DadsSolved = true;
													DadTestNum = 5;
													WinnerDad[NumWinDads][0] = i;
													WinnerDad[NumWinDads][1] = j;
													WinnerDad[NumWinDads][2] = o;
													WinnerDad[NumWinDads][3] = p;
													WinnerDad[NumWinDads][4] = q;
													NumWinDads++;
													if (NumWinDads > NumWinDadsMax)
													{
														NumWinDads = NumWinDadsMax;
														NumWinDadError = true;
													}
												} // end of if
											} // end of if
										}  // end of q loop
									}  // end of p loop
								} // end of o loop
							}// end of j loop
						}// end of i loop
					} // end of 5 dads (if statement) section

					// Can six dads explain the data?
					int r;
					if (!DadsSolved && MaxNumDadsEx > 5)
					{
						for (i = 0; i < NumGPGens; i++)
						{
							for (j = i; j < NumGPGens; j++)
							{
								for (o = j; o < NumGPGens; o++)
								{
									for (p = o; p < NumGPGens; p++)
									{
										for (q = p; q < NumGPGens; q++)
										{
											for (r = q; r < NumGPGens; r++)
											{
												if (GoodDads[i].NumberOffspring + GoodDads[j].NumberOffspring
													+ GoodDads[o].NumberOffspring + GoodDads[p].NumberOffspring
													+ GoodDads[q].NumberOffspring + GoodDads[r].NumberOffspring > NumberEmbryos - 1)
												{
													tempGoodGen = true;
													for (m = 0; m < NumberEmbryos; m++)
													{
														CompatFather[0] = true;
														for (n = 0; n < NumberLoci; n++) 
														{
															if (!compatible(GoodDads[i].Allele1[n], GoodDads[i].Allele2[n],
																EmbryoPat1[n][m], EmbryoPat2[n][m]))
																	CompatFather[0] = false;
														}

														CompatFather[1] = true;
														for (n = 0; n < NumberLoci; n++) 
														{
															if (!compatible(GoodDads[j].Allele1[n], GoodDads[j].Allele2[n],
																EmbryoPat1[n][m], EmbryoPat2[n][m]))
																	CompatFather[1] = false;
														}

														CompatFather[2] = true;
														for (n = 0; n < NumberLoci; n++) 
														{
															if (!compatible(GoodDads[o].Allele1[n], GoodDads[o].Allele2[n],
																EmbryoPat1[n][m], EmbryoPat2[n][m]))
																	CompatFather[2] = false;
														}

														CompatFather[3] = true;
														for (n = 0; n < NumberLoci; n++) 
														{
															if (!compatible(GoodDads[p].Allele1[n], GoodDads[p].Allele2[n],
																EmbryoPat1[n][m], EmbryoPat2[n][m]))
																	CompatFather[3] = false;
														}

														CompatFather[4] = true;
														for (n = 0; n < NumberLoci; n++) 
														{
															if (!compatible(GoodDads[q].Allele1[n], GoodDads[q].Allele2[n],
																EmbryoPat1[n][m], EmbryoPat2[n][m]))
																	CompatFather[4] = false;
														}

														CompatFather[5] = true;
														for (n = 0; n < NumberLoci; n++) 
														{
															if (!compatible(GoodDads[r].Allele1[n], GoodDads[r].Allele2[n],
																EmbryoPat1[n][m], EmbryoPat2[n][m]))
																	CompatFather[5] = false;
														}

														if (!CompatFather[0] && !CompatFather[1] && !CompatFather[2] && !CompatFather[3] && !CompatFather[4] && !CompatFather[5])
														{
															tempGoodGen = false;
															m = NumberEmbryos;
														}
                       
													}// end of m loop
													if (tempGoodGen == true)
													{
														DadsSolved = true;
														DadTestNum = 6;
														WinnerDad[NumWinDads][0] = i;
														WinnerDad[NumWinDads][1] = j;
														WinnerDad[NumWinDads][2] = o;
														WinnerDad[NumWinDads][3] = p;
														WinnerDad[NumWinDads][4] = q;
														WinnerDad[NumWinDads][5] = r;
														NumWinDads++;
														if (NumWinDads > NumWinDadsMax)
														{
															NumWinDads = NumWinDadsMax;
															NumWinDadError = true;
														}
													} // end of if
												} // end of if
											}  // end of r loop
										}  // end of q loop
									}  // end of p loop
								} // end of o loop
							}// end of j loop
						}// end of i loop
					} // end of 6 dads (if statement) section

					// Can seven dads explain the data?
					int r2;
					if (!DadsSolved && MaxNumDadsEx > 6)
					{
						for (i = 0; i < NumGPGens; i++)
						{
							for (j = i; j < NumGPGens; j++)
							{
								for (o = j; o < NumGPGens; o++)
								{
									for (p = o; p < NumGPGens; p++)
									{
										for (q = p; q < NumGPGens; q++)
										{
											for (r = q; r < NumGPGens; r++)
											{
											for (r2 = r; r2 < NumGPGens; r2++)
												{
													if (GoodDads[i].NumberOffspring + GoodDads[j].NumberOffspring
														+ GoodDads[o].NumberOffspring + GoodDads[p].NumberOffspring
														+ GoodDads[q].NumberOffspring + GoodDads[r].NumberOffspring 
														+ GoodDads[r2].NumberOffspring > NumberEmbryos - 1)
													{
														tempGoodGen = true;
														for (m = 0; m < NumberEmbryos; m++)
														{
															CompatFather[0] = true;
															for (n = 0; n < NumberLoci; n++) 
															{
																if (!compatible(GoodDads[i].Allele1[n], GoodDads[i].Allele2[n],
																	EmbryoPat1[n][m], EmbryoPat2[n][m]))
																		CompatFather[0] = false;
															}

															CompatFather[1] = true;
															for (n = 0; n < NumberLoci; n++) 
															{
																if (!compatible(GoodDads[j].Allele1[n], GoodDads[j].Allele2[n],
																	EmbryoPat1[n][m], EmbryoPat2[n][m]))
																		CompatFather[1] = false;
															}

															CompatFather[2] = true;
															for (n = 0; n < NumberLoci; n++) 
															{
																if (!compatible(GoodDads[o].Allele1[n], GoodDads[o].Allele2[n],
																	EmbryoPat1[n][m], EmbryoPat2[n][m]))
																		CompatFather[2] = false;
															}

															CompatFather[3] = true;
															for (n = 0; n < NumberLoci; n++) 
															{
																if (!compatible(GoodDads[p].Allele1[n], GoodDads[p].Allele2[n],
																	EmbryoPat1[n][m], EmbryoPat2[n][m]))
																		CompatFather[3] = false;
															}

															CompatFather[4] = true;
															for (n = 0; n < NumberLoci; n++) 
															{
																if (!compatible(GoodDads[q].Allele1[n], GoodDads[q].Allele2[n],
																	EmbryoPat1[n][m], EmbryoPat2[n][m]))
																		CompatFather[4] = false;
															}

															CompatFather[5] = true;
															for (n = 0; n < NumberLoci; n++) 
															{
																if (!compatible(GoodDads[r].Allele1[n], GoodDads[r].Allele2[n],
																	EmbryoPat1[n][m], EmbryoPat2[n][m]))
																		CompatFather[5] = false;
															}

															CompatFather[6] = true;
															for (n = 0; n < NumberLoci; n++) 
															{
																if (!compatible(GoodDads[r2].Allele1[n], GoodDads[r2].Allele2[n],
																	EmbryoPat1[n][m], EmbryoPat2[n][m]))
																		CompatFather[6] = false;
															}

															if (!CompatFather[0] && !CompatFather[1] && !CompatFather[2] && !CompatFather[3] && !CompatFather[4] && !CompatFather[5] && !CompatFather[6])
															{
																tempGoodGen = false;
																m = NumberEmbryos;
															}
                       
														}// end of m loop
														if (tempGoodGen == true)
														{
															DadsSolved = true;
															DadTestNum = 7;
															WinnerDad[NumWinDads][0] = i;
															WinnerDad[NumWinDads][1] = j;
															WinnerDad[NumWinDads][2] = o;
															WinnerDad[NumWinDads][3] = p;
															WinnerDad[NumWinDads][4] = q;
															WinnerDad[NumWinDads][5] = r;
															WinnerDad[NumWinDads][6] = r2;
															NumWinDads++;
															if (NumWinDads > NumWinDadsMax)
															{
																NumWinDads = NumWinDadsMax;
																NumWinDadError = true;
															}
														} // end of if
													} // end of if
												}  // end of r2 loop
											}  // end of r loop
										}  // end of q loop
									}  // end of p loop
								} // end of o loop
							}// end of j loop
						}// end of i loop
					} // end of 7 dads (if statement) section
					// end of 7 dads section

					// Can 8 dads explain the data?
					int r3;
					if (!DadsSolved && MaxNumDadsEx > 7)
					{
						for (i = 0; i < NumGPGens; i++)
						{
							for (j = i; j < NumGPGens; j++)
							{
								for (o = j; o < NumGPGens; o++)
								{
									for (p = o; p < NumGPGens; p++)
									{
										for (q = p; q < NumGPGens; q++)
										{
											for (r = q; r < NumGPGens; r++)
											{
												for (r2 = r; r2 < NumGPGens; r2++)
												{
													for(r3 = r2; r3 < NumGPGens; r3++)
													{

														if (GoodDads[i].NumberOffspring + GoodDads[j].NumberOffspring
															+ GoodDads[o].NumberOffspring + GoodDads[p].NumberOffspring
															+ GoodDads[q].NumberOffspring + GoodDads[r].NumberOffspring 
															+ GoodDads[r2].NumberOffspring + GoodDads[r3].NumberOffspring > NumberEmbryos - 1)
														{
															tempGoodGen = true;
															for (m = 0; m < NumberEmbryos; m++)
															{
																CompatFather[0] = true;
																for (n = 0; n < NumberLoci; n++) 
																{
																	if (!compatible(GoodDads[i].Allele1[n], GoodDads[i].Allele2[n],
																		EmbryoPat1[n][m], EmbryoPat2[n][m]))
																			CompatFather[0] = false;
																}

																CompatFather[1] = true;
																for (n = 0; n < NumberLoci; n++) 
																{
																	if (!compatible(GoodDads[j].Allele1[n], GoodDads[j].Allele2[n],
																		EmbryoPat1[n][m], EmbryoPat2[n][m]))
																			CompatFather[1] = false;
																}

																CompatFather[2] = true;
																for (n = 0; n < NumberLoci; n++) 
																{
																	if (!compatible(GoodDads[o].Allele1[n], GoodDads[o].Allele2[n],
																		EmbryoPat1[n][m], EmbryoPat2[n][m]))
																			CompatFather[2] = false;
																}

																CompatFather[3] = true;
																for (n = 0; n < NumberLoci; n++) 
																{
																	if (!compatible(GoodDads[p].Allele1[n], GoodDads[p].Allele2[n],
																		EmbryoPat1[n][m], EmbryoPat2[n][m]))
																			CompatFather[3] = false;
																}

																CompatFather[4] = true;
																for (n = 0; n < NumberLoci; n++) 
																{
																	if (!compatible(GoodDads[q].Allele1[n], GoodDads[q].Allele2[n],
																		EmbryoPat1[n][m], EmbryoPat2[n][m]))
																			CompatFather[4] = false;
																}

																CompatFather[5] = true;
																for (n = 0; n < NumberLoci; n++) 
																{
																	if (!compatible(GoodDads[r].Allele1[n], GoodDads[r].Allele2[n],
																		EmbryoPat1[n][m], EmbryoPat2[n][m]))
																			CompatFather[5] = false;
																}

																CompatFather[6] = true;
																for (n = 0; n < NumberLoci; n++) 
																{
																	if (!compatible(GoodDads[r2].Allele1[n], GoodDads[r2].Allele2[n],
																		EmbryoPat1[n][m], EmbryoPat2[n][m]))
																			CompatFather[6] = false;
																}

																CompatFather[7] = true;
																for (n = 0; n < NumberLoci; n++) 
																{
																	if (!compatible(GoodDads[r3].Allele1[n], GoodDads[r3].Allele2[n],
																		EmbryoPat1[n][m], EmbryoPat2[n][m]))
																			CompatFather[7] = false;
																}

																if (!CompatFather[0] && !CompatFather[1] && !CompatFather[2] && !CompatFather[3] && !CompatFather[4] && !CompatFather[5] 
																		&& !CompatFather[6] && !CompatFather[7])
																{
																	tempGoodGen = false;
																	m = NumberEmbryos;
																}
															}// end of m loop
															
															if (tempGoodGen == true)
															{
																DadsSolved = true;
																DadTestNum = 8;
																WinnerDad[NumWinDads][0] = i;
																WinnerDad[NumWinDads][1] = j;
																WinnerDad[NumWinDads][2] = o;
																WinnerDad[NumWinDads][3] = p;
																WinnerDad[NumWinDads][4] = q;
																WinnerDad[NumWinDads][5] = r;
																WinnerDad[NumWinDads][6] = r2;
																WinnerDad[NumWinDads][7] = r3;
																NumWinDads++;
																if (NumWinDads > NumWinDadsMax)
																{
																	NumWinDads = NumWinDadsMax;
																	NumWinDadError = true;
																}
															} // end of if
														} // end of if
													} // end of r3 loop
												}  // end of r2 loop
											}  // end of r loop
										}  // end of q loop
									}  // end of p loop
								} // end of o loop
							}// end of j loop
						}// end of i loop
					} // end of 8 dads (if statement) section
					
					// end of 8 dads section

					// Can 9 dads explain the data?
					int r4;
					if (!DadsSolved && MaxNumDadsEx > 8)
					{
						for (i = 0; i < NumGPGens; i++)
						{
							for (j = i; j < NumGPGens; j++)
							{
								for (o = j; o < NumGPGens; o++)
								{
									for (p = o; p < NumGPGens; p++)
									{
										for (q = p; q < NumGPGens; q++)
										{
											for (r = q; r < NumGPGens; r++)
											{
												for (r2 = r; r2 < NumGPGens; r2++)
												{
													for(r3 = r2; r3 < NumGPGens; r3++)
													{
														for (r4 = r3; r4 < NumGPGens; r4++)
														{

															if (GoodDads[i].NumberOffspring + GoodDads[j].NumberOffspring
																+ GoodDads[o].NumberOffspring + GoodDads[p].NumberOffspring
																+ GoodDads[q].NumberOffspring + GoodDads[r].NumberOffspring 
																+ GoodDads[r2].NumberOffspring + GoodDads[r3].NumberOffspring 
																+ GoodDads[r4].NumberOffspring > NumberEmbryos - 1)
															{
																tempGoodGen = true;
																for (m = 0; m < NumberEmbryos; m++)
																{
																	CompatFather[0] = true;
																	for (n = 0; n < NumberLoci; n++) 
																	{
																		if (!compatible(GoodDads[i].Allele1[n], GoodDads[i].Allele2[n],
																			EmbryoPat1[n][m], EmbryoPat2[n][m]))
																				CompatFather[0] = false;
																	}

																	CompatFather[1] = true;
																	for (n = 0; n < NumberLoci; n++) 
																	{
																		if (!compatible(GoodDads[j].Allele1[n], GoodDads[j].Allele2[n],
																			EmbryoPat1[n][m], EmbryoPat2[n][m]))
																				CompatFather[1] = false;
																	}

																	CompatFather[2] = true;
																	for (n = 0; n < NumberLoci; n++) 
																	{
																		if (!compatible(GoodDads[o].Allele1[n], GoodDads[o].Allele2[n],
																			EmbryoPat1[n][m], EmbryoPat2[n][m]))
																				CompatFather[2] = false;
																	}

																	CompatFather[3] = true;
																	for (n = 0; n < NumberLoci; n++) 
																	{
																		if (!compatible(GoodDads[p].Allele1[n], GoodDads[p].Allele2[n],
																			EmbryoPat1[n][m], EmbryoPat2[n][m]))
																				CompatFather[3] = false;
																	}

																	CompatFather[4] = true;
																	for (n = 0; n < NumberLoci; n++) 
																	{
																		if (!compatible(GoodDads[q].Allele1[n], GoodDads[q].Allele2[n],
																			EmbryoPat1[n][m], EmbryoPat2[n][m]))
																				CompatFather[4] = false;
																	}

																	CompatFather[5] = true;
																	for (n = 0; n < NumberLoci; n++) 
																	{
																		if (!compatible(GoodDads[r].Allele1[n], GoodDads[r].Allele2[n],
																			EmbryoPat1[n][m], EmbryoPat2[n][m]))
																				CompatFather[5] = false;
																	}

																	CompatFather[6] = true;
																	for (n = 0; n < NumberLoci; n++) 
																	{
																		if (!compatible(GoodDads[r2].Allele1[n], GoodDads[r2].Allele2[n],
																			EmbryoPat1[n][m], EmbryoPat2[n][m]))
																				CompatFather[6] = false;
																	}

																	CompatFather[7] = true;
																	for (n = 0; n < NumberLoci; n++) 
																	{
																		if (!compatible(GoodDads[r3].Allele1[n], GoodDads[r3].Allele2[n],
																			EmbryoPat1[n][m], EmbryoPat2[n][m]))
																				CompatFather[7] = false;
																	}

																	CompatFather[8] = true;
																	for (n = 0; n < NumberLoci; n++) 
																	{
																		if (!compatible(GoodDads[r4].Allele1[n], GoodDads[r4].Allele2[n],
																			EmbryoPat1[n][m], EmbryoPat2[n][m]))
																				CompatFather[8] = false;
																	}


																	if (!CompatFather[0] && !CompatFather[1] && !CompatFather[2] && !CompatFather[3] && !CompatFather[4] && !CompatFather[5] 
																			&& !CompatFather[6] && !CompatFather[7] & !CompatFather[8])
																	{
																		tempGoodGen = false;
																		m = NumberEmbryos;
																	}
																}// end of m loop
															
																if (tempGoodGen == true)
																{
																	DadsSolved = true;
																	DadTestNum = 9;
																	WinnerDad[NumWinDads][0] = i;
																	WinnerDad[NumWinDads][1] = j;
																	WinnerDad[NumWinDads][2] = o;
																	WinnerDad[NumWinDads][3] = p;
																	WinnerDad[NumWinDads][4] = q;
																	WinnerDad[NumWinDads][5] = r;
																	WinnerDad[NumWinDads][6] = r2;
																	WinnerDad[NumWinDads][7] = r3;
																	WinnerDad[NumWinDads][8] = r4;
																	NumWinDads++;
																	if (NumWinDads > NumWinDadsMax)
																	{
																		NumWinDads = NumWinDadsMax;
																		NumWinDadError = true;
																	}
																} // end of if
															} // end of if
														} // end of r4 loop
													} // end of r3 loop
												}  // end of r2 loop
											}  // end of r loop
										}  // end of q loop
									}  // end of p loop
								} // end of o loop
							}// end of j loop
						}// end of i loop
					} // end of 9 dads (if statement) section
					
					// end of 9 dads section

					// Can 10 dads explain the data?
					int r5;
					if (!DadsSolved && MaxNumDadsEx > 9)
					{
						for (i = 0; i < NumGPGens; i++)
						{
							for (j = i; j < NumGPGens; j++)
							{
								for (o = j; o < NumGPGens; o++)
								{
									for (p = o; p < NumGPGens; p++)
									{
										for (q = p; q < NumGPGens; q++)
										{
											for (r = q; r < NumGPGens; r++)
											{
												for (r2 = r; r2 < NumGPGens; r2++)
												{
													for(r3 = r2; r3 < NumGPGens; r3++)
													{
														for (r4 = r3; r4 < NumGPGens; r4++)
														{
															for (r5 = r4; r5 < NumGPGens; r5++)
															{

																if (GoodDads[i].NumberOffspring + GoodDads[j].NumberOffspring
																	+ GoodDads[o].NumberOffspring + GoodDads[p].NumberOffspring
																	+ GoodDads[q].NumberOffspring + GoodDads[r].NumberOffspring 
																	+ GoodDads[r2].NumberOffspring + GoodDads[r3].NumberOffspring 
																	+ GoodDads[r4].NumberOffspring + GoodDads[r5].NumberOffspring > NumberEmbryos - 1)
																{
																	tempGoodGen = true;
																	for (m = 0; m < NumberEmbryos; m++)
																	{
																		CompatFather[0] = true;
																		for (n = 0; n < NumberLoci; n++) 
																		{
																			if (!compatible(GoodDads[i].Allele1[n], GoodDads[i].Allele2[n],
																				EmbryoPat1[n][m], EmbryoPat2[n][m]))
																					CompatFather[0] = false;
																		}

																		CompatFather[1] = true;
																		for (n = 0; n < NumberLoci; n++) 
																		{
																			if (!compatible(GoodDads[j].Allele1[n], GoodDads[j].Allele2[n],
																				EmbryoPat1[n][m], EmbryoPat2[n][m]))
																					CompatFather[1] = false;
																		}

																		CompatFather[2] = true;
																		for (n = 0; n < NumberLoci; n++) 
																		{
																			if (!compatible(GoodDads[o].Allele1[n], GoodDads[o].Allele2[n],
																				EmbryoPat1[n][m], EmbryoPat2[n][m]))
																					CompatFather[2] = false;
																		}

																		CompatFather[3] = true;
																		for (n = 0; n < NumberLoci; n++) 
																		{
																			if (!compatible(GoodDads[p].Allele1[n], GoodDads[p].Allele2[n],
																				EmbryoPat1[n][m], EmbryoPat2[n][m]))
																					CompatFather[3] = false;
																		}

																		CompatFather[4] = true;
																		for (n = 0; n < NumberLoci; n++) 
																		{
																			if (!compatible(GoodDads[q].Allele1[n], GoodDads[q].Allele2[n],
																				EmbryoPat1[n][m], EmbryoPat2[n][m]))
																					CompatFather[4] = false;
																		}

																		CompatFather[5] = true;
																		for (n = 0; n < NumberLoci; n++) 
																		{
																			if (!compatible(GoodDads[r].Allele1[n], GoodDads[r].Allele2[n],
																				EmbryoPat1[n][m], EmbryoPat2[n][m]))
																					CompatFather[5] = false;
																		}

																		CompatFather[6] = true;
																		for (n = 0; n < NumberLoci; n++) 
																		{
																			if (!compatible(GoodDads[r2].Allele1[n], GoodDads[r2].Allele2[n],
																				EmbryoPat1[n][m], EmbryoPat2[n][m]))
																					CompatFather[6] = false;
																		}

																		CompatFather[7] = true;
																		for (n = 0; n < NumberLoci; n++) 
																		{
																			if (!compatible(GoodDads[r3].Allele1[n], GoodDads[r3].Allele2[n],
																				EmbryoPat1[n][m], EmbryoPat2[n][m]))
																					CompatFather[7] = false;
																		}

																		CompatFather[8] = true;
																		for (n = 0; n < NumberLoci; n++) 
																		{
																			if (!compatible(GoodDads[r4].Allele1[n], GoodDads[r4].Allele2[n],
																				EmbryoPat1[n][m], EmbryoPat2[n][m]))
																					CompatFather[8] = false;
																		}

																		CompatFather[9] = true;
																		for (n = 0; n < NumberLoci; n++) 
																		{
																			if (!compatible(GoodDads[r5].Allele1[n], GoodDads[r5].Allele2[n],
																				EmbryoPat1[n][m], EmbryoPat2[n][m]))
																					CompatFather[9] = false;
																		}


																		if (!CompatFather[0] && !CompatFather[1] && !CompatFather[2] && !CompatFather[3] && !CompatFather[4] && !CompatFather[5] 
																				&& !CompatFather[6] && !CompatFather[7] & !CompatFather[8] && !CompatFather[9])
																		{
																			tempGoodGen = false;
																			m = NumberEmbryos;
																		}
																	}// end of m loop
															
																	if (tempGoodGen == true)
																	{
																		DadsSolved = true;
																		DadTestNum = 10;
																		WinnerDad[NumWinDads][0] = i;
																		WinnerDad[NumWinDads][1] = j;
																		WinnerDad[NumWinDads][2] = o;
																		WinnerDad[NumWinDads][3] = p;
																		WinnerDad[NumWinDads][4] = q;
																		WinnerDad[NumWinDads][5] = r;
																		WinnerDad[NumWinDads][6] = r2;
																		WinnerDad[NumWinDads][7] = r3;
																		WinnerDad[NumWinDads][8] = r4;
																		WinnerDad[NumWinDads][9] = r5;
																		NumWinDads++;
																		if (NumWinDads > NumWinDadsMax)
																		{
																			NumWinDads = NumWinDadsMax;
																			NumWinDadError = true;
																		}
																	} // end of if
																} // end of if
															} // end of r5 loop
														} // end of r4 loop
													} // end of r3 loop
												}  // end of r2 loop
											}  // end of r loop
										}  // end of q loop
									}  // end of p loop
								} // end of o loop
							}// end of j loop
						}// end of i loop
					} // end of 10 dads (if statement) section
					
					// end of 10 dads section



					if (NumWinDads > 0 && DadTestNum > 0 && DadTestNum < MaxNumDadsEx)
						MaxNumDadsEx = DadTestNum;

					for (i = 0; i < NumWinDads; i++)
					{
						Solutions[TotalNumberSolutions].NumberFathers = DadTestNum;

						for (j = 0; j < NumberLoci; j++)
						{
							Solutions[TotalNumberSolutions].MoAll1[j] = MomAllele1[j];
							Solutions[TotalNumberSolutions].MoAll2[j] = MomAllele2[j];
						} // end of j

						for (j = 0; j < DadTestNum; j++)
						{
							for (l = 0; l < NumberLoci; l++)
							{
								Solutions[TotalNumberSolutions].FaAll1[l][j] = GoodDads[WinnerDad[i][j]].Allele1[l];
								Solutions[TotalNumberSolutions].FaAll2[l][j] = GoodDads[WinnerDad[i][j]].Allele2[l];
							} // end of l
							Solutions[TotalNumberSolutions].FaEmb[j] = GoodDads[WinnerDad[i][j]].NumberOffspring;
						} // end of j

						TotalNumberSolutions++;
						if (TotalNumberSolutions > TotalNumberSolutionsMax)
						{
							TotalNumberSolutions = TotalNumberSolutionsMax;
							TotalNumberSolutionsError = true;
						}

					} // end of i loop

					if (DadsSolved) BigDadsSolved = true;

					// Exhaustive Done ----------------------------------------------------------------------
					// --------------------------------------------------------------------------------------

					// MinimumAlleles -----------------------------------------------------------------------
					//---------------------------------------------------------------------------------------

					// This routine calculates the minimum number of fathers based on the single locus with
					// the largest number of alleles segregating in the progeny array.

					int tempListA[100];
					int tempNumA;
					bool onlist1, allonlist;
					bool MPdone;
					int mostalleles;
					int tempMinDadsOneLocus;

					for (j = 0; j < NumberLoci; j++)
					{
						tempNumA = 0;
						MPdone = false;

						for (i = 0; i < PatAllNum[j]; i++)
						{
							if (PatAllList[j][i] != MomAllele1[j] && PatAllList[j][i] != MomAllele2[j])
							{
								tempListA[tempNumA] = PatAllList[j][i];
								tempNumA++;
							} // end of if
						} // end of i loop

						n = 0;
						while (!MPdone)
						{
							n++;
							if (n == 2)
							{
								tempListA[tempNumA] = MomAllele1[j];
								tempNumA++;
							}

							if (n == 3)
							{
								tempListA[tempNumA-1] = MomAllele2[j];
							}

							if (n == 4)
							{
								tempListA[tempNumA] = MomAllele1[j];
								tempNumA++;
							}

							if (n > 4)
							{
								MPdone = true;
							}

							allonlist = true;
							for (k = 0; k < NumberEmbryos; k++)
							{
								onlist1 = false;
								for (m = 0; m < tempNumA; m++)
								{
									if (EmbryoPaternal1[j][k] == tempListA[m] || EmbryoPaternal2[j][k] == tempListA[m])
										onlist1 = true;
								} // end of m loop
								if (onlist1 == false)
									allonlist = false;
							} // end of k loop

							if (allonlist)
							{
								MinPatAllNum[j] = tempNumA;
								MPdone = true;
							}
						} // end of while loop
					}   // end of j loop

					mostalleles = 0;

					for (j = 0; j < NumberLoci; j++)
					{
						if (MinPatAllNum[j] > mostalleles)
							mostalleles = MinPatAllNum[j];
					}

					tempMinDadsOneLocus = mostalleles/2;
					if (mostalleles % 2 > 0)
						tempMinDadsOneLocus++;

					if (tempMinDadsOneLocus < MinDadsOneLocus)
						MinDadsOneLocus = tempMinDadsOneLocus;

					// End of MinimumAlleles ----------------------------------------------------------------
					//---------------------------------------------------------------------------------------

					} //end of BSi loop

					for (BSi = 0; BSi < TotalNumberSolutions; BSi++)
					{
						if (Solutions[BSi].NumberFathers <= MaxNumDadsEx)
						{
							Solutions[tempTNS] = Solutions[BSi];
							tempTNS++;
						}
					}
					TotalNumberSolutions = tempTNS;

				// End of Big Search --------------------------------------------------------

				// Next step: count up the results

				if (BigDadsSolved)
				{
					if (TotalNumberSolutions > 1)
					{
						// ProbRanking Routine --------------------------------------------------------------
						// ----------------------------------------------------------------------------------



						int NumDadsLeft[1000];
						bool tempflag;
						int NumDadsComp[1000];
						int numAllele1[50][10];
						int numAllele2[50][10];
						bool CompDad[10][1000];
						bool reduceNum;
						bool done;
						bool foundallele;
						bool freqError;
						bool sameloci;
						bool foundallele2;
						int ambigMomAllele[50][1000];
						int numMatAll1[50];
						int numMatAll2[50];
						long double MomChiSqr;

						long double fa, fb, expect, chival, altchi;

						int i, j, k, l, m ,n, o, p, q, r, s, t, u, v, w, x;
						int rowpos, rowpos2;

						for (i = 0; i < TotalNumberSolutions; i++)
						{
							for (q = 0; q < 1000; q++)
							{
								NumDadsComp[q] = 0;
								for (r = 0; r < MaxNumDadsEx; r++)
									CompDad[r][q] = false;
								for (l = 0; l < NumberLoci; l++)
									ambigMomAllele[l][q] = 0;
							}

							for (j = 0; j < NumberLoci; j++)
							{
								MomAllele1[j] = Solutions[i].MoAll1[j];
								MomAllele2[j] = Solutions[i].MoAll2[j];
							}

							// PaternalAlleles->Execute() --------------------------------------------------------------------
							// -----------------------------------------------------------------------------------------------

								int ma1, ma2, ea1, ea2, PAi, PAj;
								for (PAi = 0; PAi < NumberEmbryos; PAi++)
								{
									for (PAj = 0; PAj < NumberLoci; PAj++)
									{
										ma1 = MomAllele1[PAj];
										ma2 = MomAllele2[PAj];
										ea1 = EmbryoAllele1[PAj][PAi];
										ea2 = EmbryoAllele2[PAj][PAi];

										if (compatible(ma1, ma2, ea1, ea2))
										{
											if (ea1 == ma1 || ea1 == ma2) EmbryoPaternal1[PAj][PAi] = ea2;
											if (ea2 == ma1 || ea2 == ma2) EmbryoPaternal1[PAj][PAi] = ea1;
											EmbryoPaternal2[PAj][PAi] = 0;
										}
										else
										{
											EmbryoPaternal1[PAj][PAi] = 9999;
											EmbryoPaternal2[PAj][PAi] = 9999;
										}

										if (identical(ma1, ma2, ea1, ea2) && ma1 != ma2)
										{
											EmbryoPaternal1[PAj][PAi] = ea1;
											EmbryoPaternal2[PAj][PAi] = ea2;
										}
									} // end of j loop
								} // end of i loop

								// This next section will tally up all of the paternal alleles
								int ep1, ep2, PAflag, PAk;

								for (PAj = 0; PAj < NumberLoci; PAj++)
								{
									PatAllList[PAj][0] = 9999;
									PAflag = 0;

									while (PatAllList[PAj][0] == 9999)
									{
										PatAllList[PAj][0] = EmbryoPaternal1[PAj][PAflag];
										PAflag++;
									}

									PatAllNum[PAj] = 1;

									for (PAi = 0; PAi < NumberEmbryos; PAi++)
									{
										ep1 = EmbryoPaternal1[PAj][PAi];
										ep2 = EmbryoPaternal2[PAj][PAi];

										if (ep1 != 0 && ep1 != 9999)
										{
											PAflag = 0;
											for (PAk = 0; PAk < PatAllNum[PAj]; PAk++)
											{
												if (ep1 == PatAllList[PAj][PAk])
													PAflag = 1;
											} // end of k loop
											if (PAflag == 0)
											{
												PatAllList[PAj][PatAllNum[PAj]] = ep1;
												PatAllNum[PAj]++;
											}
										}// end of if

										if (ep2 != 0 && ep2 != 9999)
										{
											PAflag = 0;
											for (PAk = 0; PAk < PatAllNum[PAj]; PAk++)
											{
												if (ep2 == PatAllList[PAj][PAk])
													PAflag = 1;
											} // end of k loop
											if (PAflag == 0)
											{
												PatAllList[PAj][PatAllNum[PAj]] = ep2;
												PatAllNum[PAj]++;
											}
										}// end of if
									} // end of i loop
								} // end of j loop

  
								for (PAi = 0; PAi < 1000; PAi++)
								{
									for (PAj = 0; PAj < 50; PAj++)
									{
										EmbryoPat1[PAj][PAi] = EmbryoPaternal1[PAj][PAi];
										EmbryoPat2[PAj][PAi] = EmbryoPaternal2[PAj][PAi];
									}
								}

							// End of PaternalAlleles ------------------------------------------------------------------------
							// -----------------------------------------------------------------------------------------------

							// This loop asks how many parents the offspring is compatible with.

							for (j = 0; j < NumberEmbryos; j++)
							{
								NumDadsComp[j] = 0;
								for (k = 0; k < MaxNumDadsEx; k++)
								{
									tempflag = true;
									for (l = 0; l < NumberLoci; l++)
									{
										if (!compatible(Solutions[i].FaAll1[l][k],Solutions[i].FaAll2[l][k],
											EmbryoPaternal1[l][j],EmbryoPaternal2[l][j]))
												tempflag = false;
									}
									if (tempflag)
									{
										NumDadsComp[j]++;
										CompDad[k][j] = true;
									}
								}   // end of k
							} // end of j

							for (j = 0; j < NumberEmbryos; j++)
								NumDadsLeft[j] = NumDadsComp[j];

							for (j = 0; j < MaxNumDadsEx; j++)
							{
								for (k = 0; k < NumberLoci; k++)
								{
									numAllele1[k][j] = 0;
									numAllele2[k][j] = 0;
								}

								for (k = 0; k < NumberEmbryos; k++)
								{
									reduceNum = false;
									if (NumDadsLeft[k] > 1 && CompDad[j][k] == true)
									{
										x = randnum(10000);
										v = 10000/NumDadsComp[k];

										if (x < v)
											NumDadsLeft[k] = 1;
										else
											reduceNum = true;
									} // end of if NumDadsLeft[k]>1

									if (NumDadsLeft[k] == 1 && CompDad[j][k] == true)
									{
										for (l = 0; l < NumberLoci; l++)
										{
											done = false;

											if (EmbryoPaternal2[l][k] == 0)
											{
												if (EmbryoPaternal1[l][k] == Solutions[i].FaAll1[l][j])
												{
													numAllele1[l][j]++;
													done = true;
												}

												if (!done && EmbryoPaternal1[l][k] == Solutions[i].FaAll2[l][j])
												{
													numAllele2[l][j]++;
													done = true;
												}
											}

											if (!done && EmbryoPaternal2[l][k] > 0)
											{
												if (identical(Solutions[i].FaAll1[l][j],Solutions[i].FaAll2[l][j],
													EmbryoPaternal1[l][k],EmbryoPaternal2[l][k]))
												{
													x = randnum(10000);
													if (x < 5000)
													{
														numAllele1[l][j]++;
														ambigMomAllele[l][k] = Solutions[i].FaAll2[l][j];
													}
													else
													{
														numAllele2[l][j]++;
														ambigMomAllele[l][k] = Solutions[i].FaAll1[l][j];
													}
													done = true;
												}

												if (!done && (EmbryoPaternal1[l][k] == Solutions[i].FaAll1[l][j]
													|| EmbryoPaternal2[l][k] == Solutions[i].FaAll1[l][j]))
												{
													numAllele1[l][j]++;
													done = true;
													if (EmbryoPaternal1[l][k] == Solutions[i].FaAll1[l][j])
														ambigMomAllele[l][k] = EmbryoPaternal2[l][k];
													else
														ambigMomAllele[l][k] = EmbryoPaternal1[l][k];
												}
												if (!done && (EmbryoPaternal1[l][k] == Solutions[i].FaAll2[l][j]
													|| EmbryoPaternal2[l][k] == Solutions[i].FaAll2[l][j]))
												{
													numAllele2[l][j]++;
													done = true;
													if (EmbryoPaternal1[l][k] == Solutions[i].FaAll2[l][j])
														ambigMomAllele[l][k] = EmbryoPaternal2[l][k];
													else
														ambigMomAllele[l][k] = EmbryoPaternal1[l][k];
												}
											} // end of if != 0
										} // end of l loop

										NumDadsLeft[k] = 0;

									} // end of if NumDadsLeft[k] == 1

									if (reduceNum  == true)
										NumDadsLeft[k] = NumDadsLeft[k] - 1;

								}// end of k loop
							}// end of j loop

							// The above loops counted up the numbers of allele 1 versus allele 2 per male
							// Next, it's necessary to count up the number of allele 1 versus allele 2
							// of the maternal genotype in the progeny array.

							for (l = 0; l < NumberLoci; l++)
							{
								numMatAll1[l] = 0;
								numMatAll2[l] = 0;
							}

							for (k = 0; k < NumberEmbryos; k++)
							{
								for (l = 0; l < NumberLoci; l++)
								{
									done = false;
									if (identical(MomAllele1[l],MomAllele2[l],EmbryoAllele1[l][k],EmbryoAllele2[l][k]))
									{
										if (MomAllele1[l] == ambigMomAllele[l][k])
											numMatAll1[l]++;
										else
											numMatAll2[l]++;
										done = true;
									}

									if (!done && (MomAllele1[l] == EmbryoAllele1[l][k] || MomAllele1[l] == EmbryoAllele2[l][k]))
									{
										numMatAll1[l]++;
										done = true;
									}

									if (!done && (MomAllele2[l] == EmbryoAllele1[l][k] || MomAllele2[l] == EmbryoAllele2[l][k]))
									{
										numMatAll2[l]++;
										done = true;
									}
								}// end of l loop
							} // end of k loop

							if (!AllelesLoaded)   // Calculate probabilities on the basis of Mendelian segregation alone.
							{
								WDChiSqr[i] = 1;
								// Calculate Mendelian probability for mother

								for (m = 0; m < NumberLoci; m++)
								{
									if (MomAllele1[m] != MomAllele2[m])
									{
										//fa = factorial(numMatAll1[m]+numMatAll2[m])/(factorial(numMatAll1[m])*factorial(numMatAll2[m]));

										fa = NchooseK(numMatAll1[m]+numMatAll2[m],numMatAll1[m]);
										fb = 1;
										for (j = 0; j < numMatAll1[m] + numMatAll2[m]; j++)
										{
											fb = fb * 0.5;
										}
										chival = fa*fb;
										WDChiSqr[i] = WDChiSqr[i] * chival;
									}
									else
									{
										chival = 1;
									}
								} // end of m

								// Calculate Mendelian prob for fathers
								for (n = 0; n < MaxNumDadsEx; n++)
								{
									for (m = 0; m < NumberLoci; m++)
									{
										if (Solutions[i].FaAll1[m][n] != Solutions[i].FaAll2[m][n])
										{
											//fa = factorial(numAllele1[m][n]+numAllele2[m][n])/(factorial(numAllele1[m][n])*factorial(numAllele2[m][n]));

											fa = NchooseK(numAllele1[m][n]+numAllele2[m][n],numAllele1[m][n]);
											fb = 1;
											for (j = 0; j < numAllele1[m][n] + numAllele2[m][n]; j++)
											{
												fb = fb * 0.5;
											}
											chival = fa*fb;
											WDChiSqr[i] = WDChiSqr[i] * chival;
										}
										else
										{
											chival = 1;
										}
									} // end of m
								} // end of n
							}  // End of if !AllelesLoaded or freqError.
							w = 0;
							if (AllelesLoaded)
							{
								WDChiSqr[i] = 1;
								// Calculate joint Mendelian and allele frequency probability for mother

								for (m = 0; m < NumberLoci; m++)
								{
									w = m;

									if (MomAllele1[m] != MomAllele2[m])
									{
										//fa = factorial(numMatAll1[m]+numMatAll2[m])/(factorial(numMatAll1[m])*factorial(numMatAll2[m]));

										fa = NchooseK(numMatAll1[m]+numMatAll2[m],numMatAll1[m]);
										fb = 1;
										for (j = 0; j < numMatAll1[m] + numMatAll2[m]; j++)
										{
											fb = fb * 0.5;
										}
										chival = fa*fb;

										for (k = 0; k < Loci[w].NumberofAlleles; k++)
										{
											if (MomAllele1[m] == Loci[w].AlleleName[k])
												chival = chival * Loci[w].AlleleFreq[k];
											if (MomAllele2[m] == Loci[w].AlleleName[k])
												chival = chival * Loci[w].AlleleFreq[k];
										} // end of k
										chival = chival*2;
									} // end of if (MomAllele1[m] != MomAllele2[m])

									if (MomAllele1[m] == MomAllele2[m])
									{
										chival = 1;
										altchi = 1;  //altchi is the probability for an A/? heterozygote.
										for (j = 0; j < numMatAll1[m] + numMatAll2[m]; j++)
										{
											altchi = altchi * 0.5;
										}

										for (k = 0; k < Loci[w].NumberofAlleles; k++)
										{
											if (MomAllele1[m] == Loci[w].AlleleName[k])
											{
												chival = chival * Loci[w].AlleleFreq[k] * Loci[w].AlleleFreq[k];
												altchi = altchi * Loci[w].AlleleFreq[k] * (1 - Loci[w].AlleleFreq[k]) * 2;
											}
										} // end of k

										if (altchi > chival && !MomKnown)
										{
											chival = altchi;
											Solutions[i].MoAll2[m] = 0;
										}
									} // end of if MomAllele1[m] == MomAllele2[m]
									WDChiSqr[i] = WDChiSqr[i] * chival;
								} // end of m

								// Now calculate the probabilities for the males.

								for (m = 0; m < NumberLoci; m++)
								{
									w = m;

									for (n = 0; n < MaxNumDadsEx; n++)
									{
										if (Solutions[i].FaAll1[m][n] != Solutions[i].FaAll2[m][n])
										{
											//fa = factorial(numAllele1[m][n] + numAllele2[m][n])/(factorial(numAllele1[m][n])*factorial(numAllele2[m][n]));

											fa = NchooseK(numAllele1[m][n] + numAllele2[m][n],numAllele1[m][n]);
											fb = 1;
											for (j = 0; j < numAllele1[m][n] + numAllele2[m][n]; j++)
											{
												fb = fb*0.5;
											}
											chival = fa*fb;

											for (k = 0; k < Loci[w].NumberofAlleles; k++)
											{
												if (Solutions[i].FaAll1[m][n] == Loci[w].AlleleName[k])
													chival = chival * Loci[w].AlleleFreq[k];
												if (Solutions[i].FaAll2[m][n] == Loci[w].AlleleName[k])
													chival = chival * Loci[w].AlleleFreq[k];
											}
											chival = chival * 2;
										}  // end of if FaAll1 = FaAll2

										if (Solutions[i].FaAll1[m][n] == Solutions[i].FaAll2[m][n])
										{
											chival = 1;
											altchi = 1;
											for (j = 0; j < numAllele1[m][n] + numAllele2[m][n]; j++)
											{
												altchi = altchi*0.5;
											} // end of j
											for (k = 0; k < Loci[w].NumberofAlleles; k++)
											{
												if (Solutions[i].FaAll1[m][n] == Loci[w].AlleleName[k])
												{
													chival = chival * Loci[w].AlleleFreq[k] * Loci[w].AlleleFreq[k];
													altchi = chival * Loci[w].AlleleFreq[k] * (1 - Loci[w].AlleleFreq[k]) * 2;
												}
											} // end of k

											if (altchi > chival)
											{
												chival = altchi;
												Solutions[i].FaAll2[m][n] = 0;
											}
										} // end of if FaAll1 = FaAll2
										WDChiSqr[i] = WDChiSqr[i] * chival;
									} // end of n
								} // end of m
							} // end of if AllelesLoaded

						} // end of i loop

						for (i = 0; i < TotalNumberSolutions; i++)
							WDPriority[i] = WDChiSqr[i];

						//This next loop will keep the 50 best combinations of fathers (in order).

						NumBestDads = 50;

						if (NumBestDads > TotalNumberSolutions)
							NumBestDads = TotalNumberSolutions;

						long double lowestpriority = 1.1;
						int LPgroup;
						LPgroup = 0;

						for (i = 0; i < NumBestDads; i++)
						{
							BestDadsPriority[i] = 0;
						}

						for (i = 0; i < TotalNumberSolutions; i++)
						{
							if (WDPriority[i] < lowestpriority)
							{
								lowestpriority = WDPriority[i];
								LPgroup = i;
							}
						}

						for (k = 0; k < NumBestDads; k++)
						{
							BestDadsPriority[k] = lowestpriority;
							BestSolutions[k] = Solutions[LPgroup];
						}

						for (i = 0; i < TotalNumberSolutions; i++)
						{
							j = 0;
							while (WDPriority[i] < BestDadsPriority[j] && j < NumBestDads)
								j++;

							for (k = NumBestDads; k > j; k--)
							{
								BestDadsPriority[k] = BestDadsPriority[k-1];
								BestSolutions[k] = BestSolutions[k-1];
							}

							BestDadsPriority[j] = WDPriority[i];
							BestSolutions[j] = Solutions[i];

						}// end of i

						// This routine, not found in the actual analysis program, randomizes the order
						// of those solutions with the same priority score.  It's not necessary in the
						// Analysis program, because the user will see whether or not there are several
						// solutions with the same priority score.

						j = 1;
						while (BestDadsPriority[0] == BestDadsPriority[j])
							j++;

						x = randnum(j);
						BestSolutions[0] = BestSolutions[x];




						// End of ProbRanking ---------------------------------------------------------------
						// ----------------------------------------------------------------------------------
					}
					else
						BestSolutions[0] = Solutions[0];

					// Check Results Routine ----------------------------------------------------------------
					// --------------------------------------------------------------------------------------
					
					bool UsedRDad[10];
					bool UsedSDad[10];
					bool DadGenoMatch[10];
					bool DadExactMatch[10];
					int NumberCorrect;
					bool tFlag;
					bool tFlagSlvd;

					for (i = 0; i < 10; i++)
					{
						DadExactMatch[i] = false;
						DadGenoMatch[i] = false;
					}

					ExactlyRight = false;
					GenosRight = false;

					if (RealMaleNumber[run] == MaxNumDadsEx)
					{
						CorrectDadNumber = true;
						OverDadNumber = false;
						UnderDadNumber = false;
					}

					if (RealMaleNumber[run] < MaxNumDadsEx)
					{
						CorrectDadNumber = false;
						OverDadNumber = true;
						UnderDadNumber = false;
					}

					if (RealMaleNumber[run] > MaxNumDadsEx)
					{
						CorrectDadNumber = false;
						OverDadNumber = false;
						UnderDadNumber = true;
					}

					if (CorrectDadNumber)
					{
						for (i = 0; i < 6; i++)
							UsedRDad[i] = false;
						for (i = 0; i < 6; i++)
							UsedSDad[i] = false;

						for (i = 0; i < RealMaleNumber[run]; i++)
						{
							for (j = 0; j < MaxNumDadsEx; j++)
							{
								if (!UsedRDad[i] && !UsedSDad[j])
								{
									tFlag = true;
									for (k = 0; k < NumberLoci; k++)
									{
										if (!identical(RealDad[i].Allele1[k], RealDad[i].Allele2[k], BestSolutions[0].FaAll1[k][j], BestSolutions[0].FaAll2[k][j]))
											tFlag = false;
									}  // end of k loop

									if (tFlag)
									{
										UsedRDad[i] = true;
										UsedSDad[j] = true;
										DadGenoMatch[i] = true;

										if (EmsSampFromMale[i] == BestSolutions[0].FaEmb[j])
											DadExactMatch[i] = true;

									}
								} // end of if
							} // end of j loop
						} // end of i loop

						tFlag = true;
						for (i = 0; i < RealMaleNumber[run]; i++)
						{
							if (!DadExactMatch[i])
								tFlag = false;
						} // end of i loop
						if (tFlag)
							ExactlyRight = true;
						else
							ExactlyRight = false;


						tFlag = true;
						for (i = 0; i < RealMaleNumber[run]; i++)
						{
							if (!DadGenoMatch[i])
								tFlag = false;
						}
						if (tFlag)
							GenosRight = true;

					} // end of if correctdadnumber

					// Now figure out how many of the real dad genotypes are in the solved
					// dad pool.  And for those that aren't, try to figure out why not.

					for (i = 0; i < 10; i++)
						UsedSDad[i] = false;

					NumberDadsCorrect = 0;
					NumberDadsWrong = 0;
					NumberFiveCorrect = 0;
					NumberFiveWrong = 0;
					NumberTenCorrect = 0;
					NumberTenWrong = 0;
					NumberTwentCorrect = 0;
					NumberTwentWrong = 0;
					NumberMoreCorrect = 0;
					NumberMoreWrong = 0;
					NumberNotSampledWrong = 0;
					NumberNotSampledCorrect = 0;

					for (i = 0; i < RealMaleNumber[run]; i++)
					{
						tFlagSlvd = false;
						for (j = 0; j < MaxNumDadsEx; j++)
						{
							if (!UsedSDad[j])
							{
								tFlag = true;
								for (k = 0; k < NumberLoci; k++)
								{
									if (!identical(RealDad[i].Allele1[k], RealDad[i].Allele2[k], BestSolutions[0].FaAll1[k][j], BestSolutions[0].FaAll2[k][j]))
										tFlag = false;
								} // end of k loop

								if (tFlag)
								{
									UsedSDad[j] = true;
									tFlagSlvd = true;
								}
							}// end of if
						} // end of j loop

						if (tFlagSlvd)
						{
							NumberDadsCorrect++;
							if (EmsSampFromMale[i] == 0)
								NumberNotSampledCorrect++;
							if (EmsSampFromMale[i] > 0 && EmsSampFromMale[i] < 6)
								NumberFiveCorrect++;
							if (EmsSampFromMale[i] > 5 && EmsSampFromMale[i] < 11)
								NumberTenCorrect++;
							if (EmsSampFromMale[i] > 10 && EmsSampFromMale[i] < 21)
								NumberTwentCorrect++;
							if (EmsSampFromMale[i] > 20)
								NumberMoreCorrect++;
						}
						if (!tFlagSlvd)
						{
							NumberDadsWrong++;
							if (EmsSampFromMale[i] == 0)
								NumberNotSampledWrong++;
							if (EmsSampFromMale[i] > 0 && EmsSampFromMale[i] < 6)
								NumberFiveWrong++;
							if (EmsSampFromMale[i] > 5 && EmsSampFromMale[i] < 11)
								NumberTenWrong++;
							if (EmsSampFromMale[i] > 10 && EmsSampFromMale[i] < 21)
								NumberTwentWrong++;
							if (EmsSampFromMale[i] > 20)
								NumberMoreWrong++;

						}
					}// end of i loop

					tFlag = true;
					for (i = 0; i < NumLoci; i++)
					{
						if (!identical(RealMomAllele1[i],RealMomAllele2[i],BestSolutions[0].MoAll1[i],BestSolutions[0].MoAll2[i]))
							tFlag = false;
					}

					if (tFlag)
						MotherCorrect = true;
					else
						MotherCorrect = false;

					// End of Check Results -----------------------------------------------------------------
					// --------------------------------------------------------------------------------------

					if (CorrectDadNumber)
						CumCorrDadNum++;
					if (OverDadNumber)
					{
						CumOverDadNum++;
					}
					if (UnderDadNumber)
						CumUnderDadNum++;
					CumNumDadsCorr = NumberDadsCorrect + CumNumDadsCorr;
					CumNumDadsWrong = NumberDadsWrong + CumNumDadsWrong;
					CumNumFiveCorr = NumberFiveCorrect + CumNumFiveCorr;
					CumNumFiveWrong = NumberFiveWrong + CumNumFiveWrong;
					CumNumTenCorr = NumberTenCorrect + CumNumTenCorr;
					CumNumTenWrong = NumberTenWrong + CumNumTenWrong;
					CumNumTwentCorr = NumberTwentCorrect + CumNumTwentCorr;
					CumNumTwentWrong = NumberTwentWrong + CumNumTwentWrong;
					CumNumMoreCorr = NumberMoreCorrect + CumNumMoreCorr;
					CumNumMoreWrong = NumberMoreWrong + CumNumMoreWrong;
					CumNumNSampWrong = NumberNotSampledWrong + CumNumNSampWrong;
					CumNumNSampCorr = NumberNotSampledCorrect + CumNumNSampCorr;

					if (ExactlyRight)
					{
						CumExactlyRight++;
					}
					else
					{
						if(GenosRight)
							CumGenosRight++;
						else
							CumGenosWrong++;
					}

					if (MinDadsOneLocus == RealMaleNumber[run])
						MinDadsOL++;
					if (MinDadsOneLocus < RealMaleNumber[run])
						MinDadsOLLess++;
					if (MinDadsOneLocus > RealMaleNumber[run])
						MinDadsOLMore++;

					if (ExactlyRight || GenosRight)
					{
						if (TotalNumberSolutions == 1)
							NumOneSolnCorr++;
						if (TotalNumberSolutions > 1 && TotalNumberSolutions < 11)
							NumTenSolnsCorr++;
						if (TotalNumberSolutions > 10)
							NumMoreSolnsCorr++;
					}
					else
					{
						if (TotalNumberSolutions == 1)
							NumOneSolnWrong++;
						if (TotalNumberSolutions > 1 && TotalNumberSolutions < 11)
							NumTenSolnsWrong++;
						if (TotalNumberSolutions > 10)
							NumMoreSolnsWrong++;
					}

					if (MotherCorrect)
						CumMotherCorrect++;

				} // end of if dadssolved
				else
				{
					NumNotSolved++;
				}

			} // end of iteration loop

			// Write the results from this run to the output file

			o_file << "Results from Run " << run + 1 << ":\n\n";
			o_file << "SEARCH PARAMETERS:\n";
			if (MomKnown)
				o_file << "Mother Known -- Exhaustive Search\n";
			else
				o_file << "Unknown Mother -- Exhaustive Search\n";

			o_file << "Loci:\t";
			for (i = 0; i < NumberLoci; i++)
				o_file << Loci[i].Name << "\t";
			o_file << "\n";

			o_file << "Number of iterations = " << NumIterations << "\n"; 
			o_file << "Number of offspring assayed = " << NumberEmbryos << "\n";
			o_file << "Actual number of fathers = " << RealMaleNumber[run] << "\n";

			for (i = 0; i < RealMaleNumber[run]; i++)
				o_file << "Number sired in brood by male " << i+1 << ": " << NumFromMale[i][run] << "\n";

			o_file << "\n";
			o_file << "Number of FAMILIES with:" << "\n";
			o_file << "All fathers exactly right:    \t" << CumExactlyRight << "\n";
			o_file << "All paternal genotypes right: \t" << CumGenosRight << "\n";
			o_file << "One or more genotype wrong:   \t" << CumGenosWrong << "\n";
			o_file << "One unique solution (correct):\t" << NumOneSolnCorr << "\n";
			o_file << "One unique solution (wrong):  \t" << NumOneSolnWrong << "\n";
			o_file << "Two to ten sol'ns (correct):  \t" << NumTenSolnsCorr << "\n";
			o_file << "Two to ten sol'ns (wrong):    \t" << NumTenSolnsWrong << "\n";
			o_file << "More than ten sol'ns (right): \t" << NumMoreSolnsCorr << "\n";
			o_file << "More than ten sol'ns (wrong): \t" << NumMoreSolnsWrong << "\n";
			o_file << "\n";
			o_file << "Number of PATERNAL GENOTYPES:" << "\n";
			o_file << "Total number of paternal genotypes right:\t" << CumNumDadsCorr << "\n";
			o_file << "Total number of paternal genotypes wrong:\t" << CumNumDadsWrong << "\n";
			o_file << "Males with 0-5 offspring right:       \t" << CumNumFiveCorr << "\n";
			o_file << "Males with 0-5 offspring wrong:       \t" << CumNumFiveWrong << "\n";
			o_file << "Males with 6-10 offspring right:      \t" << CumNumTenCorr << "\n";
			o_file << "Males with 6-10 offspring wrong:      \t" << CumNumTenWrong << "\n";
			o_file << "Males with 11-20 offspring right:     \t" << CumNumTwentCorr << "\n";
			o_file << "Males with 11-20 offspring wrong:     \t" << CumNumTwentWrong << "\n";
			o_file << "Males with 21+ offspring right:       \t" << CumNumMoreCorr << "\n";
			o_file << "Males with 21+ offspring wrong:       \t" << CumNumMoreWrong << "\n";
			o_file << "Unsampled males, right by chance:     \t" << CumNumNSampCorr << "\n";
			o_file << "Unsampled males, wrong (no offspring):\t" << CumNumNSampWrong << "\n\n";

			o_file <<"Number of MATERNAL GENOTYPES:" << "\n";
			o_file << "Maternal Genotype correct:     \t" << CumMotherCorrect << "\n\n";

			o_file <<"Number of ITERATIONS in which: \n";
			o_file << "Reconstructed dad number equals real dad number:      \t" << CumCorrDadNum << "\n";
			o_file << "Reconstructed dad number greater than real dad number:\t" << CumOverDadNum << "\n";
			o_file << "Reconstructed dad number less than real dad number:   \t" << CumUnderDadNum << "\n\n";

			o_file << "One-locus minimum equals real dad number:         \t" << MinDadsOL << "\n";
			o_file << "One-locus minimum is greater than real dad number:\t" << MinDadsOLMore << "\n";
			o_file << "One-locus minimum is less than real dad number:   \t" << MinDadsOLLess << "\n\n";

	 }  // end of run loop

	 if (ErrorPP)
		o_file << "There were PP errors.\n";

	if (ErrorGP)
		o_file << "There were GP Errors.\n";

	if (TotalNumberSolutionsError)
		o_file << "There were TNS Errors.\n";
   
	if (MAnumberError)
		o_file << "There were MAnumber Errors.\n";

	if (NumWinDadError)
		o_file << "There were NumWinDad Errors.\n";

	cout << "Done!\n";
	

	o_file.close();

	delete[] Loci;
	delete[] GoodDads;
	delete[] Mothers;
	delete[] Solutions;
	delete[] BestSolutions;
	delete[] BestSolutionsOld;
	delete[] PossDads;
	delete[] PossMoms;
	delete[] WDChiSqr;
	delete[] WDFreq;
	delete[] WDPriority;

	delete[] AlleleFrequencyFilename;
	delete[] ParameterFilename;
	delete[] OutputFilename;
	delete[] OneArrayFilename;

	return 0;
}