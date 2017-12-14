#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
using namespace std;

// This file is part of GerudSim3.  This version is version 3.0, the first working version
// of GerudSim3.  Any changes should be indicated here.  This version is ANSI compliant and
// should compile on any standard C++ compiler.  This file was last altered on May 2, 2013.

class GerudSimParameters
{
public:
	bool KnownMother;
	int NumberLoci;
	string locusname[100];
	int NumberIterations;
	int NumberRuns;			  // Maximum number of runs is 100
	int EmbryosAssayed[100];  // Number of embryos assayed in the given run
	int EmbryosFromMale[100][10];  // Number of embryos in the progeny array from each of the 10 possible males
	char buffer[256];

	int LoadParameters(char* ParamFile)
	{
		int i, j, k;
		int substrbegin, substrend;
		string tempstring;
		string tempsubstr;
		ifstream p_file;
		p_file.open(ParamFile);
		if (!p_file.good())
		{
			return 0;
		}

		KnownMother = true;
		getline(p_file,tempstring);
		if (tempstring[0] == 'U' || tempstring[0] == 'u')
			KnownMother = false;

		getline(p_file,tempstring);
		substrend = tempstring.find_first_not_of("1234567890");
		tempsubstr = tempstring.substr(0,substrend);
		for (k = 0; k < static_cast<int>(tempsubstr.length()); k++)
			buffer[k] = tempsubstr[k];
		buffer[k] = '\0';
		NumberLoci = atoi(buffer);

		if (NumberLoci > 50)
			NumberLoci = 50;

		

		// This commented out code would retrieve locus names from the parameter file.
		// For now, however, we just want to use the first NumberLoci loci in the allele frequency
		// file, as the original Gerudsim1.0 does.
		//positioninstring = 0;
		//for(i = 0; i < NumberLoci; i++) 
		//{
		//	substrbegin = tempstring.find("!",positioninstring);
		//	positioninstring = substrbegin;
		//	substrend = tempstring.find("\t",positioninstring);
		//	locusname[i] = tempstring.substr(substrbegin,substrend-substrbegin);
		//	positioninstring = substrend;

		//	cout << locusname[i] << "\n";
		//}

		getline(p_file,tempstring);
		substrend = tempstring.find_first_not_of("0123456789");
		tempsubstr = tempstring.substr(0,substrend);
		for (k = 0; k < static_cast<int>(tempsubstr.length()); k++)
			buffer[k] = tempsubstr[k];
		buffer[k] = '\0';
		NumberIterations = atoi(buffer);

		if (NumberIterations > 10000)
			NumberIterations = 10000;

		

		getline(p_file,tempstring);
		substrend = tempstring.find_first_not_of("0123456789");
		tempsubstr = tempstring.substr(0,substrend);
		for (k = 0; k < static_cast<int>(tempsubstr.length()); k++)
			buffer[k] = tempsubstr[k];
		buffer[k] = '\0';
		NumberRuns = atoi(buffer);

		if (NumberRuns > 100)
			NumberRuns = 100;

		

		getline(p_file,tempstring);

		for (i = 0; i < NumberRuns; i++)
		{
			getline(p_file,tempstring);
			substrbegin = tempstring.find_first_not_of("0123456789") + 1;

			substrend = tempstring.find_first_not_of("0123456789",substrbegin);
			tempsubstr = tempstring.substr(substrbegin,substrend);
			for (k = 0; k < static_cast<int>(tempsubstr.length()); k++)
				buffer[k] = tempsubstr[k];
			buffer[k] = '\0';
			EmbryosAssayed[i] = atoi(buffer);
			substrbegin = substrend + 1;

			for (j = 0; j < 10; j++)
			{
				substrend = tempstring.find_first_not_of("0123456789",substrbegin);
				tempsubstr = tempstring.substr(substrbegin,substrend);
				for (k = 0; k < static_cast<int>(tempsubstr.length()); k++)
					buffer[k] = tempsubstr[k];
				buffer[k] = '\0';
				EmbryosFromMale[i][j] = atoi(buffer);
				substrbegin = substrend + 1;
			} // j
		} // i



		if (NumberRuns > 0 && NumberIterations > 0 && NumberLoci > 0)
			return 1;

		return 0;

		p_file.close();
	} // end of LoadParameters

};