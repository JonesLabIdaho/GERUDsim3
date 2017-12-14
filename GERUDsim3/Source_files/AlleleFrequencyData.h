#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <sstream>
using namespace std;

// This file is part of Gerud3.  This version is version 3.0, the first working version
// of Gerud3.  Any changes should be indicated here.  This version is ANSI compliant and
// should compile on any standard C++ compiler.  This file was last altered on Aug 31, 2015.


class AlleleFrequencies
{
public:
	int NumberLoci;
	string LocusName[100];
	int AlleleName[100][100];  // first index is allele and second is locus
	double AlleleFreq[100][100];  // first index is allele and second is locus
	int NumberAlleles[100];  // total number of alleles at the locus

	int LoadAlleleFreqs(char* FreqFile)
	{
		int i, j;
		string tstr;
		ifstream f_file;
		f_file.open(FreqFile);
		if (!f_file.good())
		{
			return 0;
		}
		string sstr;
		string nextlocusname;
		bool endoflocus;
		int allelectr;
		char buffer[256];
		stringstream ss;

		f_file.getline(buffer,255);
		ss.str(buffer);
		ss >> NumberLoci;

		getline(f_file,nextlocusname);

		for (i = 0; i < NumberLoci; i++)
		{
			allelectr = 0;
			endoflocus = false;
			LocusName[i] = nextlocusname;
			while(!endoflocus && !f_file.eof())
			{
				f_file.getline(buffer,255);
				if (buffer[0] != '!')
				{
					ss.clear();
					ss.str(buffer);
					ss >> AlleleName[allelectr][i];
					ss >> AlleleFreq[allelectr][i];
					allelectr++;
				}
				else
				{
					endoflocus = true;
					ss.str(buffer);
					nextlocusname = ss.str();
				}
			}
			NumberAlleles[i] = allelectr;
		} // end of i


		f_file.close();

		if (NumberLoci > 0) 
			return 1;
		return 0;
	}
};