#pragma once

// This file is part of GerudSim3.  This version is version 3.0, the first working version
// of GerudSim3.  Any changes should be indicated here.  This version is ANSI compliant and
// should compile on any standard C++ compiler.  This file was last altered on May 2, 2013.

   int ConvertChInt (char character1)
   {
      if (character1 == '1') return 1;
      if (character1 == '2') return 2;
      if (character1 == '3') return 3;
      if (character1 == '4') return 4;
      if (character1 == '5') return 5;
      if (character1 == '6') return 6;
      if (character1 == '7') return 7;
      if (character1 == '8') return 8;
      if (character1 == '9') return 9;
      return 0;
   }

   double ConvertChFloat (char character1)
   {
      if (character1 == '1') return 1.0000000000;
      if (character1 == '2') return 2.0000000000;
      if (character1 == '3') return 3.0000000000;
      if (character1 == '4') return 4.0000000000;
      if (character1 == '5') return 5.0000000000;
      if (character1 == '6') return 6.0000000000;
      if (character1 == '7') return 7.0000000000;
      if (character1 == '8') return 8.0000000000;
      if (character1 == '9') return 9.0000000000;
      return 0.0000000000;
   }

   inline bool compatible(int father1, int father2, int embryo1, int embryo2)
                {
                  if (embryo1 == father1) return true;
                  if (embryo1 == father2) return true;
                  if (embryo2 == father1) return true;
                  if (embryo2 == father2) return true;
                  return false;
                }

   inline bool identical(int father1, int father2, int embryo1, int embryo2)
                {
                  if (embryo1 == father1 && embryo2 == father2) return true;
                  if (embryo2 == father1 && embryo1 == father2) return true;
                  return false;
                }

   inline long double factorial(int factnum)
                {
                  int dbi;
                  long double dbda, dbdb;
                  dbda = 1;
                  for (dbi = 2; dbi < factnum+1; dbi++)
                  {
                    dbdb = dbi;
                    dbda = dbda * dbdb;
                  }
                  return dbda;
                }

   	inline long double NchooseK (int NNN, int KKK)
	{
		if (KKK == 0 || KKK == NNN)
			return 1;

		if (KKK < 0 || KKK > NNN)
			return 0;

		long double NKresult;
		long double NKdub1, NKdub2;

		int NKi;
		double dNNN = static_cast<double>(NNN);
		double dKKK = static_cast<double>(KKK);

		NKresult = 1;
		for (NKi = 1; NKi < NNN + 1; NKi++)
		{
			NKdub1 = NKi;
			if (NKi <= dKKK)
				NKdub2 = NKdub1;
			else
				NKdub2 = NKdub1 - dKKK;

			NKresult = NKresult * (NKdub1/NKdub2);
		}
		return NKresult;
	}


    class Locus
        {
          public:
             char Name[50];
             int AlleleName[200];
             double AlleleFreq[200];
             int NumberofAlleles;

          protected:
          private:

        };

     class Dad
        {
           public:
              int Allele1[50];
              int Allele2[50];
              int NumberOffspring;
           protected:
           private:
        };

     class ParentSet
        {
            public:
              int NumberFathers;
              int MoAll1[50];
              int MoAll2[50];
              int FaAll1[50][10];
              int FaAll2[50][10];
              int FaEmb[10];
            protected:
            private:
        };
