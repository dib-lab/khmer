#include "seqfuncs.hh"

using namespace std;

char getBase(char base)
{
   switch (base)  {
      case 'A':
         return 0;
      case 'C':
         return 1;
      case 'G':
         return 2;
      case 'T':
         return 3;
   }
}

char getComplement(char base)
{
   char retVal;

   switch (base)  {
      case 'A':
         retVal = 'T'; break;
      case 'T':
         retVal = 'A'; break;
      case 'C':
         retVal = 'G'; break;
      case 'G':
         retVal = 'C'; break;
   }

   return retVal;
}

std::string getReverseComplement(std::string seq)
{
   std::string revComp;

   for (int i = seq.length()-1; i >= 0; i--)
   {
      revComp += getComplement(seq[i]);
   }

   return revComp;
}

