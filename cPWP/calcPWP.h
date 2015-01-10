//
//  calcPWP.h
//  
//
//  Created by Evan McCartney-Melstad on 1/10/15.
//
//

#ifndef ____calcPWP__
#define ____calcPWP__

#include <vector>
#include <string>

int calcPWPfromBinaryFile (std::string binaryFile, unsigned long long int numLoci, const int numIndividuals, std::string outFile, std::string numThreads="10");

int calcPWPforRange (unsigned long long startingLocus, unsigned long long endingLocus, int numIndividuals, std::vector<unsigned char>& mainReadCountVector, std::vector<std::vector<long double>>& threadPWP, std::vector<std::vector<unsigned long long int>>& threadWeightings);

#endif /* defined(____calcPWP__) */
