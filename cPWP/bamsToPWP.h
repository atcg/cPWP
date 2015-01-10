//
//  readsToPWP.h
//  cPWP
//
//  Created by Evan McCartney-Melstad on 12/31/14.
//  Copyright (c) 2014 Evan McCartney-Melstad. All rights reserved.
//

#ifndef __cPWP__readsToPWP__
#define __cPWP__readsToPWP__

#include <vector>
#include <string>

int runANGSDforReadCounts (std::string bamlist, std::string angsdPrefix, std::string nThreads="25", std::string angsdOutputLog="angsdOutPut.log");

int convertANGSDcountsToBinary(std::string angsdPrefix, std::string binaryOutputFileName, int numIndividuals, int readDepthMax);

unsigned char readCountToUnsignedChar(int individual, std::string allele, std::vector<std::string> countsVector);

int calcPWPfromBinaryFile (std::string binaryFile, unsigned long long int numLoci, const int numIndividuals, std::string outFile, int numThreads = 1);

//int calcPWPforRange (unsigned long long startingLocus, unsigned long long endingLocus, int numIndividuals, std::vector<unsigned char> & mainReadCountVector, std::vector<std::vector<long double>> & threadPWP, std::vector<std::vector<unsigned long long int>> & threadWeightings);

int calcPWPforRange (unsigned long long startingLocus, unsigned long long endingLocus, int numIndividuals, std::vector<unsigned char>& mainReadCountVector, std::vector<std::vector<long double>>& threadPWP, std::vector<std::vector<unsigned long long int>>& threadWeightings);


#endif /* defined(__cPWP__readsToPWP__) */
