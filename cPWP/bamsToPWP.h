//
//  readsToPWP.h
//  cPWP
//
//  Created by Evan McCartney-Melstad on 12/31/14.
//  Copyright (c) 2014 Evan McCartney-Melstad. All rights reserved.
//

#ifndef __cPWP__readsToPWP__
#define __cPWP__readsToPWP__

#include <string>
#include <vector>

int runANGSDforReadCounts (std::string bamlist, std::string angsdPrefix, std::string nThreads="25", std::string angsdOutputLog="angsdOutPut.log");
int convertANGSDcountsToBinary(std::string angsdPrefix, std::string binaryOutputFileName, int numIndividuals, int readDepthMax);
unsigned char readCountToUnsignedChar(int individual, std::string allele, std::vector<std::string> countsVector);
int calcPWPfromBinaryFile (std::string binaryFile, int numLoci, const int numIndividuals, std::string outFile);

#endif /* defined(__cPWP__readsToPWP__) */
