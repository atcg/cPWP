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


#endif /* defined(__cPWP__readsToPWP__) */