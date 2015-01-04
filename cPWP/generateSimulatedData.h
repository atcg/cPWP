//
//  generateSimulatedData.h
//  cPWP
//
//  Created by Evan McCartney-Melstad on 12/31/14.
//  Copyright (c) 2014 Evan McCartney-Melstad. All rights reserved.
//

#ifndef __cPWP__generateSimulatedData__
#define __cPWP__generateSimulatedData__


#include <string>
/*
int generateReadsAndMap (int numIndividuals, double mutationRateStepSize, std::string baseErrorRate, std::string libFragmentSize, std::string stdevLibFragmentSize, std::string numReadPairs, std::string readLengths, std::string randomSeed, std::string reference, std::string threads = "10");
*/
int createReferenceGenome (int totalBases, double gcContent, std::string genomeOutFile);

int generatePerfectReads (std::string reference, unsigned int stagger, unsigned int readLengths, unsigned int fragmentLengths, std::string readPrefix);

int createMutatedGenome (std::string reference, std::string mutatedReferenceFile, float percDivergent);

int mapReads (std::string reference, std::string R1file, std::string R2file, std::string outBam, std::string threads = "25");

int generateReadsAndMap (int numIndividuals, double mutationRateStepSize, std::string libFragmentSize, std::string stdevLibFragmentSize, std::string numReadPairs, std::string readLengths, std::string randomSeed, std::string reference, std::string threads = "25");



#endif /* defined(__cPWP__generateSimulatedData__) */
