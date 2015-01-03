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
int generateReadsAndMap (int numIndividuals, double mutationRateStepSize, std::string libFragmentSize, std::string stdevLibFragmentSize, std::string depth, std::string readLengths, std::string randomSeed, std::string reference, std::string threads = "10");



#endif /* defined(__cPWP__generateSimulatedData__) */
