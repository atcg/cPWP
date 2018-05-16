//
//  main.cpp
//  cPWP
//
//  Created by Evan McCartney-Melstad on 12/31/14.
//  Copyright (c) 2014 Evan McCartney-Melstad. All rights reserved.
//

/*
 Program for calculating pairwise pi from a file full of unsigned char integers
*/
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include "calcPWPchunks.h"
#include "calcCOVARchunks.h"


int main (int argc, char *argv[]) {
    // PWP:
    calcPWPfromBinaryFile (argv[1], atoi(argv[2]), atoi(argv[3]), argv[4], atoi(argv[5]), argv[6], atoi(argv[7]));

    // Covariance:
    //calcCOVARfromBinaryFile(std::string binaryFile, long long int numLoci, const int numIndividuals, std::string outFile, int lociChunkSize, const int numThreads)
    //calcCOVARfromBinaryFile(argv[1], atoi(argv[2]), 272, argv[3], atoi(argv[4]), atoi(argv[5]));
}
