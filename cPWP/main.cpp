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
#include "generateSimulatedData.h"
#include "bamsToBin.h"
//#include "calcPWP.h"
#include "calcPWPchunks.h"


int main (int argc, char *argv[]) {
    
    
    // Throw out any loci that have an individual with at least 10 reads at that locus
    // Actually, the Perl version is way faster--use that
    // convertANGSDcountsToBinary("272torts_snp1e6_minmapq20minq30", "272torts_snp1e6_minmapq20minq30.binarycounts", 272, 10);
    
    /* Full 272 tort SNP list: calcPWPfromBinaryFile ("272torts_snp1e6_minmapq20minq30.binarycounts", 56575856, 272);
     */
     calcPWPfromBinaryFile (argv[1], atoi(argv[2]), 272, argv[3], atoi(argv[4]), atoi(argv[5]));
    //calcPWPfromBinaryFile (std::string binaryFile, unsigned long long int numLoci, const int numIndividuals, std::string outFile, int lociChunkSize, int numThreads=30);

    // First supply the binary readcounts file, then the number of loci to consider, then the number of individuals, then the output file
    // calcPWPfromBinaryFile (std::string binaryFile, unsigned long long int numLoci, const int numIndividuals, std::string outFile, int numThreads = 10);
}
