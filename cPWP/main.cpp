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
#include "bamsToPWP.h"


int main (int argc, char *argv[]) {
    
    
    // Throw out any loci that have an individual with at least 10 reads at that locus
    // Actually, the Perl version is way faster--use that
    // convertANGSDcountsToBinary("272torts_snp1e6_minmapq20minq30", "272torts_snp1e6_minmapq20minq30.binarycounts", 272, 10);
    
    /* Full 272 tort SNP list: calcPWPfromBinaryFile ("272torts_snp1e6_minmapq20minq30.binarycounts", 56575856, 272);
     */
     calcPWPfromBinaryFile ("allSNPs272torts.binary8bitunsigned", atoi(argv[1]), 272, argv[2]);
}
