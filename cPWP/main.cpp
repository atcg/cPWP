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
    convertANGSDcountsToBinary("272torts_snp1e6_minmapq20minq30", "272torts_snp1e6_minmapq20minq30.binarycounts", 272, 10);
    
    
    calcPWPfromBinaryFile ("272torts_snp1e6_minmapq20minq30.binarycounts", 56575856, 272);
}
