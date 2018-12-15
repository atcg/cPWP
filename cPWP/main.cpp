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
    if (std::string(argv[1]) == "pwp") {
        // PWP:
        calcPWPfromBinaryFile (argv[2], atoi(argv[3]), atoi(argv[4]), argv[5], atoi(argv[6]), argv[7], atoi(argv[8]));
    } else if (std::string(argv[1]) == "covar") {
        // Covariance:
        calcCOVARfromBinaryFile (argv[2], atoi(argv[3]), atoi(argv[4]), argv[5], atoi(argv[6]), argv[7], atoi(argv[8]));
    }
}
