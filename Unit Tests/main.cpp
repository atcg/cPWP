//
//  main.cpp
//  Unit Tests
//
//  Created by Evan McCartney-Melstad on 12/31/14.
//  Copyright (c) 2014 Evan McCartney-Melstad. All rights reserved.
//

#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "../cPWP/generateSimulatedData.h"
#include "../cPWP/bamsToPWP.h"
#include "catch.hpp"
#include <string>




TEST_CASE( "Simulated reads are generated", "[generateReads]" ) {
    // Make sure that the read simulation finishes
    REQUIRE( generateReadsAndMap(2, 0.01, "0.0", "300", "50", "100000", "100", "1234", "scaffold_0.fasta") == 0);

}

TEST_CASE( "Run ANGSD on simulated reads", "[runANGSD]" ) {
    REQUIRE( runANGSDforReadCounts("bamlist.txt", "angsdOut", "10", "angsdOutLog.txt") == 0);
}

TEST_CASE( "Convert ANGSD read counts to unsigned chars for major and minor counts", "[convertCountsToBinary]") {
    REQUIRE( convertANGSDcountsToBinary("angsdOut", "angsdOut.readCounts.binary", 3, 50) == 0); // 3 individuals, not 2, because generateReadsAndMap actually generates n+1 individuals, since one individual is identical to the reference genome
}

TEST_CASE( "Calculate PWP from the binary representations of the ANGSD readcounts", "[calcPWP]") {
    REQUIRE( calcPWPfromBinaryFile ("angsdOut.readCounts.binary", 7000000, 3) == 0);
}
