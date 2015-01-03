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



/*
TEST_CASE( "Simulated reads are generated", "[generateReads]" ) {
    // Make sure that the read simulation finishes
    REQUIRE( generateReadsAndMap(1, 0.01, "0.0", "300", "50", "1000000", "100", "1234", "scaffold_0.fasta") == 0);
}
 */

TEST_CASE( "Generate mutated reference genomes and simulate reads", "[genomeAndReadSim]") {
    REQUIRE( generateReadsAndMap(2, 0.01, "300", "25", "20", "100", "1234", "scaffold_0.fasta", "10") == 0);
}

TEST_CASE( "Run ANGSD on simulated reads", "[runANGSD]" ) {
    REQUIRE( runANGSDforReadCounts("bamlist.txt", "angsdOut", "10", "angsdOutLog.txt") == 0);
}

TEST_CASE( "Convert ANGSD read counts to unsigned chars for major and minor counts", "[convertCountsToBinary]") {
    REQUIRE( convertANGSDcountsToBinary("angsdOut", "angsdOut.readCounts.binary", 2, 5000) == 0); // 2 individuals, not 1, because generateReadsAndMap actually generates n+1 individuals, since one individual is identical to the reference genome. And 5000 as a max because we don't want to exclude any loci for this test
}

TEST_CASE( "Calculate PWP from the binary representations of the ANGSD readcounts", "[calcPWP]") {
    REQUIRE( calcPWPfromBinaryFile ("angsdOut.readCounts.binary", 74963, 2) == 0);
}
