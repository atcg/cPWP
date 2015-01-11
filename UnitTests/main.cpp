//
//  main.cpp
//  Unit Tests
//
//  Created by Evan McCartney-Melstad on 12/31/14.
//  Copyright (c) 2014 Evan McCartney-Melstad. All rights reserved.
//

#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "../cPWP/generateSimulatedData.h"
#include "../cPWP/bamsToBin.h"
#include "../cPWP/calcPWP.h"
#include "catch.hpp"
#include <string>
#include <iostream>
#include <fstream>


/*
TEST_CASE( "Simulated reads are generated", "[generateReads]" ) {
    // Make sure that the read simulation finishes
    REQUIRE( generateReadsAndMap(1, 0.01, "0.0", "300", "50", "1000000", "100", "1234", "scaffold_0.fasta") == 0);
}
 */

TEST_CASE( "Generate reference genome for simulation tests", "[generateReference]") {
    REQUIRE( createReferenceGenome(10000000, 0.42668722, "simulatedReferenceGenome.fasta") == 0 );
}

TEST_CASE ( "Mutate a reference genome", "[mutateRefGenome]") {
    REQUIRE( createMutatedGenome("simulatedReferenceGenome.fasta", "simulatedReferenceGenomeMutated.fasta", 0.02) == 0);
}


TEST_CASE( "Generate sequence reads", "[perfectReads]") {
    REQUIRE( generatePerfectReads ("simulatedReferenceGenome.fasta", 10, 100, 300, "normalRef") == 0);
    //generatePerfectReads (std::string reference, unsigned int stagger, unsigned int readLengths, unsigned int fragmentLengths, std::string readPrefix);
}


TEST_CASE( " Mapping first set of reads", "[mapReads]") {
    REQUIRE( mapReads("simulatedReferenceGenome.fasta", "normalRef_R1.fastq", "normalRef_R2.fastq", "normal.bam", "25") == 0);
}


TEST_CASE( "Generate sequence reads 2", "[perfectReads2]") {
    REQUIRE( generatePerfectReads ("simulatedReferenceGenomeMutated.fasta", 10, 100, 300, "mutatedRef") == 0);
}

/*
TEST_CASE( " Mapping second set of reads", "[mapReads2]") {
    REQUIRE( mapReads("simulatedReferenceGenome.fasta", "mutatedRef_R1.fastq", "mutatedRef_R2.fastq", "mutated.bam", "25") == 0);
}
*/


TEST_CASE( "Create heterozygous R1", "[createHet]") {
    REQUIRE( createHeterozygousGenome("normalRef_R1.fastq", "mutatedRef_R1.fastq", "hetRef_R1.fastq") == 0);
}

TEST_CASE( "Create heterozygous R2", "[createHet]") {
    REQUIRE( createHeterozygousGenome("normalRef_R2.fastq", "mutatedRef_R2.fastq", "hetRef_R2.fastq") == 0);
}


TEST_CASE( " Mapping het reads", "[mapReads2]") {
    REQUIRE( mapReads("simulatedReferenceGenome.fasta", "hetRef_R1.fastq", "hetRef_R2.fastq", "het.bam", "25") == 0);
}



TEST_CASE( "Run ANGSD on simulated reads", "[runANGSD]" ) {
    REQUIRE( runANGSDforReadCounts("bamlist.txt", "angsdOut", "25", "angsdOutLog.txt") == 0);
}


TEST_CASE( "Convert ANGSD read counts to unsigned chars for major and minor counts", "[convertCountsToBinary]") {
    REQUIRE( convertANGSDcountsToBinary("angsdOut", "angsdOut.readCounts.binary", 2, 5000) == 0); // 5000 as a max because we don't want to exclude any loci for this test
}

TEST_CASE( "Calculate PWP from the binary representations of the ANGSD readcounts", "[calcPWP]") {
    REQUIRE( calcPWPfromBinaryFile ("angsdOut.readCounts.binary", 0, 2, "testingOut.pwp", 30) == 0);
}


