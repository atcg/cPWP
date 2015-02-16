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
//#include "../cPWP/calcPWP.h"
#include "../cPWP/calcPWPchunks.h"
#include "../cPWP/calcCOVARchunks.h"
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
    REQUIRE( createReferenceGenome(1000000, 0.42668722, "simulatedReferenceGenome.fasta") == 0 );
    //createReferenceGenome (int totalBases, double gcContent, std::string genomeOutFile) 
}

TEST_CASE ( "Mutate a reference genome", "[mutateRefGenome]") {
    REQUIRE( createMutatedGenome("simulatedReferenceGenome.fasta", "simulatedReferenceGenomeMutated.fasta", 0.01) == 0);
    //createMutatedGenome (std::string reference, std::string mutatedReferenceFile, float percDivergent) 
}


TEST_CASE ( "Mutate anoter reference genome", "[mutateRefGenome2]") {
    REQUIRE( createMutatedGenome("simulatedReferenceGenome.fasta", "simulatedReferenceGenomeMutated2.fasta", 0.02) == 0);
    //createMutatedGenome (std::string reference, std::string mutatedReferenceFile, float percDivergent)
}


TEST_CASE( "Generate sequence reads 5", "[perfectReads]") {
    REQUIRE( generatePerfectReads ("simulatedReferenceGenome.fasta", 5, 100, 300, "normalRef5") == 0);
    //generatePerfectReads (std::string reference, unsigned int stagger, unsigned int readLengths, unsigned int fragmentLengths, std::string readPrefix);
}


TEST_CASE( "Generate sequence reads 10", "[perfectReads]") {
    REQUIRE( generatePerfectReads ("simulatedReferenceGenome.fasta", 10, 100, 300, "normalRef10") == 0);
    //generatePerfectReads (std::string reference, unsigned int stagger, unsigned int readLengths, unsigned int fragmentLengths, std::string readPrefix);
}


TEST_CASE( " Mapping first set of reads", "[mapReads]") {
    REQUIRE( mapReads("simulatedReferenceGenome.fasta", "normalRef5_R1.fastq", "normalRef5_R2.fastq", "normal.bam", "5") == 0);
    //mapReads (std::string reference, std::string R1file, std::string R2file, std::string outBam, std::string threads)
}



TEST_CASE( "Generate sequence reads 2", "[perfectReads2]") {
    REQUIRE( generatePerfectReads ("simulatedReferenceGenomeMutated.fasta", 10, 100, 300, "mutatedRef") == 0);
    //generatePerfectReads (std::string reference, unsigned int stagger, unsigned int readLengths, unsigned int fragmentLengths, std::string readPrefix);
}

TEST_CASE( "Generate sequence reads 2b", "[perfectReads2b]") {
    REQUIRE( generatePerfectReads ("simulatedReferenceGenomeMutated2.fasta", 10, 100, 300, "mutatedRef2") == 0);
    //generatePerfectReads (std::string reference, unsigned int stagger, unsigned int readLengths, unsigned int fragmentLengths, std::string readPrefix);
}





TEST_CASE( "Create heterozygous R1", "[createHetR1]") {
    REQUIRE( createHeterozygousGenome("normalRef10_R1.fastq", "mutatedRef_R1.fastq", "hetRef_R1.fastq") == 0);
    //createHeterozygousGenome(std::string firstFile, std::string secondFile, std::string outputGenome)
}

TEST_CASE( "Create heterozygous R2", "[createHetR1]") {
    REQUIRE( createHeterozygousGenome("normalRef10_R2.fastq", "mutatedRef_R2.fastq", "hetRef_R2.fastq") == 0);
    //createHeterozygousGenome(std::string firstFile, std::string secondFile, std::string outputGenome)
}



TEST_CASE( "Create heterozygous R1b", "[createHetR1b]") {
    REQUIRE( createHeterozygousGenome("normalRef10_R1.fastq", "mutatedRef2_R1.fastq", "hetRef2_R1.fastq") == 0);
    //createHeterozygousGenome(std::string firstFile, std::string secondFile, std::string outputGenome)
}

TEST_CASE( "Create heterozygous R2b", "[createHetR2b]") {
    REQUIRE( createHeterozygousGenome("normalRef10_R2.fastq", "mutatedRef2_R2.fastq", "hetRef2_R2.fastq") == 0);
    //createHeterozygousGenome(std::string firstFile, std::string secondFile, std::string outputGenome)
}






TEST_CASE( " Mapping het reads", "[mapReads2]") {
    REQUIRE( mapReads("simulatedReferenceGenome.fasta", "hetRef_R1.fastq", "hetRef_R2.fastq", "het.bam", "5") == 0);
    //mapReads (std::string reference, std::string R1file, std::string R2file, std::string outBam, std::string threads)
}


TEST_CASE( " Mapping het reads 2", "[mapReads2b]") {
    REQUIRE( mapReads("simulatedReferenceGenome.fasta", "hetRef2_R1.fastq", "hetRef2_R2.fastq", "het2.bam", "5") == 0);
    //mapReads (std::string reference, std::string R1file, std::string R2file, std::string outBam, std::string threads)
}



TEST_CASE( "Run ANGSD on simulated reads", "[runANGSD]" ) {
    REQUIRE( runANGSDforReadCounts("bamlist.txt", "angsdOut", "5", "angsdOutLog.txt") == 0);
    //runANGSDforReadCounts (std::string bamlist, std::string angsdPrefix, std::string nThreads, std::string angsdOutputLog)
}


TEST_CASE( "Convert ANGSD read counts to unsigned chars for major and minor counts", "[convertCountsToBinary]") {
    REQUIRE( convertANGSDcountsToBinary("angsdOut", "angsdOut.readCounts.binary", 3, 250) == 0);
    //convertANGSDcountsToBinary(std::string angsdPrefix, std::string binaryOutputFileName, int numIndividuals, int readDepthMax)
}

TEST_CASE( "Calculate PWP from the binary representations of the ANGSD readcounts", "[calcPWP]") {
    REQUIRE( calcPWPfromBinaryFile("angsdOut.readCounts.binary", 999990, 3, "testingOut.pwp", 1000, 5) == 0);
    //int calcPWPfromBinaryFile (std::string binaryFile, unsigned long long int numLoci, const int numIndividuals, std::string outFile, int lociChunkSize, int numThreads=30);
}


TEST_CASE( "Calculate covariances from the same representations of the ANGSD readcounts", "[calcCovar]") {
    REQUIRE(  calcCOVARfromBinaryFile("angsdOut.readCounts.binary", 999990, 3, "testingOut.covar", 1000, 5)== 0);
    //calcCOVARfromBinaryFile (std::string binaryFile, long long int numLoci, const int numIndividuals, std::string outFile, int lociChunkSize, const int numThreads)
}





















