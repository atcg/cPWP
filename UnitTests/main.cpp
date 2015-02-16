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



TEST_CASE( "Generate sequence reads 5", "[perfectReads]") {
    REQUIRE( generatePerfectReads ("simulatedReferenceGenome.fasta", 5, 100, 300, "normalRef5") == 0);
    //generatePerfectReads (std::string reference, unsigned int stagger, unsigned int readLengths, unsigned int fragmentLengths, std::string readPrefix);
}


// We need to create a set of reads that are half as dense as the first set, because when we concatenate the heterozygous chromosome reads to them we don't
// want the coverage to be twice as high.
TEST_CASE( "Generate sequence reads 10", "[perfectReads]") {
    REQUIRE( generatePerfectReads ("simulatedReferenceGenome.fasta", 10, 100, 300, "normalRef10") == 0);
    //generatePerfectReads (std::string reference, unsigned int stagger, unsigned int readLengths, unsigned int fragmentLengths, std::string readPrefix);
}


// Map the reads for the homozygous individual
TEST_CASE( " Mapping first set of reads", "[mapReads]") {
    REQUIRE( mapReads("simulatedReferenceGenome.fasta", "normalRef5_R1.fastq", "normalRef5_R2.fastq", "normal.bam", "5") == 0);
    //mapReads (std::string reference, std::string R1file, std::string R2file, std::string outBam, std::string threads)
}


// Generate the reads from the mutated reference chromosome
TEST_CASE( "Generate sequence reads 2", "[perfectReads2]") {
    REQUIRE( generatePerfectReads ("simulatedReferenceGenomeMutated.fasta", 10, 100, 300, "mutatedRef") == 0);
    //generatePerfectReads (std::string reference, unsigned int stagger, unsigned int readLengths, unsigned int fragmentLengths, std::string readPrefix);
}


// Concatenate the het-reference reads to the homozygous reads reference for R1
TEST_CASE( "Create heterozygous R1", "[createHetR1]") {
    REQUIRE( createHeterozygousGenome("normalRef10_R1.fastq", "mutatedRef_R1.fastq", "hetRef_R1.fastq") == 0);
    //createHeterozygousGenome(std::string firstFile, std::string secondFile, std::string outputGenome)
}

// Concatenate the het-reference reads to the homozygous reads reference for R2
TEST_CASE( "Create heterozygous R2", "[createHetR1]") {
    REQUIRE( createHeterozygousGenome("normalRef10_R2.fastq", "mutatedRef_R2.fastq", "hetRef_R2.fastq") == 0);
    //createHeterozygousGenome(std::string firstFile, std::string secondFile, std::string outputGenome)
}

// Map the reads for the heterozygous individual
TEST_CASE( " Mapping het reads", "[mapReads2]") {
    REQUIRE( mapReads("simulatedReferenceGenome.fasta", "hetRef_R1.fastq", "hetRef_R2.fastq", "het.bam", "5") == 0);
    //mapReads (std::string reference, std::string R1file, std::string R2file, std::string outBam, std::string threads)
}


// Use angsd to generate the read counts for the two individuals
TEST_CASE( "Run ANGSD on simulated reads", "[runANGSD]" ) {
    REQUIRE( runANGSDforReadCounts("bamlist.txt", "angsdOut", "5", "angsdOutLog.txt") == 0);
    //runANGSDforReadCounts (std::string bamlist, std::string angsdPrefix, std::string nThreads, std::string angsdOutputLog)
}

// Convert read counts to binary
TEST_CASE( "Convert ANGSD read counts to unsigned chars for major and minor counts", "[convertCountsToBinary]") {
    REQUIRE( convertANGSDcountsToBinary("angsdOut", "angsdOut.readCounts.binary", 2, 250) == 0);
    //convertANGSDcountsToBinary(std::string angsdPrefix, std::string binaryOutputFileName, int numIndividuals, int readDepthMax)
}

// Run the PWP calculation
TEST_CASE( "Calculate PWP from the binary representations of the ANGSD readcounts", "[calcPWP]") {
    REQUIRE( calcPWPfromBinaryFile("angsdOut.readCounts.binary", 999990, 2, "testingOut.pwp", 1000, 5) == 0);
    //int calcPWPfromBinaryFile (std::string binaryFile, unsigned long long int numLoci, const int numIndividuals, std::string outFile, int lociChunkSize, int numThreads=30);
}







/* 
 Covariance unit tests are below. They require different reference genomes to be used. In the example above, the homozygous reference genome
 has no variance in the major allele frequencies (they're always 1). The genomes created by this example have 
 
 
 
 
 */

TEST_CASE( "Generate reference genome for the covariance tests", "[generateReferenceCovar]") {
    REQUIRE( createReferenceGenome(10000000, 0.42668722, "simulatedReferenceGenome10mil.fasta") == 0 );
    //createReferenceGenome (int totalBases, double gcContent, std::string genomeOutFile)
}


TEST_CASE ( "Mutate a reference genome for the covariance test", "[mutateRefGenomeCovar]") {
    REQUIRE( createMutatedGenomesForCovar("simulatedReferenceGenome10mil.fasta", "simulatedReferenceGenomeMutatedRef1.fasta", "simulatedReferenceGenomeMutatedRef2.fasta", 0.01) == 0);
    //createMutatedGenomesForCovar (std::string reference, std::string mutatedReferenceFile1, std::string mutatedReferenceFile2, float percDivergent)
}


TEST_CASE( "Generate nonmutated chromosome covariance sequence reads for ind1 and ind2", "[perfectReadsCovarNoMut]") {
    REQUIRE( generatePerfectReads ("simulatedReferenceGenome10mil.fasta", 10, 100, 300, "covarRef") == 0);
    //generatePerfectReads (std::string reference, unsigned int stagger, unsigned int readLengths, unsigned int fragmentLengths, std::string readPrefix);
}



TEST_CASE( "Generate mutated chromosome covariance sequence reads for ind1", "[perfectReadsCovarInd1]") {
    REQUIRE( generatePerfectReads ("simulatedReferenceGenomeMutatedRef1.fasta", 10, 100, 300, "covarMutRef1") == 0);
    //generatePerfectReads (std::string reference, unsigned int stagger, unsigned int readLengths, unsigned int fragmentLengths, std::string readPrefix);
}


TEST_CASE( "Generate mutated chromosome covariance sequence reads for ind2", "[perfectReadsCovarInd2]") {
    REQUIRE( generatePerfectReads ("simulatedReferenceGenomeMutatedRef2.fasta", 10, 100, 300, "covarMutRef2") == 0);
    //generatePerfectReads (std::string reference, unsigned int stagger, unsigned int readLengths, unsigned int fragmentLengths, std::string readPrefix);
}



// Combine the mutated and non-mutated chromosome reads for each individual
// Concatenate the het-reference reads to the homozygous reads reference for R1
TEST_CASE( "Create heterozygous R1covar1", "[createHetR1covar1]") {
    REQUIRE( createHeterozygousGenome("covarRef_R1.fastq", "covarMutRef1_R1.fastq", "covarRef1_R1.fastq") == 0);
    //createHeterozygousGenome(std::string firstFile, std::string secondFile, std::string outputGenome)
}

// Concatenate the het-reference reads to the homozygous reads reference for R2
TEST_CASE( "Create heterozygous R2covar1", "[createHetR2covar1]") {
    REQUIRE( createHeterozygousGenome("covarRef_R2.fastq", "covarMutRef1_R2.fastq", "covarRef1_R2.fastq") == 0);
    //createHeterozygousGenome(std::string firstFile, std::string secondFile, std::string outputGenome)
}

// Concatenate the het-reference reads to the homozygous reads reference for R1
TEST_CASE( "Create heterozygous R1covar2", "[createHetR1covar2]") {
    REQUIRE( createHeterozygousGenome("covarRef_R1.fastq", "covarMutRef2_R1.fastq", "covarRef2_R1.fastq") == 0);
    //createHeterozygousGenome(std::string firstFile, std::string secondFile, std::string outputGenome)
}

// Concatenate the het-reference reads to the homozygous reads reference for R2
TEST_CASE( "Create heterozygous R2covar2", "[createHetR2covar2]") {
    REQUIRE( createHeterozygousGenome("covarRef_R2.fastq", "covarMutRef2_R2.fastq", "covarRef2_R2.fastq") == 0);
    //createHeterozygousGenome(std::string firstFile, std::string secondFile, std::string outputGenome)
}


// Map the reads for individual 1
TEST_CASE( " Mapping reads for covariance ind1", "[mapReadsCovar1]") {
    REQUIRE( mapReads("simulatedReferenceGenome10mil.fasta", "covarRef1_R1.fastq", "covarRef1_R2.fastq", "covarInd1.bam", "5") == 0);
    //mapReads (std::string reference, std::string R1file, std::string R2file, std::string outBam, std::string threads)
}


// Map the reads for individual 2
TEST_CASE( " Mapping reads for covariance ind2", "[mapReadsCovar2]") {
    REQUIRE( mapReads("simulatedReferenceGenome10mil.fasta", "covarRef2_R1.fastq", "covarRef2_R2.fastq", "covarInd2.bam", "5") == 0);
    //mapReads (std::string reference, std::string R1file, std::string R2file, std::string outBam, std::string threads)
}



// Use angsd to generate the read counts for the two individuals
TEST_CASE( "Run ANGSD on simulated reads for the covariance test", "[runANGSDCovar]" ) {
    REQUIRE( runANGSDforReadCounts("covarBamlist.txt", "covarAngsdOut", "5", "covarAngsdOutLog.txt") == 0);
    //runANGSDforReadCounts (std::string bamlist, std::string angsdPrefix, std::string nThreads, std::string angsdOutputLog)
}

// Convert read counts to binary
TEST_CASE( "Convert ANGSD read counts to unsigned chars for major and minor counts for the covariance test", "[convertCountsToBinaryCovar]") {
    REQUIRE( convertANGSDcountsToBinary("covarAngsdOut", "covarAngsdOut.readCounts.binary", 2, 250) == 0);
    //convertANGSDcountsToBinary(std::string angsdPrefix, std::string binaryOutputFileName, int numIndividuals, int readDepthMax)
}


TEST_CASE( "Calculate covariances from the same representations of the ANGSD readcounts", "[calcCovar]") {
    REQUIRE(  calcCOVARfromBinaryFile("covarAngsdOut.readCounts.binary", 8999990, 2, "testingOut.covar", 100000, 5)== 0);
    //calcCOVARfromBinaryFile (std::string binaryFile, long long int numLoci, const int numIndividuals, std::string outFile, int lociChunkSize, const int numThreads)
}


