//
//  generateSimulatedData.cpp
//  cPWP
//
//  Created by Evan McCartney-Melstad on 12/31/14.
//  Copyright (c) 2014 Evan McCartney-Melstad. All rights reserved.
//

#include "generateSimulatedData.h"
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <random>


int createReferenceGenome (int totalBases, double gcContent) {
    std::ofstream genomeOut("simulatedReferenceGenome.fasta");
    genomeOut << ">Fake_scaffold0";
    static const char bases[] = "ATCG";
    // Set up the random number generator:
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0, 1);
    
    
    std::string genomeString; // Will store the genome
    
    char pickBase() {
        float rando = dis(gen);
        if (rando < gcContent) {
            if (dis(gen) < 0.50) {
                return "C";
            } else {
                return "G";
            }
        } else {
            if (dis(gen) < 0.50) {
                return "A";
            } else {
                return "T";
            }
        }
    }
    
    for (int position; position < totalBases; position++) {
        genomeString += pickBase();
    }
    
    int baseCounter = 0;
    for(char& singleBase : genomeString) {
        if (baseCounter % 80 == 0) {
            genomeOut << "\n";
        }
        genomeOut << singleBase;
        baseCounter++;
    }
    return 0;
}



int generateReadsAndMap (int numIndividuals, double mutationRateStepSize, std::string libFragmentSize, std::string stdevLibFragmentSize, std::string depth, std::string readLengths, std::string randomSeed, std::string reference, std::string threads) {
/*// supply 1) number of individuals, 2) the mutation rate steps between them, 3) the base error rate, 4) the average library fragment size, 5) the standard deviation of the library fragment size, 6) the number of read pairs for each individual, 7) read lengths (per read), 8) random number generator seed
 
    So, for instance, if numIndividuals is 5, and mutationRateSteps is 0.01, then there will be four individuals generated with mutation rate steps of 0.01, 0.02, 0.03, and 0.04 away from the reference genome. Individual 0 is created as being identical to the reference genome (mutation rate parameter set to 0.00).
*/

    std::ofstream bamsFile;
    bamsFile.open("bamlist.txt", std::ios::out | std::ios::app);
    
    
    std::cout << "**********\nChecking if processor is available to run pirs...";
    if (system(NULL)) puts ("OK");
    else exit (EXIT_FAILURE);
    
    int pirsInd = 0;
    while(pirsInd <= numIndividuals) {
        double mutRate = pirsInd * mutationRateStepSize;
        // Covert the mutation rate into a string for the system command
        std::ostringstream mutStr;
        mutStr << mutRate;
        std::string mutRateString = mutStr.str();
        
        std::cout << "**********\nGenerating reference genome for individual " << pirsInd << " using a mutation rate of " << mutRateString << " from the reference genome\n**********\n";
        
        // Since we're identifying individuals by the looping variable (an int), we need to convert that to a string to use it in the file names
        std::ostringstream pirsIndSS;
        std::ostringstream pirsIndNumSS;
        pirsIndSS << "ind" << pirsInd;
        pirsIndNumSS << pirsInd;
        std::string pirsIndNum = pirsIndNumSS.str();
        std::string indName = pirsIndSS.str();
        std::string pirsGenomeSTDOUT = indName + "_genome.stdout";
        std::string pirsGenomeSTDERR = indName + "_genome.stderr";
    
        // Simulate the other strand of the mutated reference genome. On the first individual it should be identical to the reference (because mutStr = 0 * mutationRateStepSize
        // For this individual, we are not simulating any polymorphisms. So instead of the diploid pirs simulate command, we can use a haploid simulate command on the original reference
        std::string pirsSimSTDOUT = indName + "_reads.stdout";
        std::string pirsSimSTDERR = indName + "_reads.stderr";
        std::string indReadsPrefix = indName + "_reads";
        
        if (pirsInd == 0) {
            std::string pirsSimulateCommandToRun = "pirs simulate " + reference + " -l " + readLengths + " -x " + depth + " -m " + libFragmentSize + " -v " + stdevLibFragmentSize + " --no-substitution-errors --no-indel-errors --no-gc-content-bias -o " + indReadsPrefix + " >" + pirsSimSTDOUT + " 2>" + pirsSimSTDERR;
            if (system((pirsSimulateCommandToRun).c_str()) != 0) {
                std::cout << "**********\nFailure running the following command: " << pirsSimulateCommandToRun << "\n**********\n";
                exit(EXIT_FAILURE);
            } else {
                std::cout << "**********\nExecuted the following command: " << pirsSimulateCommandToRun << "\n**********\n";
            }
        } else {
            std::string pirsCommandToRun = "pirs diploid -s " + mutRateString + " -d 0.00 -v 0.00 -S 1234 -o " + indName + " " + reference + " >" + pirsGenomeSTDOUT + " 2>" + pirsGenomeSTDERR;
            if (system((pirsCommandToRun).c_str()) != 0) {
                std::cout << "**********\nFailure running the following command: " << pirsCommandToRun << "\n**********\n";
                exit(EXIT_FAILURE);
            } else {
                std::cout << "**********\nExecuted the following command: " << pirsCommandToRun << "\n**********\n";
            }
            // The following file is output by the pirs diploid command:
            std::string mutatedChromosome = indName + ".snp.fa";
            
            // After generating the mutated strand for each individual, we simulate reads from both strands for each individual. Parameterize the
            std::string pirsSimSTDOUT = indName + "_reads.stdout";
            std::string pirsSimSTDERR = indName + "_reads.stderr";
            std::string indReadsPrefix = indName + "_reads";
            std::string pirsSimulateCommandToRun = "pirs simulate --diploid " + reference + " " + mutatedChromosome + " -l " + readLengths + " -x " + depth + " -m " + libFragmentSize + " -v " + stdevLibFragmentSize + " --no-substitution-errors --no-indel-errors --no-gc-content-bias -o " + indReadsPrefix + " >" + pirsSimSTDOUT + " 2>" + pirsSimSTDERR;
            if (system((pirsSimulateCommandToRun).c_str()) != 0) {
                std::cout << "**********\nFailure running the following command: " << pirsSimulateCommandToRun << "\n**********\n";
                exit(EXIT_FAILURE);
            } else {
                std::cout << "**********\nExecuted the following command: " << pirsSimulateCommandToRun << "\n**********\n";
            }
        }
        // Generate the bwa mem command and then run it using a system call
        std::string R1 = indReadsPrefix + "_" + readLengths + "_" + libFragmentSize + "_1.fq";
        std::string R2 = indReadsPrefix + "_" + readLengths + "_" + libFragmentSize + "_2.fq";
        std::string bamOut = "ind" + pirsIndNum + ".bam";
        std::string bwaCommandToRun = "bwa mem -t " + threads + " " + reference + " " + R1 + " " + R2 + " | samtools view -bS - | samtools sort -T temp -o " + bamOut + " -";
        if (system((bwaCommandToRun).c_str()) != 0) {
            std::cout << "**********\nFailure running the following command: " << bwaCommandToRun << "\n**********\n";
            exit(EXIT_FAILURE);
        } else {
            std::cout << "**********\nExecuted the following command: " << bwaCommandToRun << "\n**********\n";
        }
        bamsFile << bamOut << std::endl;

        
        pirsInd++; // Move on to the next individual
    }
    
    /* Implementation with wgsim
    std::cout << "**********\nChecking if processor is available to run wgsim...";
    if (system(NULL)) puts ("Ok");
    else exit (EXIT_FAILURE);
    
    int step = 0;
    while(step <= numIndividuals) {
        double mutRate = step * mutationRateStepSize;
        // Covert the mutation rate into a string for the system command
        std::ostringstream mutStrs;
        mutStrs << mutRate;
        std::string mutRateString = mutStrs.str();
        
        std::cout << "**********\nGenerating sequence reads for individual " << step << " using a mutation rate of " << mutRateString << " from the reference genome\n**********\n";
        
        // Get the number of the individual as a string
        std::ostringstream stepString;
        stepString << step;
        std::string ind = stepString.str();
        
        // Generate the output file names as strings
        std::string R1out = "ind" + ind + "_R1.fastq";
        std::string R2out = "ind" + ind + "_R2.fastq";
        std::string polymorphismFile = "ind" + ind + "_polymorphisms.txt";
        
        // Generate the wgsim command and then run it using a system call
        std::string wgsimCommandToRun = "wgsim -N " + numReadPairs + " -r " + mutRateString + " -R 0.00 -X 0.00 -d " + libFragmentSize + " -s " + stdevLibFragmentSize +  " -1 " + readLengths + " -2 " + readLengths + " -S " + randomSeed + " -e0 " + reference + " " + R1out + " " + R2out + " > " + polymorphismFile; // No indels, no probability of indel extension, no base call error rates
        if (system((wgsimCommandToRun).c_str()) != 0) {
            std::cout << "**********\nFailure running the following command: " << wgsimCommandToRun << "\n**********\n";
            exit(EXIT_FAILURE);
        } else {
            std::cout << "**********\nExecuted the following command: " << wgsimCommandToRun << "\n**********\n";

        }
        
        // Generate the bwa mem command and then run it using a system call
        std::string bamOut = "ind" + ind + ".bam";
        std::string bwaCommandToRun = "bwa mem -t " + threads + " " + reference + " " + R1out + " " + R2out + " | samtools view -bS - | samtools sort -T temp -o " + bamOut + " -";
        if (system((bwaCommandToRun).c_str()) != 0) {
            std::cout << "**********\nFailure running the following command: " << bwaCommandToRun << "\n**********\n";
            exit(EXIT_FAILURE);
        } else {
            std::cout << "**********\nExecuted the following command: " << bwaCommandToRun << "\n**********\n";
        }
        bamsFile << bamOut << std::endl;
        step++; // Move on to the next individual
    }
    */
    bamsFile.close();
    return 0;
}


/*
 /home/evan/bin/angsd0.613/angsd -bam bamlist272.txt -out 272torts_allCounts_minmapq20minq30 -uniqueOnly 1 -only_proper_pairs 1 -remove_bads 1 -doCounts 1 -nThreads 8 -dumpCounts 4 -doMaf 1 -doMajorMinor 2 -GL 2 -minMapQ 20 -minQ 30 > 272tortsangsdMapQ30_allSites.log 2>&1
 
*/


/*
 wgsim -N 10000000 -r 0.01 -R 0.00 -X 0.00 -1 100 -2 100 -S11 -e0 lgv2.fasta r1_noerr1percDivergent_10mil.fq r2_noerr1percDivergent_10mil.fq > 1perc_polymorphisms.txt
 wgsim -N 10000000 -r 0.00 -1 100 -2 100 -S11 -e0 lgv2.fasta r1_noerrnomut_10mil.fq r2_noerrnomut_10mil.fq
 
 bwa mem -t 10 lgv2.fasta r1_noerrnomut_10mil.fq r2_noerrnomut_10mil.fq
 bwa mem -t 10 lgv2.fasta r1_noerr1percDivergent_10mil.fq r2_noerr1percDivergent_10mil.fq
 
 samtools sort -T noMut -o 1XcoverageNoErrNoMutSorted.bam 1XcoverageNoErrNoMut.bam
 samtools sort -T mut -o 1XcoverageNoErr1PercentDivergentSorted.bam 1XcoverageNoErr1PercDivergent.bam
 
 echo -e "1XcoverageNoErrNoMutSorted.bam\n1XcoverageNoErr1PercentDivergentSorted.bam" > bamlist2sim.txt
 
 /home/evan/bin/angsd0.613/angsd -bam bamlist2sim.txt -out 2tortsim1perDiv -uniqueOnly 1 -only_proper_pairs 1 -remove_bads 1 -doCounts 1 -nThreads 8 -dumpCounts 4 -doMaf 1 -doMajorMinor 2 -GL 2 -minMapQ 20 -minQ 30 > 2tortANGSD.log 2>&1 &
 
*/





/*
 Usage:   wgsim [options] <in.ref.fa> <out.read1.fq> <out.read2.fq>
 
 Options:
 -e FLOAT      base error rate [0.000]
 -d INT        outer distance between the two ends [500]
 -s INT        standard deviation [50]
 -N INT        number of read pairs [1000000]
 -1 INT        length of the first read [70]
 -2 INT        length of the second read [70]
 -r FLOAT      rate of mutations [0.0010]
 -R FLOAT      fraction of indels [0.15]
 -X FLOAT      probability an indel is extended [0.30]
 -S INT        seed for random generator [-1]
 -h            haplotype mode
 
 */

