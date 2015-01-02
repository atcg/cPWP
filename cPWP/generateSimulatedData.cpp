//
//  generateSimulatedData.cpp
//  cPWP
//
//  Created by Evan McCartney-Melstad on 12/31/14.
//  Copyright (c) 2014 Evan McCartney-Melstad. All rights reserved.
//

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "generateSimulatedData.h"




int generateReadsAndMap (int numIndividuals, double mutationRateStepSize, std::string baseErrorRate, std::string libFragmentSize, std::string stdevLibFragmentSize, std::string numReadPairs, std::string readLengths, std::string randomSeed, std::string reference, std::string threads) {
/*// supply 1) number of individuals, 2) the mutation rate steps between them, 3) the base error rate, 4) the average library fragment size, 5) the standard deviation of the library fragment size, 6) the number of read pairs for each individual, 7) read lengths (per read), 8) random number generator seed
 
    So, for instance, if numIndividuals is 5, and mutationRateSteps is 0.01, then there will be four individuals generated with mutation rate steps of 0.01, 0.02, 0.03, and 0.04 away from the reference genome. Individual 0 is created as being identical to the reference genome (mutation rate parameter set to 0.00).
*/

    std::ofstream bamsFile;
    bamsFile.open("bamlist.txt", std::ios::out | std::ios::app);
    
    std::cout << "**********\nChecking if processor is available...";
    if (system(NULL)) puts ("Ok");
    else exit (EXIT_FAILURE);
    
    std::vector<std::string> wgsimCommands;
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
