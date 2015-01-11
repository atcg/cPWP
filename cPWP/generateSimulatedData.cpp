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
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <random>
#include <algorithm>
#include <cstring>


int createReferenceGenome (int totalBases, double gcContent, std::string genomeOutFile) {
    std::ofstream genomeOut(genomeOutFile);
    genomeOut << ">Fake_scaffold0";
    // Set up the random number generator:
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0, 1);
    
    
    std::string genomeString; // Will store the genome
    

    for (int position; position < totalBases; position++) {
        float rando = dis(gen);
        if (rando < gcContent) {
            if (dis(gen) < 0.50) {
                genomeString += "C";
            } else {
                genomeString += "G";
            }
        } else {
            if (dis(gen) < 0.50) {
                genomeString += "A";
            } else {
                genomeString += "T";
            }
        }
    }
    
    int baseCounter = 0;
    for(char& singleBase : genomeString) {
        if (baseCounter % 80 == 0) {
            genomeOut << "\n";
        }
        genomeOut << singleBase;
        baseCounter++;
    }
    genomeOut << std::endl;
    
    // Now index the genome for bwa
    std::cout << "**********\nChecking if processor is available to run bwa index...";
    if (system(NULL)) puts ("OK");
    else exit (EXIT_FAILURE);
    std::string bwaIndexCommand = "bwa index " + genomeOutFile;
    if (system((bwaIndexCommand).c_str()) != 0) {
        std::cout << "**********\nFailure running the following command: " << bwaIndexCommand << "\n**********\n";
        exit(EXIT_FAILURE);
    } else {
        std::cout << "**********\nExecuted the following command: " << bwaIndexCommand << "\n**********\n";
    }
    return 0;
}


int createMutatedGenome (std::string reference, std::string mutatedReferenceFile, float percDivergent) {
    // First put the entire genome string into a single line
    std::ifstream referenceFile(reference);
    std::string line;
    std::string wholeGenome;
    std::ofstream mutGenomeOut(mutatedReferenceFile);
    std::string headerLine;
    int header = 1;
    while(std::getline(referenceFile, line)) {
        if (header == 1) {
            headerLine = line + "_mutated";
            header--;
            continue; // Skip the first line (fasta defline)
        }
        wholeGenome += line;
    }
    mutGenomeOut << headerLine;
    
    /* Since we don't want to put any mutations early in the reference (we want decent coverage at them),
     we'll subtract 100bp from each end of the genome, and evenly space the rest of the mutations after that
     */
    int numMutations = (wholeGenome.length()) * percDivergent;
    int mutationEveryNbp = (wholeGenome.length() - 200) / numMutations;
    
    /* So, let's imaging the reference is 10,000 bp long, and is 1% (0.01) divergent.
     At this stage we'd then have numMutations = 100 mutations, with mutationEveryNbp = (9800/100) every 98bp
    */
    unsigned int baseNum = 0;
    int mutationCounter = 0;
    int lastMutationPosition = 0;
    char abase = 'A';
    char tbase = 'T';
    char cbase = 'C';
    char gbase = 'G';
    std::string aMutBases [] = {"T", "C", "G"};
    std::string tMutBases [] = {"A", "C", "G"};
    std::string cMutBases [] = {"A", "T", "G"};
    std::string gMutBases [] = {"A", "T", "C"};
    while (baseNum < wholeGenome.length()) {
        if (baseNum % 80 == 0) {
            mutGenomeOut << std::endl;
        }
        
        if (baseNum < 100) {
            mutGenomeOut << wholeGenome[baseNum];
            baseNum++;
            continue; // Don't want to put mutations in the first 100 bp
        }
        if (baseNum > wholeGenome.length() - 100) {
            mutGenomeOut << wholeGenome[baseNum];
            baseNum++;
            continue; // Don't want to put any mutations in the last 100 bp
        }
        
        
        if (baseNum % mutationEveryNbp == 0) {
            int randomBaseNum = rand() % 3; // Pick which random base to assign later
            // Mutate that base of the reference and print out
            if (wholeGenome[baseNum] == abase) {
                mutGenomeOut << aMutBases[randomBaseNum];
                baseNum++;
                mutationCounter++;
                lastMutationPosition = baseNum;
                continue;
            } else if (wholeGenome[baseNum] == tbase) {
                mutGenomeOut << tMutBases[randomBaseNum];
                baseNum++;
                mutationCounter++;
                lastMutationPosition = baseNum;
                continue;
            } else if (wholeGenome[baseNum] == cbase) {
                mutGenomeOut << cMutBases[randomBaseNum];
                baseNum++;
                mutationCounter++;
                lastMutationPosition = baseNum;
                continue;
            } else if (wholeGenome[baseNum] == gbase) {
                mutGenomeOut << gMutBases[randomBaseNum];
                baseNum++;
                mutationCounter++;
                lastMutationPosition = baseNum;
                continue;
            } else {
                std::cout << "Base " << wholeGenome[baseNum] << " doesn't equal A, T, C, or G. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
        }
        mutGenomeOut << wholeGenome[baseNum];
        baseNum++;
    }
    mutGenomeOut << std::endl;
    mutGenomeOut.close();
    std::cout << "Inserted " << mutationCounter << " total mutations into " << mutatedReferenceFile << std::endl;
    std::cout << "Position of last mutation: " << lastMutationPosition << std::endl;
    return 0;
}


int generatePerfectReads (std::string reference, unsigned int stagger, unsigned int readLengths, unsigned int fragmentLengths, std::string readPrefix) {
    /* This function generates error-free paired-end sequencing reads from a fasta reference. The reference
     must be a single sequence in fasta format (only one > should be present), with only ATCG characters and
     no gaps. It also assumes that the genome file has 80 bases per line.
    */
    
    std::ifstream referenceFile(reference);
    std::string line;
    std::string wholeGenome;
    int header = 1;
    while(std::getline(referenceFile, line)) {
        if (header == 1) {
            header--;
            continue; // Skip the first line (fasta defline)
        }
        wholeGenome += line;
    }
    std::string R1out = readPrefix + "_R1.fastq";
    std::string R2out = readPrefix + "_R2.fastq";
    std::ofstream R1;
    R1.open(R1out, std::ios::out);
    std::ofstream R2;
    R2.open(R2out, std::ios::out);
    int positionCounter = 1;
    int maxPosition = wholeGenome.length() - fragmentLengths;
    while (positionCounter < maxPosition) {
        R1 << "@" << readPrefix << ":genome_position_left_read=" << positionCounter << "to" << positionCounter + readLengths - 1 << " 1:N:0:AAAAAAAA\n";
        R1 << wholeGenome.substr((positionCounter-1),readLengths) << "\n+\n" << std::string(readLengths, 'I') << "\n";
        R2 << "@" << readPrefix << ":genome_position_left_read=" << positionCounter << "to" << positionCounter + readLengths - 1 << " 2:N:0:AAAAAAAA\n";
        std::string R2Seq = wholeGenome.substr((positionCounter+fragmentLengths-readLengths-1),readLengths);
        reverse(R2Seq.begin(), R2Seq.end());
        char abase = 'A';
        char tbase = 'T';
        char cbase = 'C';
        char gbase = 'G';
        for(unsigned int i = 0; i < R2Seq.length(); i++) {
            if (R2Seq[i] == abase) {
                R2 << "T";
            } else if (R2Seq[i] == tbase) {
                R2 << "A";
            } else if (R2Seq[i] == cbase) {
                R2 << "G";
            } else if (R2Seq[i] == gbase) {
                R2 << "C";
            } else {
                std::cout << "Base pair is not an A, T, C, or G!\n" << std::endl;
                exit(EXIT_FAILURE);
            }
        }
        R2 << "\n+\n" << std::string(readLengths, 'I') << "\n";
        positionCounter = positionCounter + stagger;
    }
    return 0;
}


int mapReads (std::string reference, std::string R1file, std::string R2file, std::string outBam, std::string threads) {
   
    // Make sure the reference is bwa indexed
    std::cout << "**********\nChecking if processor is available to run bwa index...";
    if (system(NULL)) puts ("OK");
    else exit (EXIT_FAILURE);
    std::string bwaIndexCommand = "bwa index " + reference;
    if (system((bwaIndexCommand).c_str()) != 0) {
        std::cout << "**********\nFailure running the following command: " << bwaIndexCommand << "\n**********\n";
        exit(EXIT_FAILURE);
    } else {
        std::cout << "**********\nExecuted the following command: " << bwaIndexCommand << "\n**********\n";
    }
    
    std::string mapCommand = "bwa mem -t " + threads + " " + reference + " " + R1file + " " + R2file + " | samtools view -bS - | samtools sort -T blah -o " + outBam + " -";
    
    std::cout << "**********\nChecking if processor is available to run bwa mem...";
    if (system(NULL)) puts ("OK");
    else exit (EXIT_FAILURE);
    if (system((mapCommand).c_str()) != 0) {
        std::cout << "**********\nFailure running the following command: " << mapCommand << "\n**********\n";
        exit(EXIT_FAILURE);
    } else {
        std::cout << "**********\nExecuted the following command: " << mapCommand << "\n**********\n";
    }
    return 0;
}



int generateReadsAndMap (int numIndividuals, double mutationRateStepSize, std::string libFragmentSize, std::string stdevLibFragmentSize, std::string numReadPairs, std::string readLengths, std::string randomSeed, std::string reference, std::string threads) {
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
            /* pirs simulate implementation
            std::string pirsSimulateCommandToRun = "pirs simulate " + reference + " -l " + readLengths + " -x " + depth + " -m " + libFragmentSize + " -v " + stdevLibFragmentSize + " --no-substitution-errors --no-indel-errors --no-gc-content-bias -o " + indReadsPrefix + " >" + pirsSimSTDOUT + " 2>" + pirsSimSTDERR;
            if (system((pirsSimulateCommandToRun).c_str()) != 0) {
                std::cout << "**********\nFailure running the following command: " << pirsSimulateCommandToRun << "\n**********\n";
                exit(EXIT_FAILURE);
            } else {
                std::cout << "**********\nExecuted the following command: " << pirsSimulateCommandToRun << "\n**********\n";
            }
            */
            std::string R1out = indName + "_R1.fastq";
            std::string R2out = indName + "_R2.fastq";
            
            std::string wgsimCommandToRun = "wgsim -N " + numReadPairs + " -r 0 -R 0.00 -X 0.00 -d " + libFragmentSize + " -s " + stdevLibFragmentSize +  " -1 " + readLengths + " -2 " + readLengths + " -S " + randomSeed + " -e 0.00 " + reference + " " + R1out + " " + R2out + " > ind0_polymorphisms"; // No indels, no probability of indel extension, no base call error rates
            
            if (system((wgsimCommandToRun).c_str()) != 0) {
                std::cout << "**********\nFailure running the following command: " << wgsimCommandToRun << "\n**********\n";
                exit(EXIT_FAILURE);
            } else {
                std::cout << "**********\nExecuted the following command: " << wgsimCommandToRun << "\n**********" << std::endl;
            }
        } else {
            /* pirs implementation
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
            */
            
            double mutRate = pirsInd * mutationRateStepSize;
            // Covert the mutation rate into a string for the system command
            std::ostringstream mutStrs;
            mutStrs << mutRate;
            std::string mutRateString = mutStrs.str();
            
            std::cout << "**********\nGenerating sequence reads for individual " << indName << " using a mutation rate of " << mutRateString << " from the reference genome\n**********\n";
            
            // Get the number of the individual as a string

            
            // Generate the output file names as strings
            std::string R1out = indName + "_R1.fastq";
            std::string R2out = indName + "_R2.fastq";
            std::string polymorphismFile = indName + "_polymorphisms.txt";
            
            // Generate the wgsim command and then run it using a system call
            std::string wgsimCommandToRun = "wgsim -N " + numReadPairs + " -r " + mutRateString + " -R 0.00 -X 0.00 -d " + libFragmentSize + " -s " + stdevLibFragmentSize +  " -1 " + readLengths + " -2 " + readLengths + " -S " + randomSeed + " -e0 " + reference + " " + R1out + " " + R2out + " > " + polymorphismFile; // No indels, no probability of indel extension, no base call error rates
            if (system((wgsimCommandToRun).c_str()) != 0) {
                std::cout << "**********\nFailure running the following command: " << wgsimCommandToRun << "\n**********\n";
                exit(EXIT_FAILURE);
            } else {
                std::cout << "**********\nExecuted the following command: " << wgsimCommandToRun << "\n**********\n";
                
            }

            
            
        }
        // Generate the bwa mem command and then run it using a system call
        std::string R1 = indName + "_R1.fastq";
        std::string R2 = indName + "_R2.fastq";
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


int createHeterozygousGenome(std::string firstFile, std::string secondFile, std::string outputGenome) {
    std::ifstream chrom1(firstFile);
    std::ifstream chrom1a(secondFile);
    std::ofstream write(outputGenome);

    std::string line;
    std::string line2;
    while ( std::getline ( chrom1, line, '\n' ) )
    {
        write << line << std::endl;
    }
    while ( getline ( chrom1a, line2, '\n' ) )
    {
        write << line2 << std::endl;
    }
    chrom1.close();
    chrom1a.close();
    write.close();
    return 0;
}





