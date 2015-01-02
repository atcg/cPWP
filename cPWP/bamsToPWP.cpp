//
//  readsToPWP.cpp
//  cPWP
//
//  Created by Evan McCartney-Melstad on 12/31/14.
//  Copyright (c) 2014 Evan McCartney-Melstad. All rights reserved.
//

#include "bamsToPWP.h"
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/algorithm/string.hpp>




int runANGSDforReadCounts (std::string bamlist, std::string angsdPrefix, std::string nThreads, std::string angsdOutputLog) {
    std::cout << "**********\nChecking if processor is available to run ANGSD...**********\n";
    if (system(NULL)) puts ("Ok");
    else exit (EXIT_FAILURE);
    
    std::string angsdCommand = "/home/evan/bin/angsd0.613/angsd -bam " + bamlist + " -out " + angsdPrefix + " -uniqueOnly 1 -only_proper_pairs 1 -remove_bads 1 -doCounts 1 -nThreads " + nThreads + " -dumpCounts 4 -doMaf 1 -doMajorMinor 2 -GL 2 -minMapQ 20 -minQ 30 > " + angsdOutputLog + " 2>&1";
    
    if (system((angsdCommand).c_str()) != 0) {
        std::cout << "**********\nFailure running the following command: " << angsdCommand << "\n**********\n";
        exit(EXIT_FAILURE);
    } else {
        std::cout << "**********\nExecuted the following command: " << angsdCommand << "\n**********\n";
    }
    
    return 0;
}


unsigned char readCountToUnsignedChar(int individual, std::string allele, std::vector<std::string> countsVector) {
    /*    std::cout << "countsVector elements: ";
     for (int i=0; i < countsVector.size(); i++) {
     std::cout <<  countsVector.at(i) << " ";
     }
     std::cout << std::endl;
     */
    
    //    std::cout << "FunctionPassedAllele: " << allele << std::endl;
    
    unsigned char countChar;
    if (allele == "A") {
        //        std::cout << "Allele is an A" << std::endl;
        int countField = individual*4;
        int countInt = atoi(countsVector[countField].c_str());
        countChar = (unsigned char) countInt;
    }
    else if (allele == "C") {
        //	std::cout << "Allele is an C" << std::endl;
        int countField = (individual*4)+1;
        int countInt = atoi(countsVector[countField].c_str());
        countChar = (unsigned char) countInt;
    } else if (allele == "G") {
        //	std::cout << "Allele is a G" << std::endl;
        int countField = (individual*4) + 2;
        int countInt = atoi(countsVector[countField].c_str());
        countChar = (unsigned char) countInt;
    } else if (allele == "T") {
        //	std::cout << "Allele is a T" << std::endl;
        int countField = (individual*4) + 3;
        int countInt = atoi(countsVector[countField].c_str());
        countChar = (unsigned char) countInt;
    } else {
        std::cout << "Allele doesn't look right. The allele fed into readCountToUnsignedChar is: " << allele << ". Exiting now" << std::endl;
        exit(EXIT_FAILURE);
    }
    //    std::cout << std::to_string(countChar) << std::endl;
    return countChar;
}



int convertANGSDcountsToBinary(std::string angsdPrefix, std::string binaryOutputFileName, int numIndividuals, int readDepthMax)
{
    std::string countsFileName = angsdPrefix + ".counts.gz";
    std::string mafsFileName = angsdPrefix + ".mafs.gz";
    std::ifstream countsFile(countsFileName, std::ios_base::in | std::ios_base::binary);
    std::ifstream mafsFile(mafsFileName, std::ios_base::in | std::ios_base::binary);
    std::ofstream outCountsFile(binaryOutputFileName, std::ios_base::out | std::ios_base::binary);
    
    try {
        boost::iostreams::filtering_istream in;
        boost::iostreams::filtering_istream inMafs;
        in.push(boost::iostreams::gzip_decompressor());
        in.push(countsFile);
        inMafs.push(boost::iostreams::gzip_decompressor());
        inMafs.push(mafsFile);
        int counter = 0;
        //int numIndividuals = 3;
        for(std::string str; std::getline(in, str); )
        {
            /* Make sure to get the maf string before we do anything else to keep it in
             sync with the counts string */
            std::string mafStr;
            std::getline(inMafs, mafStr);
            
            // Skip the first line of the counts and mafs file, since it is a header
            if (counter == 0) {
                counter++;
                continue;
            }
            
            // Just for monitoring progress
            if (counter % 1000000 == 0) {
                std::cout << "Processed " << counter << " loci" << std::endl;
            }
            
            /* The first thing we'll do after getting the lines from the gzipped file is to split 
             the tab-separated lines into separate vectors for the counts and mafs. The vector will 
             be recreated for every line of the gzipped counts and mafs files output by angsd (each 
             line represents a locus) */
            std::vector<std::string> countsFields;
            boost::algorithm::split(countsFields, str, boost::is_any_of("\t"), boost::token_compress_on);
            
            std::vector<std::string> mafFields;
            boost::algorithm::split(mafFields, mafStr, boost::is_any_of("\t"), boost::token_compress_on);
            
            
            /* Since we're using a weighted PWP calculation, and since we don't want to overweight
             repetitive regions, we'll throw out any locus that has more than readDepthMax reads
             for any individual at that locus */
            std::vector<int> countsFieldsInts;
            for (size_t i=0; i<countsFields.size(); i++) {
                int countInt = atoi((countsFields[i]).c_str());
                countsFieldsInts.push_back(countInt);
            }
            
            
            auto countMax = std::minmax_element(countsFieldsInts.begin(),countsFieldsInts.end());
            if (*countMax.second > readDepthMax) {
                continue; // Skip to the next locus
            }
            

            // Now that the lines are parsed, pull the major and minor alleles from angsd's mafs output
            std::string majorAllele = mafFields[2];
            std::string minorAllele = mafFields[3];
            
            
            /* Next, we'll count find how many reads there are supporting the major and minor
             alleles for each individual by looking through the corresponding counts file. The
             function that does this, "readCountToUnsignedChar" is implemented above. It returns
             the readcount number as an unsigned char for binary printing. The max value of an
             unsigned char is 255, but we are skipping all loci with individual read depths 
             above readDepthMax (which will usually be much less than 255)*/
            for (int individual = 0; individual < numIndividuals; individual++) {
                unsigned char majorAlleleUnsignedChar = readCountToUnsignedChar(individual, majorAllele, countsFields);
                unsigned char minorAlleleUnsignedChar = readCountToUnsignedChar(individual, minorAllele, countsFields);
                outCountsFile << majorAlleleUnsignedChar << minorAlleleUnsignedChar;
            }
            counter++; // Keep track of how many loci have been processed
        }
    }
    catch(const boost::iostreams::gzip_error& e) {
        std::cout << e.what() << '\n';
    }
    outCountsFile.close();
    
    return 0;
}



















