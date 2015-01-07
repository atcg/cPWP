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
#include <thread>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/algorithm/string.hpp>

typedef unsigned char BYTE;



int runANGSDforReadCounts (std::string bamlist, std::string angsdPrefix, std::string nThreads, std::string angsdOutputLog) {
    std::cout << "**********\nChecking if processor is available to run ANGSD...**********\n";
    if (system(NULL)) puts ("Ok");
    else exit (EXIT_FAILURE);
    
    //std::string angsdCommand = "/home/evan/bin/angsd0.613/angsd -bam " + bamlist + " -out " + angsdPrefix + " -uniqueOnly 1 -only_proper_pairs 1 -remove_bads 1 -doCounts 1 -nThreads " + nThreads + " -dumpCounts 4 -doMaf 1 -doMajorMinor 2 -GL 2 -minMapQ 20 -minQ 30 -SNP_pval 1e-6 > " + angsdOutputLog + " 2>&1";
    
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
    //std::string countsFileName = angsdPrefix + ".counts.gz";
    std::string countsFileName = angsdPrefix + ".counts";
    //std::string mafsFileName = angsdPrefix + ".mafs.gz";
    std::string mafsFileName = angsdPrefix + ".mafs";
    //std::ifstream countsFile(countsFileName, std::ios_base::in | std::ios_base::binary);
    std::ifstream countsFile(countsFileName, std::ios_base::in);
    //std::ifstream mafsFile(mafsFileName, std::ios_base::in | std::ios_base::binary);
    std::ifstream mafsFile(mafsFileName, std::ios_base::in);
    std::ofstream outCountsFile(binaryOutputFileName, std::ios_base::out | std::ios_base::binary);
    
    try {
        /*
        boost::iostreams::filtering_istream in;
        boost::iostreams::filtering_istream inMafs;
        in.push(boost::iostreams::gzip_decompressor());
        in.push(countsFile);
        inMafs.push(boost::iostreams::gzip_decompressor());
        inMafs.push(mafsFile);
        */
        int counter = 0;
        //int numIndividuals = 3;
        for(std::string str; std::getline(countsFile, str); )
        {
            /* Make sure to get the maf string before we do anything else to keep it in
             sync with the counts string */
            std::string mafStr;
            std::getline(mafsFile, mafStr);
            
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





int calcPWPfromBinaryFile (std::string binaryFile, unsigned long long int numLoci, const int numIndividuals, std::string outFile, int numThreads) {
    typedef unsigned char BYTE;
    
    //****MODIFY THIS TO ONLY READ IN N LOCI AT A TIME, INSTEAD OF USING THE ENTIRE FILE****
    
    
    std::streampos size;
    std::ifstream file (binaryFile, std::ios::in|std::ios::binary|std::ios::ate);
    //ifstream file ("test500k.binary8bitunsigned", ios::in|ios::binary|ios::ate);
    if (file.is_open()) {
        size = file.tellg(); // Just a variable that shows position of stream--at end since ios::ate, so it's the file size. PROBABLY WON'T WORK FOR FILES LARGER THAN ~ 2GB!
        file.seekg (0, std::ios::beg); // Go back to the beginning of the file
        //file.read((char*)readCounts, size); // cast to a char* to give to file.read

        //unsigned char* readCounts;
        //readCounts = new unsigned char[size];
        std::vector<BYTE> readCounts(size);
        file.read((char*) &readCounts[0], size);
        file.close();
        
        std::cout << "the entire file content is in memory" << std::endl;
        std::cout << "the total size of the file is " << size << std::endl;
        std::cout << "the number of elements in the readCounts vector is: " << readCounts.size() << std::endl; // Will give the total size bytes divided by the size of one element--so it gives the number of elements
        
        // We now have an array of numIndividuals * 2 (major and minor allele) * 1million (loci)
        //int totalLoci = (int)size / (numIndividuals*2); // The 1 million locus file has 999,999 sites in it (because of header line)
        //int totalLoci = size/(272*2);
        std::vector< std::vector<long double> > pwp(numIndividuals, std::vector<long double>(numIndividuals,0));
        std::vector< std::vector<unsigned long long int> > weightings(numIndividuals, std::vector<unsigned long long int>(numIndividuals,0));
        //long double pwp[numIndividuals][numIndividuals] = {0.0}; // This is the matrix that will hold the pwp estimates
        //unsigned long long int weightings[numIndividuals][numIndividuals] = {0.0}; // This is the matrix that will hold the weightings--need to use a long long because the values are basically equal to the coverage squared by the end
        
        
        
        /* We are going to split the loci between numThreads threads. Each thread will modify two multidimensional
         vectors of the forms std::vector< std::vector<long double> > pwp(numIndividuals, std::vector<long double>(numIndividuals,0))    and   std::vector< std::vector<unsigned long long int> > weightings(numIndividuals, std::vector<unsigned long long int>(numIndividuals,0))
         
         First, we'll generate all of these vectors, which apparently in C++ needs to be constructed of a 
         vector of two-dimensional vectors... 
         */
        //std::vector<std::vector<std::vector<unsigned long long int>>> pwpThreads(numThreads, std::vector<std::vector<unsigned long long int>> (numIndividuals, std::vector<unsigned long long int> (numIndividuals,0) ) ); //pwpThreads[0] is the first 2D array for the first thread, etc...
        //std::vector<std::vector<std::vector<unsigned long long int>>> weightingsThreads(numThreads, std::vector<std::vector<unsigned long long int> > (numIndividuals, std::vector<unsigned long long int> (numIndividuals,0) ) );

        // Now we need to determine how many loci for each thread. If we want to use the entire binary file, instead of numLoci loci, then change this to lociPerThread = (size/(numIndividuals*2))/numThreads
        unsigned long long int lociPerThread = numLoci / numThreads;
        
        
        std::thread t;
        //std::thread t[numThreads];
        for (int threadRunning; threadRunning < numThreads; threadRunning++) {
            unsigned long long firstLocus = (unsigned long long) threadRunning * lociPerThread;
            unsigned long long finishingLocus = ((unsigned long long) threadRunning * lociPerThread) + lociPerThread - (unsigned long long)1.0;
            
            // Since we're passing the vectors in by reference, calPWPforRange can modify pwpThreads[thread] and weightingsThreads[thread] (but not mainReadCountVector because that is declared in the function as const
            //t[threadRunning]
            //std::thread t = std::thread(calcPWPforRange, firstLocus, finishingLocus, std::ref(readCounts), std::ref(pwpThreads[threadRunning]), std::ref(weightingsThreads[threadRunning]));
            calcPWPforRange(firstLocus, finishingLocus, 272, readCounts, pwp, weightings);
        }
        
        //t.join();
        /*// Wait on threads to finish
        for (int i = 0; i < numThreads; ++i) {
            t[i].join();
        }
         */
        
        // Now aggregate the results of the threads and print final results

        
        
    }
    else std::cout << "Unable to open file";
    
    return 0;
}


//int calcPWPforRange (unsigned long long startingLocus, unsigned long long endingLocus, int numIndividuals, const std::vector<BYTE>& mainReadCountVector, std::vector< std::vector<long double> > & threadPWP, std::vector< std::vector<long double> > & threadWeightings) {
int calcPWPforRange (unsigned long long startingLocus, unsigned long long endingLocus, int numIndividuals, std::vector<unsigned char> &mainReadCountVector, std::vector<std::vector<long double>> &threadPWP, std::vector<std::vector<long double>>& threadWeightings) {

    
    //usage: calcPWPforRange(0, 1000000, readCounts) // where readCounts is a vector with all the read count data
    //this function will return both the pwp calculations for the matrix as well as the weightings, which will later
    //be summed across all threads
    
    for( unsigned long long locus = startingLocus; locus < endingLocus; locus++) {
    //std::cout << "Processing locus # " << locus << std::endl;
        if (locus % 100000 == 0) {
            std::cout << locus << " loci processed through calcPWPfromBinaryFile" << std::endl;
        }
        
        int coverages[numIndividuals];
        double *majorAlleleFreqs = new double[numIndividuals]; // This will hold the major allele frequencies for that locus for each tortoise
        
        for( int tortoise = 0; tortoise <= (numIndividuals-1); tortoise++ ) {
            unsigned long long majorIndex = locus * (numIndividuals*2) + 2 * tortoise;
            unsigned long long minorIndex = locus * (numIndividuals*2) + 2 * tortoise + 1;
            
            //std::cout << "\tTrying to access readCounts[" << minorIndex << "]" << ". locus: " << locus << ". numIndividuals: " << numIndividuals << ". tortoise: " << tortoise << std::endl;
            
            coverages[tortoise] = int(mainReadCountVector[majorIndex]) + int(mainReadCountVector[minorIndex]); // Hold the coverages for each locus
            //std::cout << "\t\tCalced coverage in line 232" << std::endl;
            //std::cout << coverages[tortoise] << std::endl;
            
            //std::cout << "Total coverage for tortoise " << tortoise << " at locus " << locus+1 << ": " << coverages[tortoise] << std::endl;
            
            if ( coverages[tortoise] > 0 ) {
                //std::cout << "Made it to line 222 for locus " << locus << std::endl;
                majorAlleleFreqs[tortoise] = (double)mainReadCountVector[majorIndex] / (double)coverages[tortoise]; // Not necessarily an int, but could be 0 or 1
                //std::cout << "Major allele frequency for individual " << tortoise << " at locus " << locus << ": " << majorAlleleFreqs[tortoise] << std::endl;
                //std::cout << "\t\t\tCalced majorAlleleFreqs[" << tortoise << "] in line 239" << std::endl;
                if (coverages[tortoise] > 1) {
                    unsigned long long locusWeighting = coverages[tortoise]*(coverages[tortoise]-1);
                    threadWeightings[tortoise][tortoise] += (unsigned long long)locusWeighting; // This is an int--discrete number of reads
                    //std::cout << "\t\t\t\tCalced weightings in line 245 for tortoise[" << tortoise << "]" << std::endl;
                    
                    threadPWP[tortoise][tortoise] += double(locusWeighting) * (2.0 * majorAlleleFreqs[tortoise] * (double(coverages[tortoise]) - double(mainReadCountVector[majorIndex]))) / (double((coverages[tortoise])-1.0));
                    //std::cout << "\t\t\t\t\tCalced pwp in line 247 for tortoise[" << tortoise << "]" << std::endl;
                    //std::cout << "Locus self weighting for individual " << tortoise << " at locus: " << locus << ": " << locusWeighting << ". Locus self PWP: " << double(locusWeighting) * (2.0 * majorAlleleFreqs[tortoise] * (double(coverages[tortoise]) - double(readCounts[locus * numIndividuals * 2 + 2 * tortoise]))) / (double((coverages[tortoise])-1.0)) << std::endl;
                    //std::cout << "\tmajorAlleleFreq: " << majorAlleleFreqs[tortoise] << ". Coverages: " << double(coverages[tortoise]) << ". readCounts: " << double(readCounts[locus * numIndividuals * 2 + 2 * tortoise]) << std::endl;
                    //std::cout << "PWP for self:" << pwp[tortoise][tortoise] << std::endl;
                    //std::cout << "Made it to line 233 for locus " << locus << ". Weightings = " << weightings[tortoise][tortoise] << ". PWP = " << pwp[tortoise][tortoise] << std::endl;
                }
                
                
                for( int comparisonTortoise = 0; comparisonTortoise < tortoise; comparisonTortoise++) {
                    if (coverages[comparisonTortoise] > 0) {
                        
                        //std::cout << "Coverages for the two comparison individuals: " << std::to_string((double)coverages[tortoise]) << " and " << std::to_string((double)coverages[comparisonTortoise]) << std::endl;
                        //std::cout << "Major allele freqs for the two comparison individuals: " << std::to_string(majorAlleleFreqs[tortoise]) << " and " << std::to_string(majorAlleleFreqs[comparisonTortoise]) << std::endl;
                        double locusWeighting = (double)coverages[tortoise] * (double)coverages[comparisonTortoise];
                        threadWeightings[tortoise][comparisonTortoise] += locusWeighting;
                        //std::cout << "locusDiffPWP: " << (double)locusWeighting * ((double)majorAlleleFreqs[tortoise] * (1-(double)majorAlleleFreqs[comparisonTortoise]) + (double)majorAlleleFreqs[comparisonTortoise] * (1-(double)majorAlleleFreqs[tortoise])) << std::endl;
                        threadPWP[tortoise][comparisonTortoise] += (double)locusWeighting * (majorAlleleFreqs[tortoise] * (1.0-majorAlleleFreqs[comparisonTortoise]) + majorAlleleFreqs[comparisonTortoise] * (1.0-majorAlleleFreqs[tortoise]));
                        //std::cout << "\t\t\t\t\t\tCalced pwp for tortoise[" << tortoise << "] and comparisonTortoise[" << comparisonTortoise << "]" << std::endl;
                        //std::cout << pwp[tortoise][comparisonTortoise] << std::endl;
                        //std::cout << "Cumulative weightings = " << weightings[tortoise][comparisonTortoise] << ". Cumulative PWP = " << pwp[tortoise][comparisonTortoise] << std::endl;
                    }
                }
            }
        }
        delete[] majorAlleleFreqs; // Needed to avoid memory leaks
        //delete[] coverages; // Since that locus is done, and this variable holds per-locus coverages, we can nuke the coverages to start anew
    }
    return 0;
}


