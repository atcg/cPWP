//
//  calcCOVARchunks.cpp
//  
//
//  Created by Evan McCartney-Melstad on 2/13/15.
//
//
#include "calcCOVARchunks.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <thread>
#include <string>

int calcCOVARfromBinaryFile (std::string binaryFile, long long int numLoci, const int numIndividuals, std::string outFile, int lociChunkSize, std::string sampleNamesFile, const int numThreads) {
    // Get the list of sample names from the sampleNames file. We'll use these to make the output file readable:
    std::ifstream sampleFile(sampleNamesFile);
    std::string sampleLine;
    std::vector<std::string> sampleLines;
    while (std::getline(sampleFile, sampleLine)) {
        sampleLines.push_back(sampleLine);
    }

    std::cout << "Number of threads: " << numThreads << std::endl;

    std::ifstream file (binaryFile, std::ios::in|std::ios::binary);
    if (file.is_open()) {
        std::cout << "Calculating covariance based on " << numLoci << " total loci." << std::endl;
        
        // Determine how many bytes to read in at one time. The number of loci read in at a time will actually be numLoci*numThreads
        unsigned long long int lociChunkByteSize = (unsigned long long int)lociChunkSize * numIndividuals * 2 * numThreads;
        int numFullChunks = (numLoci*numIndividuals*2)/lociChunkByteSize; // Truncates answer to an integer
        unsigned long long int remainingBytesAfterFullChunks = (numLoci*numIndividuals*2) % lociChunkByteSize;
        
        if (remainingBytesAfterFullChunks != 0) {
            std::cout << "Total number of chunks to run: " << numFullChunks + 1 << std::endl;
        } else {
            std::cout << "Total number of chunks to run: " << numFullChunks << std::endl;
        }
        
        /* We are going to split the loci between numThreads threads. Each thread will modify two multidimensional
         vectors that will hold the weightings, the weighted sum of products, and weighted sums of the firsts.
         First, we'll generate all of these vectors of two-dimensional vectors:
         */
        std::vector<std::vector<std::vector<long long int>>> weightSumProductsThreads(numThreads, std::vector<std::vector<long long int>> (numIndividuals, std::vector<long long int> (numIndividuals,0) ) ); //covarThreads[0] is the first 2D array for the first thread, etc...
        std::vector<std::vector<std::vector<long long int>>> weightSumFirstThreads(numThreads, std::vector<std::vector<long long int>> (numIndividuals, std::vector<long long int> (numIndividuals,0) ) ); //covarThreads[0] is the first 2D array for the first thread, etc...
        std::vector<std::vector<std::vector<long long int>>> weightingsThreads(numThreads, std::vector<std::vector<long long int> > (numIndividuals, std::vector<long long int> (numIndividuals,0) ) );
        std::cout << "Initialized the 3d weighting and covar vectors" << std::endl;
        
        
        // Read in the data in chunks, and process each chunk using numThreads threads
        int chunkCounter = 0;
        while (chunkCounter < numFullChunks) {
            unsigned long long int bytesPerThread = lociChunkByteSize / numThreads;
            unsigned long long int lociPerThread = bytesPerThread / (numIndividuals*2);
            
            std::cout << "Running chunk #" << chunkCounter << std::endl;
            // Load the read counts into the vector readCounts, starting at byte 0
            std::vector<unsigned char> readCounts(lociChunkByteSize);
            file.read((char*) &readCounts[0], lociChunkByteSize);
                        
            std::vector<std::thread> threadsVec;
            for (int threadRunning = 0; threadRunning < numThreads; threadRunning++) {
                unsigned long long int firstLocus = (unsigned long long int) threadRunning * lociPerThread;
                unsigned long long int finishingLocus = ((unsigned long long int) threadRunning * lociPerThread) + lociPerThread - (unsigned long long int)1.0;
                
                threadsVec.push_back(std::thread(calcCOVARforRange, firstLocus, finishingLocus, numIndividuals, std::ref(readCounts), std::ref(weightSumProductsThreads[threadRunning]), std::ref(weightSumFirstThreads[threadRunning]), std::ref(weightingsThreads[threadRunning])));
            }

            // Wait on threads to finish
            for (int i = 0; i < numThreads; ++i) {
                threadsVec[i].join();
            }
            if (remainingBytesAfterFullChunks != 0) {
                std::cout << "All threads completed running for chunk " << chunkCounter + 1 << " of " << numFullChunks + 1 << std::endl;
            } else {
                std::cout << "All threads completed running for chunk " << chunkCounter + 1 << " of " << numFullChunks << std::endl;
            }
            chunkCounter++;
            std::cout << "Finished processing " << chunkCounter * lociPerThread * numThreads << " loci out of " << numLoci << std::endl;
        }

        //////// LAST CHUNK ////////
        // For the last chunk, we'll just run it in a single thread. Just add everything to the vectors for the first thread.
        // We can skip this if there aren't any bytes remaining (i.e. if total loci are perfectly divisibile by chunk size * number of threads).
        if (remainingBytesAfterFullChunks != 0) {
            std::vector<unsigned char> readCountsRemaining(remainingBytesAfterFullChunks);
            file.read((char*) &readCountsRemaining[0], remainingBytesAfterFullChunks);
            unsigned long long int finishingLocus = (readCountsRemaining.size()/(numIndividuals*2)) - 1;
            calcCOVARforRange(0, finishingLocus, numIndividuals, std::ref(readCountsRemaining), std::ref(weightSumProductsThreads[0]), std::ref(weightSumFirstThreads[0]), std::ref(weightingsThreads[0]));
        }
        
        // Now aggregate the results of the threads and print final results to output file
        std::vector<std::vector<long double>> weightingsSUM(numIndividuals, std::vector<long double>(numIndividuals,0));
        std::vector<std::vector<long double>> weightSumProductsSUM(numIndividuals, std::vector<long double>(numIndividuals,0));
        std::vector<std::vector<long double>> weightSumFirstSUM(numIndividuals, std::vector<long double>(numIndividuals,0));
        
        for (int tortoise = 0; tortoise < numIndividuals; tortoise++) {
            for (int comparisonTortoise = 0; comparisonTortoise <= tortoise; comparisonTortoise++) {
                for (int threadVector = 0; threadVector < numThreads; threadVector++) {
                    weightingsSUM[tortoise][comparisonTortoise] += weightingsThreads[threadVector][tortoise][comparisonTortoise];
                    weightSumProductsSUM[tortoise][comparisonTortoise] += weightSumProductsThreads[threadVector][tortoise][comparisonTortoise];
                }
            }
        }
        
        //weightSumFirstSUM is a full matrix (not half), so run comparison tortoise from 0 to numIndividuals, not 0 to tortoise
        for (int tortoise = 0; tortoise < numIndividuals; tortoise++) {
            for (int comparisonTortoise = 0; comparisonTortoise < numIndividuals; comparisonTortoise++) {
                for (int threadVector = 0; threadVector < numThreads; threadVector++) {
                    weightSumFirstSUM[tortoise][comparisonTortoise] += weightSumFirstThreads[threadVector][tortoise][comparisonTortoise];
                }
            }
        }

        std::cout << "Finished summing the threads vectors" << std::endl;\
        file.close(); // This is the binary file that holds all the read count data
        
        // Now print out the final output to the pairwise covariance file:
        std::ofstream covarOUT (outFile);
        if (!covarOUT) {
            std::cerr << "Shoot! " << outFile << "didn't open!" << std::endl;
        } else {
            for (int tortoise=0; tortoise < numIndividuals; tortoise++) {
                for (int comparisonTortoise = 0; comparisonTortoise <= tortoise; comparisonTortoise++) {
                    std::cout << "weightSumProducts of " << tortoise << " and " << comparisonTortoise << ": " << weightSumProductsSUM[tortoise][comparisonTortoise] << std::endl;
                    std::cout << "weightings of " << tortoise << " and " << comparisonTortoise << ": " << weightingsSUM[tortoise][comparisonTortoise] << std::endl;
                    std::cout << "weightSumFirst of " << tortoise << " and " << comparisonTortoise << ": " << weightSumFirstSUM[tortoise][comparisonTortoise] << std::endl;
                    std::cout << "weightSumFirst of " << comparisonTortoise << " and " << tortoise << ": " << weightSumFirstSUM[comparisonTortoise][tortoise] << std::endl;
                    
                    //(C[N,M]/W[N,M]) - (Z[N,M]/W[N,M])*(Z[M,N]/W[N,M]);
                    long double covar = ((weightSumProductsSUM[tortoise][comparisonTortoise])/(weightingsSUM[tortoise][comparisonTortoise])) - ((weightSumFirstSUM[tortoise][comparisonTortoise])/(weightingsSUM[tortoise][comparisonTortoise]))*((weightSumFirstSUM[comparisonTortoise][tortoise])/(weightingsSUM[tortoise][comparisonTortoise]));
                    covarOUT << sampleLines[tortoise] << "\t" << sampleLines[comparisonTortoise] << "\t" << covar << std::endl;
                }
            }
        }
    } else std::cout << "Unable to open file";
    
    return 0;
}

/*
 D[N] =  total coverage
 P[N] = major allele counts
 W[N,N] = weights
 C[N,N] = weighted sum of products
 Z[N,N] = weighted sum of first one (note cancelation of D[N])
 W[N,M] = weights
 C[N,M] = weighted sum of products
 Z[N,M] = weighted sum of first one (note cancelation of D[M]) : weightSumFirst
 Z[M,N] = weighted sum of second one (IN TRANSPOSE)            : weightSumFirst
 */


int calcCOVARforRange (long long int startingLocus, long long int endingLocus, const int numIndividuals, std::vector<unsigned char>& mainReadCountVector, std::vector<std::vector<long long int>>& weightSumProducts, std::vector<std::vector<long long int>>& weightSumFirst, std::vector<std::vector<long long int>>& threadWeightings) {    
    for( long long int locus = startingLocus; locus < endingLocus; locus++) {
        long long int coverages[numIndividuals];
        long long int *majorAlleleCounts = new long long int[numIndividuals]; // This will hold the major allele counts for that locus for each tortoise
        
        for( int tortoise = 0; tortoise < numIndividuals; tortoise++ ) {
            long long int majorIndex = locus * (numIndividuals*2) + 2 * tortoise;
            long long int minorIndex = locus * (numIndividuals*2) + 2 * tortoise + 1;
            coverages[tortoise] = int(mainReadCountVector[majorIndex]) + int(mainReadCountVector[minorIndex]); // Hold the coverages for each locus
            
            if ( coverages[tortoise] > 0 ) {
                majorAlleleCounts[tortoise] = (long long int)mainReadCountVector[majorIndex]; // This will be an integer value
                
                if (coverages[tortoise] > 1) {

                    long long int locusWeighting = (long long int) (coverages[tortoise]*(coverages[tortoise]-1));
                    threadWeightings[tortoise][tortoise] += (long long int)locusWeighting; // This is an integer--discrete number of reads

                    weightSumProducts[tortoise][tortoise] += (majorAlleleCounts[tortoise]) * (majorAlleleCounts[tortoise]) * (coverages[tortoise]-1)/(coverages[tortoise]);
                    weightSumFirst[tortoise][tortoise] += (majorAlleleCounts[tortoise]) * (coverages[tortoise]-1);

                }

                for( int comparisonTortoise = 0; comparisonTortoise < tortoise; comparisonTortoise++) {
                    if (coverages[comparisonTortoise] > 0) {
                        long long int locusWeighting = (long long int)coverages[tortoise] * (long long int)(coverages[comparisonTortoise]);
                        threadWeightings[tortoise][comparisonTortoise] += locusWeighting;
                        
                        weightSumProducts[tortoise][comparisonTortoise] += (majorAlleleCounts[tortoise]) * (majorAlleleCounts[comparisonTortoise]);
                        weightSumFirst[tortoise][comparisonTortoise] += (majorAlleleCounts[comparisonTortoise]) * coverages[tortoise];
                        weightSumFirst[comparisonTortoise][tortoise] += (majorAlleleCounts[tortoise]) * coverages[comparisonTortoise];
                    }
                }
            }
        }
        delete[] majorAlleleCounts; // Needed to avoid memory leaks
    }
    return 0;
}








