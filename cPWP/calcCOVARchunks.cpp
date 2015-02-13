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


int calcCOVARfromBinaryFile (std::string binaryFile, unsigned long long int numLoci, const int numIndividuals, std::string outFile, int lociChunkSize, int numThreads) {
    
    std::cout << "Number of threads: " << numThreads << std::endl;
    //std::streampos size;
    std::ifstream file (binaryFile, std::ios::in|std::ios::binary);
    std::cout << "Made it to line 22" << std::endl;
    
    if (file.is_open()) {

        unsigned long long int maxLocus = numLoci;
        
        std::cout << "Calculating divergence based on " << maxLocus << " total loci." << std::endl;
        
        // How many bytes to read in at one time (this number of loci will be split amongs numThreads threads, so it should be divisible exactly by numThreads. So the number of loci read in at a time will actually be numLoci*numThreads
        unsigned long long int lociChunkByteSize = (unsigned long long)lociChunkSize * numIndividuals * 2 * numThreads;
        int numFullChunks = (maxLocus*numIndividuals*2)/lociChunkByteSize; // Truncates answer to an integer
        unsigned long long remainingBytesAfterFullChunks = (maxLocus*numIndividuals*2) % lociChunkByteSize;
        
        std::cout << "Total number of chunks to run: " << numFullChunks + 1 << std::endl;
        
        /* We are going to split the loci between numThreads threads. Each thread will modify two multidimensional
         vectors that will hold the weightings, the weighted sum of products, and weighted sums of the firsts.
         First, we'll generate all of these vectors, which apparently in C++ needs to be constructed of a
         vector of two-dimensional vectors...
         */
        std::vector<std::vector<std::vector<unsigned long long>>> weightSumProductsThreads(numThreads, std::vector<std::vector<unsigned long long>> (numIndividuals, std::vector<unsigned long long> (numIndividuals,0) ) ); //covarThreads[0] is the first 2D array for the first thread, etc...
        std::vector<std::vector<std::vector<unsigned long long>>> weightSumFirstThreads(numThreads, std::vector<std::vector<unsigned long long>> (numIndividuals, std::vector<unsigned long long> (numIndividuals,0) ) ); //covarThreads[0] is the first 2D array for the first thread, etc...
        std::vector<std::vector<std::vector<unsigned long long int>>> weightingsThreads(numThreads, std::vector<std::vector<unsigned long long int> > (numIndividuals, std::vector<unsigned long long int> (numIndividuals,0) ) );
        std::cout << "Initialized the 3d weighting and covar vectors" << std::endl;
        
        
        // Read in the data in chunks, and process each chunk using numThreads threads
        int chunkCounter = 0;
        while (chunkCounter < numFullChunks) {
            unsigned long long bytesPerThread = lociChunkByteSize / numThreads;
            unsigned long long int lociPerThread = bytesPerThread / (numIndividuals*2);
            
            std::cout << "Running chunk #" << chunkCounter << std::endl;
            std::vector<unsigned char> readCounts(lociChunkByteSize);
            file.read((char*) &readCounts[0], lociChunkByteSize);
            
            std::cout << "Number of bytes in the chunk vector: " << readCounts.size() << std::endl;
            
            std::vector<std::thread> threadsVec;
            for (int threadRunning = 0; threadRunning < numThreads; threadRunning++) {
                
                unsigned long long int firstLocus = (unsigned long long int) threadRunning * lociPerThread;
                unsigned long long int finishingLocus = ((unsigned long long int) threadRunning * lociPerThread) + lociPerThread - (unsigned long long)1.0;
                
                std::cout << "Got to the function call in main loop. Running thread # " << threadRunning << std::endl;
                
                //threadsVec.push_back(std::thread(calcCOVARforRange, numIndividuals, lociPerThread, std::ref(readCounts), std::ref(covarThreads[threadRunning]), std::ref(weightingsThreads[threadRunning])));
                threadsVec.push_back(std::thread(calcCOVARforRange, firstLocus, finishingLocus, numIndividuals, std::ref(readCounts), std::ref(weightSumProductsThreads[threadRunning]), std::ref(weightSumFirstThreads[threadRunning]), std::ref(weightingsThreads[threadRunning])));
            }
            // Wait on threads to finish
            for (int i = 0; i < numThreads; ++i) {
                threadsVec[i].join();
                std::cout << "Joined thread " << i << std::endl;
            }
            std::cout << "All threads completed running for chunk " << chunkCounter << " of " << numFullChunks + 1 << std::endl;
            chunkCounter++;
        }
        
        
        
        // For the last chunk, we'll just run it in a single thread, so we don't have to worry about numRemainingLoci/lociPerThread remainders... Just add everything to the vectors for thread 1
        std::vector<unsigned char> readCountsRemaining(remainingBytesAfterFullChunks);
        file.read((char*) &readCountsRemaining[0], remainingBytesAfterFullChunks);
        unsigned long long int finishingLocus = (readCountsRemaining.size()/(numIndividuals*2)) - 1;
        calcCOVARforRange(0, finishingLocus, numIndividuals, std::ref(readCountsRemaining), std::ref(weightSumProductsThreads[0]), std::ref(weightSumFirstThreads[0]), std::ref(weightingsThreads[0]));
        
        
        // Now aggregate the results of the threads and print final results
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
            std::cerr << "Crap, " << outFile << "didn't open!" << std::endl;
        } else {
            for (int tortoise=0; tortoise < numIndividuals; tortoise++) {
                for (int comparisonTortoise = 0; comparisonTortoise <= tortoise; comparisonTortoise++) {

                    covar = long double(((weightSumProductsSUM[tortoise][comparisonTortoise])/(weightingsSUM[tortoise][comparisonTortoise])) - ((weightSumFirstSUM[tortoise][comparisonTortoise])/(weightingsSUM[tortoise][comparisonTortoise]))*((weightSumFirstSUM[comparisonTortoise][tortoise])/weightingsSUM[tortoise][comparisonTortoise]));
                    covarOUT << covar << std::endl;
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


int calcCOVARforRange (unsigned long long startingLocus, unsigned long long endingLocus, int numIndividuals, std::vector<unsigned char>& mainReadCountVector, std::vector<std::vector<long double>>& weightSumProducts, std::vector<std::vector<unsigned long long int>>& weightSumFirst, std::vector<std::vector<unsigned long long int>>& threadWeightings) {
    
    std::cout << "Calculating COVAR for the following locus range: " << startingLocus << " to " << endingLocus << std::endl;
    
    for( unsigned long long locus = startingLocus; locus < endingLocus; locus++) {
        if (locus % 100000 == 0) {
            std::cout << locus << " loci processed through calcCOVARfromBinaryFile" << std::endl;
        }
        
        unsigned long long coverages[numIndividuals];
        unsigned long long int *majorAlleleCounts = new unsigned long long int[numIndividuals]; // This will hold the major allele counts for that locus for each tortoise
        
        for( int tortoise = 0; tortoise < numIndividuals; tortoise++ ) {
            unsigned long long majorIndex = locus * (numIndividuals*2) + 2 * tortoise;
            unsigned long long minorIndex = locus * (numIndividuals*2) + 2 * tortoise + 1;
            coverages[tortoise] = int(mainReadCountVector[majorIndex]) + int(mainReadCountVector[minorIndex]); // Hold the coverages for each locus
            
            if ( coverages[tortoise] > 0 ) {
                majorAlleleCounts[tortoise] = (unsigned long long int)mainReadCountVector[majorIndex]; // This will be an integer value
                
                if (coverages[tortoise] > 1) {

                    unsigned long long int locusWeighting = (unsigned long long int) (coverages[tortoise]*(coverages[tortoise]-1));
                    threadWeightings[tortoise][tortoise] += (unsigned long long int)locusWeighting; // This is an integer--discrete number of reads

                    weightSumProducts[tortoise][tortoise] += (majorAlleleCounts[tortoise]) * (majorAlleleCounts[tortoise]) * (coverages[tortoise]-1)/(coverages[tortoise]);
                    weightSumFirst[tortoise][tortoise] += (majorAlleleCounts[tortoise]) * (coverages[tortoise]-1);

                }

                for( int comparisonTortoise = 0; comparisonTortoise < tortoise; comparisonTortoise++) {
                    if (coverages[comparisonTortoise] > 0) {
                        unsigned long long int locusWeighting = (unsigned long long int)coverages[tortoise] * (unsigned long long int)(coverages[comparisonTortoise]);
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
    std::cout << "Finished thread ending on locus " << endingLocus << std::endl;
    return 0;
}








