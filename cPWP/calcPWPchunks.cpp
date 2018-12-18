//
//  calcPWP.cpp
//
//
//  Created by Evan McCartney-Melstad on 1/10/15.
//
//
#include "calcPWPchunks.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <thread>
#include <string>

int calcPWPfromBinaryFile (std::string binaryFile, unsigned long long int numLoci, const int numIndividuals, std::string outFile, int lociChunkSize, std::string sampleNamesFile, int numThreads) {
    // Get the list of sample names from the sampleNames file. We'll store these in a vector so that we can make the output readable at the end
    std::ifstream sampleFile(sampleNamesFile);
    std::string sampleLine;
    std::vector<std::string> sampleLines;
    while (std::getline(sampleFile, sampleLine)) {
        sampleLines.push_back(sampleLine);
    }

    std::cout << "Number of threads: " << numThreads << std::endl;

    std::ifstream file (binaryFile, std::ios::in|std::ios::binary);
    if (file.is_open()) {
        std::cout << "Calculating divergence based on " << numLoci << " total loci." << std::endl;

        // How many bytes to read in at one time. The number of loci read in at a time will actually be numLoci*numThreads
        unsigned long long int lociChunkByteSize = (unsigned long long int)lociChunkSize * numIndividuals * 2 * numThreads;
        int numFullChunks = (numLoci*numIndividuals*2)/lociChunkByteSize; // Truncates answer to an integer
        unsigned long long int remainingBytesAfterFullChunks = (numLoci*numIndividuals*2) % lociChunkByteSize;

        if (remainingBytesAfterFullChunks != 0) {
            std::cout << "Total number of chunks to run: " << numFullChunks + 1 << std::endl;
        } else {
            std::cout << "Total number of chunks to run: " << numFullChunks << std::endl;
        }

        /* We are going to split the loci in the chunk between numThreads threads. Each thread will modify two multidimensional
         vectors of the forms std::vector< std::vector<long double> > pwp(numIndividuals, std::vector<long double>(numIndividuals,0))    
         and std::vector< std::vector<unsigned long long int> > weightings(numIndividuals, std::vector<unsigned long long int>(numIndividuals,0))
         First, we'll generate all of these vectors of two-dimensional vectors:
         */
        std::vector<std::vector<std::vector<long double>>> pwpThreads(numThreads, std::vector<std::vector<long double>> (numIndividuals, std::vector<long double> (numIndividuals,0) ) ); //pwpThreads[0] is the first 2D array for the first thread, etc...
        std::vector<std::vector<std::vector<unsigned long long int>>> weightingsThreads(numThreads, std::vector<std::vector<unsigned long long int> > (numIndividuals, std::vector<unsigned long long int> (numIndividuals,0) ) );
        std::cout << "Initialized the 3d weighting and pwp vectors" << std::endl;

        // Read in the data in chunks, and process each chunk using numThreads threads
        int chunkCounter = 0;
        while (chunkCounter < numFullChunks) {
            unsigned long long bytesPerThread = lociChunkByteSize / numThreads;
            unsigned long long int lociPerThread = bytesPerThread / (numIndividuals*2);

            std::cout << "Running chunk #" << chunkCounter << std::endl;
            // Load the read counts into the vector readCounts, starting at byte 0
            std::vector<unsigned char> readCounts(lociChunkByteSize);
            file.read((char*) &readCounts[0], lociChunkByteSize);
            
            std::vector<std::thread> threadsVec;
            for (int threadRunning = 0; threadRunning < numThreads; threadRunning++) {
                unsigned long long int firstLocus = (unsigned long long int) threadRunning * lociPerThread;
                unsigned long long int finishingLocus = ((unsigned long long int) threadRunning * lociPerThread) + lociPerThread - (unsigned long long int)1.0;

                threadsVec.push_back(std::thread(calcPWPforRange, firstLocus, finishingLocus, numIndividuals, std::ref(readCounts), std::ref(pwpThreads[threadRunning]), std::ref(weightingsThreads[threadRunning])));
            }

            // Wait on all threads to finish
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
            calcPWPforRange(0, finishingLocus, numIndividuals, std::ref(readCountsRemaining), std::ref(pwpThreads[0]), std::ref(weightingsThreads[0]));
        }

        //////// Aggregate the results of the threads and print final results into the output file ////////
        std::vector<std::vector<long double>> weightingsSum(numIndividuals, std::vector<long double>(numIndividuals,0));
        std::vector<std::vector<long double>> pwpSum(numIndividuals, std::vector<long double>(numIndividuals,0));

        for (int tortoise = 0; tortoise < numIndividuals; tortoise++) {
            for (int comparisonTortoise = 0; comparisonTortoise <= tortoise; comparisonTortoise++) {
                for (int threadVector = 0; threadVector < numThreads; threadVector++) {
                    weightingsSum[tortoise][comparisonTortoise] += weightingsThreads[threadVector][tortoise][comparisonTortoise];
                    pwpSum[tortoise][comparisonTortoise] += pwpThreads[threadVector][tortoise][comparisonTortoise];
                }
            }
        }
        std::cout << "Finished summing the threads vectors" << std::endl;\
        file.close(); // This is the binary file that holds all the read count data

        //////// Print out the final output to the pairwise pi file ////////
        std::ofstream pwpOUT (outFile);
        pwpOUT << "Sample1\tSample2\tPWP" << std::endl;
        int rowCounter = 0;
        if (!pwpOUT) {
            std::cerr << "Crap, " << outFile << "didn't open!" << std::endl;
        } else {
            for (int tortoise=0; tortoise < numIndividuals; tortoise++) {
                for (int comparisonTortoise = 0; comparisonTortoise <= tortoise; comparisonTortoise++) {
                    rowCounter++;

                    //std::cout << "Tortoise numbers: " << tortoise << " and " << comparisonTortoise << std::endl;
                    if (weightingsSum[tortoise][comparisonTortoise] > 0) {
                        /*
                        std::cout << weightings[tortoise][comparisonTortoise] << std::endl;
                        std::cout << pwp[tortoise][comparisonTortoise] / weightings[tortoise][comparisonTortoise] << std::endl;
                        std::cout << std::fixed;
                        std::cout << "Weightings for tortoise " << tortoise << " and comparisonTortoise " << comparisonTortoise << " : " << weightingsSum[tortoise][comparisonTortoise] << std::endl;
                        std::cout << "PWP for tortoise " << tortoise << " and comparisonTortoise " << comparisonTortoise << " : " << pwpSum[tortoise][comparisonTortoise] << std::endl;
                        std::cout << std::scientific;
                        */
                        pwpOUT << sampleLines[tortoise] << "\t" << sampleLines[comparisonTortoise] << "\t" << pwpSum[tortoise][comparisonTortoise] / weightingsSum[tortoise][comparisonTortoise] << std::endl;
                    } else {
                        pwpOUT << "NA" << std::endl;
                    }
                }
            }
        }
    } else {
        std::cout << "Unable to open file";
    }
    return 0;
}

int calcPWPforRange (unsigned long long startingLocus, unsigned long long endingLocus, int numIndividuals, std::vector<unsigned char>& mainReadCountVector, std::vector<std::vector<long double>>& threadPWP, std::vector<std::vector<unsigned long long int>>& threadWeightings) {
    //std::cout << "Calculating PWP for the following locus range: " << startingLocus + ((chunkCounter + 1) * numLoci * numThreads) << " to " << endingLocus + ((chunkCounter + 1) * numLoci * numThreads) << std::endl;
    for( unsigned long long locus = startingLocus; locus < endingLocus; locus++) {
        //std::cout << "Processing locus # " << locus << std::endl;
        //if (locus % 100000 == 0) {
        //    std::cout << locus << " loci processed through calcPWPfromBinaryFile" << std::endl;
        //}

        unsigned long long coverages[numIndividuals];
        long double *majorAlleleFreqs = new long double[numIndividuals]; // This will hold the major allele frequencies for that locus for each tortoise

        for( int tortoise = 0; tortoise < numIndividuals; tortoise++ ) {
            unsigned long long majorIndex = locus * (numIndividuals*2) + 2 * tortoise;
            unsigned long long minorIndex = locus * (numIndividuals*2) + 2 * tortoise + 1;

            coverages[tortoise] = int(mainReadCountVector[majorIndex]) + int(mainReadCountVector[minorIndex]); // Hold the coverages for each locus
            if ( coverages[tortoise] > 0 ) {
                majorAlleleFreqs[tortoise] = (long double)mainReadCountVector[majorIndex] / (long double)coverages[tortoise]; // Not necessarily an int, but could be 0 or 1

                if (coverages[tortoise] > 1) {
                    unsigned long long locusWeighting = (unsigned long long) (coverages[tortoise]*(coverages[tortoise]-1));
                    threadWeightings[tortoise][tortoise] += (unsigned long long)locusWeighting; // This is an integer--discrete number of reads

                    threadPWP[tortoise][tortoise] += (long double)(locusWeighting) * ((long double)2.0 * majorAlleleFreqs[tortoise] * ((long double)(coverages[tortoise]) - (long double)(mainReadCountVector[majorIndex])) / (long double)((coverages[tortoise])-(long double)1.0));
                }

                for( int comparisonTortoise = 0; comparisonTortoise < tortoise; comparisonTortoise++) {
                    if (coverages[comparisonTortoise] > 0) {
                        unsigned long long locusWeighting = (unsigned long long)coverages[tortoise] * (unsigned long long)(coverages[comparisonTortoise]);
                        threadWeightings[tortoise][comparisonTortoise] += locusWeighting;

                        threadPWP[tortoise][comparisonTortoise] += (long double)locusWeighting * (majorAlleleFreqs[tortoise] * ((long double)1.0-majorAlleleFreqs[comparisonTortoise]) + majorAlleleFreqs[comparisonTortoise] * ((long double)1.0-majorAlleleFreqs[tortoise]));
                    }
                }
            }
        }
        delete[] majorAlleleFreqs; // Needed to avoid memory leaks
    }
    //std::cout << "Finished thread ending on locus " << endingLocus + ((chunkCounter + 1) * numLoci * numThreads) << std::endl;
    return 0;
}






