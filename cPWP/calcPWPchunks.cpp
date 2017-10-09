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

    std::cout << "Number of threads: " << numThreads << std::endl;
    //std::streampos size;
    std::ifstream file (binaryFile, std::ios::in|std::ios::binary);


    std::ifstream sampleFile(sampleNamesFile);

    // Get the list of sample names from the sampleNames file
    std::string sampleLine;
    std::vector<std::string> sampleLines;
    while (std::getline(sampleFile, sampleLine)) {
        sampleLines.push_back(sampleLine);
    }

    if (file.is_open()) {
        //size = file.tellg(); // Just a variable that shows position of stream--at end since ios::ate, so it's the file size. PROBABLY WON'T WORK FOR FILES LARGER THAN ~ 2GB!
        //file.seekg (0, std::ios::beg); // Go back to the beginning of the file
        //std::cout << "The total size of the file is " << size << "bytes. This corresponds to " << size/(numIndividuals*2) << " loci" << std::endl;



        unsigned long long int maxLocus = numLoci;
        /*
        if (numLoci == 0) {
            maxLocus = (unsigned long long)(size/(numIndividuals*2));
        } else {
            maxLocus = numLoci;
        }
         */

        std::cout << "Calculating divergence based on " << maxLocus << " total loci." << std::endl;
        // How many bytes to read in at one time (this number of loci will be split amongs numThreads threads, so it should be divisible exactly by numThreads. So the number of loci read in at a time will actually be numLoci*numThreads
        unsigned long long int lociChunkByteSize = (unsigned long long)lociChunkSize * numIndividuals * 2 * numThreads;
        int numFullChunks = (maxLocus*numIndividuals*2)/lociChunkByteSize; // Truncates answer to an integer
        unsigned long long remainingBytesAfterFullChunks = (maxLocus*numIndividuals*2) % lociChunkByteSize;

        std::cout << "Total number of chunks to run: " << numFullChunks + 1 << std::endl;


        /* We are going to split the loci between numThreads threads. Each thread will modify two multidimensional
         vectors of the forms std::vector< std::vector<long double> > pwp(numIndividuals, std::vector<long double>(numIndividuals,0))    and   std::vector< std::vector<unsigned long long int> > weightings(numIndividuals, std::vector<unsigned long long int>(numIndividuals,0))
         First, we'll generate all of these vectors, which apparently in C++ needs to be constructed of a
         vector of two-dimensional vectors...
         */
        std::vector<std::vector<std::vector<long double>>> pwpThreads(numThreads, std::vector<std::vector<long double>> (numIndividuals, std::vector<long double> (numIndividuals,0) ) ); //pwpThreads[0] is the first 2D array for the first thread, etc...
        std::vector<std::vector<std::vector<unsigned long long int>>> weightingsThreads(numThreads, std::vector<std::vector<unsigned long long int> > (numIndividuals, std::vector<unsigned long long int> (numIndividuals,0) ) );
        std::cout << "Initialized the 3d weighting and pwp vectors" << std::endl;


        //file.read((char*) &readCounts[0], size);

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

                //threadsVec.push_back(std::thread(calcPWPforRange, numIndividuals, lociPerThread, std::ref(readCounts), std::ref(pwpThreads[threadRunning]), std::ref(weightingsThreads[threadRunning])));
                threadsVec.push_back(std::thread(calcPWPforRange, firstLocus, finishingLocus, numIndividuals, std::ref(readCounts), std::ref(pwpThreads[threadRunning]), std::ref(weightingsThreads[threadRunning])));

            }

            // Wait on threads to finish
            for (int i = 0; i < numThreads; ++i) {
                threadsVec[i].join();
                std::cout << "Joined thread " << i << std::endl;
            }
            std::cout << "All threads completed running for chunk " << chunkCounter << " of " << numFullChunks + 1 << std::endl;
            chunkCounter++;
            std::cout << "Finished processing " << chunkCounter * lociPerThread * numThreads << " loci out of " << maxLocus << std::endl;
        }



        // For the last chunk, we'll just run it in a single thread, so we don't have to worry about numRemainingLoci/lociPerThread remainders... Just add everything to the vectors for thread 1
        std::vector<unsigned char> readCountsRemaining(remainingBytesAfterFullChunks);
        file.read((char*) &readCountsRemaining[0], remainingBytesAfterFullChunks);
        unsigned long long int finishingLocus = (readCountsRemaining.size()/(numIndividuals*2)) - 1;
        calcPWPforRange(0, finishingLocus, numIndividuals, std::ref(readCountsRemaining), std::ref(pwpThreads[0]), std::ref(weightingsThreads[0]));



        // Now aggregate the results of the threads and print final results
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


        // Now print out the final output to the pairwise pi file:
        std::ofstream pwpOUT (outFile);
        pwpOUT << "Sample1\tSample2\tPWP" << std::endl;
        int rowCounter = 0;
        if (!pwpOUT) {
            std::cerr << "Crap, " << outFile << "didn't open!" << std::endl;
        } else {
            for (int tortoise=0; tortoise < numIndividuals; tortoise++) {
                for (int comparisonTortoise = 0; comparisonTortoise <= tortoise; comparisonTortoise++) {
                    rowCounter++;

                    //std::cout << "Made it past the beginning of the last end for loop" << std::endl;
                    //std::cout << "Tortoise numbers: " << tortoise << " and " << comparisonTortoise << std::endl;
                    if (weightingsSum[tortoise][comparisonTortoise] > 0) {
                        //std::cout << weightings[tortoise][comparisonTortoise] << std::endl;
                        //std::cout << pwp[tortoise][comparisonTortoise] / weightings[tortoise][comparisonTortoise] << std::endl;
                        std::cout << std::fixed;
                        std::cout << "Weightings for tortoise " << tortoise << " and comparisonTortoise " << comparisonTortoise << " : " << weightingsSum[tortoise][comparisonTortoise] << std::endl;
                        std::cout << "PWP for tortoise " << tortoise << " and comparisonTortoise " << comparisonTortoise << " : " << pwpSum[tortoise][comparisonTortoise] << std::endl;
                        std::cout << std::scientific;
                        //pwpOUT << pwpSum[tortoise][comparisonTortoise] / weightingsSum[tortoise][comparisonTortoise] << std::endl;
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

    //std::cout << "Calculating PWP for the following locus range: " << startingLocus << " to " << endingLocus << std::endl;
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
                    //unsigned long long locusWeighting = (unsigned long long) (coverages[tortoise]*(coverages[tortoise])); // Shift weightings to match the weightings of the inter-comparisons by removing the -1
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
    std::cout << "Finished thread ending on locus " << endingLocus << std::endl;
    return 0;
}





/*
int calcPWPforRange (int numberIndividuals, unsigned long long int lociToCalc, std::vector<unsigned char>& mainReadCountVector, std::vector<std::vector<long double>>& threadPWP, std::vector<std::vector<unsigned long long int>>& threadWeightings) {

    std::cout << "Processing " << lociToCalc << " total loci" << std::endl;

    for( unsigned long long locus = 0; locus < lociToCalc; locus++) {
        std::cout << "Working on locus " << locus << std::endl;
        //std::cout << "Processing locus # " << locus << std::endl;
        if (locus % 100000 == 0) {
            std::cout << locus << " loci processed through calcPWPfromBinaryFile" << std::endl;
        }

        std::cout << "Made it to line 169" << std::endl;

        unsigned long long coverages[numberIndividuals];
        long double *majorAlleleFreqs = new long double[numberIndividuals]; // This will hold the major allele frequencies for that locus for each tortoise
        std::cout << "Made it to line 173" << std::endl;

        for( int tortoise = 0; tortoise < numberIndividuals; tortoise++ ) {
            unsigned long long majorIndex = locus * (numberIndividuals*2) + 2 * tortoise;
            std::cout << "Major index for tortoise " << tortoise << " and locus " << locus << " : " << majorIndex << std::endl;
            unsigned long long minorIndex = locus * (numberIndividuals*2) + 2 * tortoise + 1;
            std::cout << "Minor index for tortoise " << tortoise << " and locus " << locus << " : " << minorIndex << std::endl;

            std::cout << "Made it to line 178 for tortoise " << tortoise << " and locus " << locus << std::endl;

            coverages[tortoise] = int(mainReadCountVector[majorIndex]) + int(mainReadCountVector[minorIndex]); // Hold the coverages for each locus
            std::cout << "Made it to line 182 for tortoise" << tortoise << " and locus " << locus << std::endl;
            if ( coverages[tortoise] > 0 ) {
                majorAlleleFreqs[tortoise] = (long double)mainReadCountVector[majorIndex] / (long double)coverages[tortoise]; // Not necessarily an int, but could be 0 or 1
                std::cout << "Made it to line 185 for tortoise" << tortoise << " and locus " << locus << std::endl;
                if (coverages[tortoise] > 1) {
                    unsigned long long locusWeighting = (unsigned long long) (coverages[tortoise]*(coverages[tortoise]-1));
                    //unsigned long long locusWeighting = (unsigned long long) (coverages[tortoise]*(coverages[tortoise])); // Shift weightings to match the weightings of the inter-comparisons by removing the -1
                    threadWeightings[tortoise][tortoise] += (unsigned long long)locusWeighting; // This is an integer--discrete number of reads
                    std::cout << "Made it to line 190 for tortoise" << tortoise << " and locus " << locus << std::endl;
                    threadPWP[tortoise][tortoise] += (long double)(locusWeighting) * ((long double)2.0 * majorAlleleFreqs[tortoise] * ((long double)(coverages[tortoise]) - (long double)(mainReadCountVector[majorIndex])) / (long double)((coverages[tortoise])-(long double)1.0));
                }
                std::cout << "Made it to line 193 for tortoise" << tortoise << " and locus " << locus << std::endl;
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
    return 0;
}






*/
