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


int calcPWPfromBinaryFile (std::string binaryFile, unsigned long long int numLoci, const int numIndividuals, std::string outFile, int lociChunkSize, int numThreads) {
    
    std::cout << "Number of threads: " << numThreads << std::endl;
    std::streampos size;
    std::ifstream file (binaryFile, std::ios::in|std::ios::binary|std::ios::ate);

    if (file.is_open()) {
        size = file.tellg(); // Just a variable that shows position of stream--at end since ios::ate, so it's the file size. PROBABLY WON'T WORK FOR FILES LARGER THAN ~ 2GB!
        file.seekg (0, std::ios::beg); // Go back to the beginning of the file
        std::cout << "The total size of the file is " << size << "bytes. This corresponds to " << size/(numIndividuals*2) << " loci" << std::endl;
        
        unsigned long long int maxLocus;
        if (numLoci == 0) {
            maxLocus = (unsigned long long)(size/(numIndividuals*2));
        } else {
            maxLocus = numLoci;
        }
        
        std::cout << "Calculating divergence based on " << maxLocus << " total loci." << std::endl;
        
        // How many bytes to read in at one time (this number of loci will be split amongs numThreads threads, so it should be divisible exactly by numThreads. So the number of loci read in at a time will actually be numLoci*numThreads
        unsigned long long int lociChunkByteSize = (unsigned long long)lociChunkSize * numIndividuals * 2 * numThreads;
        int numFullChunks = (maxLocus*numIndividuals*2)/lociChunkByteSize; // Truncates answer to an integer
        unsigned long long remainingBytesAfterFullChunks = (maxLocus*numIndividuals*2) % lociChunkByteSize;
        
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
            std::vector<unsigned char> readCounts(lociChunkByteSize);
            file.read((char*) &readCounts[0], lociChunkByteSize);
            
            std::vector<std::thread> threadsVec;
            for (int threadRunning = 0; threadRunning < numThreads; threadRunning++) {
                unsigned long long bytesPerThread = lociChunkByteSize / numThreads;
                std::vector<unsigned char> readCounts(bytesPerThread);
                file.read((char*) &readCounts[0], bytesPerThread);
                unsigned long long int lociPerThread = bytesPerThread / (numIndividuals*2);
                
                std::cout << "Got to the function call in main loop. Running thread # " << threadRunning << std::endl;
                
                threadsVec.push_back(std::thread(calcPWPforRange, numIndividuals, lociPerThread, std::ref(readCounts), std::ref(pwpThreads[threadRunning]), std::ref(weightingsThreads[threadRunning])));
            }
            
            // Wait on threads to finish
            for (int i = 0; i < numThreads; ++i) {
                threadsVec[i].join();
                std::cout << "Joined thread " << i << std::endl;
            }
            std::cout << "All threads completed running for chunk " << chunkCounter << " of " << numFullChunks + 1 << std::endl;
            chunkCounter++;
        }
        
        // That takes care of all the full-sized loci chunks, now deal with the remainder of loci
        std::vector<unsigned char> readCountsRemaining(remainingBytesAfterFullChunks);
        file.read((char*) &readCountsRemaining[0], remainingBytesAfterFullChunks);
        unsigned long long remainingLociAfterFullChunks = (remainingBytesAfterFullChunks/(numIndividuals*2));

        std::vector<std::thread> threadsVecRemaining;
        for (int threadRunning = 0; threadRunning < numThreads; threadRunning++) {
            std::cout << "Got to the function call in the remaining loci. Running thread # " << threadRunning << std::endl;
          
            threadsVecRemaining.push_back(std::thread(calcPWPforRange, numIndividuals, remainingLociAfterFullChunks, std::ref(readCountsRemaining), std::ref(pwpThreads[threadRunning]), std::ref(weightingsThreads[threadRunning])));
        }
        
        // Wait on threads to finish
        for (int i = 0; i < numThreads; ++i) {
            threadsVecRemaining[i].join();
            std::cout << "Joined thread " << i << std::endl;
        }
        std::cout << "All threads completed running for last chunk." << std::endl;
        
        

        
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
                        pwpOUT << pwpSum[tortoise][comparisonTortoise] / weightingsSum[tortoise][comparisonTortoise] << std::endl;
                    } else {
                        pwpOUT << "NA" << std::endl;
                    }
                }
            }
        }
    } else std::cout << "Unable to open file";
    
    return 0;
}



int calcPWPforRange (int numberIndividuals, unsigned long long int lociToCalc, std::vector<unsigned char>& mainReadCountVector, std::vector<std::vector<long double>>& threadPWP, std::vector<std::vector<unsigned long long int>>& threadWeightings) {
    
    for( unsigned long long locus = 0; locus < lociToCalc; locus++) {
        //std::cout << "Processing locus # " << locus << std::endl;
        if (locus % 100000 == 0) {
            std::cout << locus << " loci processed through calcPWPfromBinaryFile" << std::endl;
        }
        
        unsigned long long coverages[numberIndividuals];
        long double *majorAlleleFreqs = new long double[numberIndividuals]; // This will hold the major allele frequencies for that locus for each tortoise
        
        for( int tortoise = 0; tortoise < numberIndividuals; tortoise++ ) {
            unsigned long long majorIndex = locus * (numberIndividuals*2) + 2 * tortoise;
            unsigned long long minorIndex = locus * (numberIndividuals*2) + 2 * tortoise + 1;
            
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
    return 0;
}
















