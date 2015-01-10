//
//  calcPWP.cpp
//  
//
//  Created by Evan McCartney-Melstad on 1/10/15.
//
//

#include "calcPWP.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <thread>
#include <string>


int calcPWPfromBinaryFile (std::string binaryFile, unsigned long long int numLoci, const int numIndividuals, std::string outFile, int numThreads) {
    
    //****MODIFY THIS TO ONLY READ IN N LOCI AT A TIME, INSTEAD OF USING THE ENTIRE FILE****
    
    std::cout << "Number of threads: " << numThreads << std::endl;
    std::streampos size;
    std::ifstream file (binaryFile, std::ios::in|std::ios::binary|std::ios::ate);
    //ifstream file ("test500k.binary8bitunsigned", ios::in|ios::binary|ios::ate);
    if (file.is_open()) {
        size = file.tellg(); // Just a variable that shows position of stream--at end since ios::ate, so it's the file size. PROBABLY WON'T WORK FOR FILES LARGER THAN ~ 2GB!
        file.seekg (0, std::ios::beg); // Go back to the beginning of the file
        //file.read((char*)readCounts, size); // cast to a char* to give to file.read
        
        //unsigned char* readCounts;
        //readCounts = new unsigned char[size];
        std::vector<unsigned char> readCounts(size);
        file.read((char*) &readCounts[0], size);
        file.close();
        
        std::cout << "the entire file content is in memory" << std::endl;
        std::cout << "the total size of the file is " << size << std::endl;
        std::cout << "the number of elements in the readCounts vector is: " << readCounts.size() << std::endl; // Will give the total size bytes divided by the size of one element--so it gives the number of elements
        
        // We now have an array of numIndividuals * 2 (major and minor allele) * 1million (loci)
        //int totalLoci = (int)size / (numIndividuals*2); // The 1 million locus file has 999,999 sites in it (because of header line)
        //int totalLoci = size/(272*2);
        //std::vector< std::vector<long double> > pwp(numIndividuals, std::vector<long double>(numIndividuals,0));
        //std::vector< std::vector<unsigned long long int> > weightings(numIndividuals, std::vector<unsigned long long int>(numIndividuals,0));
        //long double pwp[numIndividuals][numIndividuals] = {0.0}; // This is the matrix that will hold the pwp estimates
        //unsigned long long int weightings[numIndividuals][numIndividuals] = {0.0}; // This is the matrix that will hold the weightings--need to use a long long because the values are basically equal to the coverage squared by the end
        
        
        
        /* We are going to split the loci between numThreads threads. Each thread will modify two multidimensional
         vectors of the forms std::vector< std::vector<long double> > pwp(numIndividuals, std::vector<long double>(numIndividuals,0))    and   std::vector< std::vector<unsigned long long int> > weightings(numIndividuals, std::vector<unsigned long long int>(numIndividuals,0))
         
         First, we'll generate all of these vectors, which apparently in C++ needs to be constructed of a
         vector of two-dimensional vectors...
         */
        std::vector<std::vector<std::vector<long double>>> pwpThreads(numThreads, std::vector<std::vector<long double>> (numIndividuals, std::vector<long double> (numIndividuals,0) ) ); //pwpThreads[0] is the first 2D array for the first thread, etc...
        std::vector<std::vector<std::vector<unsigned long long int>>> weightingsThreads(numThreads, std::vector<std::vector<unsigned long long int> > (numIndividuals, std::vector<unsigned long long int> (numIndividuals,0) ) );
        
        std::cout << "Initialized the 3d vectors" << std::endl;
        
        // Now we need to determine how many loci for each thread. If we want to use the entire binary file, instead of numLoci loci, then change this to lociPerThread = (size/(numIndividuals*2))/numThreads
        //unsigned long long int lociPerThread = numLoci / numThreads;
        
        //unsigned long long int lociPerThread = (readCounts.size()-1)/numThreads; // loci starts with 0, so need to subtract 1 from numLoci
        unsigned long long int lociPerThread = (unsigned long long)(numLoci-1)/(unsigned long long)numThreads;

        std::cout << "Initialized lociPerThread with " << numLoci << std::endl;
        
        std::vector<std::thread> threadsVec;
        for (int threadRunning = 0; threadRunning < numThreads; threadRunning++) {
            std::cout << "Got to the function call. Running thread # " << threadRunning << std::endl;
            unsigned long long int firstLocus = (unsigned long long int) threadRunning * lociPerThread;
            unsigned long long int finishingLocus = ((unsigned long long int) threadRunning * lociPerThread) + lociPerThread - (unsigned long long)1.0;
            std::cout << "Set firstLocus to " << firstLocus << " and finishingLocus to " << finishingLocus << std::endl;
            
            threadsVec.push_back(std::thread(calcPWPforRange, firstLocus, finishingLocus, numIndividuals, std::ref(readCounts), std::ref(pwpThreads[threadRunning]), std::ref(weightingsThreads[threadRunning])));
        }
        
        std::cout << "All threads completed launching" << std::endl;
        
        // Wait on threads to finish
        for (int i = 0; i < numThreads; ++i) {
            threadsVec[i].join();
            std::cout << "Joined thread " << i << std::endl;
        }
        std::cout << "All threads completed running" << std::endl;
        
        
        
        // Now aggregate the results of the threads and print final results
        std::vector<std::vector<long double>> weightingsSum(272, std::vector<long double>(272,0));
        std::vector<std::vector<long double>> pwpSum(272, std::vector<long double>(272,0));
        std::cout << "Created the final results vectors" << std::endl;
        
        for (int element = 0; element < 272; element++) {
            for (int comparisonElement = 0; comparisonElement <= element; comparisonElement++) {
                for (int threadVector = 0; threadVector <= numThreads; threadVector++) {
                    weightingsSum[element][comparisonElement] += weightingsThreads[threadVector][element][comparisonElement];
                    pwpSum[element][comparisonElement] += pwpThreads[threadVector][element][comparisonElement];
                }
            }
        }
        std::cout << "Finished summing the threads vectors" << std::endl;
        
        
        // Now print out the final output to the pairwise pi file:
        std::ofstream pwpOUT (outFile);
        int rowCounter = 0;
        if (!pwpOUT) {
            std::cerr << "Crap, " << outFile << "didn't open!" << std::endl;
        } else {
            for (int tortoise=0; tortoise <= (numIndividuals-1); tortoise++) {
                for (int comparisonTortoise = 0; comparisonTortoise <= tortoise; comparisonTortoise++) {
                    rowCounter++;
                    
                    //std::cout << "Made it past the beginning of the last end for loop" << std::endl;
                    //std::cout << "Tortoise numbers: " << tortoise << " and " << comparisonTortoise << std::endl;
                    if (weightingsSum[tortoise][comparisonTortoise] > 0) {
                        //std::cout << weightings[tortoise][comparisonTortoise] << std::endl;
                        //std::cout << pwp[tortoise][comparisonTortoise] / weightings[tortoise][comparisonTortoise] << std::endl;
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


//int calcPWPforRange (unsigned long long startingLocus, unsigned long long endingLocus, int numIndividuals, const std::vector<BYTE>& mainReadCountVector, std::vector< std::vector<long double> > & threadPWP, std::vector< std::vector<long double> > & threadWeightings) {
int calcPWPforRange (unsigned long long startingLocus, unsigned long long endingLocus, int numIndividuals, std::vector<unsigned char>& mainReadCountVector, std::vector<std::vector<long double>>& threadPWP, std::vector<std::vector<unsigned long long int>>& threadWeightings) {
     
     std::cout << "Calculating PWP for the following locus range: " << startingLocus << " to " << endingLocus << std::endl;
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
             
             coverages[tortoise] = int(mainReadCountVector[majorIndex]) + int(mainReadCountVector[minorIndex]); // Hold the coverages for each locus
             if ( coverages[tortoise] > 0 ) {
                 //std::cout << "Made it to line 222 for locus " << locus << std::endl;
                 majorAlleleFreqs[tortoise] = (double)mainReadCountVector[majorIndex] / (double)coverages[tortoise]; // Not necessarily an int, but could be 0 or 1

                 if (coverages[tortoise] > 1) {
                     unsigned long long locusWeighting = coverages[tortoise]*(coverages[tortoise]-1);
                     threadWeightings[tortoise][tortoise] += (unsigned long long)locusWeighting; // This is an int--discrete number of reads
     
                     threadPWP[tortoise][tortoise] += double(locusWeighting) * (2.0 * majorAlleleFreqs[tortoise] * (double(coverages[tortoise]) - double(mainReadCountVector[majorIndex]))) / (double((coverages[tortoise])-1.0));
                 }
     
                 for( int comparisonTortoise = 0; comparisonTortoise <= tortoise; comparisonTortoise++) {
                     if (coverages[comparisonTortoise] > 0) {
                         double locusWeighting = (double)coverages[tortoise] * (double)coverages[comparisonTortoise];
                         
                         threadWeightings[tortoise][comparisonTortoise] += locusWeighting;
                         threadPWP[tortoise][comparisonTortoise] += (double)locusWeighting * (majorAlleleFreqs[tortoise] * (1.0-majorAlleleFreqs[comparisonTortoise]) + majorAlleleFreqs[comparisonTortoise] * (1.0-majorAlleleFreqs[tortoise]));
                    }
                }
            }
        }
        //std::cout << "Made it to deleting majorAlleleFreqs" << std::endl;
        delete[] majorAlleleFreqs; // Needed to avoid memory leaks
    }
    std::cout << "Finished thread ending on locus " << endingLocus << std::endl;
    return 0;
}
