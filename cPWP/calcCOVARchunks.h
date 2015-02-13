//
//  calcCOVARchunks.h
//  
//
//  Created by Evan McCartney-Melstad on 2/13/15.
//
//

#ifndef ____calcCOVARchunks__
#define ____calcCOVARchunks__

#include <iostream>
#include <vector>


int calcCOVARfromBinaryFile (std::string binaryFile, unsigned long long int numLoci, const int numIndividuals, std::string outFile, int lociChunkSize, int numThreads);

int calcCOVARforRange (unsigned long long startingLocus, unsigned long long endingLocus, int numIndividuals, std::vector<unsigned char>& mainReadCountVector, std::vector<std::vector<long double>>& weightSumProducts, std::vector<std::vector<unsigned long long int>>& weightSumFirst, std::vector<std::vector<unsigned long long int>>& threadWeightings);

#endif /* defined(____calcCOVARchunks__) */
