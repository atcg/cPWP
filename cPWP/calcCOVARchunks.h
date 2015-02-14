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


int calcCOVARfromBinaryFile (std::string binaryFile, long long int numLoci, const int numIndividuals, std::string outFile, int lociChunkSize, const int numThreads);


int calcCOVARforRange (long long int startingLocus, long long int endingLocus, int numIndividuals, std::vector<unsigned char>& mainReadCountVector, std::vector<std::vector<long long int>>& weightSumProducts, std::vector<std::vector<long long int>>& weightSumFirst, std::vector<std::vector<long long int>>& threadWeightings);
#endif /* defined(____calcCOVARchunks__) */
