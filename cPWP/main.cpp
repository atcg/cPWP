//
//  main.cpp
//  cPWP
//
//  Created by Evan McCartney-Melstad on 12/31/14.
//  Copyright (c) 2014 Evan McCartney-Melstad. All rights reserved.
//

/*
 Program for calculating pairwise pi from a file full of unsigned char integers
*/


#include <iostream>
#include <fstream>
#include <vector>
#include "bamsToPWP.h"

using namespace std;



int main (int argc, char *argv[]) {
    streampos size;
    
    ifstream file ("counts.1.binary8bitunsigned", ios::in|ios::binary|ios::ate);
    //ifstream file ("test500k.binary8bitunsigned", ios::in|ios::binary|ios::ate);
    if (file.is_open()) {
        size = file.tellg(); // Just a variable that shows position of stream--at end since ios::ate, so it's the file size
        
        unsigned char* readCounts;
        readCounts = new unsigned char[size];
        
        file.seekg (0, ios::beg); // Go back to the beginning of the file
        file.read((char*)readCounts, size); // cast to a char* to give to file.read
        file.close();
        
        cout << "the entire file content is in memory" << endl;
        cout << "the total size of the file is " << size << endl;
        
        // We now have an array of 276 (number of torts) * 2 (major and minor allele) * 1million (loci) = 552 million elements, assuming 1 million loci and 276 tortoises
        // We will work with them in chunks of 552, as that is the number of integers within each locus
        int totalLoci = (int)size / (276*2); // The 1 million locus file has 999,999 sites in it (because of header line)
        //int totalLoci = 10000;
        
        long double pwp[276][276]; // This is the matrix that will hold the pwp estimates
        unsigned long long int weightings[276][276]; // This is the matrix that will hold the weightings--need to use a long long because the values are basically equal to the coverage squared by the end
        
        
        
        
        for( int locus = 0; locus < totalLoci; locus++) {
            int coverages[276];
            //int *coverages = new int[276]; // This will hold the coverages for the locus being evaluated
            double *majorAlleleFreqs = new double[276]; // This will hold the major allele frequencies for that locus for each tortoise
            
            for( int tortoise = 0; tortoise <= 275; tortoise++ ) {
                
                coverages[tortoise] = int(readCounts[locus * 552 + 2 * tortoise]) + int(readCounts[locus * 552 + 2 * tortoise + 1]); // Hold the coverages for each locus
                //cout << coverages[tortoise] << endl;
                
                //cout << "Total coverage for tortoise " << tortoise << " at locus " << locus+1 << ": " << coverages[tortoise] << endl;
                
                if ( coverages[tortoise] > 0 ) {
                    majorAlleleFreqs[tortoise] = (double)readCounts[locus * 552 + (2*tortoise)] / (double)coverages[tortoise]; // Not necessarily an int, but could be 0 or 1
                    //cout << "Major allele frequency: " << majorAlleleFreqs[tortoise] << endl;
                    
                    if (coverages[tortoise] > 1) {
                        unsigned long long locusWeighting = coverages[tortoise]*(coverages[tortoise]-1);
                        weightings[tortoise][tortoise] += (unsigned long long)locusWeighting; // This is an int--discrete number of reads
                        
                        
                        pwp[tortoise][tortoise] += double(locusWeighting) * (2.0 * majorAlleleFreqs[tortoise] * (double(coverages[tortoise]) - double(readCounts[locus * 552 + 2 * tortoise]))) / (double((coverages[tortoise])-1.0));
                        //cout << "PWP for self:" << pwp[tortoise][tortoise] << endl;
                    }
                    
                    
                    for( int comparisonTortoise = 0; comparisonTortoise < tortoise; comparisonTortoise++) {
                        if (coverages[comparisonTortoise] > 0) {
                            
                            
                            double locusWeighting = (double)coverages[tortoise] * (double)coverages[comparisonTortoise];
                            weightings[tortoise][comparisonTortoise] += locusWeighting;
                            //cout << "locusDiffPWP: " << (double)locusWeighting * ((double)majorAlleleFreqs[tortoise] * (1-(double)majorAlleleFreqs[comparisonTortoise]) + (double)majorAlleleFreqs[comparisonTortoise] * (1-(double)majorAlleleFreqs[tortoise])) << endl;
                            pwp[tortoise][comparisonTortoise] += (double)locusWeighting * (majorAlleleFreqs[tortoise] * (1.0-majorAlleleFreqs[comparisonTortoise]) + majorAlleleFreqs[comparisonTortoise] * (1.0-majorAlleleFreqs[tortoise]));
                            //cout << pwp[tortoise][comparisonTortoise] << endl;
                        }
                    }
                }
            }
            //delete[] coverages; // Since that locus is done, and this variable holds per-locus coverages, we can nuke the coverages to start anew
        }
        //cout << pwp[2][1] / weightings[2][1] << endl;
        
        // Now print out the final output to the pairwise pi file:
        ofstream pwpOUT ("pwptest.txt");
        int rowCounter = 0;
        //cout << "Made it past the ofstream call" << endl;
        if (!pwpOUT) {
            cerr << "Crap, pwptest.txt didn't open!" << endl;
        } else {
            // cout << "Made it past the check to see if !pwpOUT" << endl;
            for (int tortoise=0; tortoise <= 275; tortoise++) {
                // cout << "Made it past the beginning of the end for loop" << endl;
                for (int comparisonTortoise = 0; comparisonTortoise < tortoise; comparisonTortoise++) {
                    rowCounter++;
                    
                    if (pwp[tortoise][comparisonTortoise] / weightings[tortoise][comparisonTortoise] < 0.12 and pwp[tortoise][comparisonTortoise] / weightings[tortoise][comparisonTortoise] > 0) {
                        cout << "Possible duplicate pair: " << tortoise << " and " << comparisonTortoise << ". Look at line:" << rowCounter << endl;
                    }
                    //cout << "Made it past the beginning of the last end for loop" << endl;
                    //cout << "Tortoise numbers: " << tortoise << " and " << comparisonTortoise << endl;
                    if (weightings[tortoise][comparisonTortoise] > 0) {
                        //cout << weightings[tortoise][comparisonTortoise] << endl;
                        //cout << pwp[tortoise][comparisonTortoise] / weightings[tortoise][comparisonTortoise] << endl;
                        pwpOUT << pwp[tortoise][comparisonTortoise] / weightings[tortoise][comparisonTortoise] << endl;
                        if (pwp[tortoise][comparisonTortoise] / weightings[tortoise][comparisonTortoise] < 0) {
                            cout << "Negative value for tortoise " << tortoise << " and comparisonTortoise " << comparisonTortoise << " . Respective weightings: " << weightings[tortoise][comparisonTortoise] << " . Respective pwp: " << pwp[tortoise][comparisonTortoise] << endl;
                            
                        }
                        if (pwp[tortoise][comparisonTortoise] / weightings[tortoise][comparisonTortoise] > 1e100) {
                            cout << "Extremely high value for tortoise " << tortoise << " and comparisonTortoise " << comparisonTortoise << " . Respective weightings: " << weightings[tortoise][comparisonTortoise] << " . Respective pwp: " << pwp[tortoise][comparisonTortoise] << endl;
                        }
                    } else {
                        cout << "Nonqualifying weighting: for tortoise " << tortoise << " and comparisonTortoise " << comparisonTortoise << " : " << weightings[tortoise][comparisonTortoise] << endl;
                        pwpOUT << "NA" << endl;
                    }
                }
            }
        }
        
        delete[] readCounts; // This holds ALL of the read counts for all loci
        /*
         {
         for ( N=1; N<=(NF/2); N++ ) {
         D[N] = ( $(2*N-1) + $(2*N)  );  # total coverage
         if ( D[N] > 0 ) {
         P[N] = $(2*N-1) / D[N];  # major allele freq
         if ( D[N] > 1 ) {
         WW = D[N]*(D[N]-1) ;
         W[N,N] += WW ;   # weights
         PI[N,N] += WW * ( 2 * P[N] * ( D[N] - $(2*N-1) ) / (D[N]-1) ) ;  # prob of difference
         }
         for ( M=1; M<N; M++ ) {
         if (D[M] > 0) {
         WW = D[N]*D[M] ;
         W[N,M] += WW;
         PI[N,M] += WW * ( P[N]*(1-P[M]) + P[M]*(1-P[N]) ) ;
         }
         }
         }
         }
         }
         
         END {
         for ( N=1; N<=(NF/2); N++ ) {
         for ( M=1; M<=N; M++ ) {
         print ( W[N,M]>0 ? PI[N,M]/W[N,M] : "NA" );
         }
         }
         }
         */
        
    }
    else cout << "Unable to open file";
    
    return 0;
}


