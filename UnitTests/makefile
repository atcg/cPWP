all:
	g++ -Wall main.cpp ../cPWP/generateSimulatedData.cpp ../cPWP/bamsToBin.cpp ../cPWP/calcPWPchunks.cpp ../cPWP/calcCOVARchunks.cpp -o unitTests -lboost_iostreams -pthread -std=c++0x

clean:
	rm -rf angsdOut* ind* unitTests pwptest.txt simulatedReferenceGenome.fasta* simulatedReferenceGenomeMutated.fasta normal* mutated* testingOut.pwp het* covarRef* simulatedReferenceGenome10mil* covarAngsdOut* testingOut.covar covarMutRef* covarInd* simulatedReferenceGenomeMutatedRef*