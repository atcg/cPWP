all:
	g++ -Wall main.cpp ../cPWP/generateSimulatedData.cpp ../cPWP/bamsToBin.cpp ../cPWP/calcPWP.cpp -o cPWP -lboost_iostreams -pthread -std=c++0x

clean:
	rm -rf angsdOut* ind* cPWP pwptest.txt simulatedReferenceGenome.fasta* simulatedReferenceGenomeMutated.fasta normal* mutated* testingOut.pwp het*