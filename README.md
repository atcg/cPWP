Overview:
---------
This software can be used to calculate genetic differentiation between
individual samples based on aligned short-read sequencing data. It is
particularly useful for calculating global genome-wide differentiation
between samples that have been sequenced at very low coverage, as it
does not rely on called genotypes at any locus. There are functions to
calculate both pairwise pi and covariance among samples.


Instructions:
-------------
Compilation has been tested with g++ 4.8.2, and requires C++11.

To download and compile:
```
git clone https://github.com/atcg/cPWP.git
cd cPWP/cPWP
make
```

That will create an executable file called cPWP in the current directory which we
will use to run the divergence calculation on the binary dataset we will create below.

Data must be supplied to cPWP in the form of binary unsigned chars, two per
individual for each locus. The first unsigned char is the number of reads
that support the major allele for that individual, and the second unsigned
char is the number of reads that suppor the minor allele for that individual.
The total size of the input file should be `num_individuals * num_loci * 2` bytes.

The software [ANGSD](http://popgen.dk/wiki/index.php/ANGSD) is very efficient
in parsing BAM files and outputting read counts at every locus for every individual
in text format. The output read counts file then needs to be processed to
discern which bases are the major and minor alleles, and the counts need to be
transformed into binary for input into cPWP. The included script
[convertANGSDreadcountsToBinary.pl](https://github.com/atcg/cPWP/blob/master/convertANGSDreadcountsToBinary.pl)
does both of those things.

To use it, you must first generate the read counts for each base with ANGSD, which
corresponds to the `-dumpCounts 4` argument. One such command might be:

`angsd0.613/angsd -bam bamlist.txt -out countsPrefix -doCounts 1 -dumpCounts 4 -doMajorMinor 2 -doMaf 8`

This will generate three output files: countsPrefix.counts.gz, countsPrefix.pos.gz, and
countsPrefix.mafs.gz. We will use countsPrefix.mafs.gz to determine which allele
is major and which is minor (this is created by supplying the -doMaf 8 argument),
and with that information we will pull the appropriate columns from each line of
countsPrefix.counts.gz. We'll unzip these first:

```
gunzip countsPrefix.counts.gz
gunzip countsPrefix.mafs.gz
```

We'll then feed both of these files to [convertANGSDreadcountsToBinary.pl](https://github.com/atcg/cPWP/blob/master/convertANGSDreadcountsToBinary.pl)
for conversion to binary like so:

```
perl convertANGSDreadcountsToBinary.pl --maf countsPrefix.mafs --counts countsPrefix.counts --individuals <integer> --maxdepth <integer> --out majorMinorCounts.binary
```

You must to give it the integer number of individuals (which is the same as the number of
bam files you included in your bamfile.txt list you gave to ANGSD), and also an integer
that represents the maximum read depth you wish to allow for any particular read count
at any one locus. This is included because the divergence calculation is weighted by the
read depth at each locus, and we wish to avoid overweighting repetitive regions that
are mistakenly included only once in the reference genome (and as such may have extremely
high read depths). **Note that the largest allowed number for this argument is 255**, as
the numbers must be represented as unsigned chars, which have maximum value of 255. If
the read depth for any allele at a particular locus for an individual is higher than
the value for maxdepth, the entire locus is skipped and not output to the binary readcounts
file.

The script will output a single file, in this case called `majorMinorCounts.binary`. 

A total of eight arguments must be supplied to cPWP:
1) Either "pwp" for pairwise pi or "covar" for covariance
2) The (full or relative) path to the binary read counts file
3) The number of loci to analyze. This should be the total number of loci in your read counts file, unless you want to use a smaller subset for testing
4) The number of individuals represented in the read counts file (this information is not encoded in the read counts file itself, so you need to supply this so that cPWP knows when a new locus begins)
5) A name for the output file to be created
6) The number of loci to be analyzed in a block per thread (see below).
7) A file containing the sample names, with one name per line. For instance, the bamlist.txt file from above would work here
8) The number of threads to use for computation

To run pwp calculations based on the data in the file, you can use the following
command:

```
cPWP pwp majorMinorCounts.binary <numLociToAnalyze> <numIndividuals> <outputFileName> <blockSize> <orderedNameList> <numThreads>
```

Similarly, to run covariance calculations based on the data in the file, you can use the following command:
```
cPWP covar majorMinorCounts.binary <numLociToAnalyze> <numIndividuals> <outputFileName> <blockSize> <orderedNameList> <numThreads>
```

So, if we wanted to calculate PWP over the first 800,000 loci in the majorMinorCounts.binary file and
we had 10 total individuals represented in that file, we could do something like:

```
./cPWP/cPWP/cPWP pwp majorMinorCounts.binary 800000 10 divergenceOutput.txt 40000 names.list 3
```

This command would take the first 800,000 sites in the binary file and break them
into chunks of 40,000 sites, deploying those chunks to three independent threads for
computation. Once all chunks across all threads are finished running, results across
all threads are pooled and the final divergence numbers are output to the specified
output file. Here, names.list is an ordered list of sample names in the same order as
the list of bams that was originally passed to ANGSD (you could actually just use the 
same bam list that you passed to ANGSD, but then your output columns 1 and 2 would look
like "sample1.bam" and "sample2.bam" instead of "sample1" and "sample2".

cPWP should experience a near-linear speed increase with increasing threads. To determine the chunk size 
and number of threads to use, first consider the amount of free RAM available. RAM usage will be equal 
to slightly more than the number of threads * chunk size * number of individuals * 2, in bytes. If you have more bytes
of free RAM than the total size of the binary input file (again, in bytes), then simply
make the chunk size equal to slightly less than the total number of loci to analyze divided by
the number of free threads you have available. Otherwise, reduce the chunk size to reduce
RAM usage. 

For example, let's say you have a machine with 32 cores and 8GB of RAM, with read count data for 100,000,000 
loci and 1,000 individuals. The binary read counts file would be 100000000*1000*2 = 200 billion bytes in size. If you 
want to use all 32 cores, then you need to set the chunk size so that 32 * 2 * chunk size * number of samples 
is smaller than 8GB (about 8 billion). So the absolute maximum chunk size would be about 8 billion / (32 * 2 * 1000) = 125,000.
In practice, you shouldn't try to maximize the chunk size. In the scenario above I'd set the chunk size to 100,000, for instance,
to leave some RAM available on the machine.

