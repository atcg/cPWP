#!/usr/bin/perl


# For our unbiased estimates of pairwise divergence, we want to run calculations
# on hundreds of millions of sites, where we have read counts for the major and
# minor alleles.
#
# This script does a few things. It takes as input the read counts dumped by
# ANGSD, finds the minor and major read counts for each individual at that locus (from the mafs file),
# then packs them into binary file representations of that data. The binary file
# format is as follows:
#    -Every read depth (integer) for the major or minor allele at
#     a locus is encoded as an 8-bit unsigned integer (so values can range from 0 to 255).
#     In this experiment, we don't want to include loci with read depth greater
#     than 255, because with coverage around 1X, such read depth would be an indication
#     that something is wrong. In fact, we will cap the maximum read depth to consider
#     a locus at a much lower number (10 or 15, set with $maxReadDepth)
#    -Every individual gets two integers in the output file for each locus, with
#     the major allele read count first, and the minor allele read count second
#        -As such, there are N*2 integers per locus, where N is the number of bams
#         (individuals) input into ANGSD, in the same order as the bam file list


use strict;
use warnings;
use List::Util qw(sum max);
use Data::Dumper;
use Getopt::Long;

my $MAF;
my $countsFile;
my $numIndividuals;
my $outFile;
my $maxReadDepth;
my $help;

GetOptions  ("maf=s"         => \$MAF,
             "counts=s"      => \$countsFile,
             "individuals=i" => \$numIndividuals,
             "out=s"         => \$outFile,
             "maxdepth=i"    => \$maxReadDepth,
             "help|man"      => \$help) || die "Couldn't get options with GetOpt::Long: $!\n";

if (!$MAF or !$countsFile or !$outFile or !$numIndividuals or !$maxReadDepth or ($maxReadDepth > 255) or $help) {
    die "Must supply --maf, --counts, --individuals, --out, and --maxdepth.
    --maf designates the mafs file output by ANGSD (uncompressed)
    --counts is for the counts file output by ANGSD (uncompressed)
    --individuals is the number of individuals input into ANGSD (the number of files in the bam file list)
    --out is the desired name of the output file that holds the binary representation of the readcounts.
    --maxdepth is the highest allowed number for a major or minor allele read count (to avoid overweighting repetitive regions.
The mafs and counts files should be unzipped first.
Also note that the maximum value for --maxdepth is 255";
}

my $locusCounter = 1;
my $fileCounter = 1;

my %alleleOrderHash = ("A" => 0,
                       "C" => 1,
                       "G" => 2,
                       "T" => 3
                       );

open(my $MAFfh, "<", $MAF) or die "Couldn't open $MAF for reading: $!\n";
open(my $countsFH, "<", $countsFile) or die "Couldn't open $countsFile for reading: $!\n";
open(my $outFH, ">", $outFile) or die "Couldn't open $outFile for writing: $!\n";

while ( not eof $countsFH and not eof $MAFfh ) {
    # Each iteration of this loop processes a single locus from the ANGSD output
    
    my $MAFline = <$MAFfh>;
    my $countsLine = <$countsFH>;
    chomp($countsLine);
    if ($countsLine =~ /^ind/) { # This is the header line of the counts file, so skip
        next;
    }
    
    my @mafFields = split(/\t/, $MAFline);
    my $majorAllele = $mafFields[2];
    my $minorAllele = $mafFields[3];
    
    if ($majorAllele !~ /[ATCG]/ or $minorAllele !~ /[ATCG]/) {
        die "Either the major allele ($majorAllele) or the minor allele ($minorAllele) \
        does not equal A, C, G, or T\n";
    }
    
    # Get the information on the read counts for all the individuals
    my @countsFields = split(/\t/, $countsLine);
    if (max(@countsFields) > $maxReadDepth) { # Make sure none of the read counts are higher than the max allowed
        next;
    }
    
    my $tortCounter = 0;
    while ($tortCounter < $numIndividuals) {
        my @individualDepths = splice(@countsFields,0,4); # Remove the first 4 elements of the array    
        my $majorDepth = $individualDepths[$alleleOrderHash{$majorAllele}];
        my $minorDepth = $individualDepths[$alleleOrderHash{$minorAllele}];
        
        # Uncomment the line below if you want a text output instead of binary
        #print $outFH "$majorDepth\t$minorDepth\t";
        
        print $outFH pack("C", $majorDepth);  # "C" = 8-bit integer
        print $outFH pack("C", $minorDepth);  # "C" = 8-bit integer 
        $tortCounter++;   
    }
    $locusCounter++;
}
