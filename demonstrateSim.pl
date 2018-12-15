#!/usr/bin/perl

use strict;
use warnings;
use Math::Random qw(random_binomial random_poisson);
use Data::Dumper;
use Getopt::Long;

my $numSites;
my $fracPolymorphic;
my $fracHet;
my $outFile;
my $meanReadDepth;
my $help;

GetOptions  ("sites=i"      => \$numSites,
             "poly=f"       => \$fracPolymorphic,
             "het=s"        => \$fracHet,
             "out=s"        => \$outFile,
             "depth=f"      => \$meanReadDepth,
             "help|man"     => \$help) || die "Couldn't get options with GetOpt::Long: $!\n";

if (!$numSites or !$fracPolymorphic or !$fracHet or !$outFile or !$meanReadDepth or $help) {
    die "Must supply --sites, --poly, --het, --samps, --out, and --depth.
    --sites is the total number of sites to simulate
    --poly is the fraction of segregating sites
    --het is  per-sample heterozygosities, separated by commas (and no whitespace)
    --samps is the number of individuals simulated
    --out is the name of the binary output file that holds the readcount
    --depth is the average per-sample read depth";
}

my $adjustedHeterozygosity = $fracHet / $fracPolymorphic; # This is the fraction of polymorphic sites that are het for samples

die if ($fracHet > $fracPolymorphic);


# Generate a list of 0s and 1s. The 1s are polymorphic sites. Binomial draw with $fracPolymorphic chance of being polymorphic
my @sites = random_binomial($numSites, 1, $fracPolymorphic);
my @hetSitesA;
my @hetSitesB;

my $hetCounterA = 0;
my $hetCounterB = 0;

my $numberOfPolymorphicSites = 0;
foreach my $site (@sites) {
    if ($site == 0) {
        # Nonpolymorphic, so can't be a het
        push(@hetSitesA, 0);
        push(@hetSitesB, 0);
    } else {
        # This is a polymorphic site, so it could be het in either or both samples
        $numberOfPolymorphicSites++;
        # Draw for sample A
        if (random_binomial(1,1,$adjustedHeterozygosity) == 1) {
            push(@hetSitesA, 1);
            $hetCounterA++;
        } else {
            push(@hetSitesA, 0);
        }
        # Draw for sample B
        if (random_binomial(1,1,$adjustedHeterozygosity) == 1) {
            push(@hetSitesB, 1);
            $hetCounterB++;
        } else {
            push(@hetSitesB, 0);
        }
    }
}
print "Empirical heterozygosity for Sample 1: " . $hetCounterA / $numSites . "\n";
print "Empirical heterozygosity for Sample 2: " . $hetCounterB / $numSites . "\n";

### Generate per-site coverages
my @coveragesA;
my @coveragesB;
foreach my $site (@sites) {
    # Poisson-distribution with a mean of 1
    push(@coveragesA, random_poisson(1,$meanReadDepth));
    push(@coveragesB, random_poisson(1,$meanReadDepth));
}

### Generate per-site tracked-allele depths
my @alleleDepthsA;
my @alleleDepthsB;

# Sample 1
my $coverageIndex = 0;
foreach my $hetSiteA (@hetSitesA) {
    if ($hetSiteA == 1) {
        push(@alleleDepthsA, random_binomial(1, $coveragesA[$coverageIndex], .5))

    } else {
        push(@alleleDepthsA, 0);
    }
    $coverageIndex++;
}
# Sample 2
$coverageIndex = 0;
foreach my $hetSiteB (@hetSitesB) {
    if ($hetSiteB == 1) {
        push(@alleleDepthsB, random_binomial(1, $coveragesB[$coverageIndex], .5))

    } else {
        push(@alleleDepthsB, 0);
    }
    $coverageIndex++;
}

### Now print out the results in the right format
open(my $outFH, ">", $outFile) or die "Couldn't open $outFile for writing: $!\n";
my $counter = 0;
my $fixedDiff = 0;
my $hetOneFixOther = 0;
my $actuallyNotPolymorphic = 0;
my $hetInBoth = 0;

foreach my $site (@sites) {

    if ($site == 0) {
        print $outFH pack("C", $coveragesA[$counter]); # Sample 1 allele A depth as an 8-bit char
        print $outFH pack("C", 0); # Sample 1 allele B depth as an 8-bit char
        print $outFH pack("C", $coveragesB[$counter]); # Sample 2 allele A depth as an 8-bit char
        print $outFH pack("C", 0); # Sample 2 allele B depth as an 8-bit char
    } else {
        my $alleleSamp1;
        my $alleleSamp2;
        # This is a polymorphic site, so samples could be het. If they're NOT
        # het, then the sample needs to choose whether it's homozygous for
        # allele A or allele B
        if ($hetSitesA[$counter] == 1) {
            # This site is heterozygous for sample 1
            $alleleSamp1 = "het";
            print $outFH pack("C", $alleleDepthsA[$counter]); # Sample 1 allele A depth as an 8-bit char
            print $outFH pack("C", $coveragesA[$counter] - $alleleDepthsA[$counter]); # Sample 1 allele B depth as an 8-bit char      
        } else {
            # Choose A or B
            if (random_binomial(1, 1, .5) == 0) {
                $alleleSamp1 = "alleleA";
                print $outFH pack("C", $coveragesA[$counter]); # Sample 1 allele A depth as an 8-bit char
                print $outFH pack("C", 0); # Sample 1 allele B depth as an 8-bit char
            } else {
                $alleleSamp1 = "alleleB";
                print $outFH pack("C", 0); # Sample 1 allele A depth as an 8-bit char
                print $outFH pack("C", $coveragesA[$counter]); # Sample 1 allele B depth as an 8-bit char
            }
        }

        if ($hetSitesB[$counter] == 1) {
            # This site is heterozygous for sample 2
            $alleleSamp2 = "het";
            print $outFH pack("C", $alleleDepthsB[$counter]); # Sample 1 allele A depth as an 8-bit char
            print $outFH pack("C", $coveragesB[$counter] - $alleleDepthsB[$counter]); # Sample 1 allele B depth as an 8-bit char
        } else {
            # Choose A or B
            if (random_binomial(1, 1, .5) == 0) {
                $alleleSamp2 = "alleleA";
                print $outFH pack("C", $coveragesB[$counter]); # Sample 2 allele A depth as an 8-bit char
                print $outFH pack("C", 0); # Sample 1 allele B depth as an 8-bit char
            } else {
                $alleleSamp2 = "alleleB";
                print $outFH pack("C", 0); # Sample 1 allele A depth as an 8-bit char
                print $outFH pack("C", $coveragesB[$counter]); # Sample 2 allele B depth as an 8-bit char
            }
        }

        # We'll now count up the different categories of differences between sample 1 and sample 2.
        # These will be used to calculate the empirical pairwise pi values
        if ($alleleSamp1 eq "alleleA" and $alleleSamp2 eq "alleleB") {
            $fixedDiff++;
        }
        if ($alleleSamp1 eq "alleleB" and $alleleSamp2 eq "alleleA") {
            $fixedDiff++;
        }
        if ($alleleSamp1 eq "alleleA" and $alleleSamp2 eq "het") {
            $hetOneFixOther++;
        }
        if ($alleleSamp1 eq "alleleB" and $alleleSamp2 eq "het") {
            $hetOneFixOther++;
        }
        if ($alleleSamp1 eq "het" and $alleleSamp2 eq "alleleA") {
            $hetOneFixOther++;
        }
        if ($alleleSamp1 eq "het" and $alleleSamp2 eq "alleleB") {
            $hetOneFixOther++;
        }
        if ($alleleSamp1 eq "alleleA" and $alleleSamp2 eq "alleleA") {
            $actuallyNotPolymorphic++;
        }
        if ($alleleSamp1 eq "alleleB" and $alleleSamp2 eq "alleleB") {
            $actuallyNotPolymorphic++;
        }
        if ($alleleSamp1 eq "het" and $alleleSamp2 eq "het") {
            $hetInBoth++;
        }

    }
    $counter++;
}
print $numSites . " total sites\n";
my $numberActuallyPolymorphic = $numberOfPolymorphicSites - $actuallyNotPolymorphic;
print $numberActuallyPolymorphic . " total polymorphic sites\n";
print $hetOneFixOther . " sites that are fixed in one sample and heterozygous in the other\n";
print $fixedDiff . " sites that are fixed different between the two samples\n";
print $hetInBoth . " sites that are het in both samples\n";

print "Self PWP A (fraction of het sites / 2): " . ($hetCounterA / $numSites) / 2 . "\n";
print "Self PWP B (fraction of het sites / 2): " . ($hetCounterB / $numSites) / 2 . "\n";


my $pwp = ($fixedDiff / $numSites) + 0.5 * ($hetOneFixOther/$numSites) + 0.5 * ($hetInBoth/$numSites);
print "A vs B pwp: " . $pwp . "\n";
