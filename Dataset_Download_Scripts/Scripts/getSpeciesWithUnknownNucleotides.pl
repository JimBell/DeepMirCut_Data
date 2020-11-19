#!/usr/bin/perl -w
use strict;

my $usage = "USAGE:\n$0 <mir dataset file>\n";

my $mirDatasetFile = $ARGV[0] or die $usage;

my %species;

open(MDF,$mirDatasetFile) or die "failed to open $mirDatasetFile\n";
while (<MDF>) {
    chomp;
    unless ( /^#/ ) {
	my($speciesPrefix,$name,$id,$genomicLocation,$altGenomicLocation,$product5p,$product3p,$drosha5p,$dicer5p,$dicer3p,$drosha3p,$hpStart,$hpStop,$leftBuffer,$rightBuffer,$seqTag,$sequence)=split(/\t/);
	my(@seqChars) = split('',$sequence);
	foreach my $nuc (@seqChars) {
	    if (uc($nuc) ne 'A' && uc($nuc) ne 'C' && uc($nuc) ne 'G' && uc($nuc) ne 'T' && uc($nuc) ne 'U' && uc($nuc) ne 'N') {
		$species{$speciesPrefix} = 1;
	    }
	}
    }
}
close(MDF);

foreach my $speciesPrefix (keys %species) {
    print $speciesPrefix . "\n";
}
