#!/usr/bin/perl -w
use strict;

my $usage = "USAGE:\n$0 <mir dataset file> <output file>\n";

my $miRDataSetFile = $ARGV[0] or die $usage;
my $outputFile = $ARGV[1] or die $usage;

open(OPTF,">$outputFile") or die "failed to open $outputFile for writing\n";
open(MDSF,$miRDataSetFile) or die "failed to open $miRDataSetFile\n";
while (<MDSF>) {
    chomp;
    unless ( /^#/ ) {
	my($speciesPrefix,$name,$id,$genomicLocation,$altGenomicLocation,$product5p,$product3p,$drosha5p,$dicer5p,$dicer3p,$drosha3p,$hpStart,$hpStop,$leftBuffer,$rightBuffer,$seqTag,$sequence) = split(/\t/);
	if ($product5p ne '-' && $product3p ne '-') {
	    print OPTF "$name\n";
	}
   }
}
close(MDSF);
close(OPTF);
