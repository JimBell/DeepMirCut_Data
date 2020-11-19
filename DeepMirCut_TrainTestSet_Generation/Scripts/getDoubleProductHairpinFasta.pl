#!/usr/bin/perl -w
use strict;

my $usage = "USAGE:\n$0 <hairpins.fa> <mir dataset file> <output file>\n";

my $miRBaseHPFasta = $ARGV[0] or die $usage;
my $miRDataSetFile = $ARGV[1] or die $usage;
my $outputFasta = $ARGV[2] or die $usage;


my $hairpins = loadMiRBaseFasta($miRBaseHPFasta);

open(OPTF,">$outputFasta") or die "failed to open $outputFasta for writing\n";
open(MDSF,$miRDataSetFile) or die "failed to open $miRDataSetFile\n";
while (<MDSF>) {
    chomp;
    unless ( /^#/ ) {
	my($speciesPrefix,$name,$id,$genomicLocation,$altGenomicLocation,$product5p,$product3p,$drosha5p,$dicer5p,$dicer3p,$drosha3p,$hpStart,$hpStop,$leftBuffer,$rightBuffer,$seqTag,$sequence) = split(/\t/);
	if ($product5p ne '-' && $product3p ne '-' && $hairpins->{$name}) {
	    print OPTF ">$name\n$hairpins->{$name}\n";
	}
   }
}
close(MDSF);
close(OPTF);

sub loadMiRBaseFasta {
    my($miRBaseHPFasta) = @_;
    my %miRBaseSeqs;
    open(MBHPF,$miRBaseHPFasta) or die "failed to open $miRBaseHPFasta\n";
    my $seqId = "";
    while (<MBHPF>) {
	chomp;
	my $line = $_;
	if ( $line =~ /^>/ ) {
	    $line =~ s/>//;
	    ($seqId) = split(/\s+/,$line);
	    unless ($miRBaseSeqs{$seqId}) {
		$miRBaseSeqs{$seqId} = "";
	    } else {
		print "Warning: $seqId already used as a key for another sequence in $miRBaseHPFasta";
	    }
	} else {
	    $miRBaseSeqs{$seqId} .= $line;
	}
    }
    close(MBHPF);
    return \%miRBaseSeqs;
}
