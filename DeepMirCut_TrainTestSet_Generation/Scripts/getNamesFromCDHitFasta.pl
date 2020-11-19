#!/usr/bin/perl -w
use strict;

my $usage = "USAGE:\n$0 <CDHit fasta> <output file>\n";

my $fastaFile = $ARGV[0] or die $usage;
my $outputFile = $ARGV[1] or die $usage;

my $hairpins = loadMiRBaseFasta($fastaFile);

open(OPTF,">$outputFile") or die "failed to open $outputFile for writing\n";
foreach my $name (keys %{$hairpins}) {
    print OPTF "$name\n";
}
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
