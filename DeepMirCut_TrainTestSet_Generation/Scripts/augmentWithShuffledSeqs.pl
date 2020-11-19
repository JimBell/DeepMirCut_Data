#!/usr/bin/perl
use List::Util;
use strict;

my $parameters = {
    minBuffer => 30,
    maxBuffer => 50,
    maxSize => 250,
    maxSimilar => 1000,
    counter => 0,
};

my $usage = "USAGE:\n$0 <trainSet>\n";

my $trainSetFile = $ARGV[0] or die $usage;

my $outputFile = $trainSetFile . "\_wShuffled.txt";
if ( $trainSetFile =~ /.*\.txt$/ ) {
    ($outputFile) = $trainSetFile =~ /(.*)\.txt$/;
    $outputFile = $outputFile . "\_wShuffled.txt";
}

my($trainSet) = loadTrainSet($trainSetFile);
my $augmentedSet = augmentWithShuffledSequences($trainSet);

open(OPTF,">$outputFile") or die "failed to open $outputFile for writing\n";
foreach my $seqData (@{$augmentedSet}) {
    my($id,$name,$miRBaseId,$product5p,$product3p,$newDrosha5p,$newDicer5p,$newDicer3p,$newDrosha3p,$newHPStart,$newHPStop,$newSequence) = @{$seqData};
    print OPTF "$id\t$name\t$miRBaseId\t$product5p\t$product3p\t$newDrosha5p\t$newDicer5p\t$newDicer3p\t$newDrosha3p\t$newHPStart\t$newHPStop\t$newSequence\n";
}
close(OPTF);

sub augmentWithShuffledSequences {
    my($trainSet) = @_;
    my @augmentedSet;
    my $counter = 0;
    foreach my $seqData (@{$trainSet}) {
	my($id,$name,$miRBaseId,$product5p,$product3p,$newDrosha5p,$newDicer5p,$newDicer3p,$newDrosha3p,$newHPStart,$newHPStop,$sequence) = @{$seqData};
	push(@augmentedSet,[$id,$name,$miRBaseId,$product5p,$product3p,$newDrosha5p,$newDicer5p,$newDicer3p,$newDrosha3p,$newHPStart,$newHPStop,$sequence]);
	$sequence = prepareSeqForShuffle($sequence);
	my $shuffledSeq = `python Scripts/altschulEriksonDinuclShuffle.py $sequence 1`;
	$shuffledSeq =~ tr/tT/uU/;
	chomp($shuffledSeq);
	push(@augmentedSet,[$id."_shuffled",$name."_shuffled",$miRBaseId."_shuffled","-","-","-,-","-,-","-,-","-,-","-","-",$shuffledSeq]);
	$counter++;
    }
    return \@augmentedSet;
}

sub prepareSeqForShuffle {
    my($sequence) = @_;
    my @validChars = ('A','C','G','T');
    my $validHash = {'A' => 1,'C' => 1,'G' => 1,'T' => 1};
    my @seqChars = split('',$sequence);
    for (my $itr = 0; $itr < @seqChars; $itr++) {
	if ($seqChars[$itr] ne uc($seqChars[$itr])) {
	    $seqChars[$itr] = uc($seqChars[$itr]);
	}
	if ($seqChars[$itr] eq "U") {
	    $seqChars[$itr] = "T"
	}
	unless ($validHash->{$seqChars[$itr]}) {
	    $seqChars[$itr] = $validChars[int(rand(@validChars))];
	}
    }
    my $newSequence = join('',@seqChars);
    return $newSequence;
}

sub loadTrainSet {
    my($trainSetFile) = @_;
    my @trainSet;
    open(TSF,$trainSetFile) or die "failed to open $trainSetFile\n";
    while (<TSF>) {
	chomp;
	unless ( /^#/ ) {
	    my($id,$name,$miRBaseId,$product5p,$product3p,$drosha5p,$dicer5p,$dicer3p,$drosha3p,$hpStart,$hpStop,$sequence) = split(/\t/);
	    push(@trainSet,[$id,$name,$miRBaseId,$product5p,$product3p,$drosha5p,$dicer5p,$dicer3p,$drosha3p,$hpStart,$hpStop,$sequence]);
	}
    }
    close(TSF);
    return(\@trainSet);
}
