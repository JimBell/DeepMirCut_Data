#!/usr/bin/perl
use List::Util;
use strict;

my $parameters = {
    mutationMultiples => 5,
    mutations => 2,
    minBuffer => 30,
    maxBuffer => 50,
    maxSize => 250,
    counter => 0,
};

my $usage = "USAGE:\n$0 <trainSet>\n";

my $trainSetFile = $ARGV[0] or die $usage;

my $outputFile = $trainSetFile . "\_mut".$parameters->{mutations}.".txt";
if ( $trainSetFile =~ /.*\.txt$/ ) {
    ($outputFile) = $trainSetFile =~ /(.*)\.txt$/;
    $outputFile = $outputFile . "\_mut".$parameters->{mutations}.".txt";
}

my($trainSet) = loadTrainSet($trainSetFile);
my $augmentedSet = augmentWithRandomMutations($trainSet,$parameters);

open(OPTF,">$outputFile") or die "failed to open $outputFile for writing\n";
foreach my $seqData (@{$augmentedSet}) {
    my($id,$name,$miRBaseId,$product5p,$product3p,$newDrosha5p,$newDicer5p,$newDicer3p,$newDrosha3p,$newHPStart,$newHPStop,$newSequence) = @{$seqData};
    print OPTF "$id\t$name\t$miRBaseId\t$product5p\t$product3p\t$newDrosha5p\t$newDicer5p\t$newDicer3p\t$newDrosha3p\t$newHPStart\t$newHPStop\t$newSequence\n";
}
close(OPTF);

sub augmentWithRandomMutations {
    my($trainSet,$parameters) = @_;
    my $multiples = $parameters->{mutationMultiples};
    my $numMutations = $parameters->{mutations};
    my @newTrainSet;
    foreach my $trainSetInfo (@{$trainSet}) {
	my($id,$name,$miRBaseId,$product5p,$product3p,$drosha5p,$dicer5p,$dicer3p,$drosha3p,$hpStart,$hpStop,$sequence) = @{$trainSetInfo};
	push(@newTrainSet,[$id,$name,$miRBaseId,$product5p,$product3p,$drosha5p,$dicer5p,$dicer3p,$drosha3p,$hpStart,$hpStop,$sequence]);
	for (my $itr = 1; $itr < $multiples; $itr++) {
	    my($newSeq,$mutationString) = introduceRandomMutations($sequence,$numMutations);
	    my $newId = $id."_mut".$itr."_".$mutationString;
	    push(@newTrainSet,[$newId,$name,$miRBaseId,$product5p,$product3p,$drosha5p,$dicer5p,$dicer3p,$drosha3p,$hpStart,$hpStop,$newSeq]);
	}
    }
    return \@newTrainSet;
}

sub introduceRandomMutations {
    my($sequence,$numMutations) = @_;
    my @mutationList;
    my @seqChars = split('',$sequence);
    my @idxList = (0..length($sequence)-1);
    my @nucleotides = ('A','C','G','U');
    for (1..$numMutations) {
	my $idxListIdx = int(rand(@idxList));
	my $idx = splice(@idxList,$idxListIdx,1);
	my $nucIdx = int(rand(@nucleotides-1));
	if ($nucleotides[$nucIdx] eq $seqChars[$idx]) {
	    $nucIdx = ++$nucIdx % @nucleotides;
	}
	push(@mutationList, [$idx , $nucleotides[$nucIdx]]);
	$seqChars[$idx] = $nucleotides[$nucIdx];
    }
    my $newSeq = join('',@seqChars);
    my $mutationString = "";
    foreach my $mutationInfo (sort {$a->[0] <=> $b->[0]} @mutationList) {
	my($idx,$nucleotide) = @{$mutationInfo};
	$mutationString .= "$idx$nucleotide"
    }
    return($newSeq,$mutationString);
}

sub loadNames {
    my($mirNamesFile) = @_;
    my @names;
    open(MNF,$mirNamesFile) or die "failed to open $mirNamesFile";
    while (<MNF>) {
	chomp;
	unless ( /^#/ ) {
	    my($name) = split(/\t/);
	    push(@names,$name);
	}
    }
    close(MNF);
    return \@names;
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
