#!/usr/bin/perl
use List::Util;
use strict;

my $parameters = {
    multiples => 9,
    minBuffer => 30,
    maxBuffer => 50,
    maxSize => 250,
    maxSimilar => 1000,
    trainOutputFile => "trainSet.txt",
    validateOutputFile => "validateSet.txt",
    testOutputFile => "testSet.txt",
    counter => 0,
};

my $usage = "USAGE:\n <mir_dataset> <trainSet> <cdHit File>\n";

my $mirDataSetFile = $ARGV[0] or die $usage;
my $trainSetFile = $ARGV[1] or die $usage;
my $CDHitClusterFile = $ARGV[2] or die $usage;

my $outputFile = $trainSetFile . "\_wSimilar.txt";
if ( $trainSetFile =~ /.*\.txt$/ ) {
    ($outputFile) = $trainSetFile =~ /(.*)\.txt$/;
    $outputFile = $outputFile . "\_wSimilar.txt";
}

my $mirData = loadMiRDataSet($mirDataSetFile);
$mirData = trimPolyNEnds($mirData);
my($trainSet) = loadTrainSet($trainSetFile);
my $cluster = readClusterFile($CDHitClusterFile);
my $augmentedSet = augmentWithSimilarSequences($mirData,$trainSet,$cluster,$parameters);

open(OPTF,">$outputFile") or die "failed to open $outputFile for writing\n";
foreach my $seqData (@{$augmentedSet}) {
    my($id,$name,$miRBaseId,$product5p,$product3p,$newDrosha5p,$newDicer5p,$newDicer3p,$newDrosha3p,$newHPStart,$newHPStop,$newSequence) = @{$seqData};
    print OPTF "$id\t$name\t$miRBaseId\t$product5p\t$product3p\t$newDrosha5p\t$newDicer5p\t$newDicer3p\t$newDrosha3p\t$newHPStart\t$newHPStop\t$newSequence\n";
}
close(OPTF);

sub augmentWithSimilarSequences {
    my($mirData,$trainSet,$cluster,$parameters) = @_;
    my @newTrainSet;
    foreach my $trainSetInfo (@{$trainSet}) {
	my($id,$name,$miRBaseId,$product5p,$product3p,$drosha5p,$dicer5p,$dicer3p,$drosha3p,$hpStart,$hpStop,$sequence) = @{$trainSetInfo};
	push(@newTrainSet,[$id,$name,$miRBaseId,$product5p,$product3p,$drosha5p,$dicer5p,$dicer3p,$drosha3p,$hpStart,$hpStop,$sequence]);
	my $counter = 1;
	if ($cluster->{$name}) {
	    my @clusterNames = List::Util::shuffle(@{$cluster->{$name}});
	    my $numGet = ($parameters->{maxSimilar} > @clusterNames) ? @clusterNames : $parameters->{maxSimilar};
	    for (my $itr = 0; $itr < $numGet; $itr++) {
		my $newId = $id."_s".$counter++;
		my $randSeqData = generateRandSeqData($mirData,$newId,$clusterNames[$itr],$parameters);
		push(@newTrainSet,$randSeqData);
	    }
	}
    }
    return \@newTrainSet;
}

sub generateRandSeqData {
    my($mirData,$id,$name,$parameters) = @_;
    my @data = @{$mirData->{$name}};
    my $dataEntry = $data[int(rand(@data))];
    my($speciesPrefix,$miRBaseId,$genomicLocation,$altGenomicLocation,$product5p,$product3p,$drosha5p,$dicer5p,$dicer3p,$drosha3p,$hpStart,$hpStop,$leftBuffer,$rightBuffer,$seqTag,$sequence) = @{$dataEntry};
    my($leftRemove,$rightRemove) = getRandomTrim($dataEntry,$parameters);
    my $newDrosha5p = updateCutSite($drosha5p,$leftRemove);
    my $newDicer5p = updateCutSite($dicer5p,$leftRemove);
    my $newDrosha3p = updateCutSite($drosha3p,$leftRemove);
    my $newDicer3p = updateCutSite($dicer3p,$leftRemove);
    my $newHPStart = $hpStart - $leftRemove;
    my $newHPStop = $hpStop - $leftRemove;
    my $newLeftBuffer = $leftBuffer - $leftRemove;
    my $newRightBuffer = $rightBuffer - $rightRemove;
    my $newSeqLength = length($sequence)-($leftRemove+$rightRemove);
    my $newSequence = substr($sequence,$leftRemove,$newSeqLength);
    if ($newSeqLength > $parameters->{maxSize}) {
	die "Error: Sequence Length = $newSeqLength is greater than maxSize parameter.\n";
    }
    unless (length($newSequence) == ($newLeftBuffer + $newRightBuffer + ($hpStop - $hpStart + 1)) ) {
	die "Error: error cutting sequence";
    }
    my @randSeqData = ($id,$name,$miRBaseId,$product5p,$product3p,$newDrosha5p,$newDicer5p,$newDicer3p,$newDrosha3p,$newHPStart,$newHPStop,$newSequence);
    return \@randSeqData;
}

sub generateShuffledSequences {
    my($dataSet) = @_;
    my @shuffledSequences;
    my $counter = 0;
    foreach my $seqData (@{$dataSet}) {
	my($id,$name,$miRBaseId,$product5p,$product3p,$newDrosha5p,$newDicer5p,$newDicer3p,$newDrosha3p,$newHPStart,$newHPStop,$sequence) = @{$seqData};
	$sequence =~ tr/uU/tT/;
	my $shuffledSeq = `python Scripts/altschulEriksonDinuclShuffle.py $sequence 1`;
	$shuffledSeq =~ tr/tT/uU/;
	chomp($shuffledSeq);
	push(@shuffledSequences,[$id."_shuffled",$name."_shuffled",$miRBaseId."_shuffled","-","-","-,-","-,-","-,-","-,-","-","-",$shuffledSeq]);
	$counter++;
	if ($counter % 1000 == 0) {
	    print $counter . "\n";
	}
    }
    return \@shuffledSequences;
}

sub getRandomTrim {
    my($dataEntry,$parameters) = @_;
    my $maxBuffer = $parameters->{maxBuffer};
    my $minBuffer = $parameters->{minBuffer};
    my($speciesPrefix,$miRBaseId,$genomicLocation,$altGenomicLocation,$product5p,$product3p,$drosha5p,$dicer5p,$dicer3p,$drosha3p,$hpStart,$hpStop,$leftBuffer,$rightBuffer,$seqTag,$sequence) = @{$dataEntry};
    my $seqLen = length($sequence);
    my $hpLength = $hpStop - $hpStart + 1;
    my($leftRemove,$rightRemove) = (0,0);
    if (int(rand(2))) {
	#trim left end first
	#print "leftSide\t$leftBuffer\t$rightBuffer\t$maxBuffer\t$hpLength\n";
	if ($maxBuffer + $hpLength > $parameters->{maxSize}) {
	    #print "$maxBuffer + $hpLength=".($maxBuffer + $hpLength)."\n";
	    $maxBuffer = $parameters->{maxSize} - $hpLength;
	    $minBuffer = $maxBuffer if ($maxBuffer < $minBuffer);
	    $maxBuffer = 0 if ($maxBuffer < 0);
	    $minBuffer = 0 if ($minBuffer < 0);
	    #print "leftTrim: $minBuffer\t$maxBuffer\n";
	}
	my $leftUpperLimit = ($maxBuffer + 1 < $leftBuffer) ? $maxBuffer + 1 : $leftBuffer;
	my $leftLowerLimit = ($minBuffer < $leftBuffer) ? $minBuffer : $leftBuffer;
	$leftRemove = $leftBuffer - (int(rand($leftUpperLimit - $leftLowerLimit)) + $leftLowerLimit);
	my $leftRemaining = $leftBuffer - $leftRemove;
	#print "leftRemaining = $leftRemaining\n";
	if ($maxBuffer + $hpLength + $leftRemaining > $parameters->{maxSize}) {
	    #print "$maxBuffer + $hpLength + $leftRemaining=".($maxBuffer + $hpLength + $leftRemaining)."\n";
	    $maxBuffer = $parameters->{maxSize} - $hpLength - $leftRemaining;
	    $minBuffer = $maxBuffer if ($maxBuffer < $minBuffer);
	    $maxBuffer = 0 if ($maxBuffer < 0);
	    $minBuffer = 0 if ($minBuffer < 0);
	    #print "rightTrim: $minBuffer\t$maxBuffer\n";
	}
	my $rightUpperLimit = ($maxBuffer + 1 < $rightBuffer) ? $maxBuffer + 1 : $rightBuffer;
	my $rightLowerLimit = ($minBuffer < $rightBuffer) ? $minBuffer : $rightBuffer;
	$rightRemove = $rightBuffer - (int(rand($rightUpperLimit - $rightLowerLimit)) + $rightLowerLimit);
    } else {
	#trim right end first
	#print "rightSide\t$leftBuffer\t$rightBuffer\t$maxBuffer\t$hpLength\n";
	if ($maxBuffer + $hpLength > $parameters->{maxSize}) {
	    $maxBuffer = $parameters->{maxSize} - $hpLength;
	    $minBuffer = $maxBuffer if ($maxBuffer < $minBuffer);
	    $maxBuffer = 0 if ($maxBuffer < 0);
	    $minBuffer = 0 if ($minBuffer < 0);
	    #print "rightTrim: $minBuffer\t$maxBuffer\n";
	}
	my $rightUpperLimit = ($maxBuffer + 1 < $rightBuffer) ? $maxBuffer + 1 : $rightBuffer;
	my $rightLowerLimit = ($minBuffer < $rightBuffer) ? $minBuffer : $rightBuffer;
	$rightRemove = $rightBuffer - (int(rand($rightUpperLimit - $rightLowerLimit)) + $rightLowerLimit);
	my $rightRemaining = $rightBuffer - $rightRemove;
	#print "rightRemaining = $rightRemaining\n";
	if ($maxBuffer + $hpLength + $rightRemaining > $parameters->{maxSize}) {
	    #print "$maxBuffer + $hpLength + $rightRemaining=".($maxBuffer + $hpLength + $rightRemaining)."\n";
	    $maxBuffer = $parameters->{maxSize} - $hpLength - $rightRemaining;
	    $minBuffer = $maxBuffer if ($maxBuffer < $minBuffer);
	    $maxBuffer = 0 if ($maxBuffer < 0);
	    $minBuffer = 0 if ($minBuffer < 0);
	    #print "leftTrim: $minBuffer\t$maxBuffer\n";
	}
	my $leftUpperLimit = ($maxBuffer + 1 < $leftBuffer) ? $maxBuffer + 1 : $leftBuffer;
	my $leftLowerLimit = ($minBuffer < $leftBuffer) ? $minBuffer : $leftBuffer;
	$leftRemove = $leftBuffer - (int(rand($leftUpperLimit - $leftLowerLimit)) + $leftLowerLimit);
    }
    return($leftRemove,$rightRemove);
}

sub updateCutSite {
    my($cutSite,$shiftLeft) = @_;
    my($cutSiteStart,$cutSiteStop) = split(",",$cutSite);
    my $newCutSiteStart = $cutSiteStart - $shiftLeft;
    my $newCutSiteStop = $cutSiteStop - $shiftLeft;
    return "$newCutSiteStart,$newCutSiteStop";
}

sub loadMiRDataSet {
    my($mirDataSetFile) = @_;
    my %mirData;
    open(MDSF,$mirDataSetFile) or die "failed to open $mirDataSetFile\n";
    while (<MDSF>) {
	chomp;
	unless ( /^#/ ) {
	    my($speciesPrefix,$name,$miRBaseId,$genomicLocation,$altGenomicLocation,$product5p,$product3p,$drosha5p,$dicer5p,$dicer3p,$drosha3p,$hpStart,$hpStop,$leftBuffer,$rightBuffer,$seqTag,$sequence) = split(/\t/);
	    push(@{$mirData{$name}},[$speciesPrefix,$miRBaseId,$genomicLocation,$altGenomicLocation,$product5p,$product3p,$drosha5p,$dicer5p,$dicer3p,$drosha3p,$hpStart,$hpStop,$leftBuffer,$rightBuffer,$seqTag,$sequence]);
	}
    }
    close(MDSF);
    return \%mirData;
}

sub readClusterFile {
    my($CDHitClusterFile) = @_;
    my %clusterHash;
    my $currentCluster;
    open(CDH,$CDHitClusterFile) or die "failed to open $CDHitClusterFile\n";
    while (<CDH>) {
	chomp;
	if ( /^>/ ) {
	    s/>//;
	    my $clusterName = $_;
	    $currentCluster = $clusterName;
	} else {
	    my($num,$aa,$name) = split(/\s+/);
	    ($name) = $name =~ /^>(.*?)\.\.\./;
	    $clusterHash{$currentCluster}{$name} = 1;
	}
    }
    close(CDH);
    my %cluster;
    foreach my $clusterName (keys %clusterHash) {
	foreach my $name (keys %{$clusterHash{$clusterName}}) {
	    foreach my $name2 (keys %{$clusterHash{$clusterName}}) {
		if ($name ne $name2) {
		    push(@{$cluster{$name}},$name2);
		}
	    }
	}
    }
    return \%cluster;
}

sub trimPolyNEnds {
    my($mirData) = @_;
    my %newMirData;
    foreach my $name (keys %{$mirData}) {
	foreach my $mirDataInfo (@{$mirData->{$name}}) {
	    my($speciesPrefix,$miRBaseId,$genomicLocation,$altGenomicLocation,
	       $product5p,$product3p,$drosha5p,$dicer5p,$dicer3p,$drosha3p,
	       $hpStart,$hpStop,$leftBuffer,$rightBuffer,$seqTag,$sequence) = @{$mirDataInfo};
	    my @seqChars = split('',$sequence);
	    my($leftRemove,$rightRemove) = (0,0);
	    for (my $itr = 0; $itr < $hpStart - 1; $itr++) {
		if (uc($seqChars[$itr]) eq "N") {
		    $leftRemove += 1;
		} else {
		    last;
		}
	    }
	    for (my $itr = @seqChars - 1; $itr > $hpStop - 1; $itr--) {
		if (uc($seqChars[$itr]) eq "N") {
		    $rightRemove += 1;
		} else {
		    last;
		}
	    }
	    my $newDrosha5p = updateCutSite($drosha5p,$leftRemove);
	    my $newDicer5p = updateCutSite($dicer5p,$leftRemove);
	    my $newDrosha3p = updateCutSite($drosha3p,$leftRemove);
	    my $newDicer3p = updateCutSite($dicer3p,$leftRemove);
	    my $newHPStart = $hpStart - $leftRemove;
	    my $newHPStop = $hpStop - $leftRemove;
	    my $newLeftBuffer = $leftBuffer - $leftRemove;
	    my $newRightBuffer = $rightBuffer - $rightRemove;
	    my $newSeqLength = length($sequence)-($leftRemove+$rightRemove);
	    my $newSequence = substr($sequence,$leftRemove,$newSeqLength);
	    unless (length($newSequence) == ($newLeftBuffer + $newRightBuffer + ($hpStop - $hpStart + 1)) ) {
		die "error cutting sequence";
	    }
	    push(@{$newMirData{$name}},[$speciesPrefix,$miRBaseId,$genomicLocation,$altGenomicLocation,
					$product5p,$product3p,$newDrosha5p,$newDicer5p,$newDicer3p,$newDrosha3p,
					$newHPStart,$newHPStop,$newLeftBuffer,$newRightBuffer,$seqTag,$newSequence]);
	}
    }
    return \%newMirData;
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
