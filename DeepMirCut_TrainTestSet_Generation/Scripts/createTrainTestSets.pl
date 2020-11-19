#!/usr/bin/perl
use List::Util;
use strict;

my $parameters = {
    minBuffer => 30,
    maxBuffer => 50,
    maxSize => 250,
    idPrefix => "ex",
    trainOutputFile => "trainSet.txt",
    validateOutputFile => "validationSet.txt",
    testOutputFile => "testSet.txt",
    TVTRatio => "8:1:1",
};

my $usage = "USAGE:\n <mir_dataset> <mir names file> <train:validate:test> <cdHit File> <output prefix>\n";

my $mirDataSetFile = $ARGV[0] or die $usage;
my $mirNamesFile = $ARGV[1] or die $usage;
my $outputPrefix = $ARGV[2];

if ($outputPrefix) {
    $parameters->{trainOutputFile} = $outputPrefix ."\_trainSet.txt";
    $parameters->{validateOutputFile} = $outputPrefix ."\_validationSet.txt";
    $parameters->{testOutputFile} = $outputPrefix ."\_testSet.txt";    
}

my $mirData = loadMiRDataSet($mirDataSetFile);
$mirData = trimPolyNEnds($mirData);
my $names = loadNames($mirNamesFile);
my($trainNames,$validateNames,$testNames) = divyNames($names,$parameters->{TVTRatio});
print "set counts:\n";
print "train set = " .  @{$trainNames} . "\n";
print "validate set = ".  @{$validateNames}. "\n";
print "test set = ".  @{$testNames} . "\n\n";
my $trainSet = generateRandSequenceSet($mirData,$trainNames,$parameters);
my $validationSet = generateRandSequenceSet($mirData,$validateNames,$parameters);
my $testSet = generateRandSequenceSet($mirData,$testNames,$parameters);
printSets($trainSet,$validationSet,$testSet,$parameters);

sub printSets {
    my($trainSet,$validationSet,$testSet,$parameters) = @_;
    my $trainOutputFile = $parameters->{trainOutputFile};
    my $validateOutputFile = $parameters->{validateOutputFile};
    my $testOutputFile = $parameters->{testOutputFile};
    open(TROUT,">$trainOutputFile") or die "failed to open $trainOutputFile for writing\n";
    foreach my $seqData (@{$trainSet}) {
	my($id,$name,$miRBaseId,$product5p,$product3p,$newDrosha5p,$newDicer5p,$newDicer3p,$newDrosha3p,$newHPStart,$newHPStop,$newSequence) = @{$seqData};
	print TROUT "$id\t$name\t$miRBaseId\t$product5p\t$product3p\t$newDrosha5p\t$newDicer5p\t$newDicer3p\t$newDrosha3p\t$newHPStart\t$newHPStop\t$newSequence\n";
    }
    close(TROUT);
    open(VOUT,">$validateOutputFile") or die "failed to open $validateOutputFile for writing\n";
    foreach my $seqData (@{$validationSet}) {
	my($id,$name,$miRBaseId,$product5p,$product3p,$newDrosha5p,$newDicer5p,$newDicer3p,$newDrosha3p,$newHPStart,$newHPStop,$newSequence) = @{$seqData};
	print VOUT "$id\t$name\t$miRBaseId\t$product5p\t$product3p\t$newDrosha5p\t$newDicer5p\t$newDicer3p\t$newDrosha3p\t$newHPStart\t$newHPStop\t$newSequence\n";
    }
    close(TSTOUT);
    open(TSTOUT,">$testOutputFile") or die "failed to open $testOutputFile for writing\n";
    foreach my $seqData (@{$testSet}) {
	my($id,$name,$miRBaseId,$product5p,$product3p,$newDrosha5p,$newDicer5p,$newDicer3p,$newDrosha3p,$newHPStart,$newHPStop,$newSequence) = @{$seqData};
	print TSTOUT "$id\t$name\t$miRBaseId\t$product5p\t$product3p\t$newDrosha5p\t$newDicer5p\t$newDicer3p\t$newDrosha3p\t$newHPStart\t$newHPStop\t$newSequence\n";
    }
    close(TSTOUT);
}

sub generateRandSequenceSet {
    my($mirData,$setNames,$parameters) = @_;
    my @randSequenceSet;
    foreach my $name (@{$setNames}) {
	my $id = $parameters->{idPrefix}.$parameters->{counter}++;
	my $randSeqData = generateRandSeqData($mirData,$id,$name,$parameters);
	push(@randSequenceSet,$randSeqData);
    }
    return \@randSequenceSet;
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
    my @shuffledNames = List::Util::shuffle(@names);    
    return \@shuffledNames;
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

sub divyNames {
    my($names,$TVTRatio) = @_;
    unless ($TVTRatio =~ /^\d+:\d+:\d+$/) {
	die "train:validate:test ratio must be in the form dd:dd:dd\n";
    } 
    my($trainPart,$validatePart,$testPart) = split(":",$TVTRatio);
    my $trainPercent = $trainPart / ($trainPart + $validatePart + $testPart);
    my $validatePercent = $validatePart / ($trainPart + $validatePart + $testPart);
    my $testPercent = $testPart / ($trainPart + $validatePart + $testPart);
    my @newNames = List::Util::shuffle(@{$names});
    my(@trainNames,@validateNames,@testNames);
    for (my $itr; $itr < int($trainPercent * @{$names} + 0.5); $itr++) {
	push(@trainNames,pop(@newNames));
    }
    for (my $itr; $itr < int($validatePercent * @{$names} + 0.5); $itr++) {
	push(@validateNames,pop(@newNames));
    }
    return(\@trainNames,\@validateNames,\@newNames);
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
