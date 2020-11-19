#!/usr/bin/perl -w
use strict;

my $parameters = {
    "minRemaining" => 30,
};

my $usage = "USAGE:\n$0 <mir_sequences.fa> <mir_sequences gff> <mir_seqInfo.txt> <hairpinFolds.txt>\n";

my $outputFile = "mir_dataset.txt";

my $mirSeqFa = $ARGV[0] or die $usage;
my $mirSeqGff = $ARGV[1] or die $usage;
my $mirSeqInfoFile = $ARGV[2] or die $usage;
my $hairpinFoldsFile = $ARGV[3] or die $usage;

my $fullSeqs = loadFasta($mirSeqFa);
my($hairpins,$products) = readMirbaseGff3_withTests($mirSeqGff);
my($hairpinSeqs,$hairpinFolds) = readFoldFile($hairpinFoldsFile);
my $miRSeqInfo = loadMiRSeqInfoFile($mirSeqInfoFile);
my $cutInfo = locateDroshaDicerCuts($hairpins,$products,$hairpinFolds,$fullSeqs,$hairpinSeqs);
my $preparedSequences = prepareSequences($cutInfo,$miRSeqInfo,$fullSeqs,$parameters);
printOutput($preparedSequences,$outputFile);


sub printOutput {
    my($preparedSequences,$outputFile) = @_;
    open(MDS,">$outputFile") or die "failed to open $outputFile for writing\n";
    print MDS "\#speciesPrefix\tname\tid\tgenomicLocation\taltGenomicLocation\tproduct5p\tproduct3p\t5drosha5p\tdicer5p\tdicer3p\tdrosha3p\thpStart\thpStop\tleftBuffer\trightBuffer\tseqTag\tsequence\n";
    foreach my $data (@{$preparedSequences}) {
	my($speciesPrefix,$name,$newId,$extendedLocation,$altExtendedLocation,$drosha5p,$dicer5p,$dicer3p,$drosha3p,$hpStart,$hpStop,$leftBuffer,$rightBuffer,$seqTag,$sequence,$prodName5p,$prodName3p) = @{$data};
	print MDS "$speciesPrefix\t$name\t$newId\t$extendedLocation\t$altExtendedLocation\t$prodName5p\t$prodName3p\t$drosha5p\t$dicer5p\t$dicer3p\t$drosha3p\t$hpStart\t$hpStop\t$leftBuffer\t$rightBuffer\t$seqTag\t$sequence\n";
	#print "$speciesPrefix\t$name\t$newId\t$extendedLocation\t$altExtendedLocation\t$drosha5p\t$dicer5p\t$dicer3p\t$drosha3p\t$hpStart\t$hpStop\t$leftBuffer\t$rightBuffer\t$seqTag\t$sequence\n";
	#if ($drosha5p ne "-") {
	#    my($drosha5pStart,$drosha5pStop) = split(",",$drosha5p);
	#    my($dicer5pStart,$dicer5pStop) = split(",",$dicer5p);
	#    my($miR5pStart,$miR5pStop) = ($drosha5pStop,$dicer5pStart);
	#    if ($miR5pStart > 0) {
	#	print "".substr($sequence,$miR5pStart-1,$miR5pStop-$miR5pStart+1)."\n";
	#    } else {
	#	print "Warning: 5p miR for $name ($newId) hangs off edge of sequence\n";
	#    }
	#}
	#if ($dicer3p ne "-") {
	#    my($dicer3pStart,$dicer3pStop) = split(",",$dicer3p);
	#    my($drosha3pStart,$drosha3pStop) = split(",",$drosha3p);
	#    my($miR3pStart,$miR3pStop) = ($dicer3pStop,$drosha3pStart);
	#    if ($miR3pStart > 0) {
	#	print "".substr($sequence,$miR3pStart-1,$miR3pStop-$miR3pStart+1)."\n";
	#    } else {
	#	print "Warning: 3p miR for $name ($newId) hangs off edge of sequence\n";
	#    }
	#}     
    }
    close(MDS);
}

sub prepareSequences {
    my($cutInfo,$miRSeqInfo,$fullSeqs,$parameters) = @_;
    my @preparedSequences;
    my $minRemaining = $parameters->{"minRemaining"};
    print "Trimming away sequences with neighboring miRs:\na minimum buffer of $minRemaining will be set for each side\n";
    foreach my $seqTag (%{$fullSeqs}) {
	my $sequence = $fullSeqs->{$seqTag};
	foreach my $hpName (keys %{$cutInfo->{$seqTag}}) {
	    my($drosha5p,$dicer5p,$dicer3p,$drosha3p,$hpStart,$hpStop,$leftBuffer,$rightBuffer,$prodName5p,$prodName3p) = @{$cutInfo->{$seqTag}{$hpName}};
	    my($speciesPrefix,$name,$newId,$relHPLocation,$extendedLocation,$altExtendedLocation,$neighbors) = @{$miRSeqInfo->{$seqTag}};
	    my($trim3pEnd,$trim5pEnd) = (0,0);
	    #print "\n\n$name:\n";
	    foreach my $neighborInfo (@{$neighbors}) {
		my($neighborName,$neighborChrom,$neighborStart,$neighborStop) = @{$neighborInfo};
		if ($neighborStop < $hpStart && $neighborStop > $trim5pEnd) {
		    $trim5pEnd = $neighborStop;
		    #print "($hpStart..$hpStop): neighbor at ($neighborStart..$neighborStop).  trim5pEnd=$trim5pEnd\n";
		} elsif ($neighborStart > $hpStop && ($rightBuffer - ($neighborStart - $hpStop)) > $trim3pEnd) {
		    $trim3pEnd = $rightBuffer - ($neighborStart - $hpStop);
		    #print "($hpStart..$hpStop): neighbor at ($neighborStart..$neighborStop).  trim3pEnd=$trim3pEnd\n";
		} elsif ($neighborStart < $hpStart && $neighborStop > $hpStop) {
		    print "Warning: $neighborName ($neighborStart..$neighborStop) stradles $name ($hpStart..$hpStop). No adjustment made.\n";
		} elsif ($neighborStart > $hpStart && $neighborStop < $hpStop) {
		    print "Warning: $neighborName ($neighborStart..$neighborStop) is stradled by $name ($hpStart..$hpStop). No adjustment made.\n";
		} elsif ($neighborStart < $hpStart && $neighborStop > $hpStart) {
		    print "Warning: $neighborName ($neighborStart..$neighborStop) overlaps $name ($hpStart..$hpStop)\n";		    
		    $trim5pEnd = $leftBuffer;
		    #print "($hpStart..$hpStop): neighbor at ($neighborStart..$neighborStop).  trim5pEnd=$trim5pEnd\n";
		} elsif ($neighborStart < $hpStop && $neighborStop > $hpStop) {
		    print "Warning: $neighborName ($neighborStart..$neighborStop) overlaps $name ($hpStart..$hpStop)\n";
		    $trim3pEnd = $rightBuffer;
		    #print "($hpStart..$hpStop): neighbor at ($neighborStart..$neighborStop).  trim3pEnd=$trim3pEnd\n";
		}
		#adjusting so that at least minRemaining nucleotides are flanking the sequence
		if ($trim5pEnd > $leftBuffer - $minRemaining) {
		    $trim5pEnd = ($leftBuffer - $minRemaining > 0) ? $leftBuffer - $minRemaining : 0;
		    #print "adj for minRemaining: trim5pEnd = $trim5pEnd\n";
		}
		if ($trim3pEnd > $rightBuffer - $minRemaining) {
		    $trim3pEnd = ($rightBuffer - $minRemaining > 0) ? $rightBuffer - $minRemaining : 0;
		    #print "adj for minRemaining: trim3pEnd = $trim3pEnd\n";
		}
	    }
	    #updating Sequence
	    my $len = length($sequence);
	    my $newSeq = substr($sequence,$trim5pEnd,$len-$trim5pEnd-$trim3pEnd);
	    #updating drosha5p
	    my $newDrosha5p = "-";
	    if ($drosha5p ne "-") {
		my($drosha5pStart,$drosha5pStop) = split(",",$drosha5p);
		my($newDrosha5pStart,$newDrosha5pStop) = ($drosha5pStart-$trim5pEnd,$drosha5pStop-$trim5pEnd);
		$newDrosha5p = "$newDrosha5pStart,$newDrosha5pStop";
	    }
	    #updating dicer5p
	    my $newDicer5p = "-";
	    if ($dicer5p ne "-") {
		my($dicer5pStart,$dicer5pStop) = split(",",$dicer5p);
		my($newDicer5pStart,$newDicer5pStop) = ($dicer5pStart-$trim5pEnd,$dicer5pStop-$trim5pEnd);
		$newDicer5p = "$newDicer5pStart,$newDicer5pStop";
	    }
	    #updating dicer3p
	    my $newDicer3p = "-";
	    if ($dicer3p ne "-") {
		my($dicer3pStart,$dicer3pStop) = split(",",$dicer3p);
		my($newDicer3pStart,$newDicer3pStop) = ($dicer3pStart-$trim5pEnd,$dicer3pStop-$trim5pEnd);
		$newDicer3p = "$newDicer3pStart,$newDicer3pStop";
	    }
	    #updating drosha3p
	    my $newDrosha3p = "-";
	    if ($drosha3p ne "-") {
		my($drosha3pStart,$drosha3pStop) = split(",",$drosha3p);
		my($newDrosha3pStart,$newDrosha3pStop) = ($drosha3pStart-$trim5pEnd,$drosha3pStop-$trim5pEnd);
		$newDrosha3p = "$newDrosha3pStart,$newDrosha3pStop";
	    }
	    #updating hairpin start and stop locations
	    my $newHPStart = $hpStart - $trim5pEnd;
	    my $newHPStop = $hpStop - $trim5pEnd;
	    #updating left and right Buffers
	    my $newLeftBuffer = $leftBuffer - $trim5pEnd;
	    my $newRightBuffer = $rightBuffer - $trim3pEnd;
	    #my $newSeqLen = length($newSeq);
	    #print "---------\ndrosha5p\tdicer5p\tdicer3p\tdrosha3p\thpStart\thpStop\tleftBuffer\trightBuffer\tlen\n";
	    #print "$drosha5p\t$dicer5p\t$dicer3p\t$drosha3p\t$hpStart\t$hpStop\t$leftBuffer\t$rightBuffer\t$len\n";
	    #print "$newDrosha5p\t$newDicer5p\t$newDicer3p\t$newDrosha3p\t$newHPStart\t$newHPStop\t$newLeftBuffer\t$newRightBuffer\t$newSeqLen\n";	    
	    my $hpSeq = substr($sequence,$hpStart-1,$hpStop-$hpStart+1);
	    my $newHPSeq = substr($newSeq,$newHPStart-1,$newHPStop-$newHPStart+1);
	    #print "----------\n$hpSeq\n$newHPSeq\n";
	    if ($hpSeq ne $newHPSeq) {
		die "Error: Sequences are different\n";
	    }
	    push(@preparedSequences,[$speciesPrefix,$name,$newId,$extendedLocation,$altExtendedLocation,$newDrosha5p,$newDicer5p,$newDicer3p,$newDrosha3p,$newHPStart,$newHPStop,$newLeftBuffer,$newRightBuffer,$seqTag,$newSeq,$prodName5p,$prodName3p]);
	}
    }
    my @sortedPreparedSequences = sort {$a->[0] cmp $b->[0] || $a->[1] cmp $b->[1]} @preparedSequences;
    return \@sortedPreparedSequences;
}

sub locateDroshaDicerCuts {
    my($hairpins,$products,$hairpinFolds,$fullSeqs,$hairpinSeqs) = @_;
    my %cutInfo;
    foreach my $seqTag (keys %{$hairpins}) {
	foreach my $hairpinInfo (@{$hairpins->{$seqTag}}) {
	    my($hpStart,$hpStop,$hpStrand,$hpId,$hpName) = @{$hairpinInfo};
	    my $testHPSeq = substr($fullSeqs->{$seqTag},$hpStart-1,$hpStop-$hpStart+1);
	    #if ($testSeq ne $hairpinSeqs->{$hpName}) {
		#print "$hpName\n$testSeq\n$hairpinSeqs->{$hpName}\n";
	    #}
	    my($drosha5p,$dicer5p,$dicer3p,$drosha3p) = ('-','-','-','-');
	    my($prodName5p,$prodName3p) = ('-','-');
	    foreach my $productInfo (@{$products->{$hpId}}) {
		my($prodChrom,$prodStart,$prodStop,$prodStrand,$prodId,$prodName) = @{$productInfo};
		my($relProdStart,$relProdStop) = ($prodStart - $hpStart,$prodStop - $hpStart);
		if ($relProdStop < length($hairpinFolds->{$hpName})) {
		    my $testSeq = substr($fullSeqs->{$seqTag},$hpStart -1 + $relProdStart,$prodStop-$prodStart+1);
		    my $testSeq2 = substr($hairpinSeqs->{$hpName},$relProdStart,$relProdStop-$relProdStart+1);
		    my $foldSubString = substr($hairpinFolds->{$hpName},$relProdStart,$relProdStop-$relProdStart+1) or die "$hpName";
		    my $arm = determineArm($foldSubString);
		    if ($arm eq "5p") {
			$drosha5p = "".($prodStart-1).",".$prodStart."";
			$dicer5p = "".$prodStop.",".($prodStop+1)."";
			$prodName5p = $prodName;
		    } elsif ($arm eq "3p") {
			$dicer3p = "".($prodStart-1).",".$prodStart."";
			$drosha3p = "".$prodStop.",".($prodStop+1)."";
			$prodName3p = $prodName;
		    } else {
			print "warning: couldn't determine arm for $hpName $prodName\n";
		    }
		    #print "$testSeq\n$testSeq2\n$foldSubString\n";
		    #print "$relProdStart..$relProdStop\n";
		} else {
		    print "warning: product ends after known length for $hpName\n";
		}
	    }
	    #print "$drosha5p\t$dicer5p\t$dicer3p\t$drosha3p\n";
	    my $leftBuffer = $hpStart - 1;
	    my $rightBuffer = length($fullSeqs->{$seqTag}) - $hpStop;
	    @{$cutInfo{$seqTag}{$hpName}} = ($drosha5p,$dicer5p,$dicer3p,$drosha3p,$hpStart,$hpStop,$leftBuffer,$rightBuffer,$prodName5p,$prodName3p);
	    #print "$seqTag\t$hpName\t$drosha5p\t$dicer5p\t$dicer3p\t$drosha3p\t$hpStart\t$hpStop\t$leftBuffer\t$rightBuffer\n";
	}
    }
    return \%cutInfo;
}

sub determineArm {
    my($foldSubString) = @_;
    my $leftCount = $foldSubString =~ tr/\(//;
    my $rightCount = $foldSubString =~ tr/\)//;
    my $dotCount = $foldSubString =~ tr/\.//;
    if ($leftCount && !($rightCount)) {
	return "5p";
    } elsif ($rightCount && !($leftCount)) {
	return "3p";
    }
    return "na";
}

sub readFoldFile {
    my($foldFile) = @_;
    my %sequences;
    my %folds;
    open(FF,$foldFile) or die "failed to open $foldFile\n";
    while (<FF>) {
	chomp;
	if ( /^>/ ) {
	    s/^>//;
	    my $name = $_;
	    my $seq = <FF>;
	    chomp($seq);
	    my $foldInfo = <FF>;
	    chomp($foldInfo);
	    my($fold,@other) = split(/\s+/,$foldInfo);
	    $sequences{$name} = $seq;
	    $folds{$name} = $fold;
	}
    }
    close(FF);
    return(\%sequences,\%folds);
}

sub loadFasta {
    my($fasta) = @_;
    my %seqs;
    open(MBHPF,$fasta) or die "failed to open $fasta\n";
    my $seqId = "";
    while (<MBHPF>) {
	chomp;
	my $line = $_;
	if ( $line =~ /^>/ ) {
	    $line =~ s/>//;
	    ($seqId) = split(/\s+/,$line);
	    unless ($seqs{$seqId}) {
		$seqs{$seqId} = "";
	    } else {
		print "Warning: $seqId already used as a key for another sequence in $fasta";
	    }
	} else {
	    $seqs{$seqId} .= $line;
	}
    }
    close(MBHPF);
    return \%seqs;
}

sub loadMiRSeqInfoFile {
    my($mirSeqInfoFile) = @_;
    my %miRSeqInfo;
    open(MSIF,$mirSeqInfoFile) or die "failed to open $mirSeqInfoFile for reading\n";
    while (<MSIF>) {
	chomp;
	unless (/ ^\#/ ) {
	    my($speciesPrefix,$newName,$name,$newId,$relHPLocation,$extendedLocation,$altExtendedLocation,$neighborField) = split(/\t/);
	    if ($miRSeqInfo{$newName}) {
		die "$newName already exists in miRSeqInfoHash";
	    }
	    #aae-mir-286a=aae-mir-2944b:-17..79;aae-mir-2944a=aae-mir-2944b:438..499;aae-mir-309b=aae-mir-
	    my @neighbors;
	    my @neighborInfoArray = split(';',$neighborField);
	    for my $neighborInfo (@neighborInfoArray) {
		my($neighborName,$neighborLoc) = split('=',$neighborInfo);
		my($neighborChrom,$neighborStartStop) = split(':',$neighborLoc);
		my($neighborStart,$neighborStop) = split(/\.\./,$neighborStartStop);
		push(@neighbors,[$neighborName,$neighborChrom,$neighborStart,$neighborStop]);
	    }
	    @{$miRSeqInfo{$newName}} = ($speciesPrefix,$name,$newId,$relHPLocation,$extendedLocation,$altExtendedLocation,\@neighbors);
	}
    }
    close(MSIF);
    return \%miRSeqInfo;
}

sub readMirbaseGff3_withTests {
    my($mirbaseGff3) = @_;
    my $hairpinTotal = 0;
    my $productTotal = 0;
    my %hairpins;
    my %products;
    my $lastHPId = "";
    open(MBGFF3,$mirbaseGff3) or die "could not open $mirbaseGff3\n";
    while(<MBGFF3>) {
	unless(/^\#/) {
	    chomp;
	    my($chrom,$source,$type,$start,$stop,$score,$strand,$phase,$info) = split(/\t/);
	    my %info;
	    my @terms = split(/;/,$info);
	    foreach my $term (@terms) {
		my($key,$value) = $term =~ /(.*)=(.*)/;
		$info{$key} = $value;
	    }
	    if($type eq "miRNA_primary_transcript") {
		# hairpin region from the gff file
		my $id = $info{ID} or die "No ID found for the line:\n$_\n";
		my $name = $info{Name} or die "No Name found for the line:\n$_\n";
		push(@{$hairpins{$chrom}},[$start,$stop,$strand,$id,$name]);
		$lastHPId = $id;
		$hairpinTotal += 1;
		if ($products{$id}) {
		    die "Error: in readMirbaseGff3_withTests() there may be more than hairpin identified by $id\n";
		}
	    }
	    if($type eq "miRNA") {
		# mature product from the hairpin file
		my $id = $info{ID} or die "No ID found for the line:\n$_\n";
                my $name = $info{Name} or die "No Name found for the line:\n$_\n";
		my $parentId = $info{Derives_from} or die "No Derives_from found for the line:\n$_\n";
		if ($parentId ne $lastHPId) {
		    my($testId) = split(/\_/,$lastHPId);
		    if ($parentId eq $testId) {
			print "Warning: product is derived from $parentId but last hairpin in gff was $lastHPId.  Using $lastHPId instead.\n";
			$parentId = $lastHPId;
		    } else {
			die "Error: last hairpin id in gff was $lastHPId and this product derives from $parentId.\n";
		    }
		}
		push(@{$products{$parentId}},[$chrom,$start,$stop,$strand,$id,$name]);
		$productTotal += 1;
	    }
	}
    }
    close(MBGFF3);
    return(\%hairpins,\%products,$hairpinTotal,$productTotal);
}

