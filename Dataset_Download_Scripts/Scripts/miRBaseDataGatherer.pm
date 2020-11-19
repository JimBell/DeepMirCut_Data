package miRBaseDataGatherer;
use strict;

my $defaultParameters = {
    "indexed" => 1,
    "buffer" => 300,
};

sub getDefaultParameters {
    my($speciesPrefix,$outDir) = @_;
    unless ($outDir) {
	#outDir is null
	$outDir = ".";
    }
    if ($outDir eq "") {
	#outDir is empty string
	$outDir = ".";
    }
    $defaultParameters->{outFasta} = $outDir . "/" . $speciesPrefix . "\_sequences\.fa";
    $defaultParameters->{outGff} = $outDir . "/" . $speciesPrefix . "\_sequences\.gff";
    $defaultParameters->{seqInfoFile} = $outDir . "/" . $speciesPrefix . "\_seqInfo\.txt";
    $defaultParameters->{runInfoFile} = $outDir . "/" . $speciesPrefix . "\_runInfo\.txt";
    return $defaultParameters;
}

sub setParameters {
    my($parameters,$newParameters) = @_;
    foreach my $key (%{$newParameters}) {
	$parameters->{$key} = $newParameters->{$key};
    }
    return $parameters;
}

sub parseLocation {
    my($location)=@_;
    my($chrom,$start,$end,$strand);    
    if($location =~ /(.*)\:(-?\d+)\-(-?\d+)\:(.*)/) {
        $chrom=$1;
        $start=$2;
        $end=$3;
        $strand=$4;
    } elsif($location =~ /(.*)\:(-?\d+)\-(-?\d+)/) {
        $chrom=$1;
        $start=$2;
        $end=$3;
        $strand="+";
    } elsif($location =~ /(.*)\:(-?\d+)\.\.(-?\d+)\:(.*)/) {
	$chrom=$1;
        $start=$2;
        $end=$3;
        $strand=$4;   	
    } elsif($location =~ /(.*)\:(-?\d+)\.\.(-?\d+)/) {
	$chrom=$1;
        $start=$2;
        $end=$3;
        $strand="+";   	
    }
    return ($chrom,$start,$end,$strand);
}    

sub getOverlap {
    my($relStart1,$relStop1,$relStart2,$relStop2) = @_;
    die "relStart1 > relStop1 ($relStart1 > $relStop1) in getOverlap()\n" if ($relStart1 > $relStop1);
    die "relStart2 > relStop2 ($relStart2 > $relStop2) in getOverlap()\n" if ($relStart2 > $relStop2);
    if(($relStart1 <= $relStart2)&&($relStop2 <= $relStop1)) {
        #relstart2..relstop2 is within relstart1..relstop1
        return $relStop2 - $relStart2 + 1;
    } elsif (($relStart2 <= $relStart1)&&($relStop1 <= $relStop2)) {
        #relstart1..relstop1 is within relstart2..relstop2
        return $relStop1 - $relStart1 + 1;
    } elsif(($relStart1 <= $relStart2)&&($relStart2 <= $relStop1)) {
        #relStart2..relstop2 overlaps on right of relstart1..relstop1 
        return $relStop1 - $relStart2 + 1;
    } elsif(($relStart1 <= $relStop2)&&($relStop2 <= $relStop1)) {
        #relstart2..relStop2 overlaps on left of relstart1..relstop1
        return $relStop2 - $relStart1 + 1;
    } elsif(($relStart2 <= $relStart1)&&($relStart1 <= $relStop2)) {
        #relstart1..relstop1 overlaps on right of relstart2..relstop2
        return $relStop2 - $relStart1 + 1;
    } elsif(($relStart2 <= $relStop1)&&($relStop1<=$relStop2)) {
        #relstart1..relstop1 overlaps on left of relstart2..relstop2
        return $relStop2-$relStop1 + 1;
    }
    return 0;
}

sub extractRelSequence {
    my($seq,$relativeLocation) = @_;
    my($chrom,$start,$stop) = parseLocation($relativeLocation);
    $start = $start - 1;  #converting to zero-based
    $stop = $stop - 1;    #converting to zero-based
    my $newSeq = substr($seq,$start,$stop-$start+1);
    return $newSeq;
}

sub reverseComplement {
# Returns the reverse complement of the input sequence.
    my($seq)=@_;
    $seq =~ tr/acgturykmbdhvACGTURYKMBDHV/tgcaayrmkvhdbTGCAAYRMKVHDB/;
    $seq=reverse($seq);
    return $seq;
}

sub extractSequenceUsingIndex {
    my($fastaFile,$location) = @_;
    my($chrom,$start,$stop,$strand) = parseLocation($location);
    my $cmd = "samtools faidx $fastaFile \'$chrom:$start-$stop\'";
    my @fastaOutput = `$cmd`;
    my $locationTag;
    my $sequence = "";
    for my $line (@fastaOutput) {
	chomp($line);
	if ($line =~ /^>/) {
	    ($locationTag) = $line =~ /^>(.*)$/;
	} else {
	    $sequence .= $line;
	}
    }
    if ($strand eq '-') {
	$sequence = reverseComplement($sequence);
    }
    $sequence =~ tr/atTucg/AUUUCG/;
    if ($sequence eq "") {
	die "Error: failed to find sequence for $location\n";
    }
    return $sequence;
}

sub convertToRelativeLocation {
    #returns relative 1-based coordinates.
    my($location,$extendedLocation) = @_;
    my($chrom,$start,$stop,$strand) = parseLocation($location);
    my($extChrom,$extStart,$extStop,$extStrand) = parseLocation($extendedLocation);
    unless ($chrom eq $extChrom) {
	die "cannot convert coords to relative location $extChrom != $chrom";
    }
    my $newStart = ($strand eq '+') ? $start - $extStart : -($stop - $extStop);
    my $newStop = ($strand eq '+') ? $stop - $extStart : -($start - $extStop);
    $newStart++;  #convert to 1-based
    $newStop++;   #convert to 1-based
    return "chrom:$newStart..$newStop";
}

sub extendSequences {
    my($location,$chromSizes,$parameters) = @_;
    my($chrom,$start,$stop,$strand) = parseLocation($location);
    my $bufferLength = $parameters->{"buffer"};
    my $newStart = ($start - $bufferLength > 0) ? $start - $bufferLength : 1;
    my $newStop = ($stop + $bufferLength < $chromSizes->{$chrom}) ? $stop + $bufferLength : $chromSizes->{$chrom};
    my $newLocation = "$chrom:$newStart..$newStop:$strand";
    return $newLocation;
}

sub getNeighboringMirs {
    my($hairpins,$chrom,$chromAliases,$relHPLocation,$extendedLocation) = @_;
    my @neighbors;
    my($extChrom,$extStart,$extStop,$extStrand) = parseLocation($extendedLocation);
    unless ($extChrom eq $chromAliases->{$chrom}) {
	die "Error: chroms not the same in getNeighboringMirs(): $extChrom vs $chromAliases->{$chrom}\n";
    }
    foreach my $hairpinInfo (@{$hairpins->{$chrom}}) {
	my($start,$stop,$strand,$id,$name) = @{$hairpinInfo};
	if ($strand eq $extStrand && getOverlap($start,$stop,$extStart,$extStop)) {
	    my $location = "$extChrom:$start..$stop:$strand";
	    my $neighborHPLocation = convertToRelativeLocation($location,$extendedLocation);
	    my($neighborHPChrom,$neighborHPStart,$neighborHPStop) = parseLocation($neighborHPLocation);
	    unless ($neighborHPLocation eq $relHPLocation) {
		push(@neighbors,[$name,$id,$neighborHPStart,$neighborHPLocation,$location]);
	    } 
	}
    }
    return \@neighbors;
}

sub getExtendedSequencesFromGff {
    my($fastaFile,$hairpins,$products,$chromSizes,$chromAliases,$parameters) = @_;
    my %extSeqs;
    #my %nameHash;
    foreach my $chrom (keys %{$hairpins}) {
	if ($chromAliases->{$chrom}) {
	    my $newChrom = $chromAliases->{$chrom};
	    foreach my $hairpinInfo (@{$hairpins->{$chrom}}) {
		my($start,$stop,$strand,$id,$name) = @{$hairpinInfo};
		my $location = "$newChrom:$start..$stop:$strand";
		my $extendedLocation = extendSequences($location,$chromSizes,$parameters);
		my $extSequence = extractSequenceUsingIndex($fastaFile,$extendedLocation);
		my $relHPLocation = convertToRelativeLocation($location,$extendedLocation);
		my $neighbors = getNeighboringMirs($hairpins,$chrom,$chromAliases,$relHPLocation,$extendedLocation);
		my @extSeqProds;
		unless (@{$products->{$id}}) {
		    print "Warning: no products found for $id\n";
		}
		foreach my $productInfo (@{$products->{$id}}) {
		    my($productChrom,$productStart,$productStop,$productStrand,$productId,$productName) = @{$productInfo};
		    unless ($chrom eq $productChrom && $strand eq $productStrand && getOverlap($start,$stop,$productStart,$productStop)) {
			print "WARNING: product $productName does not overlap precursor $id $name in gff\n";
		    }
		    my $productLocation = "$newChrom:$productStart..$productStop:$productStrand";
		    my $relProdLocation = convertToRelativeLocation($productLocation,$extendedLocation);
		    my($relProdChrom,$relProdStart,$relProdStop) = parseLocation($relProdLocation);
		    my $productSeq = extractRelSequence($extSequence,$relProdLocation);
		    push(@extSeqProds,[$productId,$productName,$productSeq,$relProdStart,$relProdLocation]);
		}
		push(@{$extSeqs{$newChrom}},[$id,$name,$extSequence,$start,$relHPLocation,$location,$extendedLocation,$chrom,$neighbors,\@extSeqProds]);
	    }
	} else {
	    die "Error: can't find alias for $chrom\n";
	}
    }
    return \%extSeqs;
}

sub maskUnknownNucleotides {
    my($seq,$seq2) = @_;
    my @seqChars = split('',$seq);
    my @seqChars2 = split('',$seq2);
    my(@newSeqChars,@newSeqChars2);
    my $length = (length($seq) < length($seq2)) ? length($seq) : length($seq2);
    for (my $itr = 0; $itr < $length; $itr++) {
	if (($seqChars[$itr] =~ /[RYSKWMBDHV]/) || ($seqChars2[$itr] =~ /[RYSKWMBDHV]/)) {
	    $newSeqChars[$itr] = "N";
	    $newSeqChars2[$itr] = "N";
	} else {
	    $newSeqChars[$itr] = $seqChars[$itr];
	    $newSeqChars2[$itr] = $seqChars2[$itr];
	}
    }
    my $newSeq = join('',@newSeqChars);
    my $newSeq2 = join('',@newSeqChars2);
    return($newSeq,$newSeq2);
}


sub getSequencesFromGff {
    my($fastaFile,$hairpins,$products,$chromAliases,$parameters) = @_;
    my %sequences; 
    my %productSequences;
    my %nameHash;
    my %prodNameHash;
    foreach my $chrom (keys %{$hairpins}) {
	if ($chromAliases->{$chrom}) {
	    my $newChrom = $chromAliases->{$chrom};
	    foreach my $hairpinInfo (@{$hairpins->{$chrom}}) {
		my($start,$stop,$strand,$id,$name) = @{$hairpinInfo};
		my $location = "$newChrom:$start..$stop:$strand";
		my $seq = extractSequenceUsingIndex($fastaFile,$location);
		if (($nameHash{$name}) && $sequences{$name} ne $seq) {
		    my($tSeq1,$tSeq2) = maskUnknownNucleotides($sequences{$name},$seq);
		    if (($nameHash{$name}) && $tSeq1 eq $tSeq2) {
			print "Warning: two versions of $name has two differnt atypical characters for nucleotides at certain positions:\n";
			print "$sequences{$name}\n$seq\n";
		    } else {
			die "Error: more than one sequence with name= $name in set\n";
		    }
		} else {
		    $nameHash{$name} = 1;
		}
		$sequences{$name} = $seq;
		foreach my $productInfo (@{$products->{$id}}) {
		    my($productChrom,$productStart,$productStop,$productStrand,$productId,$productName) = @{$productInfo};
		    unless ($chrom eq $productChrom && $strand eq $productStrand && getOverlap($start,$stop,$productStart,$productStop)) {
			print "WARNING: product $productName does not overlap precursor $id $name in gff\n";
		    }
		    my $productLocation = "$newChrom:$productStart..$productStop:$productStrand";
		    my $prodSeq = extractSequenceUsingIndex($fastaFile,$productLocation);
		    if (($prodNameHash{$productName}) && $productSequences{$productName} ne $prodSeq) {
			print "WARNING: more than one sequence with name=$productName in set\n";
			print "$productLocation\n";
			print "$productSequences{$productName}\t$prodSeq\n";
		    } else {
			$prodNameHash{$productName} = 1;
		    }
		    $productSequences{$productName} = $prodSeq;
		}
	    }
	} else {
	    die "Error: can't find alias for $chrom\n";
	}
    }
    return(\%sequences,\%productSequences);
}

sub printOutputFiles {
    my($extSeqs,$parameters) = @_;
    my %nameCount;
    my %idCount;
    my $hairpinTotal = 0;
    my $productTotal = 0;
    my $miRBaseVersion = $parameters->{"miRBaseVersion"};
    my $genomeFileLocation = $parameters->{"genomeFileLocation"};
    my $assemblyInfoFileLocation = $parameters->{"assemblyInfoFileLocation"};
    my $gffFileLocation = $parameters->{"gffFileLocation"};
    my $outDir = $parameters->{"outputDir"};
    my $speciesPrefix = $parameters->{"speciesId"};
    my $genomeId = $parameters->{"genomeId"};
    my $outFasta = $defaultParameters->{outFasta};
    my $outGff = $defaultParameters->{outGff};
    my $seqInfoFile = $defaultParameters->{seqInfoFile};
    my $runInfoFile = $defaultParameters->{runInfoFile};
    open(OFSTA,">$outFasta") or die "failed to open $outFasta for writing\n";
    open(OGFF,">$outGff") or die "failed to open $outGff for writing\n";
    open(SIFLE,">$seqInfoFile") or die "failed to open $seqInfoFile for writing\n";
    open(RIFLE,">$runInfoFile") or die "failed to open $runInfoFile for writing\n";
    print OGFF "##gff-version 3\n";
    foreach my $chrom (sort {$a eq $b} keys %{$extSeqs}) {
	foreach my $extSeqInfo (sort {$a->[3] <=> $b->[3]} @{$extSeqs->{$chrom}}) {
	    my($id,$name,$extSequence,$start,$relHPLocation,$location,$extendedLocation,$altChrom,$neighbors,$extSeqProds) = @{$extSeqInfo};
	    my $newName = $name;
	    if ($nameCount{$name}) {
		$nameCount{$name} += 1;
		$newName = "$name\_seq$nameCount{$name}";
		print "Warning: $name already in output fasta. Changing name to $newName\n";
	    } else {
		$nameCount{$name} = 1;
	    }
	    my($id2) = split(/\_/,$id);
	    my $newId = $id2;
	    if ($idCount{$id2}) {
		$idCount{$id2} += 1;
		$newId = "$id2\_$idCount{$id2}";
		print "Warning: $id2 already in output gff. Changing id to $newId\n";
	    } else {
		$idCount{$id2} = 1;
	    }
	    my($relChrom,$relStart,$relStop) = parseLocation($relHPLocation);
	    my($extChrom,$extStart,$extStop,$extStrand) = parseLocation($extendedLocation);
	    my $altExtendedLocation = "$altChrom:$extStart..$extStop:$extStrand";
	    $relHPLocation = "$newName:$relStart..$relStop";
	    #FastaFile
	    print OFSTA ">$newName\n";
	    print OFSTA "$extSequence\n";
	    #seqInfo
	    my $neighborField = "";
	    foreach my $neighborInfo (sort {$a->[2] <=> $b->[2]} @{$neighbors}) {
		my($neighborName,$neighborId,$relNeighborStart,$relNeighborLocation,$neighborLocation) = @{$neighborInfo};
		my($neighborChrom,$neighborStart,$neighborStop) = parseLocation($relNeighborLocation);
		$relNeighborLocation = "$newName:$neighborStart..$neighborStop";  #in case a new sequence name was applied
		$neighborField = ($neighborField eq "") ? "$neighborName=$relNeighborLocation" : $neighborField . ";$neighborName=$relNeighborLocation";
	    }
	    print SIFLE "$speciesPrefix\t$newName\t$name\t$newId\t$relHPLocation\t$extendedLocation\t$altExtendedLocation\t$neighborField\n";
	    #gffInfo
	    my $hpInfo = "ID=$newId;Alias=$newId;Name=$name";
	    print OGFF "$newName\t.\tmiRNA_primary_transcript\t$relStart\t$relStop\t.\t+\t.\t$hpInfo\n";
	    $hairpinTotal += 1;
	    foreach my $extSeqProdInfo (sort {$a->[3] <=> $b->[3]} @{$extSeqProds}) {
		my($productId,$productName,$productSeq,$relPStart,$relProdLocation) = @{$extSeqProdInfo};
		my($relProdChrom,$relProdStart,$relProdStop) = parseLocation($relProdLocation);
		my $prodInfo = "ID=$productId;Alias=$productId;Name=$productName;Derives_from=$newId;Seq=$productSeq";
		print OGFF "$newName\t.\tmiRNA\t$relProdStart\t$relProdStop\t.\t+\t.\t$prodInfo\n";
		$productTotal += 1;
	    }
	}
    }
    #runInfo
    print RIFLE "$speciesPrefix\t$miRBaseVersion\t$genomeId\t$genomeFileLocation\t$assemblyInfoFileLocation\t$gffFileLocation\n";
    #final Checks
    if ($parameters->{expectedHPCount} != $hairpinTotal) {
	print "Warning: total number of hairpins output is not equal to input from gff ($hairpinTotal vs. $parameters->{expectedHPCount})\n";
    }
    if ($parameters->{expectedProdCount} != $productTotal) {
	print "Warning: total number of products output is not equal to input from gff ($productTotal vs. $parameters->{expectedProdCount})\n";
    }
}

sub testExtendedSequences {
    my($extSeqs,$sequences,$testSeqs,$parameters) = @_;
    foreach my $chrom (sort {$a eq $b} keys %{$extSeqs}) {
	foreach my $extSeqInfo (sort {$a->[3] <=> $b->[3]} @{$extSeqs->{$chrom}}) {
	    my($id,$name,$extSequence,$start,$relHPLocation,$location,$extendedLocation,$altChrom,$neighbors,$extSeqProds) = @{$extSeqInfo};
	    unless ($sequences->{$name}) {
		die "Error: Extended Sequence for $name found but regular sequence was not\n";
	    }
	    my $extractedSeq = extractRelSequence($extSequence,$relHPLocation);
	    unless ($extractedSeq eq $sequences->{$name}) {
		my($tSeq1,$tSeq2) = maskUnknownNucleotides($extractedSeq,$sequences->{$name});
		if ($tSeq1 eq $tSeq2) {
		    print "Warning: sequence or extracted sequence have atypical nucleotides that don't match:\n";
		    print "$extractedSeq\n$sequences->{$name}\n";
		} else {
		    die "Error: Sequences not equal ". $extractedSeq . " vs ". $sequences->{$name} . "\n";
		}
	    }
	    if ($testSeqs->{$name}) {
		unless ($extractedSeq eq $testSeqs->{$name}) {
		    my($tSeq1,$tSeq2) = maskUnknownNucleotides($extractedSeq,$testSeqs->{$name});
		    if ($tSeq1 eq $tSeq2) {
			print "Warning: sequence or test sequence have atypical nucleotides that don't match:\n";
			print "$extractedSeq\n$testSeqs->{$name}\n";
		    } else {
			print "Warning: Sequences not equal to what is in mirBase for $name ". $extractedSeq . " vs ". $testSeqs->{$name} . "\n";
			#die "Error: Sequences not equal ($name) ". $extractedSeq . " vs ". $testSeqs->{$name} . "\n";
		    }
		}
	    } else {
		die "Error: miRBase test sequence not found for $name\n";
	    }
	    foreach my $neighborInfo (@{$neighbors}) {
		my($neighborName,$neighborId,$relNeighborHPStart,$relNeighborLocation,$neighborLocation) = @{$neighborInfo};
		my($relNeighborChrom,$relNeighborStart,$relNeighborStop) = parseLocation($relNeighborLocation);
		my $newRelNeighborStart = $relNeighborStart;
		my $newRelNeighborStop = $relNeighborStop;
		my $startCorrection = 0;
		my $stopCorrection = 0;
		if ($newRelNeighborStart < 1) {
		    $newRelNeighborStart = 1;
		    $startCorrection = $newRelNeighborStart - $relNeighborStart;
		}
		if ($newRelNeighborStop > length($extSequence)) {
		    $newRelNeighborStop = length($extSequence);
		    $stopCorrection = $newRelNeighborStop - $relNeighborStop;
		}
		my $newRelNeighborLocation = "$relNeighborChrom:$newRelNeighborStart..$newRelNeighborStop";
		my $neighborSeq = extractRelSequence($extSequence,$newRelNeighborLocation);
		my $testNeighborSeq = substr($sequences->{$neighborName}, $startCorrection, length($sequences->{$neighborName}) - $startCorrection + $stopCorrection);
		unless ($neighborSeq eq $testNeighborSeq) {
		    my($tSeq1,$tSeq2) = maskUnknownNucleotides($neighborSeq,$testNeighborSeq);
		    if ($tSeq1 eq $tSeq2) {
			print "Warning: neighbor sequence and test sequence have atypical nucleotides that don't match:\n";
			print "$neighborSeq\n$testNeighborSeq\n";
		    } else {
			die "Error: neighboring seq from $name called $neighborName not equal to what is expected:\n$neighborSeq\n$testNeighborSeq\n";
		    }
		}
	    }
	}
    }
}

sub testMatureProducts {
    my($extSeqs,$productSequences,$prodTestSeqs,$parameters) = @_;
    foreach my $chrom (sort {$a eq $b} keys %{$extSeqs}) {
	foreach my $extSeqInfo (sort {$a->[3] <=> $b->[3]} @{$extSeqs->{$chrom}}) {
	    my($id,$name,$extSequence,$start,$relHPLocation,$location,$extendedLocation,$altChrom,$neighbors,$extSeqProds) = @{$extSeqInfo};
	    foreach my $extSeqProdInfo (sort {$a->[3] <=> $b->[3]} @{$extSeqProds}) {
		my($productId,$productName,$productSeq,$relProdStart,$relProdLocation) = @{$extSeqProdInfo};
		my($relChrom,$relStart,$relStop) = parseLocation($relProdLocation);
		$relStart--;  #converting to zero based
		$relStop--;   #converting to zero based
		my $seq = substr($extSequence,$relStart,$relStop - $relStart + 1);
		unless ($seq eq $productSeq) {
		    my($tSeq1,$tSeq2) = maskUnknownNucleotides($seq,$productSeq);
		    if ($tSeq1 eq $tSeq2) {
			print "Warning: product sequence or test sequence have atypical nucleotides that don't match:\n";
			print "$seq\n$productSeq\n";
		    } else {
			die "Error: sequence doesn't match up with sequence for $productName at $relProdLocation ($seq vs $productSeq)\n";
		    }
		}
		unless ($seq eq $productSequences->{$productName}) {
		    print "WARNING: sequence doesn't match the sequence directly extracted for $productName.  ($seq vs $productSequences->{$productName})\n";
		}
		if ($prodTestSeqs->{$productName}) {
		    unless ($seq eq $prodTestSeqs->{$productName}) {
			print $productName . "\n";
			print "WARNING: sequence doesn't match the sequence within miRBase mature file for $productName.  ($seq vs $prodTestSeqs->{$productName})\n";
		    }
		} else {
		    print "Warning: $productName not found in mature miRBase file\n";
		}
	    }
	}
    }
}

sub generateMiRDataFiles {
    my($fastaFile,$hairpins,$products,$chromSizes,$chromAliases,$testSeqs,$prodTestSeqs,$parameters) = @_;
    my($sequences,$productSequences) = getSequencesFromGff($fastaFile,$hairpins,$products,$chromAliases,$parameters);
    my $extSeqs = getExtendedSequencesFromGff($fastaFile,$hairpins,$products,$chromSizes,$chromAliases,$parameters);
    testExtendedSequences($extSeqs,$sequences,$testSeqs,$parameters);
    testMatureProducts($extSeqs,$productSequences,$prodTestSeqs,$parameters);
    printOutputFiles($extSeqs,$parameters);
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

sub readMirbaseGff3 {
    my($mirbaseGff3) = @_;
    my %hairpins;
    my %products;
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
	    }
	    if($type eq "miRNA") {
		# mature product from the hairpin file
		my $id = $info{ID} or die "No ID found for the line:\n$_\n";
                my $name = $info{Name} or die "No Name found for the line:\n$_\n";
		my $parentId = $info{Derives_from} or die "No Derives_from found for the line:\n$_\n";
		push(@{$products{$parentId}},[$chrom,$start,$stop,$strand,$id,$name]);
	    }
	}
    }
    close(MBGFF3);
    return(\%hairpins,\%products);
}

1;
