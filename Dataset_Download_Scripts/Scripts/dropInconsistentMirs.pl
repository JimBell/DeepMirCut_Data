#!/usr/bin/perl -w 
use strict;

my $usage = "USAGE:\n$0 <mir_sequences.gff> <mir_seqInfo.txt> <mir_sequences.fa> <inconsistentMirs.txt>\n";

my $mirSeqGff = $ARGV[0] or die $usage;
my $mirSeqInfoFile = $ARGV[1] or die $usage;
my $mirSeqFasta = $ARGV[2] or die $usage;
my $inconsistentMirsFile = $ARGV[3] or die $usage;

my $outputSeqGffFile = "new_mir_sequences.gff";
my $outputSeqInfoFile = "new_mir_seqInfo.txt";
my $outputSeqFasta = "new_mir_sequences.fa";

my $inconsistentMirs = loadInconsistentMirsFile($inconsistentMirsFile);

my($hairpins,$products,$gffLines,$hairpinTotal,$productTotal,$linesTotal) = loadMirSeqGffLines($mirSeqGff);
my $mirSeqInfoLines = loadMirSeqInfoLines($mirSeqInfoFile);
my($sequences,$fastaEntries)=loadMirSeqFastaLines($mirSeqFasta);
my($newMirSeqInfoLines,$newGffLines,$newFastaEntries) = dropInconsistentMirs($hairpins,$gffLines,$mirSeqInfoLines,$fastaEntries);
printOutput($newMirSeqInfoLines,$newGffLines,$newFastaEntries,$outputSeqGffFile,$outputSeqInfoFile,$outputSeqFasta);

sub printOutput {
    my($newMirSeqInfoLines,$newGffLines,$newFastaEntries,$outputSeqGffFile,$outputSeqInfoFile,$outputSeqFasta) = @_;
    open(SINF,">$outputSeqInfoFile") or die "failed to open $outputSeqInfoFile for writing\n";
    for my $seqInfo (@{$newMirSeqInfoLines}) {
	my($species,$seqName,$name,$id,$relativeLocation,$genomicLocation,$altGenomicLocation,$neighbors) = @{$seqInfo};
	my @neighborEntries;
	for my $neighborInfo (@{$neighbors}) {
	    my($neighborName,$neighborLoc) = @{$neighborInfo};
	    push(@neighborEntries,"$neighborName=$neighborLoc");
	}
	print SINF "$species\t$seqName\t$name\t$id\t$relativeLocation\t$genomicLocation\t$altGenomicLocation";
	if (@neighborEntries) {
	    print SINF "\t".join(";",@neighborEntries)."\n";
	} else {
	    print SINF "\t".join(";",@neighborEntries)."\n";
	}
    }
    close(SINF);
    open(SGFF,">$outputSeqGffFile") or die "failed to open $outputSeqGffFile for writing\n";
    print SGFF "##gff-version 3\n";
    for my $gffLineInfo (@{$newGffLines}) {
	my($type,$seqName,$start,$stop,$strand,$id,$name,$parentId,$line) = @{$gffLineInfo};
	print SGFF "$line\n";
    }
    close(SGFF);
    open(SFA,">$outputSeqFasta") or die "failed to open $outputSeqFasta for writing\n";
    for my $fastaEntryInfo (@{$fastaEntries}) {
	my($seqName,$seq) = @{$fastaEntryInfo};
	print SFA ">$seqName\n$seq\n";
    }    
    close(SFA);
}

sub dropInconsistentMirs {
    my($hairpins,$gffLines,$mirSeqInfoLines,$fastaEntries) = @_;
    my @newMirSeqInfoLines;
    my %droppedSeqNames;
    print "Dropping from mir_seqinfo file:\n";
    for my $seqInfo (@{$mirSeqInfoLines}) {
	my($species,$seqName,$name,$id,$relativeLocation,$genomicLocation,$altGenomicLocation,$neighbors) = @{$seqInfo};
	my $newLocation = getGenomicFromRelativeLocation($relativeLocation,$altGenomicLocation);
	if (checkInconsistentMirsSeqInfo($id,$name,$newLocation)) {
	    print "Dropping: $species,$seqName,$name,$id,$relativeLocation,$genomicLocation,$altGenomicLocation\n";
	    $droppedSeqNames{$seqName} = 1;
	} else {
	    push(@newMirSeqInfoLines,[$species,$seqName,$name,$id,$relativeLocation,$genomicLocation,$altGenomicLocation,$neighbors]);
	}
    }
    print "\n\nDropping from mir_sequences gff file:\n";
    my @newGffLines;
    for my $gffLineInfo (@{$gffLines}) {
	my($type,$seqName,$start,$stop,$strand,$id,$name,$parentId,$line) = @{$gffLineInfo};
	if ($droppedSeqNames{$seqName}) {
	    print "Dropping: $line\n";
	} else {
	    push(@newGffLines,[$type,$seqName,$start,$stop,$strand,$id,$name,$parentId,$line]);
	}
    }
    print "\n\nDropping from mir_sequences fasta file:\n";
    my @newFastaEntries;
    for my $fastaEntryInfo (@{$fastaEntries}) {
	my($seqName,$seq) = @{$fastaEntryInfo};
	if ($droppedSeqNames{$seqName}) {
	    print "Dropping: $seqName\n";
	} else {
	    push(@newFastaEntries,[$seqName,$seq]);
	}	
    }
    return(\@newMirSeqInfoLines,\@newGffLines,\@newFastaEntries);
}

sub checkInconsistentMirsSeqInfo {
    my($id,$name,$location) = @_;
    ($id) = split("_",$id);
    my $idFound = 0;
    foreach my $mbId (keys %{$inconsistentMirs}) {
	if ($id eq $mbId) {
	    $idFound = 1;
	}
	my $nameFound = 0;
	my $locationFound = 0;
	foreach my $inconsistentMirInfo (@{$inconsistentMirs->{$mbId}}) {
	    my($mbName,$mbLocation) = @{$inconsistentMirInfo};
	    my($mbChrom,$mbStart,$mbStop,$mbStrand) = parseLocation($mbLocation);
	    #if ($mbName eq "MI0015015") {
#		print "$mbId\t$mbName\t$mbChrom\t$mbStart\t$mbStop\t$mbStrand\n";
#	    }
	    if ($name eq $mbName) {
		$nameFound = 1;
	    }
#	    if ($mbId eq "MI0002382") {
#		print "$mbId\t$mbName\t$location\n$mbLocation\n\n";
#	    }
	    if ($location eq $mbLocation) {
		$locationFound = 1;
	    }
#	    if ($mbId eq "MI0002382") {
#		print "$id\t$idFound\t$nameFound\t$locationFound\n\n";
#	    }
	}
	if ($idFound || $nameFound || $locationFound) {
	    if ($idFound && $nameFound && $locationFound) {
		return 1;
	    } elsif (!($nameFound)) {
		die "Error: the names in the inconsistentMirsFile and miRSeq info file did not match for $mbId\n";
	    } elsif (!($locationFound)) {
		die "Error: the locations in the inconsistentMirsFile and miRSeq info file did not match for $mbId\n";
	    } else {
		die "Error: information for the $mbId does not match up between inconsistentMirsFile and miRSeq info file for an unknown reason.";
	    }
	}
    }
    return 0;
}

sub loadInconsistentMirsFile {
    my($inconsistentMirsFile) = @_;
    my %inconsistentMirs;
    open(DMF,$inconsistentMirsFile) or die "failed to open $inconsistentMirsFile\n";
    my $FIRST = 1;
    while (<DMF>) {
	chomp;
	unless ($FIRST) {
	    my($name,$id,$mirbaseLocation) = split(/\t/);
	    ($id) = split(/\_/,$id);
	    push(@{$inconsistentMirs{$id}},[$name,$mirbaseLocation]);
	} else {
	    $FIRST = 0;
	}
    }
    close(DMF);
    return \%inconsistentMirs;
}

sub loadMirSeqFastaLines {
    my($mirSeqFasta) = @_;
    my %sequences;
    my @fastaEntries;
    open(MSFA,$mirSeqFasta) or die "could not open $mirSeqFasta\n";
    while (<MSFA>) {
	chomp;
	my $seqName = $_;
	$seqName =~ s/>//;
	my $seq = <MSFA>;
	chomp($seq);
	$sequences{$seqName} = $seq;
	push(@fastaEntries,[$seqName,$seq]);
    }
    close(MSFA);
    return(\%sequences,\@fastaEntries);
}

sub loadMirSeqInfoLines {
    my($mirSeqInfoFile) = @_;
    my @mirSeqInfoLines;
    open(MSIF,$mirSeqInfoFile) or die "could not open $mirSeqInfoFile\n";
    while (<MSIF>) {
	chomp;
	unless ( /^#/ ) {
	    my $line = $_;
	    my($species,$seqName,$name,$id,$relativeLocation,$genomicLocation,$altGenomicLocation,$neighborLine) = split(/\t/,$line);
	    #aae-mir-286a=aae-mir-2944b:-17..79;aae-mir-2944a=aae-mir-2944b:438..499;aae-mir-309b=aae-mir-
	    my @neighbors;
	    if ($neighborLine) {
		my(@neighborData) = split(";",$neighborLine);
		foreach my $neighborInfo (@neighborData) {
		    my($neighborName,$neighborLocation) = split("=",$neighborInfo);
		    push(@neighbors,[$neighborName,$neighborLocation]);
		}
	    }
	    push(@mirSeqInfoLines,[$species,$seqName,$name,$id,$relativeLocation,$genomicLocation,$altGenomicLocation,\@neighbors]);
	}
    }
    close(MSIF);
    return \@mirSeqInfoLines;
}

sub loadMirSeqGffLines {
    my($mirSeqGff) = @_;
    my $hairpinTotal = 0;
    my $productTotal = 0;
    my $linesTotal = 0;
    my @gffLines;
    my %hairpins;
    my %products;
    my $lastHPId = "";
    open(MBGFF3,$mirSeqGff) or die "could not open $mirSeqGff\n";
    while(<MBGFF3>) {
	unless(/^\#/) {
	    chomp;
	    my $line = $_;
	    my($chrom,$source,$type,$start,$stop,$score,$strand,$phase,$info) = split(/\t/,$line);
	    $linesTotal += 1;
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
		push(@{$hairpins{$chrom}},[$start,$stop,$strand,$id,$name,$line]);
		push(@gffLines,["HP",$chrom,$start,$stop,$strand,$id,$name,0,$line]);
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
			die "Error: product is derived from $parentId but last hairpin in gff was $lastHPId.\n";
			#$parentId = $lastHPId;
		    } else {
			die "Error: last hairpin id in gff was $lastHPId and this product derives from $parentId.\n";
		    }
		}
		push(@{$products{$parentId}},[$chrom,$start,$stop,$strand,$id,$name,$line]);
		push(@gffLines,["PROD",$chrom,$start,$stop,$strand,$id,$name,$parentId,$line]);
		$productTotal += 1;
	    }
	}
    }
    close(MBGFF3);
    return(\%hairpins,\%products,\@gffLines,$hairpinTotal,$productTotal,$linesTotal);
}

sub getGenomicFromRelativeLocation {
    my($relLocation,$genomicLocation) = @_;
    my($seqTag,$relStart,$relStop,$relStrand) = parseLocation($relLocation);
    my($gChrom,$gStart,$gStop,$gStrand) = parseLocation($genomicLocation);
    my $newLocation;
    if ($relStrand eq "-") {
	#we're not expecting negative strands from this dataset
	die "Error: function not written to handle negative relative strands";
    } else {
	if ($gStrand eq "+") {
	    $newLocation = "$gChrom:".($gStart + ($relStart - 1))."..".($gStart + ($relStop - 1)).":$gStrand";
	} elsif ($gStrand eq "-") {
	    $newLocation = "$gChrom:".($gStop - ($relStop - 1))."..".($gStop - ($relStart - 1)).":$gStrand";
	} else {
	    die "Error: strand for genomic coord mut be + or -\n";
	}
    }
    return $newLocation;
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
