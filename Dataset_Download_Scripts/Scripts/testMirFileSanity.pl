#!/usr/bin/perl -w
use miRBaseDataGatherer;
use strict;

my $usage = "USAGE:\n <mir_sequences.fa> <mir_sequences.gff> <mir_seqInfo.txt> <mirbase hairpin.fa> <mirbase mature.fa>\n";

my $mirSeqFasta = $ARGV[0] or die $usage;
my $mirSeqGff = $ARGV[1] or die $usage;
my $mirSeqInfoFile = $ARGV[2] or die $usage;
my $mirbaseHairpinFasta = $ARGV[3] or die $usage;
my $mirbaseMatureFasta = $ARGV[4] or die $usage;

my $mirbaseHairpins = loadMiRBaseFasta($mirbaseHairpinFasta);
my $mirbaseProducts = loadMiRBaseFasta($mirbaseMatureFasta);
my $mirSeqs = loadMiRBaseFasta($mirSeqFasta);

testGffSanity($mirSeqGff);
testGffSequenceSanity($mirSeqGff,$mirSeqs,$mirbaseHairpins,$mirbaseProducts);
testSeqInfoFileSanity($mirSeqInfoFile,$mirSeqs,$mirbaseHairpins);

sub testSeqInfoFileSanity {
    my($mirSeqInfoFile,$mirSeqs,$mirbaseHairpins) = @_;
    open(SIFLE,$mirSeqInfoFile) or die "faield to open $mirSeqInfoFile\n";
    while (<SIFLE>) {
	chomp;
	my $line = $_;
	unless ( $line =~ /^\#/ ) {
	    my($speciesPrefix,$seqTag,$hpName,$hpId,$relLocation,$genomicLoc,$altGenomicLoc,$neighborEntry) = split(/\t/,$line);
	    my($relChrom,$relStart,$relStop) = miRBaseDataGatherer::parseLocation($relLocation);
	    $relStart--; #converting to zero-based
	    $relStop--;  #converting to zero-based
	    unless ($mirSeqs->{$seqTag}) {
		die "$seqTag not found in mir sequences fasta file\n"
	    }
	    my $seq = substr($mirSeqs->{$seqTag},$relStart,$relStop - $relStart + 1);
	    if ($mirbaseHairpins->{$hpName}) {
		unless ($seq eq $mirbaseHairpins->{$hpName}) {
		    die "sequence in mirbase hairpins.fa not equal to sequence in mir_sequences.fa\n$seq vs $mirbaseHairpins->{$hpName}\n";
		}
	    } else {
		print "Warning $hpName not found in mirbase hairpins.fa\n";
	    }
	    #checking neighbors
	    my $mirSeqLength = length($mirSeqs->{$seqTag});
	    if (($neighborEntry) && $neighborEntry ne "") {
		my @neighbors = split(';',$neighborEntry);
		foreach my $neighborInfo (@neighbors) {
		    my($neighborName,$relNeighborLocation) = split("=",$neighborInfo);
		    my($relNeighborChrom,$relNeighborStart,$relNeighborStop) = miRBaseDataGatherer::parseLocation($relNeighborLocation);
		    my $newRelNeighborStart = $relNeighborStart;
		    my $newRelNeighborStop = $relNeighborStop;
		    my $startCorrection = 0;
		    my $stopCorrection = 0;
		    if ($newRelNeighborStart < 1) {
			$newRelNeighborStart = 1;
			$startCorrection = $newRelNeighborStart - $relNeighborStart;
		    }
		    if ($newRelNeighborStop > $mirSeqLength) {
			$newRelNeighborStop = $mirSeqLength;
			$stopCorrection = $newRelNeighborStop - $relNeighborStop;
		    }
		    my $newRelNeighborLocation = "$relNeighborChrom:$newRelNeighborStart..$newRelNeighborStop";
		    $newRelNeighborStart--;
		    $newRelNeighborStop--;
		    my $neighborSeq = substr($mirSeqs->{$seqTag},$newRelNeighborStart,$newRelNeighborStop-$newRelNeighborStart+1);
		    my $testNeighborSeq = substr($mirbaseHairpins->{$neighborName}, $startCorrection, length($mirbaseHairpins->{$neighborName}) - $startCorrection + $stopCorrection);
		    unless ($neighborSeq eq $testNeighborSeq) {
			die "Error: neighboring seq from $hpName called $neighborName not equal to what is expected:\n$neighborSeq\n$testNeighborSeq\n";
		    }
		}
	    }
	}
    }
    close(SIFLE);
}

sub testGffSequenceSanity {
    my($mirSeqGff,$mirSeqs,$mirbaseHairpins,$mirbaseProducts) = @_;
    my($hairpins,$products) = miRBaseDataGatherer::readMirbaseGff3($mirSeqGff);
    my($tagCount,$hpCount) = (0,0);
    foreach my $tag (keys %{$hairpins}) {
	$tagCount++;
	foreach my $hairpinInfo (@{$hairpins->{$tag}}) {
	    $hpCount++;
	    my($hpStart,$hpStop,$hpStrand,$hpId,$hpName) = @{$hairpinInfo};
	    $hpStart--; #convert to zero-based
	    $hpStop--;  #convert to zero-based
	    unless ($mirSeqs->{$tag}) {
		die "$tag not found in mir sequences fasta file\n"
	    }
	    my $seq = substr($mirSeqs->{$tag},$hpStart,$hpStop - $hpStart + 1);
	    if ($mirbaseHairpins->{$hpName}) {
		unless ($seq eq $mirbaseHairpins->{$hpName}) {
		    die "sequence in mirbase hairpins.fa not equal to sequence in mir_sequences.fa\n$seq vs $mirbaseHairpins->{$hpName}\n";
		}
	    } else {
		print "Warning $hpName not found in mirbase hairpins.fa\n";
	    }
	    foreach my $productInfo (@{$products->{$hpId}}) {
		my($prodChrom,$prodStart,$prodStop,$prodStrand,$prodId,$prodName) = @{$productInfo};
		$prodStart--; #convert to zero-based
		$prodStop--;  #convert to zero-based
		my $prodSeq = substr($mirSeqs->{$tag},$prodStart,$prodStop - $prodStart + 1);
		if ($mirbaseProducts->{$prodName}) {
		    unless ($prodSeq eq $mirbaseProducts->{$prodName}) {
			die "sequence in mirbase hairpins.fa not equal to sequence in mir_sequences.fa\n$prodSeq vs $mirbaseProducts->{$prodName}\n";
		    }
		} else {
		    print "Warning $prodName not found in mirbase hairpins.fa\n";
		}	
	    }
	}
    }
    unless ($tagCount == $hpCount) {
	die "Error: tagCount and hpCount expected to be equal.\ntagCount = $tagCount\thpCount = $hpCount\n";
    }
}

sub testGffSanity {
    my($mirSeqGff) = @_;
    my($hairpins,$products,$hairpinTotal,$productTotal) = miRBaseDataGatherer::readMirbaseGff3_withTests($mirSeqGff);
    my($hairpins2,$products2) = miRBaseDataGatherer::readMirbaseGff3($mirSeqGff);
    my($hairpinTotal2,$productTotal2) = (0,0);
    foreach my $chrom (keys %{$hairpins}) {
	unless ($hairpins2->{$chrom}) {
	    die "$chrom in readMirbaseGff3_withTests output but not in readMirbaseGff3 output\n";
	}
	foreach my $hairpinInfo (@{$hairpins->{$chrom}}) {
	    my($hpStart,$hpStop,$hpStrand,$hpId,$hpName) = @{$hairpinInfo};
	    my $hpFound = 0;
	    foreach my $hairpinInfo2 (@{$hairpins2->{$chrom}}) {
		my($hpStart2,$hpStop2,$hpStrand2,$hpId2,$hpName2) = @{$hairpinInfo2};
		if ($hpStart == $hpStart2 && $hpStop == $hpStop2 && $hpStrand eq $hpStrand2 && $hpId eq $hpId2 && $hpName eq $hpName2) {
		    $hpFound = 1;
		}
	    }
	    unless ($hpFound) {
		die "$hpId - $hpName found in readMirbaseGff3_withTests output but not in readMirbaseGff3 output\n";
	    }
	}
    }
    foreach my $chrom (keys %{$hairpins2}) {
	unless ($hairpins->{$chrom}) {
	    die "$chrom in readMirbaseGff3 output but not in readMirbaseGff3_withTests output\n";
	}
	foreach my $hairpinInfo2 (@{$hairpins2->{$chrom}}) {
	    my($hpStart2,$hpStop2,$hpStrand2,$hpId2,$hpName2) = @{$hairpinInfo2};
	    my $hpFound = 0;
	    foreach my $hairpinInfo (@{$hairpins->{$chrom}}) {
		my($hpStart,$hpStop,$hpStrand,$hpId,$hpName) = @{$hairpinInfo};
		if ($hpStart == $hpStart2 && $hpStop == $hpStop2 && $hpStrand eq $hpStrand2 && $hpId eq $hpId2 && $hpName eq $hpName2) {
		    $hpFound = 1;
		}
	    }
	    unless ($hpFound) {
		die "$hpId2 - $hpName2 found in readMirbaseGff3 output but not in readMirbaseGff3_withTests output\n";
	    }
	    $hairpinTotal2++;
	}
    }
    foreach my $hpId (keys %{$products}) {
	unless ($products2->{$hpId}) {
	    die "products from readMirbaseGff3_withTests keyed by $hpId but not products from readMirbaseGff3\n";
	}
	foreach my $prodInfo (@{$products->{$hpId}}) {
	    my($prodChrom,$prodStart,$prodStop,$prodStrand,$prodId,$prodName) = @{$prodInfo};
	    my $prodFound = 0;
	    #print "$prodId\t$prodName\n";
	    foreach my $prodInfo2 (@{$products2->{$hpId}}) {
		my($prodChrom2,$prodStart2,$prodStop2,$prodStrand2,$prodId2,$prodName2) = @{$prodInfo2};
		if ($prodChrom eq $prodChrom2 && $prodStart == $prodStart2 && $prodStop == $prodStop2 && $prodId eq $prodId2 && $prodName eq $prodName2) {
		    $prodFound = 1;
		}
	    }
	    unless ($prodFound) {
		die "$prodId $prodName found in readMirbaseGff3_withTests output but not in readMirbaseGff3 output\n";
	    }
	}

    }
    foreach my $hpId2 (keys %{$products2}) {
	unless ($products->{$hpId2}) {
	    die "products from readMirbaseGff3 keyed by $hpId2 but not products from readMirbaseGff3_withTests\n";
	}
	foreach my $prodInfo2 (@{$products2->{$hpId2}}) {
	    my($prodChrom2,$prodStart2,$prodStop2,$prodStrand2,$prodId2,$prodName2) = @{$prodInfo2};
	    my $prodFound = 0;
	    foreach my $prodInfo (@{$products->{$hpId2}}) {
		my($prodChrom,$prodStart,$prodStop,$prodStrand,$prodId,$prodName) = @{$prodInfo};
		if ($prodChrom eq $prodChrom2 && $prodStart == $prodStart2 && $prodStop == $prodStop2 && $prodId eq $prodId2 && $prodName eq $prodName2) {
		    $prodFound = 1;
		}
	    }
	    unless ($prodFound) {
		die "$prodId2 - $prodName2 found in readMirbaseGff3 output but not in readMirbaseGff3_withTests output\n";
	    }
	    $productTotal2++;
	}
    }
    unless ($hairpinTotal == $hairpinTotal2 && $productTotal == $productTotal2) {
	die "error reading gff hairpins: $hairpinTotal vs $hairpinTotal2 ..... products: $productTotal vs $productTotal2\n";
    }
}

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

