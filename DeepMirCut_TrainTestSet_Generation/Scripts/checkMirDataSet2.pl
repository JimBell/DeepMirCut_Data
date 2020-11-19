#!/usr/bin/perl -w
use strict;

my $usage = "USAGE:\n$0 <mir dataset file> <hairpin.fa> <mirbase mature.fa>";

my $miRDatasetFile = $ARGV[0] or die $usage;
my $hairpinsFasta = $ARGV[1] or die $usage;
my $matureFasta = $ARGV[2] or die $usage;

my $hairpinSeqs = loadMiRBaseFasta($hairpinsFasta);
my $matureSeqs = loadMiRBaseFasta($matureFasta);
my $miRDataset = loadDatasetFile($miRDatasetFile);

for my $miRInfo (@{$miRDataset}) {
    my($csId,$name,$id,$product5p,$product3p,$drosha5p,$dicer5p,$dicer3p,$drosha3p,$hpStart,$hpStop,$sequence) = @{$miRInfo};
    my $hpSeq = substr($sequence,$hpStart-1,$hpStop-$hpStart+1);
    if ($hairpinSeqs->{$name}) {
	unless ($hairpinSeqs->{$name} eq $hpSeq) {
	    print "Warning: $name  Sequences are not the same\n$hairpinSeqs->{$name}\n$hpSeq\n\n"
	}
    } else {
	print "Warning:$name not found in $hairpinsFasta\n";
    }
    if ($drosha5p ne '-' && $dicer5p ne '-') {
	my($drosha5pStart,$drosha5pStop) = split(",",$drosha5p);
	my($dicer5pStart,$dicer5pStop) = split(",",$dicer5p);
	my $miRSeq = substr($sequence,$drosha5pStop-1,$dicer5pStart-$drosha5pStop+1);
	if ($matureSeqs->{$product5p}) {
	    unless ($matureSeqs->{$product5p} eq $miRSeq) {
		print "Warning: $product5p ($name) Sequences are not the same\n$matureSeqs->{$product5p}\n$miRSeq\n\n";
	    }
	} else {
	    print "Warning: $product5p ($name) not found in $matureFasta\n";
	}
    } elsif ($drosha5p ne '-' ||  $dicer5p ne '-') {
	die "Error: $product5p cut in one spot but not the other. $drosha5p and $dicer5p";
    }

    if ($dicer3p ne '-' && $drosha3p ne '-') {
	my($dicer3pStart,$dicer3pStop) = split(",",$dicer3p);
	my($drosha3pStart,$drosha3pStop) = split(",",$drosha3p);
	my $miRSeq = substr($sequence,$dicer3pStop-1,$drosha3pStart-$dicer3pStop+1);
	if ($matureSeqs->{$product3p}) {
	    unless ($matureSeqs->{$product3p} eq $miRSeq) {
		print "Warning: $product3p ($name) Sequences are not the same\n$matureSeqs->{$product3p}\n$miRSeq\n\n";
	    }
	} else {
	    print "Warning: $product3p ($name) not found in $matureFasta\n";
	}
    } elsif ($dicer3p ne '-' ||  $drosha3p ne '-') {
	die "Error: $product3p cut in one spot but not the other. $dicer3p and $drosha3p";
    }

    
}


sub loadDatasetFile {
    my($miRDataSetFile) = @_;
    my @miRDataSet;
    open(MDSF,$miRDataSetFile) or die "failed to open $miRDataSetFile\n";
    while (<MDSF>) {
	chomp;
	unless ( /^#/ ) {
	    my($csId,$name,$id,$product5p,$product3p,$drosha5p,$dicer5p,$dicer3p,$drosha3p,$hpStart,$hpStop,$sequence) = split(/\t/);
	    push(@miRDataSet,[$csId,$name,$id,$product5p,$product3p,$drosha5p,$dicer5p,$dicer3p,$drosha3p,$hpStart,$hpStop,$sequence]);
	}
    }
    close(MDSF);
    return \@miRDataSet;
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
			die "Error: product is derived from $parentId but last hairpin in gff was $lastHPId.\n";
			#$parentId = $lastHPId;
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
