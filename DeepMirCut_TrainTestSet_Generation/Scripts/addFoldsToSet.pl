#!/usr/bin/perl -w 
do "Scripts/bpRNA.pl";
use strict;

my $usage = "USAGE:\n$0 <data set>";

my $dataSetFile = $ARGV[0] or die $usage;

my $outputFile = $dataSetFile . "\_wFolds.txt";
if ( $dataSetFile =~ /.*\.txt$/ ) {
    ($outputFile) = $dataSetFile =~ /(.*)\.txt$/;
    $outputFile = $outputFile . "\_wFolds.txt";
}

my $fastaFile = $dataSetFile . ".fa";
my $foldFile = $dataSetFile . ".folds";

my $dataSet = readDataSet($dataSetFile);
createFastaFromDataset($dataSet,$fastaFile);

my $out = `RNAfold -i $fastaFile --noLP --noPS > $foldFile`;
print $out;
my($folds,$seqs) = readFoldFile($foldFile);
my $context = getContextFromBPRNA($folds,$seqs);
printOutputFile($folds,$context,$dataSetFile,$outputFile);

sub printOutputFile {
    my($folds,$context,$dataSetFile,$outputFile) = @_;
    open(OPTF,">$outputFile") or die "failed to open $outputFile for writing\n";
    open(DS,$dataSetFile) or die "failed to open $dataSetFile\n";
    while (<DS>) {
	chomp;
	my $line = $_;
	unless (/^\#/) {
	    my($id,$name,$mirbaseId,$product5p,$product3p,$drosha5p,$dicer5p,$dicer3p,$drosha3p,$hpStart,$hpStop,$sequence) = split(/\t/,$line);
	    print OPTF "$line\t$folds->{$id}\t$context->{$id}\n";
	}
    }
    close(DS);
    close(OPTF);
}

sub readFoldFile {
    my($foldFile) = @_;
    my %folds;
    my %seqs;
    open(FF,$foldFile) or die "failed to open $foldFile\n";
    while (<FF>) {
	chomp;
	$_ =~ s/>//;
	my $id = $_;
	my $sequence = <FF>;
	chomp($sequence);
	$seqs{$id} = $sequence;
	my $foldLine = <FF>;
	chomp($foldLine);
	my($fold,$mfe) = split(/\s+/,$foldLine);
	$folds{$id} = $fold;
    }
    close(FF);
    return(\%folds,\%seqs);
}

sub createFastaFromDataset {
    my($dataSet,$fastaFile) = @_;
    open(FA,">$fastaFile") or die "failed to open $fastaFile for writing\n";
    foreach my $id (keys %{$dataSet}) {
	my($name,$mirbaseId,$product5p,$product3p,$drosha5p,$dicer5p,$dicer3p,$drosha3p,$hpStart,$hpStop,$sequence) = @{$dataSet->{$id}};
	print FA ">$id\n$sequence\n";
    }
    close(FA);
}

sub readDataSet {
    my($dataSetFile) = @_;
    my %dataSet;
    open(DS,$dataSetFile) or die "failed to open $dataSetFile\n";
    while (<DS>) {
	chomp;
	unless (/^\#/) {
	    my($id,$name,$mirbaseId,$product5p,$product3p,$drosha5p,$dicer5p,$dicer3p,$drosha3p,$hpStart,$hpStop,$sequence) = split(/\t/);
	    $dataSet{$id} = [$name,$mirbaseId,$product5p,$product3p,$drosha5p,$dicer5p,$dicer3p,$drosha3p,$hpStart,$hpStop,$sequence];
	}
    }
    close(DS);
    return \%dataSet;
}

sub getContextFromBPRNA {
    my($folds,$seqs) = @_;
    my %context;
    foreach my $id (sort keys %{$folds}) {
	my $sequence = $seqs->{$id};
	my $fold = $folds->{$id};
	$context{$id} = dotBracketToStructureArray($sequence,$fold);
    }
    return \%context;
}
