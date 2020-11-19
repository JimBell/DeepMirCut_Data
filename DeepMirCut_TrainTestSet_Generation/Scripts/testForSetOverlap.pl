#!/usr/bin/perl -w
use strict;

my $usage = "USAGE:\n$0 <dataset1> <dataset2>\n";

my $dataSetFile1 = $ARGV[0] or die $usage;
my $dataSetFile2 = $ARGV[1] or die $usage;

my $dataSetNameHash1 = readNames($dataSetFile1);
my $dataSetNameHash2 = readNames($dataSetFile2);

foreach my $name (keys %{$dataSetNameHash1}) {
    if ($dataSetNameHash2->{$name}) {
	print "Warning: $name in both datasets\n";
    }
}

sub readNames {
    my($dataSetFile) = @_;
    my %dataSetNameHash;
    open(DSF,$dataSetFile) or die "failed to open $dataSetFile\n";
    while (<DSF>) {
	chomp;
	my($id,$name) = split(/\t/);
	$dataSetNameHash{$name} = 1;
    }
    close(DSF);
    return \%dataSetNameHash;
}

