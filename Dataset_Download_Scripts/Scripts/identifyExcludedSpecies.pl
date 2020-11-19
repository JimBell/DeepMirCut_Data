#!/usr/bin/perl -w
use strict;

my $usage = "USAGE:\n$0 <mirbase file list> <mir_dataset.txt>\n";

my $mirbaseFileList = $ARGV[0] or die $usage;
my $mirDatasetFile = $ARGV[1] or die $usage;

my $speciesHash = getOrgsFromMirDatasetFile($mirDatasetFile);
my $speciesList = getOrgsFromMirbaseFileList($mirbaseFileList);

foreach my $species (@{$speciesList}) {
    unless ($speciesHash->{$species}) {
	print "$species\n";
    }
}

sub getOrgsFromMirDatasetFile {
    my($mirDatasetFile) = @_;
    my %speciesHash;
    open(MDF,$mirDatasetFile) or die "failed to open $mirDatasetFile";
    while (<MDF>) {
	chomp;
	unless ( /^#/ ) {
	    my($species) = split(/\t/);
	    $speciesHash{$species} = 1;
	}
    }
    close(MDF);
    return \%speciesHash;
}

sub getOrgsFromMirbaseFileList {
    my($mirbaseFileList) = @_;
    my @speciesList;
    open(MBFL,$mirbaseFileList) or die "failed to open $mirbaseFileList\n";
    while (<MBFL>) {
	chomp;
	unless ( /^#/ ) {
	    my($file) = $_;
	    chomp($file);
	    my($filePrefix) = split(/\./);
	    push(@speciesList,$filePrefix);
	}
    }
    close(MBFL);
    return \@speciesList;
}
