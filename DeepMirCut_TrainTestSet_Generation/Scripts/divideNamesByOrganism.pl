#!/usr/bin/perl -w
use strict;

my $usage = "USAGE:\n <mirdataset file> <names file> <organisms file>\n";

my $mirDatasetFile = $ARGV[0] or die $usage;
my $namesFile = $ARGV[1] or die $usage;
my $organismFile = $ARGV[2] or die $usage;

my $outputPrefix = "names";

my $mirSpecies = getSpeciesPrefixesFromMiRDataset($mirDatasetFile);
my $organismLineage = readOrganismFile($organismFile);
my $names = readNamesFile($namesFile);
my $splitNames = splitNamesByOrg($names,$mirSpecies,$organismLineage);
printOutput($splitNames,$outputPrefix);

sub printOutput {
    my($splitNames,$outputPrefix) = @_;
    foreach my $orgType (keys %{$splitNames}) {
	my $outputFile = $outputPrefix . "\_$orgType.txt";
	open(OPTF,">$outputFile") or die "failed to open $outputFile for writing\n";
	foreach my $name (@{$splitNames->{$orgType}}) {
	    print OPTF "$name\n";
	}
	close(OPTF);
    }
}

sub splitNamesByOrg {
    my($names,$mirSpecies,$organismLineage) = @_;
    my %splitNames;
    foreach my $name (@{$names}) {
	my $speciesPrefix = $mirSpecies->{$name};
	my($orgType) = @{$organismLineage->{$speciesPrefix}};
	push(@{$splitNames{$orgType}},$name);
    }
    return \%splitNames;
}

sub getSpeciesPrefixesFromMiRDataset {
    my($mirDatasetFile) = @_;
    my %mirSpecies;
    open(MDF,$mirDatasetFile) or die "failed to open $mirDatasetFile\n";
    while (<MDF>) {
	chomp;
	unless ( /^\#/ ) {
	    my($speciesPrefix,$name) = split(/\t/);
	    $mirSpecies{$name} = $speciesPrefix;
	}
    }
    close(MDF);
    return \%mirSpecies;
}

sub readNamesFile {
    my($namesFile) = @_;
    my @names;
    open(NMF,$namesFile) or die "failed to open $namesFile\n";
    while (<NMF>) {
	chomp;
	unless ( /^\#/ ) {
	    push(@names,$_);
	}
    }
    close(NMF);
    return \@names;
}

sub readOrganismFile {
    my($organismFile) = @_;
    my %organismLineage;
    open(ORGF,$organismFile) or die "failed to open $organismFile\n";
    while (<ORGF>) {
	chomp;
	unless ( /^\#/ ) {
	    my($organism,$division,$name,$tree) = split(/\t/);
	    my @lineage = split(';',$tree); 
	    @{$organismLineage{$organism}} = @lineage;
	}
    }
    close(ORGF);
    return \%organismLineage;
}
