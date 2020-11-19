#!/usr/bin/perl -w
use strict;

my $usage = "USAGE:\n$0 <2 product fasta> <mir dataset> <organisms file>";

my $twoProductFasta = $ARGV[0] or die $usage;
my $mirDatasetFile = $ARGV[1] or die $usage;
my $organismFile = $ARGV[2] or die $usage;


my $twoProductNames = getNamesFromFasta($twoProductFasta);
my $mirSpecies = getSpeciesPrefixesFromMiRDataset($mirDatasetFile);
my $organismLineage = readOrganismFile($organismFile);

my $count = 0;
foreach my $name (@{$twoProductNames}) {
    my $speciesPrefix = $mirSpecies->{$name};
    my($orgType) = @{$organismLineage->{$speciesPrefix}};
    if ($orgType eq "Metazoa") {
	print $orgType . " " . $name . "\n";
	$count++;
    }
}
print "Metazoan = $count\n";


sub getNamesFromFasta {
    my($twoProductFasta) = @_;
    my @twoProductNames;
    open(TPF,$twoProductFasta) or die "failed to open $twoProductFasta\n";
    while (<TPF>) {
	chomp;
	if ( /^>/ ) {
	    s/>//;
	    my($name) = $_;
	    push(@twoProductNames,$name);
	}
    }
    close(TPF);
    return \@twoProductNames;
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
