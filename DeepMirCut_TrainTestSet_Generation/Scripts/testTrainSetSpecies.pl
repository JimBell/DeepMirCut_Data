#!/usr/bin/perl -w
use strict;

my $usage = "USAGE:\n$0 <train set> <organism list>\n";

my $trainSetFile = $ARGV[0] or die $usage;
my $organismFile = $ARGV[1] or die $usage;

my $organismLineage = readOrganismFile($organismFile);
my $names = getNamesFromTrainset($trainSetFile);

my %speciesSet;
foreach my $name (@{$names}) {
    my($speciesPrefix) = split("-",$name);
    $speciesSet{$speciesPrefix} = 1;
}

my %speciesGroup;
foreach my $speciesPrefix (keys %speciesSet) {
    $speciesGroup{${$organismLineage->{$speciesPrefix}}[0]} = 1;
}

foreach my $orgType (keys %speciesGroup) {
    print $orgType . "\n"
}

sub getNamesFromTrainset {
    my($trainSetFile) = @_;
    my @names;
    open(TSF,$trainSetFile) or die "failed to open $trainSetFile\n";
    while (<TSF>) {
	chomp;
	unless (/^#/) {
	    my($id,$name) = split(/\t/);
	    push(@names,$name);
	}
    }
    close(TSF);
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
