#!/usr/bin/perl -w
use strict;

my $usage = "USAGE:\n$0 <sequence directory>";

my $miRSeqDir = $ARGV[0] or die $usage;

my $runInfoPostfix = "\_runInfo.txt";
my $seqInfoPostfix = "\_seqInfo.txt";
my $seqFaPostfix = "\_sequences.fa";
my $seqGffPostfix = "\_sequences.gff";
my $warningsPostfix = "\_warnings.txt";

my $runInfoOutputFile = "mir$runInfoPostfix";
my $seqInfoOutputFile = "mir$seqInfoPostfix";
my $seqFaOutputFile = "mir$seqFaPostfix";
my $seqGffOutputFile = "mir$seqGffPostfix";
my $warningsOutputFile = "mir$warningsPostfix";

my @runInfoFiles = `ls $miRSeqDir/*$runInfoPostfix`;
my @seqInfoFiles = `ls $miRSeqDir/*$seqInfoPostfix`;
my @seqFaFiles = `ls $miRSeqDir/*$seqFaPostfix`;
my @seqGffFiles = `ls $miRSeqDir/*$seqGffPostfix`;
my @warningsFiles = `ls $miRSeqDir/*$warningsPostfix`;

print `rm $runInfoOutputFile`;
print `rm $seqInfoOutputFile`;
print `rm $seqFaOutputFile`;
print `rm $seqGffOutputFile`;
print `rm $warningsOutputFile`;
print `echo "##gff-version 3" > $seqGffOutputFile`;

foreach my $file (@runInfoFiles) {
    my $outputFile = $runInfoOutputFile;
    chomp($file);
    open(OPTF,">>$outputFile") or die "failed to open $runInfoOutputFile for appending";
    open(FLE,$file) or die "failed to open $file";
    print "appending $file to $outputFile\n";
    while (<FLE>) {
	chomp;
	print OPTF $_ . "\n";
    }
    close(FLE);
    close(OPTF);
};

foreach my $file (@seqInfoFiles) {
    my $outputFile = $seqInfoOutputFile;
    chomp($file);
    open(OPTF,">>$outputFile") or die "failed to open $runInfoOutputFile for appending";
    open(FLE,$file) or die "failed to open $file";
    print "appending $file to $outputFile\n";
    while (<FLE>) {
	chomp;
	print OPTF $_ . "\n";
    }
    close(FLE);
    close(OPTF);
};

foreach my $file (@seqFaFiles) {
    my $outputFile = $seqFaOutputFile;
    chomp($file);
    open(OPTF,">>$outputFile") or die "failed to open $runInfoOutputFile for appending";
    open(FLE,$file) or die "failed to open $file";
    print "appending $file to $outputFile\n";
    while (<FLE>) {
	chomp;
	print OPTF $_ . "\n";
    }
    close(FLE);
    close(OPTF);
};

foreach my $file (@seqGffFiles) {
    my $outputFile = $seqGffOutputFile;
    chomp($file);
    open(OPTF,">>$outputFile") or die "failed to open $runInfoOutputFile for appending";
    open(FLE,$file) or die "failed to open $file";
    print "appending $file to $outputFile\n";
    while (<FLE>) {
	chomp;
	unless ( /^#/ ) {
	    print OPTF $_ . "\n";
	}
    }
    close(FLE);
    close(OPTF);
};

foreach my $file (@warningsFiles) {
    my $outputFile = $warningsOutputFile;
    chomp($file);
    open(OPTF,">>$outputFile") or die "failed to open $runInfoOutputFile for appending";
    open(FLE,$file) or die "failed to open $file";
    print "appending $file to $outputFile\n";
    my($filePrefix) = $file =~ /([^\/]*)\Q$warningsPostfix\E/;
    print OPTF "\n" . $filePrefix . ":\n";
    while (<FLE>) {
	chomp;
	print OPTF $_ . "\n";
    }
    close(FLE);
    close(OPTF);
};


