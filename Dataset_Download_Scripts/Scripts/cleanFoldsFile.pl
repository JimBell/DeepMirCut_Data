#!/usr/bin/perl -w
use strict;

my $usage = "USAGE:\n$0 <folds file>\n";

my $foldsFile = $ARGV[0] or die $usage;

my $newLines = "";
open(FF,$foldsFile) or die "failed to open $foldsFile\n";
while (<FF>) {
    chomp;
    my $line = $_;
    if ( /^>/ ) {
	($line) = split;
    }
    $newLines .= $line . "\n";
}
close(FF);

open(NFF,">$foldsFile") or die "fialed to open $foldsFile for writing\n";
print NFF $newLines;
close(NFF);
