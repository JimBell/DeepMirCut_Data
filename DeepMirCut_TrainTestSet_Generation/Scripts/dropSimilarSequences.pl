#!/usr/bin/perl
use List::Util;
use strict;

my $parameters = {
    multiples => 9,
    minBuffer => 30,
    maxBuffer => 50,
    maxSize => 250,
    maxSimilar => ~0,
    trainOutputFile => "trainSet.txt",
    validateOutputFile => "validateSet.txt",
    testOutputFile => "testSet.txt",
    counter => 0,
};

my $usage = "USAGE:\n$0 <trainSet> <validationSet> <testSet> <hairpins.fa> <max identity score>\n";

