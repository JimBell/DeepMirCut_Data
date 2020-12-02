#!/usr/bin/perl -w
use miRBaseDataGatherer;
use strict;


my $mirBaseVersion = "22.1";
my $genomeId = "GCA_000004015.1";
my $speciesPrefix = "aae";
my $miRBaseGff = "$speciesPrefix\.gff3";
my $genomeFileLocation = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/004/015/GCA_000004015.1_Aedes_aegypti/GCA_000004015.1_Aedes_aegypti_genomic.fna.gz";
my $assemblyInfoFileLocation = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/004/015/GCA_000004015.1_Aedes_aegypti/GCA_000004015.1_Aedes_aegypti_assembly_report.txt";
my $gffFileLocation = "ftp://mirbase.org/pub/mirbase/$mirBaseVersion/genomes/$miRBaseGff";
my $precursorFastaLocation = "hairpin.fa";
my $matureFastaLocation = "mature.fa";

my $gDir = "./genomes";
my $outputDir = "./miRSequenceOutput";

unless ($genomeFileLocation =~ /\Q$genomeId\E/) {
    print "Warning: Double check that $genomeFileLocation is for $genomeId\n";
}

unless ($assemblyInfoFileLocation =~ /\Q$genomeId\E/) {
    print "Warning: Double check that $assemblyInfoFileLocation is for $genomeId\n";
}

#prepareing input Directory
prepareDirectory($gDir);

#prepareing output Directory
prepareDirectory($outputDir);

#downloading Genome
my $fastaFile = prepareFasta($genomeFileLocation,$gDir);
my $faiFile = $fastaFile ."\.fai";
unless (-e $faiFile) {
    indexFasta($fastaFile)
}

#retrieving chrom aliases
my($chromNameCol,$chromAliasCol) = (0,4);
my $chromAliases = getChromAliases($assemblyInfoFileLocation,$gDir,$chromNameCol,$chromAliasCol);

#retrieving sequences from mirbase hairpin fasta for testing purposes
my $miRBaseHPFasta = prepareFasta($precursorFastaLocation,".");
my $testSeqs = loadMiRBaseFasta($miRBaseHPFasta);

#retrieving sequences from mirbase hairpin fasta for testing purposes
my $miRBaseProductFasta = prepareFasta($matureFastaLocation,".");
my $prodTestSeqs = loadMiRBaseFasta($miRBaseProductFasta);

#retrieving chrom sizes
my $chromSizes = getChromSizesFromFAI($faiFile);

#retrieving hairpins and products from gff
my $gffFile = prepareGff($gDir,$gffFileLocation);
my($hairpins,$products,$expectedHPCount,$expectedProdCount) = miRBaseDataGatherer::readMirbaseGff3_withTests($gffFile);

#generate data for deepMirCut
my $newParams = {
    "miRBaseVersion" => $mirBaseVersion,
    "genomeId" => $genomeId,
    "speciesId" => $speciesPrefix,
    "genomeFileLocation" => $genomeFileLocation,
    "assemblyInfoFileLocation" => $assemblyInfoFileLocation,
    "gffFileLocation" => $gffFileLocation,
    "testFastaLocation" => $precursorFastaLocation,
    "processingDir" => $gDir,
    "expectedHPCount" => $expectedHPCount,
    "expectedProdCount" => $expectedProdCount,
    "outputDir" => $outputDir
};
my $parameters = miRBaseDataGatherer::getDefaultParameters($speciesPrefix,$outputDir);
$parameters = miRBaseDataGatherer::setParameters($parameters,$newParams);
miRBaseDataGatherer::generateMiRDataFiles($fastaFile,$hairpins,$products,$chromSizes,$chromAliases,$testSeqs,$prodTestSeqs,$parameters);

sub getChromAliases {
    my($assemblyInfoFileLocation,$gDir,$chromNameCol,$chromAliasCol) = @_;
    my $assemblyInfoFile = $assemblyInfoFileLocation;
    if ((uc($assemblyInfoFileLocation) =~ /^HTTP:/) || (uc($assemblyInfoFileLocation) =~ /^FTP:/)) {
	($assemblyInfoFile) = $assemblyInfoFileLocation =~ /([^\/]*$)/;
	$assemblyInfoFile = "$gDir/$assemblyInfoFile";
	unless ((-e $assemblyInfoFile)) {
	    my $cmd = "wget $assemblyInfoFileLocation --directory-prefix=$gDir";
	    my $out = `$cmd`;
	    print $out;
	}
    }
    my %chromAliases;
    open(AFLE,$assemblyInfoFile) or die "failed to open $assemblyInfoFile\n";
    while (<AFLE>) {
	chomp;
	unless ( /^\#/ ) {
	    my(@rowData) = split(/\t/);
	    unless ($chromAliases{$rowData[$chromNameCol]}) {
		$chromAliases{$rowData[$chromNameCol]} = $rowData[$chromAliasCol];
	    } else {
		print "Warning: Cannot key more than one chrom alias by " . $rowData[$chromNameCol] . "\n";
	    }
	}
    }
    close(AFLE);
    return \%chromAliases;
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

sub getChromSizesFromFAI {
    my($faiFile) = @_;
    my %chromSizes;
    open(FAI,$faiFile) or die "failed to open $faiFile";
    while (<FAI>) {
	chomp;
	my($chrom,$chromSize) = split(/\t/);
	$chromSizes{$chrom} = $chromSize;
    }
    close(FAI);
    return \%chromSizes;
}

sub indexFasta {
    my $cmd = "samtools faidx $fastaFile";
    my $out = `$cmd`;
    print $out;
}

sub prepareGff {
    my($gDir,$gffFileLocation) = @_;
    my $gffFile = $gffFileLocation;
    if ((uc($gffFileLocation) =~ /^HTTP:/) || (uc($gffFileLocation) =~ /^FTP:/)) {
	($gffFile) = $gffFileLocation =~ /([^\/]*$)/;
	$gffFile = "$gDir/$gffFile";
	unless ((-e $gffFile)) {
	    my $cmd = "wget $gffFileLocation --directory-prefix=$gDir";
	    my $out = `$cmd`;
	    print $out;
	}
    }
    return $gffFile;
}

sub prepareFasta {
    my($genomeFileLocation,$gDir) = @_;
    my $fastaFile = $genomeFileLocation;
    if ((uc($genomeFileLocation) =~ /^HTTP:/) || (uc($genomeFileLocation) =~ /^FTP:/)) {
	my($fileName) = $genomeFileLocation =~ /([^\/]*\.f.*$)/;
	if ($fileName =~ /\.gz$/) {
	    ($fastaFile) = $fileName =~ /(^.*?)\.gz/;
	    $fastaFile = "$gDir/$fastaFile";
	}
	my $newGenomeFileLocation = "$gDir/$fileName";
	unless ((-e $newGenomeFileLocation) || (-e $fastaFile)) {
	    my $cmd = "wget $genomeFileLocation --directory-prefix=$gDir";
	    my $out = `$cmd`;
	    print $out;
	}
	if ($fileName =~ /\.gz$/ && !(-e $fastaFile)) {
	    my $cmd = "gunzip $newGenomeFileLocation";
	    my $out = `$cmd`;
	    print $out;
	}
    }
    return $fastaFile;
}

sub prepareDirectory {
    my($directory) = @_;
    unless ($directory eq "" || $directory eq "." || $directory eq "./") {
	unless (-d $directory) {
	    unless (mkdir($directory)) {
		die "failed to make $directory";
	    }
	}
    }
}
