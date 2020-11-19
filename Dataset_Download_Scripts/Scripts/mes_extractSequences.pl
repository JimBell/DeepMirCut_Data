#!/usr/bin/perl -w
use miRBaseDataGatherer;
use strict;

my $mirBaseVersion = "22.1";
my $genomeId = "Cassava4";
my $speciesPrefix = "mes";
my $miRBaseGff = "$speciesPrefix\.gff3";
my $genomeFileLocation = "ftp://ftp.psb.ugent.be/pub/plaza/plaza_public_dicots_03//Genomes/mes.con.gz";
my $assemblyInfoFileLocation = "";
my $gffFileLocation = "ftp://mirbase.org/pub/mirbase/$mirBaseVersion/genomes/$miRBaseGff";
my $precursorFastaLocation = "/nfs0/BB/Hendrix_Lab/MicroRNAs/miRBase_v22/hairpin.fa";
my $matureFastaLocation = "/nfs0/BB/Hendrix_Lab/MicroRNAs/miRBase_v22/mature.fa";

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
my $fastaFile = prepareFasta_MES($genomeFileLocation,$gDir);
my $faiFile = $fastaFile ."\.fai";
unless (-e $faiFile) {
    indexFasta($fastaFile)
}

#retrieving chrom aliases
my($chromNameCol,$chromAliasCol) = (0,0);
my $chromAliases = getChromAliases($faiFile,$gDir,$chromNameCol,$chromAliasCol);

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
	    my $cmd2 = "sed -i \'s/\r//g\' $assemblyInfoFile";
	    my $out2 = `$cmd2`;
	    print $out2;
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

sub getChromAliases_replacePrefix {
    my($assemblyInfoFileLocation,$gDir,$chromNameCol,$chromAliasCol,$prefix,$replacementPrefix) = @_;
    my $assemblyInfoFile = $assemblyInfoFileLocation;
    if ((uc($assemblyInfoFileLocation) =~ /^HTTP:/) || (uc($assemblyInfoFileLocation) =~ /^FTP:/)) {
	($assemblyInfoFile) = $assemblyInfoFileLocation =~ /([^\/]*$)/;
	$assemblyInfoFile = "$gDir/$assemblyInfoFile";
	unless ((-e $assemblyInfoFile)) {
	    my $cmd = "wget $assemblyInfoFileLocation --directory-prefix=$gDir";
	    my $out = `$cmd`;
	    print $out;
	    my $cmd2 = "sed -i \'s/\r//g\' $assemblyInfoFile";
	    my $out2 = `$cmd2`;
	    print $out2;
	}
    }
    my %chromAliases;
    open(AFLE,$assemblyInfoFile) or die "failed to open $assemblyInfoFile\n";
    while (<AFLE>) {
	chomp;
	unless ( /^\#/ ) {
	    my(@rowData) = split(/\t/);
	    my $chromName = $rowData[$chromNameCol];
	    $chromName =~ s/\Q$prefix\E/\Q$replacementPrefix\E/;
	    unless ($chromAliases{$chromName}) {
		$chromAliases{$chromName} = $rowData[$chromAliasCol];
	    } else {
		print "Warning: Cannot key more than one chrom alias by " . $chromName . "\n";
	    }
	}
    }
    close(AFLE);
    return \%chromAliases;
}

sub getChromAliases_chromOnly_appendPrefix {
    my($assemblyInfoFileLocation,$gDir,$chromNameCol,$chromAliasCol,$prefix) = @_;
    my $assemblyInfoFile = $assemblyInfoFileLocation;
    if ((uc($assemblyInfoFileLocation) =~ /^HTTP:/) || (uc($assemblyInfoFileLocation) =~ /^FTP:/)) {
	($assemblyInfoFile) = $assemblyInfoFileLocation =~ /([^\/]*$)/;
	$assemblyInfoFile = "$gDir/$assemblyInfoFile";
	unless ((-e $assemblyInfoFile)) {
	    my $cmd = "wget $assemblyInfoFileLocation --directory-prefix=$gDir";
	    my $out = `$cmd`;
	    print $out;
	    my $cmd2 = "sed -i \'s/\r//g\' $assemblyInfoFile";
	    my $out2 = `$cmd2`;
	    print $out2;
	}
    }
    my %chromAliases;
    open(AFLE,$assemblyInfoFile) or die "failed to open $assemblyInfoFile\n";
    while (<AFLE>) {
	chomp;
	unless ( /^\#/ ) {
	    my(@rowData) = split(/\t/);
	    if ($rowData[3] eq "Chromosome") {
		unless ($chromAliases{"$prefix$rowData[$chromNameCol]"}) {
		    $chromAliases{"$prefix$rowData[$chromNameCol]"} = $rowData[$chromAliasCol];
		} else {
		    print "Warning: Cannot key more than one chrom alias by $prefix$rowData[$chromNameCol]\n";
		}
	    }
	}
    }
    close(AFLE);
    return \%chromAliases;
}

sub getChromAliases_withChromNamePrefix {
    my($assemblyInfoFileLocation,$gDir,$chromNameCol,$chromAliasCol,$chromNamePrefix) = @_;
    my $assemblyInfoFile = $assemblyInfoFileLocation;
    if ((uc($assemblyInfoFileLocation) =~ /^HTTP:/) || (uc($assemblyInfoFileLocation) =~ /^FTP:/)) {
	($assemblyInfoFile) = $assemblyInfoFileLocation =~ /([^\/]*$)/;
	$assemblyInfoFile = "$gDir/$assemblyInfoFile";
	unless ((-e $assemblyInfoFile)) {
	    my $cmd = "wget $assemblyInfoFileLocation --directory-prefix=$gDir";
	    my $out = `$cmd`;
	    print $out;
	    my $cmd2 = "sed -i \'s/\r//g\' $assemblyInfoFile";
	    my $out2 = `$cmd2`;
	    print $out2;
	}
    }
    my %chromAliases;
    open(AFLE,$assemblyInfoFile) or die "failed to open $assemblyInfoFile\n";
    while (<AFLE>) {
	chomp;
	unless ( /^\#/ ) {
	    my(@rowData) = split(/\t/);
	    unless ($chromAliases{"$chromNamePrefix$rowData[$chromNameCol]"}) {
		$chromAliases{"$chromNamePrefix$rowData[$chromNameCol]"} = $rowData[$chromAliasCol];
	    } else {
		print "Warning: Cannot key more than one chrom alias by " . $rowData[$chromNameCol] . "\n";
	    }
	}
    }
    close(AFLE);
    return \%chromAliases;
}

sub getChromAliases_dropChromNamePrefix {
    my($assemblyInfoFileLocation,$gDir,$chromNameCol,$chromAliasCol,$chromNamePrefix) = @_;
    my $assemblyInfoFile = $assemblyInfoFileLocation;
    if ((uc($assemblyInfoFileLocation) =~ /^HTTP:/) || (uc($assemblyInfoFileLocation) =~ /^FTP:/)) {
	($assemblyInfoFile) = $assemblyInfoFileLocation =~ /([^\/]*$)/;
	$assemblyInfoFile = "$gDir/$assemblyInfoFile";
	unless ((-e $assemblyInfoFile)) {
	    my $cmd = "wget $assemblyInfoFileLocation --directory-prefix=$gDir";
	    my $out = `$cmd`;
	    print $out;
	    my $cmd2 = "sed -i \'s/\r//g\' $assemblyInfoFile";
	    my $out2 = `$cmd2`;
	    print $out2;
	}
    }
    my %chromAliases;
    open(AFLE,$assemblyInfoFile) or die "failed to open $assemblyInfoFile\n";
    while (<AFLE>) {
	chomp;
	unless ( /^\#/ ) {
	    my(@rowData) = split(/\t/);
	    my($chromName) = $rowData[$chromNameCol] =~ /^\Q$chromNamePrefix\E(.*)/;
	    unless ($chromName) {
		die "Error: failed to drop $chromNamePrefix from $rowData[$chromNameCol]\n";
	    }
	    unless ($chromAliases{$chromName}) {
		$chromAliases{$chromName} = $rowData[$chromAliasCol];
	    } else {
		print "Warning: Cannot key more than one chrom alias by " . $chromName . "\n";
	    }
	}
    }
    close(AFLE);
    return \%chromAliases;
}

sub getChromAliases_withNormalChromNames {
    my($assemblyInfoFileLocation,$gDir,$chromNameCol,$chromAliasCol,$normalChromNameCol) = @_;
    my $assemblyInfoFile = $assemblyInfoFileLocation;
    if ((uc($assemblyInfoFileLocation) =~ /^HTTP:/) || (uc($assemblyInfoFileLocation) =~ /^FTP:/)) {
	($assemblyInfoFile) = $assemblyInfoFileLocation =~ /([^\/]*$)/;
	$assemblyInfoFile = "$gDir/$assemblyInfoFile";
	unless ((-e $assemblyInfoFile)) {
	    my $cmd = "wget $assemblyInfoFileLocation --directory-prefix=$gDir";
	    my $out = `$cmd`;
	    print $out;
	    my $cmd2 = "sed -i \'s/\r//g\' $assemblyInfoFile";
	    my $out2 = `$cmd2`;
	    print $out2;
	}
    }
    my %chromAliases;
    open(AFLE,$assemblyInfoFile) or die "failed to open $assemblyInfoFile\n";
    while (<AFLE>) {
	chomp;
	unless ( /^\#/ ) {
	    my(@rowData) = split(/\t/);
	    if ($rowData[3] eq "Chromosome") {
		#Assigned-molecule is labeled as a chromosome
		unless ($chromAliases{$rowData[$normalChromNameCol]}) {
		    $chromAliases{$rowData[$normalChromNameCol]} = $rowData[$chromAliasCol];
		} else {
		    print "Warning: Cannot key more than one chrom alias by " . $rowData[$normalChromNameCol] . "\n";
		}
	    } else {
		unless ($chromAliases{$rowData[$chromNameCol]}) {
		    $chromAliases{$rowData[$chromNameCol]} = $rowData[$chromAliasCol];
		} else {
		    print "Warning: Cannot key more than one chrom alias by " . $rowData[$chromNameCol] . "\n";
		}
	    }
	}
    }
    close(AFLE);
    return \%chromAliases;
}

sub getChromAliases_withNormalChromNames_pref {
    my($assemblyInfoFileLocation,$gDir,$chromNameCol,$chromAliasCol,$normalChromNameCol,$pref) = @_;
    my $assemblyInfoFile = $assemblyInfoFileLocation;
    if ((uc($assemblyInfoFileLocation) =~ /^HTTP:/) || (uc($assemblyInfoFileLocation) =~ /^FTP:/)) {
	($assemblyInfoFile) = $assemblyInfoFileLocation =~ /([^\/]*$)/;
	$assemblyInfoFile = "$gDir/$assemblyInfoFile";
	unless ((-e $assemblyInfoFile)) {
	    my $cmd = "wget $assemblyInfoFileLocation --directory-prefix=$gDir";
	    my $out = `$cmd`;
	    print $out;
	    my $cmd2 = "sed -i \'s/\r//g\' $assemblyInfoFile";
	    my $out2 = `$cmd2`;
	    print $out2;
	}
    }
    my %chromAliases;
    open(AFLE,$assemblyInfoFile) or die "failed to open $assemblyInfoFile\n";
    while (<AFLE>) {
	chomp;
	unless ( /^\#/ ) {
	    my(@rowData) = split(/\t/);
	    if ($rowData[3] eq "Chromosome" && ($rowData[$chromNameCol] =~ /^\Q$pref\E/)) {
		#Assigned-molecule is labeled as a chromosome
		unless ($chromAliases{$rowData[$normalChromNameCol]}) {
		    $chromAliases{$rowData[$normalChromNameCol]} = $rowData[$chromAliasCol];
		} else {
		    print "Warning: Cannot key more than one chrom alias by " . $rowData[$normalChromNameCol] . "\n";
		}
	    } else {
		unless ($chromAliases{$rowData[$chromNameCol]}) {
		    $chromAliases{$rowData[$chromNameCol]} = $rowData[$chromAliasCol];
		} else {
		    print "Warning: Cannot key more than one chrom alias by " . $rowData[$chromNameCol] . "\n";
		}
	    }
	}
    }
    close(AFLE);
    return \%chromAliases;
}

sub getChromAliases_withNormalChromNames_withNormalChromPrefix {
    my($assemblyInfoFileLocation,$gDir,$chromNameCol,$chromAliasCol,$normalChromNameCol,$normalChromNamePrefix) = @_;
    my $assemblyInfoFile = $assemblyInfoFileLocation;
    if ((uc($assemblyInfoFileLocation) =~ /^HTTP:/) || (uc($assemblyInfoFileLocation) =~ /^FTP:/)) {
	($assemblyInfoFile) = $assemblyInfoFileLocation =~ /([^\/]*$)/;
	$assemblyInfoFile = "$gDir/$assemblyInfoFile";
	unless ((-e $assemblyInfoFile)) {
	    my $cmd = "wget $assemblyInfoFileLocation --directory-prefix=$gDir";
	    my $out = `$cmd`;
	    print $out;
	    my $cmd2 = "sed -i \'s/\r//g\' $assemblyInfoFile";
	    my $out2 = `$cmd2`;
	    print $out2;
	}
    }
    my %chromAliases;
    open(AFLE,$assemblyInfoFile) or die "failed to open $assemblyInfoFile\n";
    while (<AFLE>) {
	chomp;
	unless ( /^\#/ ) {
	    my(@rowData) = split(/\t/);
	    if ($rowData[3] eq "Chromosome") {
		#Assigned-molecule is labeled as a chromosome
		unless ($chromAliases{"$normalChromNamePrefix$rowData[$normalChromNameCol]"}) {
		    $chromAliases{"$normalChromNamePrefix$rowData[$normalChromNameCol]"} = $rowData[$chromAliasCol];
		} else {
		    print "Warning: Cannot key more than one chrom alias by " . $rowData[$normalChromNameCol] . "\n";
		}
	    } else {
		unless ($chromAliases{$rowData[$chromNameCol]}) {
		    $chromAliases{$rowData[$chromNameCol]} = $rowData[$chromAliasCol];
		} else {
		    print "Warning: Cannot key more than one chrom alias by " . $rowData[$chromNameCol] . "\n";
		}
	    }
	}
    }
    close(AFLE);
    return \%chromAliases;
}

sub getChromAliases_withNormalChromNames_withNormalChromPrefix_pref {
    my($assemblyInfoFileLocation,$gDir,$chromNameCol,$chromAliasCol,$normalChromNameCol,$normalChromNamePrefix,$pref) = @_;
    my $assemblyInfoFile = $assemblyInfoFileLocation;
    if ((uc($assemblyInfoFileLocation) =~ /^HTTP:/) || (uc($assemblyInfoFileLocation) =~ /^FTP:/)) {
	($assemblyInfoFile) = $assemblyInfoFileLocation =~ /([^\/]*$)/;
	$assemblyInfoFile = "$gDir/$assemblyInfoFile";
	unless ((-e $assemblyInfoFile)) {
	    my $cmd = "wget $assemblyInfoFileLocation --directory-prefix=$gDir";
	    my $out = `$cmd`;
	    print $out;
	    my $cmd2 = "sed -i \'s/\r//g\' $assemblyInfoFile";
	    my $out2 = `$cmd2`;
	    print $out2;
	}
    }
    my %chromAliases;
    open(AFLE,$assemblyInfoFile) or die "failed to open $assemblyInfoFile\n";
    while (<AFLE>) {
	chomp;
	unless ( /^\#/ ) {
	    my(@rowData) = split(/\t/);
	    if ($rowData[3] eq "Chromosome" && ($rowData[$chromAliasCol] =~ /^\Q$pref\E/)) {
		#Assigned-molecule is labeled as a chromosome
		unless ($chromAliases{"$normalChromNamePrefix$rowData[$normalChromNameCol]"}) {
		    $chromAliases{"$normalChromNamePrefix$rowData[$normalChromNameCol]"} = $rowData[$chromAliasCol];
		} else {
		    print "Warning: Cannot key more than one chrom alias by " . $rowData[$normalChromNameCol] . "\n";
		}
	    } else {
		unless ($chromAliases{$rowData[$chromNameCol]}) {
		    $chromAliases{$rowData[$chromNameCol]} = $rowData[$chromAliasCol];
		} else {
		    print "Warning: Cannot key more than one chrom alias by " . $rowData[$chromNameCol] . "\n";
		}
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

sub prepareFasta_MES {
    my($genomeFileLocation,$gDir) = @_;
    my $fastaFile = $genomeFileLocation;
    if ((uc($genomeFileLocation) =~ /^HTTP:/) || (uc($genomeFileLocation) =~ /^HTTPS:/) || (uc($genomeFileLocation) =~ /^FTP:/)) {
	my($fileName) = $genomeFileLocation =~ /([^\/]*mes\.con.*$)/;
	if ($fileName =~ /\.gz$/) {
	    ($fastaFile) = $fileName =~ /(^.*?)\.gz/;
	    $fastaFile = "$gDir/$fastaFile";
	} else {
	    $fastaFile = "$gDir/$fileName";
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

sub prepareFasta {
    my($genomeFileLocation,$gDir) = @_;
    my $fastaFile = $genomeFileLocation;
    if ((uc($genomeFileLocation) =~ /^HTTP:/) || (uc($genomeFileLocation) =~ /^HTTPS:/) || (uc($genomeFileLocation) =~ /^FTP:/)) {
	my($fileName) = $genomeFileLocation =~ /([^\/]*\.f.*$)/;
	if ($fileName =~ /\.gz$/) {
	    ($fastaFile) = $fileName =~ /(^.*?)\.gz/;
	    $fastaFile = "$gDir/$fastaFile";
	} else {
	    $fastaFile = "$gDir/$fileName";
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
