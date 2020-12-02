#!/usr/bin/perl -w
use miRBaseDataGatherer;
use strict;

my $mirBaseVersion = "22.1";
my $genomeId = "v1.1";
my $speciesPrefix = "fve";
my $miRBaseGff = "$speciesPrefix\.gff3";
my $genomeFileLocation = "ftp://ftp.bioinfo.wsu.edu/species/Fragaria_vesca/Fvesca-genome.v1.1/assembly/*";
my $assemblyInfoFileLocation = "";
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
my $fastaFile = prepareFasta_fve($genomeFileLocation,$gDir);
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

sub prepareFasta {
    my($genomeFileLocation,$gDir) = @_;
    my $fastaFile = $genomeFileLocation;
    if ((uc($genomeFileLocation) =~ /^HTTP:/) || (uc($genomeFileLocation) =~ /^FTP:/)) {
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

sub prepareFasta_fve {
    my($genomeFileLocation,$gDir) = @_;
    my $fastaFile = "$gDir/fve.fa";
    unless (-e "$gDir/fve.fa") {
	my $cmd = "wget -O $gDir/fvesca_v1.1_LG1.fna.gz ftp://ftp.bioinfo.wsu.edu/species/Fragaria_vesca/Fvesca-genome.v1.1/assembly/fvesca_v1.1_LG1.fna.gz";
	my $out = `$cmd`;
	print $out;
	my $cmd1a = "gunzip $gDir/fvesca_v1.1_LG1.fna.gz";
	my $out1a = `$cmd1a`;
	print $out1a;
	my $cmd2 = "wget -O $gDir/fvesca_v1.1_LG2.fna.gz ftp://ftp.bioinfo.wsu.edu/species/Fragaria_vesca/Fvesca-genome.v1.1/assembly/fvesca_v1.1_LG2.fna.gz";
	my $out2 = `$cmd2`;
	print $out2;
	my $cmd2a = "gunzip $gDir/fvesca_v1.1_LG2.fna.gz";
	my $out2a = `$cmd2a`;
	print $out2a;
	my $cmd3 = "wget -O $gDir/fvesca_v1.1_LG3.fna.gz ftp://ftp.bioinfo.wsu.edu/species/Fragaria_vesca/Fvesca-genome.v1.1/assembly/fvesca_v1.1_LG3.fna.gz";
	my $out3 = `$cmd3`;
	print $out3;
	my $cmd3a = "gunzip $gDir/fvesca_v1.1_LG3.fna.gz";
	my $out3a = `$cmd3a`;
	print $out3a;
	my $cmd4 = "wget -O $gDir/fvesca_v1.1_LG4.fna.gz ftp://ftp.bioinfo.wsu.edu/species/Fragaria_vesca/Fvesca-genome.v1.1/assembly/fvesca_v1.1_LG4.fna.gz";
	my $out4 = `$cmd4`;
	print $out4;
	my $cmd4a = "gunzip $gDir/fvesca_v1.1_LG4.fna.gz";
	my $out4a = `$cmd4a`;
	print $out4a;
	my $cmd5 = "wget -O $gDir/fvesca_v1.1_LG5.fna.gz ftp://ftp.bioinfo.wsu.edu/species/Fragaria_vesca/Fvesca-genome.v1.1/assembly/fvesca_v1.1_LG5.fna.gz";
	my $out5 = `$cmd5`;
	print $out5;
	my $cmd5a = "gunzip $gDir/fvesca_v1.1_LG5.fna.gz";
	my $out5a = `$cmd5a`;
	print $out5a;
	my $cmd6 = "wget -O $gDir/fvesca_v1.1_LG6.fna.gz ftp://ftp.bioinfo.wsu.edu/species/Fragaria_vesca/Fvesca-genome.v1.1/assembly/fvesca_v1.1_LG6.fna.gz";
	my $out6 = `$cmd6`;
	print $out6;
	my $cmd6a = "gunzip $gDir/fvesca_v1.1_LG6.fna.gz";
	my $out6a = `$cmd6a`;
	print $out6a;
	my $cmd7 = "wget -O $gDir/fvesca_v1.1_LG7.fna.gz ftp://ftp.bioinfo.wsu.edu/species/Fragaria_vesca/Fvesca-genome.v1.1/assembly/fvesca_v1.1_LG7.fna.gz";
	my $out7 = `$cmd7`;
	print $out7;
	my $cmd7a = "gunzip $gDir/fvesca_v1.1_LG7.fna.gz";
	my $out7a = `$cmd7a`;
	print $out7a;
	my $cmd8 = "wget -O $gDir/fvesca_v1.1_unanchored.fna.gz ftp://ftp.bioinfo.wsu.edu/species/Fragaria_vesca/Fvesca-genome.v1.1/assembly/fvesca_v1.1_unanchored.fna.gz";
	my $out8 = `$cmd8`;
	print $out8;
	my $cmd8a = "gunzip $gDir/fvesca_v1.1_unanchored.fna";
	my $out8a = `$cmd8a`;
	print $out8a;
	my $cmd9 = "cat $gDir/fvesca_v1.1_unanchored.fna $gDir/fvesca_v1.1_LG1.fna $gDir/fvesca_v1.1_LG2.fna $gDir/fvesca_v1.1_LG3.fna $gDir/fvesca_v1.1_LG4.fna $gDir/fvesca_v1.1_LG5.fna $gDir/fvesca_v1.1_LG6.fna $gDir/fvesca_v1.1_LG7.fna > $gDir/fve.fa";
	my $out9 = `$cmd9`;
	print $out9;
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
