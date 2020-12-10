
# Dataset Download and Preparation

If you plan on running these scripts, please add the PerlModules/miRBaseDataGatherer.pm to your PERL5LIB path.

The following dependencies are required for these scripts to work properly:
* [The ViennaRNA Package](https://www.tbi.univie.ac.at/RNA/)
* [samtools](http://www.htslib.org/download/)
* [Emboss](ftp://emboss.open-bio.org/pub/EMBOSS/)
* [cd-hit](https://github.com/weizhongli/cdhit)


The downloadSequenceData.sh script will download all the neccessary files and extract information about each organism microRNAs (download includes entire genome for each organsim so it is around 100GB in size.)

```sh
$ bash downloadSequenceData.sh
```

Once the download is complete, prepareInitialDataset.sh can be used to prepare the deep-mir-cut dataset.

```sh
$ bash prepareInitialDataset.sh
```

# Details of Run Scripts

## downloadSequenceData.sh

Sequences from miRBase are used to verify that the correct sequences are being extracted from each species genome.  The downloadSequenceData.sh script downloads these sequences from the miRBase server using the following commands:

```sh
$ wget ftp://mirbase.org/pub/mirbase/22/hairpin.fa.gz
$ gunzip hairpin.fa.gz
$ wget ftp://mirbase.org/pub/mirbase/22/mature.fa.gz
$ gunzip mature.fa.gz
```

Next downloadSequenceData.sh runs a set of scripts to download and extract data for microRNAs from each species into the miRSequenceOutput folder.  Each script has the abbreviated species name as it's prefix.

```sh
$ mkdir miRSequenceOutput
$ perl Scripts/aae_extractSequences.pl | tee miRSequenceOutput/aae_warnings.txt
$ perl Scripts/abu_extractSequences.pl | tee miRSequenceOutput/abu_warnings.txt
$ perl Scripts/aca_extractSequences.pl | tee miRSequenceOutput/aca_warnings.txt
$ perl Scripts/aga_extractSequences.pl | tee miRSequenceOutput/aga_warnings.txt
$ perl Scripts/aly_extractSequences.pl | tee miRSequenceOutput/aly_warnings.txt
$ perl Scripts/ame_extractSequences.pl | tee miRSequenceOutput/ame_warnings.txt
...................see downloadSequenceData.sh.................................
```


## prepareInitialDataset.sh

microRNA data that was placed in the miRSequenceOutput directory is combined using the following command:
```sh
$ perl Scripts/combineMiRSequenceOutput.pl miRSequenceOutput
```

Next, hairpin precursor sequences from miRBase are folded using RNAfold.  These folds will be used to determine if microRNA are on the 5p or 3p arm.  A script called cleanFoldsFile.pl was added to the pipeline because some versions of RNAfold allow additional information on the defline which may cause other scripts to throw errors.

```sh
$ RNAfold hairpin.fa --noLP --noPS > hairpin.folds
$ perl Scripts/cleanFoldsFile.pl hairpin.folds
```

MicroRNAs with inconsistencies between files found on miRBase were dropped from the dataset.

```sh
$ perl Scripts/dropInconsistentMirs.pl mir_sequences.gff mir_seqInfo.txt mir_sequences.fa inconsistentMirs.txt > dropInconsistentMirs.out
```

All the information is put toghether into one dataset using a script called createDicerCutsiteSet.pl.  The folded hairpin.fa file is used to determine the arm of the hairpin that each mature microRNA is on.  Dicer and drosha cutsites are identified using each microRNAs arm and start and stop location.  In rare cases where the hairpin is folded such that a microRNA is in the loop, no cutsites will be recorded since the hairpin arm is ambiguous.


```sh
$ perl Scripts/createDicerCutsiteSet.pl new_mir_sequences.fa new_mir_sequences.gff new_mir_seqInfo.txt hairpin.folds > createDicerCutsiteSet.out
```
