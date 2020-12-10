
# deep-mir-cut Dataset Download and Preparation

If you plan on running these scripts to download and prepare the data yourself, please add the PerlModules/miRBaseDataGatherer.pm to your PERL5LIB path.

The downloadSequenceData.sh script may be run to download all the neccessary files and extract information about each organism microRNAs (download includes entire genome for each organsim so it is around 100GB in size.)

```sh
$ bash downloadSequenceData.sh
```

Once the download is complete, the prepareInitialDataset.sh script can be used to prepare the deep-mir-cut dataset.

```sh
$ bash prepareInitialDataset.sh
```

# Details of Run Scripts

## downloadSequenceData.sh

Sequences from miRBase are used to verify that the genomic regions where sequences are being extracted contain the microRNAs assumed to be at that loci.  The downloadSequenceData.sh script downloads these sequences from the miRBase server using the following commands:

```sh
$ wget ftp://mirbase.org/pub/mirbase/22/hairpin.fa.gz
$ gunzip hairpin.fa.gz
$ wget ftp://mirbase.org/pub/mirbase/22/mature.fa.gz
$ gunzip mature.fa.gz
```

Next downloadSequenceData.sh runs a set of scripts to download and extract data for miRs from each species into the miRSequenceOutput folder.  Each script has the abbreviated species name as it's prefix.

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

Next, RNAfold is used to fold the hairpin precursors that were downloaded from miRBase.  These folds will be used to determine if microRNA are on the 5p or 3p arm.  Because some versions of RNAfold let additional information on the defline of the fasta file to be included in the output, a script called cleanFoldsFile.pl was added to prevent other scripts in the pipeline from throwing errors.

```sh
$ RNAfold hairpin.fa --noLP --noPS > hairpin.folds
$ perl Scripts/cleanFoldsFile.pl hairpin.folds
```

Some microRNAs were inconsistent between files found on miRBase.  These microRNAs were dropped from the dataset.

```sh
$ perl Scripts/dropInconsistentMirs.pl mir_sequences.gff mir_seqInfo.txt mir_sequences.fa inconsistentMirs.txt > dropInconsistentMirs.out
```

A script called createDicerCutsiteSet.pl puts toghether all the information into one dataset (mir_dataset.txt).  The folded hairpin.fa file is used to determine the arm of the hairpin that each mature microRNA is on.  Dicer and drosha cutsites are identified using each microRNAs arm and start and stop location.  In rare cases where the hairpin is folded such that a microRNA is in the loop, no cutsites will be recorded because the arm is ambiguous.

