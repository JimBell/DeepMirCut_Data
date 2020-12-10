
# deep-mir-cut Train, Validation, and Test Set Generation

The generateTrainTestSets.sh script may be used to create training, validation, and testing sets.  

First make sure you have all dependencies installed:
* [The ViennaRNA Package](https://www.tbi.univie.ac.at/RNA/)
* [bpRNA](https://github.com/hendrixlab/bpRNA)
* [Emboss](ftp://emboss.open-bio.org/pub/EMBOSS/)
* [cd-hit](https://github.com/weizhongli/cdhit)
* [Graph.pm perl module](https://metacpan.org/pod/distribution/Graph/lib/Graph.pod)

Place a copy of the bpRNA script into the scripts directory.

Place the microRNA dataset file (mir_dataset.txt) into this directory and run the generateTrainTestSets.sh script to produce the train, validation, and test sets.

```sh
$ bash generateTrainTestSets.sh
```

# Details of generateTrainTestSets.sh Script

generateTrainTestSets.sh downloads the hairpin fasta from miRBase and the organisms file and runs scripts to find all microRNAs with 2 products.  The cd-hit command is used to filter these down to just microRNA within 80% identity with one another.  The precursors identified by cd-hit are divided up by the type of organism.  Precursor names from metazoan species are passed into a script called createTrainTestSets.pl which divides them up into the train, test, and validation sets.

```sh
$ wget ftp://mirbase.org/pub/mirbase/22/hairpin.fa.gz
$ wget ftp://mirbase.org/pub/mirbase/22/organisms.txt.gz
$ gunzip hairpin.fa.gz
$ gunzip organisms.txt.gz
$ perl Scripts/getDoubleProductHairpinFasta.pl hairpin.fa mir_dataset.txt hairpin_2prod.fa
$ cd-hit -c 0.8 -b 1000 -i hairpin_2prod.fa -o hairpin_cdHit.txt
$ perl Scripts/getNamesFromCDHitFasta.pl hairpin_cdHit.txt names.txt
$ perl Scripts/divideNamesByOrganism.pl mir_dataset.txt names.txt organisms.txt
$ perl Scripts/createTrainTestSets.pl mir_dataset.txt names_Metazoa.txt Metazoa
```

A script called augmentWithSimilarSequences2.pl adds precursors that were filtered out by cd-hit back into the train set as long as they have less than 80% identity with the validation and testing sets.  It also adds a random 30 to 50 nt buffer region around the precursor.  Further data augmentation is done with augmentWithRandomCuts.pl which adds 9x additional examples with random 30 to 50 nt buffer regions around the precursor regions. 

```sh
$perl Scripts/augmentWithSimilarSequences2.pl mir_dataset.txt Metazoa_trainSet.txt Metazoa_validationSet.txt Metazoa_testSet.txt hairpin_cdHit.txt.clstr hairpin.fa 
$perl Scripts/augmentWithRandomCuts.pl mir_dataset.txt Metazoa_trainSet_wSimilar.txt
$perl Scripts/augmentWithRandomCuts.pl mir_dataset.txt Metazoa_validationSet.txt
```

Dot-bracket folds and the bpRNA structure arrays are added to each of the datasets using addFoldsToSet.pl. 

```sh
perl Scripts/addFoldsToSet.pl Metazoa_trainSet_wSimilar_mult.txt
perl Scripts/addFoldsToSet.pl Metazoa_validationSet_mult.txt
perl Scripts/addFoldsToSet.pl Metazoa_testSet.txt
```
