#!/bin/bash
wget ftp://mirbase.org/pub/mirbase/22/hairpin.fa.gz
wget ftp://mirbase.org/pub/mirbase/22/organisms.txt.gz
gunzip hairpin.fa.gz
gunzip organisms.txt.gz
perl Scripts/getDoubleProductHairpinFasta.pl hairpin.fa mir_dataset.txt hairpin_2prod.fa
cd-hit -c 0.8 -b 1000 -i hairpin_2prod.fa -o hairpin_cdHit.txt
perl Scripts/getNamesFromCDHitFasta.pl hairpin_cdHit.txt names.txt
perl Scripts/divideNamesByOrganism.pl mir_dataset.txt names.txt organisms.txt
perl Scripts/createTrainTestSets.pl mir_dataset.txt names_Metazoa.txt Metazoa

#Augmenting Data
perl Scripts/augmentWithSimilarSequences2.pl mir_dataset.txt Metazoa_trainSet.txt Metazoa_validationSet.txt Metazoa_testSet.txt hairpin_cdHit.txt.clstr hairpin.fa 
perl Scripts/augmentWithRandomCuts.pl mir_dataset.txt Metazoa_trainSet_wSimilar.txt
perl Scripts/augmentWithRandomCuts.pl mir_dataset.txt Metazoa_validationSet.txt

#Adding Folds
perl Scripts/addFoldsToSet.pl Metazoa_trainSet_wSimilar_mult.txt
perl Scripts/addFoldsToSet.pl Metazoa_validationSet_mult.txt
perl Scripts/addFoldsToSet.pl Metazoa_testSet.txt
