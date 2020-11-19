#!/bin/bash
perl Scripts/combineMiRSequenceOutput.pl miRSequenceOutput
wget ftp://mirbase.org/pub/mirbase/22/hairpin.fa.gz
gunzip hairpin.fa.gz
RNAfold hairpin.fa --noLP --noPS > hairpin.folds
perl Scripts/dropInconsistentMirs.pl mir_sequences.gff mir_seqInfo.txt mir_sequences.fa inconsistentMirs.txt > dropInconsistentMirs.out
perl Scripts/createDicerCutsiteSet.pl new_mir_sequences.fa new_mir_sequences.gff new_mir_seqInfo.txt hairpin.folds > createDicerCutsiteSet.out
