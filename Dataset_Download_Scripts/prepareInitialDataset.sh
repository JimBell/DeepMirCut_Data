#!/bin/bash
perl Scripts/combineMiRSequenceOutput.pl miRSequenceOutput
RNAfold hairpin.fa --noLP --noPS > hairpin.folds
perl Scripts/dropInconsistentMirs.pl mir_sequences.gff mir_seqInfo.txt mir_sequences.fa inconsistentMirs.txt > dropInconsistentMirs.out
perl Scripts/createDicerCutsiteSet.pl new_mir_sequences.fa new_mir_sequences.gff new_mir_seqInfo.txt hairpin.folds > createDicerCutsiteSet.out
