#!/bin/bash

#ancient DNA sequence pre-processing script 
#trim and merge all samples from a NovaSeq Run
#adapted from Jonas
#edits Kim Ballare 21-12-14

ids=/redser4/projects/kim/TRBL_list.txt #list of raw sequence filenames ending with S number ex. sample1_S1
rawdir=/redser4/raw/211026_Duke_NovaSeq_bison_trbb/Johnstone_7226_211018B6 #raw sequence directory

N=2

#SeqPrep2 Variables
SEQPREP_MIN_LENGTH=30    #Min length to merge <-VARIABLE -L
SEQPREP_OVERLAP=15       #Overlap variable for merging <-VARIABLE -o
SEQPREP_MINQUAL=15       #Minumum quality score
MISMATCH_FRAC=0.05

while read -r id; do

        s=${id%S*}
    samp=${s%_*}

    ((i=i%N)); ((i++==0)) && wait

        mkdir -p ${samp} && cd ${samp} || exit
        for lane in {L001,L002,L003,L004}; do
                /redser4/personal/jonas/bin/SeqPrep2 -f ${rawdir}/${id}_${lane}_R1_001.perfect.fastq.gz -r ${rawdir}/${id}_${lane}_R2_001.perfect.fastq.gz \
                -1 ${samp}_${lane}_R1_unmerged.fastq.gz -2 ${samp}_${lane}_R2_unmerged.fastq.gz \
                -q ${SEQPREP_MINQUAL} -L ${SEQPREP_MIN_LENGTH} \
                -A AGATCGGAAGAGCACACGTC -B AGATCGGAAGAGCGTCGTGT \
                -s ${samp}_${lane}_merged.fastq.gz \
                -o ${SEQPREP_OVERLAP} -m ${MISMATCH_FRAC} \
                -C ATCTCGTATGCCGTCTTCTGCTTG -D GATCTCGGTGGTCGCCGTATCATT -S >& ${samp}_${lane}.seqprep_output.txt &
        done
    cd ..

done < $ids
