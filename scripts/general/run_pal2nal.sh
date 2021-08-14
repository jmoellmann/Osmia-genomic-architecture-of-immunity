#!/bin/bash

PROJECT_DIR=/media/jannik/Samsung_T5/masters_project/ApisTranscriptomics
PHYLO_DIR=${PROJECT_DIR}/data/reference/immune_gene_identification
SINGLE_COPY_ORTHOGROUPS=${PROJECT_DIR}/data/reference/single_copy_orthologues.txt
CDS_ALIGNMENTS_DIR=${PHYLO_DIR}/fasta/cds_alignments

mkdir -p ${CDS_ALIGNMENTS_DIR}

for ORTHOGROUP in $(cat ${SINGLE_COPY_ORTHOGROUPS});
do
	perl ${PROJECT_DIR}/src/pal2nal.v14/pal2nal.pl ${PHYLO_DIR}/fasta/single_copy_orthologues_alignments/${ORTHOGROUP}.fa ${PHYLO_DIR}/fasta/cds_for_alignments/${ORTHOGROUP}.fa -output fasta > ${CDS_ALIGNMENTS_DIR}/${ORTHOGROUP}.fa &
done

wait
