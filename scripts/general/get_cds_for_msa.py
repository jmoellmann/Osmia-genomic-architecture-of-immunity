import os
import HTSeq
import pandas as pd
import re 

os.chdir("/media/jannik/Samsung_T5/masters_project/ApisTranscriptomics/")

def get_cds_for_msa(msa_path, coding_sequences):
   out_dir = "data/reference/immune_gene_identification/fasta/cds_for_alignments/"
   out_fasta_file = open(out_dir + os.path.split(msa_path)[-1], "w")
   sequence_ids = [seq.name for seq in HTSeq.FastaReader(msa_path)]
   for sequence_id in sequence_ids:
      try:
         coding_sequences[sequence_id].write_to_fasta_file(out_fasta_file)
      except Exception:
         out_fasta_file.write("ERROR\n")
   out_fasta_file.close()
   
def get_cds_for_msa_with_mrna(msa_path, coding_sequences):
   out_dir = "data/reference/immune_gene_identification/fasta/cds_for_alignments/"
   out_fasta_file = open(out_dir + os.path.split(msa_path)[-1], "w")
   sequence_ids = [seq.name for seq in HTSeq.FastaReader(msa_path)]
   for sequence_id in sequence_ids:
      accession = prot2transcript[sequence_id][1]
      try:
         coding_sequences[accession].write_to_fasta_file(out_fasta_file)
      except Exception:
         out_fasta_file.write("ERROR\n")
   out_fasta_file.close()
   
def gather_cds_files(cds_file_paths):
   coding_sequences = dict()
   for cds_file_path in cds_file_paths:
      species_name = os.path.split(cds_file_path)[-1][:-3]
      temp_dict_v1 = dict( (seq.name, seq) for seq in HTSeq.FastaReader(cds_file_path) )
      temp_dict_v2 = dict()
      for key, value in temp_dict_v1.items():
         temp_key = re.search("[XNY]P_[0-9]+\\.[0-9]", key)[0]
         temp_seq = value
         temp_seq.name = species_name + "_" + temp_seq.name[:-2] + "_" + temp_seq.name[:-1]
         temp_dict_v2[temp_key] = temp_seq
      coding_sequences.update(temp_dict_v2)
   return coding_sequences

def gather_cds_files_from_mrna(cds_file_paths):
   coding_sequences = dict()
   for cds_file_path in cds_file_paths:
      species_name = os.path.split(cds_file_path)[-1][:-3]
      temp_dict = dict( (seq.name, seq) for seq in HTSeq.FastaReader(cds_file_path) )
      for key, value in temp_dict.items():
         temp_seq = value
         temp_seq.name = species_name + "_" + temp_seq.name[:-2] + "_" + temp_seq.name[:-1]
         temp_dict[key] = temp_seq
      coding_sequences.update(temp_dict)
   return coding_sequences

# cds_from_mrna_dir = "data/reference/immune_gene_identification/fasta/transcripts/"
# cds_from_mrna_file_paths = [cds_from_mrna_dir + i for i in os.listdir(cds_from_mrna_dir)]
# coding_sequences_from_mrna = gather_cds_files(cds_from_mrna_file_paths)

cds_dir = "data/reference/immune_gene_identification/fasta/cds_from_genomic/"
cds_file_paths = [cds_dir + i for i in os.listdir(cds_dir)]
coding_sequences = gather_cds_files(cds_file_paths)

msa_dir = "data/reference/immune_gene_identification/fasta/single_copy_orthologues_alignments/"

prot2transcript_df = pd.read_csv("data/reference/immune_gene_identification/prot2transcript.tsv", 
   sep = "\t", header = None)
prot2transcript_df = prot2transcript_df.set_index(0)
prot2transcript = prot2transcript_df.to_dict(orient = "index")

try:
   os.mkdir("data/reference/immune_gene_identification/fasta/cds_for_alignments/")
except:
   None

for msa in os.listdir(msa_dir):
   get_cds_for_msa(msa_dir + msa, coding_sequences)

