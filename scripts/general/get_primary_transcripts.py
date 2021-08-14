import os
import HTSeq
import pandas as pd

os.chdir("/media/jannik/Samsung_T5/masters_project/ApisTranscriptomics/data/reference/immune_gene_identification")
try:
	os.mkdir("fasta/primary_transcripts")
except:
	None
	
def get_primary_transcripts(in_fasta, out_fasta, gene2prot):
	gene2prot = pd.read_csv(gene2prot, names = ["GeneID", "ProteinID"])
	gene2prot['GeneID'] = gene2prot['GeneID'].str.replace('GeneID:', '')
	gene2prot['ProteinID'] = gene2prot['ProteinID'].str.replace('Genbank:', '')
	
	sequences = dict( (seq.name, seq) for seq in HTSeq.FastaReader(in_fasta) )
	out_fasta_file = open(out_fasta, "w")
	
	for geneID in gene2prot['GeneID'].unique():
		proteins = gene2prot.loc[gene2prot['GeneID'] == geneID, 'ProteinID'].tolist()
		longest_protein = proteins[0]
		longest_protein_len = len(sequences[proteins[0]])
		for protein in proteins[1:]:
			if len(sequences[protein]) > longest_protein_len:
				longest_protein = protein
				longest_protein_len = len(sequences[protein])
		sequences[longest_protein].write_to_fasta_file(out_fasta_file)
	out_fasta_file.close()
	
	out_sequences = [seq.name for seq in HTSeq.FastaReader(out_fasta)]

	
	print(os.path.split(in_fasta)[1], ":", len(sequences), "proteins read in;", len(gene2prot['ProteinID'].unique()), "protein accessions in gff;", len(out_sequences), "primary transcripts written out;", len(gene2prot['GeneID'].unique()), "gene accessions in gff")


get_primary_transcripts("fasta/proteins/Acyrthosiphon_pisum.faa", "fasta/primary_transcripts/Acyrthosiphon_pisum.faa", "gene2prot/Acyrthosiphon_pisum.txt")
get_primary_transcripts("fasta/proteins/Aedes_aegypti.faa", "fasta/primary_transcripts/Aedes_aegypti.faa", "gene2prot/Aedes_aegypti.txt")
get_primary_transcripts("fasta/proteins/Anopheles_gambiae.faa", "fasta/primary_transcripts/Anopheles_gambiae.faa", "gene2prot/Anopheles_gambiae.txt")
get_primary_transcripts("fasta/proteins/Apis_mellifera.faa", "fasta/primary_transcripts/Apis_mellifera.faa", "gene2prot/Apis_mellifera.txt")
get_primary_transcripts("fasta/proteins/Bombus_terrestris.faa", "fasta/primary_transcripts/Bombus_terrestris.faa", "gene2prot/Bombus_terrestris.txt")
get_primary_transcripts("fasta/proteins/Bombyx_mori.faa", "fasta/primary_transcripts/Bombyx_mori.faa", "gene2prot/Bombyx_mori.txt")
get_primary_transcripts("fasta/proteins/Ceratina_calcarata.faa", "fasta/primary_transcripts/Ceratina_calcarata.faa", "gene2prot/Ceratina_calcarata.txt")
get_primary_transcripts("fasta/proteins/Drosophila_melanogaster.faa", "fasta/primary_transcripts/Drosophila_melanogaster.faa", "gene2prot/Drosophila_melanogaster.txt")
get_primary_transcripts("fasta/proteins/Dufourea_novaeangliae.faa", "fasta/primary_transcripts/Dufourea_novaeangliae.faa", "gene2prot/Dufourea_novaeangliae.txt")
get_primary_transcripts("fasta/proteins/Eufriesea_mexicana.faa", "fasta/primary_transcripts/Eufriesea_mexicana.faa", "gene2prot/Eufriesea_mexicana.txt")
get_primary_transcripts("fasta/proteins/Habropoda_laboriosa.faa", "fasta/primary_transcripts/Habropoda_laboriosa.faa", "gene2prot/Habropoda_laboriosa.txt")
get_primary_transcripts("fasta/proteins/Megachile_rotundata.faa", "fasta/primary_transcripts/Megachile_rotundata.faa", "gene2prot/Megachile_rotundata.txt")
get_primary_transcripts("fasta/proteins/Megalopta_genalis.faa", "fasta/primary_transcripts/Megalopta_genalis.faa", "gene2prot/Megalopta_genalis.txt")
get_primary_transcripts("fasta/proteins/Nasonia_vitripennis.faa", "fasta/primary_transcripts/Nasonia_vitripennis.faa", "gene2prot/Nasonia_vitripennis.txt")
get_primary_transcripts("fasta/proteins/Nomia_melanderi.faa", "fasta/primary_transcripts/Nomia_melanderi.faa", "gene2prot/Nomia_melanderi.txt")
get_primary_transcripts("fasta/proteins/Osmia_bicornis.faa", "fasta/primary_transcripts/Osmia_bicornis.faa", "gene2prot/Osmia_bicornis.txt")
get_primary_transcripts("fasta/proteins/Osmia_lignaria.faa", "fasta/primary_transcripts/Osmia_lignaria.faa", "gene2prot/Osmia_lignaria.txt")
get_primary_transcripts("fasta/proteins/Polistes_dominula.faa", "fasta/primary_transcripts/Polistes_dominula.faa", "gene2prot/Polistes_dominula.txt")
get_primary_transcripts("fasta/proteins/Solenopsis_invicta.faa", "fasta/primary_transcripts/Solenopsis_invicta.faa", "gene2prot/Solenopsis_invicta.txt")
get_primary_transcripts("fasta/proteins/Tribolium_castaneum.faa", "fasta/primary_transcripts/Tribolium_castaneum.faa", "gene2prot/Tribolium_castaneum.txt")
get_primary_transcripts("fasta/proteins/Vespa_mandarinia.faa", "fasta/primary_transcripts/Vespa_mandarinia.faa", "gene2prot/Vespa_mandarinia.txt")
