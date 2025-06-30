# Before starting this, make sure all the gene IDs across the two CDS fastas and orthologues match

from Bio import SeqIO

# Filepaths
lyrata_fasta = "your_cds_fromGTF_matchingOrthoIDs.fasta"
thaliana_fasta = "your_cds_matchingOrthoIDs.fasta"
ortholog_file = "your_2species_orthologs.txt"

# Load all CDS into dictionaries for quick access
lyrata_cds = SeqIO.to_dict(SeqIO.parse(lyrata_fasta, "fasta"))
thaliana_cds = SeqIO.to_dict(SeqIO.parse(thaliana_fasta, "fasta"))

with open(ortholog_file) as orths:
    for line in orths:
        lyrata_id, thaliana_id = line.strip().split()
        if lyrata_id not in lyrata_cds or thaliana_id not in thaliana_cds:
            continue  # Skip missing genes
        # Output FASTA for this pair
        with open(f"{lyrata_id}_vs_{thaliana_id}.fa", "w") as out:
            SeqIO.write([lyrata_cds[lyrata_id], thaliana_cds[thaliana_id]], out, "fasta")



