#!/usr/bin/env python3
from Bio import SeqIO
import gffutils
import os

genome_fasta = "your_input.fasta"
gtf_file     = "your_input.gtf"
output_fasta = "your_input_CDS.fasta"

# 1. Load genome
genome = SeqIO.to_dict(SeqIO.parse(genome_fasta, "fasta"))

# 2. Build an in‐memory database from the GTF
db = gffutils.create_db(
    gtf_file,
    dbfn=":memory:",
    force=True,
    keep_order=True,
    merge_strategy="merge",
    sort_attribute_values=True,
    disable_infer_transcripts=False,
    disable_infer_genes=False
)

# 3. Extract CDS per transcript
with open(output_fasta, "w") as out_fasta:
    # In GTF the parent feature for CDS is usually 'transcript'
    for transcript in db.features_of_type("transcript"):
        cds_parts = list(db.children(transcript, featuretype="CDS", order_by="start"))
        if not cds_parts:
            continue
        # concatenate exons in genomic order, fixing strand
        seq_fragments = []
        for cds in cds_parts:
            seq = genome[cds.seqid].seq[cds.start - 1 : cds.end]
            if cds.strand == "-":
                seq = seq.reverse_complement()
            seq_fragments.append(str(seq))
        cds_seq = "".join(seq_fragments)

        # write FASTA: use transcript_id attribute (fallback to transcript.id)
        tid = transcript.attributes.get("transcript_id", [transcript.id])[0]
        out_fasta.write(f">{tid}\n")
        for i in range(0, len(cds_seq), 60):
            out_fasta.write(cds_seq[i : i + 60] + "\n")
