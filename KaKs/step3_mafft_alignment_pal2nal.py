import sys, os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import subprocess

if len(sys.argv) < 2:
    print("Usage: python run_align_pal2nal.py <input_fasta>")
    sys.exit(1)

cds_path = sys.argv[1]
output_dir = "/your/output/dir" # REPLACE it here
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

fname = os.path.basename(cds_path)
base = fname[:-3]

# 1. Translate CDS to protein
records = list(SeqIO.parse(cds_path, "fasta"))
prot_records = []
for rec in records:
    prot_seq = Seq(str(rec.seq)).translate(to_stop=True)
    prot_rec = SeqRecord(prot_seq, id=rec.id, description="")
    prot_records.append(prot_rec)
prot_path = os.path.join(output_dir, base + ".prot.fa")
SeqIO.write(prot_records, prot_path, "fasta")

# 2. Align protein
aligned_prot_path = os.path.join(output_dir, base + ".prot.aln.fa")
with open(aligned_prot_path, "w") as fout:
    subprocess.run(["mafft", "--auto", prot_path], stdout=fout)

# 3. Codon alignment with PAL2NAL
aligned_cds_path = os.path.join(output_dir, base + ".codon.aln.fa")
with open(aligned_cds_path, "w") as fout:
    subprocess.run([
        "/your/path/to/pal2nal.pl", # REPLACE it here
        aligned_prot_path,
        cds_path,
        "-output", "fasta"
    ], stdout=fout)
