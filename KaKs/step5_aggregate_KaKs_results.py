import os
import glob

kaks_dir = "/your/path/to/kaks_results"
output_table = "/your/path/to/kaks_results_aggregated.tsv"

# New header with P-value column
header = "Lyrata_ID\tThaliana_ID\tKa\tKs\tKa/Ks\tP_value\n"

with open(output_table, "w") as fout:
    fout.write(header)
    for fname in glob.glob(os.path.join(kaks_dir, "*.kaks.txt")):
        with open(fname) as fin:
            for line in fin:
                # Skip header lines and blanks
                if line.startswith("Sequence") or not line.strip():
                    continue
                fields = line.strip().split('\t')
                # We expect at least: SeqA, SeqB, Ka, Ks, Ka/Ks, P-value
                if len(fields) < 6:
                    continue
                # Derive gene IDs from filename
                base = os.path.basename(fname).replace(".kaks.txt", "")
                try:
                    lyrata_id, thaliana_id = base.split("_vs_")
                except ValueError:
                    lyrata_id, thaliana_id = fields[0], fields[1]
                ka, ks, ratio = fields[2], fields[3], fields[4]
                pval = fields[5]
                fout.write(f"{lyrata_id}\t{thaliana_id}\t{ka}\t{ks}\t{ratio}\t{pval}\n")


