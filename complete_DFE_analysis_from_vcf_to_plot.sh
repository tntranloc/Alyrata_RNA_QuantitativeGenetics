

##### Required softwares/packages #####
# snpEff, angsd, python with pandas and fastdfe installed




###### Find synonymous and non-synonymous variant sites ######

# install snpeff if needed
java -jar snpEff.jar download <reference-genome>
# annotate the vcf with synonymous and non-synonymous sites
java -jar snpEff.jar -v <reference-genome> input.vcf > annotated.vcf

# filter 
grep "synonymous_variant" annotated.vcf > synonymous.txt
grep "missense_variant" annotated.vcf > nonsynonymous.txt

###### SFS calculation ######
# Generate site allele frequency (saf) file
angsd -bam bamlist.txt -doSaf 1 -anc reference.fa -out output -GL 1
# Estimate SFS 
realSFS output.saf.idx -P 4 > output.sfs
  # The output.sfs file contains the estimated SFS for the population

#### To calculate selected (non-synonymous) and neutral (synonymous) SFS:

# For synonymous (neutral) SFS
angsd -bam bamlist.txt -doSaf 1 -anc reference.fa -GL 1 -sites synonymous.txt -out sfs_neutral

# For nonsynonymous (selected) SFS
angsd -bam bamlist.txt -doSaf 1 -anc reference.fa -GL 1 -sites nonsynonymous.txt -out sfs_selected

# Estimate SFS 
realSFS sfs_neutral.saf.idx > sfs_neutral.sfs
realSFS sfs_selected.saf.idx > sfs_selected.sfs

###### Ancestral allele annotation ######
angsd -doAncError 1 -doMaf 1 -anc outgroup.fasta -bam bamlist.txt -out ancestral_out
## Output files are
	#1.	MAF file (*.mafs): This file contains allele frequencies at each site.
	#2.	Ancestral state error file (*.ancError): This file provides information on potential errors in ancestral state inference.
	#3.	SAF file (*.saf): This file stores the Site Allele Frequency likelihoods.


# To create an ancestral allele annotation file for fastDFE using the output files from ANGSD

	# 1.	Extract the ancestral alleles from the *.mafs file:
	  # The MAF file contains the reference and alternative alleles at each site, along with their frequencies.
	  # You can infer the ancestral allele based on which allele matches the outgroup sequence (reference.fa).
	# 2.	Cross-reference with the *.ancError file:
	  #	Use this file to correct potential errors in ancestral state inference.
	# 3.	Format the output:
	  #	The output file should include columns for chromosome, position, reference allele, alternative alleles, and inferred ancestral allele. This can be saved as a text file suitable for fastDFE.
#CHROM   POS     REF     ALT     Ancestral
#chr1    1000    A       G       A
#chr1    1050    C       T       C

#### A Python script that extracts the ancestral alleles from the MAF and ancestral state error files, 
  #### and creates an ancestral allele annotation file suitable for fastDFE:

  import pandas as pd

# Input files
maf_file = "ancestral_out.mafs"
anc_error_file = "ancestral_out.ancError"
output_file = "ancestral_annotation.txt"

# Load the MAF file into a DataFrame
maf_df = pd.read_csv(maf_file, delim_whitespace=True)

# Load the ancestral error file
anc_error_df = pd.read_csv(anc_error_file, delim_whitespace=True, header=None, names=["chrom", "pos", "ref", "anc", "error"])

# Open output file to write the ancestral annotation
with open(output_file, "w") as out_file:
    out_file.write("CHROM\tPOS\tREF\tALT\tAncestral\n")
    
    # Iterate over the MAF DataFrame
    for idx, row in maf_df.iterrows():
        chrom = row['chromo']
        pos = row['position']
        ref_allele = row['ref']
        alt_alleles = row['knownEM'].split(",") if isinstance(row['knownEM'], str) else ['N']  # split alt alleles if they exist
        
        # Find the corresponding row in the ancError DataFrame
        anc_error_row = anc_error_df[(anc_error_df['chrom'] == chrom) & (anc_error_df['pos'] == pos)]
        
        # If ancestral allele is available in the error file and has low error, assign it as ancestral
        if not anc_error_row.empty and anc_error_row['error'].values[0] < 0.5:
            ancestral_allele = anc_error_row['anc'].values[0]
        else:
            ancestral_allele = "N"  # uncertain
        
        # Write the result to output file
        out_file.write(f"{chrom}\t{pos}\t{ref_allele}\t{','.join(alt_alleles)}\t{ancestral_allele}\n")

print("Ancestral allele annotation file created.")

  ## Explanation:
    
	#1.	MAF file parsing: The MAF file is read, and reference and alternative alleles are extracted.
	#2.	Ancestral error file: This file is cross-referenced to check for any errors in the ancestral allele estimation. If the error is low (here set arbitrarily to 0.5), the inferred ancestral allele is used.
	#3.	Output: The ancestral allele annotation is written in a tab-delimited format suitable for fastDFE.


from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

with open("cds.fa") as cds_file, open("protein.fa", "w") as protein_file:
    for record in SeqIO.parse(cds_file, "fasta"):
        # Translate CDS sequence to protein
        protein_seq = record.seq.translate(to_stop=True)
        protein_record = SeqRecord(protein_seq, id=record.id, description="translated protein")
        SeqIO.write(protein_record, protein_file, "fasta")

##################################
 from Bio import SeqIO

# Input files
fasta_file = "reference.fa"  # Your FASTA file with scaffold-based names
gff_file = "annotations.gff"  # Your GFF file with gene IDs
output_fasta = "matched_sequences.fa"  # Output file for matched sequences

# Load FASTA sequences
scaffold_sequences = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

# Function to parse GFF and extract regions
def extract_gene_sequences_from_gff(gff_file, scaffold_sequences):
    gene_to_sequence = {}

    with open(gff_file, "r") as gff:
        for line in gff:
            if not line.startswith("#"):
                fields = line.strip().split("\t")
                chrom = fields[0]  # Chromosome or scaffold name
                feature_type = fields[2]  # Type of the feature (e.g., gene, exon)
                start = int(fields[3]) - 1  # Start position (0-based)
                end = int(fields[4])  # End position
                attributes = fields[8]  # GFF attribute field
                gene_id = ""

                # Extract gene ID from attributes
                for attribute in attributes.split(";"):
                    if attribute.startswith("ID="):
                        gene_id = attribute.split("=")[1]

                # Get sequence from FASTA using scaffold/coordinate
                if chrom in scaffold_sequences and feature_type == "gene":
                    sequence = scaffold_sequences[chrom].seq[start:end]
                    gene_to_sequence[gene_id] = sequence

    return gene_to_sequence

# Extract the gene sequences from GFF and match with scaffolds
matched_sequences = extract_gene_sequences_from_gff(gff_file, scaffold_sequences)

# Write the matched sequences to a FASTA file
with open(output_fasta, "w") as output_handle:
    for gene_id, sequence in matched_sequences.items():
        output_handle.write(f">{gene_id}\n{sequence}\n")

print(f"Matched sequences have been written to {output_fasta}")

#######################
#####################
from Bio import SeqIO

# Input files
fasta_file = "gene_sequences.fa"  # Your FASTA file with gene names
gff_file = "annotations.gff"  # Your GFF file with scaffold coordinates
output_fasta = "scaffold_mapped_sequences.fa"  # Output FASTA with scaffold names

# Dictionary to map gene names to scaffold and coordinates
gene_to_scaffold = {}

# Parse GFF file to map genes to scaffolds
with open(gff_file, "r") as gff:
    for line in gff:
        if not line.startswith("#"):
            fields = line.strip().split("\t")
            chrom = fields[0]  # Scaffold name
            feature_type = fields[2]  # Feature type (e.g., gene, CDS)
            start = fields[3]
            end = fields[4]
            attributes = fields[8]

            # Extract gene ID from GFF attributes
            gene_id = None
            for attr in attributes.split(";"):
                if attr.startswith("ID="):
                    gene_id = attr.split("=")[1]

            # Map gene ID to scaffold and coordinates
            if gene_id and feature_type == "gene":
                gene_to_scaffold[gene_id] = f"{chrom}:{start}-{end}"

# Parse the FASTA file and rename sequences with scaffold information
with open(output_fasta, "w") as out_handle:
    for record in SeqIO.parse(fasta_file, "fasta"):
        gene_id = record.id
        if gene_id in gene_to_scaffold:
            # Rename sequence with scaffold information
            new_id = gene_to_scaffold[gene_id]
            record.id = new_id
            record.description = ""
            SeqIO.write(record, out_handle, "fasta")

print(f"FASTA file with scaffold names written to {output_fasta}")

######################
######## FAST DFE #######

