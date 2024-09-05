

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
######## FAST DFE #######

