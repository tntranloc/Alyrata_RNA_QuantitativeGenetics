## Prepare bimbam file for TWAS
  # first column: gene id
  # second column: a nucleotide (A f.e.)
  # third column: a corresponding nucleotide (T f.e.)
  # fourth column onwards: sample names
  # save the data with write.table(data, "dataname.txt", sep = ",", quote = F, row.names = F, col.names = F)

  
# Make kinship
gemma -g bimbamfile.txt -p trait.txt -gk 2 -o kinshipname 
  # gk 2 is for making kinship with bimbam input
# Run TWAS
gemma -g bimbamfile.txt -p traits.txt -k kinshipname.sXX.txt -lmm 4 -o twas_outputname
