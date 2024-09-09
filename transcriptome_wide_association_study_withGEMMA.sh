## Prepare bimbam file for TWAS
# Make kinship
gemma -g bimbamfile.txt -p trait.txt -gk 2 -o kinshipname 
  # gk 2 is for making kinship with bimbam input
# Run TWAS
gemma -g bimbamfile.txt -p traits.txt -k kinshipname.sXX.txt -lmm 4 -o twas_outputname
