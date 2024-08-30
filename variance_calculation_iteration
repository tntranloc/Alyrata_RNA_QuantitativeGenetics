## This script runs in chunks, given by start and end in the function
## Last modified 30 August 2024
## Nhu LT Tran
## on RNA count data

## Data strcuture
    
    ##### Main data matrix #####
    ## rownames are sample name
    ## colnames are gene name
    ## values are normalised read count (I used DESeq2 median-to-mean ratios method)
   
    ##### Random effects #####
    ## two columns: one is "animal" with sample ID, one is name of the random effect
    ## f.e.: batch_info has column "animal" and "batch"
    ## make sure it's "factor" when needed

## Package required 
    ## brms to run the model
    ## nadiv to create Additive and Dominance matrices before hand
        ## Note that Amat can be saved as txt by write.table() in R
        ## but Dmat must be written by writeMM() from R's Matrix package, and read by readMM function
        #### SEE the preparing_animal_model_input.R script regarding how to create the matrices  

## This script can be ran on cluster 
## It can skip existing file, and output individual files for each gene (memory efficient)

#####################################
####### Beginning of R script #######

## Create a temporary directory ## do NOT skip this!
temp_dir = “/path/to/your/main/directory/cheopsusername/tmp”
Sys.setenv(TMPDIR = temp_dir)

## Load libraries 
library(brms)
library(Matrix) # to load Dominance matrix

setwd(“/path/to/your/working/directory/”)
data = read.csv(“/path/to/your/directory/your_count_table.csv”, header= T)
## cleaning up the matrix, as rownames become first column "X" when loading in
rownames(data) = data$X
data$X= as.numeric(data$X)
# example script for extracting intra-population, with ID less than 400
data_intra = subset(data, data$X < 400)
data_intra = data_intra[,-1] # remove "X" column

# Load additive and dominance matrices
Amat_intra = read.csv(“/path/to/your/directory/your_Amatrix.csv”, header = T)
Dmat_intra = readMM(“/path/to/your/directory/your_Dmatrix.csv”)

# same cleaning up as earlier
a = Amat_intra$X 
Amat_intra = Amat_intra[,-1]
colnames(Amat_intra) = a
rownames(Amat_intra) = a
rm(a)
colnames(Dmat_intra) = colnames(Amat_intra)
rownames(Dmat_intra) = rownames(Amat_intra)

# Load random effects
    # here I have so called sowing batch, dam (maternal), and tray ID as random effects
batch_info_intra = read.csv(“/path/to/your/directory/your_batch_info.csv”, header = T)
dam_info_intra = read.csv(“/path/to/your/directory/your_dam_info.csv”, header = T)
tray_info_intra = read.csv(“/path/to/your/directory/your_tray_info.csv”, header = T)

# renaming columns 
colnames(batch_info_intra)[1] = “animal”
colnames(dam_info_intra)[1] = “animal”
colnames(tray_info_intra)[1] = “animal”
print(head(batch_info_intra))
print(head(dam_info_intra))


# set the global Amat and Dmat variables which are used in the function below
    # if not doing this, remember to provide Amat and Dmat args in the function below
Amat = Amat_intra
Dmat = Dmat_intra


## The function ##
variancesCalc = function(myfile, st, end, batch_info, dam_info, tray_info, outdir, individual_files = FALSE) {
  out1 = list()
  counts = myfile[, st:end]

  for (iter in 1:length(colnames(counts))) {
    temp = as.data.frame(counts[, iter])
    gene = colnames(counts)[iter]
    output_gene_file = paste0(outdir, “/output_intra_“, gene, “.csv”)

    # Skip processing if the output file already exists
    if (individual_files && file.exists(output_gene_file)) {
      next
    }

    colnames(temp) = c(“counts”)
    temp$animal = as.character(rownames(counts))
    temp$dom = as.character(rownames(counts))
    temp$counts = as.numeric(temp$counts)
    temp$batch = batch_info$batch[match(temp$animal, batch_info$animal)]
    temp$batch = factor(temp$batch)
    temp$dam = dam_info$dam[match(temp$animal, dam_info$animal)]
    temp$dam = factor(temp$dam)
    temp$tray = tray_info$tray[match(temp$animal, tray_info$animal)]
    temp$tray = factor(temp$tray)
    print(head(temp))

    myformula = brmsformula(counts ~ 1 + (1 | gr(animal, cov = Amat)) + (1 | gr(dom, cov = Dmat)) + (1 | dam) + (1 | batch)) + (1 | tray))
    fit = brm(myformula,
              data = temp, family = skew_normal(),
              data2 = list(Amat = Amat, Dmat = Dmat),
              chains = 4, threads = threading(4), cores = 4, iter = 4000, warmup = 1000, thin = 2
    )
    Va = median(unlist(VarCorr(fit, summary = FALSE)$animal))
    Vd = median(unlist(VarCorr(fit, summary = FALSE)$dom))
    Vr = median(unlist(VarCorr(fit, summary = FALSE)$residual))
    Vm = median(unlist(VarCorr(fit, summary = FALSE)$dam))
    Vb = median(unlist(VarCorr(fit, summary = FALSE)$batch))
    ha = round(Va / (Va + Vm + Vr + Vd + Vb), digits = 2)
    hd = round(Vd / (Va + Vm + Vr + Vd + Vb), digits = 2)
    hr = round(Vr / (Va + Vm + Vr + Vd + Vb), digits = 2)
    hm = round(Vm / (Va + Vm + Vr + Vd + Vb), digits = 2)
    hb = round(Vb / (Va + Vm + Vr + Vd + Vb), digits = 2)
    sdA = sd(unlist(VarCorr(fit, summary = FALSE)$animal))
    sdD = sd(unlist(VarCorr(fit, summary = FALSE)$dom))
    sdR = sd(unlist(VarCorr(fit, summary = FALSE)$residual))
    sdM = sd(unlist(VarCorr(fit, summary = FALSE)$dam))
    sdB = sd(unlist(VarCorr(fit, summary = FALSE)$batch))
    RA = summary(fit)$random$animal$Rhat
    RD = summary(fit)$random$dom$Rhat
    RM = summary(fit)$random$dam$Rhat
    RR = summary(fit)$spec_pars$Rhat
    RB = summary(fit)$random$batch$Rhat

    outp = c(as.character(gene), Va, Vd, Vm, Vr, Vb, ha, hd, hm, hr, hb, sdA, sdD, sdM, sdR, sdB, RA, RD, RM, RR, RB)
    out1[[iter]] = outp

    if (individual_files) {
      individual_output = data.frame(t(outp))
      write.table(individual_output, file = output_gene_file, sep = ”	“, row.names = FALSE, col.names = TRUE)
    }
  }
  out1_df = do.call(rbind, out1)
  write.table(out1_df, file = output_file, sep = ”	“, row.names = FALSE, col.names = TRUE)

  return(out1_df)
}


## Apply the function
variancesCalc(data_intra, 801, 1000, batch_info_intra, dam_info_intra, tray_info_intra, “/path/to/your/output/directory”, individual_files = TRUE)
    ## The function takes data matrix (data_intra), processing column number 801 to 1000 (genes), taking batch info, dam info, and tray info as random effects
    ## This will create 200 csv files, for gene index 801 and 1000 respectively. 

## Simply merge the csv for downstream analysis


