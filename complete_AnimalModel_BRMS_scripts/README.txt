# To run animal model using R brms package
# check dependencies ins script
# Written to run on a HPC cluster (Slurm)  
# ..._script_generator will generate jobs in chunk, 200 genes each (this can be changed in the script)
# ...._wrapped_function is the masterfile, it takes each gene as a phenotype to fit the model, output individual csv file, as well as a merged csv file of 200 genes when the job is done

# "animal_model_5phenotypes_....": modified the wrapped function to run locally for small number of phenotypes 


