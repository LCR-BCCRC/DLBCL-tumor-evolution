#! /bin/bash

if [[ $CONDA_PREFIX != "/projects/rmorin_scratch/conda_environments/GAMBLR" ]]; then
    echo "ERROR: You must have the GAMBLR conda env loaded to run this script. "
    echo "Your current conda env is $CONDA_PREFIX. "
    echo "Run conda activate /projects/rmorin_scratch/conda_environments/GAMBLR first. "
    exit 1
fi

# Generates a stable table of trios metadata from the GAMBL metadata. 
# Should not be re-run regularly. 
# Rscript src/create_sequencing_md.R

# Generates the clinical metadata. 
# Depends on files that have been extensively manually curated. 
Rscript src/create_clinical_md.R

# Create the Swimmer plot for figure 1. 
Rscript src/swimmer_plot.R