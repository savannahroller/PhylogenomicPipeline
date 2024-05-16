#!/bin/bash
#SBATCH --job-name=phylogenomicpipeline
#SBATCH --ntasks-per-node=12
#SBATCH --time=24:0:0
#SBATCH --output=phylogenomicpipeline.out
#SBATCH --error=phylogenomicpipeline.err
#SBATCH --mail-user=rollers@oregonstate.edu
#SBATCH --mail-type=END

#run python code 
python3 phylogenomic_pipeline.py 