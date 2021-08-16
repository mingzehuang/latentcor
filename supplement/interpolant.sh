#!/bin/bash
##ENVIRONMENT SETTINGS; CHANGE WITH CAUTION
#SBATCH --export=NONE                #Do not propagate environment
#SBATCH --get-user-env=L             #Replicate login environment
#SBATCH --partition=bigmem
##NECESSARY JOB SPECIFICATIONS
#SBATCH --job-name=interpolation      #Set the job name to
#SBATCH --time=01:00:00               #Set the wall clock limit to 6hr and 30min
#SBATCH --nodes=1                    #Request 1 node
#SBATCH --ntasks-per-node=80         #Request 8 tasks/cores per node
#SBATCH --mem=320GB                     #Request 320GB per node
#SBATCH --output=interpolant_Out.%j            #Send stdout/err to "Example2Out.[jobID]"

##OPTIONAL JOB SPECIFICATIONS
#SBATCH --mail-type=ALL                     #Send email on all job events
#SBATCH --mail-user=sharkmanhmz@tamu.edu    #Send all emails

#First Executable Line
module load iccifort/2020.4.304  impi/2019.9.304 R/4.1.0
R CMD BATCH --no-save --no-restore --slave interpolant.R interpolant.Rout

