#!/bin/sh

#SBATCH --partition=t2standard
#SBATCH --nodes=16
#SBATCH --mail-user=memuscarella@alaska.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --output=mothur-batch.%j

ulimit -s unlimited
ulimit -l unlimited

# Load any desired modules, usually the same as loaded to compile
source /center1/MICROBES/Mothur/mothur-1.46.1/setEnvironment.sh

cd /center1/MICROBES/memuscarella/MicroEuks/analysis/

# Run Mothur Batch
mothur tadpole_MicroEuks.batch 