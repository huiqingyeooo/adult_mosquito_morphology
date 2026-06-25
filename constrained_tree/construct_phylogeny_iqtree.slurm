#!/bin/bash
####### Reserve computing resources #############
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --time=24:00:00
#SBATCH --mem=100G

####### Set environment variables ###############
export PATH=/work/soghigian_lab/apps/iqtree-latest/bin:$PATH
module load perl/5.32.0
source ~/software/init-conda

directory=/path/to/working/directory
inputFasta=concatenated.fasta
constraintTree=S2017.S2023.genbank.monophyly.treefile
paritionFile=partitions.txt
prefixName=S2017S2023genbank_monophyly


echo "starting run at: `date`"

####### Script ##################################

cd ${directory}

# Running iqtree
iqtree2 -s ${inputFasta} -bb 1000 -g ${constraintTree} \
-p ${paritionFile} -m TEST -T 12 --prefix ${prefixName} -keep_empty_seq

echo "Job finised with exit code $? at: `date`"
