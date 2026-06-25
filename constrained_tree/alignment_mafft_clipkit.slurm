#!/bin/bash
####### Reserve computing resources #############
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=10:00:00
#SBATCH --mem=10G

echo "starting run at: `date`"

####### Set environment variables ###############
source ~/software/init-conda
PATH=/work/soghigian_lab/apps/AMAS/amas:$PATH

cd /path/to/working/directory

####### Script ###################################
# Alignment
#conda activate mafft
mafft --reorder 18S.fasta > 18S_aligned.fasta
mafft --reorder 28S.fasta > 28s_aligned.fasta
mafft --reorder ak.fasta > ak_aligned.fasta
mafft --reorder coi.fasta > coi_aligned.fasta
mafft --reorder coii.fasta > coii_aligned.fasta
mafft --reorder eno.fasta > eno_aligned.fasta
einsi --reorder its2.fasta > its2_aligned.fasta

# Alignment was visually inspected and minor adjustments were made in Geneious

# Standardized sequence headers
cut -d ':' -f 1 18S_aligned.fasta > 18S.out.fasta
cut -d ':' -f 1 28S_aligned.fasta > 28S.out.fasta
cut -d ':' -f 1 ak_aligned.fasta > ak.out.fasta
cut -d ':' -f 1 coi_aligned.fasta > coi.out.fasta
cut -d ':' -f 1 coii_aligned.fasta > coii.out.fasta
cut -d ':' -f 1 eno_aligned.fasta > eno.out.fasta
cut -d ':' -f 1 its2_aligned.fasta > its2.forClipkit.fasta

conda deactivate
conda activate /work/soghigian_lab/apps/conda/envs/clipkit

# Recommended to run multiple trimming strategies in clipkit and examine resulting Average Bipartition Support (ABS) for trees
# Default gappyness is 90%
# Very high gappyness thresholds are recommended (e.g. retain sites with <95% gaps), -g 0.95
# Low thresholds may remove too many sites and worsen phylogenetic inference

clipkit -m kpic-gappy -l its2.forClipkit.fasta -o its2.out.fasta

# Concatenate genes
python3 /work/soghigian_lab/apps/AMAS/amas/AMAS.py concat -f fasta -d dna -i *.out.fasta -u fasta
mv concatenated.out concatenated.fasta

echo "Job finised with exit code $? at: `date`"
