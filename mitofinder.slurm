#!/bin/bash
####### Reserve computing resources #############
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --mem=50G

echo "starting run at: `date`"

####### Set environment variables ###############
PATH=/work/soghigian_lab/apps/MitoFinder/:$PATH

referenceDB=/path/to/mitofinder/culicidae_mitoGenome.gb
fasta=/path/to/fasta/files
output=/path/to/mitofinder/output

#for file in $(ls ${fasta}/*.fasta);
#do
#taxa=$(basename $file .fasta)
#sed "s/TEMPLATE/${taxa}/g" mitofinder_fasta_TEMPLATE.slurm > mitofinder_fasta_${taxa}.slurm
#sbatch mitofinder_fasta_${taxa}.slurm
#done

###### Script ###################################
cd /scratch/${SLURM_JOB_ID}
echo results in ${SLURM_JOB_ID}

# Mitofinder
# -o: genetic code 5 for the invertebrate mt code
# -j: name of the process
# -r: reference mt genome in genbank format
# -m: max memory
# --min-contig-size: minimum size of contig to be considered. default=1000
# --blast-size: % of overlap in best blast hit to be retained. default=30
# --blast-identity-nucl: Nucleotide identity percentage for a hit to be retained. Default = 50

mitofinder --megahit -j TEMPLATE -r ${referenceDB} \
-a ${fasta}/TEMPLATE.fasta \
--min-contig-size 500 --blast-size 20 -o 5 -p 8 -m 50

cp /scratch/${SLURM_JOB_ID}/*.log ${output}
cp /scratch/${SLURM_JOB_ID}/TEMPLATE/TEMPLATE_*Results/*_NT.fasta ${output}

rm /work/soghigian_lab/huiqing.yeo/aedes_coi/mitofinder/mitofinder_fasta_TEMPLATE.slurm

echo "Job finised with exit code $? at: `date`"
