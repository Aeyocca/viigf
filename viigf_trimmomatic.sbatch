#!/bin/sh -login
#Alan E Yocca
#04-24-18
#run trimmomatic
#06-26-18
#add as	part of	viigf pipeline
#01-15-19
#slurm conversion

#SBATCH --job-name=viigf_trimmomatic
#SBATCH --mem=50Gb
#SBATCH --ntasks=30
#SBATCH --time=03:00:00
#SBATCH --output=/mnt/research/edgerpat_lab/AlanY/Error_files/%x-%j.SLURMout


set -e

source ${CONFIG}

mkdir -p ${WKDIR}
cd ${WKDIR}

#make sure EO exists
mkdir -p ${WKDIR}/Error_files

date

LINE=`/bin/sed -n ${SLURM_ARRAY_TASK_ID}p ${INFILE}`
OUTDIR="04_trimmed"

module load Trimmomatic/0.38-Java-1.8.0_162

mkdir -p ${OUTDIR}

java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.38.jar PE \
-threads 19 \
-phred33 \
01_raw_reads/${LINE}_1.fastq \
01_raw_reads/${LINE}_2.fastq \
${OUTDIR}/${LINE}_1_paired.fastq \
${OUTDIR}/${LINE}_1_unpaired.fastq \
${OUTDIR}/${LINE}_2_paired.fastq \
${OUTDIR}/${LINE}_2_unpaired.fastq \
ILLUMINACLIP:/mnt/research/common-data/Bio/Trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 TRAILING:30 MINLEN:41

date

#fastqc on raw reads
#qsub -t ${PBS_ARRAYID} -v INFILE=${INFILE},CONFIG=${CONFIG} /mnt/research/edgerpat_lab/AlanY/01_athal_cns/qsub/viigf_fastqc_raw.sh
#sbatch --array=${SLURM_ARRAY_TASK_ID} \
#--export=INFILE=${INFILE},CONFIG=${CONFIG} \
#${SCRIPT_DIR}/viigf_fastqc_raw.sbatch


#fastqc on cleaned reads
#qsub -t ${PBS_ARRAYID} -v INFILE=${INFILE},CONFIG=${CONFIG} /mnt/research/edgerpat_lab/AlanY/01_athal_cns/qsub/viigf_fastqc_cleaned.sh
#sbatch --array=${SLURM_ARRAY_TASK_ID} \
#--export=INFILE=${INFILE},CONFIG=${CONFIG} \
#${SCRIPT_DIR}/viigf_fastqc_cleaned.sbatch

#submit var_inc_pipe
if [[ ${CONTINUE} =~ ^[Yy] ]]; then
	#qsub -t ${PBS_ARRAYID} -v CONFIG="/mnt/research/egerpat_lab/AlanY/01_athal_cns/qsub/viigf_config.conf" -l walltime=${VAR_INC_WALLTIME},nodes=40:ppn=1,mem=40gb /mnt/research/edgerpat_lab/AlanY/01_athal_cns/qsub/viigf_var_inc.sh
	sbatch --array=${SLURM_ARRAY_TASK_ID} \
	--export=CONFIG=${CONFIG} \
	--time=${VAR_INC_WALLTIME} \
	--ntasks=40 \
	--mem=40Gb \
	--output=${WKDIR}/Error_files/%x-%A_%a.SLURMout \
	${SCRIPT_DIR}/viigf_var_inc.sbatch
fi

#${TRIM_OUT_STAT} ${PBS_JOBID}

#qstat -f ${PBS_JOBID}
echo "Finished"; date
