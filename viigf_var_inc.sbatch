#!/bin/sh -login
#Alan E Yocca
#05-17-18
#Variant incorporation pipe, script should take about 4 hours to complete each round, try and get it to run 2 rounds overnight
#06-27-18 added to viigf pipeline
#01-15-19 slurm conversion

#SBATCH --job-name=viigf_var_inc
#SBATCH --mem=4Gb
#SBATCH --ntasks=1
#SBATCH --time=04:00:00
#SBATCH --output=/mnt/research/edgerpat_lab/AlanY/Error_files/%x-%j.SLURMout

#set -e

source ${CONFIG}

mkdir -p ${WKDIR}
cd ${WKDIR}

#Make sure EO exists:
mkdir -p ${EO_DIR}

echo ${WKDIR}

date

#mass module load
#module load bcftools
#module load GNU/6.4.0-2.28  OpenMPI/2.1.1
#module load BWA/0.7.17
#module load SAMtools/1.7
#module load picard/2.18.1-Java-1.8.0_152
#module load GATK/4.0.5.1-Python-3.6.4
#module load BamTools/2.5.1

#Save bwa for later because needs different version of openmpi


##************TEST************
#echo "Modules loaded"
#echo $?
#echo ""

#bwa_mem

#LINE=`/bin/sed -n ${PBS_ARRAYID}p ${INFILE}`
#echo ${LINE}
#echo ${START}
#echo ${STOP}

#define i, keeping that as my curent iteration variable because it is everywhere and really don't want to have to replace it, I think it is fine
#if i is not defined (ie this script was submitted via the trimmomatic script or by command line and want to start at the first iteration),
##set i to first iteration
#otherwise, i will have been defined by submitting it at the bottom of this script from the previous iteration 
if [ -z ${I+x} ]; then
	I=${VAR_START}
	echo "i unset, setting to var start"
else
	echo "i is set to ${I}"
fi

echo "Loop set to: ${I}"

##************TEST************
#echo "Loop variable checked"
#echo $?
#echo ""

PAST_REF=$(expr ${I} - 1)
#CONSENSUS_PATH="08_consensus/${LINE}_bwa_alt_${PAST_REF}.fasta"

#if run 1, used TAIR10 as reference
#if bwa index already exists for reference, don't make another, otherwise, index it
INDEXED=false

#Define outside loop
CONSENSUS_PATH=""
if [ ${I} == 1 ]; then
        CONSENSUS_PATH=${REF_GENOME}
else
    	CONSENSUS_PATH="08_consensus/${LINE}_bwa_alt_${PAST_REF}.fasta"
fi

#samtools to index
module load SAMtools/1.7
samtools faidx ${CONSENSUS_PATH}

module load GNU/6.4.0-2.28  OpenMPI/2.1.1
module load BWA/0.7.17

if [ ${I} == 1 ]; then
	CONSENSUS_PATH=${REF_GENOME}
	#check if bwa index exists
	BWA_INDEX=$(echo ${CONSENSUS_PATH} | sed "s/\.fasta/_bwa_index/")
	bwa_index_check=$(ls ${BWA_INDEX}* 2> /dev/null | wc -l)
	if [ ${bwa_index_check} != 0 ]; then
		echo "BWA index for ${CONSENSUS_PATH} exists, not making"
		INDEXED=true
	else
		echo "BWA index for ${CONSENSUS_PATH} does not exist, creating"
		bwa index -p ${BWA_INDEX} ${CONSENSUS_PATH}
		INDEXED=true
	fi	
fi

TAG="${LINE}_bwa_alt_${I}"
CONSENSUS_PATH_BASE=$(echo $CONSENSUS_PATH | sed 's/\.fasta$//')
BWA_INDEX=$(echo ${CONSENSUS_PATH} | sed "s/\.fasta/_bwa_index/")
TAG_DIR=$(dirname ${TAG})

echo "started on ${LINE}"
echo "past_ref: ${PAST_REF}"
echo "consensus path: ${CONSENSUS_PATH}"
echo "tag: ${TAG}"
echo "tag dir: ${TAG_DIR}"

#module load bwa

####*********BWA INDEXING********

#if first round, need to index TAIR10, but also need to check first if TAIR10 has already been indexed

if ${INDEXED};
then
	echo "${BWA_INDEX} files for reference exist, not indexing again"
	echo "if these files are corrupt, remove and rerun this script or on command line"
else
	echo "creating bwa index for non-reference genome: ${CONSENSUS_PATH}"
	bwa index -p ${BWA_INDEX} ${CONSENSUS_PATH}
fi

mkdir -p 06_SAM/${TAG_DIR}
#rm -f 06_SAM/${TAG}.sam

#Check if output already exists / if force is specified
#Done throughout

file_check=$(ls 06_SAM/${TAG}.sam 2> /dev/null | wc -l )

if [ ${file_check} != 0 ] \
&& [ ${FORCE} == "false" ]; then
	echo "BWA output exists (${file_check}), and FORCE set to false, NOT running BWA MEM"
else
	echo "Either bwa output doesnt exist or force not set to false (${FORCE})"
	echo "Running BWA MEM"; echo ""
	bwa mem -R'@RG\tID:group1\tSM:sample1\tPL:illumina\tLB:lib1\tPU:unit1' \
	-t 39 \
	${BWA_INDEX} \
	${READ_DIR}/${LINE}_1_paired.fastq \
	${READ_DIR}/${LINE}_2_paired.fastq \
	> 06_SAM/${TAG}.sam
	echo ""; echo "Finished BWA MEM"; date; echo "";
fi


#mass module load
module load OpenMPI/2.1.2
module load bcftools
#module load GNU/6.4.0-2.28  OpenMPI/2.1.1
#module load BWA/0.7.17
module load SAMtools/1.7
module load picard/2.18.1-Java-1.8.0_152
#module load GATK/4.0.5.1-Python-3.6.4
#Keeping version 3 because thats how others were assembled:
module load GATK/3.8-1-0-gf15c1c3ef-Java-1.8.0_112
module load BamTools/2.5.1


####################### BAM Processing ##############################

mkdir -p 07_BAM/${TAG_DIR}

file_check=$(ls 07_BAM/${TAG}.bam 2> /dev/null | wc -l )

if [ ${file_check} != 0 ] \
&& [ ${FORCE} == "false" ]; then
	echo "Sorted BAM exists (${file_check}), and FORCE set to false, NOT running Picard SortSam"
else
	echo "Either Sorted BAM doesnt exist, or force not set to false (${FORCE})"
	echo "Running Picard SortSam"; echo ""
	java -jar $EBROOTPICARD/picard.jar \
	SortSam \
	INPUT=06_SAM/${TAG}.sam \
	OUTPUT=07_BAM/${TAG}.bam \
	SORT_ORDER=coordinate
	echo ""; echo "Finished Picard SortSam"; date; echo ""
fi

#index
samtools index 07_BAM/${TAG}.bam

#mark duplicates
file_check=$(ls 07_BAM/${TAG}_markdup.bam 2> /dev/null | wc -l )

if [ ${file_check} != 0 ] \
&& [ ${FORCE} == "false" ]; then
	echo "BAM Markdup exists (${file_check}), and FORCE set to false, NOT running Picard MarkDuplicates"
else
	echo "Either Markdup BAM doesnt exist, or force not set to false (${FORCE})"
	echo "Running Picard MarkDuplicates"; echo ""
	java -jar $EBROOTPICARD/picard.jar \
	MarkDuplicates \
	M=07_BAM/${TAG}_dup_stats.txt \
	I=07_BAM/${TAG}.bam \
	REMOVE_DUPLICATES=true \
	O=07_BAM/${TAG}_markdup.bam
	echo ""; echo "Finished Picard MarkDuplicates"; date; echo ""
fi


####################### Variant Calling ##############################

#index
samtools index 07_BAM/${TAG}_markdup.bam

CONSENSUS_PATH_BASE_DIR=$(dirname ${CONSENSUS_PATH_BASE})

file_check=$(ls ${CONSENSUS_PATH_BASE}.dict 2> /dev/null | wc -l)

if [ ${file_check} == 1 ]; then
	echo "PicardTools Sequence Dictionary exists already, not creating a new one."
	echo "If you suspect the existing dictionary is corrupt, delete this file and run this script again"
	echo "${CONSENSUS_PATH_BASE}.dict"
else
	echo "making sequence dictionary"
	java -jar $EBROOTPICARD/picard.jar \
	CreateSequenceDictionary \
	R= ${CONSENSUS_PATH} \
	O= ${CONSENSUS_PATH_BASE}.dict
fi

#indelrealigner
file_check=$(ls 07_BAM/${TAG}.intervals 2> /dev/null | wc -l )

if [ ${file_check} != 0 ] \
&& [ ${FORCE} == "false" ]; then
	echo "IndelRealigner intervals file exists (${file_check}), and FORCE set to false, NOT running GATK RealignerTargetCreator"
else
	echo "Either IndelRealigner intervals file doesnt exist, or force not set to false (${FORCE})"
	echo "Running GATK RealignerTargetCreator"; echo ""
	java -jar $EBROOTGATK/GenomeAnalysisTK.jar \
	-T RealignerTargetCreator \
	-R ${CONSENSUS_PATH} \
	-I 07_BAM/${TAG}_markdup.bam \
	-o 07_BAM/${TAG}.intervals
	echo ""; echo "Finished running GATK RealignerTargetCreator"; date; echo ""
fi

file_check=$(ls 07_BAM/${TAG}_realign.bam 2> /dev/null | wc -l )

if [ ${file_check} != 0 ] \
&& [ ${FORCE} == "false" ]; then
	echo "IndelRealigner output file exists (${file_check}), and FORCE set to false, NOT running GATK IndelRealigner"
else
	echo "Either InderRealigner output file doesnt exist, or force not set to false (${FORCE})"
	echo "Running GATK IndelRealigner"; echo ""
	java -jar $EBROOTGATK/GenomeAnalysisTK.jar \
	-T IndelRealigner \
	-R ${CONSENSUS_PATH} \
	-I 07_BAM/${TAG}_markdup.bam \
	-targetIntervals 07_BAM/${TAG}.intervals \
	-o 07_BAM/${TAG}_realign.bam
	echo ""; echo "Finished running GATK IndelRealigner"; date; echo ""
fi


samtools index 07_BAM/${TAG}_realign.bam

#module load GATK

mkdir ${LINE}_tmp

samtools faidx ${CONSENSUS_PATH}

mkdir -p 12_VCF/${TAG_DIR}

file_check=$(ls 12_VCF/${TAG}.vcf 2> /dev/null | wc -l )

if [ ${file_check} != 0 ] \
&& [ ${FORCE} == "false" ]; then
	echo "UnifiedGenotyper output file exists (${file_check}), and FORCE set to false, NOT running GATK IndelRealigner"
else
	echo "Either UnifiedGenotyper output file doesnt exist, or force not set to false (${FORCE})"
	echo "Running GATK UnifiedGenotyper"; echo ""
	java -Xmx10g -Djava.io.tmpdir=./${LINE}_tmp \
	-jar $EBROOTGATK/GenomeAnalysisTK.jar \
	-R ${CONSENSUS_PATH} \
	-T UnifiedGenotyper \
	-I 07_BAM/${TAG}_realign.bam \
	--out 12_VCF/${TAG}.vcf \
	--num_threads 39 \
	--genotype_likelihoods_model BOTH \
	--sample_ploidy 2 \
	--defaultBaseQualities 30
	echo ""; echo "Finished running GATK UnifiedGenotyper"; date; echo ""
fi

rm -rf ${LINE}_tmp

#vcfutils varFilter
file_check=$(ls 12_VCF/${TAG}_vf.vcf 2> /dev/null | wc -l )

if [ ${file_check} != 0 ] \
&& [ ${FORCE} == "false" ]; then
	echo "Varfilter output file exists (${file_check}), and FORCE set to false, NOT running VarFilter"
else
	echo "Either VarFilter output files doesnt exist, or force not set to false (${FORCE})"
	echo "Running VarFilter"; echo ""
	rm -f 12_VCF/${TAG}_vf.vcf
	
	${SCRIPT_DIR}/vcfutils.pl \
	varFilter 12_VCF/${TAG}.vcf \
	> 12_VCF/${TAG}_vf.vcf
	
	echo ""; echo "Finished runing VarFilter"; date; echo ""
fi

mkdir -p 08_consensus/${TAG_DIR}


file_check=$(ls 08_consensus/${TAG}.fasta 2> /dev/null | wc -l )

if [ ${file_check} != 0 ] \
&& [ ${FORCE} == "false" ]; then
	echo "FARM output file exists (${file_check}), and FORCE set to false, NOT running GATK FastaAlternativeReferenceMaker"
else
	echo "Either FARM output file doesnt exist, or force not set to false (${FORCE})"
	echo "Running GATK FastaAlternativeReferenceMaker"; echo ""

	java -jar $EBROOTGATK/GenomeAnalysisTK.jar \
	-T FastaAlternateReferenceMaker \
	-R ${CONSENSUS_PATH} \
	-o 08_consensus/${TAG}.fasta \
	-V 12_VCF/${TAG}_vf.vcf

	echo ""; echo "Finished running GATK FastaAlternativeReferenceMaker"; date; echo ""
fi

	
#sed -i "s/1 Chr1:1/Chr1/g" 08_consensus/${TAG}.fasta
#sed -i "s/2 Chr2:1/Chr2/g" 08_consensus/${TAG}.fasta
#sed -i "s/3 Chr3:1/Chr3/g" 08_consensus/${TAG}.fasta
#sed -i "s/4 Chr4:1/Chr4/g" 08_consensus/${TAG}.fasta
#sed -i "s/5 Chr5:1/Chr5/g" 08_consensus/${TAG}.fasta
#sed -i "s/6 mitochondria:1/mitochondria/g" 08_consensus/${TAG}.fasta
#sed -i "s/7 chloroplast:1/chloroplast/g" 08_consensus/${TAG}.fasta

LAST_VAR_INC="08_consensus/${TAG}.fasta"

#DO NOT	CONTINUE IF CONSENSUS NOT MADE SUCCESSFULLY
consensus_check=$(ls 08_consensus/${TAG}.fasta 2> /dev/null | wc -l)
if [ ${consensus_check} == 0 ]; then
	echo "Round ${I} Failed, consensus not created, exiting"
	exit 1
fi

date

#submit the next round
#check for what walltime to use:
if [[ ${MANUAL_WALLTIME} =~ ^[nN] ]]; then
        READ_FILE_SIZE=$(stat ${READ_DIR}/${LINE}_1_paired.fastq | grep "Size:" | cut -f4 -d" " | cut -f1 -d"	")
        if [[ ${READ_FILE_SIZE} -gt 12000000000 ]]; then
                #more than 4 hours for var_inc
                VAR_INC_WALLTIME=12:00:00
                if [[ ${READ_FILE_SIZE} -gt 35000000000 ]]; then
                        PBJ_WALLTIME=6:00:00
                fi
        else
                VAR_INC_WALLTIME=4:00:00
                PBJ_WALLTIME=4:00:00
        fi
else
    	echo "using walltime specified in configuration file"
fi

if [ ${I} -lt ${VAR_LAST} ]; then
	#increase for next round
	I=$(expr ${I} + 1)

	#qsub -t ${PBS_ARRAYID} -v CONFIG=${CONFIG},I=${I} -l walltime=${VAR_INC_WALLTIME},nodes=40:ppn=1,mem=40gb /mnt/research/edgerpat_lab/AlanY/01_athal_cns/qsub/viigf_var_inc.sh
	sbatch --array=${SLURM_ARRAY_TASK_ID} \
	--export=CONFIG=${CONFIG},I=${I} \
	--time=${VAR_INC_WALLTIME} \
	--ntasks=40 \
	--mem=40Gb \
	--output=${WKDIR}/Error_files/%x-%A_%a.SLURMout \
	${SCRIPT_DIR}/viigf_var_inc.sbatch

	echo "submitted next round var_inc with i=${I}"
	echo "submitted within the var_inc script"
	echo "slrum array: ${SLURM_ARRAY_TASK_ID}"
	echo "line: ${LINE}"
elif [[ ${CONTINUE} =~ ^[Yy] ]]
then
	#submit PBJelly
	#if this was just the last iteration, submit pbj
	#qsub -t ${PBS_ARRAYID} -v CONFIG=${CONFIG},CONSENSUS_PATH=${LAST_VAR_INC} -l walltime=${PBJ_WALLTIME},nodes=50:ppn=1,mem=40gb /mnt/research/edgerpat_lab/AlanY/01_athal_cns/qsub/viigf_pbj.sh
	sbatch --array=${SLURM_ARRAY_TASK_ID} \
	--export=CONFIG=${CONFIG},CONSENSUS_PATH=${LAST_VAR_INC} \
	--time=${PBJ_WALLTIME} \
	--ntasks=50 \
	--mem=40Gb \
	${SCRIPT_DIR}/viigf_viigf_pbj.sbatch

	echo "submitted pbj"
	echo "slurm array: ${SLURM_ARRAY_TASK_ID}"
        echo "line: ${LINE}"
	echo "with conensus path: ${LAST_VAR_INC}"
else
	echo "NOT continuing on to PBJ"
	echo "if you would like to continue, start up from PBJ"
	echo "if you believe it should have continued, make sure CONTINUE variable in configuration file is set to Yes (or anything that matches this regex: =~ ^[Yy])"
	if [[ ${CONFIG} =~ *launch_edited.conf$ ]]
	then
		rm -f ${CONFIG}
	fi
fi

#run the cleaner
if [[ ${RUN_CLEANERS} =~ ^[Yy] ]]; then
	echo "running the cleaner"
	source ${CLEANER_VAR_INC}
else
	echo "NOT running cleaner"
fi

#qstat -f ${PBS_JOBID}
