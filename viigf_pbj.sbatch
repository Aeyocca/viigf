a#!/bin/sh -login
#Alan E Yocca
#06-08-18
#gap_fill_pbj_test
#iterative gap fill using pbjelly
#06-28-18
#added to viigf pipeline

#! /usr/bin/env bash
#PBS -N viigf_pbj
#PBS -j oe
#PBS -o /mnt/research/edgerpat_lab/AlanY/Error_files
#PBS -e /mnt/research/edgerpat_lab/AlanY/Error_files
#PBS -m n

#set -e

source ${CONFIG}

mkdir -p ${WKDIR}
cd ${WKDIR}

#Make sure EO exists
mkdir -p ${EO_DIR}

date

#mass module load
module load SAMTools
module load BAMTools
module load bcftools
module load picardTools
module load bwa
module load GATK
module load BEDTools
module load seqtk/1.0

#start with a consensus from variant incorporation. 
#map all reads to it, 
#filter bam file,
#collect unmapped reads, de novo assemble these unmapped reads, 
#use assembled unmapped reads to gapfill using jellyfish

#again here, ${I} used as loop variable, fed in from config file
##set i to first iteration
#otherwise, i will have been defined by submitting it at the bottom of this script from the previous iteration 
if [ -z ${I+x} ]; then
	I=${PBJ_START}
	echo "i unset, setting to ${PBJ_START}"
else
	echo "Starting round ${I}!!!"
fi

echo "Checked for existence of i"
echo ""
echo "STARTED ROUND: ${I}"
#I think I got rid of this variable
#LAST_GF=$( expr ${I} - 1 )
#echo "last_gf: ${LAST_GF}"

CONSENSUS_PATH_BASE=$(echo ${CONSENSUS_PATH} | sed 's/\.fasta//')
echo "cpb: ${CONSENSUS_PATH_BASE}"

BWA_INDEX=$(echo ${CONSENSUS_PATH} | sed 's/\.fasta/_bwa_index/')
echo "bwa_index: ${BWA_INDEX}"

#CONSENSUS_PATH_RELATIVE_BASE=$(basename $CONSENSUS_PATH_BASE )
#echo "cprb: ${CONSENSUS_PATH_RELATIVE_BASE}"

CONSENSUS_NO_PATH_BASE=$(basename ${CONSENSUS_PATH} | sed 's/\.fasta//' )

TAG_LAST="14_PBJelly/${CONSENSUS_NO_PATH_BASE}" 
echo "tag_last: ${TAG_LAST}"

#TAG="04_farm/01_trimmed/${LINE}_bwa_alt_${I}"

TAG_DIR=$(dirname ${TAG_LAST})
echo "tag_dir: ${TAG_DIR}"

echo "variables made"

echo ${CONSENSUS_PATH}
echo ${TAG_LAST}


#******************************START MAPPING*****************************************

#module load OpenMPI/2.1.2
module load GNU/6.4.0-2.28  OpenMPI/2.1.1
module load BWA/0.7.17

file_check=$(ls 07_BAM/${TAG}_markdup.bam 2> /dev/null | wc -l )

if [ ${file_check} != 0 ] \
&& [ ${FORCE} == "false" ]; then
        echo "BAM Markdup exists, and FORCE set to false, NOT running Picard MarkDuplicates"
else
        echo "Either Markdup BAM doesnt exist, or force not set to false (${FORCE})"
        echo "Running Picard MarkDuplicates"; echo ""
        java -jar $PICARD/MarkDuplicates.jar \
        M=07_BAM/${TAG}_dup_stats.txt \
        I=07_BAM/${TAG}.bam \
        REMOVE_DUPLICATES=true \
        O=07_BAM/${TAG}_markdup.bam
        echo ""; echo "Finished Picard MarkDuplicates"; date; echo ""
fi


samtools faidx ${CONSENSUS_PATH}
bwa index -p ${BWA_INDEX} ${CONSENSUS_PATH}
mkdir -p 06_SAM/${TAG_DIR}

file_check=$(ls 06_SAM/${TAG_LAST}.sam 2> /dev/null | wc -l )

if [ ${file_check} != 0 ] \
&& [ ${FORCE} == "false" ]; then
	echo "SAM mapping to last consensus exists, and FORCE set to false, NOT running BWA MEM"
else
	echo "Either SAM mapping to last consensus doesnt exist, or force not set to false (${FORCE})"
	echo "Running BWA MEM"; echo ""

	rm -f 06_SAM/${TAG_LAST}.sam

	#bwa mem reads on the last consensus made

	bwa mem -R'@RG\tID:group1\tSM:sample1\tPL:illumina\tLB:lib1\tPU:unit1' \
	-t 39 \
	${BWA_INDEX} \
	${READ_DIR}/${LINE}_1_paired.fastq \
	${READ_DIR}/${LINE}_2_paired.fastq \
	> 06_SAM/${TAG_LAST}.sam

	echo ""; echo "Finished running BWA MEM"; date; echo ""
fi


module load OpenMPI/2.1.2


mkdir -p 07_BAM/${TAG_DIR}

#sam_to_bam
file_check=$(ls 07_BAM/${TAG_LAST}.bam 2> /dev/null | wc -l )

if [ ${file_check} != 0 ] \
&& [ ${FORCE} == "false" ]; then
	echo "Picard SortSam output exists, and FORCE set to false, NOT running	Picard SortSam"
else
	echo "Either Picard SortSam doesnt exist, or force not set to false (${FORCE})"
	echo "Running Picard SortSam"; echo ""
	java -jar $PICARD/SortSam.jar \
	INPUT=06_SAM/${TAG_LAST}.sam \
	OUTPUT=07_BAM/${TAG_LAST}.bam \
	SORT_ORDER=coordinate

	echo ""; echo "Finished running Pircard SortSam"; date; echo ""
fi

##******************************END MAPPING*****************************************

#******************************START MASKING*****************************************
echo "start mask:"; date; echo ""

#module load BEDTools

mkdir -p 13_BED/${TAG_DIR}

#just checking first one, rerun second no matter what since simple grep
file_check=$(ls 13_BED/${TAG_LAST}_all.bed 2> /dev/null | wc -l )

if [ ${file_check} != 0 ] \
&& [ ${FORCE} == "false" ]; then
	echo "Bed file exists, and FORCE set to false, NOT running bedGenomeCov"
else
 	echo "Bed file doesnt exist, or force not set to false (${FORCE})"
	echo "Running bedtools genomecov"; echo ""

	rm -f 13_BED/${TAG_LAST}_all.bed

	bedtools genomecov -bga -ibam 07_BAM/${TAG_LAST}.bam > 13_BED/${TAG_LAST}_all.bed

	echo ""; echo "Finished running bedtools genomecov"; date; echo ""
fi

rm -f 13_BED/${TAG_LAST}_empty.bed
grep -w 0$ 13_BED/${TAG_LAST}_all.bed > 13_BED/${TAG_LAST}_empty.bed

#This step also quick so just rerun anyway
bedtools maskfasta -fi ${CONSENSUS_PATH} -bed 13_BED/${TAG_LAST}_empty.bed -fo ${CONSENSUS_PATH_BASE}_masked.fasta
echo "end mask:"; date; echo ""

##******************************END MASKING*****************************************

##******************************START MASURCA*****************************************

#collect unmapped reads
	#maybe shoulders too

rm -f 09_UNMAPPED/${TAG_LAST}_um.sam
rm -f 09_UNMAPPED/${TAG_LAST}_um_ids.lst

mkdir -p 09_UNMAPPED/${TAG_DIR}
mkdir -p 09_UNMAPPED/00_shoulder/${TAG_DIR}

samtools view -S -f4 07_BAM/${TAG_LAST}.bam > 09_UNMAPPED/${TAG_LAST}_um.sam

cut -f1 09_UNMAPPED/${TAG_LAST}_um.sam | sort | uniq > 09_UNMAPPED/${TAG_LAST}_um_ids.lst

#module load seqtk/1.0

TAG_DOT=$(echo ${TAG_LAST} | sed 's|.*/||' | sed 's/_/\./g' )
#TAG_DIR=$(dirname ${TAG_LAST})

rm -f 09_UNMAPPED/${TAG_DIR}/${TAG_DOT}.1.um.fastq
rm -f 09_UNMAPPED/${TAG_DIR}/${TAG_DOT}.2.um.fastq

seqtk subseq ${READ_DIR}/${LINE}_1_paired.fastq 09_UNMAPPED/${TAG_LAST}_um_ids.lst > 09_UNMAPPED/${TAG_DIR}/${TAG_DOT}.1.um.fastq
seqtk subseq ${READ_DIR}/${LINE}_2_paired.fastq 09_UNMAPPED/${TAG_LAST}_um_ids.lst > 09_UNMAPPED/${TAG_DIR}/${TAG_DOT}.2.um.fastq

echo "finished creating unmapped fastas:"
date

#collect shoulders and add to unmapped read pool:

#to one line because scripts currently can't handle wrapped fasta
/mnt/research/edgerpat_lab/AlanY/01_athal_cns/perl/fasta_one_line.pl -f ${CONSENSUS_PATH_BASE}_masked.fasta -o ${CONSENSUS_PATH_BASE}_one_line.fasta
/mnt/research/edgerpat_lab/AlanY/01_athal_cns/perl/split_fasta_by_Ns.pl -f ${CONSENSUS_PATH_BASE}_one_line.fasta -o ${CONSENSUS_PATH_BASE}_split.fasta
/mnt/research/edgerpat_lab/AlanY/01_athal_cns/perl/fasta_shoulder.pl -f ${CONSENSUS_PATH_BASE}_split.fasta -o ${CONSENSUS_PATH_BASE}_split_shoulder.fasta --min 50 --length 200
/mnt/research/edgerpat_lab/AlanY/01_athal_cns/perl/fasta_to_fastq.pl -f ${CONSENSUS_PATH_BASE}_split_shoulder.fasta -o ${CONSENSUS_PATH_BASE}_split_shoulder.fastq

rm -f 09_UNMAPPED/00_shoulder/${TAG_DIR}/${TAG_DOT}.1.um.200s.fastq
rm -f 09_UNMAPPED/00_shoulder/${TAG_DIR}/${TAG_DOT}.2.um.200s.fastq
#200s for 200 bp shoulders
cat 09_UNMAPPED/${TAG_DIR}/${TAG_DOT}.1.um.fastq ${CONSENSUS_PATH_BASE}_split_shoulder.fastq > 09_UNMAPPED/00_shoulder/${TAG_DIR}/${TAG_DOT}.1.um.200s.fastq
cat 09_UNMAPPED/${TAG_DIR}/${TAG_DOT}.2.um.fastq ${CONSENSUS_PATH_BASE}_split_shoulder.fastq > 09_UNMAPPED/00_shoulder/${TAG_DIR}/${TAG_DOT}.2.um.200s.fastq

#for now since running many times over:
rm -rf 10_masurca/${TAG_LAST}
mkdir -p 10_masurca/${TAG_LAST}

UNMAPPED1="${WKDIR}/09_UNMAPPED/00_shoulder/${TAG_DIR}/${TAG_DOT}.1.um.200s.fastq"
UNMAPPED2="${WKDIR}/09_UNMAPPED/00_shoulder/${TAG_DIR}/${TAG_DOT}.2.um.200s.fastq"

cd /mnt/research/edgerpat_lab/AlanY/MaSuRCA-3.2.6_6.2

rm -f ${WKDIR}/10_masurca/${TAG_LAST}/temp_config.txt
rm -f ${WKDIR}/10_masurca/${TAG_LAST}/${TAG_DOT}_sr_config.txt

sed "s|/mnt/research/edgerpat_lab/AlanY/01_athal_cns/assembly_golicz16/01_raw_data/test.1.fastq|${UNMAPPED1}|" sr_config_test.txt > ${WKDIR}/10_masurca/${TAG_LAST}/temp_config.txt
sed "s|/mnt/research/edgerpat_lab/AlanY/01_athal_cns/assembly_golicz16/01_raw_data/test.2.fastq|${UNMAPPED2}|" ${WKDIR}/10_masurca/${TAG_LAST}/temp_config.txt > ${WKDIR}/10_masurca/${TAG_LAST}/${TAG_DOT}_sr_config.txt

export PATH=$PATH:/mnt/research/edgerpat_lab/AlanY/MaSuRCA-3.2.6_6.2/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/mnt/research/edgerpat_lab/AlanY/MaSuRCA-3.2.6_6.2/lib

module load GNU/6.2

cd ${WKDIR}/10_masurca/${TAG_LAST}

#TAG_BASE=$(echo ${TAG} | sed 's|*/||g')

/mnt/research/edgerpat_lab/AlanY/MaSuRCA-3.2.6_6.2/bin/masurca ${TAG_DOT}_sr_config.txt

sed -i.bak 's/rm -rf CA//' assemble.sh

./assemble.sh

echo "Finished MaSuRCA assembly"
date
##******************************END MASURCA*****************************************

##******************************START PBJELLY*****************************************
cd ${WKDIR}

#mass module load, back to GNU/4.4.5
module load GNU/4.4.5
module load OpenMPI/1.4.3
module load PBSuite

#run PBJ:
#need to set up some loops to make sure steps finish
#some substitutions to specify things in the template
#CONSENSUS_NO_PATH_BASE=$(basename ${CONSENSUS_PATH} | sed 's/\.fasta//' )

#if first one, tag doesn't have pbj in it yet
if [ ${I} -eq 1 ]; then
	TAG_PBJ="${CONSENSUS_NO_PATH_BASE}_pbj"
else
	TAG_PBJ=$(echo ${CONSENSUS_NO_PATH_BASE} | sed 's/_[0-9]\+$//' )
fi

echo "tag_pbj: ${TAG_PBJ}"

#fake quality file:
rm -rf ${WKDIR}/10_masurca/${TAG_LAST}/CA/final.genome.scf.qual

#check masurca worked, if not, and output file doesn't exist, a few code blocks later will be infinite loop, not great so shutting down here if that happens

if [ ! -f ${WKDIR}/10_masurca/${TAG_LAST}/CA/final.genome.scf.fasta ]; then
	echo "Could not find tha MaSuRCA output:"
	echo "${WKDIR}/10_masurca/${TAG_LAST}/CA/final.genome.scf.fasta"
	echo "exiting"
	date
	exit 1;
fi

fakeQuals.py ${WKDIR}/10_masurca/${TAG_LAST}/CA/final.genome.scf.fasta ${WKDIR}/10_masurca/${TAG_LAST}/CA/final.genome.scf.qual

#make directories needed
#rm -rf step incase rerunning, best way to overwrite old stuff
rm -rf 14_PBJelly/${LINE}_gf/02_RUN/${I}
mkdir -p 14_PBJelly/${LINE}_gf/02_RUN/${I}
#for final output
#**********DO NOT DO THIS, IT STORES ALL THE OTHER RUNS
##rm -rf 14_PBJelly/${LINE}_gf/01_OUT/
mkdir -p 14_PBJelly/${LINE}_gf/01_OUT/

#toss in template
#if fails, then need to go to this line of code and specify where copy of default template is, going to use one in hpcc build of PBJelly
cp /opt/software/PBSuite/14.7.14--GCC-4.4.5/docs/TemplateProtocol.xml 14_PBJelly/${LINE}_gf/02_RUN/${I}

#EDITING TEMPLATE
#specify reference,
#if not using a genome gnerated by a round of PBJelly, use the masked version... hmm maybe should always used the masked version, yea because have zero reads mapping there
sed -i "s|<reference>/FULL/PATH/TO__/PBJelly/data/reference/lambda.fasta</reference>|<reference>${CONSENSUS_PATH_BASE}_masked.fasta</reference>|" 14_PBJelly/${LINE}_gf/02_RUN/${I}/TemplateProtocol.xml

#specify output dir,
sed -i "s|<outputDir>/FULL/PATH/TO__/PBJelly/lambdaExample/</outputDir>|<outputDir>${WKDIR}/14_PBJelly/${LINE}_gf/02_RUN/${I}</outputDir>|" 14_PBJelly/${LINE}_gf/02_RUN/${I}/TemplateProtocol.xml

#specify resources
sed -i "s|<command\ notes=\"For\ PBS/Moab\">echo\ '\${CMD}'\ \|\ msub\ -N\ \"\${JOBNAME}\"\ -o\ \${STDOUT}\ -e\ \${STDERR}\ -l\ nodes=1:ppn=8,mem=48000mb</command>|<command>\${CMD}\ \${JOBNAME}\ 2>\ \${STDERR}\ 1>\ \${STDOUT}\ \&amp;</command>|" 14_PBJelly/${LINE}_gf/02_RUN/${I}/TemplateProtocol.xml
sed -i "s|<command\ notes=\"For\ single\ node,\ multi-core\ machines\"\ >\${CMD}\ \${JOBNAME}\ 2>\ \${STDERR}\ 1>\ \${STDOUT}\ &amp;</command>||" 14_PBJelly/${LINE}_gf/02_RUN/${I}/TemplateProtocol.xml

#specify assembled unmapped reads
sed -i "s|<input baseDir=\"/FULL/PATH/TO__/PBJelly/lambdaExample/data/reads/\">|<input baseDir=\"${WKDIR}/10_masurca/${TAG_LAST}/CA\">|" 14_PBJelly/${LINE}_gf/02_RUN/${I}/TemplateProtocol.xml
sed -i "s|<job>pacbioReads.fasta</job>|<job>final.genome.scf.fasta</job>|" 14_PBJelly/${LINE}_gf/02_RUN/${I}/TemplateProtocol.xml

#time to run
#SETUP
Jelly.py setup 14_PBJelly/${LINE}_gf/02_RUN/${I}/TemplateProtocol.xml

if [[ $? -ne 0 ]]; then
	echo "Failed Setup"
	exit
fi

##wait loop, check every depends how long run takes:
HOW_LONG=0
WAIT=1
CONSENSUS_DIR=$(dirname ${CONSENSUS_PATH})
while [[ "$WAIT" -ne "0" ]]
do 
        sleep 1
        HOW_LONG=$(expr ${HOW_LONG} + 1)
        grep "Finished!" ${CONSENSUS_DIR}/setup_chunk0.err > /dev/null
        if [[ $? -eq 0 ]]; then
                #break out of the loop
                WAIT=0
        fi
done
echo "${HOW_LONG}:Setup (roughly) seconds"

#MAPPING
Jelly.py mapping 14_PBJelly/${LINE}_gf/02_RUN/${I}/TemplateProtocol.xml

if [[ $? -ne 0 ]]; then
        echo "Failed Mapping"
        exit
fi

HOW_LONG=0
WAIT=1
while [[ "$WAIT" -ne "0" ]]
do
  	sleep 1
        HOW_LONG=$(expr ${HOW_LONG} + 1)
        grep "tails\ mapped" 14_PBJelly/${LINE}_gf/02_RUN/${I}/mapping/mapping_chunk0.err > /dev/null
        if [[ $? -eq 0 ]]; then
                #break out of the loop
                WAIT=0
        fi
	grep "exiting" 14_PBJelly/${LINE}_gf/02_RUN/${I}/mapping/mapping_chunk0.err > /dev/null
	if [[ $? -eq 0 ]]; then
		#kill it
		echo "No tails mapped in PBJ mapping stage"
		echo "Round ${I}"
		echo "on ${LINE}"
		echo "most often caused by shallow read depth"
		echo "consider increasing shoulder size, or generating more data"
		exit
	fi
done
echo "${HOW_LONG}:Mapping (roughly) seconds"

#SUPPORT
Jelly.py support 14_PBJelly/${LINE}_gf/02_RUN/${I}/TemplateProtocol.xml

if [[ $? -ne 0 ]]; then
        echo "Failed Support"
        exit
fi

HOW_LONG=0
WAIT=1
while [[ "$WAIT" -ne "0" ]]
do
  	sleep 1
        HOW_LONG=$(expr ${HOW_LONG} + 1)
        grep "Finished" 14_PBJelly/${LINE}_gf/02_RUN/${I}/support/support_chunk0.err > /dev/null
        if [[ $? -eq 0 ]]; then
                #break out of the loop
                WAIT=0
        fi
done
echo "${HOW_LONG}:Support (roughly) seconds"

#EXTRACTION
Jelly.py extraction 14_PBJelly/${LINE}_gf/02_RUN/${I}/TemplateProtocol.xml

if [[ $? -ne 0 ]]; then
        echo "Failed Extraction"
        exit
fi

HOW_LONG=0
WAIT=1
while [[ "$WAIT" -ne "0" ]]
do
  	sleep 1
        HOW_LONG=$(expr ${HOW_LONG} + 1)
        grep "Finished" 14_PBJelly/${LINE}_gf/02_RUN/${I}/assembly/extraction_chunk0.err > /dev/null
        if [[ $? -eq 0 ]]; then
                #break out of the loop
                WAIT=0
        fi
done
echo "${HOW_LONG}:Extraction (roughly) seconds"

#ASSEMBLY
Jelly.py assembly 14_PBJelly/${LINE}_gf/02_RUN/${I}/TemplateProtocol.xml

if [[ $? -ne 0 ]]; then
        echo "Failed Assembly"
        exit
fi

IN_ROUNDS=$(ls -l 14_PBJelly/${LINE}_gf/02_RUN/${I}/assembly/ | grep -c "^d")
HOW_LONG=0
WAIT=1
while [[ "$WAIT" -ne "0" ]]
do
  	sleep 1
        HOW_LONG=$(expr ${HOW_LONG} + 1)
        OUT_ROUNDS=$(wc -l 14_PBJelly/${LINE}_gf/02_RUN/${I}/assembly/assembly_chunk0.out | cut -f1 -d" ")
        if [[ "$OUT_ROUNDS" -eq "$IN_ROUNDS" ]]; then
                #break out of the loop
                WAIT=0
        fi
done
echo "${HOW_LONG}:Assembly (roughly) seconds"

#OUTPUT
Jelly.py output 14_PBJelly/${LINE}_gf/02_RUN/${I}/TemplateProtocol.xml

if [[ $? -ne 0 ]]; then
        echo "Failed Output"
        exit
fi

HOW_LONG=0
WAIT=1
while [[ "$WAIT" -ne "0" ]]
do
  	sleep 1
        HOW_LONG=$(expr ${HOW_LONG} + 1)
        grep "Finished" 14_PBJelly/${LINE}_gf/02_RUN/${I}/assembly/output_chunk0.err
        if [[ $? -eq 0 ]]; then
                #break out of the loop
                WAIT=0
        fi
done
echo "${HOW_LONG}:Output (roughly) seconds"

#remove previous run if rerunning
rm -rf 14_PBJelly/${LINE}_gf/01_OUT/${TAG_PBJ}_${I}.fasta
cp 14_PBJelly/${LINE}_gf/02_RUN/${I}/jelly.out.fasta 14_PBJelly/${LINE}_gf/01_OUT/${TAG_PBJ}_${I}.fasta

#rename chromosomes to what they were
#executable bash script
#source /mnt/research/edgerpat_lab/AlanY/01_athal_cns/qsub/viigf_rename_chr.sh

echo "I did it!! well atleast for the ${I} round"

date

#check for what walltime to use:
if [[ ${MANUAL_WALLTIME} =~ ^[nN] ]]; then
        READ_FILE_SIZE=$(stat ${READ_DIR}/${LINE}_1_paired.fastq | grep "Size:" | cut -f4 -d" " | cut -f1 -d"	")
        if [[ ${READ_FILE_SIZE} -gt 12000000000 ]]; then
                #more than 4 hours for var_inc
                VAR_INC_WALLTIME=12:00:00
		PBJ_WALLTIME=4:00:00
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


#define consensus for next round
NEXT_ROUND_CONSENSUS="${WKDIR}/14_PBJelly/${LINE}_gf/01_OUT/${TAG_PBJ}_${I}.fasta"

#next round
if [ "$I" -lt "$PBJ_FINISH" ]; then
	#next step
	I=$( expr ${I} + 1 )
	qsub -t ${PBS_ARRAYID} -v CONFIG=${CONFIG},I=${I},CONSENSUS_PATH=${NEXT_ROUND_CONSENSUS} -l walltime=${PBJ_WALLTIME},nodes=50:ppn=1,mem=50gb /mnt/research/edgerpat_lab/AlanY/01_athal_cns/qsub/viigf_pbj.sh
	echo "submitted within the pbj script"
	echo "pbs_array: ${PBS_ARRAYID}"
        echo "line: ${LINE}"
	echo "with consensus: ${NEXT_ROUND_CONSENSUS}"
	echo "with loop variable: ${I}"
else
	echo "Success!! Visit ${OUTPUT_STATS} for summary statistics of the pipeline"
	echo "I hope this helped, remember, always check the quality yourself to see if anything funky happened"
	echo "If any files need revisited, check $SCRATCH"
	echo "Any major problems, you probably know where to find me"
	if [[ ${CONFIG}	=~ *launch_edited.conf$	]]
	then
                rm -f ${CONFIG}
        fi
fi

#run the cleaner
if [[ ${RUN_CLEANERS} =~ ^[Yy] ]]; then
        echo "running the cleaner"
        source ${CLEANER_PBJ}
else
        echo "NOT running cleaner"
fi

qstat -f ${PBS_JOBID}
