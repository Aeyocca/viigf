#!/bin/sh
#Alan E Yocca
#07-05-18
#launch viigf pipeline from anywhere
#01-15-19
#slurm conversion

set -e

#making things full path so this doesn't matter
#cd ${WKDIR}

#read in cmd line args 
#launch at trimmomatic, var_inc, or pbj
#I say still feed in diff config yea that good idead
#read in which line of infile

print_usage() {
        echo "Usage:"
        echo "$0"
        echo "  --step:"
        echo "          <which part of the viigf pipe to begin from>"
	echo "			<options: trimmomatic, var_inc, pbj>"
	echo "			<******if pbj specified, --consensus flag required*******>"
	echo "  --consensus:"
	echo "		<full path of consensus genome to gap fill with PBJ>"
	echo "		<only specify if pbj step selected>"
        echo "  --array:"
        echo "          <which lines of INFILE to run the pipeline on>"
	echo "			<e.g. --array 1-4>"
        echo "  --config:"
        echo "          <full path of the configuration file>"
	echo "			<e.g. --config /mnt/research/edgerpat_lab/AlanY/01_athal_cns/qsub/viigf_config.conf>"
	echo "		<*****ASSUMES ALL VIIGF FILE executables are in same directory, if they aren't they should be>"
	echo "  --continue:"
	echo "		<continue past the specified step or not>"
	echo "			<def: yes, or whatever is in config file under variable CONTINUE>"
	echo "  --var_start:"
	echo "		<specify manually the var_start, will ignore whatever is in the config file>"
	echo "  --var_last:"
	echo "		<specify manually the var_last, will ignore whatever is in the config file>"
	echo "  --pbj_start:"
	echo "		<specify manually the pbj_start, will ignore whatever is in the config file>"
        echo "  --pbj_finish:"
        echo "	        <specify manually the pbj_finish, will ignore whatever is in the config file>"
	echo " --script_dir:"
	echo "		<specify script directory if different than>"
        }

# goddamn this is annoying to have to put at the beginning, o whale
# credit to stack overflow mcoolive
# Transform long options to short ones
for arg in "$@"; do
  shift
  case "$arg" in
    "--step") set -- "$@" "-s" ;;
    "--consensus") set -- "$@" "-g" ;;
    "--array") set -- "$@" "-a" ;;
    "--config") set -- "$@" "-c" ;;
    "--continue") set -- "$@" "-m" ;;
    "--pbj_start") set -- "$@" "-p" ;;
    "--pbj_finish") set -- "$@" "-f" ;;
    "--var_start") set -- "$@" "-v" ;;
    "--var_last") set -- "$@" "-l" ;;
    "--script_dir") set -- "$@" "-d" ;;
    *)        set -- "$@" "$arg"
  esac
done

# Default behavior
dflag=false
sflag=false
gflag=false	#only mandatory if --step pbj specified
aflag=false
cflag=false
mflag=false
pflag=false
fflag=false
vflag=false
lflag=false

# Parse short options
OPTIND=1
while getopts ":s:g:a:c:m:p:f:v:l:d:" opt
do
  case "$opt" in
    "s") STEP=$OPTARG; sflag=true ;;
    "g") CONSENSUS_PATH=$OPTARG; gflag=true ;;
    "a") ARRAYID=$OPTARG; aflag=true ;;
    "c") CONFIG=$OPTARG; cflag=true ;;
    "m") CONTINUE_CMD=$OPTARG; mflag=true ;;
    "p") PBJ_START=$OPTARG; pflag=true ;;
    "f") PBJ_FINISH=$OPTARG; fflag=true ;;
    "v") VAR_START=$OPTARG; vflag=true ;;
    "d") SCRIPT_DIR=$OPTARG; dflag=true ;;
    "l") VAR_LAST=$OPTARG; lflag=true ;;
    "?") print_usage >&2; exit 1 ;;
  esac
done
shift $(expr $OPTIND - 1) # remove options from positional parameter

#check for some mandatory options:
if ! ${sflag} || ! ${aflag} || ! ${cflag}
then
	print_usage
	echo "MISSING SOME MANDATORY ARGUMENTS!!"
	exit 1
fi


#if continue specified fix in config file
#if [ ${CONTINUE} =~ ^[Yy] ]; then
#	sed -i 's/CONTINUE*/CONTINUE="Yawwwwwww"/' ${CONFIG}
#elif [ ${CONTINUE} =~ ^[Nn] ] then
#	sed -i 's/CONTINUE*/CONTINUE="Nokle dokle artle chokle"/' ${CONFIG} 
#fi
if [ -z ${CONTINUE_CMD+x} ]; then
        #not defined, stick with default in config file
        echo "Continue is not defined, sticking with default"
else
        sed -i "s/CONTINUE.*/CONTINUE=\"${CONTINUE_CMD}\"/" ${CONFIG}
fi

#add in the pbj_start / other thingies
if ${pflag} || ${fflag} || ${vflag} || ${lflag}
then
	DATE=$(date)
	DATE_SCORE=$(echo ${DATE} | sed "s/\s/_/g" )
	cp ${CONFIG} ${CONFIG}_${DATE_SCORE}
	CONFIG="${CONFIG}_${DATE_SCORE}"
	if ${pflag}
	then
		sed -i "s/PBJ_START=[0-9]+/PBJ_START=${PBJ_START}/" ${CONFIG}
	fi
	if ${fflag}
        then	
                sed -i "s/PBJ_FINISH=[0-9]+/PBJ_FINISH=${PBJ_FINISH}/" ${CONFIG}
        fi
        if ${vflag}
	then
            	sed -i "s/VAR_START=[0-9]+/VAR_START=${VAR_START}/" ${CONFIG}
        fi
        if ${lflag}
	then
            	sed -i "s/VAR_LAST=[0-9]+/VAR_LAST=${VAR_LAST}/" ${CONFIG}
        fi
fi

#test
#echo "${STEP}"
#echo "${ARRAYID}"
#echo "${CONFIG}"

source ${CONFIG}


if ! ${dflag}
then
	SCRIPT_DIR=$(dirname ${CONFIG})
fi


#for now, considering adding READ_DIR to config so can borrow reads for other project and make separate assembly space
#READ_DIR="${WKDIR}"
#depreciated 07-06-18

#echo "${VIIGF_DIR}"

#Make output directory
mkdir -p ${EO_DIR}


#split up the --array into an array to loop through
array_1=(`echo ${ARRAYID} | sed 's/,/\n/g'`)
for I in "${array_1[@]}"
do
	if [[ ${I} =~ - ]]; then
		array_2=(`echo ${I} | sed 's/-/\n/g'`)
		array_3=( $(seq ${array_2[0]} ${array_2[1]}) )
		for j in "${array_3[@]}"
		do
			sub_array[${#sub_array[@]}]="${j}"
#			echo "${j}"
		done
	else
		sub_array[${#sub_array[@]}]="${I}"
	fi
done

#for k in "${sub_array[@]}"
#do
#	echo "${k}"
#	echo "next"
#done
#echo "Killed"
#exit

for each_array in "${sub_array[@]}"
do
	SLURM_ARRAY_TASK_ID=${each_array}
	source ${CONFIG}

	if [[ ${MANUAL_WALLTIME} =~ ^[nN] ]]; then
		#SLURM_ARRAY_TASK_ID=${each_array}
		#source ${CONFIG}
		#check for consensus path if pbj specified
		if [[ ${STEP} == "pbj" ]]; then
        		if ! ${gflag}
        		then
#   		            print_usage
#       		        echo "MISSING SOME MANDATORY ARGUMENTS!!"
        		        echo "Consensus path not given, using default consensus path based on PBJ_start"
#       		        exit 1
        		        if [[ ${PBJ_START} == "1" ]]; then
                		        echo "Setting consensus path to: "
                		        source ${CONFIG}
                		        CONSENSUS_PATH="${WKDIR}/08_consensus/04_farm/01_trimmed/${LINE}_bwa_alt_3.fasta"
                		        echo "  ${CONSENSUS_PATH}"
               		 else
                		    	echo "Setting consensus path to: "
                		        LAST=$(expr ${PBJ_START} - 1)
                		        source ${CONFIG}
                		        CONSENSUS_PATH="${WKDIR}/14_PBJelly/${LINE}_gf/01_OUT/${LINE}_bwa_alt_3_pbj_${LAST}.fasta"
				fi
			fi
		fi
		if [ ${STEP} == "trimmomatic" ]; then
			#qsub -t ${PBS_ARRAYID} -l walltime=4:00:00,nodes=40:ppn=1,mem=20gb -v CONFIG=${CONFIG},CONSENSUS_PATH=${CONSENSUS_PATH} ${VIIGF_DIR}/viigf_${STEP}.sh
			sbatch --array=${SLURM_ARRAY_TASK_ID} \
			--time=4:00:00 \
			--ntasks=40 \
			--mem=20gb \
			--export=CONFIG=${CONFIG},CONSENSUS_PATH=${CONSENSUS_PATH} \
			--output=${EO_DIR}/%x-%A_%a.SLURMout \
			${SCRIPT_DIR}/viigf_${STEP}.sbatch

			echo "Started pipe from: ${STEP}"
			date
			echo "Running on line(s) ${SLUR_ARRAY_TASK_ID} (${LINE}) from file: ${INFILE}"
			continue 1
		fi
		READ_FILE_SIZE=$(stat ${READ_DIR}/${LINE}_1_paired.fastq | grep "Size:" | cut -f4 -d" " | cut -f1 -d"	" )
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
		if [ ${STEP} == "var_inc" ]; then
		        WALLTIME=${VAR_INC_WALLTIME}
			#qsub -t ${PBS_ARRAYID} -l walltime=${WALLTIME},nodes=40:ppn=1,mem=40gb -v CONFIG=${CONFIG} ${VIIGF_DIR}/viigf_${STEP}.sh
#			sub_var_inc=true
			sbatch --array=${SLURM_ARRAY_TASK_ID} \
			--time=${WALLTIME} \
			--ntasks=40 \
			--mem=40gb \
			--export=CONFIG=${CONFIG} \
			--output=${EO_DIR}/%x-%A_%a.SLURMout \
			${SCRIPT_DIR}/viigf_${STEP}.sbatch

			DATE=$(date)
			echo "Started pipe from: ${STEP}"
			echo "Started on round ${VAR_START}"
			echo "${DATE}"
			echo "Running on line ${SLURM_ARRAY_TASK_ID} (${LINE}) from file: ${INFILE}"
			echo "with walltime: ${WALLTIME}"
		elif [ ${STEP} == "pbj" ]; then
		        WALLTIME=${PBJ_WALLTIME}
			#qsub -t ${PBS_ARRAYID} -l walltime=${WALLTIME},nodes=40:ppn=1,mem=100gb -v CONFIG=${CONFIG},CONSENSUS_PATH=${CONSENSUS_PATH} ${VIIGF_DIR}/viigf_${STEP}.sh
#			sub_pbj=true
			sbatch --array=${SLURM_ARRAY_TASK_ID} \
			--time=${WALLTIME} \
			--ntasks=40 \
			--mem=100gb \
			--export=CONFIG=${CONFIG},CONSENSUS_PATH=${CONSENSUS_PATH} \
			--output=${EO_DIR}/%x-%A_%a.SLURMout \
			${SCRIPT_DIR}/viigf_${STEP}.sbatch

			DATE=$(date)
                        echo "Started pipe from: ${STEP}"
			echo "Started on round: ${PBJ_START}"
                        echo "${DATE}"
                        echo "Running on line ${SLURM_ARRAY_TASK_ID} (${LINE}) from file: ${INFILE}"
                        echo "with walltime: ${WALLTIME}"
			echo "with consensus:"
			echo "${CONSENSUS_PATH}"
		fi
	else
		echo "using walltime specified in configuration	file"
		if [ ${STEP} == "var_inc" ]; then
			WALLTIME=${VAR_INC_WALLTIME}
			#qsub -t ${PBS_ARRAYID} -l walltime=${WALLTIME},nodes=40:ppn=1,mem=50gb -v CONFIG=${CONFIG},CONSENSUS_PATH=${CONSENSUS_PATH} ${VIIGF_DIR}/viigf_${STEP}.sh
			sbatch --array=${SLURM_ARRAY_TASK_ID} \
			--time=${WALLTIME} \
			--ntasks=40 \
			--mem=50gb \
			--export=CONFIG=${CONFIG},CONSENSUS_PATH=${CONSENSUS_PATH} \
			--output=${EO_DIR}/%x-%A_%a.SLURMout \
			${SCRIPT_DIR}/viigf_${STEP}.sbatch

			echo "Started pipe from: ${STEP}"
			date
			echo "Running on line(s) ${SLURM_ARRAY_TASK_ID} (${LINE}) from file: ${INFILE}"
			echo "with consensus:"
			echo "${CONSENSUS_PATH}"
		elif [ ${STEP} == "pbj" ]; then
		        WALLTIME=${PBJ_WALLTIME}
			#qsub -t ${PBS_ARRAYID} -l walltime=${WALLTIME},nodes=40:ppn=1,mem=50gb -v CONFIG=${CONFIG},CONSENSUS_PATH=${CONSENSUS_PATH} ${VIIGF_DIR}/viigf_${STEP}.sh
			sbatch --array=${SLURM_ARRAY_TASK_ID} \
			--time=${WALLTIME} \
			--ntasks=40 \
			--mem=50gb \
			--export=CONFIG=${CONFIG},CONSENSUS_PATH=${CONSENSUS_PATH} \
			--output=${EO_DIR}/%x-%A_%a.SLURMout \
			${SCRIPT_DIR}/viigf_${STEP}.sbatch

			echo "Started pipe from: ${STEP}"
			date
			echo "Running on line(s) ${SLURM_ARRAY_TASK_ID} (${LINE}) from file: ${INFILE}"
			echo "with consensus:"
			echo "${CONSENSUS_PATH}"
		elif [ ${STEP} == "trimmomatic" ]; then
			#qsub -t ${PBS_ARRAYID} -l walltime=4:00:00,nodes=40:ppn=1,mem=20gb -v CONFIG=${CONFIG},CONSENSUS_PATH=${CONSENSUS_PATH} ${VIIGF_DIR}/viigf_${STEP}.sh
			sbatch --array=${SLURM_ARRAY_TASK_ID} \
			--time=${WALLTIME} \
			--ntasks=40 \
			--mem=50gb \
			--export=CONFIG=${CONFIG},CONSENSUS_PATH=${CONSENSUS_PATH} \
			--output=${EO_DIR}/%x-%A_%a.SLURMout \
			${SCRIPT_DIR}/viigf_${STEP}.sbatch

			echo "Started pipe from: ${STEP}"
			date
			echo "Running on line(s) ${SLURM_ARRAY_TASK_ID} (${LINE}) from file: ${INFILE}"
		fi
	fi
done

#if [ ${STEP} == "trimmomatic" ]; then
#	WALLTIME="4:00:00"
#fi


#leave resource request static outside walltime for now, think about changing later

#if [ -z "$sub_var_inc" ] && [ -z "$sub_pbj" ] ; then
#	# leave CONSENSUS_PATH in here regardless of what step, will be reset in them, but required for pbj
#	qsub -t ${ARRAYID} -l walltime=${WALLTIME},nodes=40:ppn=1,mem=100gb -v CONFIG=${CONFIG},CONSENSUS_PATH=${CONSENSUS_PATH} ${VIIGF_DIR}/viigf_${STEP}.sh
#	DATE=$(date)
#	echo "Started pipe from: ${STEP}"
#	echo "${DATE}"
#	echo "Running on line(s) ${ARRAYID} from file: ${INFILE}"
#	echo "with consensus (if specified):"
#	echo "${CONSENSUS_PATH}"
#fi

#thing about adding later
#echo "All error/output directed to ${ERROR_OUTPUT}"

