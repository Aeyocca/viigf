# VIIGF
Alan E. Yocca
aeyap42@gmail.com
@Aeyocca (twitter)
Uploaded 03-27-19

This is a horrible acronym but it keeps the file names nice and concise

## Built for use on University's cluster. Compatible with SLURM submission system
- PBJelly part is still untested, working on that

# Inputs:
- Illumina paired-end reads
- Reference genome

- Run just the viigf_var_inc step to get a consensus genome from incorporating SNPs into the reference genome
- I found marginal improvement after 3 rounds
  - test your own though, compare number of variants that still exist after 3 rounds

# Configuration file:
- set variables to what you need them to be
- example configuration file: viigf_config_ragoo.conf

# To set it all off!
- make sure the configuration file is set how you would like it to be set
- make sure reads are in either 01_raw_reads/ within ${WKDIR} or 04_trimmed
  - depending if they are trimmed reads or not
- run ./viigf_launch.sh
  - type: ```$ ./viigf_launch.sh ``` without any arguments to see the command line options

# viigf_pbj.sbatch
- untested as of 03-27-19
- additionally, code written to run MaSuRCA, prefer other assembly methods though
  - since my MaSuRCA build is fragile





