#!/bin/bash

# Define sample list
samples=$(cat "/home/projects/dtu_00032/analysis/milk-cohort/meta/non-stool.txt")

# Define base directory and paths
base_dir="/home/projects/dtu_00032/analysis/milk-cohort"
scriptpath="${base_dir}/meta/job_scripts/mp"
logpath="${base_dir}/meta/logs/mp"
mkdir -p $scriptpath $logpath

# Iterate through each sample and create a PBS script
for N in $samples; do
    script_name="${scriptpath}/${N}_mp.sh"
    log_name="${logpath}/${N}_mp.log"
    err_name="${logpath}/${N}_mp.err"

    # Create PBS script
    cat <<EOL > $script_name
#!/bin/sh
#PBS -W group_list=dtu_00032 -A dtu_00032
#PBS -N ${N}-Bif-MP
#PBS -e ${err_name}
#PBS -o ${log_name}
#PBS -l nodes=1:thinnode:ppn=40
#PBS -l mem=188gb
#PBS -l walltime=24:00:00

module load tools
module load metaphlan/4.1.0


cd ${base_dir}

tf="data/trimmomatic/${N}_trimmed_1.fq.gz"
tr="data/trimmomatic/${N}_trimmed_2.fq.gz"
mkdir -p meta/bif-metaphlan meta/bif-metaphlan/bowtie2 meta/bif-metaphlan/sams meta/bif-metaphlan/original meta/bif-metaphlan/convGTDB
d="/home/projects/dtu_00032/db/metaphlan_bif_db_202307"

if [[ -f \$tf && -f \$tr ]]; then
    sam="meta/bif-metaphlan/sams/${N}.sam"
    bowtie="meta/bif-metaphlan/bowtie2/${N}.bowtie2.bz2"
    mpfile="meta/bif-metaphlan/original/${N}_profiled.txt"
    convfile="meta/bif-metaphlan/convGTDB/${N}_conv_profiled.txt"
    
    echo "Running Metaphlan with ${N}"

    metaphlan "\${tf},\${tr}" \
        --nproc 40 \
        --input_type fastq \
        --bowtie2db \$d \
        -x "mpa_vJun23_CHOCOPhlAnSGB_202403_lon_subsp" \
        -s \$sam \
        --bowtie2out \$bowtie \
        -o \$mpfile;

    #python /services/tools/metaphlan/4.1.0/bin/sgb_to_gtdb_profile.py \
    #    -i \$mpfile \
    #    -o  \$convfile \
    #    -d "/home/projects/dtu_00032/db/metaphlan_db_2307/mpa_vJun23_CHOCOPhlAnSGB_202307.pkl";
else
    echo "Error: Input files for sample ${N} not found"
    exit 1
fi

EOL

    # Make the script executable and submit it
    chmod +x $script_name
    qsub $script_name
done




