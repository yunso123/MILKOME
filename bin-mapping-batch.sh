#!/bin/bash

cd /home/projects/dtu_00032/analysis/milk-cohort
mkdir -p bins/new-ani98/bins-mapping
mkdir -p bins/new-ani98/bin-mapping-scripts
scriptpath="bins/new-ani98/bin-mapping-scripts"

for f in data/host_removed/*.1.fastq.gz; do
    b=$(basename "$f" '.1.fastq.gz')
    r="data/host_removed/${b}.2.fastq.gz"
    o="bins/new-ani98/bins-mapping/HR_${b}.txt"
    
    script_name="${scriptpath}/run_HR_${b}.sh"
    log_name="${scriptpath}/HR_${b}.log"
    err_name="${scriptpath}/HR_${b}.err"

    echo "#!/bin/bash" > $script_name
    echo "#PBS -W group_list=dtu_00032 -A dtu_00032" >> $script_name
    echo "#PBS -N bin-map-${b}" >> $script_name
    echo "#PBS -e ${err_name}" >> $script_name
    echo "#PBS -o ${log_name}" >> $script_name
    echo "#PBS -l nodes=1:thinnode:ppn=40" >> $script_name
    echo "#PBS -l mem=188gb" >> $script_name
    echo "#PBS -l walltime=24:00:00" >> $script_name
    echo "" >> $script_name
    echo "module load tools" >> $script_name
    echo "module load gcc/12.2.0" >> $script_name
    echo "module load fastani/1.34" >> $script_name
    echo "module load dashing/0.4.0" >> $script_name
    echo "module load minimap2/2.24r1122" >> $script_name
    echo "module load samtools/1.18" >> $script_name
    echo "module load coverm/0.7.0" >> $script_name
    echo "cd /home/projects/dtu_00032/analysis/milk-cohort" >> $script_name
    echo "" >> $script_name
    echo "coverm genome \\" >> $script_name
    echo "    -t 40 \\" >> $script_name
    echo "    -1 $f \\" >> $script_name
    echo "    -2 $r \\" >> $script_name
    echo "    -o $o \\" >> $script_name
    echo "    --genome-fasta-list bins/new-ani98/bins-75-5.txt \\" >> $script_name
    echo "    -x fna;" >> $script_name

    chmod +x $script_name
    qsub $script_name
done
