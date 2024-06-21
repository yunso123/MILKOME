#!/bin/sh
#PBS -W group_list=dtu_00032 -A dtu_00032
#PBS -N stool-mp
#PBS -e stool-mp.err
#PBS -o stool-mp.log
#PBS -l nodes=1:thinnode:ppn=40
#PBS -l mem=188gb
#PBS -l walltime=48:00:00


module load tools
module load bowtie2/2.5.2

cd /home/projects/dtu_00032/analysis/milk-cohort
#mkdir data/host_removed
#mkdir data/bowtie

file=$(cat logs/stools.txt)
#for N in $file; do
#    f="data/trimmomatic/${N}_trimmed_1.fq.gz"; 
#    r="data/trimmomatic/${N}_trimmed_2.fq.gz"
#    echo 
#    echo 'processing' $N
#    bowtie2 -p 40 -x /home/projects/dtu_00032/GRCh38_noalt_as/GRCh38_noalt_as -1 $f -2 $r --un-conc-gz "data/host_removed/${N}" > "data/bowtie/${N}_mapped_unmapped.sam"
#done


module purge
module load tools
module load megahit/1.2.9

#N="E150024856_L01_UDB-64"
#for N in $file; do
#f=$"/home/projects/dtu_00032/analysis/milk-cohort/data/host_removed/${N}.1"
#r=$"/home/projects/dtu_00032/analysis/milk-cohort/data/host_removed/${N}.2"
#echo 'processing' ${N}
#rm -r "bins/megahit/${N}"
#megahit -1 $f -2 $r -o "bins/megahit/${N}"



module purge
module load tools
module load metaphlan/4.1.0

d="/home/projects/dtu_00032/db/metaphlan_db_2307"

file=$(cat logs/stools.txt)
for N in $file; do
    f=$"/home/projects/dtu_00032/analysis/milk-cohort/data/host_removed/${N}.1"
    r=$"/home/projects/dtu_00032/analysis/milk-cohort/data/host_removed/${N}.2"
    echo 
    echo "Running metaphlan with $N"
    echo
    metaphlan "$f","$r" \
        --nproc 40 \
        --input_type fastq \
        --bowtie2db "$d" \
        --bowtie2out "meta/metaphlan/bowtie2/${N}.bowtie2.bz2"  \
        -s "meta/metaphlan/sams/${N}.sam.bz2" \
        -o "/home/projects/dtu_00032/analysis/milk-cohort/meta/metaphlan/original/${N}_profiled.txt"; 

    python /services/tools/metaphlan/4.1.0/bin/sgb_to_gtdb_profile.py  \
            -i "/home/projects/dtu_00032/analysis/milk-cohort/meta/metaphlan/original/${N}_profiled.txt" \
            -o "/home/projects/dtu_00032/analysis/milk-cohort/meta/metaphlan/convGTDB/${N}_profiled.txt" \
            -d "/home/projects/dtu_00032/db/metaphlan_db_2307/mpa_vJun23_CHOCOPhlAnSGB_202307.pkl";
 
done


mail -s 'Stool MP is done' yunso@dtu.dk

exit


#bowtie2-build ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna GRCh38