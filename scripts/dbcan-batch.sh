#!/bin/bash

cd /home/projects/dtu_00032/analysis/milk-cohort
scriptpath="bins/new-ani98/dbcan-scripts"
dbcanpath="bins/new-ani98/bakta-dbcan"
mkdir -p $scriptpath
mkdir -p $dbcanpath


for bin in bins/new-ani98/bakta/*; do
    b=$(basename $bin)
    script_name="${scriptpath}/run_${b}.sh"
    log_name="${scriptpath}/${b}.log"
    err_name="${scriptpath}/${b}.err"

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
    echo "module load signalp/6.0h-slow_sequential  " >> $script_name
    echo "module load dbcan/4.1.4" >> $script_name
    echo "cd /home/projects/dtu_00032/analysis/milk-cohort" >> $script_name
    echo "" >> $script_name
    echo "rm -r $dbcanpath/$b" >> $script_name
    echo "mkdir -p $dbcanpath/$b" >> $script_name
    echo "run_dbcan \\" >> $script_name
    echo "    $bin/$b.faa protein \\" >> $script_name
    echo "    -c $bin/$b.gff3 \\" >> $script_name
    echo "    --tools all \\" >> $script_name
    echo "    --dia_cpu 40 \\" >> $script_name
    echo "    --hmm_cpu 40 \\" >> $script_name
    echo "    --tf_cpu 40 \\" >> $script_name
    echo "    --stp_cpu 40 \\" >> $script_name
    echo "    --cgc_substrate \\" >> $script_name
    echo "    --db_dir /home/projects/dtu_00032/db/dbcan_db_09072024 \\" >> $script_name
    echo "    --out_dir $dbcanpath/$b; \\" >> $script_name


    chmod +x $script_name
    qsub $script_name
done
