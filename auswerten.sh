#!/usr/bin/env bash
result="../../$(echo "${1}" | grep -o "/[a-zA-Z0-9\.]\{1,\}/\?$" | grep -o "[a-zA-Z0-9\.]\{1,\}" ).csv"
cd $1/Step_Component_2um_Ahtermal1_diffNuc_Cath1.5_Boundary0
echo "Folder;JobID;Hostname;Start;End;Diff;Error;Cores;Mem_req;Mem_used;Elementsize;Cache;Success;Additional Info;" > ${result}
for dir in ./*/; do
    #Print Foldername
    line="$(echo ${dir} | grep -o "[0-9_]*\+")"

    outputfile="${dir}$(ls ${dir} | grep -o "output\.[0-9]\{5,\}\.txt")"

    if [[ -f ${outputfile} ]]; then
        #print JobID
        line="${line};$(echo ${outputfile} | grep -o "output\.[0-9]\{5,\}\.txt" | grep -o "[0-9]\{5,\}")"
        #print hostname
        line="${line};$(cat ${outputfile} | grep "hpc\.itc\.rwth-aachen\.de" | grep -o "^[a-z]\{2,4\}[0-9]\{3,5\}" | tail -n 1)"
        #print Start and End Time
        start="$(cat ${outputfile} | grep "\(CEST 2019\)\|\(CET 2019\)" | sed '1q;d')"
        end="$(cat ${outputfile} | grep "\(CEST 2019\)\|\(CET 2019\)" | sed '2q;d')"
        line="${line};${start};${end}"
        timeformat="%Y-%m-%d-%T"
        if [[ ${end} != "" ]]; then
            line="${line};$(ddiff -S -f "%dd %Hh %Mm %Ss" -i "${timeformat}" $(date -d "${start}" "+${timeformat}") $(date -d "${end}" "+${timeformat}"))"
        else
            line="${line};"
        fi
        #Check for Errors
        line="${line};"
        line="${line}$(cat ${outputfile} | grep "killed by signal")"
        line="${line}$(cat ${outputfile} | grep "slurmstepd: error:")"
        #Print #CPUS and mem
        line="${line};$(cat ${outputfile} | grep "#SBATCH --cpus-per-task" | grep -o "[0-9]\+\+")"
        line="${line};$(cat ${outputfile} | grep "#SBATCH --mem" | grep -o "[0-9]\+\+[A-Z]\?")"
        line="${line};$(cat ${outputfile} | grep "memusage: peak usage:" | grep -o "[0-9]\+\+\.[0-9]\+\+ [A-Z]\+\+")"
    else
        #If the outputfile is not found
        line="${line};;;;;;;;;"
    fi
    inputfile="${dir}$(ls ${dir} | grep -o "\([0-9]*\.\)\?InputFile[A-Za-z_]*\.xml")"
    #Print Cellsize and Cache value
    line="${line};$(cat ${inputfile} | grep -o "cell_size_m=\"[0-9]\+e-[0-9]\+\"" | grep -o "[0-9]\+\+e-[0-9]\+\+")"
    line="${line};$(cat ${inputfile} | grep -o "<cache_number_entries>[-0-9]\+</cache_number_entries>" | grep -o "[-0-9]\+\+")"
    #Print Success

    if [[ ! -d "${dir}1/" ]]; then
        line="${line};no"
    else
        numFiles=$(ls -l ${dir}1/ | wc -l)
        numFiles=$((numFiles-1))
        if (( ${numFiles} < 450 )); then
            line="${line};no"
        elif (( ${numFiles} < 500 )); then
            line="${line};probably"
        else
            line="${line};yes"
        fi
    fi
    #Print additional Info
    if [[ -f ${outputfile} ]]; then
        line="${line};$(cat ${outputfile} | grep "####" | grep -o "[^#]*\+")"
        line="${line/\$\{SLURM_ARRAY_TASK_ID\}/$(head -n 3 ${outputfile})}"
        line="${line//$'\n'/ }"
    else
        line="${line};"
    fi
    echo "${line}" >> ${result}
done