#!/usr/local_rwth/bin/zsh

#SBATCH --job-name=Sphaerosim
#SBATCH --mem=0G
#SBATCH --time=0-02:30:00

#SBATCH --array=1-14

#OpenMP settings:
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --ntasks-per-node=1

#SBATCH -v
#SBATCH --hint=nomultithread

#SBATCH --export=NONE
#SBATCH --output=/work/xu077052/Projekte/sphaerosim/out/test/output.%J.txt

#export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

### Change to the work directory
cd /work/xu077052/Projekte/sphaerosim/bin
mkdir ${SLURM_JOB_ID}
cd ${SLURM_JOB_ID}
echo "Array ID:${SLURM_ARRAY_TASK_ID}"

#combinations=(
#"intel;"
#"intel;-inline-forceinline"
#"intel;-xHOST"
#"intel;-ipo"
#"intel;-fast"
#"intel;-inline-forceinline -fast"
#"intel;-inline-forceinline -xHost"
#"gcc/8;"
#"gcc/8;-Ofast"
#"gcc/8;-march=native"
#"gcc/8;-Ofast -march=native"
#"clang;"
#)
#IFS=';' read -rA OUT <<< "${combinations[SLURM_ARRAY_TASK_ID]}"
#echo "${OUT[1]}"
#echo "${OUT[2]}"
#module switch intel "${OUT[1]}"
export OMP_NUM_THREADS="${combinations[$SLURM_ARRAY_TASK_ID]}"
echo "Threads: ${OMP_NUM_THREADS}"
echo ""

RunType="release"
export OUTPUT_DIR="/rwthfs/rz/cluster/hpcwork/xu07705a/Projekte/sphaerosim/out/test"


if [[ "$RunType" = "likwid" ]]; then
    CmakeFlags="-DLIKWID_PERFMON=1 -DCMAKE_BUILD_TYPE=release"
    RunPrefix="likwid-perfctr -m -C S0:0-23 -g DATA"
    module load likwid
elif [[ "$RunType" = "release" ]]; then
    CmakeFlags="-DCMAKE_BUILD_TYPE=release"
    RunPrefix=""
    export OMP_PLACES=cores
    export OMP_PROC_BIND=close
elif [[ "$RunType" = "debug" ]]; then
    CmakeFlags="-DCMAKE_BUILD_TYPE=debug"
    RunPrefix=""
    export OMP_PLACES=cores
    export OMP_PROC_BIND=close
elif [[ "$RunType" = "memusage" ]]; then
    CmakeFlags="-DCMAKE_BUILD_TYPE=release"
    RunPrefix="cgroup_memusage -i=10 --step=6"
    #export OMP_PLACES="{$(numactl --hardware | grep "node [0-9] cpus:" | grep -o ":\( [0-9]\{1,3\}\)*" | grep -o "[0-9][ 0-9]*" | sed -r 's/[ ]+/,/g' | paste -sd "#" | sed -r 's/[#]+/},{/g')}"
    export OMP_PLACES=cores
    export OMP_PROC_BIND=close
elif [[ "$RunType" = "numamem" ]]; then
    CmakeFlags="-DCMAKE_BUILD_TYPE=release"
    RunPrefix="no_numa_balancing numamem -s 60"
    export OMP_PLACES=cores
    export OMP_PROC_BIND=close
elif [[ "$RunType" = "PerfReport" ]]; then
    CmakeFlags="-DCMAKE_BUILD_TYPE=release"
    RunPrefix="perf-report --nompi"
    export OMP_PLACES=cores
    export OMP_PROC_BIND=close
    module load reports
elif [[ "$RunType" = "VTune tr" ]]; then
    CmakeFlags="-DDebug=1 -DCMAKE_BUILD_TYPE=release"
    RunPrefix="/rwthfs/rz/SW/intel/vtune/XE2019-u02/vtune_amplifier_2019.2.0.584348/bin64/amplxe-cl -collect threading -data-limit=1024 -app-working-dir $(pwd) -result-dir /work/xu077052/Projekte/sphaerosim/VTune/Sphaerosim/${SLURM_JOB_ID}tr/ --"
    export OMP_PLACES=cores
    export OMP_PROC_BIND=close
elif [[ "$RunType" = "VTune hs" ]]; then
    CmakeFlags="-DDebug=1 -DCMAKE_BUILD_TYPE=release"
    RunPrefix="/rwthfs/rz/SW/intel/vtune/XE2019-u02/vtune_amplifier_2019.2.0.584348/bin64/amplxe-cl -collect hotspots -data-limit=1024 -app-working-dir $(pwd) -result-dir /work/xu077052/Projekte/sphaerosim/VTune/Sphaerosim/${SLURM_JOB_ID}hs/ --"
    export OMP_PLACES=cores
    export OMP_PROC_BIND=close
elif [[ "$RunType" = "Inspector" ]]; then
    CmakeFlags="-DDebug=1 -DCMAKE_BUILD_TYPE=release"
    RunPrefix="inspxe-cl -collect mi2 -knob detect-uninit-read=false -knob revert-uninit=false -knob detect-leaks-on-exit=true -knob enable-memory-growth-detection=true -knob enable-on-demand-leak-detection=true -knob still-allocated-memory=true -knob stack-depth=8 -module-filter-mode=include -appdebug=off -no-trace-mpi -app-working-dir $(pwd) -result-dir /work/xu077052/Projekte/sphaerosim/Inspector/Sphaerosim/${SLURM_JOB_ID}/ --"
    export OMP_PLACES=cores
    export OMP_PROC_BIND=close
    module load intelixe
elif [[ "$RunType" = "ScoreP" ]]; then
    CmakeFlags="-DScoreP=1 -DCMAKE_BUILD_TYPE=release"
    RunPrefix=""
    export OMP_PLACES=cores
    export OMP_PROC_BIND=close
    module load DEV_TOOLS scorep
    export SCOREP_TOTAL_MEMORY=4GB
    export SCOREP_ENABLE_TRACING=false
    export SCOREP_PROFILING_MAX_CALLPATH_DEPTH=20
elif [[ "$RunType" = "ThreadSanitizer" ]]; then
    CmakeFlags="-DTHREAD_SANITIZER=1 -DCMAKE_BUILD_TYPE=release"
    RunPrefix=""
    export OMP_PLACES=cores
    export OMP_PROC_BIND=close
    module switch intel clang/7.0
else
    CmakeFlags=""
    RunPrefix=""
    export OMP_PLACES=cores
    export OMP_PROC_BIND=close
fi

### Execute your application
cat ../../run_SLURM_array.sh
echo "\n\n"
~/Dokumente/Übungen/Laufzeittest/testUlimitPin.sh
echo "\n"
source ~/Dokumente/Übungen/Laufzeittest/.printGeneralInformation.sh

~/Downloads/cmake-3.13.1-Linux-x86_64/bin/cmake ../.. ${CmakeFlags}
make -j
date
${RunPrefix} ./Sphaerosim -f /work/xu077052/Projekte/sphaerosim/InputFile_Linux.xml
date
#### Skalierbarkeitstest ${SLURM_ARRAY_TASK_ID}
