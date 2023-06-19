#!/usr/bin/env zsh

date
OMP_NUM_THREADS=48 OMP_PLACES=cores OMP_PROC_BIND=close /work/xu077052/Projekte/sphaerosim/bin/test/Sphaerosim -f /work/xu077052/Projekte/sphaerosim/InputFile_Linux.xml
#OMP_NUM_THREADS=24 OMP_PLACES=cores OMP_PROC_BIND=close SCOREP_PROFILING_MAX_CALLPATH_DEPTH=20 SCOREP_ENABLE_TRACING=false SCOREP_TOTAL_MEMORY=10G cgroup_memusage -i 1 --step=1 /work/xu077052/Projekte/sphaerosim/bin/test/src/Sphaerosim -f /work/xu077052/Projekte/sphaerosim/InputFile_Linux.xml
#OMP_NUM_THREADS=24 likwid-perfctr -m -C S0:0-23 -g DATA /work/xu077052/Projekte/sphaerosim/bin/test/src/Sphaerosim -f /work/xu077052/Projekte/sphaerosim/InputFile_Linux.xml
date
