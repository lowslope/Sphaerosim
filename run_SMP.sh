#!/usr/bin/env zsh
 
### Job name
#BSUB -J sphaerosim_ErsterLauf
 
### File / path where STDOUT & STDERR will be written
###    %J is the job ID, %I is the array ID
#BSUB -o sphaerosim_ErsterLauf.%J.%I
 
### Request the time you need for execution in minutes
### The format for the parameter is: [hour:]minute,
### that means for 80 minutes you could also use this: 1:20
#BSUB -W 5:0
 
### Request memory you need for your job in TOTAL in MB
#BSUB -M 307200
 
### Request the number of compute slots you want to use
#BSUB -n 144
 
### Use esub for OpenMP/shared memeory jobs
#BSUB -a openmp
#BSUB -P jara0001
 
### Change to the work directory
cd /work/xu077052/Projekte/sphaerosim/test
 
### Execute your application
~/Dokumente/Übungen/Laufzeittest/.printGeneralInformation.sh

rm * -r
~/Downloads/cmake-3.13.1-Linux-x86_64/bin/cmake ..
make
r_memusage -i=10 --step=6 ../run.sh