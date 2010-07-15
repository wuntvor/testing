#!/usr/local/bin/bash
#PBS -q test                         # queue: dq  for <=8 CPUs, multi for more
#PBS -T mpisx                           # job type: mpisx for MPI
#PBS -l cpunum_job=1                   # cpus per node
#PBS -b 1                               # number of nodes
#PBS -l elapstim_req=01:00:00         # max elapse time
#PBS -l cputim_job=01:00:00             # max accumulated cputime
#PBS -l memsz_job=16gb                 # max accumulated memory
#PBS -A hpc43598                  # account code
#PBS -j o                               # join stdout/stderr
#PBS -N lbc                      # job name
#PBS -M hasert@hlrs.de               # email address

# variables for monitoring:
# -------------------------
export F_PROGINF=DETAIL
export F_FTRACE=YES
export MPIPROGINF=YES
export MPIMULTITASKMIX=YES

# export variables for MPI:
# -------------------------
MPIEXPORT="F_RSVTASK F_PROGINF F_FTRACE MPIPROGINF MPIMULTITASKMIX
C_SETBUF F_SETBUF C_SETBUF_VERBOSE"
export MPIEXPORT

NNODES=1
PES_PER_NODE=1
MPI_PROCESSES=1

# Make sure that the stacksize is unlimited in all relevant shells
ulimit -s unlimited

# set your correct enviroment if necessary
#. /opt/modules/3.1.63.1.6/init/bash
#module load mpich-126-pgi/6.2.5 pgi/6.2.5

EXECUTABLE=lbc #vec #trats_sx9_ftrace

#cd /nfs/nas/homeB/home2/HLRS/hlrs/hpcmhase/lbbench/utl/benchmark
cd $PBS_O_WORKDIR


for j in 1 #2 #3

do

for k in 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 30 32 34 40 44 48 56 60 62 64 72 80 88 90 92 96 100 107 114 128 192 256 512 1024

do

  if(($k > 128)); then
LY=128
  else
LY=$k
  fi 
  echo 0.005 > lbc.params
  echo 10 >> lbc.params
  echo 50,50 >> lbc.params
  echo 2 >> lbc.params
  echo $k,$LY,$LY >> lbc.params
  echo f >> lbc.params
  echo 1,1,1 >> lbc.params
  echo f >> lbc.params

#mpirun -np $MPI_PROCESSES ./lbc_ftrace
./lbc
done

done

