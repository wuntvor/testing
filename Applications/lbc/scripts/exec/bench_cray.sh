#!/bin/bash
#PBS -l mppwidth=8,mppnppn=8
#PBS -l walltime=1:00:00
#PBS -m abe
#PBS -M hasert@hlrs.de
#PBS -N LBC

NPROCX=1
NPROCY=1
NPROCZ=1
NNODEX=1
NNODEY=1
NNODEZ=1
NNODES=$(($NNODEX*$NNODEY*$NNODEZ))
PES_PER_NODE=$(($NPROCX*$NPROCY*$NPROCZ))
MPI_PROCESSES=$(($NNODES*$NPROCX*$NPROCY*$NPROCZ))

export XT_SYMMETRIC_HEAP_SIZE=2G

# Make sure that the stacksize is unlimited in all relevant shells
ulimit -s unlimited

#module unload PrgEnv-pgi
#module load PrgEnv-cray

# set your correct enviroment if necessary
#. /opt/modules/3.1.63.1.6/init/bash
#module load mpich-126-pgi/6.2.5 pgi/6.2.5

EXECUTABLE=lbc

#cd /nfs/HOME/HLRS/hlrs/hpcmhase/lbbench/trunk/utl/benchmark
cd $PBS_O_WORKDIR
#cat $PBS_NODEFILE > nodefile.txt
loop=0
loopend=$(( $PES_PER_NODE - 1 ))
LIST_NODE_FILES=$PBS_NODEFILE
while [ $loop -ne $loopend ]
do
    loop=$(( $loop + 1 ))
    LIST_NODE_FILES="$LIST_NODE_FILES"" $PBS_NODEFILE"
done
#cat $LIST_NODE_FILES | sort > machines
#cat $PBS_NODEFILE | sort > hf


for LX in 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 30 32 34 40 44 48 56 60 62 64 72 80 88 90 92 96 100 107 114 128 192 256 512 1024

do

if(($LX > 128)); then
LY=128
LZ=128
else
LY=$LX
LZ=$LX

fi

  echo 0.005 > lbc.params
  echo 10 >> lbc.params
  echo 50,50 >> lbc.params
  echo 2 >> lbc.params
  echo $(($LX*$NPROCX*$NNODEX)),$(($LY*$NPROCY*$NNODEY)),$(($LZ*$NPROCZ*$NNODEZ)) >> lbc.params
  echo f >> lbc.params
  echo $(($NPROCX*$NNODEX)),$(($NPROCY*$NNODEY)),$(($NPROCZ*$NNODEZ)) >> lbc.params
  echo f >> lbc.params
sleep 2
aprun -n $MPI_PROCESSES -N $(($NPROCX*$NPROCY*$NPROCZ)) ./lbc
done
