#!/bin/bash
#PBS -l nodes=1:ppn=8:nehalem
#PBS -l walltime=1:00:00
#PBS -m abe
#PBS -M hasert@hlrs.de
#PBS -N LBC

NPROCX=2
NPROCY=2
NPROCZ=2
NNODES=1
PES_PER_NODE=$(($NPROCX*$NPROCY*$NPROCZ))
MPI_PROCESSES=$(($NNODES*$NPROCX*$NPROCY*$NPROCZ))

# Make sure that the stacksize is unlimited in all relevant shells
ulimit -s unlimited

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

for j in 1 #2 #3

do

for k in 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 30 32 34 40 44 48 56 60 62 64 72 80 88 90 92 96 100 107 114 128 #192 #256 512 1024

do

  if(($k > 128)); then
LY=128
LZ=128
  else
LY=$k
LZ=$k

  fi 
k=$(($k*$NPROCX))
LY=$(($LY*$NPROCY))
LZ=$(($LZ*$NPROCZ))
  echo 0.05 > lbc.params
  echo 10 >> lbc.params
  echo 50,50 >> lbc.params
  echo 2 >> lbc.params
  echo $k,$LY,$LZ >> lbc.params
  echo f >> lbc.params
  echo $NPROCX,$NPROCY,$NPROCZ >> lbc.params
  echo f >> lbc.params

 # cp param.cfg_$k param.cfg
#  mpirun -np $MPI_PROCESSES -machinefile machines ./$EXECUTABLE --pe_affinity pes_per_node $PES_PER_NODE which_pe 1
# tut nicht mpirun -nn $NNODES -nnp $PES_PER_NODE ./$EXECUTABLE > output 2>&1


#mpirun -np $MPI_PROCESSES -hostfile $PBS_NODEFILE ./$EXECUTABLE"_"$i
mpirun -np $MPI_PROCESSES ./lbc
done
done

mv performance.res ${MPI_PROCESSES}_performance.res

#done

#done

#done
