#!/bin/bash
#PBS -l nodes=64:ppn=8:nehalem
#PBS -l walltime=1:00:00
#PBS -m abe
#PBS -M hasert@hlrs.de
#PBS -N LBC

# Make sure that the stacksize is unlimited in all relevant shells
ulimit -s unlimited

# Set domain size 
# Weak scaling:   domain size per process
# Strong scaling: total domain size
NX=100

# Set scaling
SCALING=strong
SCALING=weak

EXECUTABLE=lbc

PARAM_FILE="lbc.params"

cd $PBS_O_WORKDIR

loop=0
loopend=$(( $PES_PER_NODE - 1 ))
LIST_NODE_FILES=$PBS_NODEFILE
while [ $loop -ne $loopend ]
do
    loop=$(( $loop + 1 ))
    LIST_NODE_FILES="$LIST_NODE_FILES"" $PBS_NODEFILE"
done

for NPROCZ in 1 2 4 8
do
for NPROCY in 1 2 4 8
do
for NPROCX in 1 2 4 8
do

if [ $SCALING="weak"   ]; then 
  LX=$(($NX*$NPROCX))
  LY=$(($NX*$NPROCY))
  LZ=$(($NX*$NPROCZ))
else
  LX=$(($NX))
  LY=$(($NX))
  LZ=$(($NX))
fi

  echo 0.001 > $PARAM_FILE            
  echo 1.80  >> $PARAM_FILE            
  echo 50,50,1,1,-10  >> $PARAM_FILE            
  echo 2   >> $PARAM_FILE                           
  echo $LX,$LY,$LZ   >> $PARAM_FILE                   
  echo f   >> $PARAM_FILE                           
  echo $NPROCX,$NPROCY,$NPROCZ   >> $PARAM_FILE                       
  echo f    >> $PARAM_FILE                          
  echo 1.0   >> $PARAM_FILE                         
  echo 0.01   >> $PARAM_FILE                        

# cp param.cfg_$k param.cfg
#  mpirun -np $MPI_PROCESSES -machinefile machines ./$EXECUTABLE --pe_affinity pes_per_node $PES_PER_NODE which_pe 1
# tut nicht mpirun -nn $NNODES -nnp $PES_PER_NODE ./$EXECUTABLE > output 2>&1

MPI_PROCESSES=$(($NPROCX*$NPROCY*$NPROCZ))

echo "# Scaling $SCALING   " >> performance.res
mpirun -np $MPI_PROCESSES ./lbc

done
done
done

