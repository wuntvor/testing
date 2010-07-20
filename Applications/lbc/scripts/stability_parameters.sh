#!/bin/bash
#PBS -l nodes=1:ppn=8:nehalem
#PBS -l walltime=6:00:00
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

 cd $PBS_O_WORKDIR
cat $PBS_NODEFILE > mpd.hosts 

for executable in lbc_mrt_les lbc_bgk_les; do

for bulk in 0.001 0.01 0.1 0.5 1.0 2.0 ;do

for umax in 0.0001 0.001 0.01 0.05;do

for LX in 10 30 60 120;do

LY=$LX
LZ=$LX

  echo $umax  > lbc.params            
  echo $omega >> lbc.params            
  echo 10000,10000,1,1,-10  >> lbc.params            
  echo 2      >> lbc.params                           
  echo $LX,$LY,$LZ   >> lbc.params                   
  echo f      >> lbc.params                           
  echo $NPROCX,$NPROCY,$NPROCZ   >> lbc.params                       
  echo f      >> lbc.params                          
  echo 1.0    >> lbc.params                         
  echo $bulk  >> lbc.params                        

mpirun -np $MPI_PROCESSES -r ssh ./$executable

done
done
done
