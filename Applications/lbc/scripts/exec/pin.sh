#!/bin/bash
#
# The execution of this script is as follows:
# mpirun -np # pin.sh /path/to/execution/file
# This script ensures, that all mpi processes are pinned to a physical core
#
#

export KMP_AFFINITY=verbose,scatter
export OMP_NUM_THREADS=4
RANK=${OMPI_COMM_WORLD_RANK:=$PMI_RANK}
if [ $(expr $RANK % 2) = 0  ]
then
     export GOMP_CPU_AFFINITY=0-3
     numactl --preferred=0 --cpunodebind=0 $@
else
     export GOMP_CPU_AFFINITY=4-7
     numactl --preferred=1 --cpunodebind=1 $@
fi
