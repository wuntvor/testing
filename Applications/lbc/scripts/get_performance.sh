#!/bin/bash

comm="no-adcl-sr no-adcl-isir no-adcl-subarr brute hypo  "
comm="$comm IsIr_aao IsIr_pair IsIr_aao_pack IsIr_pair_pack SIr_aao SIr_pair SIr_aao_pack SIr_pair_pack "
comm="$comm S_R_pair Sr_pair S_R_pair_pack Sr_pair_pack"
comm="$comm IsIr_aao_buf SIr_aao_buf IsIr_pair_buf SIr_pair_buf S_R_pair_buf Sr_pair_buf"
comm="$comm WinFencePut_aao WinFenceGet_aao PostStartPut_aao PostStartGet_aao"
comm="$comm WinFencePut_pair WinFenceGet_pair PostStartPut_pair PostStartGet_pair"


rundirs=`find . -name "*run*"`
#irun=1
#while [[ irun -le $nrun ]]
#do
#    run="$mpi_lib"_"$proc"_run"$irun"
#
#    cd $run
#    if [[ "X$?" == "X1" ]];  then
#       echo "Could not enter directory $run. Exiting ..."    && exit
#    fi
#    echo "entering $run ..."

    fixeddir=`echo $rundirs | cut -d" " -f1`/no-adcl-sr
    echo n 
    string=`cat $fixeddir/performance.res | grep ijk`
    echo $string | cut -d" " -f7


    echo lx ly lz 
    echo $string | cut -d" " -f4-6

for method in `echo $comm`
do
    printf "%s," $method

    for dir in `echo $rundirs`
    do
       string=`cat $dir/$method/performance.res | grep ijk`
       t_tot=`echo $string | cut -d" " -f8`
       t_comm=`echo $string | cut -d" " -f9`   
       perc=`echo $string | cut -d" " -f10`   
       printf "%s %s %s ", $t_tot, $t_comm, $perc 
    done
    printf "\n"
done 

#   irun=$(($irun+1))
#done

