#!/bin/bash
# This is for Bluegene queues
# @ job_name = LBC_ADCL_Bench
# @ comment = "LBC benchmark"
# @ error = $(job_name).$(jobid).out
# @ output = $(job_name).$(jobid).out
# @ environment = COPY_ALL
# @ wall_clock_limit = 00:30:00
# @ notification = error
# @ notify_user = m.hasert@grs-sim.de
# @ job_type = bluegene
# @ bg_size = 8 
# @ queue

#
# nehalem
#x#PBS -l nodes=4:ppn=8:nehalem,walltime=6:00:00
#x#PBS -m abe
#x#PBS -M benkert@hlrs.de
#
# sgi
#x#PBS -l select=32 -l walltime=6:00:00 -M benkert@hlrs.de

# cray
#x#PBS -l mppwidth=32,mppnppn=8,walltime=6:00:00
#x#PBS -M m.hasert@grs-sim.de

# -----------------------------------------------------------------------------
# purpose: creates directory structure necessary for performance measurements 
#          with ADCL
#          Copy this file to the location, where the ADCL Benchmark runs 
#          should be run and execute with mode=setup
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# set variables
# -----------------------------------------------------------------------------

mode=setup
#mode=run

# number of MPI processes
proc=32

# number of runs
nrun=3

#onesided=_true_
onesided=_false_


#hpc_type=nehalem
#hpc_type=sgi
#hpc_type=bwgrid
#hpc_type=cray
hpc_type=bluegene



# Intel MPI or OpenMPI? if ompi, use pinning script
#mpi_lib=altix-mpi
#mpi_lib=ompi
mpi_lib=impi
#mpi_lib=mvapich
#mpi_lib=cray

export MPI_TYPE_DEPTH=9

# adcl_home is path with ADCL sources
#adcl_home=/nfs/home2/HLRS/hlrs/hpcmhase/adcl
#adcl_home=~/projects/adcl/adcl_lib_svn_edgar/trunk/      # pcglap10
#adcl_home=/lustre/ws1/ws/hpcbenka-adcl_lbc-0/adcl-svn    # nehalem
adcl_home=/homea/grs100/grs10004/work/adcl/nn_ext         # jugene
#adcl_home=/home/hlrb2/pr47te/lu43pay/adcl-nnext-r3221  # sgi

# WORKDIR is path where directory structure should be build
workdir=$PBS_O_WORKDIR
#workdir=$PWD

# LBC files
#lbc_home=~/projects/skalb/lbc/lbc-svn/                   # pcglap10
#lbc_home=/lustre/ws1/ws/hpcbenka-adcl_lbc-0/lbc-svn      # nehalem
lbc_home=/homea/grs100/grs10004/work/lbc                  # jugene
#lbc_home=/home/hlrb2/pr47te/lu43pay/lbc-r218              # sgi

lbc_exe_with_adcl="$lbc_home/"lbc-with-adcl-$mpi_lib
lbc_exe_no_adcl_subarr="$lbc_home"/lbc-no-adcl-subarr-$mpi_lib
lbc_exe_no_adcl_isir="$lbc_home"/lbc-no-adcl-isir-$mpi_lib
lbc_exe_no_adcl_sr="$lbc_home"/lbc-no-adcl-sr-$mpi_lib

lbc_params="$lbc_home"/lbc.params.$proc

pinning_script="$lbc_home"/scripts/pin.sh


# -----------------------------------------------------------------------------
# no changes below necessary
# -----------------------------------------------------------------------------
# check input parameters
# -----------------------------------------------------------------------------
if [[ "X$workdir" == "X" ]]; then
    echo "No work directory given. Using $PWD"
    workdir=$PWD
fi
if [[ ! -d $workdir ]]; then
    echo "Directory $workdir does not exist. Exiting ..."
    exit
fi

if [[ "X$mode" == "Xsetup" ]]; then
   if [[ ! -d $adcl_home ]]; then
       echo "Directory $adcl_home does not exist. Exiting ..."
       exit
   fi

   if [[ ! -f $lbc_exe_with_adcl ]]; then
       echo "File $lbc_exe_with_adcl does not exist. Exiting ..."
       exit
   fi

   if [[ ! -f $lbc_exe_no_adcl_subarr ]]; then
       echo "File $lbc_no_adcl_subarr does not exist. Exiting ..."
       exit
   fi

   if [[ ! -f $lbc_exe_no_adcl_isir ]]; then
       echo "File $lbc_exe_no_adcl_isir does not exist. Exiting ..."
       exit
   fi

   if [[ ! -f $lbc_exe_no_adcl_sr ]]; then
       echo "File $lbc_exe_no_adcl_sr does not exist. Exiting ..."
       exit
   fi

   if [[ ! -f $lbc_params ]]; then
       echo "File $lbc_params does not exist. Exiting ..."
       exit
   fi

         mkdir exe 
         mkdir config

fi

# -----------------------------------------------------------------------------
# copy necessary ADCL files to ADCL_PATH
# -----------------------------------------------------------------------------
if [[ "X$mode" == "Xsetup" ]]; then
   cp "$adcl_home"/config/config.adcl.brute "$adcl_home"/config/config.adcl.hypo    \
      "$adcl_home"/config/config.adcl.*pair*  "$adcl_home"/config/config.adcl.*aao* \
      config/
   cp  $lbc_exe_with_adcl $lbc_exe_no_adcl_subarr $lbc_exe_no_adcl_isir $lbc_exe_no_adcl_sr exe/
   cp $lbc_params ./exe/lbc.params
fi

comm="no-adcl-isir no-adcl-sr no-adcl-subarr brute hypo  "
comm="$comm IsIr_aao SIr_aao IsIr_aao_pack SIr_aao_pack IsIr_pair SIr_pair"
comm="$comm IsIr_pair_pack SIr_pair_pack S_R_pair Sr_pair S_R_pair_pack Sr_pair_pack"
comm="$comm IsIr_aao_buf SIr_aao_buf IsIr_pair_buf SIr_pair_buf S_R_pair_buf Sr_pair_buf"
#------------------------------------
# Add one-sided Communication?

if [[ "X$onesided" == "X_true_" ]]; then
   echo "One-sided Communication is also performed"
   comm="$comm PostStartGet_aao PostStartGet_pair PostStartPut_aao PostStartPut_pair"
   comm="$comm WinFenceGet_aao WinFenceGet_pair WinFencePut_aao WinFencePut_pair"
fi

date

cd $workdir

if [[ "X$mode" == "Xrun" ]]; then
    if [[ "X$hpc_type" == "Xnehalem" ]]; then
        module unload mpi
        module unload compiler
        # module load mpi/openmpi/1.4-gnu-4.1.2
        ulimit -s unlimited
        module load compiler/intel/11.0
        if [[ "$mpi_lib" == "impi"  ]]; then
           module load mpi/impi/intel-11.0.074-impi-4.0.0.017
        else
           module load mpi/openmpi/1.4-intel-11.0
        fi
    elif [[ "X$hpc_type" == "Xsgi" ]]; then
        . /etc/profile.d/modules.sh

        module unload mpi.altix
        module unload ccomp
        module unload fortran
        module load ccomp/intel/11.0
        module load fortran/intel/11.0
        if [[ "$mpi_lib" == "altix-mpi"  ]]; then
           module load mpi.altix
        elif [[ "$mpi_lib" == "impi"  ]]; then
           module load mpi.intel/3.2
        else
           module load mpi.ompi
        fi
    elif [[ "X$hpc_type" == "Xbluegene" ]]; then
           module load UNITE
     fi
fi

irun=1
while [[ irun -le $nrun ]]
do
    run="$mpi_lib"_"$proc"_run"$irun"

    if [[ "X$mode" == "Xsetup" ]]; then
         mkdir $run
    fi
    #if [[ "X$?" == "X1" ]];  then
    #   echo "Creation of directory $run failed. Exiting ..." && exit
    #fi
    cd $run
    if [[ "X$?" == "X1" ]];  then
       echo "Could not enter directory $run. Exiting ..."    && exit
    fi
    echo "entering $run ..."

    for method in `echo $comm`
    do
        if [[ "X$mode" == "Xsetup" ]]; then
            mkdir $method
        fi
        #if [[ "X$?" == "X1" ]];  then
        #   echo "Creation of directory $method failed. Exiting ..." && exit
        #fi
        cd "$method"
        if [[ "X$?" == "X1" ]];  then
           echo "Could not enter directory $method. Exiting ..."    && exit
        fi
        echo "entering $method ..."

        if [[ "X$method" = "Xno-adcl-subarr" ]]; then
                lbc_exe=./../../exe/`basename $lbc_exe_no_adcl_subarr`
        elif [[ "X$method" = "Xno-adcl-isir" ]]; then
                lbc_exe=./../../exe/`basename $lbc_exe_no_adcl_isir`
        elif [[ "X$method" = "Xno-adcl-sr" ]]; then
                lbc_exe=./../../exe/`basename $lbc_exe_no_adcl_sr`
        else
                lbc_exe=./../../exe/`basename $lbc_exe_with_adcl`
        fi

        # set execution file
        if [[ "X$mode" == "Xsetup" ]]; then
            # create links
            # - executable
            if [[ ! -f $lbc_exe ]]
            then
                echo "File $lbc_exe does not exist. Exiting ..."
                exit
            fi
            ln -s $lbc_exe

            # - adcl config file
            if [[ "X$method" != "Xno-adcl-subarr" && "X$method" != "Xno-adcl-isir" && \
                  "X$method" != "Xno-adcl-sr" ]]; then
                adclconfig_file=./../../config/config.adcl."$method"
                if [[ ! -f $adclconfig_file ]]; then
                    echo "File $adclconfig_file does not exist. Exiting ..."
                    exit
                fi
                ln -s $adclconfig_file ./config.adcl
            fi
            # - lbc param file
            lbcparam_file=./../../exe/lbc.params
            if [[ ! -f $lbcparam_file ]]; then
                echo "File $lbcparam_file does not exist. Exiting ..."
                exit
            fi
            ln -s $lbcparam_file ./lbc.params


        elif [[ "X$mode" == "Xrun" ]]; then
           # execute program

           if [[ "X$hpc_type" == "Xnehalem" ]]; then
               if [[ "X$mpi_lib" == "Ximpi"  ]]; then
                  cat $PBS_NODEFILE > mpd.hosts
                  mpirun -np $proc -r ssh $lbc_exe
                else
                  mpirun -np $proc $pinning_script $lbc_exe
                fi
            elif [[ "X$hpc_type" == "Xsgi" ]]; then
                pwd
                echo mpiexec -n $proc $lbc_exe
                mpiexec -n $proc $lbc_exe
            elif [[ "X$hpc_type" == "Xbluegene" ]]; then
                echo mpiexec -n $proc $lbc_exe
                mpirun -mode VN -np $proc -exe $lbc_exe 
            elif [[ "X$hpc_type" == "Xcray" ]]; then
                echo aprun -n $proc $lbc_exe
                aprun -n $proc $lbc_exe 
            else
                echo "Machine type missing for correct execution."
                exit
            fi
            cp tstep* tstep
        fi

        cd ..
        echo "leaving $method..."
    done

    cd ..
    echo "leaving $run ..."
    irun=$(($irun+1))
done

date


