#!/bin/bash
# -----------------------------------------------------------------------------
# purpose: parse directories for gnuplot output 
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# set variables
# -----------------------------------------------------------------------------

#mode=setup
echo $1
if [ X$1 = X ]; then
echo "Error: No directory to parse given"
echo "Call with directory in argument"
exit
else
echo "Parsing through the directories $1"
      if [ -d $1 ]; then
         echo "Found directory"
         run=$1 
      else
         echo "Error: directory not found. exiting" 
         exit
      fi
fi
mode=run
gplfile=plot_adcl_tstep.gpl
echo "# this is a gnuplot file created by parsing runs of adcl-benchmarks" > $gplfile 
echo "set key under" >> $gplfile
echo "set xlabel \"Timestep\"  " >> $gplfile

countdir=0
date


    # count how many subdirs
    for subdir in `ls $run`
    do
      countdir=$(($countdir + 1))
    done
    echo "I have located $countdir subdirs"

    # parse through all subdirs  
    for plot in single total; do
   counter=0
    echo "this is for plot $plot "
    if [ $plot = "single" ]; then
      echo "set ylabel \"Duration for one timestep [ms]\"  " >> $gplfile
      echo "set title  \"LBC Calculation time for each Timestep for $run\"  " >> $gplfile
      col=2
    else
      echo "set ylabel \"Accumulated duration [ms]\"  " >> $gplfile
      echo "set title  \"LBC Accumulated Calculation time for $run\"  " >> $gplfile
      col=3
    fi
    for subdir in `ls $run`; do
      counter=$(($counter+1))
    if [ $subdir = "no-adcl" ]; then
      linetype="with linespoints"
    else
    if [[ $subdir = *pair* ]]; then
      linetype="with points"
    else
      linetype="with lines"
    fi
    fi
      if [ $counter -eq 1 ]; then
         echo "plot \"$run/$subdir/tstep\" using 1:$col $linetype title \"$subdir\",\\" >> $gplfile 
      else 
      if [ $counter -lt $countdir ]; then 
      if [ -e $run/$subdir/tstep ]; then
         echo "\"$run/$subdir/tstep\" using 1:$col $linetype title \"$subdir\",\\" >> $gplfile 
      fi
      else
      if [ -e $run/$subdir/tstep ]; then
         echo "\"$run/$subdir/tstep\" using 1:$col $linetype  title \"$subdir\"" >> $gplfile 
      fi
      fi
      fi
    done

   filename="${run}_$plot"
   echo "pause -1"  >> $gplfile 
   echo ""  >> $gplfile 
   echo "# Save to disk: png"  >> $gplfile 
   echo "set terminal png" >> $gplfile 
   echo "set output \"$filename.png\"" >> $gplfile 
   echo "replot" >> $gplfile 

   echo "# Save to disk: ps "  >> $gplfile 
   echo "set size 1.0, 0.6" >> $gplfile 
   echo " set terminal postscript portrait enhanced mono dashed lw 1 \"Helvetica\" 14" >> $gplfile 
   echo " set output \"$filename.ps\"" >> $gplfile 
   echo " replot " >> $gplfile 
   echo " set terminal x11" >> $gplfile 
   echo " set size 1,1 " >> $gplfile 
   echo " " >> $gplfile 
   echo " " >> $gplfile 

    done
    exit



