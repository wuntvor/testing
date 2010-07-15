#! /bin/bash
echo "This is the output Script that generates the environment variables and runs tecplot"
basefilename=$1  #"d3q19_Re10.0_Lx_30_u0.005_"
if [ "X${basefilename}" = "X" ]
then
echo "Error: No basefilename given. Please pass the Basefile-Name as an argument"
exit
fi

max_loops=10000
i=0
while [ $i -lt $max_loops ]
do 
testfile="${basefilename}${i}.plt"
if [ -f $testfile ]
then
j=1
else
break
fi
i=$(($i+1))
done
if [ $i -eq  0 ]; then
   echo "No files found"
exit
fi

echo "Found $i files! Setting environment variables..."
export Tec_Num_Files=$i
export Tec_Output_File=$basefilename
echo $Tec_Num_Files
echo $Tec_Output_File

tec360 -p tec_out.mcr
