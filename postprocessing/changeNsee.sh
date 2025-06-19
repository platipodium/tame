#!/bin/bash
# tame minimal parameter variation and visualization
# automated workflow script
# takes param_name and value as command line input (e.g. "changeNsee.sh rmax 4.5")
# by kai wirtz  2025/06/17
param_name=$1
val=$2
shift
# sets arguments if not given in command line
# parameter name
if [ -z "$param_name" ]
then
 param_name='hydrolysis'
fi
# new value
if [ -z "$val" ]
then
  val='0.21'
fi
#thisdir=\`pwd\`
thisdir=`pwd`  
#echo $thisdir
# requires fabm-executable ./fabm0d and assumes the ref state to be stored in output_0.nc
# checks whether binary fabm0d and ref result file exist
# cp output.yaml output.yaml.0

if ! [[ -f fabm0d && -f output_0.nc ]]; 
then
    echo "reset output.yaml"
    sed -r 's/output( )*[0-9._-]+/output_0/g' output.yaml > output.yaml.tmp
    mv output.yaml.tmp output.yaml
    echo "create new fabm0d"
    make 
fi
#   ./fabm0d > /dev/null 2 > tmp.log

# replace old value of parameter in fabm.yaml
../../postprocessing/replace_y.sh fabm.yaml $param_name $val

# replace FABM result file name in output.yaml
sed -r 's/output( )*[0-9._-]+/output_1/g' output.yaml > output.yaml.tmp
mv output.yaml.tmp output.yaml

# TODO: entry in control file
#  echo $p1 ' ' $dip $din $dil $par $co2 >>  $varfile

# run fabm0d, sending shell output to log file
echo "Executing ./fabm0d"
./fabm0d > /dev/null 2>tmp.log

# visualize by sending matlab script to opened pipe

# Named Pipe erstellen
cd ../../postprocessing/
if ! [[ -f matlab_pipe ]];
then
    mkfifo matlab_pipe
    # MATLAB mit Pipe starten
    matlab -nodisplay -nosplash < matlab_pipe &
fi
# Befehle senden
echo "run('cmp_0Dres.m');" > matlab_pipe
# echo "../../postprocessing/cmp_0Dres;" > matlab_pipe
cd $thisdir

# matlab -nodisplay -nosplash -r "while true; if exist('tmp.m','file'); run('tmp.m'); delete('tmp.m'); end; pause(0.1); end" &
# cp ../../postprocessing/cmp_0Dres.m tmp.m
# generates an update of simres.png that can be viewed by, e.g., okular or preview
 
# restore old fabm.yaml
cp fabm.yaml.backup fabm.

# remove pipe if wanted
# rm ../../matlab_pipe