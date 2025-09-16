#!/bin/bash
#
# tame minimal parameter variation and visualization automated workflow script
# takes param_name and value as command line input (e.g. "changeNsee.sh rmax 4.5")
# by kai wirtz  2025/06/17
# Assign default values to positional arguments if not provided
# param_name is the name of the parameter to be changed in fabm.yaml
# val is the new value for that parameter
# If no parameters are provided, it defaults to "hydrolysis" and "0.21"
param_name=${1:-hydrolysis}
val=${2:-0.21}

thisdir="$PWD"

cp output.yaml output.yaml.backup
#echo $thisdir
# requires fabm-executable ./fabm0d and assumes the ref state to be stored in output_0.nc
# checks whether binary fabm0d and ref result file exist
if ! [[ -x fabm0d && -f output_0.nc ]];
then
    echo "reset output.yaml"
    sed -r 's/output( )*[0-9._-]+/output_0/g' output.yaml > output.yaml.tmp
    mv output.yaml.tmp output.yaml
    echo "create new fabm0d"
    make
fi

# replace old value of parameter in fabm.yaml
replace_y.sh fabm.yaml $param_name $val

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
#MDIR=$(realpath ../../postprocessing)
#cd ${MDIR}
if [[ ! -p matlab_pipe ]]; then
    rm -f matlab_pipe
    mkfifo matlab_pipe
    # MATLAB mit Pipe starten
    echo "Starting MATLAB with pipe matlab_pipe"
    #matlab -nodisplay -nosplash  < matlab_pipe > matlab_output.log 2> matlab_errors.log &
    matlab -nodisplay -nosplash  < matlab_pipe 2>&1 | tee matlab_output.log &
    #matlab -nodisplay -nosplash  -r "try;cmp_0Dres ; catch ME; fprintf('Error: %s\n', ME.message); end; exit" 2>&1 | tee matlab_output.log
    MATLAB_PID=$!
    echo "Matlab process ID is $MATLAB_PID"
    #printf "cd ${MDIR};\n" > matlab_pipe
fi

# Befehle senden
##printf "run('cmp_0Dres.m');\n" > matlab_pipe
echo "cmp_0Dres" > matlab_pipe

# echo "../../postprocessing/cmp_0Dres;" > matlab_pipe
#cd $thisdir

# exit
# matlab -nodisplay -nosplash -r "while true; if exist('tmp.m','file'); run('tmp.m'); delete('tmp.m'); end; pause(0.1); end" &
# cp ../../postprocessing/cmp_0Dres.m tmp.m
# generates an update of simres.png that can be viewed by, e.g., okular or preview

# restore old fabm.yaml
cp fabm.yaml.backup fabm.yaml
cp output.yaml.backup output.yaml

# remove pipe if wanted
# rm ${MDIR}/matlab_pipe
echo "Done. Check simres.png for results."

# Kill matlab if wanted
# killall -q matlab MATLAB
