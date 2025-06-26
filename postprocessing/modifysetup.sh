## Feedback
# - Call it and plot the figures in R
# - Code it so it adds the module, instead of removing it
# - sed does not work on other OS
# - can just turn true / false in the fabm.yaml 'use' part in a module
# - chek jq to work with yaml files, more adapted than bash

## Script to run a reference setup (saved as the fabm.yaml) and a modified setup, with one of the modules removed.
## To use, call the script and give :
## - The path to the fabm.yaml within tame (e.g. setup/0d/fabm.yaml)
## - The module to remove for the modified setup (e.g. zoo)
# --> bash modifysetup.sh setup/0d/fabm.yaml zoo
# The bash script directory should be postprocessing
echo '>> START'

set -e # Exit on any error

YAML_FILE="$1"  # YAML file to modify
MODULE="$2"     # Paragraph to remove
DIR=`pwd`
#echo $DIR

TAME="$(dirname "$(pwd)")" # TAME main directory
# --> This script assumes that the ./fabm0d is in /setup/0d/
#echo $TAME

cd $TAME
#echo `pwd`

echo '> Copying a backup of '${YAML_FILE}' into '${YAML_FILE}'.backup'
cp $YAML_FILE "${YAML_FILE}.backup" # Saves the full coupling setup .yaml 

## Modify the setup by removing a module from the yaml file

# Removes the paragraph correspondig to $MODULE
rmmodule () {
	sed -n -i '
	    /^[[:space:]]\+\'$MODULE':$/ {
		n
		:c
		    /^[[:space:]]\+\//! {
		        n
		        bc
		    }
	    }
	    p
	' $YAML_FILE
}

echo '> Removing the module '${MODULE}''
rmmodule # Updates the fabm.yaml with the modified setup 

## Makes sure that all the variables will be in the output file
echo '> Copying a backup of the output.yaml into output.yaml.backup'
cp "${TAME}/$(dirname "${YAML_FILE}")/output.yaml" "${TAME}/$(dirname "${YAML_FILE}")/output.yaml.backup" # Saves the full coupling setup .yaml 
sed -i '/      - source: /d' "${TAME}/$(dirname "${YAML_FILE}")/output.yaml"     # To rewrite the output
echo "      - source: *" >> "${TAME}/$(dirname "${YAML_FILE}")/output.yaml"      # Record all the variables
sed -i 's|output|output_modif|g' "${TAME}/$(dirname "${YAML_FILE}")/output.yaml" # To rewrite the output namefile

echo "> Running the model with removed ${MODULE}"
cd setup/0d
./fabm0d # Run the model

## Restore the backup to the original yaml file after run
cd ${TAME}
#echo `pwd`

echo "> Running the model with reference setup"
cp "${YAML_FILE}.backup" $YAML_FILE
echo '> '"${YAML_FILE}"' restored'

# Rename the output file for the reference run to output.nc
sed -i 's|output_modif|output|g' "${TAME}/$(dirname "${YAML_FILE}")/output.yaml"  # To rewrite the output

cd setup/0d
./fabm0d # Run the model

cd $TAME
cp "${TAME}/$(dirname "${YAML_FILE}")/output.yaml.backup" "${TAME}/$(dirname "${YAML_FILE}")/output.yaml" # Restore the original output file
echo '> Initial '"${TAME}/$(dirname "${YAML_FILE}")/output.yaml"' restored'

echo '>> DONE'

## Read the output files with an external script


