# SPDX-FileCopyrightText: 2025 Helmholtz-Zentrum hereon GmbH
# SPDX-License-Identifier: CC0-1.0
# SPDX-FileContributor Ovidio Garcia-Oliva <ovidio.garcia@hereon.de>

export GOTMDIR=$HOME/tools/gotm6
export FABMDIR=$HOME/tools/fabm/fabm
export TAMEDIR=$HOME/tools/generalized-aquatic-ecosystem-model/fortran

mkdir -p ./build
cd ./build
cmake -S $GOTMDIR -DGOTM_USE_FABM=ON -DFABM_BASE=$FABMDIR -DFABM_INSTITUTES="gotm;tame" -DFABM_TAME_BASE=$TAMEDIR    
make 
make install

cd ../

ln -f ./build/gotm ./gotm
ln -sf fabm_cuxhaven.yaml fabm.yaml

if [ ! -f ./data/meteofile.csv ]; then
    bash get_data.sh
fi


./gotm