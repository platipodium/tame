#! /bin/bash
# SPDX-FileCopyrightText: 2025 Helmholtz-Zentrum hereon GmbH
# SPDX-License-Identifier: CC0-1.0
# SPDX-FileContributor Ovidio Garcia-Oliva <ovidio.garcia@hereon.de>

export GOTM_BASE="${GOTM_BASE:-${HOME}/tools/gotm6}"
export FABM_BASE="${FABM_BASE:-${HOME}/tools/fabm/fabm}"
export FABM_TAME_BASE="${FABM_TAME_BASE:-${HOME}/tools/tame}"

mkdir -p ./build && cd ./build || exit
cmake -S "${GOTM_BASE}" -DGOTM_USE_FABM=ON \
  -DFABM_BASE="${FABM_BASE}" -DFABM_INSTITUTES="gotm;tame" \
  -DFABM_TAME_BASE="${FABM_TAME_BASE}/tame" \
  -DCMAKE_POLICY_VERSION_MINIMUM=3.5 || exit
make && make install  || exit

cd ../ || exit

ln -f ./build/gotm ./gotm
ln -sf fabm_cuxhaven.yaml fabm.yaml

if [ ! -f ./data/meteofile.csv ]; then
    cd ./data || exit
    bash get_data.sh
    cd ../ || exit
fi

./gotm
