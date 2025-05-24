#! /bin/bash
# SPDX-FileCopyrightText: 2025 Helmholtz-Zentrum hereon GmbH
# SPDX-License-Identifier: CC0-1.0
# SPDX-FileContributor Ovidio Garcia-Oliva <ovidio.garcia@hereon.de>

NAMES=( "Cuxhaven_DWD!Lufttemperatur"
        "Cuxhaven_DWD!Windgeschwindigkeit"
        "Cuxhaven_DWD!Windrichtung"
)

for NAME in "${NAMES[@]}"; do
    FILENAME="$NAME.zip"
    test -f "${FILENAME}" && continue
    URL="https://www.kuestendaten.de/DE/dynamisch/appl/data/daten_prod/prodNeuProd/direct_download/$FILENAME"
    wget "${URL}" -O "${FILENAME}"
    unzip "${FILENAME}"
    TXTFILE="${NAME}.txt"
    sed -i'' '/^[^0-9]/ s/^/# /' "${TXTFILE}" #Comment lines that are not numerical
done

Rscript setup_data.R && (rm -f ./*.zip ./*.txt)
