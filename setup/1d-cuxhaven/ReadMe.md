<!---
SPDX-FileCopyrightText: 2025 Helmholtz-Zentrum hereon GmbH
SPDX-License-Identifier: CC0-1.0
SPDX-FileContributor Ovidio Garcia-Oliva <ovidio.garcia@hereon.de>
-->

# Testcase 1d-cuxhaven

This directory includes a testcase for validation of TAME in a station in the Elbe estuary (Cuxhaven).
It includes:

* configuration files
* external data acquisition and management

Just run `bash test_gotm.sh`, if successful, `output.nc` should be produced.

Typically, no editing is needed to any file as long as the enviroment variables `GOTMDIR`, `FABMDIR`, and `TAMEDIR` are correctly defined in .

## Configuration files

Three configuration files are required to set-up the simulation:

* `gotm.yaml`: physical configuration (coordinates, tidal components, depth, etc.) and location of forcing files.

* `fabm_cuxhaven.yaml`: configuration (parameterization and initial conditions) of TAME models.

* `output.yaml`: configuration for model output.

## External data acquisition and management

Realistic forcing data (`.\data\meteofile.csv`) inspired in the Cuxhaven sampling station in the Elbe estuary is included.
The complete dataset can be downloaded runing the script `./data/get_data.sh`.
`get_data.sh` prepares the forcing data from [kuestendaten.de](https://www.kuestendaten.de) and executes
`setup_data.R`, which formats the data to be read by `GOTM`.
The downloaded files are used under license [DL-DE->Zero-2.0](https://www.govdata.de/dl-de/zero-2-0).

## Licenses

Apache-2.0
CC0-1.0

## To do

Implementing visualization scripts.
