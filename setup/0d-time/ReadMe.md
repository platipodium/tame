<!--
SPDX-FileCopyrightText: 2024-2025 Helmholtz-Zentrum hereon GmbH
SPDX-FileContributor: Carsten Lemmen <carsten.lemmen@hereon.de>
SPDX-License-Identifier: CC0-1.0
-->

# Demonstration 0D setup for the time module

The `tame_time` model is a demonstration model to showcase the capability of FABM models using diagnostics and self-dependencies to
keep track of

1. the time step since the last invocation of FABM.
2. the past value of a 3D variable and its change to the current time step.

Simply invoke `make` to compile and run the model and output to `stdout` tabular values of time and concentration.
