<!--
SPDX-FileCopyrightText: 2025 Helmholtz-Zentrum hereon GmbH
SPDX-FileContributor: Carsten Lemmen <carsten.lemmen@hereon.de>
SPDX-License-Identifier: CC0-1.0
-->

# Dockers for TAME

TAME can be built in a Docker container, which is useful for testing and development. The Dockerfile is located hereis called [./Dockerfile-fabm](./Dockerfile-fabm). It is based on a minimal alpine image and installs all necessary dependencies to build TAME, FABM, and GOTM.

The `[Makefile](./Makefile)` in this directory can be used to build the Docker image and run it. The image is built with the command:

```bash
make fabm
```

Which internally creates a Docker image tagged `platipodium/tame:latest` as a multi-layer image with build-stages `gotm`, `gotm-fabm`, `fabm`, and `runtime`. The four layers are

- `gotm`: installs GOTM and its dependencies, creates a runnable `/opt/bin/gotm` executable
- `gotm-fabm`: installs FABM and its dependencies, creates a runnable `/opt/bin/gotm` executable with a 1D GOTM host and a `/opt/bin/fabm0d` executable.
- `fabm`: installs TAME and its dependencies, creates a runnable `/opt/bin/gotm` executable with a 1D GOTM host and a `/opt/bin/fabm0d` executable, but this time containing all TAME models.

You can run the image interactively with the command:

```bash
docker run -it platipodium/tame:latest
```

Or, you can run a GOTM TAME simulation with the command:

```bash
docker run -t platipodium/tame:latest /opt/bin/gotm
```
