# -*- coding: utf-8 -*-
"""
This script creates from a yaml file the structure for a
new FABM model using a Jinja2 template

SPDX-FileCopyrightText: 2024 Helmholtz-Zentrum hereon GmbH
SPDX-License-Identifier: Apache-2.0
SPDX-FileContributor: Carsten Lemmen <carsten.lemmen@hereon.de>
"""

import argparse
from pathlib import Path

import yaml
from jinja2 import Environment, FileSystemLoader


def create_fabm(yaml_path, output_path):
    # Load YAML configuration
    with open(yaml_path, "r") as yaml_file:
        config = yaml.safe_load(yaml_file)

    # Setup Jinja2 environment
    env = Environment(
        loader=FileSystemLoader("."), trim_blocks=True, lstrip_blocks=True
    )

    # Load template
    template = env.get_template(
        "model.F90.in"
    )  # Assuming your template file is named 'template.f90'

    # Render template with configuration data
    rendered_code = template.render(config)

    # Write rendered code to output file
    with open(output_path, "w") as output_file:
        output_file.write(rendered_code)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Generate Fortran code from a YAML configuration file and a template."
    )
    parser.add_argument(
        "yaml_path",
        type=Path,
        help="Path to the YAML configuration file",
        default="entities.yaml",
        nargs="?",
    )
    parser.add_argument(
        "output_path",
        type=Path,
        help="Path to the output Fortran file",
        default="model.F90",
        nargs="?",
    )

    args = parser.parse_args()

    create_fabm(args.yaml_path, args.output_path)
