# SPDX-FileContributor Carsten Lemmen <carsten.lemmen@hereon.de>
# SPDX-License-Identifier: GPL-2.0
# SPDX-FileCoprightText: 2024 Helmholtz-Zentrum hereon GmbH
# SPDX-FileCoprightText: 2014 Bolding and Bruggemann ApS

add_library(fabm_models_tame OBJECT
            tame_model_library.F90
            library/tame_library.F90
            library/tame_types.F90
           )

add_dependencies(fabm_models_tame fabm_base)
