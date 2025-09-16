# SPDX-License-Identifier: CC0-1.0
# SPDX-FileCopyrightText: 2024-2025 Helmholtz-Zentrum hereon GmbH
# SPDX-FileContributor Carsten Lemmen <carsten.lemmen@hereon.de>
# SPDX-FileContributor Kai Wirtz <kai.wirtz@hereon.de>

add_library(fabm_models_tame OBJECT
            bgc.F90
            phy.F90
            zoo.F90
            tame_model_library.F90
            #temperature.F90
            #library/tame_library.F90
            library/tame_stoich_functions.F90
            library/tame_functions.F90
            library/tame_types.F90
            library/chemistry_types.F90
            chemistry/chemistry.F90
            chemistry/detritus.F90
            chemistry/nutrient.F90
            #library/chemistry_types.F90
            time.F90
           )

add_dependencies(fabm_models_tame fabm_base)
