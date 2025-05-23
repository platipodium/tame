add_library(fabm_models_tame OBJECT
           tame_model_library.F90
           npzdep.F90
           )

#target_link_libraries(fabm_models_tame PRIVATE fabm_base)
add_dependencies(fabm_models_tame fabm_base)

# Copyright (C) 2014 - Karsten Bolding and Jorn Bruggeman (BB)         !
# SPDX-FileCopyrightText: 2025 Helmholtz-Zentrum hereon GmbH
#
# SPDX-License-Identifier: Apache-2.0
