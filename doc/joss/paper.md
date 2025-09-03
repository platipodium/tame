---
title: "OxyPOM and DiaMO: simple FABM models for dissolved oxygen and biogeochemistry"
tags:
    - FABM
    - Biogeochemistry
    - Ecology
    - Food web
    - Fortran
    - Water quality
    - Aquatic
    - Trait-based
authors:
    - name: Carsten Lemmen
      orcid: 0000-0003-3483-6036
      affiliation: 1
    - name: Ovidio García-Oliva
      orcid: 0000-0001-6060-2001
      affiliation: 1
    - name: Thomas Imbert
      orcid:
      affiliation: 1
    - name: Kai W. Wirtz
      orcid:
      affiliation: 1
affiliations:
    - name: Helmholtz-Zentrum Hereon, Institute of Coastal Systems - Modeling and Analysis, Germany, ovidio.garcia@hereon.de
      index: 1
      ror: 03qjp1d79
date: 3 Sep 2025
year: 2025
bibliography: paper.bib
SPDX-FileCopyrightText: 2025 Helmholtz-Zentrum hereon GmbH
SPDX-License-Identifier: CC-BY-4.0
---

# Summary

The Trait-based Adaptive Ecosystem Model (TAME) is a flexible and modular model system representing the material flows and trophic interactions in an aquatic environment.  It is space-agnostic and can therefore be applied in zero-dimensional to three-dimensional hydrodynamic contexts via the Framework for Aquatic Biogeochemical Models (FABM @Bruggeman2014).

As components of TAME, biogeochemical and ecological processes mediated by bacteria, by phytoplankton, and by zooplankton and higher trophic levels are represented.  Material flows of main carbon and a user-configurable constellation of major nutrients like phosphorous and nitrogen, or minor nutrients like silicate or iron are represented, along with their variable stoichiometry in different biological compartments.
Stoichiometry and body size of biological compartments are adaptive traits that evolve towards optimizing the fitness of the respective compartments.

# Statement of need

The rich spatial and temporal complexity of aquatic food webs and material flows is currently not represented by state-of-the art biogeochemical models, that often lack the food web interaction, nor by ecological (food-web) models, that often are not suitable to be included in dynamic simulation models.  On the other hand, trait-based models have so far been mostly theoretic in nature and rarely tested in the field.

We therefore take on a perspective of rigorously applying adaptive trait dynamics with the goal to
* represent the complex time and space-dynamic behaviour of ecosystems
* provide minimally complex models that can be efficiently used even in Earth System Models.

A predecessor to TAME was the Model for Aquatic EcoSystems (MAECS, @Wirtz2016).  TAME builds on the stoichiometric and physiological trait formulations from MAECS, but adds body size as an important trait, flexible coupling between its submodels and choice of chemical species resolution.


# Model documentation and license

The models are documented in short form in the `ReadMe.md` section of the repository.
Open access data from third parties are not included with the model, and scripts for their download are included.
Our own models, scripts and documentations are are released under open source licenses, foremost Apache 2.0, CC0-1.0, and CC-BY-SA-4.0; a comprehensive documentation of all licenses is provided via REUSE Software.

# Acknowledgements

The development of this model was made possible by the grant no. 03F0954D of the BMBF as part of the DAM mission ‘mareXtreme’, project ElbeXtreme. We are grateful for the open source community that facilitated this research, amongst them the developers of and contributors to FABM, GOTM, R, pandoc, and LaTeX.

# References
