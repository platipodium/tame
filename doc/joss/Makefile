# SPDX-FileCopyrightText: 2023-2025 Helmholtz-Zentrum hereon GmbH
# SPDX-License-Identifier: CC0-1.0
# SPDX-FileContributor: Carsten Lemmen <carsten.lemmen@hereon.de>

MEXT = md
SRC=paper.$(MEXT)

PDF=$(SRC:.md=.pdf)
XREF=$(SRC:.md=.crossref)
HTML=$(SRC:.md=.html)

.PHONY: all cff test pdf local-pdf clean
.FORCE: pdf

all:
cff pdf: $(PDF)

%.pdf:	%.md %.bib Makefile
	docker run -it -v $(PWD):/data \
    -u $(id -u):$(id -g) \
    openjournals/inara \
    -o pdf,cff $(SRC)

test:
	docker run -it -v $(PWD):/data \
    -u $(id -u):$(id -g) \
    openjournals/inara \
    -o pdf,cff $(SRC) /bin/bash

local-pdf:
	pandoc --output=$(PDF) --template=joss-template.tex --bibliography=paper.bib  $(SRC)

clean:
	@rm -f $(PDF) $(XREF) $(HTML) CITATION.cff paper.preprint.* *aux *bcf *.xml *.out *.log $(TEX) *.html
	@rm -rf jats
