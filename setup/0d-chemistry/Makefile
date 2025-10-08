# SPDX-FileCopyrightText: 2025 Helmholtz-Zentrum hereon GmbH
# SPDX-FileContributor: Kai Wirtz <kai.wirtz@hereon.de>
# SPDX-FileContributor: Carsten Lemmen <carsten.lemmen@hereon.de>
# SPDX-License-Identifier: CC0-1.0


FABM_BASE := $(or $(FABM_BASE), $(FABMDIR), $(HOME)/prog/fabm)
ifneq ($(wildcard $(FABM_BASE)),$(FABM_BASE))
  $(error Directory FABM_BASE=$(FABM_BASE) does not exist)
else
  $(info Using FABM_BASE=$(FABM_BASE))
endif

FABM_TAME_BASE := $(or $(FABM_TAME_BASE), $(FABM_TAME_BASE), $(HOME)/prog/tame/tame)
ifneq ($(wildcard $(FABM_TAME_BASE)),$(FABM_TAME_BASE))
  $(error Directory FABM_TAME_BASE=$(FABM_TAME_BASE) does not exist)
else
  $(info Using FABM_TAME_BASE=$(FABM_TAME_BASE))	
endif

GOTM_BASE := $(or $(GOTM_BASE), $(GOTMDIR), $(HOME)/prog/gotm/v6.0)
ifneq ($(wildcard $(GOTM_BASE)),$(GOTM_BASE))
  $(error Directory GOTM_BASE=$(GOTM_BASE) does not exist)
else
  $(info Using GOTM_BASE=$(GOTM_BASE))
endif

DATE    =   $(shell date +%d%m%y)
APPNAME =   $(shell basename $(PWD))
CDIR    =   $(CURDIR)

.PHONY: all build exec clean backup run build-fortran build-python

all: exec

exec: build
	@echo $(CDIR)
	ln -sf $(FABM_BASE)/build-0d/fabm0d .
	./fabm0d

build: build-fortran # build-python

build-fortran:
	cmake -B $(FABM_BASE)/build-0d -S $(FABM_BASE)/src/drivers/0d \
	  -DFABM_HOST=0d -DFABM_INSTITUTES=tame -DFABM_TAME_BASE=$(FABM_TAME_BASE)/tame \
	  -DGOTM_BASE=$(GOTM_BASE) \
	  -DCMAKE_POLICY_VERSION_MINIMUM=3.5
	cmake --build $(FABM_BASE)/build-0d --parallel 4 --target fabm0d

build-python:
	python -m pip install --no-build-isolation --no-deps -v -C--build-option=build_ext \
	  -C--build-option="--cmake-opts=-DFABM_INSTITUTES='tame' -DFABM_FABM_TAME_BASE=${FABM_TAME_BASE}/tame" $(FABM_BASE)

clean:
	@echo -n "Making all clean ... "; \
	rm -f *.o *~ *% tmp* *.log "#"*"#"; \
	rm -i *.tgz; \
	echo "Done."

backup:	clean
	@echo -n "Creating Archiv of $(APPNAME) .. \n"; \
	tar zcvf $(APPNAME)$(DATE).tgz *.*  Makefile par/* mlab/*.m data/*.dat ; \
	echo "Done.";

run:
	@echo  "Executing fabm0d\n"; \
	./fabm0d
