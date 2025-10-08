# SPDX-FileCopyrightText: 2025 Helmholtz-Zentrum hereon GmbH
# SPDX-FileContributor: Kai Wirtz <kai.wirtz@hereon.de>
# SPDX-FileContributor: Carsten Lemmen <carsten.lemmen@hereon.de>
# SPDX-License-Identifier: CC0-1.0

# Created originally by Kai Wirtz, HZG, Nov 2013

FABM_BASE := $(or $(FABM_BASE), $(FABMDIR), $(HOME)/prog/fabm)
ifneq ($(wildcard $(FABM_BASE)),$(FABM_BASE))
  $(error Directory FABM_BASE=$(FABM_BASE) does not exist)
else
  $(info Using FABM_BASE=$(FABM_BASE))
endif

FABM_TAME_BASE := $(or $(FABM_TAME_BASE), $(HOME)/prog/tame/tame)
ifneq ($(wildcard $(FABM_TAME_BASE)),$(FABM_TAME_BASE))
  $(error Directory FABM_TAME_BASE=$(FABM_TAME_BASE) does not exist)
else
  $(info Using FABM_TAME_BASE=$(FABM_TAME_BASE))
endif

GOTM_BASE := $(or $(GOTM_BASE), $(GOTMDIR), $(HOME)/prog/gotm)
ifneq ($(wildcard $(GOTM_BASE)),$(GOTM_BASE))
  $(error Directory GOTM_BASE=$(GOTM_BASE) does not exist)
else
  $(info Using GOTM_BASE=$(GOTM_BASE))
endif

DATE    =   $(shell date +%d%m%y)
APPNAME =   $(shell basename $(PWD))
CDIR    =   $(CURDIR)
FABM_BUILD_DIR :=  $(or $(FABM_BUILD_DIR), $(HOME)/build/fabm)

.PHONY: all build exec clean backup run

all: exec

exec: build
	@echo $(CDIR)
	ln -sf $(FABM_BUILD_DIR)/fabm0d .
	./fabm0d

build:
	cmake -B $(FABM_BUILD_DIR) -S $(FABM_BASE)/src/drivers/0d \
	  -DFABM_HOST=0d -DFABM_INSTITUTES=tame -DFABM_TAME_BASE=$(FABM_TAME_BASE) \
	  -DGOTM_BASE=$(GOTM_BASE) \
	  -DCMAKE_POLICY_VERSION_MINIMUM=3.5
	cmake --build $(FABM_BUILD_DIR) --parallel 4 --target fabm0d

clean:
	@echo -n "Making all clean ... " ; \
	rm -f *.o *~ *% tmp* *.log "#"*"#" ; \
	if ls *.tgz 1> /dev/null 2>&1; then rm -i *.tgz; fi ; \
	echo "Done."

backup:	clean
	@echo -n "Creating Archiv of $(APPNAME) .. \n"; \
	tar zcvf $(APPNAME)$(DATE).tgz *.*  Makefile par/* mlab/*.m data/*.dat ; \
	echo "Done.";

run:
	@echo  "Executing fabm0d\n"; \
	./fabm0d
