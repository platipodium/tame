# SPDX-FileCopyrightText: 2025 Helmholtz-Zentrum hereon GmbH
# SPDX-FileContributor: Kai Wirtz <kai.wirtz@hereon.de>
# SPDX-FileContributor: Carsten Lemmen <carsten.lemmen@hereon.de>
# SPDX-License-Identifier: CC0-1.0

# Created originally by Kai Wirtz, HZG, Nov 2013

FABMDIR := $(or $(FABMDIR), $(FABM_BASE), $(HOME)/prog/fabm)
ifneq ($(wildcard $(FABMDIR)),$(FABMDIR))
  $(error Directory FABMDIR=$(FABMDIR) does not exist)
endif

FABM_TAME_BASE := $(or $(FABM_TAME_BASE), $(HOME)/prog/tame/tame)
ifneq ($(wildcard $(FABM_TAME_BASE)),$(FABM_TAME_BASE))
  $(error Directory FABM_TAME_BASE=$(FABM_TAME_BASE) does not exist)
endif

DATE    =   $(shell date +%d%m%y)
APPNAME =   $(shell basename $(PWD))
CDIR    =   $(CURDIR)

.PHONY: all build exec clean backup run

all: exec

exec: build
	@echo $(CDIR)
	ln -sf $(HOME)/build/fabm/fabm0d .
	./fabm0d

build:
	cmake -B $(HOME)/build/fabm -S $(FABMDIR)/src/drivers/0d \
	  -DFABM_HOST=0d -DFABM_INSTITUTES=tame -DFABM_TAME_BASE=$(FABM_TAME_BASE) \
	  -DGOTM_BASE=$(GOTM_BASE) \
	  -DCMAKE_POLICY_VERSION_MINIMUM=3.5
	cmake --build $(HOME)/build/fabm --parallel 4 --target fabm0d

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
