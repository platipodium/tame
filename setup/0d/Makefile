# SPDX-FileCopyrightText: 2025 Helmholtz-Zentrum hereon GmbH
#
# SPDX-License-Identifier: Apache-2.0

# Makefile created by kai wirtz, HZG  Nov 2013
#
ifndef FABMDIR
   FABMDIR = $(HOME)/prog/fabm
endif

DATE    =   `date +%d%m%y`
APPNAME =   `basename $(PWD)`
CDIR    =   $(CURDIR)

all: exec

exec:
	( cmake -B $(HOME)/build/fabm -S $(FABMDIR)/src/drivers/0d -DFABM_HOST=0d -DFABM_INSTITUTES=tame -DFABM_TAME_BASE=~/prog/tame/tame && make -C  $(HOME)/build/fabm install)
	@echo $(CDIR)
	( cp $(HOME)/local/fabm/0d/bin/fabm0d .)
#	( export FABMHOST=0d ; $(MAKE) -C $(FABMDIR)/src/ )
#	( export FABMHOST=0d ; $(MAKE) -C $(FABMDIR)/src/drivers/0d fabm0d )
#	( mv $(FABMDIR)/src/drivers/0d/fabm0d . )
#cd /home/build/fabm/0d && cmake $FABMDIR/src/drivers/0d -DCMAKE_INSTALL_PREFIX=$HOME/local/fabm/0d -DGOTM_BASE=$GOTMDIR -DCMAKE_Fortran_COMPILER=gfortran -DFABM_USE_NetCDF=OFF
# sed -i '.tmp' '1 s/.*/output2:/g' output.yaml
	./fabm0d

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
