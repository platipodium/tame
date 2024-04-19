# Makefile created by kai wirtz, HZG  Nov 2013
#
ifndef FABMDIR
   FABMDIR = $(HOME)/mossco/code/external/fabm/code
endif

DATE    =   `date +%d%m%y`
APPNAME =   `basename $(PWD)`
CDIR    =   $(CURDIR)

all: exec

exec:
	(  export FABMHOST=0d ; cd $(HOME)/build/fabm0 && $(MAKE) install )
	@echo $(CDIR)
	( cd $(CDIR) && mv $(HOME)/local/fabm0/0d/bin/fabm0d .)
#	( export FABMHOST=0d ; $(MAKE) -C $(FABMDIR)/src/ )
#	( export FABMHOST=0d ; $(MAKE) -C $(FABMDIR)/src/drivers/0d fabm0d )
#	( mv $(FABMDIR)/src/drivers/0d/fabm0d . )
#cd /home/build/fabm/0d && cmake $FABMDIR/src/drivers/0d -DCMAKE_INSTALL_PREFIX=$HOME/local/fabm/0d -DGOTM_BASE=$GOTMDIR -DCMAKE_Fortran_COMPILER=gfortran -DFABM_USE_NetCDF=OFF
#
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
