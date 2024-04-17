/*************************************************************/
/*  Simple yaml Parser for creating FABM models in  F2003    */
/*                                                           */
/*   Kai W. Wirtz  (Hereon)                       10/03/24   */
/*************************************************************/
#include<stdio.h>
#include<stdlib.h>
#include<ctype.h>
#include<string.h>
#define MAXP 49 	/*   elements (pars, vars) of FABM model, should increase if you create highly complicated models */
#define NAML 49 	/*   model and directory environments     */
#define BLKN 3 	  /*   number of types in FABM yaml     */
#define NDEPV 3 	/*   number of dep variables types  */
#define SCALEFAC 0  /* 1: transforms all rate parameters already during get(.. */
                    /* 0: insert UNIT into _SET_ODE_(.. */
/* --------------------------------------------------------------------------- */
char *fil(char *,char *, int ),*insul(char *),*trim_space(char *);
int count_space(char *);

int main(int argn, char **argv )
{
/* ----------------------------------------------------------------------- */
/*  configuration variables for the model structure; you SHOULD edit here: */
/* ----------------------------------------------------------------------- */
char modname[NAML] ="selma";/* mmaecs odel name; this string will be also used for building the nml-filenames */
char yamlname[NAML]="fabm-selma";/* hzg_ model-tree name for the type definition; FABM-convention includes the    */
	/* "institutional" origin, corresponding to folders within fabm/src/model/... */
char elements[7]   = "N";	/* elements: "N", "CNP", or "CNPSF" with "S" for Si/silicon or "F" for Fe/iron */
char longelemnam[7][NAML]={"carbon","nitrogen","phosphorus","silicate"};
char dnam[17][NAML]={"par","temp","taub","chla","GPP","NPP"};
char dnam2[17][NAML]={"downwelling_photosynthetic_radiative_flux","temperature","bottom_stress","'mg chl a/m3', 'chlorophyll concentration'"," 'mmol/m3/d',   'gross primary production'"," 'mmol/m3/d',   'net primary production'"};
int dnami[17] =     {    0,     0,     1,     2,    2,    2, -1},dn0[NDEPV],dnn[NDEPV]; //index to depname
char dir_setup[NAML]	= ""; //"selma/";	/* directory where all the nml files  reside ./*/
char dirn_f90[2*NAML] 	= ""; //"/Users/wirtz/mossco/code/external/fabm/code/src/models/";	/* directory where all the input sources (model.F90,...) reside */

/* ---------------------------------------------------------------- */
/*            do not edit below ...                                 */
/* ---------------------------------------------------------------- */
int eoi,mi,mi0,pi,i,j,nvar[BLKN],nvart[BLKN],scalf[BLKN][MAXP],pv,nmli=0,out=1,tti,nls,indent,indent0,indenti,pc,nbi[2][MAXP],npb[2],pti[3][MAXP],ptn[3],ipb,IsBottom;
char line[256],tmp[256],*lr,c,*cp,*cp1,coli[21]="--------------------";
FILE *sp,*so,*spt;
char ptkey[4][17]={"real(rk)","integer","logical","character(len=64)"}, keys[5][4]={"RHS","ODE","GET"};
char bsn[2][NAML]={"","bottom_"},depname[NDEPV][NAML]={"dependency","horizontal_dependency","diagnostic_variable"},vardepname[3][NAML]={"register_state_variable","register_state_dependency","get_parameter"};
char scalname[3][NAML]={"",", scale_factor=1.0_rk/secs_per_day",""},wsnam0[NAML]=", vertical_movement=",wsname[NAML]="ws";
char lstname0[BLKN][NAML] = {"deps","diags","aux"};
/* names of the files for 1) external forcing and 2) diagnostics and, optional (AUX=1), 3) derived parameters */
		/* order is fixed; ""deps" creates the full filename "modname_deps.lst"  */
// char traitpre[5]="phy%";/* structure name for functional group variables */
char insname[NAML],modulname[NAML],modulname0[NAML]="qwetz",yamlfname[NAML],outn[3*NAML];
char blockname[BLKN][NAML] = {"initialization","coupling","parameters"}; // 3 groups: parameters, init=states, coupling between modules
char yesno[2][NAML] = {"false","true"};
char dirn[NAML],modn[NAML], varname[BLKN][MAXP][NAML],varval[BLKN][MAXP][NAML],unit[BLKN][MAXP][NAML],pcom[BLKN][MAXP][3*NAML];
int nmltype,partype[MAXP],nump,ChemSpec[MAXP], ni0,nis,nelements;
/* printf("argn=%d\t",argn);printf("argv=%s\t",argv[1]);
if (argn>=2)   strcpy(simfile,argv[1]); else  exit(0);
if (argn>=3)   strcpy(cas,argv[2]);printf("AS=%s\n",cas);*/
nelements = strlen(elements);
tti=-1;
for (j=0;j<NDEPV;j++) dn0[j]=-1,dnn[j]=0;
i=j=0;
while(dnami[i]>=0)
  {
	while(dnami[i]>j) j++;
	if(dn0[j]<0) dn0[j]=i,dnn[j]=1;
	else dnn[j]++;
	i++;
	}
for (j=0;j<NDEPV;j++) printf("%d %d\t",dn0[j],dnn[j]);

// ----------------------------------------------------------
//  read all info from single yaml file
// ----------------------------------------------------------
sprintf(yamlfname,"%s%s.yaml",dir_setup,yamlname);
if(out) printf("extracting from %s \n",yamlfname);
// ----------------------------------------------------------
//  open yaml file for reading
// ----------------------------------------------------------
sp=fopen(yamlfname,"r");
if(sp==NULL) {printf("Problem while opening %s !\n",yamlfname),exit(0);}

// ----------------------------------------------------------
// trailing comments, if any
lr=(char *)1; i=0;
// ----  reads  lines until "instances:" is found
while(lr!=NULL && i==0)
	  {
    lr=fgets(line,256,sp);
    if(strstr(line,"instances:")>=0) i=1;
		}
// instances name
lr=fgets(line,256,sp);
// ----------------------------------------------------------
// main parser loop
//     retrieve info line per line
// ----------------------------------------------------------
pi=0; mi=0;mi0=0;
while(lr!=NULL && eoi==0)
  {
	for (j=0;j<BLKN;j++) nvar[j]=nvart[j]=0;

	indenti=count_space(line);
	line[strlen(line)-1]='\0'; // cut final :
  strcpy(insname,trim_space(line));
	// module name
  lr=fgets(line,256,sp);
	sscanf(line," %s %s",tmp,modulname);

	if(out) printf("instance \t%s %s\n",insname,modulname);
// ----  reads  lines until blockname is found
	indent0=count_space(line);
// ----------------------------------------------------------
// new block
// ----------------------------------------------------------
//  coupling
  lr=fgets(line,256,sp);indent=indenti+1;
	while(lr!=NULL && indent>indenti)
    {
		//printf("%d > %d %d ? \nline=%s\n",indent,indenti,indent0,line);
		//  getchar();

		//for (j=0;j<BLKN;j++) printf("%s %s\n",blockname[j],strstr(line,blockname[j]));;
		for (j=0;j<BLKN;j++) if(strstr(line,blockname[j]) != NULL) break;

		if(j==BLKN)
   		if(strchr(line,'#')==NULL) {
	   		printf("line %s contains neither blockname nor comment\n",line);
		  	exit(0);
		    }
// ----------------------------------------------------------
// ----  retrieve parameter/variable/coupling info line per line
    i=0;
// ----  check for indent
    lr=fgets(line,256,sp);
    indent=count_space(line);
		while(lr!=NULL && indent>indent0)
			{
		// parameters name
			sscanf(line," %s %s",tmp,varval[j][i]);
			tmp[strlen(tmp)-1]='\0'; // cut final :
			strcpy(varname[j][i],tmp);
			cp=strchr(line,'#');
			strcpy(tmp,(++cp));
			strcpy(line,trim_space(tmp));
     // checks for rates
		  strcat(tmp,line);
			if(strstr(tmp,"vertical velocity")!=NULL && j==2)
			   strcpy(wsname,varname[j][i]);

		  if(strstr(tmp," rate")==NULL && strstr(tmp,"velocity")==NULL )
			  if(strstr(tmp,"concentration")==NULL || j==2 )
			    scalf[j][i]=0;
			  else
			    scalf[j][i]=2;
			else
			  scalf[j][i]=1;

				// retrieve units from comment in parenthesis
			cp=strchr(line,'(');
			cp1=strchr(line,')');
			*unit[j][i]='\0';
			if(cp!=NULL && cp1>cp+1)
				 {
				 *pcom[j][i]='\0';
				 strncat(pcom[j][i],line,cp-line-1);
				 strcat(pcom[j][i],cp1+1);
				 strncat(unit[j][i],cp+1,cp1-cp-1);
				 }
			else
				 strcpy(pcom[j][i],line);
      //cp=strstr(pcom[j][i],',default');
			if(out) printf("%d %d\t%s %s\t%s\n",j,i,varname[j][i],varval[j][i],pcom[j][i]);
			lr=fgets(line,256,sp);
			indent=count_space(line);
			i++;
      if(i>=MAXP) {printf("Too many entries found %d>=%d !\n",i,MAXP),exit(0);}
			}
		nvar[j]=i;
		//printf("indent=%d  %s %d %d\n",indent,line,lr,j);
	  } // end module

  // ------------------------------------------------------------
  // ------------------------------------------------------------

//printf("----------------\n%s %s\n\n",insname,modulname);
	// output of real and not yet used models
	//printf("modulname\t%s %s\t%d %d\n",modulname,modulname0,strstr(modulname,"constant_"),strstr(modulname,modulname0));

  if( strstr(modulname,"constant_") == NULL && strstr(modulname,modulname0) == NULL)
	  {
		// ------------------------------------------------
		//   organize FABM model dir and name
    strcpy(modulname0,modulname);
    strcpy(tmp,modulname);
		cp=strchr(tmp,'/');
		*dirn = '\0';
    if(cp!=NULL) strncat(dirn,tmp,cp-tmp),strcpy(modn,cp+1);
		else strcpy(modn,tmp);
    strcpy(varval[1][MAXP-1],"");
		// ------------------------------------------------
		//   open file for FABM model code
		// ------------------------------------------------
		sprintf(outn,"%s%s_gen.F90",dirn_f90,modulname);
	  so=fopen(outn,"w");
	  printf("Opening %s !\n",outn);
	  if(so==NULL) {printf("Problem while opening %s !\n",outn),exit(0);}
		if(out) printf("writing FABM model code into %s...\n",outn);
		// ------------------------------------------------
		//   write header info
    fprintf(so,"#include \"fabm_driver.h\"\n!%s%s\n!\n!\t%s\n!\n! here some comments\n!%s%s\n!\n!%s%s\n",coli,coli,modulname,coli,coli,coli,coli);
		fprintf(so,"! !INTERFACE:\n MODULE %s_%s\n!\n",dirn,modn);
		fprintf(so,"! !USES:\n use fabm_types\n\n"); //use fabm_driver
		fprintf(so," implicit none\n\n private\n!\n ! !PUBLIC_DERIVED_TYPES:\n");
		fprintf(so," type,extends(type_base_model),public :: type_%s_%s\n",dirn,modn);
		// types for state variables
		//printf("\n\n nvar = %d %d\n",nvar[0],nvar[1]);
    strcpy(scalname[2],wsnam0); strcat(scalname[2],wsname);

		for (j=0,IsBottom=0;j<2;j++)
			if (nvar[j]>0)
			  {
				// distinguish pelagic and bottom states
				npb[0]=npb[1]=0;
				for(i=0;i<nvar[j];i++)
				  {
					ipb=0;
				  if(strstr(unit[j][i],"m2") != NULL ||  strstr(pcom[j][i],"benthic") != NULL) ipb=1;
					nbi[ipb][npb[ipb]++]=i;
          }
				// write declaration with list for each case (pelagic/bottom)
				printf("npb=%d %d\n",npb[0],npb[1]);

				for (ipb=0;ipb<2;ipb++)
	        if(npb[ipb]>0)
					  {
						if(ipb==1) IsBottom=1;
		        fprintf(so,"\ttype (type_%sstate_variable_id) :: ",bsn[ipb]);
				  	for(i=0;i<npb[ipb]-1;i++)
					  	fprintf(so,"id_%s,",varname[j][nbi[ipb][i]]);
					  fprintf(so,"id_%s\n",varname[j][nbi[ipb][i]]);
					  }
	      } // for (j=0;j<2  state variables
		// ------------------------------------------------
		//   declaration of dependency and diag ids
		// ------------------------------------------------
		for (j=0;j<NDEPV;j++)
		  if(dnn[j]>0)
			  {
			 	fprintf(so,"\ttype (type_%s_id) :: ",depname[j]);
	      for(i=0;i<dnn[j]-1;i++)
				  fprintf(so,"id_%s,",dnam[dn0[j]+i]);
				fprintf(so,"id_%s\n",dnam[dn0[j]+i]);
	      }
		// ------------------------------------------------
		//   declaration of parameters
		// ------------------------------------------------
		printf("\nnpar = %d\n",nvar[2]);
		j=2;
	  if(nvar[j]>0)
			{
			for (i=0;i<3;i++) ptn[i]=0;

			for(i=0;i<nvar[j];i++)
			  {
				ipb=0;
			  if(strstr(varval[j][i],".") == NULL ||  strstr(pcom[j][i],"switch") != NULL) ipb=1;  // integer
				if(strstr(varval[j][i],"true") != NULL || strstr(varval[j][i],"false") != NULL) ipb=2;  // bool
				//printf("%d %s %s\t%d %d\n",i,varname[j][i],varval[j][i],ipb,ptn[ipb]);
				pti[ipb][ptn[ipb]++]=i;
	      }
		// write declaration with list for each parameter type
			for (ipb=0;ipb<3;ipb++)
				if(ptn[ipb]>0)
					{
					fprintf(so,"\t%s :: ",ptkey[ipb]);
			  	for(i=0;i<ptn[ipb]-1;i++)
				  	fprintf(so,"id_%s,",varname[j][pti[ipb][i]]);
				  fprintf(so,"id_%s\n",varname[j][pti[ipb][i]]);
					}
		  // wrap-up
		  fprintf(so,"\n contains\n\tprocedure :: initialize\n\tprocedure :: do\n");
		  if(IsBottom) fprintf(so,"\tprocedure :: do_bottom\n");
		  fprintf(so," end type\n!EOP\n!%s%s\n",coli,coli);
		  fprintf(so," CONTAINS\n\n!%s%s!BOP\n!\n",coli,coli);
			fprintf(so,"! !IROUTINE: Initialise the %s model\n!\n",modulname);
		  fprintf(so,"! !INTERFACE:\n subroutine initialize(self,configunit)\n!\n");
		  fprintf(so,"! !DESCRIPTION:\n! Reading namelist and registration of variables with FABM\n!\n! !INPUT PARAMETERS:");
			fprintf(so," class(type_%s_%s),intent(inout),target :: self\n integer,\t\tintent(in)\t\t:: configunit\n\n",dirn,modn);
		  fprintf(so," real(rk),parameter :: secs_per_day = 86400._rk\n\n");
		// ------------------------------------------------
    // declaration of all elements (vars, pars, ..)
			for (j=0;j<3;j++)
			  {
				if (j==1) strcpy(tmp,"");
				else strcpy(tmp,", default=");
				if(nvar[j]>0)
					for(i=0;i<nvar[j];i++)
						fprintf(so," call self%c%s(self%cid_%s, '%s','%s','%s'%s%s %s)\n",'%',vardepname[j],'%',varname[j][i],varname[j][i],unit[j][i],pcom[j][i],tmp,varval[j][i*(j!=1)+(MAXP-1)*(j==1)],scalname[scalf[j][i]]);
				}

	// ------------------------------------------------
	//   declaration of dependency and diag ids
	// ------------------------------------------------
			for (j=0;j<NDEPV;j++)
				for(i=0;i<dnn[j];i++)
					{
					if(j<2)
					  sprintf(tmp,"standard_variables%c%s",'%',dnam2[dn0[j]+i]);
					else
						sprintf(tmp,"'%s', %s",dnam[dn0[j]+i],dnam2[dn0[j]+i]);
					fprintf(so," call self%cregister_%s(self%cid_%s, %s)\n",'%',depname[j-(j==1)],'%',dnam[dn0[j]+i],tmp);
					}
			fprintf(so,"\n end subroutine initialize\n!EOC\n\n!%s%s!BOP\n!\n\n",coli,coli);
			fprintf(so,"! !IROUTINE: Right hand sides of %s model\n!\n",modulname);
			fprintf(so,"! !INTERFACE:\n subroutine do(self,_ARGUMENTS_DO_)\n!\n");
			fprintf(so,"! !LOCAL VARIABLES:\n	real(rk) :: ");
// write declaration with list for each parameter type
      for(i=0,j=0;i<dnn[j]-1;i++)
				fprintf(so,"%s, ",dnam[dn0[j]+i]);
      if(dnn[j]>0) fprintf(so,"%s\n",dnam[dn0[j]+dnn[j]-1]);
			for (j=0;j<2;j++)
				if(nvar[j]>0)
					{
					fprintf(so,"  real(rk) :: ");
					for(i=0;i<nvar[j]-1;i++)
						fprintf(so,"%s,",varname[j][i]);
					fprintf(so,"%s\n",varname[j][i]);
					}
		  fprintf(so,"! real(rk) :: \n  real(rk),parameter :: secs_per_day = 86400._rk\n\n");
	    fprintf(so,"! Enter spatial_loops (if any)\n _LOOP_BEGIN_\n\n");

			for (j=0;j<2;j++)
				if(nvar[j]>0)
					for(i=0;i<nvar[j]-1;i++)
						 fprintf(so,"  _GET_(self%cid_%s, %s)\t\t! %s\n",'%',varname[j][i],varname[j][i],pcom[j][i]);
			j=0;
		  for(i=0;i<dnn[j];i++)
			  fprintf(so,"  _GET_(self%cid_%s, %s)\t\t! %s\n",'%',dnam[dn0[j]+i],dnam[dn0[j]+i],dnam2[dn0[j]+i]);
			fprintf(so,"\n!%s%s\n\n\n!%s%s\n",coli,coli,coli,coli);
			for (j=0;j<2;j++)
			  for(i=0;i<dnn[j];i++)
          fprintf(so,"! _SET_ODE_(self%cid_%s,  )\n",'%',varname[j][i]);
		  j=2;
			for(i=0;i<dnn[j];i++)
				fprintf(so,"!  _SET_DIAGNOSTIC_(self%cid_%s, )\t\t! %s\n",'%',dnam[dn0[j]+i],dnam2[dn0[j]+i]);

	    fprintf(so,"\n! Leave spatial loops (if any)\n _LOOP_END_\n\n END subroutine do\n!EOC\n\n!%s%s",coli,coli);
			printf("IsBottom %d npb=%d\n",IsBottom,npb[1]);
			if(IsBottom)
			  {
				fprintf(so,"\n!BOP\n!\n\n! !IROUTINE: Right hand sides of benthic %s model\n!\n",modulname);
				fprintf(so,"! !INTERFACE:\n subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)\n!\n");
				fprintf(so," class(type_%s_%s),intent(in) :: self\n _DECLARE_ARGUMENTS_DO_BOTTOM_\n!\n",dirn,modn);
				fprintf(so,"! !LOCAL VARIABLES:\n real(rk) :: ");
	// write declaration with list for each parameter type
	      ipb=1;
  		  for(i=0;i<npb[ipb]-1;i++)
				  fprintf(so,"%s,",varname[j][nbi[ipb][i]]);
			  fprintf(so,"%s\n",varname[j][nbi[ipb][i]]);
				j=1;
				if(dnn[j]>0)
				  {
				  fprintf(so," real(rk) :: ");
		      for(i=0;i<dnn[j]-1;i++)
					  fprintf(so,"%s,",dnam[dn0[j]+i]);
					fprintf(so,"%s\n",dnam[dn0[j]+i]);
          }
			  fprintf(so,"\n!EOP\n!%s%s\n!BOC\n",coli,coli);
			  fprintf(so," CONTAINS\n\n!%s%s!BOP\n!\n! Enter spatial loops over the horizontal domain (if any).\n _HORIZONTAL_LOOP_BEGIN_\n",coli,coli);

				for(i=0;i<dnn[j];i++)
			    fprintf(so," _GET_HORIZONTAL_(self%cid_%s, %s)\t\t! %s\n",'%',dnam[dn0[j]+i],dnam[dn0[j]+i],dnam2[dn0[j]+i]);

			  for(i=0;i<npb[ipb];i++)
				  fprintf(so,"! _SET_BOTTOM_ODE_(self%cid_%s,  )\n",'%',varname[j][nbi[ipb][i]]);
        fprintf(so,"\n! Leave spatial loops over the horizontal domain (if any).\n_HORIZONTAL_LOOP_END_\n\nend subroutine do_bottom\n!EOC\n!%s%s\n",coli,coli);
        }
		  fprintf(so,"\n!%s%s\nEND  MODULE %s_%s\n\n!%s%s\n",coli,coli,dirn,modn,coli,coli);

		  fclose(so);
			for (j=0;j<BLKN;j++)
			  if(nvar[j]>0)
			    {
			    printf("---------\n%s: (%d)\n",blockname[j],nvar[j]);
			    for(i=0;i<nvar[j];i++)
			      printf("%s = %s\t%s\n",varname[j][i],varval[j][i],pcom[j][i]);
			    }
			} // if(nvar[j]>0) parameter output
		}	// if model out
	for (j=0;j<BLKN;j++) nvart[j]+=nvar[j];
  } // end all
for (j=0;j<BLKN;j++)
  printf("%s: %d\n",blockname[j],nvart[j]);
mi0=mi;
}

char * trim_space(char *str) {
    char *end;
    /* skip leading whitespace */
    while (isspace(*str)) {
        str = str + 1;
    }
    /* remove trailing whitespace */
    end = str + strlen(str) - 1;
  /*  while (end > str && isspace(*end)) {
        end = end - 1;
    }*/
    /* write null character */
    *(end) = '\0';
    return str;
}

int count_space(char *str) {
    int n=0;
    /* skip leading whitespace */
    while (isspace(*str)) str++,n++;
    return n;
}
