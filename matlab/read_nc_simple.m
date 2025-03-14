out =0 ;
if ~exist(datf,'file') warning('File %s does not exist, skipped',datf); 
else
  ncid=netcdf.open(datf,'NC_NOWRITE');
  [ndim ndvar natt udimid] = netcdf.inq(ncid);
  iv=1; vars={};
  for id=0:ndvar-1
    [datvarname,xtype,dimids,natts]=netcdf.inqVar(ncid,id);
    unit = netcdf.getAtt(ncid,id,'units');
    if strcmp(datvarname,'time')
      datime=netcdf.getVar(ncid,id); %"seconds since 2002-01-01" ;
      tvec=datevec(datime);
      dayear= tvec(:,1); dayears=unique(dayear);
      datimeunits=unit;
      if out, fprintf(' %d data times found:\t%s\n',length(datime),datimeunits); end
    else
      varid = netcdf.inqVarID(ncid,datvarname);
      vars{iv} = datvarname;
      units{iv} = unit;
      tmp = squeeze(netcdf.getVar(ncid,id));
      nt=length(find(~isnan(tmp)));
      if out, fprintf('%d %s found %d values\n',iv,(vars{iv}),nt); end
      if nt>= length(datime)
        data(is,iv,:)=tmp;
        iv=iv+1;
      end
    end
  end
end