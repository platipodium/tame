% 
% matlab script for visualizing and comparing 0D model results (from netcdf files)
%
% kai wirtz (hereon 2024-2025)
%
% name of variables to plot
varn={'bgc_NO3';'bgc_PO4';'bgc_rate';'bgc_det_C';'bgc_dom_C';'bgc_dom_N';'bgc_dom_P';'temp';'phytoplankton_phytoplankton';'phytoplankton_nut'};
% 'bgc_NH4';'bgc_det_N';'bgc_det_P';
% FABM prefix for (sub)model

% settings
yl=365.25; dayl=24*3600; fs=22;
col=[[0.9 0.6 0.25];[0.65 0. 0.3];[0 0 0];[0.7 0.1 1];[0.2 0.7 0.3];[0.1 0.4 0.8];[0.2 0.52 0.95];];%

clear data;
ns=2;   % number of scenarios
% read series of netcdf result files to compare
for is=1:ns
  datf=['~/prog/tame/setup/0d/output' num2str(is) '.nc'];
  read_nc_simple
end
tim=datime/dayl;
% ----------------------------------------
% calculate plot partitioning
totn=length(varn);
nrow = round(sqrt(totn)); ncol=ceil(totn/nrow);
dxp = 1.01/(ncol+0.05); dyp = 1./(nrow +0.05);
x00 = 0.045; y00=0.05;
eps = 1E-3; % y-margin around min-max

% open figure
gcf=figure(1);
set(gcf,'Position',[400 00 160+ncol*360 100+nrow*270],'Visible','on','Color','w');clf;

% find position of first env variables
%%idxS = strfind(varn,'temp');   % find the index of your string
%%idx = min(find(not(cellfun('isempty',idxS))));

% loop over variables to plot
for i=1:totn
  % subplot position
  ix = mod(i-1,ncol);
  iy = floor((i-1)/ncol);

  gca=subplot('Position',[x00+ix*dxp y00+iy*dyp 0.8*dxp 0.86*dyp]);
  hold on;
  set(gca,'Box','on','YScale','Lin','FontSize',fs); %'Xlim',[25 52]*yl,'Ylim',[0 4.8],'YTick',0:4,

  % retrieve data from name 
  %%if i<idx, mname=modelname;
  %%else mname=''; end
  j=find(strcmp(vars,[varn{i}])); %mname 
  ymin=0;ymax=0;

  if ~isempty(j) % if name is found 
    for is=1:ns
      y=squeeze(data(is,j,:));
      lel(is)=plot(tim,y,'-','Color',col(is,:),'LineWidth',3);
    %  ym=mean(y); ys=std(y);
      ymin=min(ymin,min(y));
      ymax=max(ymax,max(y));
    end
  else
     fprintf('Error: variable %s not found in netcdf file!\n',[mname varn{i}])
  end

  % add title
  tmp=varn{i};
  ip=strfind(tmp,'_');
  if ip, tmp=tmp(ip+1:end); end
  tmp=strrep(tmp,'_',' ');
  annotation('textbox',[x00+(ix-0.05)*dxp y00+(iy+0.65)*dyp 0.2 0.05],'String',tmp,'Color','k','Fontweight','bold','FontSize',fs,'LineStyle','none','HorizontalAlignment','center');

  set(gca,'Box','on','Xlim',[min(tim) max(tim)],'Ylim',[ymin-eps ymax+eps]);%,'YTick',0:4,
  if ymax/(ymin+eps) > 1E4, set(gca,'YScale','log','Ylim',[ymin+eps ymax]);
  else set(gca,'Ylim',[ymin-eps ymax+eps]);
  end
  yl=ylabel(units{j});
  yp=get(yl,'Position'); yp(1)=yp(1)+max(tim)/35;
  set(yl,'Position',yp);
end %i

% add legend
le=legend(lel,num2str([1:ns]'),'location','northwest');%'
set(le,'Box','off','FontSize',fs);

%dtim=datime(end)-datime(1);
% output to original setup folder 
ii=findstr(datf,'/');
fnam=[datf(1:ii(end)) 'simres.png'];
fprintf('save PNG in %s ...\n',fnam);
print(gcf,'-dpng',fnam);
