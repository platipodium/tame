




% 
% matlab script for visualizing and comparing 0D model results (from netcdf files)
%
% kai wirtz (hereon 2024-2025)
%
clear all; close all
% name of variables to plot
varn={'bgc_NO3';'bgc_din';'bgc_PO4';'phyto_rate';'phyto_Q_N'; 'bgc_dom_P'; 'phyto_dQ_dt_N'; ...
  'phyto_phytoplankton_C'; 'phyto_nut1'; 'phyto_nut2';'phyto_phy_P';'phyto_dQ_dt_P';...
  'phyto_Q_P';'total_phosphorus_calculator_result'; 'phyto_phy_N'; ...
 'total_nitrogen_calculator_result';};
% 'phyto_Q_N';'bgc_det_N';'bgc_dom_C';'temp'; 'par' 'bgc_rate';'bgc_dom_N';'bgc_det_C';'bgc_RHS_NO3';
% 'bgc_det_P'; 'bgc_RHS_NH4';'bgc_dom_N';'bgc_NH4';'bgc_RHS_PO4';
links={'phyto_nut1';};  % variables to plot together with subsequent entry
% FABM prefix for (sub)model

% settings
yl=365.25; dayl=24*3600; fs=22;
col=[[0.9 0.6 0.25];[0.65 0. 0.3];[0 0 0];[0.7 0.1 1];[0.2 0.7 0.3];[0.1 0.4 0.8];[0.2 0.52 0.95];];%
lsty={'-';'--';':';':-'};
clear data;
ns=2;   % number of scenarios
% read series of netcdf result files to compare
for is=1:ns
  datf=['~/prog/tame/setup/0d/output' num2str(is-0) '.nc'];
  read_nc_simple
end
tim=datime/dayl;
% ----------------------------------------
% calculate plot partitioning
totn=length(varn);
nrow = round(sqrt(totn)); ncol=ceil(totn/nrow);
dxp = 1.01/(ncol+0.05); dyp = 0.98/(nrow +0.05);
x00 = 0.045; y00=0.04;
eps = 1E-3; % y-margin around min-max
dt =tim(2)-tim(1);
% open figure
gcf=figure(1);
set(gcf,'Position',[180 05 160+ncol*360 100+nrow*270],'Visible','on','Color','w');clf;

% find position of first env variables
%%idxS = strfind(varn,'temp');   % find the index of your string
%%idx = min(find(not(cellfun('isempty',idxS))));
link0=-1;
% loop over variables to plot
ig=1;
for i=1:totn
  if i~= link0+1
    % subplot position
    ix = mod(ig-1,ncol);
    iy = floor((ig-1)/ncol);
    gca=subplot('Position',[x00+ix*dxp y00+iy*dyp 0.8*dxp 0.88*dyp]);
    hold on;
    set(gca,'Box','on','YScale','Lin','FontSize',fs); 
    st=1;
    ig=ig+1;
  else
    st=st+1;
  end
  j=find(strcmp(vars,varn{i})); %mname 

  ymin=9E9;ymax=0;
  if ~isempty(j) % if name is found 
    jl=find(strcmp(links,varn{i})); %linked to subsequent variable 
    if ~isempty(jl), link0=i; end

    for is=1:ns
      y=squeeze(data(is,j,:));
      lel(is)=plot(tim,y,lsty{st},'Color',col(is,:),'LineWidth',3);
    %  ym=mean(y); ys=std(y);
      ymin=min(ymin,min(y));
      ymax=max(ymax,max(y));
      if strcmp('phyto_Q_N',varn{i})
         dy=diff(y)/dt;
         time2=tim(2:end)-dt;
      end
      if strcmp('bgc_PO4',varn{i})
         ddin=diff(y)/dt;
         time2=tim(2:end)-dt;
      end
      if strcmp('phyto_dQ_dt_N',varn{i})
        plot(time2,dy,'-','Color','k','LineWidth',2);
      end    
      if strcmp('phyto_nut2',varn{i})
        plot(time2,ddin,'-','Color','k','LineWidth',2);
      end      
    end
  else
     fprintf('Error: variable %s not found in netcdf file!\n',[varn{i}])
  end

  tval=[];%[22.925 22.93];
  for t=tval
    plot(t*ones(2,1),[ymin-eps ymax+eps],'k:','LineWidth',1)
  end
  % add title
  tmp=varn{i};
  ip=strfind(tmp,'_');
  if ip, tmp=tmp(ip(1)+1:end); end
  tmp=strrep(tmp,'_calculator_result','');
  tmp=strrep(tmp,'_',' ');
  annotation('textbox',[x00+(ix-0.05)*dxp y00+(iy+0.65-(st-1)*0.2)*dyp 0.2 0.05],'String',tmp,'Color','k','Fontweight','bold','FontSize',fs,'LineStyle','none','HorizontalAlignment','center');
%  set(gca,'Box','on','Xlim',[min(tim) max(tim)],'Ylim',[ymin-eps ymax+eps]);%,'YTick',0:4,
  set(gca,'Box','on','Xlim',[7 11],'Ylim',[ymin-eps ymax+eps]);%,'YTick',0:4,
  %if(ymin<0 & isempty(strfind(varn{i},'rate'))) ymin=1.5*eps; end
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
