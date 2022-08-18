% Plot the forcasted seismicity for the Aftershock Zone
% written by: Hossein Ebrahimian 
% Last update: 06.2020 by Atefe Darzi

   
%% Input
T_start_figure = hour_ref ; 
ni= 1 ;% find(vec_m==Ml_cut)
forecast_data =  N_robust_p98_xy(:,:,ni); vec_m(ni) %  N_robust_p50_xy
% Prob of earthquake for the entire aftershock zone
Plb = P_mean ; % P_p98 , P_p84 , P_p50 , P_p16 , P_p02
 
caxis_max = 0.05; % max(max(forecast_data)) 0.05 / 0.2
Dxx = [0.0,-0.1,0.1,-0.1,0.1,-0.1];   
Dyy = [-0.15,0.10,-0.1,0.10,-0.1,0.10]; 
dxx = [0.02,0.01,0.01,0.01,0.01,0.01];
dyy = [0.00,0.01,0.01,0.01,0.01,0.01]; 
X_upperText = lonMin+2*0.01;% 13.22;
Y_upperText = latMax+3*0.01;% 42.50-0.05/2;

%% Plot Forecasts N
gcf=figure;
set(gcf,'Position',[200   80   1400   700])

pcolor(Xcgrid,Ycgrid,forecast_data);
labels = ['Forecasted Seismicity for Earthquakes with M', num2str(vec_m(ni)),'^{+}'];
hcb = colorbar; % t=get(hcb,'Limits'); set(hcb,'Ticks',linspace(t(1),t(2),5)) 
caxis([0.0 caxis_max])
% colorbar('Ticks',[0.01,0.02,0.06,0.08,0.1],'TickLabels',{'0.01','0.02','0.04','0.06','0.08','0.1'}) 
cmap = colormap(hot(256)); % jet(5)
cmap = flipud(cmap); %Flip array up to down
colormap(cmap)
ylabel(hcb,labels,'FontSize',22,'Rotation',-270,'FontWeight','normal'); hcb.Label.Position(1) = 5;
shading interp
hold on

%% Plot Region
Regions = shaperead('SIce_shp/S.IceTowns_prj.shp'); % geographic data structure array in projected map coordinates
vec_Regions   = 1:42; % number of polygones of the shapefile
% vec_Locations = [13,42.75;13.30,42.80;13.05,42.55;13.30,42.50]; % province capital:[lon, lat (j); ... ; lon,lat (j)]
for j=1:length(vec_Regions)
    jj = vec_Regions(j);
    plot(Regions(jj).X,Regions(jj).Y,'b','linewidth',1.3)
    hold on
end
                   
%% plot observations in the forecasting interval
indx_catalog = find(time >= T_start & time < T_end & M >= Ml_cat); 
if Ml_cat < 2.0
    [h,nameh] = plot_ASzone(indx_catalog,longitude,latitude,M,1.0,max(M),1,...
    {[0.7,0.9,0.3],'c','m','g','b','r','y'},{[0.5,0.5,0.5],[0.2,0.2,0.2],'k','k','k','k'},{'^','o','s','^','o','p','p'},[10,12,15,17,20,28,30],0);% 
else
    [h,nameh] = plot_ASzone(indx_catalog,longitude,latitude,M,Ml_cat,max(M),1,...
    {'c','m','g','b','r','y'},{[0.5,0.5,0.5],[0.2,0.2,0.2],'k','k','k','k'},{'o','s','^','o','p','p'},[12,15,17,20,28,30],0);% [0.8,0.8,0.8]
end
hold on

% add mag and date of large events on the map: "M6.5 and 30/10/2016" 
indxMg6 = find(M(indx_catalog)>=6.0);
% h1 = []; nameh1 = {};
longitude_i = longitude(indx_catalog); latitude_i = latitude(indx_catalog); M_i = M(indx_catalog);
day_i = day(indx_catalog); month_i = month(indx_catalog); year_i = year(indx_catalog);
if ~isempty(indxMg6)
   for k=1:length(indxMg6)
        j=indxMg6(k);
%         h1 = plot(longitude(j),latitude(j),'h','MarkerEdgeColor','b','MarkerFaceColor','none','MarkerSize',24);% AD: or do_plotMS=1 in plot_ASzone
%         nameh1 = ['M \geq ',num2str(6.0)];
        hh = annotation('arrow');
        set(hh,'parent',gca,'position',[longitude_i(j),latitude_i(j),Dxx(k),Dyy(k)],'linewidth',2,'color','r');
        text(longitude_i(j)+Dxx(k)+dxx(k),latitude_i(j)+Dyy(k)+dyy(k),...
            {['M',num2str(M_i(j),'%3.1f')];[num2str(day_i(j)),'/',num2str(month_i(j)),'/',num2str(year_i(j))]},...
            'FontName','Times New Roman','FontSize',25,'FontWeight','normal','Color','r');
   end
end

% learning and forecasting interval info
if T_end-T_start >= 1
    deltaT_figure = T_end-T_start;
    deltaT_figure_scale = '[day]';
else
    deltaT_figure = (T_end-T_start)*24;
    deltaT_figure_scale = '[hours]';
end 
  
for k= 2:length(vec_m)
    if Plb(k) >= 0.10
        text(lonMin+2*deltaGrid_X,latMax-.02-0.03*(k-2),['Pr(M\geq',num2str(vec_m(k)),') = ',num2str(Plb(k),'%2.2f')],'FontName','Times New Roman','FontSize',22,'FontWeight','normal')
    elseif (Plb(k) >= 0.010 && Plb(k) < 0.10)
        text(lonMin+2*deltaGrid_X,latMax-.02-0.03*(k-2),['Pr(M\geq',num2str(vec_m(k)),') = ',num2str(Plb(k),'%3.2f')],'FontName','Times New Roman','FontSize',22,'FontWeight','normal')
    else
        text(lonMin+2*deltaGrid_X,latMax-.02-0.03*(k-2),['Pr(M\geq',num2str(vec_m(k)),') = ',num2str(Plb(k),'%10.1e')],'FontName','Times New Roman','FontSize',22,'FontWeight','normal')
    end
end

h1=rectangle('Position',[lonMin,latMin+.02,lonMax-lonMin-.3,latMax-latMin-.02],'FaceColor','none','EdgeColor',[.5 .5 .5],'LineWidth',4);
% h1=rectangle('Position',[lonMin,latMin,lonMax-lonMin,latMax-latMin],'FaceColor',[.5 .5 .5],'EdgeColor',[1 1 1]);
% set(h1,'facealpha',.9)

text(X_upperText,Y_upperText,{['{\it T}_{start}= ',num2str(T_start_figure,'%4.2f'),'UTC ',...
    datestr(addtodate(datenum(Date_T0,'mm/dd/yy'),T_start,'day'),'dd-mmm-yy')...
    ' Forecasting interval= ',num2str(deltaT_figure,3),' ',deltaT_figure_scale]},'FontName','Times New Roman','FontSize',22,'FontWeight','normal')

% text(0.5*(lonMax+lonMin)-50*deltaGrid_X,latMin-.08,{['\itT\rm_{start}= ',num2str(T_start_figure,'%4.2f'),' UTC ',...
%     datestr(addtodate(datenum(refer_Date,'mm/dd/yyyy'),T_start,'day'),'dd-mmm-yyyy')];...
%     ['Forecasting interval= ',num2str(deltaT_figure,3),' ',deltaT_figure_scale]},'FontName','Times New Roman','FontSize',24,'FontWeight','normal','FontAngle','Normal','Color','k')

xlabel('Longitude','fontsize',24,'FontWeight','normal')
ylabel('Latitude','fontsize',24,'FontWeight','normal')
xlim([lonMin-0.1 lonMax])
ylim([latMin-0.02 latMax+0.07])

set(gca,'fontsize',24)
set(gca,'XTick',str2double(num2str(min(xlim),'%4.1f')):0.4:str2double(num2str(max(xlim),'%4.1f')),...
    'YTick',str2double(num2str(min(ylim),'%4.1f')):0.1:str2double(num2str(max(ylim),'%4.1f')),'GridLineStyle',':')
if ~isempty(h)
    legend(h,nameh,'fontsize',23)%,'FontWeight','bold'
end
% box on
% grid on

saveas(gcf,[output_Dir,'/','forecast_NMag',num2str(vec_m(ni)),'-N98Pmean.tif']) 

