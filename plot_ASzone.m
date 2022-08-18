% plot Mainshock + Aftershocks by different symbols on the geographic map
% plotted by plot_AS_spatial.m
%%% AD: [h,nameh] = plot_ASzone(1:length(M),longitude,latitude,M,3.0,max(M),1.0,{[0.80,0.80,0.80],'b','m','r'},{[0.50,0.50,0.50],'k','k','k'},{'o','o','^','p'},[4,10,12,20],0);

function [h,nameh] = plot_ASzone(index,longitude,latitude,M,Ml,Mms,dM,markerColor,markerEdgeColor,markerType,markerSize,do_plotMS)

longitude = longitude(index);
latitude  = latitude(index);
M         = M(index);

vecM = Ml:dM:Mms;

h = [];
nameh = {};
count = 0;
for j=1:length(vecM)-1
    index_vecMjplus  = find (M >= vecM(j) & M < vecM(j+1));
    if ~isempty(index_vecMjplus)
        count = count+1;
        h(count) = plot(longitude(index_vecMjplus),latitude(index_vecMjplus),markerType{1,j},...
            'MarkerEdgeColor',markerEdgeColor{1,j},'MarkerFaceColor',markerColor{1,j},'MarkerSize',markerSize(j));
        nameh{count} = [num2str(vecM(j)),' \leq M < ',num2str(vecM(j+1))];
        hold on
    end
end

%% plot the Mmax using yellow star symbol
% in case vecM(end)< Mms   6<6.5
j=length(vecM);
index_vecMjplus  = find (M >= vecM(end));
if ~isempty(index_vecMjplus)
    count = count+1;
%     for ii=1:length(index_vecMjplus)
        h(count) = plot(longitude(index_vecMjplus),latitude(index_vecMjplus),markerType{1,j},...
            'MarkerEdgeColor',markerEdgeColor{1,j},'MarkerFaceColor',markerColor{1,j},'MarkerSize',markerSize(j));
        nameh{count} = ['M \geq ',num2str(vecM(end))];
        hold on
%     end
end

% Note: catalog is sorted in time-order started by mainshock, followed by aftershocks
% MS=main shock
longitude_MS = longitude(1);
latitude_MS  = latitude(1);
m_MS         = M(1); 

if do_plotMS  %plot the first lat,long in catalog which is the mainshock
    hold on
    h(count+1) = plot(longitude_MS,latitude_MS,'h','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',20);
    nameh{count+1} = ['Main-shock, M=',num2str(m_MS)];
end