%% function for plotting histograms

function [Xbins,PMF,h] = plot_hist(Data,Xmin,Xmax,nbins,barWidth)

Xbins = linspace(Xmin,Xmax,nbins);
counts = hist(Data,Xbins);
PMF = counts/sum(counts);
if strcmp(barWidth,'default')
    h = bar(Xbins,counts/sum(counts));
    set(get(gca,'Children'),'FaceColor',[.95 .95 1])
else
    h = bar(Xbins,counts/sum(counts),str2num(barWidth)); 
end

