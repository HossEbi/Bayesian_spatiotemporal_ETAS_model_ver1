%   This m-file is for Ploting MCMC Samples
%   Final Revision by: Hossein Ebrahimian 06/2016
%   Modified by Atefe Darzi- July 2022


plotXlabel = 1;    % do plot X labels    
plotPrior  = 1;    % do plot Priors
nbins = 100;        % number of bins for the histogram

gcf=figure;

set(0,'DefaultAxesFontName', 'Times New Roman')
set(gcf,'Position',[50,500,2000,850/3])

MCMC_statistics = cell(1,numUP);

count=0;
for n = 1:numUP
    subplot(1,numUP+1,n)
    [Xbins,~,h1] = plot_hist(samples(n,:),min(vec_uncerParameter{1,n}),max(vec_uncerParameter{1,n}),nbins,'default'); % sample posterior by MCMC
    hold on
    if plotPrior
        if use_adaptive_prior == 0
            PDFprior = priorPDFfunction{1,n}(Xbins,priorPDF_par{1,n});% func(inputs), func:@(x,par)normpdf(x,par(1),par(2))-->normpdf(x,mu,sigma) 
             h2 = plot(Xbins,PDFprior/sum(PDFprior),'-ob','markerSize',2,'LineWidth',2);  % PMF = P/sum(P)
        elseif use_adaptive_prior == 1
            PDFprior = priorPDFfunction{1,n}(:,2);    
             h2 = plot(priorPDFfunction{1,n}(:,1),PDFprior/max(PDFprior)*max(PMF),'-ob','markerSize',2,'LineWidth',2); 
        elseif use_adaptive_prior == 2
            PDFprior = priorPDFfunction{1,n}(Xbins,priorPDF_par{1,n});
             h2 = plot(Xbins,PDFprior/sum(PDFprior),'-ob','markerSize',2,'LineWidth',2);  % PMF = P/sum(P)
         end
    end
        
%     xlim([min(min(samples(n,:)),min(vec_uncerParameter{1,n})),max(max(vec_uncerParameter{1,n}),max(samples(n,:)))]) %
     xlim([min(vec_uncerParameter{1,n}),max(vec_uncerParameter{1,n})])
     
   set(gca,'fontsize',20)
   
   if plotXlabel
        xlabel(Name_uncerParameter{1,n},'fontsize',30)%,'fontweight','bold'
   end
    
%     if n==1
%         ylabel('PMF','fontsize',30) 
%     elseif (n==2 && plotPrior==1)    
%          legend([h1 h2],{'Posterior','Prior'},'fontsize',9,'location','northeast')
%     end
    
    set(gca,'YTick',[])
     
    disp(' ') 
    display( ['Mean(',Name_uncerParameter{1,n},') = ',num2str(mean(samples(n,:)),'%3.3f')]); 
    display( ['CV (',Name_uncerParameter{1,n},') = ',num2str(std(samples(n,:))/mean(samples(n,:)),'%3.4f')]); 
    
    [~,~,x16,x84,~,~] = ordered_statistic(samples(n,:)); %% [xmean,x50,x16,x84,x02,x98] 
    display( ['CV (',Name_uncerParameter{1,n},') = ',num2str(0.5*(x84-x16)/mean(samples(n,:)),'%3.4f')]); 
       
    [~,~,~,~,x025,x975] = stat(samples(n,:)); 
    display( ['95% CI (',Name_uncerParameter{1,n},') = [',num2str(x025,'%3.3f'),',',num2str(x975,'%3.3f'),']']); 

    count=count+1;
    MCMC_statistics{1,count} = ['mean ',Name_uncerParameter{1,n}];
    MCMC_statistics{2,count} = str2double(num2str(mean(samples(n,:)),'%3.4f'));
    MCMC_statistics{3,count} = str2double(num2str(mean(samples(n,:)),'%3.4f'));

    count=count+1;
    MCMC_statistics{1,count} = ['CV ',Name_uncerParameter{1,n}];
    MCMC_statistics{2,count} = str2double(num2str(std(samples(n,:))/mean(samples(n,:)),'%3.4f'));
    MCMC_statistics{3,count} = str2double(num2str(0.5*(x84-x16)/mean(samples(n,:)),'%3.4f'));

    count=count+1;
    MCMC_statistics{1,count} = ['95% CI ',Name_uncerParameter{1,n}];
    MCMC_statistics{2,count} = str2double(num2str(x025,'%3.4f')); %LL 
    MCMC_statistics{3,count} = str2double(num2str(x975,'%3.4f')); %UL 
end

subplot(1,numUP+1,numUP+1)
[Xbins,~,~] = plot_hist(K,min(vec_K),max(vec_K),nbins,'default');
xlim([min(vec_K),max(vec_K)])
set(gca,'fontsize',20)
legend([h1 h2],{'Posterior','Prior'},'fontsize',13,'location','northeast')
if plotXlabel
    xlabel('K','fontsize',30)
end
set(gca,'YTick',[])

saveas(gcf,[output_Dir,'\','Thetta.tiff']) 

%% K
disp(' ') 
display( ['Mean(K) = ',num2str(mean(K),'%3.3f')]); 
display( ['CV (K) = ',num2str(std(K)/mean(K),'%3.3f')]); 
[~,~,x16,x84,~,~] = ordered_statistic(K);
display( ['CV (K) = ',num2str(0.5*(x84-x16)/mean(K),'%3.3f')]); 
[~,~,~,~,x025,x975] = stat(K); 
display( ['95% CI (K)= [',num2str(x025,'%3.3f'),',',num2str(x975,'%3.3f'),']']); 
    
count=count+1;
MCMC_statistics{1,count} = 'mean K';
MCMC_statistics{2,count} = str2double(num2str(mean(K),'%3.3f'));
MCMC_statistics{3,count} = str2double(num2str(mean(K),'%3.3f'));
count=count+1;
MCMC_statistics{1,count} = 'CV K';
MCMC_statistics{2,count} = str2double(num2str(std(K)/mean(K),'%3.3f'));
MCMC_statistics{3,count} = str2double(num2str(0.5*(x84-x16)/mean(K),'%3.3f'));
    count=count+1;
    MCMC_statistics{1,count} = '95% CI K';
    MCMC_statistics{2,count} = str2double(num2str(x025,'%3.4f')); %LL 
    MCMC_statistics{3,count} = str2double(num2str(x975,'%3.4f')); %UL 

%% Kt
disp(' ') 
display( ['Mean(Kt) = ',num2str(mean(Kt),'%3.3f')]); 
display( ['CV (Kt) = ',num2str(std(Kt)/mean(Kt),'%3.3f')]); 
[~,~,x16,x84,~,~] = ordered_statistic(Kt);
display( ['CV (Kt) = ',num2str(0.5*(x84-x16)/mean(Kt),'%3.3f')]); 
[~,~,~,~,x025,x975] = stat(Kt); 
display( ['95% CI (Kt)= [',num2str(x025,'%3.3f'),',',num2str(x975,'%3.3f'),']']); 
    
count=count+1;
MCMC_statistics{1,count} = 'mean Kt';
MCMC_statistics{2,count} = str2double(num2str(mean(Kt),'%3.3f'));
MCMC_statistics{3,count} = str2double(num2str(mean(Kt),'%3.3f'));
count=count+1;
MCMC_statistics{1,count} = 'CV Kt';
MCMC_statistics{2,count} = str2double(num2str(std(Kt)/mean(Kt),'%3.3f'));
MCMC_statistics{3,count} = str2double(num2str(0.5*(x84-x16)/mean(Kt),'%3.3f'));
    count=count+1;
    MCMC_statistics{1,count} = '95% CI Kt';
    MCMC_statistics{2,count} = str2double(num2str(x025,'%3.4f')); %LL 
    MCMC_statistics{3,count} = str2double(num2str(x975,'%3.4f')); %UL 

%% Kr
disp(' ') 
display( ['Mean(Kr) = ',num2str(mean(Kr),'%3.3f')]); 
display( ['CV (Kr) = ',num2str(std(Kr)/mean(Kr),'%3.3f')]);
[~,~,x16,x84,~,~] = ordered_statistic(Kr);
display( ['CV (Kr) = ',num2str(0.5*(x84-x16)/mean(Kr),'%3.3f')]); 
[~,~,~,~,x025,x975] = stat(Kr); 
display( ['95% CI (Kr)= [',num2str(x025,'%3.3f'),',',num2str(x975,'%3.3f'),']']); 
    
count=count+1;
MCMC_statistics{1,count} = 'mean Kr';
MCMC_statistics{2,count} = str2double(num2str(mean(Kr),'%3.3f'));
MCMC_statistics{3,count} = str2double(num2str(mean(Kr),'%3.3f'));
count=count+1;
MCMC_statistics{1,count} = 'CV Kr';
MCMC_statistics{2,count} = str2double(num2str(std(Kr)/mean(Kr),'%3.3f'));
MCMC_statistics{3,count} = str2double(num2str(0.5*(x84-x16)/mean(Kr),'%3.3f'));
    count=count+1;
    MCMC_statistics{1,count} = '95% CI Kr';
    MCMC_statistics{2,count} = str2double(num2str(x025,'%3.4f')); %LL 
    MCMC_statistics{3,count} = str2double(num2str(x975,'%3.4f')); %UL 

save([output_Dir,'\posterior_stat.mat'],'MCMC_statistics')
