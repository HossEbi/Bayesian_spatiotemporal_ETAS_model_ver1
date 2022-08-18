%%% This scripts uses MCMC simulation with Metropolis-Hastings procedure
%%% to estimate the probability P(Data|Model) following the
%%% Supplementary Figure S5 in Ebrahimian & Jalayer (2017) 
%%% Outputs:
%%% 1- The posterior samples for five main ETAS model parameters 
%%% 2- Plots of sampled posterior distributions for each chain and for each ETAS parameter 
%%% 3- Acceptance ratio of sampling

%% Initialize the Metropolis-Hastings sampler

numUP         = 5;           % number of Uncertain Parameters
BandWidth_threshold = NaN;   % for Kernel Density
Name_uncerParameter = {'\beta','c','p','d','q'};  % name of Uncertain Parameters
vec_uncerParameter = {vec_beta,vec_c,vec_p,vec_d,vec_q};
nbins = 80;                  % number of bins for plotting the histogram, can be larger

%% Define the PDF for Prior and Proposal Distributions 

if use_adaptive_prior == 0
%%% Prior    
priorPDF     = {'normal','normal','normal','normal','normal'};
priorPDF_par = {[beta_gen,beta_STR],[c_gen,c_STR],[p_gen,p_STR],[d_ini,d_STR],[q_ini,q_STR]};
priorPDFfunction = {@(x,par) normpdf(x,par(1),par(2)),...
                    @(x,par) normpdf(x,par(1),par(2)),...
                    @(x,par) normpdf(x,par(1),par(2)),...
                    @(x,par) normpdf(x,par(1),par(2)),...
                    @(x,par) normpdf(x,par(1),par(2)),...
                    };
                
%%% Proposal Distributions
proposalPDF     = {'lognormal','lognormal','lognormal','lognormal','lognormal'};
proposalPDF_par = {0.10,0.10,0.10,0.10,0.10};

%%% Initial Start of THETA
THETA = [beta_gen;c_gen;p_gen;d_ini;q_ini];                

elseif use_adaptive_prior == 2
%%% Prior    
priorPDF     = {'uniform','uniform','uniform','uniform','uniform'};
priorPDF_par = {[0,2],[0,0.1],[1,3],[0,2],[1,4]}; % (thetamin,thetamax)
priorPDFfunction = {@(x,par) unifpdf(x,par(1),par(2)),...
                    @(x,par) unifpdf(x,par(1),par(2)),...
                    @(x,par) unifpdf(x,par(1),par(2)),...
                    @(x,par) unifpdf(x,par(1),par(2)),...
                    @(x,par) unifpdf(x,par(1),par(2)),...
                    };
 
%%% Proposal Distributions
proposalPDF     = {'lognormal','lognormal','lognormal','lognormal','lognormal'};
proposalPDF_par = {0.10,0.10,0.10,0.10,0.10};

%%% Initial Start of THETA
THETA = [beta_gen;c_gen;p_gen;d_ini;q_ini];                


%%% fit to posterior distribution of the preceding forecasting interval
elseif use_adaptive_prior == 1
load([output_Dir,'\samples for prior.mat'],'samples')

%%% Prior 
priorPDF     = {'kernel','kernel','kernel','kernel','kernel'};
priorPDF_par = cell(1,numUP);
priorPDFfunction = cell(1,numUP);

%%% Proposal Distributions
proposalPDF     = {'kernel','kernel','kernel','kernel','kernel'};
proposalPDF_par = cell(1,numUP);

%%% Initial Start of THETA
THETA = zeros(numUP,1);

%%% Assignments
for n = 1:numUP
    priorPDF_par{1,n} = fitdist((samples(n,:))','kernel','support','positive');
    [fn,xn] = ksdensity(priorPDF_par{1,n}.InputData.data,'support','positive','function','pdf');
    priorPDFfunction{1,n} = [xn',fn'];
    proposalPDF_par{1,n} = priorPDF_par{1,n};
    THETA(n,1) = priorPDF_par{1,n}.random;
end
    
end

%% likelihoodFunction definition
likelihoodFunction = @LikelihoodFunction_spatialModel;
DATA = {M(indexSeq),time_T0(indexSeq),tstart,Ml_cat,r};

%% Start sampling
nchain = 1;    
samples_NC = zeros(numChain,numUP,maxIterations);accept_NC = zeros(numChain,numUP,maxIterations);
samples_Nch= zeros(numChain,numUP,maxIterations);

while nchain <= numChain
    
    display(['-------------- Chain Number = ',num2str(nchain)]) 
    
    state  = zeros(numUP,maxIterations);  % Storage space for our samples
    accept = zeros(numUP,maxIterations);  % Storage space for accept/reject decisions

    for iter = 1:maxIterations
       for n = 1:numUP
            if nchain==1
                [THETA,accept(n,iter)] = sample_posterior_MCMC(THETA,n,proposalPDF{1,n},proposalPDF_par{1,n},priorPDF{1,n},priorPDF_par{1,n},likelihoodFunction,DATA, dA); % MH algorithm
            accept_NC(nchain,n,iter)= accept(n,iter);
            else
                [THETA,accept(n,iter)] = sample_posterior_MCMC(THETA,n,'kernel',psample{1,n},priorPDF{1,n},priorPDF_par{1,n},likelihoodFunction,DATA, dA);
            accept_NC(nchain,n,iter)= accept(n,iter);
            end
        end        
        state(:,iter) = THETA; % storage samples of theta  
        samples_NC(nchain,:,iter)= THETA; % storage samples of theta in all N (herein=5) number of chains
    end
    
 %%% Samples we take for further analysis excluding burnin to consider the initial transient effect
    samples = state(:,burnin+1:thin:maxIterations); % discarding first few samples
   
    psample = cell(1,numUP);
    for n = 1:numUP
        psample{1,n} = fitdist((samples(n,:))','kernel','support','positive');
        if ~isnan(BandWidth_threshold)
            psample{1,n}.BandWidth = min(BandWidth_threshold,psample{1,n}.BandWidth);
        end
        THETA(n,1) = psample{1,n}.random;
    end
    samples_Nch(nchain,:,:)= THETA(n,1);
    samples_NchNB = samples_Nch(:,:,burnin+1:thin:maxIterations);
   
%% Plot sampled posterior distributions of ETAS parameters for each chain 
for n = 1:numUP 
    gcf=figure;
    h = [];
    nameh = {};
    [Xbins,PMF,h(1)] = plot_hist(samples(n,:),min(vec_uncerParameter{1,n}),max(vec_uncerParameter{1,n}),nbins,'default');
    nameh{1} = 'Sample Posterior by MCMC';
    hold on
    
    if use_adaptive_prior == 0
        PDFprior = priorPDFfunction{1,n}(Xbins,priorPDF_par{1,n});% func(inputs), func:@(x,par)normpdf(x,par(1),par(2))-->normpdf(x,mu,sigma) 
        h(2) = plot(Xbins,PDFprior/sum(PDFprior),'-ob','LineWidth',2);
    elseif use_adaptive_prior == 1
        PDFprior = priorPDFfunction{1,n}(:,2);    
        h(2) = plot(priorPDFfunction{1,n}(:,1),PDFprior/max(PDFprior)*max(PMF),'-ob','LineWidth',2); 
    elseif use_adaptive_prior == 2
        PDFprior = priorPDFfunction{1,n}(Xbins,priorPDF_par{1,n});
        h(2) = plot(Xbins,PDFprior/sum(PDFprior),'-ob','LineWidth',2);
    end        
    nameh{2} = 'prior';
    
    if nchain > 1 
        hold on
        [fn,xn] = ksdensity(psample{1,n}.InputData.data,'support','positive','function','pdf');
        h(3) = plot(xn,fn/max(fn)*max(PMF),'-','color',[0.50,0.50,0.50],'LineWidth',2);  
        nameh{3} = 'Proposal PDF';
    end
     hold off
    xlim([min(vec_uncerParameter{1,n}),max(vec_uncerParameter{1,n})])
    xlabel(Name_uncerParameter{1,n},'fontsize',30)
    ylabel('Probability Distribution','fontsize',30)
    legend(h,nameh,'fontsize',15)
    set(gca,'fontsize',22)    
    
% saveas(gcf,[output_Dir,'/',strcat(num2str(n),'-nP-',num2str(nchain),'-nC-sample.tiff')]) 
end
    if nchain < numChain
        close all
    end
  nchain = nchain + 1;
end

%% Acceptance ratio
disp(' ')
disp('MCMC outputs: ')
for ii=1:length(THETA)
    fprintf( ['Acceptance ratio for THETA(',num2str(ii),') = ',num2str(mean(accept(ii,:)),'%3.3f'),'\n']);
    AccR_m(ii)= mean(accept(ii,:))'; 
end
% save([output_Dir,'\meanAccR for Tstart_',num2str(T_start,3),'_day and dt_',num2str(T_end-T_start,3),'_day.mat'],'AccR_m')
