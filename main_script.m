% Simulation-based Bayesian Epidemic Type Aftershock Sequence (ETAS) model 
% Written by: Hossein Ebrahimian and Fatemeh Jalayer (2017)
% Updated by: Atefe Darzi to be applied to South Iceland seismicity

%% Authors
% Hossein Ebrahimian, PhD, University of Naples Federico II, hossein.ebrahimianchelehkhaneh@unina.it 
% Fatemeh Jalayer, PhD, University of Naples Federico II, fatemeh.jalayer@unina.it 
% Edited by: Atefe Darzi, PhD, University of Iceland, atefe@hi.is 
% DATE: 11/08/2022

clear all ; close all; clc

%% Control Parameters

set(0,'DefaultAxesFontName', 'Times New Roman')
use_BkGd           = 0;   % background seismicity 
use_adaptive_prior = 0;   % define PDF for Prior and Proposal Distribution in MCMC 

do_MCMC_updating             = 1;
do_find_Robust_estimate      = 1;
do_plot_ASzone_daily         = 1;

%% ************************* Inputs  ************************

%% Aftershock Zone, SISZ+RPOR in SW-Iceland

lonMin = -22.9; lonMax = -19.5;
latMin = 63.74; latMax = 64.23;
Mmax   = 7.0; % upper-bound magnitude according to the region

%% ETAS prior distributions
% Define the generic values of the mean prior distribution
beta_gen = 1.25; b_gen= beta_gen/log(10);
c_gen = 0.005; 
p_gen = 1.5; 
d_ini = 1; 
q_ini = 1.5; 

% Define CV and STD of prior distribution
% NOTE: we recommend that the initial prior distribution varies on a broad
% enough range of values.
CV_beta_gen  = 0.3;  beta_STR = CV_beta_gen *beta_gen;
CV_c_gen     = 0.3;  c_STR = CV_c_gen*c_gen;
CV_p_gen     = 0.3;  p_STR = CV_p_gen*p_gen;
CV_d_ini     = 0.3;  d_STR = CV_d_ini*d_ini;
CV_q_ini     = 0.3;  q_STR = CV_q_ini*q_ini;

% Vector of ETAS parameters for posterior plots  
% Note: the boundaries of the ETAS parameter vectors are recommended to be changed based on
% the sequence/region under study
vec_beta  = 0.01:0.01:2.2;
vec_c     = 0.001:0.001:0.012; 
vec_p     = 0.01:0.01:2.5;
vec_d     = 0.01:0.01:2;
vec_q     = 0.01:0.01:3;
vec_K     = 0.01:0.01:4.0;

% Inputs for MCMC simulation with Metropolis-Hastings procedure 
maxIterations = 200; 
burnin        = 0.1*maxIterations;          % burnin period
numChain      = 5;           % number of Markov chains
thin          = 1;           % thinning parameter (burnin+1:thin:maxIterations) 

%% Short-term Forecasting Interval [T_start,T_end]
hour_ref  = 0.0;    
T_start = 1; % from T0 
T_end   = 2; % one-day forecasting interval with Tstart
Date_T0='06/17/2000';

%% Read the catalog of an ongoing seismic sequence
% Note: the first earthquake of the sequence happened at T0
Catalog = load('Cat_2000Seq.txt');
  latitude  = Catalog(:,1); longitude = Catalog(:,2); M = Catalog(:,3);
  time      = Catalog(:,4); % origin time 
  time_T0   = Catalog(:,5); % time from T0

% Define Ml_cat (lower cut-off magnitude) 
   Ml_cat = 2; % Note that Ml_cat>= Mc
% magnitude vector for forecasting number of events 
   vec_m = [Ml_cat 3.5:1.0:6.5];
 
% organize the catalog : sort catalog in time order
    [time_T0,indx_sort_time] = sort(time_T0);
    time  = time(indx_sort_time); latitude  = latitude(indx_sort_time);
    longitude = longitude(indx_sort_time); M   = M(indx_sort_time);
    tstart = T_start- time(1); tend   = T_end- time(1);

% Use only earthquake data within the aftershock zone
    filter =  find(latitude>=latMin & latitude<=latMax & longitude>=lonMin & longitude<=lonMax & M>=Ml_cat); 
    latitude  = latitude(filter); longitude = longitude(filter);
    M = M(filter); time_T0 = time_T0(filter); time = time(filter);

%% Converting Aftershock Area Sources to ETAS grids
% to grid the aftershock zone into square cells 
ratio = 0.0151/0.00664; % ratio=longitude/latitude, for Iceland's geographic coordinate
deltaGrid_Y =  0.01; % latitude in degree
deltaGrid_X =  ratio*deltaGrid_Y; % longitude in degree

    % Area of each grid cell
    [dxgrid,AZ_x] = dis_az(latMin,lonMin,latMin,lonMin+deltaGrid_X);
    [dygrid,AZ_y] = dis_az(latMin,lonMin,latMin+deltaGrid_Y,lonMin);
    dA = dxgrid*dygrid; 
    % Define grid cells
    Xgrid = lonMin:deltaGrid_X:lonMax;
    Ygrid = latMin:deltaGrid_Y:latMax;
    % Define centre of each grid cell 
    Xcgrid = Xgrid(1:end-1)+deltaGrid_X/2; % long
    Ycgrid = Ygrid(1:end-1)+deltaGrid_Y/2; % lat
    
    num_grid = length(Xcgrid)*length(Ycgrid); %number of grids
    Ggrid = zeros(num_grid,2); count=1;
    for j=1:length(Ycgrid)
        for k=1:length(Xcgrid)
            Ggrid(count,:) = [Ycgrid(j),Xcgrid(k)];
            count=count+1;
        end
    end

%% Distance of each grid from Seqi events
   rxy = calculate_rxy(latitude,longitude,Ggrid(:,1:2)); 
   % size rxy = [length(Seq:event(lon,lat)),length(Ggrid)];

%%  ***********  Start Short-term Seismicity Forecasting *************
% output address
   output_Dir = 'DailyForecast_18Jun';
    if exist(output_Dir,'dir')==0 
        mkdir(output_Dir) % Make new folder
    end

%% MH-MCMC to sample posteriors of ETAS model parameters 
if do_MCMC_updating == 1  
   
indexSeq = find(time_T0 >= 0 & time <= T_start & M >= Ml_cat); % find seqi (the observation histroty)
    
    r = cell(length(indexSeq)-1,1);
    for j=2:length(indexSeq)
        Grdata  = [latitude(indexSeq(j)),longitude(indexSeq(j))]*pi/180;           % Read the aftershock latitude and longitude, gradi a radian
        Grpdata = [latitude(indexSeq(1:j-1)),longitude(indexSeq(1:j-1))]*pi/180;   % Read the aftershock latitude and longitude, gradi a radian    
        xv = topgeo(Grpdata(:,1),Grpdata(:,2),Grdata(1),Grdata(2));            
        r{j-1,1} = sqrt(xv(:,1).^2+xv(:,2).^2);
    end
   
 sample_posterior_STETAS_main

  % Calculate K,Kt,Kr
    K  = zeros(1,size(samples,2));
    c = samples(2,:); p = samples(3,:); d = samples(4,:); q = samples(5,:);
    for j=1:size(samples,2)
      K(j)  = calculate_Kseq (M(indexSeq),time_T0(indexSeq), tstart, Ml_cat, samples(:,j)); % first part of Eq7 in EJ2017
    end
      Kr = (samples(5,:)-1)/pi.*(samples(4,:).^(2*(samples(5,:)-1))); % Eq.5 in EJ2017
      Kt = (samples(3,:)-1).*(samples(2,:).^(samples(3,:)-1)); % Eq.4 in EJ2017
       
save([output_Dir,'\posterior_samples.mat'],'samples','Ml_cat','indexSeq','K','Kr','Kt')

% plot MCMC posterior samples for the last chain + posterior statistics 
plot_posterior_samples

    else
load([output_Dir,'\posterior_samples.mat'])
end    

%% Obtain Robust Estimates using an stochastic simulation procedure
if do_find_Robust_estimate == 1
      disp(' ')     
      disp('---------- Robust estimate for the number of events----------') 
      disp(' ')     
%%% load posterior samples obtained from the MH-MCMC algorithm
    load([output_Dir,'\posterior_samples.mat'])
    indexSeq = find(time_T0 >= 0 & time < T_start & M >= Ml_cat); % seg_i 

sampleN_mgrM = zeros(num_grid,size(samples,2),length(vec_m));
M_seqgen    = cell(1,size(samples,2)); time_seqgen = cell(1,size(samples,2));
Lon_seqgen  = cell(1,size(samples,2)); Lat_seqgen  = cell(1,size(samples,2));
r_seqgen    = cell(1,size(samples,2)); sum_intLambda_seqgen  = cell(1,size(samples,2));

for j = 1:size(samples,2) 
        display(['--------------Generated Seq no. ',num2str(j)]) 
    %%% Generate seqgen
        [M_seqgen{1,j},time_seqgen{1,j},Lon_seqgen{1,j},Lat_seqgen{1,j},r_seqgen{1,j},sum_intLambda_seqgen{1,j}] = ...
            generateSEQ (M(indexSeq),time_T0(indexSeq),rxy(indexSeq,:),tstart,tend,Ml_cat,K(j),samples(:,j),Xcgrid,Ycgrid,Ggrid,Mmax, dA);

    %%% Calculate Robust N
    for k=1:length(vec_m)
        for ngrid=1:num_grid            
            if ~isempty(M_seqgen{1,j})
                % Calculation of integral of rate of events with M=m in the interval [tstart,tend]
                sampleN_mgrM(ngrid,j,k) = calculate_N([M(indexSeq);M_seqgen{1,j}],[time_T0(indexSeq);time_seqgen{1,j}],...
                 [rxy(indexSeq,ngrid);r_seqgen{1,j}(:,ngrid)],vec_m(k),tstart,tend,Ml_cat,K(j),samples(:,j),sum_intLambda_seqgen{1,j}(:,ngrid), dA);
            else  % seqi 
                sampleN_mgrM(ngrid,j,k) = calculate_N(M(indexSeq),time_T0(indexSeq),rxy(indexSeq,ngrid),vec_m(k),tstart,tend,Ml_cat,K(j),samples(:,j),[], dA);
            end 
        end
    end
    
end

save([output_Dir,'\robust_estimate_N for Tstart_',num2str(T_start,3),'_day and dt_',num2str(T_end-T_start,3),'_day.mat'],...
        'M_seqgen','time_seqgen','Lon_seqgen','Lat_seqgen','r_seqgen','sum_intLambda_seqgen','sampleN_mgrM')
end

%% Post-Processing (Seismicity Forecasting Map)
 if do_plot_ASzone_daily == 1    
  
        if use_BkGd == 1
            load([output_Dir,'\robust_estimate_No.mat'],'sampleN0_mgrM')
            sampleNo_mgrM = sampleN0_mgrM*(T_end-T_start);
            load([output_Dir,'\robust_estimate_N for Tstart_',num2str(T_start,3),'_day and dt_',num2str(T_end-T_start,3),'_day.mat'])    
        else
            load([output_Dir,'\robust_estimate_N for Tstart_',num2str(T_start,3),'_day and dt_',num2str(T_end-T_start,3),'_day.mat'])    
            sampleNo_mgrM = zeros(size(sampleN_mgrM));
        end

% Calculate the Total Number of Events
    N_robust_mean = zeros(length(vec_m),1);
    N_robust_p50  = zeros(length(vec_m),1);
    N_robust_p16  = zeros(length(vec_m),1);
    N_robust_p84  = zeros(length(vec_m),1);
    N_robust_p02  = zeros(length(vec_m),1);
    N_robust_p98  = zeros(length(vec_m),1);
    N_observed    = zeros(length(vec_m),1);
% Calculate the Number of Events in each grid
    N_robust_mean_xy = zeros(length(Ycgrid),length(Xcgrid),length(vec_m));
    N_robust_p50_xy  = zeros(length(Ycgrid),length(Xcgrid),length(vec_m));
    N_robust_p16_xy  = zeros(length(Ycgrid),length(Xcgrid),length(vec_m));
    N_robust_p84_xy  = zeros(length(Ycgrid),length(Xcgrid),length(vec_m));
    N_robust_p02_xy  = zeros(length(Ycgrid),length(Xcgrid),length(vec_m));
    N_robust_p98_xy  = zeros(length(Ycgrid),length(Xcgrid),length(vec_m));

    P_mean = zeros(length(vec_m),1);    P_p98  = zeros(length(vec_m),1);
    P_p84  = zeros(length(vec_m),1);    P_p50  = zeros(length(vec_m),1);
    P_p16  = zeros(length(vec_m),1);    P_p02  = zeros(length(vec_m),1);

for k=1:length(vec_m)
    
 % Calculate N_robust for the whole aftershock zone
    [N_robust_mean(k),N_robust_p50(k),N_robust_p16(k),N_robust_p84(k),N_robust_p02(k),N_robust_p98(k)] = ordered_statistic(sum(sampleN_mgrM(:,:,k)+sampleNo_mgrM(:,:,k))); 
 % Calculate N_robust for each grid for spatial forecast
    [N_robust_mean_grid,N_robust_p50_grid,N_robust_p16_grid,N_robust_p84_grid,N_robust_p02_grid,N_robust_p98_grid] = ordered_statistic(sampleN_mgrM(:,:,k)+sampleNo_mgrM(:,:,k));
   
% spatial distribution of forecast_N_xy
        N_robust_mean_xy(:,:,k) = (reshape(N_robust_mean_grid,length(Xcgrid),length(Ycgrid)))';
        N_robust_p50_xy(:,:,k)  = (reshape(N_robust_p50_grid ,length(Xcgrid),length(Ycgrid)))';
        N_robust_p16_xy(:,:,k)  = (reshape(N_robust_p16_grid ,length(Xcgrid),length(Ycgrid)))';
        N_robust_p84_xy(:,:,k)  = (reshape(N_robust_p84_grid ,length(Xcgrid),length(Ycgrid)))';
        N_robust_p02_xy(:,:,k)  = (reshape(N_robust_p02_grid ,length(Xcgrid),length(Ycgrid)))';
        N_robust_p98_xy(:,:,k)  = (reshape(N_robust_p98_grid ,length(Xcgrid),length(Ycgrid)))';
         
% Calculate the probability over the entire aftershock zone
        P_mean(k) = 1-exp(-N_robust_mean(k));
        P_p98(k)  = 1-exp(-N_robust_p98(k));
        P_p84(k)  = 1-exp(-N_robust_p84(k));
        P_p50(k)  = 1-exp(-N_robust_p50(k));
        P_p16(k)  = 1-exp(-N_robust_p16(k));
        P_p02(k)  = 1-exp(-N_robust_p02(k));

    % number of observed earthquakes occurred during the forecasting interval of interest 
       N_observed(k) = length(find(time >= T_start & time < T_end & M >= vec_m(k)));

end 
 
save([output_Dir,'\N_robust_stat for Tstart_',num2str(T_start,3),'_day and dt_',num2str(T_end-T_start,3),'_day.mat'],...
        'P_mean','P_p98','P_p84','P_p50','P_p16','P_p02','N_robust_mean','N_robust_p50','N_robust_p16','N_robust_p84','N_robust_p02','N_robust_p98')
 
%%% Plot the Short-term Aftershock Forecasting Map with an Errorbar and their comparison with the observed data 
    plot_ForecastMap
    plot_N  % error bar

 end

%% Updating prior knowledge adaptively (day-by-day) during the ongoing sequence
% to inform the priors of the next forecasting interval
 for ii=1:5
    THETA_m(ii)= mean(samples(ii,:))'; 
 end

beta_gen  = THETA_m(1);  b_gen = beta_gen / log(10); 
c_gen     = THETA_m(2);
p_gen     = THETA_m(3);
d_ini     = THETA_m(4);
q_ini     = THETA_m(5);

% NOTE: we recommend choosing large enough standard deviation enabling the
% prior distribution covering a reasonable range of values
beta_STR  = str2double(num2str(std(samples(1,:)),'%3.4f'));  
c_STR     = str2double(num2str(std(samples(2,:)),'%3.5f'));
p_STR     = str2double(num2str(std(samples(3,:)),'%3.4f'));
d_STR     = str2double(num2str(std(samples(4,:)),'%3.4f'));
q_STR     = str2double(num2str(std(samples(5,:)),'%3.4f')); 
    



