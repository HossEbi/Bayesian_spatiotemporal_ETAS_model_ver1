# Bayesian_spatiotemporal_ETAS_model_ver1
 Bayesian spatio-temporal ETAS model Version 1

Guide to running the code:
This guide explains how to perform the code related to the seismicity forecasting based on a Bayesian epidemiological spatio-temporal aftershock clustering (ETAS) model. This code has been originally developed in MATLAB by Hossein Ebrahimian and Fatemeh Jalayer (University of Naples Federico II) in the following publication (the methodology is demonstrated by retrospective early forecasting of seismicity associated with the 2016 Central Italy seismic sequence activities):
---------------------------------------------
If you use this code, please cite:
Ebrahimian, H. & Jalayer, F. Robust seismicity forecasting based on Bayesian parameter estimation for epidemiological spatio-temporal aftershock clustering models. Scientific Reports 7, 9803 (2017), https://doi.org/10.1038/s41598-017-09962-z.
---------------------------------------------
The code has been recently updated by Atefe Darzi (University of Iceland) to perform retrospective early forecasting of seismicity for June 2000 South Iceland seismic sequence within the context of the TURNkey project (earthquake-turnkey.eu), supported by the European Union’s Horizon 2020 research and innovation programme under grant agreement No. 821046. 
The main script is called “main_script.mat”. The main input parameters, defined in this script, are as follows: 
•	lonMin, lonMax, latMin, latMax: The minimum and maximum longitude and latitude boundaries of the Aftershock zone.
•	Catalog: The text file containing the catalog of the ongoing seismic sequence which consists of 5 column vectors including [longitude, latitude, Magnitude, time, time_T0]. It is to note “time” is the time from the origin and “time_T0” is the time from the first event in the sequence.  
•	Mmax: Upper-bound Magnitude of the aftershock zone.
•	Ml_cat: Lower cut-off magnitude (which should be equal to or greater than the completeness magnitude).
•	vec_m: The vector containing the magnitude values m used for forecasting number of events with M≥m. It starts with Ml_cat.
•	T_start, T_end: Start and end time of the forecasting time interval of interest.
•	hour_ref: Time (hour) of issuing forecast automatically (e.g. 0.0).
•	beta_gen, c_gen, p_gen, d_ini, and q_ini: Initial mean values for prior distribution of 5 ETAS parameters [, c, p, d, q]. It is to note that the priors have a normal distribution.
•	CV_beta_gen, CV_c_gen, CV_p_gen, CV_d_gen, and CV_q_gen: CV (coefficient of variation) of the prior of 5 ETAS parameters (herein, we set it equal to 0.30).              
•	vec_beta, vec_c, vec_p, vec_d, vec_q, vec_K: the vector defining the boundaries of the ETAS parameters [, c, p, d, q, K] (to show the histogram of the posteriors and the PDF of the priors)
•	Regions: The shapefile (geographic data structure array), see for example the file 'S.IceTowns_prj.shp' located in “SIce_shp” folder for the Iceland test-case. 
•	ratio:  is equal to longitude/latitude associated with a square cell grid within the aftershock zone. It is noted that ratio varies for different geographic coordinate of the aftershock zone.
•	use_BkGd: background seismicity consideration; if use_BkGd=0, no Background seismicity is considered. If use_BkGd=1, the Background seismicity is considered by loading the file “sampleN_mgrM.mat” with the size equal to the 3D matrix [numer of grid cells, number of samples for the vector of model parameters, length of vector vec_m]. This file should be created by the user. The data should be the daily number of events with M≥m due to the background seismicity. In case of uniform background seismicity and for each value of m, the daily number of events in the column of the above matrix are equal (representing number of events with M≥m for each sample of model parameters). For more info, please see the main reference (Ebrahimian and Jalayer 2017).
•	use_adaptive_prior: define probability density function (PDF) for Prior and Proposal Distribution in MCMC, which can be 0, 1 or 2. if =0, normal prior with known/generic mean & STD and lognormal proposal are used (see [beta_gen,c_gen,p_gen,d_ini,q_ini] for mean priors and [beta_STR,c_STR,p_STR,d_STR, q_STR] for their STD); if =1, uniform priors and lognormal proposal (non-informative prior) are employed; if =2, use posterior samples from the preceding forecasting interval. To this end, the user needs to load the file “samples for prior.mat” (informative priors)
•	maxIterations: Number of posterior samples to draw
•	burnin: Burnin period, assumption: 0.1*maxIterations
•	numChain: Number of Markov chains, default =5 
•	thin: Thinning parameter, default =1. NOTE: The selected posterior samples are [burnin+1:thin:maxIterations]
•	output_Dir: Output directory address
•	Date_T0: Time of occurrence of the Main event in the sequence (dd/mm/yyyy)
The outputs are as follows:
(1)	Posterior samples for main ETAS model parameters [, c, p, d, q]
(2)	Plots of sampled posterior distributions for each chain and for each ETAS parameter 
(3)	The 2nd, 16th, 50th, 84th and 98th percentiles of the forecasted number of events with M≥m over each spatial grid cell (seismicity map) /and across the whole aftershock zone (error bar) 
(4)	Probability of experiencing an event with M≥m in the entire aftershock zone

---------------------------------------------
Other direct references that have implemented this Bayesian workflow are: 
1.	Azarbakht, A., Ebrahimian, H., Jalayer, F., & Douglas, J. Variations in hazard during earthquake sequences between 1995 and 2018 in western Greece as evaluated by a Bayesian ETAS model. Geophys. J. Int. 231, 1, 27-46 (2022), https://doi.org/10.1093/gji/ggac177 
2.	Darzi, A., Halldorsson, B., Hrafnkelsson, B., Ebrahimian, H., Jalayer, F., & Vogfjörð, K. S. Calibration of a Bayesian Spatio-temporal ETAS model to the June 2000 South Iceland seismic sequence. Geophys. J. Int. (2022), GJI-S-21-0737 (revision in review)
3.	Darzi, A., Halldorsson, B., Hrafnkelsson, B., & Vogfjörð, K. S. Short-term Bayesian ETAS spatiotemporal forecasting of the Ölfus 2008 earthquake sequence in Iceland. Tectonophysics, 229522 (2022), https://doi.org/10.1016/j.tecto.2022.229522.
---------------------------------------------


