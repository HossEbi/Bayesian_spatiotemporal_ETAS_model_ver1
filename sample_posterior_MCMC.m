%% Metropolis-Hastings algorithm
%%% Supplementary Figure S5 in  Ebrahimian & Jalayer (2017)
%%% Fig. 4 in Darzi et al. (2022)-GJI

function [THETA,accept] = sample_posterior_MCMC(THETA,rank,proposalPDF,proposalPDF_par,priorPDF,priorPDF_par,likelihoodFunction,data, dA)

theta = THETA(rank);

%% sampling new_theta from proposal PDF q(theta) / calculate Proposal Ratio q(old_theta)/q(new_theta)

if strcmp(proposalPDF,'normal')
    
    sigma_theta    = proposalPDF_par(1);
    new_theta      = normrnd(theta,sigma_theta);  
    proposal_ratio = 1.0;
    
elseif strcmp(proposalPDF,'lognormal')
    
    beta_theta     = proposalPDF_par(1);
    new_theta      = lognrnd(log(theta),beta_theta);   
    proposal_ratio = lognpdf(theta, log(new_theta), beta_theta)/lognpdf(new_theta, log(theta), beta_theta);
    
elseif strcmp(proposalPDF,'uniform')
    
    thetamin  = proposalPDF_par(1);
    thetamax  = proposalPDF_par(2);
    new_theta = unifrnd(thetamin,thetamax); 
    proposal_ratio = 1.0;
    
elseif strcmp(proposalPDF,'kernel')     
    
    psample   = proposalPDF_par;
    new_theta = psample.random;
    proposal_ratio = ksdensity(psample.InputData.data,theta,'function','pdf')/ksdensity(psample.InputData.data,new_theta,'function','pdf');

end

%% calculate the ratio of priors

if strcmp(priorPDF,'normal')
    meanpriorPDF   = priorPDF_par(1);
    sigmapriorPDF  = priorPDF_par(2);
    ratioPrior     = normpdf(new_theta, meanpriorPDF, sigmapriorPDF)/normpdf(theta, meanpriorPDF, sigmapriorPDF);
    
elseif strcmp(priorPDF,'lognormal')
    medianpriorPDF = priorPDF_par(1);
    betapriorPDF   = priorPDF_par(2);
    ratioPrior     = lognpdf(new_theta, log(medianpriorPDF), betapriorPDF)/lognpdf(theta, log(medianpriorPDF), betapriorPDF); 
    
elseif strcmp(priorPDF,'uniform')
    ratioPrior     = 1.0;
    
elseif strcmp(priorPDF,'kernel')
    ratioPrior     = ksdensity(priorPDF_par.InputData.data,new_theta,'function','pdf')/ksdensity(priorPDF_par.InputData.data,theta,'function','pdf');
    
end

%% calculate the likelihood function 

THETA(rank) = new_theta;
Likelihoodnew = likelihoodFunction(data,THETA, dA);
THETA(rank) = theta;
Likelihood    = likelihoodFunction(data,THETA, dA);

%% calculate the acceptance probability

pratio =  (Likelihoodnew/Likelihood)*ratioPrior;
alpha = min([1  pratio*proposal_ratio]); 
u = rand;                                    
if u <= alpha
    if (rank==3 && new_theta<1) || (rank==5 && new_theta<1)
        accept = 0;
    else
        theta = new_theta; 
        accept = 1;
    end
else
    accept = 0;
end
        
THETA(rank) = theta;
        