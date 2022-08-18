% Generate Mi, Ti for SEQgen (for Ti, using THINNING METHOD)

function [Mgen,tgen,numRepeat] = generate_m_time (Mi, Ti, ri, tstart, tend, Ml, K, theta, Mmax)
%K=sampleK; theta=sampletheta; 

Mgen = generate_M (Ml, Mmax, theta(1));

%%% Estimate lambda_max
% lambda_max = 0;
% for ngrid=1:size(ri,2)
%     lambda_max = lambda_max+lambdaETAS_xy(Mi,Ti,ri(:,ngrid),Mgen,tstart,Ml,K,theta);
% end
lambda_max = lambdaETAS_txy(Mi,Ti,Mgen,tstart,Ml,K,theta);

accept    = 0;
numRepeat = 0;

tgen = tstart;

while (accept == 0 && tgen <= tend)
    
    numRepeat = numRepeat + 1;
   
    %%% generate time
    tgen = generate_time (tgen, lambda_max);
    
    %%% Estimate lambdaxy
%     lambda = 0;
%     for ngrid=1:size(ri,2)
%         lambda = lambda+lambdaETAS_xy(Mi,Ti,ri(:,ngrid),Mgen,tgen,Ml,K,theta);
%     end    
    lambda = lambdaETAS_txy(Mi,Ti,Mgen,tgen,Ml,K,theta);
    
    %%% Thinning Algorithm
    ratio = lambda/lambda_max;
    u = rand;
    if u <= ratio
        accept = 1;
    end
    
end

end
