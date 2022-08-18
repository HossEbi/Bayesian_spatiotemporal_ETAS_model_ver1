% Calculation of rate of events with M=m at t=time for the whole area using ETAS model (mu_ETAS) for Sequence Generation
% Mi & Ti  : the vectors of events
% written by: Hossein Ebrahimian 
% Last update: 01/2017

function lambda = lambdaETAS_txy (Mi, Ti, m, time, Ml, K, theta)

beta  = theta(1);
c     = theta(2);
p     = theta(3);

index = find(Ti < time);    

Kt = (p-1)*c^(p-1);

lambda = Kt*beta*exp(-beta*(m-Ml))*sum(K*exp(beta*(Mi(index)-Ml))./((time-Ti(index)+c).^p));

end






        




    



















        