% Calculation of rate of events with M=m at t=time using ETAS model (mu_ETAS) for Likelihood Function 
% Mi & Ti  : the vectors of events
% written by: Hossein Ebrahimian   
% Last update: 01/2017

function lambda = lambdaETAS (Mi, Ti, ri, m, time, Ml, K, theta, dA)

beta  = theta(1);
c     = theta(2);
p     = theta(3);
d     = theta(4);
q     = theta(5);

index = find(Ti < time);    

Krt = (q-1)/pi*d^(2*(q-1))*(p-1)*c^(p-1);

lambda = beta*exp(-beta*(m-Ml))*sum((K*exp(beta*(Mi(index)-Ml))./((time-Ti(index)+c).^p)).*dA.*(Krt./(ri{index(end),1}.^2+d^2).^q));

end






        




    



















        