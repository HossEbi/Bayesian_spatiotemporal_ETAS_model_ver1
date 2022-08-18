% Calculation of integral of rate of events with M=m in the interval [tstart,tend] for Likelihood Function 
% Mi & Ti  : the vectors of events
% written by: Hossein Ebrahimian   
% Last update: 01/2017

function intLambda = intLambdaETAS (Mi, Ti, ri, m, tstart, tend, Ml, K, theta, dA)

beta  = theta(1);
c     = theta(2);
p     = theta(3);
d     = theta(4);
q     = theta(5);

Krt = (q-1)/pi*d^(2*(q-1))*(p-1)*c^(p-1);

index = find(Ti < tend); 

if p == 1    
    Io = log((tend-Ti(index)+c)./(tstart-Ti(index)+c));     
else    
    Io = ((tend-Ti(index)+c).^(1-p)-(tstart-Ti(index)+c).^(1-p))/(1-p);
end

intLambda = beta*exp(-beta*(m-Ml))*sum((K*exp(beta*(Mi(index)-Ml)).*Io).*dA.*(Krt./(ri{index(end),1}.^2+d^2).^q));

end












        




    



















        