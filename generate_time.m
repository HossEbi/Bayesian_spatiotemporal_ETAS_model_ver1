% Generate Ti for SEQgen (THINNING METHOD)

function tgen = generate_time (tstart, lambdaMax)
%tstart = tgen; lambdaMax = lambda_max;
u = rand;

iat = -1/lambdaMax*log(1-u);

%iat = exprnd(1/lambdaMax);

tgen = tstart + iat;

end

