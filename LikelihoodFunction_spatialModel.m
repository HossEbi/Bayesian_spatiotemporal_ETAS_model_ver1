% Likelihood Function

function p_teta = LikelihoodFunction_spatialModel (data, theta, dA)

M      = data{1,1};
T      = data{1,2};
tstart = data{1,3};
Ml     = data{1,4};
r      = data{1,5};

time = T(2:end);
m    = M(2:end);

K = calculate_Kseq (M, T, tstart, Ml, theta);

lambda    = zeros(length(time),1);
intlambda = zeros(length(time),1);

for j=1:length(time)
    lambda(j) = lambdaETAS (M, T, r, m(j), time(j), Ml, K, theta, dA);
    if j == 1
        intlambda(j) = intLambdaETAS (M, T, r, m(j), 0.0, time(j), Ml, K, theta, dA);        
    else
        intlambda(j) = intLambdaETAS (M, T, r, m(j), time(j-1), time(j), Ml, K, theta, dA);
    end
end

intlambdaEnd = intLambdaETAS_mgrMl (M, T, time(end), tstart, Ml, K, theta);

p_teta = prod(lambda.*exp(-intlambda))*exp(-intlambdaEnd);

end


        




    



















        