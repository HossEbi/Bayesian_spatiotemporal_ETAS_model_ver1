% Calculation of integral of rate of events with M=m in the interval [tstart,tend]
% Mi & Ti  : the vectors of events
% written by: Hossein Ebrahimian   
% Last update: 01/2017

function K = calculate_Kseq (Mi, Ti, tstart, Ml, theta)
%Mi=M(indexSeq);Ti=time_MS(indexSeq); tstart=tstart; Ml=Mc; theta=samples(:,j);


time = Ti(2:end);

intLambda = zeros(length(Ti),1);

intLambda(1) = intLambdaETAS_mgrMl (Mi, Ti, Ti(1), time(1), Ml, 1, theta); % Ti(1) was 0.0 in the original version

for j=2:length(time)
    intLambda(j) = intLambdaETAS_mgrMl (Mi, Ti, time(j-1), time(j), Ml, 1, theta);
end

intLambda(end) = intLambdaETAS_mgrMl (Mi, Ti, time(end), tstart, Ml, 1, theta);

K = length(Ti)/sum(intLambda);

end