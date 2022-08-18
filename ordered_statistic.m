function [xmean,x50,x16,x84,x02,x98] = ordered_statistic(x)

xmean = mean(x,2);

x50 = median(x,2);

N = size(x,2);

x02 = zeros(size(x50));
x16 = zeros(size(x50));
x84 = zeros(size(x50));
x98 = zeros(size(x50));

for i=1:size(x,1)
    x02(i) = interp1(1:N,sort(x(i,:)),0.02*N,'nearest');
    x16(i) = interp1(1:N,sort(x(i,:)),0.16*N,'nearest');
    x84(i) = interp1(1:N,sort(x(i,:)),0.84*N,'nearest');
    x98(i) = interp1(1:N,sort(x(i,:)),0.98*N,'nearest');   
end




