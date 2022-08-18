function [xmean,x50,x16,x84,x025,x975] = stat(x)

xmean = mean(x,2);

x50 = median(x,2);

N = size(x,2);

x025 = zeros(size(x50));
x16 = zeros(size(x50));
x84 = zeros(size(x50));
x975 = zeros(size(x50));

for i=1:size(x,1)
    x025(i) = interp1(1:N,sort(x(i,:)),0.025*N,'nearest');
    x16(i) = interp1(1:N,sort(x(i,:)),0.16*N,'nearest');
    x84(i) = interp1(1:N,sort(x(i,:)),0.84*N,'nearest');
    x975(i) = interp1(1:N,sort(x(i,:)),0.975*N,'nearest');   
end




