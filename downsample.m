function Y = downsample(X,factor)

N = size(X,1);
R = reshape(X,factor,N/factor,factor,N/factor);
t = sum(R, 1);
disp(size(t))
S = sum(sum(R, 1), 3) / (N/factor)^2;
disp(size(S))
Y = reshape(S, N/factor, N/factor);