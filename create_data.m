function artificial_data = create_data(N,prob,wood_mean,wood_std,steel_mean,steel_std)
% Default parameters:
if nargin < 6 || isempty(steel_std)
    steel_std = 3.5*10^(-5);
end
if nargin < 5 || isempty(steel_mean)
    steel_mean = 0.0029;
end
if nargin < 4 || isempty(wood_std)
    wood_std = 2.3*10^(-19);
end
if nargin < 3 || isempty(wood_mean)
    wood_mean = 8.5*10^(-4);
end

artificial_data = wood_mean + wood_std*randn(N,N);
bullet = steel_mean + steel_std*randn();

for i = 1:N
    for j = 1:N
        if (i - (N+1)/2)^2 + (j - (N+1)/2)^2 > (N/2)^2
            artificial_data(i,j) = 0;
        end
    end
end

if rand() < prob
    radius = rand()*N/2;
    angle = rand()*2*pi;
    [x,y] = deal(radius*cos(angle),radius*sin(angle));
    artificial_data(ceil(N/2 + x),ceil(N/2 + y)) = bullet;
end
