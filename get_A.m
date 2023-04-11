function A = get_A(N,ntheta,p,d)
% Computes the coefficient matrix A (Radon transformation)
%
% inputs:
%   N - number of cells in a row (or in a column) of the region
%   p - number of parallel rays simulated for each angle
%   d - distance between the first and last parallel rays
%   ntheta - number of angles to shoot rays with.
% output: coefficient matrix A

% Create two vectors, rho and theta, with respectively p and ntheta
% different values but both of length p*ntheta, such that taking an 
% element from each with the same index always gives a different 
% combination.
rholinspace = linspace(-d/2,d/2,p);
thetalinspace = linspace(0,pi,ntheta);

k = 1:ntheta*p;
[kt,kr] = quorem(sym(k-1),sym(p));
rho(k) = rholinspace(kr+1);
theta = zeros(1,ntheta*p);
theta(k) = thetalinspace(kt+1);

% Run through the indexes i and j of the cells, and k of the simulated rays
A = zeros(ntheta*p,N^2);
for k = 1:ntheta*p
    for i = 1:N
        for j = 1:N
            if intersect_cell(i,j,1,N,rho(k),theta(k)) % only call the function get_length when the ray passes through the cell
                A(k,(i-1)*N+j) = get_length(i,j,1,N,rho(k),theta(k));
            end
        end
    end
end

