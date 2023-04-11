function bool_intersect = intersect_cell(i,j,a,N,rho,theta)
% Checks if the line l(rho,theta) passes through the cell C_{ij}
% inputs:
%   i - first cell index
%   j - first cell index
%   a - parameter describing the size of the region B
%   N - number of cells in a row of the region B
%   rho - minimal distance between the origin and the ray
%   theta - angle between the ray and the y-axis
% output:
%   bool_intersect = 1 if the ray passes through the cell C_{ij}, 0
%   otherwise.

bool_intersect = 0;

% Compute the corners of the cell
[x0,y0,x1,y1] = deal(-a+(2*a/N)*(j-1), -a+(2*a/N)*(i-1), -a+(2*a/N)*(j), -a+(2*a/N)*(i)); 

% The first case is when theta is in [0,pi/2]
if theta <= pi/2
    rho0 = cos(theta)*x0 + sin(theta)*y0;
    rho1 = cos(theta)*x1 + sin(theta)*y1;

    % The corresponding conditions
    if rho0 < rho && rho < rho1
        bool_intersect = 1;
    end
    if (theta == 0 || theta == pi/2) && rho == rho1
        bool_intersect = 1;
    end
    
end

% The second case is when theta is in (pi/2,pi]
if theta > pi/2
    rho0 = cos(theta)*x1 + sin(theta)*y0;
    rho1 = cos(theta)*x0 + sin(theta)*y1;
    
    % The corresponding conditions
    if rho0 < rho && rho < rho1
        bool_intersect = 1;
    end
    if theta == pi && rho == rho0
        bool_intersect = 1;
    end
end
