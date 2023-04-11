function sijk = get_length(i,j,a,N,rho,theta)
% Calculates the length of the line segment of the ray l(rho,theta) lying
% inside the cell C_{ij}.
% inputs:
%   i - first cell index
%   j - first cell index
%   a - parameter describing the size of the region B
%   N - number of cells in a row of the region B
%   rho - minimal distance between the origin and the ray
%   theta - angle between the ray and the y-axis
% output:
%   sijk - the length of the segment s_{ij,k}

% Compute the corners of the cell
[x0,y0,x1,y1] = deal(-a+(2*a/N)*(j-1),-a+(2*a/N)*(i-1),-a+(2*a/N)*(j),-a+(2*a/N)*(i));

% For vertical and horizontal lines, the length equals the side of the cell
if theta == 0 || (theta == pi || theta == pi/2)
    sijk = 2*a/N;
else
    Points = zeros(2); % Points will store the coordinates of the points where the line enters and exits the cell
    liney0 = (rho-x0*cos(theta))/sin(theta); % y-coordinate of the line at x = x0
    if liney0 < y0
        Points(1,2) = y0;
        Points(1,1) = (rho-y0*sin(theta))/cos(theta);
    elseif liney0 > y1
        Points(1,2) = y1;
        Points(1,1) = (rho-y1*sin(theta))/cos(theta);
    else
        Points(1,:) = [x0;liney0];
    end

    liney1 = (rho-x1*cos(theta))/sin(theta);
    if liney1 < y0
        Points(2,2) = y0;
        Points(2,1) = (rho-y0*sin(theta))/cos(theta);
    elseif liney1 > y1
        Points(2,2) = y1;
        Points(2,1) = (rho-y1*sin(theta))/cos(theta);
    else
        Points(2,:) = [x1;liney1];
    end

    sijk = abs(norm(Points(1,:)-Points(2,:)));
end

