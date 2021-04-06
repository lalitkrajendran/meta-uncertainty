function [d, x, y] = calculate_shortest_distance_point_line(p1, p2, p3)
% Function to calculate distance between a point and a line.
% Equation taken from this link: http://paulbourke.net/geometry/pointlineplane/
% INPUTS:
% p1, p2: (x,y) co-ordinates of the endpoints of the line segment
% p3: (x, y) co-ordinates of the point to which the distance is to be computed
%
% OUTPUTS:
% d: shortest distance between the point and the line
% x, y: point on the line that corresponds to the shortest distance
%
% AUTHORS:
% Lalit Rajendran (lrajendr@purdue.edu)
%
% DATE:
% 05/01/2020

    % extract x, y co-ordinates
    x1 = p1(1); y1 = p1(2);
    x2 = p2(1); y2 = p2(2);
    x3 = p3(1); y3 = p3(2);

    % calculate relative distance of the shortest point along the line
    u = ((x3 - x1) * (x2 - x1) + (y3 - y1) * (y2 - y1)) / ((x2 - x1)^2 + (y2 - y1)^2);

    % calculate co-ordinates of the shortest daistnce point
    x = x1 + u * (x2 - x1);
    y = y1 + u * (y2 - y1);

    % calculate shortest distance
    d = sqrt((x - x3)^2 + (y - y3)^2);
end