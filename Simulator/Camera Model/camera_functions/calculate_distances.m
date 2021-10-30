function intersect_distances = calculate_distances(cam_inner, cam_outer, line_of_sight)
% This function calculates the distances between the left-most point of the
% line of sight and the left-most intersect, the distance between the
% left-most intersect and the right-most intersect, and the distance
% between the right-most intersect and the right-most point. 
% ========================================================================
% INPUTS:
% - Intersects between line of sight and closest track sections.
% - Line of Sight
% OUTPUTS:
% - Array with the three distance values mentioned. (If any intersect does 
%   not exist, the corresponding distance becomes 0)
function d = distance(p1,p2)
    d = sqrt(((p2(1)-p1(1))^2)+((p2(2)-p1(2))^2));
end

% If no intersects exist, set all distances to 0
if isempty(cam_inner) && isempty(cam_outer)
    distanceL = 0;
    distanceR = 0;
    center = 0;
% If the left intersect does not exist, the left most distance also
% cannot exist.
elseif isempty(cam_inner)
    distanceL = 0;
    distanceR = distance(cam_outer, line_of_sight{1,2});
    center = distance(line_of_sight{1,1},cam_outer);
% If the right intersect does not exist, the right most distance also
% cannot exist.
elseif isempty(cam_outer)
    distanceR = 0;
    distanceL = distance(cam_inner, line_of_sight{1,1}); 
    center = distance(line_of_sight{1,2},cam_inner);
% If all intersects exist, caclulate distances with distance formula.
else
    distanceL = distance(cam_inner ,line_of_sight{1,1});
    distanceR = distance(cam_outer,line_of_sight{1,2});
    center = distance(cam_inner,cam_outer);
end
% Concatenate into array
intersect_distances = [distanceL, center, distanceR];
end
