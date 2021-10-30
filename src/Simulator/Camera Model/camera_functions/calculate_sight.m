function line_of_sight = calculate_sight(position, car_angle, cam_height, aov, cam_angle)
% This function finds the left most, right most and center point of the
% virtual camera's line of sight.
% ========================================================================
% INPUTS:
% - Car XY Position
% - Car Angle
% - Camera height
% - Camera XY Plane Angle
% - Camera YZ Plane Angle of View
% OUTPUTS:
% - Cell Array with XY coordinates left most, right most and center point 
%   of the virtual camera's line of sight.
d = cam_height/cos(cam_angle);
n = cam_height*tan(cam_angle);
change_mid_x = n*cos(car_angle);
change_mid_y = n*sin(car_angle);
theta = aov/2;
change_left = d*tan(theta);
change_left_y = change_left*cos(car_angle);
change_left_x = change_left*sin(car_angle);

cam_mid_x = position(1) + change_mid_x;
cam_mid_y = position(2) + change_mid_y;
x1 = cam_mid_x - change_left_x;
y1 = cam_mid_y + change_left_y;
x2 = cam_mid_x + change_left_x;
y2 = cam_mid_y - change_left_y;

%left most point, right most point, and center point.
line_of_sight = {[x1,y1], [x2,y2], [cam_mid_x,cam_mid_y]};
end
