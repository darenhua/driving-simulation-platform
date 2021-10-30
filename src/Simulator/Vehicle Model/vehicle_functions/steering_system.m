function [v, steering_rate, steering_angle] = steering_system(servo, driveL, driveR)
% This function converts the motor values to the car's predicted velocity,
% rate of steering, and steering angle.
% ========================================================================
% INPUTS:
% - Servo motor strength
% - Drive motor strength (left and right drive motor)
% OUTPUTS:
% - Predicted velocity 
% - Predicted steering_rate 
% - Predicted steering_angle 

driveAvg = (driveL + driveR)/2;
v = driveAvg * 2;
steering_rate = servo * .1;
steering_angle = servo * 1; %in radians

end
