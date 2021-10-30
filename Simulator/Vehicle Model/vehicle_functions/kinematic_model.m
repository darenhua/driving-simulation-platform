function next_state = kinematic_model(state, v, steering_rate, steering_angle)
% This function calculates the next position and angle of the virtual car
% after applying a given velocity, steering rate, and steering angle over a
% given time-step. Uses MATLAB's builtin bicycle kinematics model.
% ========================================================================
% INPUTS:
% - The car's initial position and angle (state)
% - Velocity
% - Steering Rate
% - Angle
% OUTPUTS:
% - The car's final position and angle after .05 seconds

%initialize the model
kinematicModel = bicycleKinematics('WheelBase', 1);
tspan = 0:0.0025:.05;
inputs = [v steering_rate];
%This line solves the kinematic model differential equations over the given
%timespan to turn velocity and steering rate into displacement and angle.
[t,next_state] = ode45(@(t,state)derivative(kinematicModel,state,inputs),tspan, state);    
end

