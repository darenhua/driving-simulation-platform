% motor controller
% 
% input:    1) waveform information as txt
%           2) gain values
% output:   1) motor values as txt

clc, clf, clear;
servo = 0;
r_drive = 0;
l_drive = 0;
steer_gain = 1;
motor_gain = 1;
%can abstract above gain even more by starting from potentiometer input
waveform_str = fileread('cam_11-12-2020-21-18.txt');
waveform = str2num(waveform_str);
center = CenterIdx(waveform);
%so much potential with motor control gui:
%3d model of potentiometer knobs control gain...
steering_angle = (center - 64)*steer_gain;
%coded 1d lookup table
%steering_angle = fixpt_interp1([-30 0 0 30], [-30 -2 2 30], center, sfix(8), 2^-3., sfix(16), 2^-14, 'Floor');
%drive_motor = fixpt_interp1([1 1 0.5], [0 10 30], abs(steering_angle), sfix(8), 2^-3., sfix(16), 2^-14, 'Floor');
drive_motor = 1;
drive_motor = drive_motor*motor_gain;
driveL = drive_motor;
driveR = drive_motor;
output_motor(driveL, driveR, steering_angle)

kinematicModel = bicycleKinematics('WheelBase', 1);
tspan = 0:0.025:.5;
[v, steering_rate, steering_angle] = motor2kinematics(servo, driveL, driveR);
initialState = [0 0 steering_angle]; %x y angle

inputs = [v steering_rate];
[t,state] = ode45(@(t,position)derivative(kinematicModel,position,inputs),tspan,initialState);
output_position(initialState, inputs, state);


function [v, steering_rate, steering_angle] = motor2kinematics(servo, driveL, driveR)
    driveAvg = (driveL + driveR)/2;
    v = driveAvg * 2;
    steering_rate = servo * .1;
    steering_angle = servo * 1; %in radians
end

function output_position(initialState, inputs, state)
output_filename = strcat('pos',  datestr(now, '_mm-dd-yyyy-HH-MM')); 

track_filename = strcat(output_filename, '.txt');
outfile = fopen(track_filename,'w');
fprintf(outfile, '%s\r\n', mat2str(initialState));
fprintf(outfile, '%s\r\n', mat2str(inputs));
fprintf(outfile, '%s\r\n', mat2str(state));

end
function output_motor(driveL, driveR, steering_angle)
output_filename = strcat('motor',  datestr(now, '_mm-dd-yyyy-HH-MM')); 

track_filename = strcat(output_filename, '.txt');
outfile = fopen(track_filename,'w');
output = [driveL, driveR, steering_angle];
fprintf(outfile, mat2str(output));

end

function centerIdx = CenterIdx(waveform)
    first_edge = single(1);
    second_edge = single(128);
    threshold = .75*mean(waveform(40:80));
    for i=64:-1:2
        if waveform(i)<threshold
            first_edge = single(i);
            break;
        end
    end
    for i=64:1:127
        if waveform(i)<threshold
            second_edge = single(i);
            break;
        end
    end
    if first_edge>1 && second_edge<128
        centerIdx = (second_edge+first_edge)/2;
    else
        centerIdx = single(255);
    end
end
