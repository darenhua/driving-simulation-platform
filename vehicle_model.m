% simulation
% 
% input:    1) when the car recieved steering power input
%           2) how much steering power input recieved over time
%           3) when the car recieves drive motor input
%           4) how much drive motor input recieved over time
%           5) track pieces csv file
%           6) track sequence txt file
%
%%% Can make a nice GUI with slider inputs for time/power: like a timeline
%
% output:   1) display a figure of the track
%           2) display racing line
%           3) display car's path
%           4) inputs as txt file/hardcode as txt
%           5) figure as a png file 
%           6) video of car's motion

clc, clear;

%initialize while loop
true = 1;
index = 1;
%motor input structure
inputs = struct('driveL', 1, 'driveR', 1, 'servo', 30, 'tspan', [0:0.05:1]);
%loop
while true
    fprintf('enter ''done'' at any time to finish input');
    formatSpec = 'input No.%d | Enter %s';
    inputVar = 'driveL value ([-1,1])';
    fprintf(formatSpec,index,inputVar);
    driveL = input('');
    %future update: make it so that user can only input int [-1,1], not
    %just 'not done'
    if string(driveL) ~= 'done'
        inputs(index).driveL = driveL;
    else
        break
    end
    inputVar = 'driveR value ([-1,1])';
    fprintf(formatSpec,index,inputVar);
    driveR = input('');
    if string(driveR) ~= 'done'
        inputs(index).driveR = driveR;
    else
        break
    end
    inputVar = 'servo value ([-30,30])';
    fprintf(formatSpec,index,inputVar);
    servo = input('');
    if string(servo) ~= 'done'
        inputs(index).servo = servo;
    else
        break
    end    
    inputVar = 'tspan value ([timeStart:timeStep:timeEnd])';
    fprintf(formatSpec,index,inputVar);
    tspan = input('');
    if string(tspan) ~= 'done'
        inputs(index).tspan = tspan;
    else
        break
    end
    index = index + 1;
end

%preparation: importing and organizing data
import_track = readtable('track_12_09-18-2020-00-11.csv');
s = table2struct(import_track);
track_length = length(s);
track = fileread('track_12_09-18-2020-00-11.txt');
%potential future update: use multiple figures with figure(H) to have
%multiple events going on simultaneously
f1 = figure;
x0 = 10;
y0 = 10;
%zoom = min(max(1,log10(target_length)), 1.6);
%width = round(600*zoom);
width = 960;
%height = round(480*zoom);
height = 768;
set(f1,'position',[x0,y0,width,height])
hold on;

%fixing up track re-imported data
for i = 1:track_length
    s(i).dir = [s(i).dir_1, s(i).dir_2];
    s(i).loc = [s(i).loc_1, s(i).loc_2];
    s(i).tile = [s(i).tile_1, s(i).tile_2];
end

%main code here
%recreate track
full_plot = plot_track(track, track_length, s, f1);
%initialize model
kinematicModel = bicycleKinematics('WheelBase', 1);
initialState = [0 0 0];

%initialize simulation
%TODO: Change into a function
tspan = 0:0.05:1;

%inputs [v psiDot] where
%v is the vehicle velocity in the direction of motion in meters per second.
% psiDot is the steering angle rate in radians per second.

inputs = [2 pi/4];

%solve
[t,position] = ode45(@(t,position)derivative(kinematicModel,position,inputs),tspan,initialState);

%plot
plot(position(:,1),position(:,2))

%function position = calculate_pos(input)
%end

function [v, psiDot] = motorSystem(servo, driveL, driveR)
    driveAvg = (driveL + drive R)/2;
    v = driveAvg * 2;
    psiDot = servo * .1;
end

%recreate track function
function full_plot = plot_track(track, track_length, s, f1)
fp1 = [];
fp2 = [];
radius_inner = 1;  
radius_outer = 3;

title({['\fontsize{16}Track Auto-Generation (Length = ', num2str(track_length), ')'] ['\fontsize{11}Track Sequence = ', track]});

for i = 1:track_length  
    %Tile_Center = s(i).tile
    if track(i) == 'F'
        if s(i).dir(1) ~= 0
            p1 = line([s(i).loc(1),s(i).loc(1)-s(i).dir(1)*4],[s(i).loc(2)+1,s(i).loc(2)+1], 'Color', 'black','LineWidth',4);
            p2 = line([s(i).loc(1),s(i).loc(1)-s(i).dir(1)*4],[s(i).loc(2)-1,s(i).loc(2)-1], 'Color', 'black','LineWidth',4);
        else
            p1 = line([s(i).loc(1)+1,s(i).loc(1)+1],[s(i).loc(2),s(i).loc(2)-s(i).dir(2)*4], 'Color', 'black','LineWidth',4);
            p2 = line([s(i).loc(1)-1,s(i).loc(1)-1],[s(i).loc(2),s(i).loc(2)-s(i).dir(2)*4], 'Color', 'black','LineWidth',4);            
        end
    else
        angle_index = atan2(s(i).dir(2), s(i).dir(1))/(pi/2)*2 + 3;
        th = linspace((angle_index+sin(s(i).angle*pi/180))*pi/4, (angle_index+2+sin(s(i).angle*pi/180))*pi/4, 100);
        x = radius_inner*cos(th) + s(i).loc(1) - 2*s(i).dir(2)*sin(s(i).angle*pi/180);
        y = radius_inner*sin(th) + s(i).loc(2) + 2*s(i).dir(1)*sin(s(i).angle*pi/180);
        p1 = plot(x,y,'Color', 'black','LineWidth',4); 
        cir_x = s(i).loc(1) - 2*s(i).dir(2)*sin(s(i).angle*pi/180);
        cir_y = s(i).loc(2) + 2*s(i).dir(1)*sin(s(i).angle*pi/180);
        x = radius_outer*cos(th) + cir_x;
        y = radius_outer*sin(th) + cir_y;
        %cir = [cir_x cir_y];
        p2 = plot(x,y,'Color', 'black','LineWidth',4); 
        %viscircles(cir, 0.05, 'Color', 'black');
        axis equal;
    end
    fp1 = horzcat(fp1, p1);
    fp2 = horzcat(fp2, p2);
end
    full_plot = [fp1, fp2];
drawnow;
end
