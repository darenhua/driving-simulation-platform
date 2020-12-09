% HSAVC simulator
% 
% input:    1) track data
%           2) initial state of car
%           3) physical state of camera
%           4) camera noise ratio
%           5) motor gain values
% 
% % output:   
%           1) display a figure of the car on the track
%           2) autonomously update position of the car and camera. 

clf, clear, clc;

%initialization, inputs

import_track = readtable('track_8_12-09-2020-00-42.csv');
s = table2struct(import_track);

show_grid = 1;
box on;
initial_state = [0 0 0];

aov = pi/3;
cam_angle = pi/6;
cam_height = 2;

snr = 33;

drive_gain = 1;
steering_gain = 1;

f1 = figure;
x0 = 10;

y0 = 10;
%zoom = min(max(1,log10(target_length)), 1.6);
%width = round(600*zoom);
width = 960;
%height = round(480*zoom);
height = 1000;
set(f1,'position',[x0,y0,width,height])

track_axes = subplot(2,1,1);
cam_axes = subplot(2,1,2);

cam_axes.XLim = [1 128];
cam_axes.YLim = [0 7];

pos1 = [0.1 0.3 0.3 0.3];
track_axes.Position = pos1;
pos2 = [0.5 0.15 0.4 0.7];
cam_axes.Position = pos2;
axes(track_axes);
title('Vehicle Sim');
axes(cam_axes);
title('Camera Sim');

hold on;

track_length = length(s);

%fixing up track re-imported data
for i = 1:track_length
    s(i).dir = [s(i).dir_1, s(i).dir_2];
    s(i).loc = [s(i).loc_1, s(i).loc_2];
    s(i).tile = [s(i).tile_1, s(i).tile_2];
end

minY = s(1).tile(2);
minX = s(1).tile(1);

for i=1:length(s)
    if s(i).tile(1)<minX
        minX = s(i).tile(1);
    end
end

for i=1:length(s)
    if s(i).tile(2)<minY
        minY = s(i).tile(2);
    end
end
%minY and minX should be bot left corner tile center
minX = minX-2;
minY = minY-2;
maxX = int16(10);
maxY = int16(10);
%want max to be number so guarentees clean multiples
multipleX = idivide(maxX, 4, 'floor');
multipleY = idivide(maxY, 4, 'floor');
maxX = (1 + multipleX) * 4;
maxY = (1 + multipleY) * 4;
track_axes.XLim = [minX maxX];
track_axes.YLim = [minY maxY];
gridlines(minX, maxX, minY, maxY, 4, 4, show_grid, track_axes);


i = 1;

axes(track_axes);
hold on;
[recreated_lines, inner, outer] = recreate_track(initial_state, s);
initial_dot = plot(initial_state(1), initial_state(2), '.', 'MarkerSize', 20, 'Color', 'm');
[intersects, cam_bounds, bound_lines, cam_line, intersect_lines] = generate_camera(aov, cam_angle, cam_height, initial_state, recreated_lines, track_axes);
[waveform, waveform_line] = generate_waveform(snr, intersects, cam_bounds, cam_axes);
pause(2);
centerIdx = CenterIdx(waveform);
[drive_motor, servo_motor] = motor_controller(centerIdx, drive_gain, steering_gain);
state = vehicle_model(servo_motor, drive_motor, initial_state);
delete(initial_dot);
delete(bound_lines);
delete(cam_line);
delete(intersect_lines);
delete(waveform_line);
axes(track_axes);
hold on;
dot = plot(state(end,1), state(end,2), '.', 'MarkerSize', 20, 'Color', 'm');
pause(2);

% simulator

while i<100
    car_state = state(end, 1:3);
    delete(inner);
    delete(outer);
    [recreated_lines, inner, outer] = recreate_track(car_state, s);
    [intersects, cam_bounds, bound_lines, cam_line, intersect_lines] = generate_camera(aov, cam_angle, cam_height, car_state, recreated_lines, track_axes);
    [waveform, waveform_line] = generate_waveform(snr, intersects, cam_bounds, cam_axes);
    pause(.5);
    centerIdx = CenterIdx(waveform);
    [drive_motor, servo_motor] = motor_controller(centerIdx, drive_gain, steering_gain);
    state = vehicle_model(servo_motor, drive_motor, car_state);
    delete(dot);
    delete(bound_lines);
    delete(cam_line);
    delete(intersect_lines);
    delete(waveform_line);
    axes(track_axes);
    hold on;
    dot = plot(state(end,1), state(end,2), '.', 'MarkerSize', 20, 'Color', 'm');
    pause(1);
    i = i+1;
end

function gridlines(minX, maxX, minY, maxY, intervalX, intervalY, show_grid, track_axes)
    axes(track_axes);
    if show_grid
        for row = minY: intervalY : maxY
            line([minX, maxY], [row, row], 'Color', 'r');
        end
        for col = minX : intervalX : maxX
            line([col, col], [minY, maxX], 'Color', 'r');
        end
    end
end

function [recreated_lines, inner, outer] = recreate_track(state, s)
% line re-generation tool
% 
% input:    1) section information structure
%           2) track sequence string
%           3) car state (position, camera position, angle)
% output:   1) display a figure of the track section closest to the car
%           2) coordinate array of recreated lines closest to car


car_pos = state(end, 1:2);

%identify track piece before, current, and future
%locate which part of track car is in
for i=1:length(s)
    tile_center = s(i).tile;
    tile_range = [tile_center(1)-2 tile_center(1)+2 ; tile_center(2)-2 tile_center(2)+2];
    if (tile_range(1,1) <= car_pos(1)) && (car_pos(1) <= tile_range(1,2)) && (tile_range(2,1) <= car_pos(2)) && (car_pos(2) <= tile_range(2,2)) 
        ct_index = i;
    end
end

current_tile = s(ct_index);

if ct_index+1 <= length(s)
    %check negative condition, under first
    next_tile = s(ct_index+1);
else
    next_tile = s(1);
end

if ct_index-1 > 0
    %check negative condition, positive first
    prev_tile = s(ct_index-1);
else
    prev_tile = s(end);
end

if ct_index+2 <= length(s)
    next_tile2 = s(ct_index+2);
elseif ct_index+2 == length(s)+1
    next_tile(2) = s(1);
else
    next_tile2 = s(2);
end

if ct_index-2 > 0
    prev_tile2 = s(ct_index-2);
elseif ct_index-2 == 0
    prev_tile2 = s(end);
else
    prev_tile2 = s(end-1);
end

recreated_lines = generate_joined_line(prev_tile2, prev_tile, current_tile, next_tile, next_tile2);
%inner
inner = plot(recreated_lines{1,1}, recreated_lines{1,2}, 'Color', 'b', 'LineWidth', 4);
%outer
outer = plot(recreated_lines{2,1}, recreated_lines{2,2}, 'Color', 'r', 'LineWidth', 4);


function joined_line = generate_joined_line(prev_tile2, prev_tile, current_tile, next_tile, next_tile2)
    %check what type (F/L/R) each section is
    tile_arr = [prev_tile2 prev_tile current_tile next_tile next_tile2];
    reconstructed = struct();
    for i=2:length(tile_arr)-1
        %don't want to actually iterate over prev_tile2, next_tile2... bad
        %solution though
        if tile_arr(i).weight == 4
            %F
            reconstructed(i-1).type = 'F';
            if tile_arr(i).dir(1) ~= 0
                if tile_arr(i).dir(1) > 0
                    x1_s1 = tile_arr(i-1).loc(1);
                    x2_s1 = tile_arr(i).loc(1);
                    y1_s1 = tile_arr(i).loc(2)+1;
                    y1_s2 = tile_arr(i).loc(2)-1;
                    x_s1 = linspace(x1_s1,x2_s1);
                    x_length = length(x_s1);
                    y_s1 = zeros(1,x_length) + y1_s1;
                    x_s2 = linspace(x1_s1,x2_s1);
                    y_s2 = zeros(1,x_length) + y1_s2;
                else
                    x1_s1 = tile_arr(i-1).loc(1);
                    x2_s1 = tile_arr(i).loc(1);
                    y1_s1 = tile_arr(i).loc(2)+1;
                    y1_s2 = tile_arr(i).loc(2)-1;
                    x_s1 = linspace(x1_s1,x2_s1);
                    x_length = length(x_s1);
                    y_s1 = zeros(1,x_length) + y1_s1;
                    x_s2 = linspace(x1_s1,x2_s1);
                    y_s2 = zeros(1,x_length) + y1_s2;                    
                end
            else
                if tile_arr(i).dir(2) > 0
                    x1_s1 = tile_arr(i).loc(1)-1;
                    x1_s2 = tile_arr(i).loc(1)+1;
                    y1_s1 = tile_arr(i-1).loc(2);
                    y2_s1 = tile_arr(i).loc(2);
                    y_s1 = linspace(y1_s1,y2_s1);
                    y_length = length(y_s1);
                    x_s1 = zeros(1,y_length) + x1_s1;
                    y_s2 = linspace(y1_s1,y2_s1);
                    x_s2 = zeros(1,y_length) + x1_s2;
                else
                    x1_s1 = tile_arr(i).loc(1)+1;
                    x1_s2 = tile_arr(i).loc(1)-1;
                    y1_s1 = tile_arr(i-1).loc(2);
                    y2_s1 = tile_arr(i).loc(2);
                    y_s1 = linspace(y1_s1,y2_s1);
                    y_length = length(y_s1);
                    x_s1 = zeros(1,y_length) + x1_s1;
                    y_s2 = linspace(y1_s1,y2_s1);
                    x_s2 = zeros(1,y_length) + x1_s2;
                end
            end
            reconstructed(i-1).xRangeS1 = x_s1;
            reconstructed(i-1).xRangeS2 = x_s2;
            reconstructed(i-1).yRangeS1 = y_s1;
            reconstructed(i-1).yRangeS2 = y_s2;
        elseif tile_arr(i).angle == 90
            %L
            reconstructed(i-1).type = 'L';
            if tile_arr(i).dir(1) ~= 0
                centerX = tile_arr(i).loc(1);
                if tile_arr(i).dir(1) > 0
                    centerY = tile_arr(i).loc(2) + 2;
                    x1_s1 = tile_arr(i-1).loc(1) + 1;
                    x2_s1 = centerX;
                    x1_s2 = tile_arr(i-1).loc(1) - 1;
                    x2_s2 = centerX;
                    x_s1 = linspace(x1_s1, x2_s1);
                    x_s2 = linspace(x1_s2, x2_s2);
                    th_s1 = -acos(x_s1-centerX);
                    th_s2 = -acos((x_s2-centerX)/3);
                    y_s1 = sin(th_s1) + centerY;
                    y_s2 = 3*sin(th_s2) + centerY;
                else
                    centerY = tile_arr(i).loc(2) - 2;
                    x1_s1 = centerX;
                    x2_s1 = tile_arr(i-1).loc(1) - 1;
                    x1_s2 = centerX;
                    x2_s2 = tile_arr(i-1).loc(1) + 1;
                    x_s1 = linspace(x1_s1, x2_s1);
                    x_s2 = linspace(x1_s2, x2_s2);
                    th_s1 = acos(x_s1-centerX);
                    th_s2 = acos((x_s2-centerX)/3);
                    y_s1 = sin(th_s1) + centerY;
                    y_s2 = 3*sin(th_s2) + centerY;
                end
            else
                centerY = tile_arr(i).loc(2);
                if tile_arr(i).dir(2) > 0
                    centerX = tile_arr(i).loc(1) - 2;
                    x1_s1 = centerX;
                    x2_s1 = tile_arr(i).loc(1) - 1;
                    x1_s2 = centerX;
                    x2_s2 = tile_arr(i).loc(1) + 1;
                    x_s1 = linspace(x1_s1, x2_s1);
                    x_s2 = linspace(x1_s2, x2_s2);
                    th_s1 = -acos(x_s1-centerX);
                    th_s2 = -acos((x_s2-centerX)/3);
                    y_s1 = sin(th_s1) + centerY;
                    y_s2 = 3*sin(th_s2) + centerY;
                else
                    centerX = tile_arr(i).loc(1) + 2;
                    x1_s1 = tile_arr(i).loc(1) + 1;
                    x2_s1 = centerX;
                    x1_s2 = tile_arr(i).loc(1) - 1;
                    x2_s2 = centerX;
                    x_s1 = linspace(x1_s1, x2_s1);
                    x_s2 = linspace(x1_s2, x2_s2);
                    th_s1 = acos(x_s1-centerX);
                    th_s2 = acos((x_s2-centerX)/3);
                    y_s1 = sin(th_s1) + centerY;
                    y_s2 = 3*sin(th_s2) + centerY;
                end
            end
            reconstructed(i-1).xRangeS1 = x_s1;
            reconstructed(i-1).xRangeS2 = x_s2;
            reconstructed(i-1).yRangeS1 = y_s1;
            reconstructed(i-1).yRangeS2 = y_s2;
        elseif tile_arr(i).angle == -90
            %R
            reconstructed(i).type = 'R';
            if tile_arr(i).dir(1) ~= 0
                centerX = tile_arr(i).loc(1);
                if tile_arr(i).dir(1) > 0
                    centerY = tile_arr(i).loc(2) - 2;
                    %duplicates, should probably just include outside of
                    %loop
                    x1_s1 = tile_arr(i-1).loc(1) + 1;
                    x2_s1 = centerX;
                    x_s1 = linspace(x1_s1, x2_s1);
                    x1_s2 = tile_arr(i-1).loc(1) - 1;
                    x2_s2 = centerX;
                    x_s2 = linspace(x1_s2, x2_s2);
                    th_s1 = acos(x_s1-centerX);
                    th_s2 = acos((x_s2-centerX)/3);
                    y_s1 = sin(th_s1) + centerY;
                    y_s2 = 3*sin(th_s2) + centerY;
                else
                    centerY = tile_arr(i).loc(2) + 2;
                    x1_s1 = centerX;
                    x2_s1 = tile_arr(i-1).loc(1) - 1;
                    x_s1 = linspace(x1_s1, x2_s1);
                    x1_s2 = centerX;
                    x2_s2 = tile_arr(i-1).loc(1) + 1;
                    x_s2 = linspace(x1_s2, x2_s2);
                    th_s1 = -acos(x_s1-centerX);
                    th_s2 = -acos((x_s2-centerX)/3);
                    y_s1 = sin(th_s1) + centerY;
                    y_s2 = 3*sin(th_s2) + centerY;
                end
            else
                centerY = tile_arr(i).loc(2);
                if tile_arr(i).dir(2) > 0
                    centerX = tile_arr(i).loc(1) + 2;
                    x1_s1 = tile_arr(i).loc(1) + 1;
                    x2_s1 = centerX;
                    x_s1 = linspace(x1_s1, x2_s1);
                    x1_s2 = tile_arr(i).loc(1) - 1;
                    x2_s2 = centerX;
                    x_s2 = linspace(x1_s2, x2_s2);
                    th_s1 = -acos(x_s1-centerX);
                    th_s2 = -acos((x_s2-centerX)/3);
                    y_s1 = sin(th_s1) + centerY;
                    y_s2 = 3*sin(th_s2) + centerY;
                else
                    centerX = tile_arr(i).loc(1) - 2;
                    x1_s1 = centerX;
                    x2_s1 = tile_arr(i).loc(1) - 1;
                    x_s1 = linspace(x1_s1, x2_s1);
                    x1_s2 = centerX;
                    x2_s2 = tile_arr(i).loc(1) + 1;
                    x_s2 = linspace(x1_s2, x2_s2);
                    th_s1 = acos(x_s1-centerX);
                    th_s2 = acos((x_s2-centerX)/3);
                    y_s1 = sin(th_s1) + centerY;
                    y_s2 = 3*sin(th_s2) + centerY;
                end
            end
            reconstructed(i-1).xRangeS1 = x_s1;
            reconstructed(i-1).xRangeS2 = x_s2;
            reconstructed(i-1).yRangeS1 = y_s1;
            reconstructed(i-1).yRangeS2 = y_s2;
        end
    end
    joined_lineX1 = [];
    joined_lineY1 = [];
    joined_lineX2 = [];
    joined_lineY2 = [];

    for i=1:length(reconstructed)
        joined_lineX1 = [joined_lineX1 reconstructed(i).xRangeS1];
        joined_lineX2 = [joined_lineX2 reconstructed(i).xRangeS2];
        joined_lineY1 = [joined_lineY1 reconstructed(i).yRangeS1];
        joined_lineY2 = [joined_lineY2 reconstructed(i).yRangeS2];
    end
    joined_line = {joined_lineX1, joined_lineY1; joined_lineX2, joined_lineY2};
end

end
function [intersects, cam_bounds, bound_lines, cam_line, intersect_lines] = generate_camera(aov, cam_angle, cam_height, state, recreated_lines, track_axes)

% camera line-of-sight generation tool
% 
% input:    1) recreated lines closest to car (prev, current, next)
%           2) track sequence as txt
%           3) car information (position, camera position, angle)
% output:   1) display a line of sight line 
%           2) intersects, if any, between recreated lines and sight line

% figure initialization


%track created
%car_angle = pi/6;

car_state = state(end, 1:3);
carX = car_state(1);
carY = car_state(2);
car_angle = state(end,3);

cam_bounds = find_cam_boundary(car_state, car_angle, cam_height, aov, cam_angle);
axes(track_axes)
bound1 = plot(cam_bounds{1,1}(1), cam_bounds{1,1}(2), '.', 'MarkerSize', 20, 'Color', 'm');
bound2 = plot(cam_bounds{1,2}(1), cam_bounds{1,2}(2), '.', 'MarkerSize', 20, 'Color', 'm');
bound3 = plot(cam_bounds{1,3}(1), cam_bounds{1,3}(2), '.', 'MarkerSize', 20, 'Color', 'm');
bound_lines = [bound1 bound2 bound3];
linescanX = linspace(cam_bounds{1,1}(1), cam_bounds{1,2}(1));
linescanY = linspace(cam_bounds{1,1}(2), cam_bounds{1,2}(2));
cam_line = plot(linescanX, linescanY, 'Color', 'm');
cam_inner = InterX([linescanX;linescanY],[recreated_lines{1,1};recreated_lines{1,2}]);
cam_outer = InterX([linescanX;linescanY],[recreated_lines{2,1};recreated_lines{2,2}]);
if ~isempty(cam_inner)
    intersect1 = plot(cam_inner(1), cam_inner(2), 'x', 'MarkerSize', 20, 'Color', 'k');
end
if ~isempty(cam_outer)
    intersect2 = plot(cam_outer(1), cam_outer(2), 'x', 'MarkerSize', 20, 'Color', 'k');
end    

intersects = {cam_inner, cam_outer};
if (exist('intersect1')) && (exist('intersect2'))
    intersect_lines=[intersect1, intersect2];
elseif exist('intersect1')
    intersect_lines = intersect1;
elseif exist('intersect2')
    intersect_lines = intersect2;
end 

function P = InterX(L1,varargin)
    %...Argument checks and assignment of L2
    error(nargchk(1,2,nargin));
    if nargin == 1,
        L2 = L1;    hF = @lt;   %...Avoid the inclusion of common points
    else
        L2 = varargin{1}; hF = @le;
    end
       
    %...Preliminary stuff
    x1  = L1(1,:)';  x2 = L2(1,:);
    y1  = L1(2,:)';  y2 = L2(2,:);
    dx1 = diff(x1); dy1 = diff(y1);
    dx2 = diff(x2); dy2 = diff(y2);
    
    %...Determine 'signed distances'   
    S1 = dx1.*y1(1:end-1) - dy1.*x1(1:end-1);
    S2 = dx2.*y2(1:end-1) - dy2.*x2(1:end-1);
    
    C1 = feval(hF,D(bsxfun(@times,dx1,y2)-bsxfun(@times,dy1,x2),S1),0);
    C2 = feval(hF,D((bsxfun(@times,y1,dx2)-bsxfun(@times,x1,dy2))',S2'),0)';

    %...Obtain the segments where an intersection is expected
    [i,j] = find(C1 & C2); 
    if isempty(i),P = zeros(2,0);return; end;
    
    %...Transpose and prepare for output
    i=i'; dx2=dx2'; dy2=dy2'; S2 = S2';
    L = dy2(j).*dx1(i) - dy1(i).*dx2(j);
    i = i(L~=0); j=j(L~=0); L=L(L~=0);  %...Avoid divisions by 0
    
    %...Solve system of eqs to get the common points
    P = unique([dx2(j).*S1(i) - dx1(i).*S2(j), ...
                dy2(j).*S1(i) - dy1(i).*S2(j)]./[L L],'rows')';
              
    function u = D(x,y)
        u = bsxfun(@minus,x(:,1:end-1),y).*bsxfun(@minus,x(:,2:end),y);
    end
end

function cam_bounds = find_cam_boundary(position, car_angle, cam_height, aov, cam_angle)
    %need the car + experimentation to find aov  
    d = cam_height/cos(cam_angle);
    s = cam_height*tan(cam_angle);
    change_mid_x = s*cos(car_angle);
    change_mid_y = s*sin(car_angle);
    theta = aov/2;
    change_left = d*tan(theta);
    change_left_y = change_left*cos(car_angle);
    change_left_x = change_left*sin(car_angle);
    %the right change from the center is the same as left, but negative
    if (0<=car_angle) && (car_angle<90)
        car_dir = [0, 1];
    elseif (90<=car_angle) && (car_angle<180)
        car_dir = [1,0];   
    elseif (180<=car_angle) && (car_angle<270)
        car_dir = [-1,0];   
    elseif (270<=car_angle) && (car_angle<360)
        car_dir = [1,0];   
    end
    cam_mid_x = position(1) + change_mid_x;
    cam_mid_y = position(2) + change_mid_y;
    x1 = cam_mid_x - change_left_x;
    y1 = cam_mid_y + change_left_y;
    x2 = cam_mid_x + change_left_x;
    y2 = cam_mid_y - change_left_y;
    cam_bounds = {[x1,y1], [x2,y2], [cam_mid_x,cam_mid_y]};
end
end
function [waveform, waveform_line] = generate_waveform(snr, intersects, cam_bounds, cam_axes)


% camera waveform generation tool
% 
% input:    1) camera line of sight information
%           2) intersects, if any, between recreated lines and sight line
% output:   1) camera waveform values
%           2) display camera waveform as figure

cam_inner = intersects{1,1};
cam_outer = intersects{1,2};

%right now, the program has the full camera line and the intersection
%points. We want the magnitude of the vector from each edge of the camera
%line to the closest intersection point. 

camera_lengths = calculate_distances(cam_inner, cam_outer, cam_bounds);
waveform = generateWaveform(camera_lengths);
noise = generateNoise(waveform, snr);
waveform = noise+waveform;
axes(cam_axes);
cam_pixels = linspace(1, 128, 128);
waveform_line = plot(cam_pixels, waveform);
function noise = generateNoise(waveform, snr)
    %the lower signal to noise ratio, the more noisy function
    waveLength = length(waveform);
    signalPower = (sum(abs(waveform).^2))/waveLength;
    noisePower =  signalPower/(10^(snr/10));
    %we want random amount of additon/subtraction: range (-1,1)
    randomArr = (2 .* rand(1,waveLength)) - 1;
    noise = sqrt(noisePower) * randomArr;
end

function camera_lengths = calculate_distances(cam_inner, cam_outer, cam_bounds)
    %future fix: need to have this function be able to handle if the camera
    %has no intersect but only sees track, not floor
    if isempty(cam_inner) && isempty(cam_outer)
        distanceL = 0;
        distanceR = 0;
        center = 0;
    elseif isempty(cam_inner)
        distanceL = 0;
        distanceR = distance(cam_outer, cam_bounds{1,2});
        %norm(cam_outer - cam_bounds{1,2});
        center = distance(cam_bounds{1,1},cam_outer);
        %norm(cam_bounds{1,1}-cam_outer);
    elseif isempty(cam_outer)
        distanceR = 0;
        distanceL = distance(cam_inner, cam_bounds{1,1}); 
        %norm(cam_inner - cam_bounds{1,1});
        center = distance(cam_bounds{1,2},cam_inner);
        %norm(cam_bounds{1,2}-cam_inner);
    else
        distanceL = distance(cam_inner ,cam_bounds{1,1});
        %norm(cam_inner - cam_bounds{1,1});
        distanceR = distance(cam_outer,cam_bounds{1,2});
        %norm(cam_outer - cam_bounds{1,2});
        center = distance(cam_inner,cam_outer);
        %norm(cam_inner - cam_outer);
    end
    camera_lengths = [distanceL, center, distanceR];
end

function d = distance(p1,p2)
    d = sqrt(((p2(1)-p1(1))^2)+((p2(2)-p1(2))^2));
end

% output track info
%output_track(track, track_length, s, show_loc, f1);
function [waveform, center] = generateWaveform(camera_lengths)
    waveform = zeros(1,128);
    cam_length = 0;
    %cam_px_length = 0;
    %cam_length is length of camera sight line in feet
    cam_px = [];
    for i = 1:length(camera_lengths)
        if ~isempty(camera_lengths(i))
            cam_length = cam_length + camera_lengths(i);
            %cam_px_length = cam_px_length + 1; 
        end
    end
    %cam_px = zeros(cam_px_length);
    for i = 1:length(camera_lengths)
        if ~isempty(camera_lengths(i))
            %future: preallocate
            cam_px = [cam_px (camera_lengths(i)*128)/cam_length];
        end
    end
    
    
    if camera_lengths(1) == 0 && camera_lengths(2) == 0 && camera_lengths(3) == 0
        waveform = waveform + 3.5;
    elseif camera_lengths(1) == 0
        right_px = int16((camera_lengths(3)*128)/cam_length);
        center_px = int16((camera_lengths(2)*128)/cam_length);
        waveform(1:center_px) = waveform(1:center_px) + 8.5;
        waveform(center_px+1:center_px+right_px) = waveform(center_px+1:center_px+right_px) + 3.5;
    elseif camera_lengths(3) == 0 
        left_px = int16((camera_lengths(1)*128)/cam_length);
        center_px = int16((camera_lengths(2)*128)/cam_length);
        waveform(1:left_px-1) = waveform(1:left_px-1) + 3.5;
        waveform(left_px:center_px+left_px) = waveform(left_px:center_px+left_px) + 8.5;
    else
        %no other scenarios, as this is based off of outputs from
        %calculate_distances function
        left_px = int16((camera_lengths(1)*128)/cam_length);
        right_px = int16((camera_lengths(3)*128)/cam_length);
        center_px = int16((camera_lengths(2)*128)/cam_length);
        waveform(1:left_px-1) = waveform(1:left_px-1) + 3.5;
        waveform(left_px:center_px+left_px) = waveform(left_px:center_px+left_px) + 6.5;
        waveform(center_px+left_px+1:center_px+left_px+right_px) = waveform(center_px+left_px+1:center_px+left_px+right_px) + 3.5;
        %start = int16(.1.*left_px);
        %stop = .1.*right_px;
        %stop = int16(center_px + left_px + stop);
        %waveform(start:left_px-1) = waveform(start:left_px-1) + 2.5;
        %waveform(left_px:center_px+left_px) = waveform(left_px:center_px+left_px) + 5.5;
        %waveform(center_px+left_px+1:stop) = waveform(center_px+left_px+1:stop) + 2.5;
        %waveform = waveform+1;
    end

end

end
function [centerIdx] = CenterIdx(waveform)
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
function [drive_motor, servo_motor] = motor_controller(centerIdx, drive_gain, steering_gain)
    drive_motor = drive_gain;
    servo_motor = (centerIdx - 64)*.4*steering_gain;
end
function state = vehicle_model(servo_motor, drive_motor, state)

% vehicle model
% 
% input:    1) motor values
% output:   1) the position and angle of a car after .5 seconds

driveL = drive_motor;
driveR = drive_motor;

kinematicModel = bicycleKinematics('WheelBase', 1);
tspan = 0:0.025:.5;
[v, steering_rate, steering_angle] = motor2kinematics(servo_motor, driveL, driveR);
inputs = [v steering_rate];
[t,state] = ode45(@(t,state)derivative(kinematicModel,state,inputs),tspan, state);    



function [v, steering_rate, steering_angle] = motor2kinematics(servo, driveL, driveR)
    driveAvg = (driveL + driveR)/2;
    v = driveAvg * 2;
    steering_rate = servo * .1;
    steering_angle = servo * 1; %in radians

end
end
