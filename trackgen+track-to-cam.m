% track generation tool
% 
% input:    1) selection of mode of track generation (manual or auto mode)
%           2) section sequence (F-straight, L-left, R-right) (manual mode)
%           3) no. of sections in the track (auto mode)
% output:   1) display a figure of the track
%           2) track sequence (auto mode) as a txt file
%           3)  section info (types/locations/directions) as a csv file
%           4) section info (types/locations/directions) as a mat file
%           5) track figure as a png file 

clc, clear;
frequency = 0; % value range [0,5]
show_process = 1;  % 1: show section locations in the figure, 0: don't show
show_loc = 1;  % 1: show section locations in the figure, 0: don't show
show_grid = 1;
full_plot = [];
defaultAx = get(gca);

% ask user to input the target track length
%target_length = input('Input Target Track Length (Even Number, Greater than 4): ');
target_length = 8;
% define the initial 4-segment track loop
track = 'LLLL';
track_length = 4;
s(1) = struct('dir', [ 0  1], 'weight', 2, 'angle', 90, 'loc', [ 2 2], 'tile', [ 2 0]);
s(2) = struct('dir', [-1  0], 'weight', 2, 'angle', 90, 'loc', [ 0 4], 'tile', [ 2 4]);
s(3) = struct('dir', [ 0 -1], 'weight', 2, 'angle', 90, 'loc', [-2 2], 'tile', [-2 4]);
s(4) = struct('dir', [ 1  0], 'weight', 2, 'angle', 90, 'loc', [ 0 0], 'tile', [-2 0]);

% define patterns for track growth & mutation
pair = ["LL"   "RR"   "LF"   "RF"   "FL"   "FR"   "FF"   "FF"  ];
quad = ["FLLF" "FRRF" "FLLR" "FRRL" "RLLF" "LRRF" "LRRL" "RLLR"];
mutate = ["FLF" "FRF" "FRL" "FLR" "LRL" "RLR" "RLF" "LRF"];

% initialization
done = 0;
overlap = 0;

% figure initialization
%want to specify axes and center, otherwise axes are randomized every time
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

if show_process == 1
    disp(track);
    full_plot = plot_track(track, track_length, s, show_loc, f1);
    delay(0.5);
end

% track growth/mutation loop
while track_length ~= target_length
    % track growth
    idx1 = randi([1, track_length]);    
    idx2 = mod(idx1,track_length)+1;
    two = [track(idx1) track(idx2)];
    for n = 1:7
        pattern = convertStringsToChars(pair(n));
        if strcmp(two, pattern)
            % backup the current track sequence & structure
            track_temp = track;
            s_temp = s;

            % replace original pair by quad sections
            quad_str = convertStringsToChars(quad(n)); 
            if idx1 == track_length
                track_segment = [quad_str(4) track(1:idx1-1) quad_str(1:3)];
                track = strcat(quad_str(4), track(1:idx1-1), quad_str(1:3));
            elseif idx1 == (track_length-1)
                rack_segment = [track(1:idx1-1) quad_str];
                track = strcat(track(1:idx1-1), quad_str);
            elseif idx1 == 1
                rack_segment = [track(1:idx1-1) quad_str];
                track = strcat(quad_str, track(idx1+2:track_length));
            else
                rack_segment = [track(1:idx1-1) quad_str  track(idx1+2:track_length)];
                track = strcat(track(1:idx1-1),quad_str, track(idx1+2:track_length));
            end
            track_length = track_length + 2;
            
            % update section structure
            for i = 1:track_length
                if track(i) == 'F'
                    s(i).dir = [1 0];
                    s(i).weight = 4;
                    s(i).angle = 0;
                    s(i).loc = [4 0];
                    s(i).tile = [2, 0];
                elseif track(i) == 'L'
                    s(i).dir = [0 1];
                    s(i).weight = 2;
                    s(i).angle = 90;
                    s(i).loc = [2 2];
                    s(i).tile = [2, 0];
                elseif track(i) == 'R'
                    s(i).dir = [0 -1];
                    s(i).weight = 2;
                    s(i).angle = -90;
                    s(i).loc = [2 -2];
                    s(i).tile = [2, 0];
                end 

                if i ~= 1
                    theta = s(i).angle; % to rotate 90 counterclockwise
                    R = [cosd(theta) -sind(theta); 
                         sind(theta)  cosd(theta)];
                    direction = transpose(s(i-1).dir); 
                    s(i).dir = transpose(R * direction); % new direction
                    disp_x = s(i-1).dir(1)*s(i).weight - 2*s(i-1).dir(2)*sind(theta);
                    diff = abs(s(i-1).dir(2)) - abs(s(i-1).dir(1));
                    disp_y = (s(i-1).dir(2)*s(i).weight - 2*s(i-1).dir(1)*sind(theta)) * diff;
                    displacement = [disp_x disp_y];
                    s(i).loc = s(i-1).loc + displacement; % new section location
                    s(i).tile = s(i).loc - s(i).dir*2;
                end
            end

            % check overlap
            overlap = checkOverlap(s);
            if overlap == 1
                % restore original track
                track = track_temp;
                s = s_temp;
                track_length = track_length - 2;
            else
                if show_process == 1
                    delete(full_plot);
                    disp(track);
                    full_plot = plot_track(track, track_length, s, show_loc, f1);
                    delay(0.2);
                end
                break; 
            end
        end
    end
    
    % track mutation 
    track_length = length(track);
    for m = 1:min(round(track_length*frequency/10), 10)    
        idx = randi([1, track_length]);
        idx1 = mod(idx-1,track_length)+1;
        idx2 = mod(idx,track_length)+1;
        idx3 = mod(idx+1,track_length)+1;
        three = [track(idx1) track(idx2) track(idx3)];
        for n = 1:8
            pattern = convertStringsToChars(mutate(n));
            if strcmp(three, pattern)
                track_temp = track;
                s_temp = s;
                new_pattern = convertStringsToChars(mutate(mod(n+3,8)+1));
                track(idx1) = new_pattern(1);
                track(idx2) = new_pattern(2);
                track(idx3) = new_pattern(3);
                for i = 1:track_length
                    if track(i) == 'F'
                        s(i).dir = [1 0];
                        s(i).weight = 4;
                        s(i).angle = 0;
                        s(i).loc = [4 0];
                        s(i).tile = [2, 0];
                    elseif track(i) == 'L'
                        s(i).dir = [0 1];
                        s(i).weight = 2;
                        s(i).angle = 90;
                        s(i).loc = [2 2];
                        s(i).tile = [2, 0];
                    elseif track(i) == 'R'
                        s(i).dir = [0 -1];
                        s(i).weight = 2;
                        s(i).angle = -90;
                        s(i).loc = [2 -2];
                        s(i).tile = [2, 0];
                    end 

                    if i ~= 1
                        theta = s(i).angle; % to rotate 90 counterclockwise
                        R = [cosd(theta) -sind(theta); 
                             sind(theta)  cosd(theta)];
                        direction = transpose(s(i-1).dir); 
                        s(i).dir = transpose(R * direction); % new direction
                        disp_x = s(i-1).dir(1)*s(i).weight - 2*s(i-1).dir(2)*sind(theta);
                        diff = abs(s(i-1).dir(2)) - abs(s(i-1).dir(1));
                        disp_y = (s(i-1).dir(2)*s(i).weight - 2*s(i-1).dir(1)*sind(theta)) * diff;
                        displacement = [disp_x disp_y];
                        s(i).loc = s(i-1).loc + displacement; % new section location
                        s(i).tile = s(i).loc - s(i).dir*2;
                    end
                end
                overlap = checkOverlap(s);
                if overlap == 1
                    track = track_temp;
                    s = s_temp;
                else
                    if show_process == 1
                        delete(full_plot);
                        disp(track);
                        full_plot = plot_track(track, track_length, s, show_loc, f1);
                    end
                    break; 
                end
            end
        end
    end
end

%track created
car_angle = pi/6;
%car_angle = 0;

cam_height = 2;
aov = pi/3;
cam_angle = pi/6;
car_pos = [0,0];
carX = car_pos(1);
carY = car_pos(2);
plot(carX, carY, '.', 'MarkerSize', 30, 'Color', 'g')
%identify track piece before, current, and future
%at the moment, all track works in relative manner- based on previous
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
gridlines(minX, maxX, minY, maxY, 4, 4, show_grid);

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

coord_arr = generate_joined_line(prev_tile2, prev_tile, current_tile, next_tile, next_tile2);
%inner
plot(coord_arr{1,1}, coord_arr{1,2}, 'Color', 'b', 'LineWidth', 4);
%outer
plot(coord_arr{1,3}, coord_arr{1,4}, 'Color', 'r', 'LineWidth', 4);

cam_bounds = find_cam_boundary(car_pos, car_angle, cam_height, aov, cam_angle);
plot(cam_bounds{1,1}(1), cam_bounds{1,1}(2), '.', 'MarkerSize', 20, 'Color', 'm');
plot(cam_bounds{1,2}(1), cam_bounds{1,2}(2), '.', 'MarkerSize', 20, 'Color', 'm');
plot(cam_bounds{1,3}(1), cam_bounds{1,3}(2), '.', 'MarkerSize', 20, 'Color', 'm');

linescanX = linspace(cam_bounds{1,1}(1), cam_bounds{1,2}(1));
linescanY = linspace(cam_bounds{1,1}(2), cam_bounds{1,2}(2));
plot(linescanX, linescanY, 'Color', 'm');
cam_inner = InterX([linescanX;linescanY],[coord_arr{1,1};coord_arr{1,2}]);
cam_outer = InterX([linescanX;linescanY],[coord_arr{1,3};coord_arr{1,4}]);
plot_intersect(cam_inner, cam_outer);
%right now, the program has the full camera line and the intersection
%points. We want the magnitude of the vector from each edge of the camera
%line to the closest intersection point. 
camera_lengths = calculate_distances(cam_inner, cam_outer, cam_bounds);
waveform = generateWaveform(camera_lengths);
noise = generateNoise(waveform, 33);
wav_x = linspace(1, 128, 128);
wav_y = waveform + noise;
waveformAxes = axes(f1);
waveformAxes.Units = 'pixels';
waveformAxes.Position = [400 200 150 200];
waveformAxes.Title.String = "Linescan Camera Waveform";
waveformAxes.XLabel.String = "Pixel";
waveformAxes.YLabel.String = "Grayscale Value";
waveformAxes.XLim = [1 128];
waveformAxes.YLim = [0 9];
plot(wav_x, wav_y);

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
        distanceR = norm(cam_outer - cam_bounds{1,2});
        center = norm(cam_bounds{1,1}-cam_outer);
    elseif isempty(cam_outer)
        distanceR = 0;
        distanceL = norm(cam_inner - cam_bounds{1,1});
        center = norm(cam_bounds{1,2}-cam_inner);
    else
        distanceL = norm(cam_inner - cam_bounds{1,1});
        distanceR = norm(cam_outer - cam_bounds{1,2});
        center = norm(cam_inner - cam_outer);
    end
    camera_lengths = [distanceL, center, distanceR];
end
% output track info
%output_track(track, track_length, s, show_loc, f1);
function waveform = generateWaveform(camera_lengths)
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
        waveform(left_px:center_px+left_px) = waveform(left_px:center_px+left_px) + 8.5;
        waveform(center_px+left_px+1:center_px+left_px+right_px) = waveform(center_px+left_px+1:center_px+left_px+right_px) + 3.5;
    end

end

function plot_intersect(cam_inner, cam_outer)
    if ~isempty(cam_inner)
        plot(cam_inner(1), cam_inner(2), 'x', 'MarkerSize', 20, 'Color', 'k');
    end
    if ~isempty(cam_outer)
        plot(cam_outer(1), cam_outer(2), 'x', 'MarkerSize', 20, 'Color', 'k');
    end    
end

function gridlines(minX, maxX, minY, maxY, intervalX, intervalY, show_grid)
    if show_grid
        for row = minY: intervalY : maxY
            line([minX, maxY], [row, row], 'Color', 'r');
        end
        for col = minX : intervalX : maxX
            line([col, col], [minY, maxX], 'Color', 'r');
        end
    end
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

function coord_array = generate_joined_line(prev_tile2, prev_tile, current_tile, next_tile, next_tile2)
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
    %diagonal line issue happens around reconstruction
    for i=1:length(reconstructed)
        joined_lineX1 = [joined_lineX1 reconstructed(i).xRangeS1];
        joined_lineX2 = [joined_lineX2 reconstructed(i).xRangeS2];
        joined_lineY1 = [joined_lineY1 reconstructed(i).yRangeS1];
        joined_lineY2 = [joined_lineY2 reconstructed(i).yRangeS2];
    end
    coord_array = {joined_lineX1, joined_lineY1, joined_lineX2, joined_lineY2};
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

% ckeckOverlap function
function overlap = checkOverlap(s)
    overlap = 0; 
    track_length = length(s);
    
    for i = 1:track_length
        for j = i+1:track_length
            if s(i).tile == s(j).tile
                overlap = 1;
                break;
            end
        end
    end
end

% output_track function
function output_track(track, track_length, s, show_loc, f1);
% plot track
plot_track(track, track_length, s, show_loc, f1);
% save figure as a png file
output_filename = strcat('track_', num2str(track_length), datestr(now, '_mm-dd-yyyy-HH-MM')); 
figure_filename = strcat(output_filename, '.png');
saveas(gcf,figure_filename); 

% save track as a txt file
track_filename = strcat(output_filename, '.txt');
outfile = fopen(track_filename,'w');
fprintf(outfile, track);

% save section structure as a csv file
struct_filename = strcat(output_filename, '.csv');
writetable(struct2table(s), struct_filename)

% save section structure as a mat file
struct_mat = strcat(output_filename, '.mat');
save(struct_mat, 's');
end

function delay(seconds)
    % function pause the program
    % seconds = delay time in seconds
    tic;
    while toc < seconds
    end
end

% plot track function
function full_plot = plot_track(track, track_length, s, show_loc, f1)
fp1 = [];
fp2 = [];
if show_loc == 1
    fc1 = [];
end
radius_inner = 1;  
radius_outer = 3;

title({['\fontsize{16}Track Auto-Generation (Length = ', num2str(track_length), ')'] ['\fontsize{11}Track Sequence = ', track]});

for i = 1:track_length  
    %Tile_Center = s(i).tile
    if track(i) == 'F'
        if show_loc == 1
            c1 = viscircles(s(i).loc, 0.1, 'Color', 'red');
        end
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
        if show_loc == 1
            c1 = viscircles(s(i).loc, 0.1, 'Color', 'red');
        end
        %viscircles(cir, 0.05, 'Color', 'black');
        axis equal;
    end
    fp1 = horzcat(fp1, p1);
    fp2 = horzcat(fp2, p2);
    if show_loc == 1
        fc1 = horzcat(fc1, c1);
    end
end
if show_loc == 1
    full_plot = [fp1, fp2, fc1];
else
    full_plot = [fp1, fp2];
end
drawnow;
end
