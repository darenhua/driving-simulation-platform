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
full_plot = [];

% ask user to input the target track length
target_length = input('Input Target Track Length (Even Number, Greater than 4): ');

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



% output track info
output_track(track, track_length, s, show_loc, f1);

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
