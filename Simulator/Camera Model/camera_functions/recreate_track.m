function nearest_line = recreate_track(tile_arr)
% This function generates a 300 element array of the car's closest 5 track
% sections.
% ========================================================================
% INPUTS:
% - Array of the car's closest 5 track sections.
% OUTPUTS:
% - Cell Array with XY coordinates of inner and outer closest track lines.

reconstructed = struct();
for i=2:length(tile_arr)-1
    if tile_arr(i).weight == 4
        %F
        reconstructed(i-1).type = 'F';
        if tile_arr(i).dir_1 ~= 0
            if tile_arr(i).dir_1 > 0 
                x1_s1 = tile_arr(i-1).loc_1;
                x2_s1 = tile_arr(i).loc_1;
                y1_s1 = tile_arr(i).loc_2+1;
                y1_s2 = tile_arr(i).loc_2-1;
                x_s1 = linspace(x1_s1,x2_s1); % Creates 100 x-coordinates
                x_length = length(x_s1);
                y_s1 = zeros(1,x_length) + y1_s1; % Creates 100 y-coordinates
                x_s2 = linspace(x1_s1,x2_s1); 
                y_s2 = zeros(1,x_length) + y1_s2;
            else
                x1_s1 = tile_arr(i-1).loc_1;
                x2_s1 = tile_arr(i).loc_1;
                y1_s1 = tile_arr(i).loc_2+1;
                y1_s2 = tile_arr(i).loc_2-1;
                x_s1 = linspace(x1_s1,x2_s1);
                x_length = length(x_s1);
                y_s1 = zeros(1,x_length) + y1_s1;
                x_s2 = linspace(x1_s1,x2_s1);
                y_s2 = zeros(1,x_length) + y1_s2;
            end
        else
            if tile_arr(i).dir_2 > 0
                x1_s1 = tile_arr(i).loc_1-1;
                x1_s2 = tile_arr(i).loc_1+1;
                y1_s1 = tile_arr(i-1).loc_2;
                y2_s1 = tile_arr(i).loc_2;
                y_s1 = linspace(y1_s1,y2_s1);
                y_length = length(y_s1);
                x_s1 = zeros(1,y_length) + x1_s1;
                y_s2 = linspace(y1_s1,y2_s1);
                x_s2 = zeros(1,y_length) + x1_s2;
            else
                x1_s1 = tile_arr(i).loc_1+1;
                x1_s2 = tile_arr(i).loc_1-1;
                y1_s1 = tile_arr(i-1).loc_2;
                y2_s1 = tile_arr(i).loc_2;
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
        if tile_arr(i).dir_1 ~= 0
            centerX = tile_arr(i).loc_1;
            if tile_arr(i).dir_1 > 0
                centerY = tile_arr(i).loc_2 + 2;
                x1_s1 = tile_arr(i-1).loc_1 + 1;
                x2_s1 = centerX;
                x1_s2 = tile_arr(i-1).loc_1 - 1;
                x2_s2 = centerX;
                x_s1 = linspace(x1_s1, x2_s1);
                x_s2 = linspace(x1_s2, x2_s2);
                th_s1 = -acos(x_s1-centerX);
                th_s2 = -acos((x_s2-centerX)/3);
                y_s1 = sin(th_s1) + centerY;
                y_s2 = 3*sin(th_s2) + centerY;

            else
                centerY = tile_arr(i).loc_2 - 2;
                x1_s1 = centerX;
                x2_s1 = tile_arr(i-1).loc_1 - 1;
                x1_s2 = centerX;
                x2_s2 = tile_arr(i-1).loc_1 + 1;
                x_s1 = linspace(x2_s1, x1_s1);
                x_s2 = linspace(x2_s2, x1_s2);
                th_s1 = acos(x_s1-centerX);
                th_s2 = acos((x_s2-centerX)/3);
                y_s1 = sin(th_s1) + centerY;
                y_s2 = 3*sin(th_s2) + centerY;

            end
        else
            centerY = tile_arr(i).loc_2;
            if tile_arr(i).dir_2 > 0
                centerX = tile_arr(i).loc_1 - 2;
                x1_s1 = centerX;
                x2_s1 = tile_arr(i).loc_1 - 1;
                x1_s2 = centerX;
                x2_s2 = tile_arr(i).loc_1 + 1;
                x_s1 = linspace(x1_s1, x2_s1);
                x_s2 = linspace(x1_s2, x2_s2);
                th_s1 = -acos(x_s1-centerX);
                th_s2 = -acos((x_s2-centerX)/3);
                y_s1 = sin(th_s1) + centerY;
                y_s2 = 3*sin(th_s2) + centerY;

            else
                centerX = tile_arr(i).loc_1 + 2;
                x1_s1 = tile_arr(i).loc_1 + 1;
                x2_s1 = centerX;
                x1_s2 = tile_arr(i).loc_1 - 1;
                x2_s2 = centerX;
                x_s1 = linspace(x2_s1, x1_s1);
                x_s2 = linspace(x2_s2, x1_s2);
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
        reconstructed(i-1).type = 'R';
        if tile_arr(i).dir_1 ~= 0
            centerX = tile_arr(i).loc_1;
            if tile_arr(i).dir_1 > 0
                centerY = tile_arr(i).loc_2 - 2;
                %duplicates, should probably just include outside of
                %loop
                x1_s1 = tile_arr(i-1).loc_1 + 1;
                x2_s1 = centerX;
                x_s1 = linspace(x1_s1, x2_s1);
                x1_s2 = tile_arr(i-1).loc_1 - 1;
                x2_s2 = centerX;
                x_s2 = linspace(x1_s2, x2_s2);
                th_s1 = acos(x_s1-centerX);
                th_s2 = acos((x_s2-centerX)/3);
                y_s1 = sin(th_s1) + centerY;
                y_s2 = 3*sin(th_s2) + centerY;
            else
                centerY = tile_arr(i).loc_2 + 2;
                x1_s1 = centerX;
                x2_s1 = tile_arr(i-1).loc_1 - 1;
                x_s1 = linspace(x2_s1, x1_s1);
                x1_s2 = centerX;
                x2_s2 = tile_arr(i-1).loc_1 + 1;
                x_s2 = linspace(x2_s2, x1_s2);
                th_s1 = -acos(x_s1-centerX);
                th_s2 = -acos((x_s2-centerX)/3);
                y_s1 = sin(th_s1) + centerY;
                y_s2 = 3*sin(th_s2) + centerY;
            end
        else
            centerY = tile_arr(i).loc_2;
            if tile_arr(i).dir_2 > 0
                centerX = tile_arr(i).loc_1 + 2;
                x1_s1 = tile_arr(i).loc_1 + 1;
                x2_s1 = centerX;
                x_s1 = linspace(x2_s1, x1_s1);
                x1_s2 = tile_arr(i).loc_1 - 1;
                x2_s2 = centerX;
                x_s2 = linspace(x2_s2, x1_s2);
                th_s1 = -acos(x_s1-centerX);
                th_s2 = -acos((x_s2-centerX)/3);
                y_s1 = sin(th_s1) + centerY;
                y_s2 = 3*sin(th_s2) + centerY;
            else
                centerX = tile_arr(i).loc_1 - 2;
                x1_s1 = centerX;
                x2_s1 = tile_arr(i).loc_1 - 1;
                x_s1 = linspace(x1_s1, x2_s1);
                x1_s2 = centerX;
                x2_s2 = tile_arr(i).loc_1 + 1;
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
    if i+1<=length(reconstructed)
        lastX = int16(reconstructed(i).xRangeS1(end));
        firstX = int16(reconstructed(i+1).xRangeS1(1));
        lastY = int16(reconstructed(i).yRangeS1(end));
        firstY = int16(reconstructed(i+1).yRangeS1(1));
        if (lastX ~= firstX) || (lastY ~= firstY)
            if (((reconstructed(i).type == 'L') && (reconstructed(i+1).type == 'R')) ||...
                   ((reconstructed(i).type == 'R') && (reconstructed(i+1).type == 'L')))
                  if (i+1 == 3) &&...
                          ((int16(reconstructed(i-1).xRangeS1(end)) == int16(reconstructed(i).xRangeS1(1))) ||...
                          (int16(reconstructed(i-1).xRangeS1(end)) == int16(reconstructed(i).xRangeS1(1))))
                      temp = reconstructed(i+1).xRangeS1;
                      temp2 = reconstructed(i+1).xRangeS2;
                      reconstructed(i+1).xRangeS1 = temp2;
                      reconstructed(i+1).xRangeS2 = temp;
                      temp = reconstructed(i+1).yRangeS1;
                      temp2 = reconstructed(i+1).yRangeS2;
                      reconstructed(i+1).yRangeS1 = temp2;
                      reconstructed(i+1).yRangeS2 = temp;
                  elseif (i+1 == 2) &&...
                          ((int16(reconstructed(i+1).xRangeS1(end)) == int16(reconstructed(i+2).xRangeS1(1))) ||...
                          (int16(reconstructed(i+1).xRangeS1(end)) == int16(reconstructed(i+2).xRangeS1(1))))
                      temp = reconstructed(i).xRangeS1;
                      temp2 = reconstructed(i).xRangeS2;
                      reconstructed(i).xRangeS1 = temp2;
                      reconstructed(i).xRangeS2 = temp;
                      temp = reconstructed(i).yRangeS1;
                      temp2 = reconstructed(i).yRangeS2;
                      reconstructed(i).yRangeS1 = temp2;
                      reconstructed(i).yRangeS2 = temp;
                  else 
                      temp = reconstructed(i+1).xRangeS1;
                      temp2 = reconstructed(i+1).xRangeS2;
                      reconstructed(i+1).xRangeS1 = temp2;
                      reconstructed(i+1).xRangeS2 = temp;
                      temp = reconstructed(i+1).yRangeS1;
                      temp2 = reconstructed(i+1).yRangeS2;
                      reconstructed(i+1).yRangeS1 = temp2;
                      reconstructed(i+1).yRangeS2 = temp;
                  end
            end
        end
    end
end

for i=1:length(reconstructed)
    joined_lineX1 = [joined_lineX1 reconstructed(i).xRangeS1];
    joined_lineX2 = [joined_lineX2 reconstructed(i).xRangeS2];
    joined_lineY1 = [joined_lineY1 reconstructed(i).yRangeS1];
    joined_lineY2 = [joined_lineY2 reconstructed(i).yRangeS2];
end
nearest_line = {joined_lineX1, joined_lineY1; joined_lineX2, joined_lineY2};
end
