function [intersects, cam_line, waveform, recreated_lines] = camera_model(aov, cam_angle, cam_height, state, s)
    position = state(end, 1:3);
    angle = state(end,3);
    %Prepare: recreate the car's closest track pieces
    closest = locate_car(state, s);
    recreated_lines = recreate_track(closest);

    %1:
    %   make line of sight
    
    line_of_sight = calculate_sight(position, angle, cam_height, aov, cam_angle);
    %create line of sight line by connecting the start and the end.
    linescanX = linspace(line_of_sight{1,1}(1), line_of_sight{1,2}(1));
    linescanY = linspace(line_of_sight{1,1}(2), line_of_sight{1,2}(2));
    cam_line = {linescanX, linescanY}; %output line of sight for plotting
    
    %2: 
    %   find intersect between line of sight and nearest track
    
    cam_inner = find_intersect([linescanX;linescanY],[recreated_lines{1,1};recreated_lines{1,2}]); %left intersect
    cam_outer = find_intersect([linescanX;linescanY],[recreated_lines{2,1};recreated_lines{2,2}]); %right intersect
    intersects = {cam_inner, cam_outer}; %output for plotting
    
    %3:
    %   make waveform
    
    intersect_distances = calculate_distances(cam_inner, cam_outer, line_of_sight);
    waveform = make_waveform(intersect_distances);
    %create and add white noise to the linescan waveform
    noise = generate_noise(waveform, 33);
    waveform = noise+waveform;
end


