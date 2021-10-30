function waveform = make_waveform(intersect_distances)
% This function creates the 128 pixel wide waveform outputted by the
% virtual linescan camera.
% ========================================================================
% INPUTS:
% - Distances between intersects and line of sight
% OUTPUTS:
% - 128 element array with constant high or low values, serves as the
%   foundation of the waveform before applying noise.

waveform = zeros(1,128); %initialize
cam_length = 0;
%re-calculate the total length of the camera line of sight by summing the
%intersect distances if they exist.
for i = 1:length(intersect_distances)
    if ~isempty(intersect_distances(i))
        cam_length = cam_length + intersect_distances(i);
    end
end
%ratio used to convert length of the line of sight in inches to 128 pixel
%wide waveform. These are the amount of pixels for each distance.
left_px = abs(int16((intersect_distances(1)*128)/cam_length));
right_px = abs(int16((intersect_distances(3)*128)/cam_length));
center_px = abs(int16((intersect_distances(2)*128)/cam_length));

%if no intersects are found, the camera outputs all low waveform
if intersect_distances(1) == 0 && intersect_distances(2) == 0 && intersect_distances(3) == 0
    waveform = waveform + 3.5;
%if the left intersect does not exist, only the right side of the waveform
%is visible.( --|__ )
elseif (intersect_distances(1) == 0) || (left_px == 0)
    waveform(1:center_px) = waveform(1:center_px) + 6.5;
    waveform(center_px+1:128) = waveform(center_px+1:128) + 3.5;
%if the right intersect does not exist, only the left side of the waveform
%is visible. ( __|-- )
elseif intersect_distances(3) == 0 || (right_px == 0)
    waveform(1:left_px-1) = waveform(1:left_px-1) + 3.5;
    waveform(left_px:128) = waveform(left_px:128) + 6.5;
%if both intersects exist, the waveform appears like: ( _|-|_ )
else
    waveform(1:left_px-1) = waveform(1:left_px-1) + 3.5;
    waveform(left_px:center_px+left_px) = waveform(left_px:center_px+left_px) + 6.5;
    waveform(center_px+left_px+1:128) = waveform(center_px+left_px+1:128) + 3.5;
end
end
