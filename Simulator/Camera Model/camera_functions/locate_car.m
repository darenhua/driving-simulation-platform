function closest = locate_car(state, s)
% This function locates the car's closest 5 track sections.
% ========================================================================
% INPUTS:
% - Track sequence structure.
% - Car state
% OUTPUTS:
% - Array of the car's closest 5 track sections.

% Search through track pieces and check if the car's position is in it.
% Return the track piece the car is inside. (ct_index)
for index=1:length(s)
    tile_range = [s(index).tile_1-2 s(index).tile_1+2 ; s(index).tile_2-2 s(index).tile_2+2];
    if (tile_range(1,1) <= state(1)) && (state(1) <= tile_range(1,2)) && (tile_range(2,1) <= state(2)) && (state(2) <= tile_range(2,2)) 
        ct_index = index;
    end
end

% The closest 5 track sections are used in case the camera line of sight
% intersects with track pieces outside of the car's current tile.
current_tile = s(ct_index);

%if the current array index is the end of the array, the next index is the
%first one.
if ct_index+1 <= length(s)
    next_tile = s(ct_index+1);
else
    next_tile = s(1);
end

%if the current array index is the start of the array, the previous index 
%is the last one.
if ct_index-1 > 0
    prev_tile = s(ct_index-1);
else
    prev_tile = s(end);
end

if ct_index+2 <= length(s)
    next_tile2 = s(ct_index+2);
elseif ct_index+2 == length(s)+1
    next_tile2 = s(1);
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

%group into array
closest = [prev_tile2, prev_tile, current_tile, next_tile, next_tile2];

end                
