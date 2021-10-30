function center = locate_center(waveform)
% This function predicts the center of the HSAVC track by analyzing the
% virtual linescan camera waveform.
% ========================================================================
% INPUTS:
% - Virtual camera waveform array
% OUTPUTS:
% - Predicted track center pixel index (1-128) of the waveform. 

%COPIED AND PASTED FROM HSAVC EXERCISE BOOK.

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
    center = (second_edge+first_edge)/2;
else
    center = single(255);
end
end


