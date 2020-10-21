%VEHICLE MODEL TO CAMERA
%input position of car, angle
%output waveform the linescan camera would output at that location
%ease of life todo: implement Preset instead of Postset for addlistener
%callback, so that don't need to let go of slider


clc, clf, clear;

f2 = figure('position', [360 500 400 400]);
defaultAx = get(gca);
defaultAx.YAxis.Visible = 'off'; % remove y-axis
defaultAx.XAxis.Visible = 'off'; % remove x-axis

hold on;
posSlider(f2);
edgeLText = uicontrol('style','text', 'position', [50 295 50 30]);
edgeRText = uicontrol('style','text', 'position', [300 295 50 30]);
edgeLText.String = "Left Edge";
edgeRText.String = "Right Edge";
function posSlider(f2)
    %axes initialization
    waveformAxes = axes(f2);
    waveformAxes.Units = 'pixels';
    waveformAxes.Position = [100 50 200 200];
    waveformAxes.NextPlot = 'add';
    waveformAxes.Title.String = "Linescan Camera Waveform";
    waveformAxes.XLabel.String = "Pixel";
    waveformAxes.YLabel.String = "Grayscale Value";
    waveformAxes.XLim = [1 128];
    waveformAxes.YLim = [0 7];
    
    x = linspace(1, 128, 128);
    %ideally y should be a function of x, could change later
    edge = [0,24];
    lineInitial = generateWaveform(edge);
    noiseInitial = generateNoise(lineInitial);
    y = lineInitial + noiseInitial;
    lineH = plot(x,y);
    sliderVal = uicontrol('style','text', 'position', [140 340 120 15]);
    sliderPos = uicontrol('style','slider','position', [100 300 200 20], 'min', 0, 'max', 24);
    %car will have a position of 0-24 inches on the 24 inch track
    addlistener(sliderPos, 'Value', 'PostSet', @callbackFunction);
    movegui(f2, 'center')
    function callbackFunction(source, event)
        num = get(event.AffectedObject, 'Value');
        %a value 0-24
        edge = pos2edge(num);
        line = generateWaveform(edge);
        noise = generateNoise(line);
        out = line + noise;
        lineH.YData  = out;
        sliderVal.String = strcat(num2str(num), " inches");
    end
end

function noise = generateNoise(waveform)
    %the lower signal to noise ratio, the more noisy function
    regsnr = 33;
    waveLength = length(waveform);
    signalPower = (sum(abs(waveform).^2))/waveLength;
    noisePower =  signalPower/(10^(regsnr/10));
    %we want random amount of additon/subtraction: range (-1,1)
    randomArr = (2 .* rand(1,waveLength)) - 1;
    noise = sqrt(noisePower) * randomArr;
end

function edge = pos2edge(position)
    %turn the 0-24 inch pos input into an array with [x y] where x is  
    %inches from left edge of track and y is inches from right edge
    %constant: how many inches is black edge??? TO UPDATE TODO
    edgeL = position - 0;
    edgeR = 24 - position;
    edge = [edgeL edgeR];
    
    %potential issue: does not consider the thickness of edge. TO CHECK 
end

function shiftedWaveform = generateWaveform(edge)
    %input: d from track edge in array [x y]
    %where x is d from left edge, y is d from right
    % 0<x=y<24, 0 and 24 is cam on edge
    %future: need to account for width of car in relation
    waveform = zeros(1,128);
    edgeL = edge(1);
    edgeR = edge(2);
    %edge is in physical inch value, waveform is pixel
    %no matter values of edgeL, edgeR, track waveform same.
    edgeStrength = int16(20); %~in px, could use rand to randomize
    trackCenter = int16(80); %~in px, could use rand to randomize
    
    center = int16(64);
    centerLStart = center - (trackCenter/2);
    centerRStart = center + (trackCenter/2);
    edgeLStart = centerLStart - edgeStrength;
    edgeRStart = centerRStart + edgeStrength;
    
    waveform(edgeLStart:centerLStart-1) = waveform(edgeLStart:centerLStart-1) + 2.5;
    waveform(centerLStart:centerRStart) = waveform(centerLStart:centerRStart) + 5.5;
    waveform(centerRStart+1:edgeRStart) = waveform(centerRStart+1:edgeRStart) + 2.5;
    %^ IS DEFAULT P the only thing that edgeL, edgeR changes is the shifting of graph
    %TEMPORARY SHIFT METHOD: Subtraction. Obviously wrong. 
    shiftN = round(edgeR - edgeL);
    shiftedWaveform = zeros(size(waveform));
    if shiftN > 0
       shiftedWaveform(1+shiftN:end) = waveform(1:end-shiftN);
    elseif shiftN < 0 
       shiftN = abs(shiftN);
       shiftedWaveform(1:end-shiftN) = waveform(1+shiftN:end);
    elseif shiftN == 0
        shiftedWaveform = waveform;
    end
    
    shiftedWaveform = shiftedWaveform + 1;
end
