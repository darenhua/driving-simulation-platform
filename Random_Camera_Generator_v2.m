clc, clf, clear;

edge = [4, 4];
f2 = figure('position', [360 500 400 400]);
hold on;
plot(x,y);
title("Linescan Camera Waveform");
xlabel("Pixel");
ylabel("Greyscale Value");
xlim([1, 128]);
ylim([0,7]);


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
    edge = [4,4];
    
    y = generateWaveform(edge);
    lineH = plot(x,y);
    sliderVal = uicontrol('style','text', 'position', [170 340 40 15]);
    sliderH = uicontrol('style','slider','position', [100 280 200 20], 'min', 0, 'max', 24);
    %car will have a position of 0-24 inches on the 24 inch track
    addlistener(sliderH, 'Value', 'PreSet', @callbackFunction);
    movegui(f2, 'center')
    function callbackFunction(source, event)
        num = get(event.AffectedObject, 'Value');
        %a value 0-24
        lineH.YData  = sin(num * x);
        sliderVal.String = num2str(num);
    end
end

function edge = pos2edge(position)
    %turn the 0-24 inch pos input into an array with [x y] where x is  
    %inches from left edge of track and y is inches from right edge
    edgeStrength = 2;
    %constant: how many inches is black edge??? TO UPDATE TODO
    
end

function waveform = generateWaveform(edge)
    %input: d from track edge in array [x y]
    %where x is d from left edge, y is d from right
    % 0<x=y<24, 0 and 24 is cam on edge
    %future: need to account for width of car in relation
    waveform = zeros(1,128);
    edgeL = edge(1);
    edgeR = edge(2);
    %edge is in physical inch value, waveform is pixel

end
