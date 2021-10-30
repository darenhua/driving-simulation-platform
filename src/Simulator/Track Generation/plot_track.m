function plot_track(axes, s)
% This function re-plots a previously recreated track
% ========================================================================
% INPUTS:
% - Axes to plot onto.
% - Track sequence structure.
% OUTPUTS:
% - Track plotted on provided Axes.

radius_inner = 1;
radius_outer = 3;
for i = 1:length(s)
    % For each track section, check if it is Forwards or Curve, then plot.
    
    if s(i).weight == 4 % Forwards section:
        if s(i).dir_1 ~= 0 % Horizontal pointing:
            p1 = line(axes, [s(i).loc_1,s(i).loc_1-s(i).dir_1*4],[s(i).loc_2+1,s(i).loc_2+1], 'Color', 'black','LineWidth',4);
            p2 = line(axes, [s(i).loc_1,s(i).loc_1-s(i).dir_1*4],[s(i).loc_2-1,s(i).loc_2-1], 'Color', 'black','LineWidth',4);
        else % Vertical pointing:
            p1 = line(axes, [s(i).loc_1+1,s(i).loc_1+1],[s(i).loc_2,s(i).loc_2-s(i).dir_2*4], 'Color', 'black','LineWidth',4);
            p2 = line(axes, [s(i).loc_1-1,s(i).loc_1-1],[s(i).loc_2,s(i).loc_2-s(i).dir_2*4], 'Color', 'black','LineWidth',4);            
        end
    else % Curve section (Right or Left):
        angle_index = atan2(s(i).dir_2, s(i).dir_1)/(pi/2)*2 + 3;
        th = linspace((angle_index+sin(s(i).angle*pi/180))*pi/4, (angle_index+2+sin(s(i).angle*pi/180))*pi/4, 100);
        x = radius_inner*cos(th) + s(i).loc_1 - 2*s(i).dir_2*sin(s(i).angle*pi/180);
        y = radius_inner*sin(th) + s(i).loc_2 + 2*s(i).dir_1*sin(s(i).angle*pi/180);
        p1 = plot(axes, x,y,'Color', 'black','LineWidth',4); 
        cir_x = s(i).loc_1 - 2*s(i).dir_2*sin(s(i).angle*pi/180);
        cir_y = s(i).loc_2 + 2*s(i).dir_1*sin(s(i).angle*pi/180);
        x = radius_outer*cos(th) + cir_x;
        y = radius_outer*sin(th) + cir_y;
        p2 = plot(axes, x,y,'Color', 'black','LineWidth',4);
        axis(axes, "equal");        
    end
end
end

