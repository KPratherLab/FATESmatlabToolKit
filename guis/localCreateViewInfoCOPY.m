function viewinfo = localCreateViewInfoCOPY(hAxes)
%this is a copy of a subfunction found in resetplotview (internal matlab
%function) that is used to create a structure to hold info on an axis
axes_properties = {'DataAspectRatio',...
    'CameraViewAngle',...
    'PlotBoxAspectRatio',...
    'CameraPosition',...
    'CameraTarget',...
    'CameraUpVector',...
    'XLim',...
    'YLim',...
    'ZLim'};

% Save the value of each axes property and its mode
for i = 1:numel(axes_properties)
    current_prop = axes_properties{i};
    current_mode = [axes_properties{i} 'Mode'];
    viewinfo.(current_mode) = get(hAxes,current_mode);
    % Only get properties in manual since getting them in auto can trigger
    % an auto-calc
    if strcmp(get(hAxes,current_mode),'manual')
        viewinfo.(current_prop) = get(hAxes,current_prop);
    end
end

[az, el] = view(hAxes);
viewinfo.View = [az, el];
end