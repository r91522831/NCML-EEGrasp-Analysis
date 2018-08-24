close all; clearvars; clc

%%
% import data from .csv files
sub_dir = uigetdir;
file_list = dir(fullfile(sub_dir, '*.csv'));

data = cell(size(file_list));
for i = 1%:length(file_list)
    tmp = importfileEEGrasp(fullfile(file_list(i).folder,file_list(i).name ));
    tmp = tmp(3:end-2, :); % remove first and last 2 data to remove spline interpolation effects
    data{i} = tmp;
    
    %%
    % get Table coordinate
    % Wu, G., Cavanagh, P.R., 1995. Recommendations for standardization in
    % the reporting of kinematic data. Journal of Biomechanics 28 (10), 
    % 1257–1260.
    table_marker0 = mean(data{i}{:, {'x11', 'y11', 'z11'}}, 1);
    table_marker1 = mean(data{i}{:, {'x12', 'y12', 'z12'}}, 1);
    table_marker2 = mean(data{i}{:, {'x13', 'y13', 'z13'}}, 1);
    
    table_vector1 = table_marker1 - table_marker0;
    table_vector2 = table_marker2 - table_marker0;
    
    coord_table_origin = table_marker0;
    coord_table_x = table_vector1;
    coord_table_y = cross(table_vector1, table_vector2);
    coord_table_z = cross(coord_table_x, coord_table_y);
    
    %%
    % get trigger indices
    audio_trigger = data{i}{:, {'trigger'}};
    
    %%
    % get PS Marker6: the center of the object
    % object coordinate before audio go cue and will keep until object
    % lift*
    ind_b4go = (audio_trigger == 1);
    obj_marker0 = mean(data{i}{ind_b4go, {'x0', 'y0', 'z0'}}, 1);
    obj_marker1 = mean(data{i}{ind_b4go, {'x1', 'y1', 'z1'}}, 1);
    obj_marker2 = mean(data{i}{ind_b4go, {'x2', 'y2', 'z2'}}, 1);
    obj_marker4 = mean(data{i}{ind_b4go, {'x4', 'y4', 'z4'}}, 1);
    obj_marker5 = mean(data{i}{ind_b4go, {'x5', 'y5', 'z5'}}, 1);
    obj_marker6 = mean(data{i}{ind_b4go, {'x6', 'y6', 'z6'}}, 1);
    obj_marker7 = mean(data{i}{ind_b4go, {'x7', 'y7', 'z7'}}, 1);
    
    % coordinate fixed on the object for each frame!!!
    coord_obj_origin = obj_marker6;
    
    obj_vector1 = obj_marker1 - obj_marker4;
    obj_vector2 = obj_marker0 - obj_marker4;
    coord_obj_x = cross(obj_vector1, obj_vector2);
    
    obj_vector3 = obj_marker5 - obj_marker2;
    obj_vector4 = obj_marker7 - obj_marker2;
    % coord_obj_y = cross(obj_vector3, obj_vector4);
    
    % coord_obj_z = cross(coord_obj_x, coord_obj_y);
    
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    coord_obj_y = cross(obj_vector1, coord_obj_x);
    coord_obj_z = obj_vector1; % for now!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    %%
    % compute object tilt
    angTilt = zeros(height(data{i}), 1);
    coord_obj_z = zeros(height(data{i}), 3);
    for j = 1:height(data{i})
        obj_marker0 = data{i}{j, {'x0', 'y0', 'z0'}};
        obj_marker1 = data{i}{j, {'x1', 'y1', 'z1'}};
        obj_marker2 = data{i}{j, {'x2', 'y2', 'z2'}};
        obj_marker4 = data{i}{j, {'x4', 'y4', 'z4'}};
        obj_marker5 = data{i}{j, {'x5', 'y5', 'z5'}};
        obj_marker7 = data{i}{j, {'x7', 'y7', 'z7'}};
        
        obj_vector1 = obj_marker1 - obj_marker4;
        obj_vector2 = obj_marker0 - obj_marker4;
        coord_obj_x = cross(obj_vector1, obj_vector2);
        
        obj_vector3 = obj_marker5 - obj_marker2;
        obj_vector4 = obj_marker7 - obj_marker2;
        % coord_obj_y = cross(obj_vector3, obj_vector4);
        % coord_obj_z = cross(coord_obj_x, coord_obj_y);
        
        %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        coord_obj_y = cross(obj_vector1, coord_obj_x);
        coord_obj_z(j, :) = obj_vector1; % for now!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
% %         coord_obj_z = filtmat_class( dt, 30, data);
        %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        % Angle b/w line (object z) and plane (table xz)
        
        angTilt(j, 1) = asind( abs(dot(coord_table_y, coord_obj_z(j, :))) / (sqrt(sum(coord_table_y.^2)) * sqrt(sum(coord_obj_z(j, :).^2))) );
    end
    
    %%
    % compute finger tip coordinate without missing frames
    
    
    
    
end


