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
    table_marker0 = mean(data{i}{:, {'x11', 'y11', 'z11'}}, 1);
    table_marker1 = mean(data{i}{:, {'x12', 'y12', 'z12'}}, 1);
    table_marker2 = mean(data{i}{:, {'x13', 'y13', 'z13'}}, 1);
    
    table_vector1 = table_marker1 - table_marker0;
    table_vector2 = table_marker2 - table_marker0;
    
    coord_table_origin = table_marker0;
    coord_table_x = table_vector1;
    coord_table_z = cross(table_vector1, table_vector2);
    coord_table_y = cross(coord_table_x, coord_table_z);
    
    %%
    %
    audio_trigger = data{i}{:, {'trigger'}};
    % get trigger indices
    a = find(gradient(audio_trigger) > 0);
    
    %%
    % get PS Marker6: the center of the object
    coord_object_origin = data{i}{:, {'x6', 'y6', 'z6'}};
    
    %%
    % object coordinate before audio go cue 
    ind_b4go = (audio_trigger == 1);
    obj_marker0 = mean(data{i}{ind_b4go, {'x0', 'y0', 'z0'}}, 1);
    obj_marker1 = mean(data{i}{ind_b4go, {'x1', 'y1', 'z1'}}, 1);
    obj_marker4 = mean(data{i}{ind_b4go, {'x4', 'y4', 'z4'}}, 1);
    
    obj_vector1 = obj_marker4 - obj_marker1;
    obj_vector2 = obj_marker0 - obj_marker1;
    coord_obj_x = obj_vector1;
    coord_obj_y = cross(obj_vector1, obj_vector2);
    coord_obj_z = cross(coord_obj_x, coord_obj_y);
    
    
    
    
    
    
end


