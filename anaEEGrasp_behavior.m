close all; clearvars; clc

%%
% import data from .csv files
sub_dir = uigetdir;
file_list = dir(fullfile(sub_dir, '*.csv'));

%%
data = cell(size(file_list));
resultant_mx = zeros(size(file_list));
angTilt = zeros(length(file_list), 1);
for i = 1:length(file_list)
    tmp = importfileEEGrasp( fullfile( file_list(i).folder, file_list(i).name ) );
    
    % All ATI Nano 25 are in N and N-m
    % rotate the ATI coordinate to handle coordinate
    tmp_th0 = table2array( tmp(:, {'fx0', 'fy0', 'fz0', 'mx0', 'my0', 'mz0'}) );
    tmp_th1 = table2array( tmp(:, {'fx1', 'fy1', 'fz1', 'mx1', 'my1', 'mz1'}) );
    tmp_v2 = table2array( tmp(:, {'fx2', 'fy2', 'fz2', 'mx2', 'my2', 'mz2'}) );
    tmp_v3 = table2array( tmp(:, {'fx3', 'fy3', 'fz3', 'mx3', 'my3', 'mz3'}) );
    % both thumb ATI rotate about z for 180 deg then about y for 180 deg to
    % the object coordinate
    RTH = rotz(90) * roty(180);
    tmp_th0 = [(RTH * tmp_th0(:, 1:3)')', (RTH * tmp_th0(:, 4:end)')'];
    tmp_th1 = [(RTH * tmp_th1(:, 1:3)')', (RTH * tmp_th1(:, 4:end)')'];
    % both virtual finger ATI rotate about z for 90 deg to the object
    % coordinate
    RV = rotz(-90);
    tmp_v2 = [(RV * tmp_v2(:, 1:3)')', (RV * tmp_v2(:, 4:end)')'];
    tmp_v3 = [(RV * tmp_v3(:, 1:3)')', (RV * tmp_v3(:, 4:end)')'];
    % replace with rotated data
    tmp{:, 3:26} = [tmp_th0, tmp_th1, tmp_v2, tmp_v3];
    % remove first and last 2 data to remove spline interpolation effects
    tmp = tmp(3:end-2, :);
    data{i} = tmp;
    
    %% get Table coordinate
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
    
    %% get time
    time = data{i}{:, {'time_ms'}};
    dt = 0.001 * (time(2) - time(1)); % in second
    % get trigger indices
    audio_trigger = data{i}{:, {'trigger'}};
    
    %% get PS Marker6: the center of the object
    % object coordinate before audio go cue and will keep until object
    % lift*
    ind_b4go = (audio_trigger == 1);
    obj_marker0_b4go = mean(data{i}{ind_b4go, {'x0', 'y0', 'z0'}}, 1);
    obj_marker1_b4go = mean(data{i}{ind_b4go, {'x1', 'y1', 'z1'}}, 1);
    obj_marker2_b4go = mean(data{i}{ind_b4go, {'x2', 'y2', 'z2'}}, 1);
    obj_marker4_b4go = mean(data{i}{ind_b4go, {'x4', 'y4', 'z4'}}, 1);
    obj_marker5_b4go = mean(data{i}{ind_b4go, {'x5', 'y5', 'z5'}}, 1);
    obj_marker6_b4go = mean(data{i}{ind_b4go, {'x6', 'y6', 'z6'}}, 1);
    obj_marker7_b4go = mean(data{i}{ind_b4go, {'x7', 'y7', 'z7'}}, 1);
    
    % coordinate fixed on the object for each frame!!!
    coord_obj_origin = obj_marker6_b4go;
    
    obj_vector1_b4go = obj_marker1_b4go - obj_marker4_b4go;
    obj_vector2_b4go = obj_marker0_b4go - obj_marker4_b4go;
    coord_obj_x_b4go = cross(obj_vector1_b4go, obj_vector2_b4go);
    
    obj_vector3_b4go = obj_marker5_b4go - obj_marker2_b4go;
    obj_vector4_b4go = obj_marker7_b4go - obj_marker2_b4go;
    % coord_obj_y_b4go = cross(obj_vector3_b4go, obj_vector4_b4go);
    
    % coord_obj_z_b4go = cross(coord_obj_x_b4go, coord_obj_y_b4go);
    
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    coord_obj_y_b4go = cross(obj_vector1_b4go, coord_obj_x_b4go);
    coord_obj_z_b4go = obj_vector1_b4go; % for now!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    %% find lift onset
    lift_marker0 = zeros(height(data{i}), 1);
    for j = 1:height(data{i})
        obj_marker0 = data{i}{j, {'x0', 'y0', 'z0'}};
        lift_marker0(j, 1) = sqrt(sum((obj_marker0 - obj_marker0_b4go).^2));
    end
    avg_lft = mean(lift_marker0(ind_b4go));
    std_lft = std(lift_marker0(ind_b4go));
%     tmp = find(abs(lift_marker0 - avg_lft) > 5 * std_lft);
    tmp = find(abs(lift_marker0 - avg_lft) > 15); % larger than 15 mm
    ind_lft_onset = tmp(1) - 1;
%     b = lift_marker0(ind_lft_onset, 1);
    
    %% compute Torque compensation at lif onset
%     ind_lft_onset = 1:height(data{i});
    obj_side = file_list(i).name(end-4);
    switch obj_side
        case 'R'
            finger_Th = data{i}{ind_lft_onset, {'fx0', 'fy0', 'fz0', 'mx0', 'my0', 'mz0'}};
            finger_V  = data{i}{ind_lft_onset, {'fx2', 'fy2', 'fz2', 'mx2', 'my2', 'mz2'}};
        case 'L'
            finger_Th = data{i}{ind_lft_onset, {'fx1', 'fy1', 'fz1', 'mx1', 'my1', 'mz1'}};
            finger_V  = data{i}{ind_lft_onset, {'fx3', 'fy3', 'fz3', 'mx3', 'my3', 'mz3'}};
        otherwise
            disp('something wrong with the file name!!')
    end
    
    % find th_cop_y and v_cop_y
    % since the ATI reference origin is on tool surface: dy = mx / fz
    th_cop_y = finger_Th(:, 4) ./ finger_Th(:, 3);
    v_cop_y = finger_V(:, 4) ./ finger_V(:, 3);
    
    % width b/w two nano 25 tool surface in meter
    width_obj = 2 * (21.6 + 3) * 0.001;
    
    diff_finger = finger_Th - finger_V;
    resultant_mx(i, 1) = diff_finger(:, 2) .* width_obj - (diff_finger(:, 3) .* (th_cop_y - v_cop_y));
    disp(resultant_mx(i, 1));
        %%
        b = abs(lift_marker0 - avg_lft);

% %     plotyy(1:length(audio_trigger), resultant_max, 1:length(audio_trigger), b)
% %     disp(i)
% %     disp(obj_side)
% %     pause



    %% compute object tilt
    coord_obj_z = zeros(height(data{i}), 3);
    for j = ind_lft_onset% 1:height(data{i})
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
        %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    end
    
%     %% filter
%     cutoff = 30; % in Hz
%     coord_obj_z_filtered = filtmat_class( dt, cutoff, coord_obj_z);
    

    for j = ind_lft_onset%1:height(data{i})
        % Angle b/w line (object z) and plane (table xz)
        angTilt(i, 1) = asind( abs(dot(coord_table_y, coord_obj_z(j, :))) / (sqrt(sum(coord_table_y.^2)) * sqrt(sum(coord_obj_z(j, :).^2))) );
    end

    %% compute finger tip coordinate without missing frames
    
end
[path, ~, ~] = fileparts(sub_dir);

save([path, 'temp_data.mat'], 'resultant_mx', 'file_list', 'angTilt')

%% plot
name = {file_list.name}';

tmp = char(name);
ind_IL = (tmp(:, 11) == 'I');
ind_TR = (tmp(:, 11) == 'T');
ind_PT = (tmp(:, 11) == 'P');

trial_id = 1:100;
mx = resultant_mx * 1000;
subplot 211
plot(trial_id(ind_IL), mx(ind_IL, 1), '-*r', trial_id(ind_TR), mx(ind_TR, 1), 'ob', trial_id(ind_PT), mx(ind_PT, 1), 'xk')
subplot 212
plot(trial_id(ind_IL), angTilt(ind_IL, 1), '-*r', trial_id(ind_TR), angTilt(ind_TR, 1), 'ob', trial_id(ind_PT), angTilt(ind_PT, 1), 'xk')
legend({'IL', 'TR', 'PT'})



