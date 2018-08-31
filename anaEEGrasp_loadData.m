close all; clearvars; clc

%%
% import data from .csv files
sub_dir = uigetdir;
file_list = dir(fullfile(sub_dir, '*.csv'));

%% build variable names
var_ATI = { {'fx0', 'fy0', 'fz0', 'mx0', 'my0', 'mz0'}; ...
            {'fx1', 'fy1', 'fz1', 'mx1', 'my1', 'mz1'}; ...
            {'fx2', 'fy2', 'fz2', 'mx2', 'my2', 'mz2'}; ...
            {'fx3', 'fy3', 'fz3', 'mx3', 'my3', 'mz3'} };
var_PS = cell(14, 1);
var_PS_cond = cell(14, 1);
for i = 1:14
    tmp_id = i - 1;
    var_PS{i} = {['x', num2str(tmp_id)], ['y', num2str(tmp_id)], ['z', num2str(tmp_id)]};
    var_PS_cond{i} = ['cond', num2str(tmp_id)];
end

%%
data = cell(size(file_list));
coord_table_origin = cell(size(file_list));
for i = 1%:length(file_list)
    data_raw = importfileEEGraspWithCond( fullfile( file_list(i).folder, file_list(i).name ) );
    %% replace first and last 2 data to remove spline interpolation effects
    data_raw(1:2, 3:end) = [data_raw(3, 3:end); data_raw(3, 3:end)];
    data_raw(end-1:end, 3:end) = [data_raw(end-2, 3:end); data_raw(end-2, 3:end)];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% All ATI Nano 25 are in N and N-m
    % 1 and 3 are on Left handle; 0 and 2 are on the Right handle
    % 1 and 0 are for thumb; 2 and 3 are for virtual finger
    % align all kinetic (ATI) coordinate to handle coordinate
    % and convert moments into N-mm
    tmp_th0 = table2array( data_raw(:, var_ATI{1}) );
    tmp_th1 = table2array( data_raw(:, var_ATI{2}) );
    tmp_v2 = table2array( data_raw(:, var_ATI{3}) );
    tmp_v3 = table2array( data_raw(:, var_ATI{4}) );
    % all the angles are defined in the handle coordinate
    % both thumb ATI axes rotate about z for -90 deg then about y for -180 deg to
    % the object coordinate
    RTH = rotz(90) * roty(180); % axes rotate angle = value rotate -angle
    tmp_th0_rotated = [(RTH * tmp_th0(:, 1:3)')', ((RTH * tmp_th0(:, 4:end)')' .* 1000)];
    tmp_th1_rotated = [(RTH * tmp_th1(:, 1:3)')', ((RTH * tmp_th1(:, 4:end)')' .* 1000)];
    % both virtual finger ATI axes rotate about z for 90 deg to the object
    % coordinate
    RV = rotz(-90); % axes rotate angle = value rotate -angle
    tmp_v2_rotated = [(RV * tmp_v2(:, 1:3)')', ((RV * tmp_v2(:, 4:end)')' .* 1000)];
    tmp_v3_rotated = [(RV * tmp_v3(:, 1:3)')', ((RV * tmp_v3(:, 4:end)')' .* 1000)];
    
    % To translate all ATI coordinates along the handle z axis
    % will change only mx and my values
    d_center2nano25 = 3 + 21.6; % in mm, Nano25 base to center (3 mm) plus Nano25 height (21.6 mm)
    % translate all rotated thumb ATI to their respective handle center
    % [1, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0; 0, 0, 1, 0, 0, 0; 0, d, 0, 1, 0, 0; d, 0, 0, 0, 1, 0; 0, 0, 0, 0, 0, 0]
    TrTH = eye(6, 6);
    TrTH(4, 2) = -d_center2nano25; % in the handle coordinate
    TrTH(5, 1) = -d_center2nano25;
    tmp_th0_rotated_translated = (TrTH * tmp_th0_rotated')';
    tmp_th1_rotated_translated = (TrTH * tmp_th1_rotated')';
    
    TrV = eye(6, 6);
    TrV(4, 2) = d_center2nano25; % in the handle coordinate
    TrV(5, 1) = d_center2nano25;
    tmp_v2_rotated_translated = (TrV * tmp_v2_rotated')';
    tmp_v3_rotated_translated = (TrV * tmp_v3_rotated')';
    
    % all translated and rotated ATI to its own hanlde center
    data_ATI_aligned = [tmp_th0_rotated_translated, tmp_th1_rotated_translated, tmp_v2_rotated_translated, tmp_v3_rotated_translated];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% get Table coordinate and convert all kinematic data from PhaseSpace (PS)
    % Wu, G., Cavanagh, P.R., 1995. Recommendations for standardization in
    % the reporting of kinematic data. Journal of Biomechanics 28 (10), pp. 1257-1261.
    table_marker0 = mean(data_raw{:, var_PS{12}}, 1); % PS ID 11
    table_marker1 = mean(data_raw{:, var_PS{13}}, 1); % PS ID 12
    table_marker2 = mean(data_raw{:, var_PS{14}}, 1); % PS ID 13
    
    table_vector1 = table_marker1 - table_marker0;
    table_vector2 = table_marker2 - table_marker0;
    
    % get unit vectors for the coordinate axes
    coord_table_origin{i, 1} = table_marker0;
    coord_table_x = table_vector1 / norm(table_vector1);
    tmp = cross(table_vector1, table_vector2); % this is the gravity direction
    coord_table_y = tmp / norm(tmp);
    tmp = cross(coord_table_x, coord_table_y);
    coord_table_z = tmp / norm(tmp);
    
    %% convert all kinematics into table coordinate
    % first translate PS origin to table origin
    tmp = cell(1, length(var_PS));
    for j = 1:length(var_PS)
        tmp{1, j} = bsxfun(@minus, data_raw{:, var_PS{j}}, coord_table_origin{i, 1});
    end
    %----------------------------------------------------------------------
    % In the future:
    % rotate the coordinate system from PS to table.
    % [1, 0, 0] -> coord_table_x; [0, 1, 0] -> coord_table_y; % [0, 0, 1] -> coord_table_z
    %----------------------------------------------------------------------
    % rotate the coordinate from PS to table
%    R = [coord_table_x; coord_table_y; coord_table_z];

    
    
    % use PS coordinate with table origin for now
    data_PS_aligned = [tmp{:}];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% replace with rotated data
    tmp = data_raw;
    tmp{:, [var_ATI{:}]} = data_ATI_aligned;
    tmp{:, [var_PS{:}]} = data_PS_aligned;
    data{i} = tmp;
end
[path, subID, ~] = fileparts(sub_dir);
save(fullfile(path, [subID, '_aligned_data.mat']), 'data', 'file_list', 'sub_dir', 'var_ATI', 'var_PS');
clear tmp* i j
