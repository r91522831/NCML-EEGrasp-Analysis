close all; clearvars; clc

%%
% import data from .csv files
sub_dir = uigetdir;
file_list = dir(fullfile(sub_dir, '*.csv'));

%% build variable names
% % % var_ATI = { {'fx0', 'fy0', 'fz0', 'mx0', 'my0', 'mz0'}; ...
% % %             {'fx1', 'fy1', 'fz1', 'mx1', 'my1', 'mz1'} };
var_ATI = cell(2, 1);
for i = 1:2
    tmp_id = i - 1;
    var_ATI{i} = { ['fx', num2str(tmp_id)], ['fy', num2str(tmp_id)], ['fz', num2str(tmp_id)],...
                   ['mx', num2str(tmp_id)], ['my', num2str(tmp_id)], ['mz', num2str(tmp_id)]};
end        
        
var_PS = cell(14, 1);
var_PS_cond = cell(14, 1);
for i = 1:14
    tmp_id = i - 1;
    var_PS{i} = {['x', num2str(tmp_id)], ['y', num2str(tmp_id)], ['z', num2str(tmp_id)]};
    var_PS_cond{i} = ['cond', num2str(tmp_id)];
end

%%
data = cell(size(file_list));
data_aligned2surface = cell(size(file_list));
for i = 1:length(file_list)
    data_raw = TShape_importfileEEGraspWithCond( fullfile( file_list(i).folder, file_list(i).name ) );
    %% replace first and last 2 data to remove spline interpolation effects
    data_raw(1:2, 3:end) = [data_raw(3, 3:end); data_raw(3, 3:end)];
    data_raw(end-1:end, 3:end) = [data_raw(end-2, 3:end); data_raw(end-2, 3:end)];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% All ATI Nano 25 are in N and N-m
    % 1 and 0 are for thumb or virtual finger
    % Thumb: ATI0 for fx0 in IL and PT; ATI0 for fx1 in TR
    % align all kinetic (ATI) coordinate to handle coordinate
    % and convert moments into N-mm
    tmp_ntrial = str2double(file_list(i).name(7:9));
    %======================================================================
    % only true for TShape sub-38
    if tmp_ntrial > 15 && rem(tmp_ntrial - 15, 5) == 1
        tmp_ATI0 = table2array( data_raw(:, var_ATI{1}) ); 
        tmp_ATI1 = table2array( data_raw(:, var_ATI{2}) ); % fx1, fy1, fz1
    else
        tmp_ATI0 = table2array( data_raw(:, var_ATI{2}) ); % fx1, fy1, fz1
        tmp_ATI1 = table2array( data_raw(:, var_ATI{1}) );
    end
    % use below lines for Tshape other subjects
    %{
    if tmp_ntrial > 15 && rem(tmp_ntrial - 15, 5) == 1
        tmp_ATI0 = table2array( data_raw(:, var_ATI{2}) );
        tmp_ATI1 = table2array( data_raw(:, var_ATI{1}) ); % fx0, fy0, fz0
    else
        tmp_ATI0 = table2array( data_raw(:, var_ATI{1}) ); % fx0, fy0, fz0
        tmp_ATI1 = table2array( data_raw(:, var_ATI{2}) );
    end
    %}
    %======================================================================
    
    % ---------------------------------------------------------------------
    % all angles are defined in the handle coordinate
    % ATI0 axes rotate about z for -90 deg then about y for -180 deg to the object coordinate
    RTH = rotz(90) * roty(180); % axes rotate angle = value rotate -angle
    tmp_ATI0_rotated = [(RTH * tmp_ATI0(:, 1:3)')', ((RTH * tmp_ATI0(:, 4:end)')' .* 1000)]; 
    % ATI1 axes rotate about z for 90 deg to the object coordinate
    RV = rotz(-90); % axes rotate angle = value rotate -angle
    tmp_ATI1_rotated = [(RV * tmp_ATI1(:, 1:3)')', ((RV * tmp_ATI1(:, 4:end)')' .* 1000)];
    
    % ---------------------------------------------------------------------
    % To translate all ATI coordinates along the handle z axis will change
    % only mx and my values
    % aligned to the respective handle center
    d_center2nano25 = 3 + 21.6; % in mm, Nano25 base to center (3 mm) plus Nano25 height (21.6 mm)
    
    new_z_L_ATI = d_center2nano25;
    tmp_ATI0_aligned2center = translateATI(tmp_ATI0_rotated, new_z_L_ATI);

    new_z_R_ATI = -d_center2nano25;
    tmp_ATI1_aligned2center = translateATI(tmp_ATI1_rotated, new_z_R_ATI);
    
    % all translated and rotated ATI to its own hanlde center
    data_ATI_aligned2handle_center = [tmp_ATI0_aligned2center, tmp_ATI1_aligned2center];
    
    % ---------------------------------------------------------------------
    % aligned to the respective finger surface
    d_ATI2surface = 3 + 2; % in mm, Nano25 surface to mounting (3 mm) plus mounting to cover (2 mm)
    
    new_z_L_ATI = -d_ATI2surface;
    tmp_ATI0_aligned2surface = translateATI(tmp_ATI0_rotated, new_z_L_ATI);
    
    new_z_R_ATI = d_ATI2surface;
    tmp_ATI1_aligned2surface = translateATI(tmp_ATI1_rotated, new_z_R_ATI);
        
    % all translated and rotated ATI to its own hanlde center
    data_ATI_aligned2surface = [tmp_ATI0_aligned2surface, tmp_ATI1_aligned2surface];
    
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
    coord_table_origin = table_marker0;
    coord_table_x = table_vector1 / norm(table_vector1);
    tmp = cross(table_vector1, table_vector2); % this is the gravity direction
    coord_table_y = tmp / norm(tmp);
    tmp = cross(coord_table_x, coord_table_y);
    coord_table_z = tmp / norm(tmp);
    
    %% convert all kinematics into table coordinate
    % first translate PS origin to table origin
    tmp = cell(1, length(var_PS));
    for j = 1:length(var_PS)
        tmp{1, j} = bsxfun(@minus, data_raw{:, var_PS{j}}, coord_table_origin);
    end

    %----------------------------------------------------------------------
    % rotate the coordinate system from PS to table.
    % [1, 0, 0] -> coord_table_x; [0, 1, 0] -> coord_table_y; % [0, 0, 1] -> coord_table_z
    R = [coord_table_x; coord_table_y; coord_table_z];
    for j = 1:length(var_PS)
        tmp{1, j} = (R * tmp{1, j}')';
    end

    data_PS_aligned = [tmp{:}];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% replace with rotated data
    tmp = data_raw;
    tmp{:, [var_ATI{:}]} = data_ATI_aligned2handle_center;
    tmp{:, [var_PS{:}]} = data_PS_aligned;
    data{i} = tmp;
    
    tmp = data_raw;
    tmp{:, [var_ATI{:}]} = data_ATI_aligned2surface;
    tmp(:, [var_PS{:}, var_PS_cond{:}]) = [];
    data_aligned2surface{i} = tmp;
    
end

[path, ~, ~] = fileparts(sub_dir);
subID = file_list(1).name(1:4);
save(fullfile(path, [subID, '_aligned_data.mat']), 'data', 'data_aligned2surface', 'file_list', 'sub_dir', 'var_ATI', 'var_PS', 'var_PS_cond');
clear tmp* i j
disp("Data alignment completed!")