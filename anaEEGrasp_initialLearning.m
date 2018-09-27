close all; clearvars; clc

%% load aligned data
[filename, pathname, ~] = uigetfile;
filelist = dir(fullfile(pathname, [filename(1:4), '*.mat']));
for i = 1:length(filelist)
    load(fullfile(pathname, filelist(i).name));
end

%% lowpass filter all data
dt = diff(data{1}{1:2, 1}) * 0.001; % in second
% % % cutoff = 30; % in Hz
data_filtered = data;
% % % for i = 1:length(file_list)
% % %     data_filtered{i}{:, 3:end} = filtmat_class( dt, cutoff, data{i}{:, 3:end} );
% % % end


%%
obj = cell(length(file_list), 1);
fgr_all = cell(length(file_list), 1);

input = data_filtered; % or data
for i = 1:length(file_list)
    % get only trials from the initial learing session
    if file_list(i).name(:, 11) ~= 'I'
        continue;
    end
    
    % object coordinate before audio go cue and should keep the same until object lift onset*
    ind_b4go = (info_time_trigger{i, 2} == 1);
    tmp = cellfun(@table2array, coord_obj{i}(ind_b4go, 1),'UniformOutput',false);
    coord_obj_b4go = array2table(mean(cat(3, tmp{:}), 3), 'RowNames', {'x', 'y', 'z'}, 'VariableNames', {'x_axis', 'y_axis', 'z_axis', 'origin', 'RCenter', 'LCenter'});
    
    r_matrix_b4go = coord_obj_b4go{:, {'x_axis', 'y_axis', 'z_axis'}};
    eul_b4go = rotm2eul(r_matrix_b4go, 'ZYX');
    %% frame by frame
    eul = zeros(height(input{i}), 6);
    for time_id = 1:height(input{i})
        r_matrix = coord_obj{1, 1}{time_id, 1}{:, {'x_axis', 'y_axis', 'z_axis'}};
        tmp_eul = rotm2eul(r_matrix, 'ZYX') - eul_b4go; % ZYX -> pitch, yaw, roll
        eul(time_id, :) = [coord_obj{i, 1}{time_id, 1}{:, {'origin'}}', tmp_eul];
    end
    obj{i, 1} = array2table(eul, 'VariableNames', {'x', 'y', 'z', 'pitch', 'yaw', 'roll'});
    
    % individual fingers
    fgr_all{i, 1}.Th = [finger_Th{i, 1}, fgr_on_obj{i}.Th{:}]; % fx, fy, fz, mx, my, mz, COPx, COPy, COPz, x, y, z, cond
    fgr_all{i, 1}.Vi = finger_V{i, 1}; % fx, fy, fz, mx, my, mz, COPx, COPy, COPz
    fgr_all{i, 1}.In = fgr_on_obj{i}.In{:}; % x, y, z, cond
    fgr_all{i, 1}.Mi = fgr_on_obj{i}.Mi{:}; % x, y, z, cond
    
    
% % %     finger_Vi
% % %     fgr_on_obj % x, y, z, cond
    
    
    disp(i)
end