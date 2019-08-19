close all; clearvars; clc

%% load aligned data
[filename, pathname, ~] = uigetfile;
file_suffix = {'_aligned_data', '_fingerInfo', '_info_onset_time', '_temp_result'};
for i = 1:length(file_suffix)
    load(fullfile(pathname, [filename(1:4), file_suffix{i}, '.mat']));
end

%%
obj = cell(length(file_list), 1);
fgr_all = cell(length(file_list), 1);

input = data;
for i = 1:length(file_list)
    % object coordinate before audio go cue and should keep the same until object lift onset*
    ind_b4go = (info_time_trigger{i, 2} == 1);
    tmp = cellfun(@table2array, coord_obj{i}(ind_b4go, 1),'UniformOutput',false);
    coord_obj_b4go = array2table(mean(cat(3, tmp{:}), 3), 'RowNames', {'x', 'y', 'z'}, 'VariableNames', {'x_axis', 'y_axis', 'z_axis', 'origin', 'RCenter', 'LCenter'});
    
    r_matrix_b4go = coord_obj_b4go{:, {'x_axis', 'y_axis', 'z_axis'}};
    eul_b4go = rotm2eul(r_matrix_b4go, 'ZYX');
    %% frame by frame
    eul = nan(height(input{i}), 6);
    for time_id = 1:height(input{i})
        if isempty(coord_obj{i, 1}{time_id, 1})
            continue;
        end
        r_matrix = coord_obj{i, 1}{time_id, 1}{:, {'x_axis', 'y_axis', 'z_axis'}};
        tmp_eul = rotm2eul(r_matrix, 'ZYX') - eul_b4go; % ZYX -> pitch, yaw, roll; angles are in radian
        eul(time_id, :) = [coord_obj{i, 1}{time_id, 1}{:, {'origin'}}', tmp_eul];
    end
    obj{i, 1} = array2table(eul, 'VariableNames', {'x', 'y', 'z', 'pitch', 'yaw', 'roll'});
    
    % individual fingers
    tmp_resultant = resultantF{i, 1};
    tmp_th = [finger_Th{i, 1}, fgr_on_obj{i}.Th{:}]; % fx, fy, fz, mx, my, mz, COPx, COPy, COPz, x, y, z, cond
    tmp_vi = finger_V{i, 1}; % fx, fy, fz, mx, my, mz, COPx, COPy, COPz
    tmp_in = fgr_on_obj{i}.In{:}; % x, y, z, cond
    tmp_mi = fgr_on_obj{i}.Mi{:}; % x, y, z, cond
    fgr_all{i, 1} = cell2table({tmp_th, tmp_in, tmp_mi, tmp_vi, tmp_resultant}, 'VariableNames', {'Th', 'In', 'Mi', 'Vi', 'Re'});
    
    disp(i)
end

save(fullfile(pathname, [filename(1:4), '_for_plot.mat']), 'obj', 'fgr_all', 'info_time_trigger', 'ind_lft_onset', 'file_list', '-v7');