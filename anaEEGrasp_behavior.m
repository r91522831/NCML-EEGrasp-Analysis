close all; clearvars; clc

%% load aligned data
[filename, pathname, ~] = uigetfile;
load(fullfile(pathname, filename));

% % % coord_table_origin = [0, 0, 0];
% % % coord_table_x = [1, 0, 0];
% % % coord_table_y = [0, 1, 0];
% % % coord_table_z = [0, 0, 1];

% % % n_PSonObj = 8;
% % % n_PSonFgr = 3;
% % %
% % % kinetic_channels = {'fx', 'fy', 'fz', 'mx', 'my', 'mz'};

%%
angTilt2R = cell(size(file_list));
ind_lft_onset = nan(length(file_list), 6);
info_onset_time = nan(length(file_list), 1);

obj_height = cell(length(file_list), 3);
obj_weight = nan(length(file_list), 1);
peak_roll = table(nan(size(file_list)), nan(size(file_list)), 'VariableNames', {'peakRoll', 'index'});
peak_mx = table(nan(size(file_list)), nan(size(file_list)), 'VariableNames', {'peakMx', 'index'});
mx_onset = nan(length(file_list), 1);

%%
input = data;
for i = 1:length(file_list)
    tmp_audio = info_time_trigger{i, 2};
    dt = 0.001 * (info_time_trigger{i, 1}(2, 1) - info_time_trigger{i, 1}(1, 1)); % in second
    %% Get object coordinate before reaching initiation
    % object coordinate before audio go cue and should keep the same until object lift onset*
    ind_b4go = (tmp_audio == 1);
    tmp = cellfun(@table2array, coord_obj{i}(ind_b4go, 1),'UniformOutput',false);
    coord_obj_b4go = array2table(mean(cat(3, tmp{:}), 3), 'RowNames', {'x', 'y', 'z'}, 'VariableNames', {'x_axis', 'y_axis', 'z_axis', 'origin', 'RCenter', 'LCenter'});
    
    % the object center is a virtual center at the middle of right and left
    % handle centers
    tmp_obj_Rcenter_b4go = mean(coord_obj_b4go{:, {'RCenter'}}, 2)';
    tmp_obj_Lcenter_b4go = mean(coord_obj_b4go{:, {'LCenter'}}, 2)';
    tmp_obj_center_b4go = mean([tmp_obj_Rcenter_b4go; tmp_obj_Lcenter_b4go], 1);
    
    % Angle b/w object axis z and table z, positive represents lower right obj wrt subjects
    tmp_angTilt2R_b4go = atan2d(coord_obj_b4go{'z', 'z_axis'}, coord_obj_b4go{'y', 'z_axis'}); % range -180 < x < 180 degrees
    % Angle b/w object z axis and table xz plane
% % %     tmp_angTilt_b4go = asind( abs(dot(coord_table_y, coord_obj_b4go{:, 'z_axis'})) / (sqrt(sum(coord_table_y.^2)) * sqrt(sum(coord_obj_b4go{:, 'z_axis'} .^2))) );
    
    %% frame by frame
    obj_Rcenter_v = zeros(height(input{i}), 3);
    obj_Lcenter_v = zeros(height(input{i}), 3);
    obj_center_v = zeros(height(input{i}), 3);
    obj_Rcenter_d = zeros(height(input{i}), 1);
    obj_Lcenter_d = zeros(height(input{i}), 1);
    obj_center_d = zeros(height(input{i}), 1);
    for time_id = 1:height(input{i})
        % compute object height
        if isempty(coord_obj{i}{time_id, 1})
            continue;
        end
        
        obj_Rcenter_v(time_id, :) = (coord_obj{i}{time_id, 1}{:, {'RCenter'}}' - tmp_obj_Rcenter_b4go);
        obj_Lcenter_v(time_id, :) = (coord_obj{i}{time_id, 1}{:, {'LCenter'}}' - tmp_obj_Rcenter_b4go);
        obj_center_v(time_id, :) = mean([obj_Rcenter_v(time_id, :); obj_Lcenter_v(time_id, :)], 1);
        
        obj_Rcenter_d(time_id, :) = sqrt(sum((obj_Rcenter_v(time_id, :)).^2));
        obj_Lcenter_d(time_id, :) = sqrt(sum((obj_Lcenter_v(time_id, :)).^2));
        obj_center_d(time_id, :) = mean([obj_Rcenter_d(time_id, :), obj_Lcenter_d(time_id, :)]);
        
        % Angle b/w object axis z and table z
        angTilt2R{i}(time_id, 1) = atan2d(coord_obj{i}{time_id, 1}{'z', 'z_axis'}, coord_obj{i}{time_id, 1}{'y', 'z_axis'}) - tmp_angTilt2R_b4go;
        % Angle b/w object z axis and table xz plane
% % %         angTilt{i}(time_id, 1) = abs(tmp_angTilt_b4go - asind( abs(dot(coord_table_y, coord_obj{i}{time_id, 1}{:, 'z_axis'})) / (sqrt(sum(coord_table_y.^2)) * sqrt(sum(coord_obj{i}{time_id, 1}{:, 'z_axis'} .^2))) ));
    end
    
    %% define lift onset
    obj_height{i, 1} = min([obj_Rcenter_v(:, 2), obj_Lcenter_v(:, 2)], [], 2);
    obj_height{i, 2} = min([obj_Rcenter_d, obj_Lcenter_d], [], 2);
    obj_height{i, 3} = obj_center_d;
    
    % based on kinematics
    ind_hold = find(tmp_audio == 3);
    for j = ind_hold(end, 1):-1:1 % 10 mm
        if obj_height{i, 1}(j, 1) < 10
            ind_lft_onset(i, 1) = j;
            break;
        end
    end
    for j = ind_hold(end, 1):-1:1 %  5 mm
        if obj_height{i, 1}(j, 1) < 5
            ind_lft_onset(i, 2) = j;
            break;
        end
    end
    for j = ind_hold(end, 1):-1:1 %  3 mm
        if obj_height{i, 1}(j, 1) < 3
            ind_lft_onset(i, 3) = j;
            break;
        end
    end
    for j = ind_hold(end, 1):-1:1 %  2 mm
        if obj_height{i, 1}(j, 1) < 2
            ind_lft_onset(i, 4) = j;
            break;
        end
    end
    for j = ind_hold(end, 1):-1:1 %  1 mm
        if obj_height{i, 1}(j, 1) < 1
            ind_lft_onset(i, 5) = j;
            break;
        end
    end
    
    % based on kinetics
    % compute object weight
    tmp_stable_window = 100; % in ms
    obj_weight(i, 1) = mean( sqrt(sum(resultantF{i, :}{(ind_hold(end, 1) - tmp_stable_window):(ind_hold(end, 1) + tmp_stable_window), {'fx', 'fy', 'fz'}}.^2, 2)) );
    
    % onset when total force larger than the object weight
    for j = ind_lft_onset(i, 1):-1:1 % start from 10 mm backward
        tmp_rf = sqrt(sum(resultantF{i, :}{j, {'fx', 'fy', 'fz'}}.^2, 2));
        if tmp_rf < obj_weight(i, 1) % load force equal to object weight
            ind_lft_onset(i, 6) = j;
            break;
        end
    end
    

    %% choosen onset
    tmp_ind_onset = 3; % 3 mm    
    info_onset_time(i, 1) = 0.001 * info_time_trigger{i, 1}(ind_lft_onset(i, tmp_ind_onset)); % in seconds
    
    %% find peak roll after lift onset
    roll_win = 250; % in ms
    ind_roll_win = floor(roll_win ./ (dt * 1000));
    [~, tmp_ind] = max( abs(angTilt2R{i, 1}(ind_lft_onset(i, tmp_ind_onset):(ind_lft_onset(i, tmp_ind_onset) + ind_roll_win), 1)) );
    tmp_roll = angTilt2R{i, 1}(ind_lft_onset(i, tmp_ind_onset) + tmp_ind, 1);
    peak_roll{i, {'peakRoll', 'index'}} = [tmp_roll, ind_lft_onset(i, tmp_ind_onset) + tmp_ind];
    
    %% find peak mx around lift onset
    % this might need to be rewrite
    pmx_win = 50; % in ms
    ind_pmx_win = floor(pmx_win ./ (dt * 1000));
    [~, tmp_ind] = max( abs(resultantF{i, 1}{(ind_lft_onset(i, tmp_ind_onset) - ind_pmx_win):ind_lft_onset(i, tmp_ind_onset), 'mx'}) );
    tmp_ind = (ind_lft_onset(i, tmp_ind_onset) - ind_pmx_win) + tmp_ind;
    peak_mx{i, {'peakMx', 'index'}} = [resultantF{i, 1}{tmp_ind, 'mx'}, tmp_ind];

    
    %% mx at lift onset
    mx_onset(i) = resultantF{i, 1}{ind_lft_onset(i, tmp_ind_onset), 'mx'};
    
    %% compute finger tip coordinate without missing frames
    
    
    %% compute the finger average coordinate within the holding phase
    
    
    
    
    
    %%
    disp(i);
end

%%
ind_lft_onset = array2table(ind_lft_onset, 'VariableNames', {'h10_mm', 'h5_mm', 'h3_mm', 'h2_mm', 'h1_mm', 'fy_obj_w'});

save(fullfile(pathname, [filename(1:4), '_temp_result.mat']), 'resultantF', 'finger_Th', 'finger_V', 'angTilt2R', 'ind_lft_onset', 'info_onset_time', 'file_list', 'obj_height', 'obj_weight', 'peak_roll', 'peak_mx', 'info_time_trigger');
save(fullfile(pathname, [filename(1:4), '_info_onset_time.mat']), 'info_onset_time');