close all; clearvars; clc

%% load aligned data
[filename, pathname, ~] = uigetfile;
load(fullfile(pathname, filename));

coord_table_origin = [0, 0, 0];
coord_table_x = [1, 0, 0];
coord_table_y = [0, 1, 0];
coord_table_z = [0, 0, 1];

% % % n_PSonObj = 8;
% % % n_PSonFgr = 3;
% % %
% % % kinetic_channels = {'fx', 'fy', 'fz', 'mx', 'my', 'mz'};


%% lowpass filter all data
dt = diff(data{1}{1:2, 1}) * 0.001; % in second
% % % cutoff = 30; % in Hz
data_filtered = data;
% % % for i = 1:length(file_list)
% % %     data_filtered{i}{:, 3:end} = filtmat_class( dt, cutoff, data{i}{:, 3:end} );
% % % end

%%
angTilt = cell(size(file_list));
ind_lft_onset = zeros(length(file_list), 6);

obj_height = cell(size(file_list));
obj_weight = zeros(length(file_list), 1);
peak_roll = table(zeros(size(file_list)), zeros(size(file_list)), 'VariableNames', {'peakRoll', 'index'});
peak_mx = table(zeros(size(file_list)), zeros(size(file_list)), 'VariableNames', {'peakMx', 'index'});

fingr_th_height = cell(size(file_list));
fingr_in_height = cell(size(file_list));
fingr_mi_height = cell(size(file_list));

%%
input = data_filtered; % or data
for i = 1%:length(file_list)
    tmp_audio = info_time_trigger{i, 2};
    %% Get object coordinate before reaching initiation
    % object coordinate before audio go cue and should keep the same until object lift onset*
    ind_b4go = (tmp_audio == 1);
    tmp = cellfun(@table2array, coord_obj{i}(ind_b4go, 1),'UniformOutput',false);
    coord_obj_b4go = array2table(mean(cat(3, tmp{:}), 3), 'RowNames', {'x', 'y', 'z'}, 'VariableNames', {'x_axis', 'y_axis', 'z_axis', 'origin'});
    
    % Angle b/w line (object z) and plane (table xz)
    angTilt_b4go = asind( abs(dot(coord_table_y, coord_obj_b4go{:, 'z_axis'})) / (sqrt(sum(coord_table_y.^2)) * sqrt(sum(coord_obj_b4go{:, 'z_axis'} .^2))) );
    
    %% frame by frame
    obj_center = zeros(height(input{i}), 1);
    for time_id = 1:height(input{i})
        %% compute object height
        obj_center(time_id, 1) = sqrt(sum((markers{7, 1} - markers_b4go{7, 1}).^2));
        
        
        %% compute object tilt
        markers = cell(n_PSonObj, 1);
        for m = 1:n_PSonObj
            markers{m, 1} = input{i}{time_id, var_PS{m}};
        end
        coord_obj{i}{time_id, 1} = coordOnObj(markers, obj_side);
        
        % Angle b/w line (object z) and plane (table xz)
        angTilt{i}(time_id, 1) = abs(angTilt_b4go - asind( abs(dot(coord_table_y, coord_obj{i}{time_id, 1}{'z_axis', :})) / (sqrt(sum(coord_table_y.^2)) * sqrt(sum(coord_obj{i}{time_id, 1}{'z_axis', :} .^2))) ));
        
        
    end
    
    
    
    
    
    
    
    %{
    %% define lift onset
    avg_lft = mean(lift_marker0(ind_b4go));
    std_lft = std(lift_marker0(ind_b4go));
    %     tmp = find(abs(lift_marker0 - avg_lft) > 5 * std_lft);
    obj_height{i} = abs(lift_marker0 - avg_lft);
    
    % based on kinematics
    ind_hold = find(tmp_audio == 3);
    for j = ind_hold(end, 1):-1:1
        if obj_height{i, :}(j, 1) < 10 % 10 mm
            ind_lft_onset(i, 1) = j;
            break;
        end
    end
    for j = ind_hold(end, 1):-1:1
        if obj_height{i, :}(j, 1) < 5 % 5 mm
            ind_lft_onset(i, 2) = j;
            break;
        end
    end
    for j = ind_hold(end, 1):-1:1
        if obj_height{i, :}(j, 1) < 4 % 4 mm
            ind_lft_onset(i, 3) = j;
            break;
        end
    end
    for j = ind_hold(end, 1):-1:1
        if obj_height{i, :}(j, 1) < 3 % 3 mm
            ind_lft_onset(i, 4) = j;
            break;
        end
    end
    for j = ind_hold(end, 1):-1:1
        if obj_height{i, :}(j, 1) < 1 % 1 mm
            ind_lft_onset(i, 5) = j;
            break;
        end
    end
    
    % based on kinetics
    % compute object weight
    ind_hold = find(tmp_audio == 3);
    tmp_stable_window = 50; % in ms
    obj_weight(i, 1) = mean( sqrt(sum(resultantF{i, :}{(ind_hold(end, 1) - tmp_stable_window):ind_hold(end, 1), {'fx', 'fy', 'fz'}}.^2, 2)) );
    
    for j = ind_lft_onset(i, 1):-1:1 % start from 10 mm backward
        tmp_rf = sqrt(sum(resultantF{i, :}{j, {'fx', 'fy', 'fz'}}.^2, 2));
        if tmp_rf < obj_weight(i, 1) % load force equal to object weight
            ind_lft_onset(i, 6) = j;
            break;
        end
    end
    
    
    %% choosen onset
    tmp_ind_onset = 2; % 5mm
    
    %% find peak roll after lift onset
    roll_win = 250; % in ms
    ind_roll_win = floor(roll_win ./ (dt * 1000));
    [tmp_roll, tmp_ind] = max(angTilt{i, 1}(ind_lft_onset(i, tmp_ind_onset):(ind_lft_onset(i, tmp_ind_onset) + ind_roll_win), 1));
    peak_roll{i, {'peakRoll', 'index'}} = [tmp_roll, ind_lft_onset(i, tmp_ind_onset) + tmp_ind];
    
    
    
    %% find peak mx around lift onset
    %{
    % this might need to be rewrite
    pmx_win = 50; % in ms
    ind_pmx_win = floor(pmx_win ./ (dt * 1000));
    [~, tmp_ind] = max( abs(resultantF{i, 1}{(ind_lft_onset(i, tmp_ind_onset) - ind_pmx_win):ind_lft_onset(i, tmp_ind_onset), 'mx'}) );
    tmp_ind = (ind_lft_onset(i, tmp_ind_onset) - ind_pmx_win) + tmp_ind;
    peak_mx{i, {'peakMx', 'index'}} = [resultantF{i, 1}{tmp_ind, 'mx'}, tmp_ind];
    %}
    %%
    
    %}
    
    %% compute finger tip coordinate without missing frames
    %%
    disp(i);
end

% save(fullfile(pathname, [filename(1:4), '_temp_result.mat']), 'resultantF', 'finger_Th*', 'finger_V*', 'angTilt', 'ind_lft_onset', 'file_list', 'pathname', 'obj_height', 'obj_weight', 'peak_roll', 'peak_mx', 'info_time_trigger');