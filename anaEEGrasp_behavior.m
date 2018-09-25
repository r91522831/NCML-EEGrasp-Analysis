close all; clearvars; clc

%% load aligned data
[filename, pathname, ~] = uigetfile;
load(fullfile(pathname, filename));

coord_table_origin = [0, 0, 0];
coord_table_x = [1, 0, 0];
coord_table_y = [0, 1, 0];
coord_table_z = [0, 0, 1];

n_PSonObj = 8;
n_PSonFgr = 3;

kinetic_channels = {'fx', 'fy', 'fz', 'mx', 'my', 'mz'};

%% lowpass filter all data
dt = diff(data{1}{1:2, 1}) * 0.001; % in second
% % % cutoff = 30; % in Hz
data_filtered = data;
% % % for i = 1:length(file_list)
% % %     data_filtered{i}{:, 3:end} = filtmat_class( dt, cutoff, data{i}{:, 3:end} );
% % % end

%% 
resultantF = cell(size(file_list));
finger_Th = cell(size(file_list));
finger_V = cell(size(file_list));
finger_Th_surface = cell(size(file_list));
finger_V_surface = cell(size(file_list));
angTilt = cell(size(file_list));
coord_obj = cell(size(file_list));
fgr_on_obj = cell(size(file_list));
ind_lft_onset = zeros(length(file_list), 6);
info_time_trigger = cell(length(file_list), 2);
missing_info = cell(length(file_list), 1);
obj_height = cell(size(file_list));
obj_weight = zeros(length(file_list), 1);
peak_roll = table(zeros(size(file_list)), zeros(size(file_list)), 'VariableNames', {'peakRoll', 'index'});
peak_mx = table(zeros(size(file_list)), zeros(size(file_list)), 'VariableNames', {'peakMx', 'index'});

fingr_th_height = cell(size(file_list));
fingr_in_height = cell(size(file_list));
fingr_mi_height = cell(size(file_list));


%--------------------------------------------------------------------------
% assign filtered or raw data to analyze
input = data_filtered; % or data
input_surface = data_aligned2surface;

for i = 1%:length(file_list)
    % get which side of the handle on the object is grasped
    obj_side = file_list(i).name(end-4);
    % get time stamp
    tmp_tstamp = input{i}{:, {'time_ms'}};
    info_time_trigger{i, 1} = tmp_tstamp;
    % get trigger indices
    tmp_audio = input{i}{:, {'trigger'}};
    info_time_trigger{i, 2} = tmp_audio;
    
    %% Check the PS sensor in each frame
    markers = cell(n_PSonObj, 1);
    for m = 1:n_PSonObj
        markers{m, 1} = input{i}{:, [var_PS{m}, var_PS_cond{m}]};
    end
    markers_fgr = cell(n_PSonFgr, 1);
    for m = 1:n_PSonFgr
        markers_fgr{m, 1} = input{i}{:, [var_PS{m + 8}, var_PS_cond{m + 8}]};
    end
    
    tmp_missing = 1;
    for time_id = 1:height(input{i})
        % for markers on the object
        tmp_markers = cell(n_PSonObj, 1);
        tmp_cond = zeros(n_PSonObj, 1);
        for m = 1:n_PSonObj
            tmp_markers{m, 1} = markers{m, 1}(time_id, 1:3);
            tmp_cond(m, 1) = markers{m, 1}(time_id, 4);
        end
        tmp_yz_markers = [2, 4, 5]; % default markers
        tmp_xz_markers = [3, 6, 8];
        tmp_origin_markers = 1:8;
        
        % for markers on three finger tips
        tmp_markers_fgr = cell(n_PSonFgr, 1);
        tmp_cond_fgr = zeros(n_PSonFgr, 1);
        for m = 1:n_PSonFgr
            tmp_markers_fgr{m, 1} = markers_fgr{m, 1}(time_id, 1:3);
            tmp_cond_fgr(m, 1) = markers_fgr{m, 1}(time_id, 4);
        end
        
        % if there is any marker missing in this frame
        if any(tmp_cond < 0) || any(tmp_cond_fgr < 0)
            tmp_missing_markers = find(tmp_cond < 0);
            tmp_missing_markers_fgr = find(tmp_cond_fgr < 0);
            % check if there is at least 3 within [1, 2, 4, 5]
            tmp_critical = [1, 2, 4, 5];
            tmp_missed = ismember(tmp_critical, tmp_missing_markers);
            if sum(tmp_missed) > 1
                disp(["trial ", i, "missed too many markers for yz plan!"]);
                continue;
            elseif sum(tmp_missed) > 0
                tmp_yz_markers = tmp_critical( ~tmp_missed );
            end
            % check if there is all within [3, 6, 8]
            tmp_missed = ismember(tmp_xz_markers, tmp_missing_markers);
            if sum(tmp_missed) > 0
                disp(["trial ", i, "missed too many markers for xz plan!"]);
                continue;
            end
            % remove missing marker from markers used to get origin
            tmp_origin_markers = tmp_origin_markers(tmp_cond >= 0);
            
            % ------------------------------------------------------------
            % log missing marker information
            missing_info{i, 1}{tmp_missing, 1} = time_id;
            missing_info{i, 1}{tmp_missing, 2} = find(tmp_cond < 0);
            missing_info{i, 1}{tmp_missing, 3} = (tmp_cond < 0);
            missing_info{i, 1}{tmp_missing, 4} = find(tmp_cond_fgr < 0);
            missing_info{i, 1}{tmp_missing, 5} = (tmp_cond_fgr < 0);
            tmp_missing = tmp_missing + 1;
        end
        
        % compute coordinate bases from yz and xz plan vectors
        tmp_coordBasis = coordBasisOnObj(tmp_markers, tmp_yz_markers, tmp_xz_markers);
        % ----------------------------------------------------------------
        % reconstruct the origins (L and R) using the object rigid dimensions
        tmp_coordOrigin = coordOriginOnObj(tmp_markers, tmp_origin_markers, tmp_coordBasis, obj_side);
        coord_obj{i}{time_id, 1} = [tmp_coordBasis, tmp_coordOrigin];
        
        % compute the finger location in the coordinate frame fixed on the object
        tmp_fgr_on_obj = cell(n_PSonFgr, 1);
        for m = 1:n_PSonFgr
            tmp_fgr_on_obj{m, 1} = (tmp_markers_fgr{m, 1} - transpose( table2array(tmp_coordOrigin) )) * table2array(tmp_coordBasis);
            tmp_cond_fgr(m, 1) = markers_fgr{m, 1}(time_id, 4);
        end
        
        fgr_on_obj{i}{time_id, 1} = tmp_fgr_on_obj;
        fgr_on_obj{i}{time_id, 2} = tmp_cond_fgr;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Get object coordinate before reaching initiation
    % object coordinate before audio go cue and should keep the same until object lift onset*
    ind_b4go = (tmp_audio == 1);
    tmp = cellfun(@table2array, coord_obj{i}(ind_b4go, 1),'UniformOutput',false);
    coord_obj_b4go = array2table(mean(cat(3, tmp{:}), 3), 'RowNames', {'x', 'y', 'z'}, 'VariableNames', {'x_axis', 'y_axis', 'z_axis', 'origin'});
    
    % Angle b/w line (object z) and plane (table xz)
    angTilt_b4go = asind( abs(dot(coord_table_y, coord_obj_b4go{:, 'z_axis'})) / (sqrt(sum(coord_table_y.^2)) * sqrt(sum(coord_obj_b4go{:, 'z_axis'} .^2))) );
    
    
    
    
    %{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% get kinetic data
    switch obj_side
        case 'R'
            tmp_fTh = input{i}{:, var_ATI{1}};
            tmp_fV  = input{i}{:, var_ATI{3}};
            
            tmp_fTh_surface = input_surface{i}{:, var_ATI{1}};
            tmp_fV_surface  = input_surface{i}{:, var_ATI{3}};
        case 'L'
            tmp_fTh = input{i}{:, var_ATI{2}};
            tmp_fV  = input{i}{:, var_ATI{4}};
            
            tmp_fTh_surface = input_surface{i}{:, var_ATI{2}};
            tmp_fV_surface  = input_surface{i}{:, var_ATI{4}};
        otherwise
            disp('something wrong with the file name!!')
    end
    finger_Th{i} = array2table(tmp_fTh, 'VariableNames', kinetic_channels);
    finger_V{i} = array2table(tmp_fV, 'VariableNames', kinetic_channels);
    
    resultantF{i} = array2table(finger_Th{i}{:, kinetic_channels} + finger_V{i}{:, kinetic_channels}, 'VariableNames', kinetic_channels);
    
    finger_Th_surface{i} = array2table(tmp_fTh_surface, 'VariableNames', kinetic_channels);
    finger_V_surface{i} = array2table(tmp_fV_surface, 'VariableNames', kinetic_channels);
    
    %% find th_cop_y and v_cop_y
    %{
    % the ATI data coordinates have been translated and rotated to the
    % center of each respective handle
    z_center2finger = 3 + 21.6 + 3 + 2; % in mm; center piece + ATI Nano25 + mounting + surface
    cop_z = z_center2finger * ones(height(finger_Th), 1); % in mm
    finger_Th.COPz = -cop_z;
    finger_V.COPz = cop_z;
    % COPy = (fy * COPz - mx) / fz
    finger_Th.COPy = ((finger_Th{:, 'fy'} .* finger_Th{:, 'COPz'}) - finger_Th{:, 'mx'})./ finger_Th{:, 'fz'};
    finger_V.COPy = ((finger_V{:, 'fy'} .* finger_V{:, 'COPz'}) - finger_V{:, 'mx'})./ finger_V{:, 'fz'};
    % COPx = (my + fx * COPz) / fz
    finger_Th.COPx = (finger_Th{:, 'my'} + finger_Th{:, 'fx'} .* finger_Th{:, 'COPz'} )./ finger_Th{:, 'fz'};
    finger_V.COPx = (finger_V{:, 'my'} + finger_V{:, 'fx'} .* finger_V{:, 'COPz'} )./ finger_V{:, 'fz'};
    %}
    
    
    %% frame by frame
    lift_marker0 = zeros(height(input{i}), 1);
    for time_id = 1:height(input{i})
        %% compute object tilt
        markers = cell(n_PSonObj, 1);
        for m = 1:n_PSonObj
            markers{m, 1} = input{i}{time_id, var_PS{m}};
        end
        coord_obj{i}{time_id, 1} = coordOnObj(markers, obj_side);
        % Angle b/w line (object z) and plane (table xz)
        angTilt{i}(time_id, 1) = abs(angTilt_b4go - asind( abs(dot(coord_table_y, coord_obj{i}{time_id, 1}{'z_axis', :})) / (sqrt(sum(coord_table_y.^2)) * sqrt(sum(coord_obj{i}{time_id, 1}{'z_axis', :} .^2))) ));
        
        %% compute object height
        lift_marker0(time_id, 1) = sqrt(sum((markers{7, 1} - markers_b4go{7, 1}).^2));
    end
    
    
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
    disp(i);

    %% compute finger tip coordinate without missing frames
    
end

% save(fullfile(pathname, [filename(1:4), '_temp_result.mat']), 'resultantF', 'finger_Th*', 'finger_V*', 'angTilt', 'ind_lft_onset', 'file_list', 'pathname', 'obj_height', 'obj_weight', 'peak_roll', 'peak_mx', 'info_time_trigger');








for i = 88%1:length(missing_info)
    if ~isempty(missing_info{i})
        subplot 211
        plot(cell2mat(missing_info{i, 1}(:, 2)), 'x')
        subplot 212
        plot(cell2mat(missing_info{i, 1}(:, 4)), 'x')
        disp(i)
        pause
    end
end

%%

for j = 1:10%:length(coord_obj)
    tmp = cell(length(coord_obj{j, 1}), 1);
    for i = 1:length(coord_obj{j, 1})
        tmp{i, :} = coord_obj{j, 1}{i, 1}{'origin', :};
    end
    tmp = cell2mat(tmp);
    plot(1:length(tmp), tmp(:, 1), 'r', 1:length(tmp), tmp(:, 2), 'b', 1:length(tmp), tmp(:, 3), 'k');
    disp(j)
    pause
end

%%
for j = 1:10%:length(fgr_on_obj)
    tmp = cell(length(fgr_on_obj{j, 1}), 1);
    for i = 1:length(fgr_on_obj{j, 1})
        tmp{i, :} = fgr_on_obj{j, 1}{i, 1}{:};
    end
    tmp = cell2mat(tmp);
    plot(1:length(tmp), tmp(:, 1), 'r', 1:length(tmp), tmp(:, 2), 'b', 1:length(tmp), tmp(:, 3), 'k');
    disp(j)
    pause
end


