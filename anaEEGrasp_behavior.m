close all; clearvars; clc

%% load aligned data
[filename, pathname, ~] = uigetfile;
load(fullfile(pathname, filename));

coord_table_origin = [0, 0, 0];
coord_table_x = [1, 0, 0];
coord_table_y = [0, 1, 0];
coord_table_z = [0, 0, 1];

n_PSonObj = 8;

%% lowpass filter all data
dt = diff(data{1}{1:2, 1}) * 0.001; % in second
% % % cutoff = 30; % in Hz
data_filtered = data;
% % % for i = 1:length(file_list)
% % %     data_filtered{i}{:, 3:end} = filtmat_class( dt, cutoff, data{i}{:, 3:end} );
% % % end

%%
resultant_mx = cell(size(file_list));
angTilt = cell(size(file_list));
coord_obj = cell(size(file_list));
ind_lft_onset = zeros(size(file_list));
obj_height = cell(size(file_list));

% assign filtered or raw data to analyze
input = data_filtered; % or data

for i = 1:length(file_list)
    % get which side of the handle on the object is grasped
    obj_side = file_list(i).name(end-4);
    % get time and triggers
    time = input{i}{:, {'time_ms'}};
    % get trigger indices
    audio_trigger = input{i}{:, {'trigger'}};
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Get object coordinate before reaching initiation
    % object coordinate before audio go cue and should keep the same until object lift onset*
    ind_b4go = (audio_trigger == 1);
    markers_b4go = cell(n_PSonObj, 1);
    for m = 1:8
        markers_b4go{m, 1} = mean(input{i}{ind_b4go, var_PS{m}}, 1);
    end
    
    coord_obj_b4go = coordOnObj(markers_b4go, obj_side);
    
    % Angle b/w line (object z) and plane (table xz)
    angTilt_b4go = asind( abs(dot(coord_table_y, coord_obj_b4go{'z_axis', :})) / (sqrt(sum(coord_table_y.^2)) * sqrt(sum(coord_obj_b4go{'z_axis', :} .^2))) );
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% compute frame by frame
    lift_marker0 = zeros(height(input{i}), 1);
    for time_id = 1:height(input{i})
        %% compute Torque compensation
        switch obj_side
            case 'R'
                finger_Th = input{i}{time_id, var_ATI{1}};
                finger_V  = input{i}{time_id, var_ATI{3}};
            case 'L'
                finger_Th = input{i}{time_id, var_ATI{2}};
                finger_V  = input{i}{time_id, var_ATI{4}};
            otherwise
                disp('something wrong with the file name!!')
        end

        % find th_cop_y and v_cop_y
        % since the ATI reference origin is on tool surface: dy = (mx - fy * dz) / fz
        th_cop_z2surface = -5; % in mm
        v_cop_z2surface = 5; % in mm
        th_cop_y = (finger_Th(:, 4) - finger_Th(:, 2) .* th_cop_z2surface )./ finger_Th(:, 3);
        v_cop_y = (finger_V(:, 4) - finger_V(:, 2) .* v_cop_z2surface )./ finger_V(:, 3);

        % width b/w two nano 25 tool surface in mm
        width_obj = 2 * (21.6 + 3);

        diff_finger = finger_Th - finger_V;
        resultant_mx{i}(time_id, 1) = diff_finger(:, 2) .* width_obj - (diff_finger(:, 3) .* (th_cop_y - v_cop_y));
    
        %% compute object tilt
        markers = cell(n_PSonObj, 1);
        for m = 1:8
            markers{m, 1} = mean(input{i}{time_id, var_PS{m}}, 1);
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
    
    ind_hold = find(audio_trigger == 3);
    for j = ind_hold(end, 1):-1:1
        if obj_height{i, :}(j, 1) < 10 % 10 mm
            ind_lft_onset(i, 1) = j;
            break;
        end
    end

    disp(i);

    %% compute finger tip coordinate without missing frames
    
end

save(fullfile(pathname, [filename(1:4), '_temp_result.mat']), 'resultant_mx', 'angTilt', 'ind_lft_onset', 'file_list', 'pathname', 'obj_height')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check lift onset

for i = 1:length(file_list)
    plot(obj_height{i});
    vline(ind_lft_onset(i));
    disp(i)
    pause
end



%{
%% plot
mx = resultant_mx;
% % % cutoff_plot = 5; % in Hz
mx_filtered = mx;
obj_height_filtered = obj_height;
% % % for i = 1:length(mx)
% % %     mx_filtered{i} = filtmat_class( dt, cutoff_plot, mx{i} );
% % %     obj_height_filtered{i} = filtmat_class( dt, cutoff_plot, obj_height{i} );
% % % end


%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
%% find peak roll after lift onset
ind_peak_roll = zeros(size(file_list));
for i = 1:length(file_list)
    roll_win = 250; % in ms
    ind_roll_win = floor(roll_win ./ (dt * 1000));
    
    [~, tmp_ind] = max(angTilt{i, 1}(ind_lft_onset(i, 1):(ind_lft_onset(i, 1) + ind_roll_win), 1));
    ind_peak_roll(i, 1) = tmp_ind;
end

for i = 1:length(file_list)
    figure(2)
    subplot 211
    plot(obj_height_filtered{i, :})
    hold on
    vline(ind_lft_onset(i, 1));
    hold off
    subplot 212
    plot(mx_filtered{i, :})
    hold on
    vline(ind_lft_onset(i, 1));
    hold off
    
%     disp(i)
%     pause
end
%}




%%

time_of_plot = ind_lft_onset(1, 1);
tmp = [str2double(file_list(1).name(7:9)), mx_filtered{1}(time_of_plot, 1), angTilt{1}(ind_peak_roll(1, 1), 1)];
session = cell(33, 2);
j = 1;
for i = 2:length(file_list)
    time_of_plot = ind_lft_onset(i, 1);
    if (file_list(i).name(:, 11) ~= file_list(i - 1).name(:, 11))
        session(j, :) = {file_list(i - 1).name(:, 11), tmp};
        j = j + 1;
        tmp = [];
    end
    
    tmp = [tmp; str2double(file_list(i).name(7:9)), mx_filtered{i}(time_of_plot, 1), angTilt{i}(ind_peak_roll(i, 1), 1)];
end
session(j, :) = {file_list(end).name(:, 11), tmp};

figure(1)
subplot 211
hold on
for i = 1:length(session)
    switch session{i, 1}
        case 'I'
            line_spec = '-*r';
        case 'T'
            line_spec = '-ob';
        case 'P'
            line_spec = '-xk';
        otherwise
    end
    plot(session{i, 2}(:, 1), session{i, 2}(:, 2), line_spec)
end
hold off
subplot 212
hold on
for i = 1:length(session)
    switch session{i, 1}
        case 'I'
            line_spec = '-*r';
        case 'T'
            line_spec = '-ob';
        case 'P'
            line_spec = '-xk';
        otherwise
    end
    plot(session{i, 2}(:, 1), session{i, 2}(:, 3), line_spec)
end
hold off
% ylim([0, 15])
legend({'IL', 'TR', 'PT'})
%}
