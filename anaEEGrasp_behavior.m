close all; clearvars; clc

%% load aligned data
[filename, pathname, ~] = uigetfile;
load(fullfile(pathname, filename));

coord_table_origin = [0, 0, 0];
coord_table_x = [1, 0, 0];
coord_table_y = [0, 1, 0];
coord_table_z = [0, 0, 1];

n_PSonObj = 8;

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
ind_lft_onset = zeros(length(file_list), 6);
audio_trigger = cell(size(file_list));
obj_height = cell(size(file_list));
obj_weight = zeros(length(file_list), 1);
peak_roll = table(zeros(size(file_list)), zeros(size(file_list)), 'VariableNames', {'peakRoll', 'index'});
peak_mx = table(zeros(size(file_list)), zeros(size(file_list)), 'VariableNames', {'peakMx', 'index'});

fingr_th_height = cell(size(file_list));
fingr_in_height = cell(size(file_list));
fingr_mi_height = cell(size(file_list));

% assign filtered or raw data to analyze
input = data_filtered; % or data
input_surface = data_aligned2surface;
for i = 1:length(file_list)
    % get which side of the handle on the object is grasped
    obj_side = file_list(i).name(end-4);
    % get time and triggers
    time = input{i}{:, {'time_ms'}};
    % get trigger indices
    audio_trigger{i, 1} = input{i}{:, {'trigger'}};
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Get object coordinate before reaching initiation
    % object coordinate before audio go cue and should keep the same until object lift onset*
    ind_b4go = (audio_trigger{i, 1} == 1);
    markers_b4go = cell(n_PSonObj, 1);
    for m = 1:n_PSonObj
        markers_b4go{m, 1} = mean(input{i}{ind_b4go, var_PS{m}}, 1);
    end
    
    coord_obj_b4go = coordOnObj(markers_b4go, obj_side);
    
    % 
    
    % Angle b/w line (object z) and plane (table xz)
    angTilt_b4go = asind( abs(dot(coord_table_y, coord_obj_b4go{'z_axis', :})) / (sqrt(sum(coord_table_y.^2)) * sqrt(sum(coord_obj_b4go{'z_axis', :} .^2))) );
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% compute frame by frame
    % get kinetic data
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
    
    lift_marker0 = zeros(height(input{i}), 1);
    for time_id = 1:height(input{i})
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
    
    % based on kinematics
    ind_hold = find(audio_trigger{i, 1} == 3);
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
    ind_hold = find(audio_trigger{i, 1} == 3);
    tmp_stable_window = 50; % in ms
    obj_weight(i, 1) = mean( sqrt(sum(resultantF{i, :}{(ind_hold(end, 1) - tmp_stable_window):ind_hold(end, 1), {'fx', 'fy', 'fz'}}.^2, 2)) );
    
    for j = ind_lft_onset(i, 1):-1:1 % start from 10 mm backward
        tmp_rf = sqrt(sum(resultantF{i, :}{j, {'fx', 'fy', 'fz'}}.^2, 2));
        if tmp_rf < obj_weight(i, 1) % load force equal to object weight
            ind_lft_onset(i, 6) = j;
            break;
        end
    end
    
    
    %% find peak roll after lift onset
    roll_win = 250; % in ms
    ind_roll_win = floor(roll_win ./ (dt * 1000));
    [tmp_roll, tmp_ind] = max(angTilt{i, 1}(ind_lft_onset(i, 1):(ind_lft_onset(i, 1) + ind_roll_win), 1));
    peak_roll{i, {'peakRoll', 'index'}} = [tmp_roll, ind_lft_onset(i, 1) + tmp_ind];
    
    %% find peak mx around lift onset
    pmx_win = 50; % in ms
    ind_pmx_win = floor(pmx_win ./ (dt * 1000));
    [~, tmp_ind] = max( abs(resultantF{i, 1}{(ind_lft_onset(i, 1) - ind_pmx_win):ind_lft_onset(i, 1), 'mx'}) );
    tmp_ind = (ind_lft_onset(i, 1) - ind_pmx_win) + tmp_ind;
    peak_mx{i, {'peakMx', 'index'}} = [resultantF{i, 1}{tmp_ind, 'mx'}, tmp_ind];
    %%
    disp(i);

    %% compute finger tip coordinate without missing frames
    
end

save(fullfile(pathname, [filename(1:4), '_temp_result.mat']), 'resultantF', 'finger_Th*', 'finger_V*', 'angTilt', 'ind_lft_onset', 'file_list', 'pathname', 'obj_height', 'obj_weight', 'peak_roll', 'peak_mx', 'audio_trigger');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check lift onset
tmp_time_real = zeros(length(file_list), size(ind_lft_onset, 2));
for i = 1:length(file_list)
        tmp_time_real(i, j) = data{i}{ind_lft_onset(i, j), 1};
    end
end
delta_t = tmp_time_real - tmp_time_real(:, 6);

str_criterion = {'10', '5', '4', '3', '1'};
figure(1)
for i = 1:5
    subplot(5, 1, i)
    histogram(delta_t(2:end,i))
    ylabel('trial count')
    ntitle(['   onset criterion ', str_criterion{i}, ' mm'], 'location', 'northwest')
    if i ~= 5 
        xlim([-150, 140]);
    end
end
xlabel('time differences b/w height criteria and load equal to weight (ms)')


% % % for i = 1:length(file_list)
% % %     ind_lft_onset
% % % end

%{
for i = 1:length(file_list)
    plotyy(1:length(obj_height{i}), obj_height{i}, 1:length(obj_height{i}), angTilt{i, 1});
    vline(ind_lft_onset(i));
    vline(peak_roll{i, 'index'})
    disp(i)
    pause
end
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check mx
%{
for i = 1:length(file_list)    
    figure(2)
    subplot 211
    plot(obj_height_filtered{i, :})
    hold on
    vline(ind_lft_onset(i, 1));
    hold off
    subplot 212
    plot(mx_filtered{i}{:, 'mx'})
    hold on
    vline(ind_lft_onset(i, 1));
    hold off
    
    disp(i)
    pause
end
%}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot

% Get Torque compensation
mx = resultantF;
% % % cutoff_plot = 5; % in Hz
mx_filtered = mx;
obj_height_filtered = obj_height;
% % % for i = 1:length(mx)
% % %     mx_filtered{i} = filtmat_class( dt, cutoff_plot, mx{i} );
% % %     obj_height_filtered{i} = filtmat_class( dt, cutoff_plot, obj_height{i} );
% % % end

d_center2surface = 3 + 21.6 + 3 + 2; % in mm center piece + Nano 25 + mounting + cover
for i = 1:length(file_list)
    % remove the fy * COPx - fx * COPy effects on mz, 
    % (mz_pure_Th should be close to 0.5 * mz_pure_V) <- this statement is
    % not necessarily true, since individual finger's fy * COPx - fx * COPy
    % can cancel out each other's mz
% % %     COPx_Th = (-finger_Th_surface{i}{:, 'fx'} .* d_center2surface - finger_Th_surface{i}{:, 'my'}) ./ finger_Th_surface{i}{:, 'fz'};
% % %     COPy_Th = (finger_Th_surface{i}{:, 'mx'} - finger_Th_surface{i}{:, 'fy'} .* d_center2surface ) ./ finger_Th_surface{i}{:, 'fz'};
% % %     
% % %     COPx_V = (finger_V_surface{i}{:, 'fx'} .* d_center2surface - finger_V_surface{i}{:, 'my'}) ./ finger_V_surface{i}{:, 'fz'};
% % %     COPy_V = (finger_V_surface{i}{:, 'mx'} + finger_V_surface{i}{:, 'fy'} .* d_center2surface ) ./ finger_V_surface{i}{:, 'fz'};
% % %     
% % %     tmp_mz_pure_Th = finger_Th_surface{i}{:, 'mz'} - (finger_Th_surface{i}{:, 'fy'} .* COPx_Th - finger_Th_surface{i}{:, 'fx'} .* COPy_Th);
% % %     tmp_mz_pure_V = finger_V_surface{i}{:, 'mz'} - (finger_V_surface{i}{:, 'fy'} .* COPx_V - finger_V_surface{i}{:, 'fx'} .* COPy_V);
    

    figure(2)
    subplot 211
<<<<<<< HEAD
    plotyy(1:length(obj_height_filtered{i, :}), data{i, :}{:, 'y10'}, 1:length(obj_height_filtered{i, :}), audio_trigger{i, :})
% % %     vline(ind_lft_onset(i, 1), '-or');
%     plot(obj_height_filtered{i, :})
    plotyy(1:length(obj_height_filtered{i, :}), obj_height_filtered{i, :}, 1:length(obj_height_filtered{i, :}), audio_trigger{i, :})
=======
    plot(obj_height_filtered{i, :})
%     plotyy(1:length(obj_height_filtered{i, :}), obj_height_filtered{i, :}, 1:length(obj_height_filtered{i, :}), audio_trigger{i, :})
>>>>>>> 74cc1390dd5fd2a3911b6cfe48b2d920543474c0
    hold on
    vline(ind_lft_onset(i, 1), '-or');
    vline(ind_lft_onset(i, 2), ':or');
    vline(ind_lft_onset(i, 3), '-ok');
    vline(ind_lft_onset(i, 4), '-ob');
    vline(ind_lft_onset(i, 5), ':ob');
    vline(ind_lft_onset(i, 6), '-k');
    hold off
    xlim([ind_lft_onset(i, 4) - 100, ind_lft_onset(i, 1) + 100])
    ylim([0, 15])
    subplot 212
    plot(resultantF{i, :}{:, 'fy'})
%     plotyy(1:length(resultantF{i, :}{:, 'fy'}), resultantF{i, :}{:, 'fy'}, 1:length(resultantF{i, :}{:, 'fy'}), audio_trigger{i, :})
    hold on
    vline(ind_lft_onset(i, 1), '-or');
    vline(ind_lft_onset(i, 2), ':or');
    vline(ind_lft_onset(i, 3), '-ok');
    vline(ind_lft_onset(i, 4), '-ob');
    vline(ind_lft_onset(i, 5), ':ob');
    vline(ind_lft_onset(i, 6), '-k');
    hline(obj_weight(i, 1));
    hold off
    xlim([ind_lft_onset(i, 4) - 100, ind_lft_onset(i, 1) + 100])
    
% % %     
% % %     plot(1:length(tmp_mz_pure_Th), tmp_mz_pure_Th, '-r', 1:length(tmp_mz_pure_Th), tmp_mz_pure_V, '-b')
% % %     hold on
% % %     vline(ind_lft_onset(i, 1));
% % %     hold off
% % %     ylim([-100, 100])
    
    disp(i)
    pause
end




%% Plot peak mx and peak roll around lift onset in all trials
%{
tmp = [str2double(file_list(1).name(7:9)), peak_mx{1, 'peakMx'}, peak_roll{1, 'peakRoll'}];
session = cell(33, 2);
j = 1;
for i = 2:length(file_list)
    time_of_plot = ind_lft_onset(i, 1);
    if (file_list(i).name(:, 11) ~= file_list(i - 1).name(:, 11))
        session(j, :) = {file_list(i - 1).name(:, 11), tmp};
        j = j + 1;
        tmp = [];
    end
    
    tmp = [tmp; str2double(file_list(i).name(7:9)), peak_mx{i, 'peakMx'}, peak_roll{i, 'peakRoll'}];
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
