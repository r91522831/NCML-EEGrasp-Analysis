close all; clearvars; clc

%% load aligned data
[filename, pathname, ~] = uigetfile;
filelist = dir(fullfile(pathname, [filename(1:4), '*.mat']));
for i = 1:length(filelist)
    load(fullfile(pathname, filelist(i).name));
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check lift onset
tmp_time_real = zeros(length(file_list), size(ind_lft_onset, 2));
for i = 1:length(file_list)
    for j = 1:size(ind_lft_onset, 2)
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
% % %     vline(ind_lft_onset(i, 1), '-or');
    plot(obj_height_filtered{i, :})
%     plotyy(1:length(obj_height_filtered{i, :}), obj_height_filtered{i, :}, 1:length(obj_height_filtered{i, :}), info_time_trigger{i, 2}{:})

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
%     plotyy(1:length(resultantF{i, :}{:, 'fy'}), resultantF{i, :}{:, 'fy'}, 1:length(resultantF{i, :}{:, 'fy'}), info_time_trigger{i, 2}{:})
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