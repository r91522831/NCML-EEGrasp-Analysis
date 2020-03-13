close all; clearvars; clc

%% load aligned data
All_path = uigetdir;
All_filelist = dir(fullfile(All_path, 'sub*'));

disp([num2cell((1:length(All_filelist))'), {All_filelist.name}']);
All_selected_sub = input('Which subject(s) to process? ');

if isempty(All_selected_sub)
    All_selected_sub = 1:length(All_filelist);
end

%%
nsub = length(All_selected_sub);
All_info_trial = cell(nsub, 3);
All_sub_i = 1;
for All_i = All_selected_sub
    clearvars -except All_*; close all;
    sub_id = All_filelist(All_i).name(end-1:end);
    disp(['Start processing sub-', sub_id, ' ...']);
    
    beh = load(fullfile(All_path, All_filelist(All_i).name, 'beh', 'mat', ['S0', sub_id, '_temp_result.mat']));
    nep = length(beh.file_list);
    info_trial = cell(nep, 8);
    for ep = 1:length(beh.file_list)
        info_trial{ep, 1} = str2double(beh.file_list(ep).name(7:9)); % trial id
        info_trial{ep, 2} = beh.file_list(ep).name(11:12); % condition: IL, TR, PT, TR
        info_trial{ep, 3} = str2double(beh.file_list(ep).name(23:25)); % added weight in gram
        info_trial{ep, 4} = str2double(beh.file_list(ep).name(16:18)); % trial number in block
        
        f_mag_th = sqrt(beh.finger_Th{ep, 1}.fy .^2 + beh.finger_Th{ep, 1}.fz .^2); % F_mag_TH
        f_ang_th = atan2(beh.finger_Th{ep, 1}.fy, beh.finger_Th{ep, 1}.fz); % F_ang_TH
        f_mag_vf = sqrt(beh.finger_V{ep, 1}.fy .^2 + beh.finger_V{ep, 1}.fz .^2); % F_mag_VF
        f_ang_vf = atan2(beh.finger_V{ep, 1}.fy, beh.finger_V{ep, 1}.fz); % F_ang_VF
        
        d_copy_thvf = beh.finger_Th{ep, 1}.COPy - beh.finger_V{ep, 1}.COPy; % delta COPy TH - VF
        
        info_trial{ep, 5} = [f_mag_th, f_ang_th, f_mag_vf, f_ang_vf, d_copy_thvf];
        info_trial{ep, 6} = beh.info_onset_time(ep, :);
    end
    
    info_trial(:, 7) = mat2cell([beh.peak_roll{:, 'peakRoll'}], ones(nep, 1));
    info_trial(:, 8) = mat2cell([beh.peak_mx{:, 'peakMx'}], ones(nep, 1));
    
    
    All_info_trial{All_sub_i, 1} = ['sub-', sub_id];
    All_info_trial{All_sub_i, 2} = info_trial;
    All_info_trial{All_sub_i, 3} = diff(beh.info_time_trigger{1, 1}(1:2, 1)); % dt in ms
    All_sub_i = All_sub_i + 1;
end

%% plot force magnitude, force angle, and delta COPy for TH and VF in IL, TR, PT conditions
sub_i = 1;
sub_id = All_info_trial{sub_i, 1};
data = All_info_trial{sub_i, 2};
dt = All_info_trial{sub_i, 3}; % in ms

tmp_win = round(-1000/dt):round(2000/dt);
tch_time = nan(length(data), 1);
f_mag_th = nan(length(data), length(tmp_win));
f_mag_vf = f_mag_th; f_ang_th = f_mag_th; f_ang_vf = f_mag_th; d_copy_thvf = f_mag_th;
for ep = 1:length(data)
    lft_ind = data{ep, 6}{1, 'lft_ind'};
    lft_time = data{ep, 6}{1, 'lft_time'};
    tch_time(ep, 1) = data{ep, 6}{1, 'tch_time'} - lft_time;
    lft_win = ( lft_ind + tmp_win )';
    lft_win_time = tmp_win * dt * 0.001; % in s
    
    f_mag_th(ep, :) = data{ep, 5}(lft_win, 1)';
    f_mag_vf(ep, :) = data{ep, 5}(lft_win, 3)';
    f_ang_th(ep, :) = data{ep, 5}(lft_win, 2)';
    f_ang_vf(ep, :) = data{ep, 5}(lft_win, 4)';
    d_copy_thvf(ep, :) = data{ep, 5}(lft_win, 5)';
end

%%
ep_ind = [strcmpi(data(:, 2), 'IL'), strcmpi(data(:, 2), 'TR'), strcmpi(data(:, 2), 'PT')];
tmp_lft_win_time = repmat(lft_win_time, length(data), 1);
tt = {'IL', 'TR', 'PT'};

ep_ind = ep_ind(1:18, :);

figure('DefaultAxesFontSize', 18)
for cond = 1:3
    nep = sum(ep_ind(:, cond));
    carray = flip([logspace(0, 1, nep), logspace(0, 1, nep)])./10;
    
    subplot(2, 3, cond)
    p = polarplot(f_ang_th(ep_ind(:, cond), :)', f_mag_th(ep_ind(:, cond), :)', f_ang_vf(ep_ind(:, cond), :)', f_mag_vf(ep_ind(:, cond), :)', '--');
    for i = 1:length(p)
        p(i).Color = [carray(i), carray(i), carray(i)];
    end
    title(tt{cond})
    set(gca, 'FontSize', 18)
    
    subplot(2, 3, cond + 3)
    p = plot(tmp_lft_win_time(ep_ind(:, cond), :)', d_copy_thvf(ep_ind(:, cond), :)');
    for i = 1:length(p)
        p(i).Color = [carray(i), carray(i), carray(i)];
    end
    ylim([-50, 50])
% % %     vline(tch_time, ':k')
    xlim([-1, 2])
    vline(0, '--r', 'lift')
end
mtit(sub_id)

%%
ep_plot = 1:6; % 19:24; % 7:12; %1:6; % 7:12; % 13:18;
nep = length(ep_plot);
figure('DefaultAxesFontSize', 18, 'units', 'normalized', 'outerposition', [0, 0, 1, 1])
id_plot = 1;
for ep = ep_plot
    subplot(3, nep, id_plot)
    plot(lft_win_time, f_mag_th(ep, :), lft_win_time, f_mag_vf(ep, :))
    ylim([0, 50])
    xlim([-1, 2])
    vline([0, tch_time(ep, 1)], {'--r', '--k'}, {'lift', 'touch'})
    if id_plot == 1
        ylabel('force (N)')
    end
    title(['t', num2str(data{ep, 1}), ' cond ', data{ep, 2}, ' ', num2str(data{ep, 4})])
    
    subplot(3, nep, nep + id_plot )
    plot(lft_win_time, rad2deg(f_ang_th(ep, :)), lft_win_time, rad2deg(f_ang_vf(ep, :)))
    ylim([-180, 180])
    xlim([-1, 2])
    vline([0, tch_time(ep, 1)], {'--r', '--k'}, {'lift', 'touch'})
    if id_plot == 1
        ylabel('angle ({\circ})')
        legend({'TH', 'VF'})
    end
    
    subplot(3, nep, 2 * nep + id_plot)
    plot(lft_win_time, d_copy_thvf(ep, :));
    ylim([-30, 30])
    xlim([-1, 2])
    vline([0, tch_time(ep, 1)], {'--r', '--k'}, {'lift', 'touch'})
    ylabel('time (s)')
    id_plot = id_plot + 1;
end
suptitle(sub_id)

%%
ep_plot = 1:30; % 19:24; % 7:12; %1:6; % 7:12; % 13:18;
nep = length(ep_plot);
figure('DefaultAxesFontSize', 18, 'units', 'normalized', 'outerposition', [0, 0, 1, 1])

for ep = ep_plot
    subplot(5, 6, ep)
    plot(lft_win_time, f_mag_th(ep, :), lft_win_time, f_mag_vf(ep, :))
    ylim([0, 50])
    xlim([-1, 2])
    
    if ep == 1
        ylabel('force (N)')
        vline([0, tch_time(ep, 1)], {'--r', '--k'}, {'lift', 'touch'})
    else
        vline([0, tch_time(ep, 1)], {'--r', '--k'})
    end
    if ep <= 6
        title(['t', num2str(data{ep, 1}), ' cond ', data{ep, 2}, ' ', num2str(data{ep, 4})])
    end
end
suptitle(sub_id)