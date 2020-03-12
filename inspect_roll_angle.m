close all; clearvars; clc

%% load aligned data
All_path = uigetdir;
All_filelist = dir(fullfile(All_path, '*_temp_result.mat'));

disp([num2cell((1:length(All_filelist))'), {All_filelist.name}']);

All_selected_sub = input('Which subject(s) to process? ');
if isempty(All_selected_sub)
    All_selected_sub = 1:length(All_filelist);
end

All_roll_ang = cell(length(All_selected_sub), 1);
All_sub_count = 1;
%%
for All_i = All_selected_sub
    clearvars -except All_*;
    sub_id = All_filelist(All_i).name(1:4);
    disp(['Start processing ', sub_id, ' ...']);

    beh = load(fullfile(All_filelist(All_i).folder, All_filelist(All_i).name));
    
    trial_list = {beh.file_list.name}';
    trial_id = cellfun(@(x) x(7:9), trial_list, 'UniformOutput', false);
    trial_cond = cellfun(@(x) x(11:12), trial_list, 'UniformOutput', false);

    trial_cond_count = nan(length(trial_list), 1);
    for i = 1:length(trial_list)
        if strcmpi(trial_cond{i}, 'IL')
            trial_cond_count(i) = str2double(trial_list{i, 1}(16:18));
        elseif strcmpi(trial_cond{i}, 'TR')
            trial_cond_count(i) = str2double(trial_list{i, 1}(13:14)) / 2;
        else
            trial_cond_count(i) = (floor(str2double(trial_list{i, 1}(13:14)) / 2) - 1) * 3 + str2double(trial_list{i, 1}(16:18));
        end
    end
    
    onset_id = nan(length(beh.file_list), 1);
    roll_ang = cell(length(beh.file_list), 3);
    for ep = 1:length(beh.file_list)
        [~, onset_id(ep)] = min( abs(beh.info_onset_time(ep, 1) * 1000 - beh.info_time_trigger{ep, 1}) );
        id_win = (onset_id(ep) - 1000/4):(onset_id(ep) + 2500/4); % -1000 to 2500 ms
        roll_ang{ep, 3} = beh.info_time_trigger{ep, 1}(id_win, :) - beh.info_onset_time(ep, 1) * 1000;
        roll_ang{ep, 1} = beh.angTilt2R{ep, 1}(id_win, 1);
        roll_ang{ep, 2} = trial_cond{ep, 1};
    end
    
    %%
    All_roll_ang{All_sub_count, 1} = sub_id;
    All_roll_ang{All_sub_count, 2} = roll_ang{1, 3};
    % calculate roll angle for IL1, IL2-19, TR1, TR2-19, and PT
    if mod(str2double(sub_id(2:end)), 2) ~= 1
        flip_side = 1;
    else
        flip_side = -1;
    end
    
    
    % IL1
    All_roll_ang{All_sub_count, 2 + 1} = flip_side .* cell2mat(roll_ang(find(strcmpi(roll_ang(:, 2), 'IL'), 1), 1));
    % IL2-19
    tmp_ind = strcmpi(roll_ang(:, 2), 'IL');
    All_roll_ang{All_sub_count, 2 + 2} = flip_side .* cell2mat(roll_ang(tmp_ind(2:end, :), 1)');
    % TR1
    All_roll_ang{All_sub_count, 2 + 3} = flip_side .* cell2mat(roll_ang(find(strcmpi(roll_ang(:, 2), 'TR'), 1), 1));
    % TR2-19
    tmp_ind = strcmpi(roll_ang(:, 2), 'TR');
    All_roll_ang{All_sub_count, 2 + 4} = flip_side .* cell2mat(roll_ang(tmp_ind(2:end, :), 1)');
    % PT
    All_roll_ang{All_sub_count, 2 + 5} = flip_side .* cell2mat(roll_ang(strcmpi(roll_ang(:, 2), 'PT')'));


    All_sub_count = All_sub_count + 1;
end

ncond = size(All_roll_ang, 2) - 2;
ntime = length(All_roll_ang{1, 2});

sig_roll = nan(ncond, ntime);
mu_roll = nan(ncond, ntime);
for c = 1:ncond
    mu_roll(c, :) = mean(cell2mat(All_roll_ang(:, 2 + c)'), 2)';
    sig_roll(c, :) = var(cell2mat(All_roll_ang(:, 2 + c)'), 0, 2)';
end

%% plot
figure
timerstamps = All_roll_ang{1, 2};
h1 = shadedErrorBar(timerstamps, mu_roll(1, :), sqrt(sig_roll(1, :)), '-r', 1);
hold on
h2 = shadedErrorBar(timerstamps, mu_roll(2, :), sqrt(sig_roll(2, :)), '--m', 1);
h3 = shadedErrorBar(timerstamps, mu_roll(3, :), sqrt(sig_roll(3, :)), '-b', 1);
h4 = shadedErrorBar(timerstamps, mu_roll(4, :), sqrt(sig_roll(4, :)), '--c', 1);
h5 = shadedErrorBar(timerstamps, mu_roll(5, :), sqrt(sig_roll(5, :)), '--w', 1);
ylim([-10, 20])
xlim([-1000, 2000])
vline(0, '--r')
xlabel('time (ms)')
ylabel('Roll angle ({\circ})')
hold off
lgnd = legend([h1.mainLine, h2.mainLine, h3.mainLine, h4.mainLine, h5.mainLine], 'IL1', 'IL', 'TR1', 'TR', 'PT');
set(lgnd,'color','none');
set(gca, 'Color', [.8, .8, .8])
title(['Object roll trajectories average across ', num2str(length(All_selected_sub)), ' subjects'])
set(gca, 'FontSize', 24)