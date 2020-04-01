close all; clear; clc;
%%
All_path = uigetdir;
All_filelist = dir(fullfile(All_path, '*_timefreq.mat'));

disp([num2cell((1:length(All_filelist))'), {All_filelist.name}']);
All_selected_sub = input('Which subject(s) to run correlation? ');
if isempty(All_selected_sub)
    All_selected_sub = 1:length(All_filelist);
end

% channel id:
% {'Fz', 'FCz', 'C3', 'CP3', 'C1'} = {6, 41, 15, 47, 44}
Fz = 6; FCz = 41; C3 = 15; CP3 = 47; C1 = 44;
electrodes_name = {'Fz'; 'FCz'; 'C3'; 'CP3'; 'C1'};
electrodes = {Fz; FCz; C3; CP3; C1};
disp([electrodes(:), electrodes_name(:)]);
All_selected_elec = input('Which electrode(s) to process? ');
if isempty(All_selected_elec)
    All_selected_elec = 1:63; % 63 electrodes
end

%%
All_summary = load('/Users/yenhsunw/Dropbox (ASU)/BIDS_format/NCML-EEGrasp/behavior_resultSummary_w_dist2ideal.mat');

All_power_dist = cell(length(All_selected_sub), 3);
for All_sub_i = All_selected_sub
    clearvars -except All_*;
    subID = All_filelist(All_sub_i).name(1:6);
    disp(['start correlating ', subID])
    
    % tf = load('/Users/yenhsunw/Dropbox (ASU)/NCML-EEGrasp/EEG/eeglab/004 RemoveEye90/linear/RAW/sub-02_timefreq.mat');
    tf = load(fullfile(All_path, All_filelist(All_sub_i).name));

    subID_inSummary = strcmpi(All_summary.All_info_trial.subID, subID);
    if ~any(subID_inSummary)
        disp([subID, ' has no behavior data'])
        continue;
    end
    bdata = All_summary.All_info_trial.data{ subID_inSummary, 1 }; % get the right 
        
    times = tf.tf_times';
    freqz = tf.tf_freqs';
    elect = tf.tf_ersp.Properties.RowNames;
    
    ntime = length(times);
    nfreq = length(freqz);
    nelec = length(elect);
    nep   = size(tf.tf_ersp{1, 1}{:}, 3);
    
    % theta band: 4 to 8 Hz; alpha: 8 ~ 13 Hz; low beta: 13 ~ 19 Hz; high beta: 20 ~ 30 Hz
    rg_freq = {find(freqz > 4 & freqz <= 8), find(freqz > 8 & freqz <= 13), find(freqz > 13 & freqz <= 20), find(freqz > 20 & freqz <= 30) };
    rg_freq_name = {'\theta', '\alpha', 'low\beta', 'high\beta'};
    nfreqband = length(rg_freq);
    
    % define a window size of ~300 ms moving across the range from 0 to 5000 ms
    rg_all_time = times(times >= -3000 & times < 5000);
    rg_ind_dt = round(50 / mean(diff(rg_all_time))); % find number of indices corresponding to 50 ms window
    n_mv_win = length(rg_all_time) - rg_ind_dt;
    rg_50_time = cell(n_mv_win, 1);
    for i_win = 1:n_mv_win
        [~, rg_50_time{i_win, 1}] = ismember(rg_all_time(i_win:(i_win + rg_ind_dt)), times); % find the indices in the timerstamps
    end
    
    % correlation for each electrode
    sub_power_dist = cell2table(cell(nelec, nfreqband), 'VariableNames', {'theta', 'alpha', 'lowBeta', 'highBeta'}, 'RowNames', tf.tf_ersp.Properties.RowNames);
    for i_elect = All_selected_elec
        for i_freqband = 1:nfreqband
            tmp_power_mw = cell(n_mv_win, 2);
            % mean with a window size of ~300 ms moving across the range from 0 to 5000 ms for a freqency band
            for i_win = 1:n_mv_win
                for i_ep = 1:nep
                    tmp_p_mw = abs(tf.tf_ersp.Var1{elect(i_elect), 1}(rg_freq{i_freqband}, rg_50_time{i_win, 1}, i_ep));
                    [~, tmp_power_mw{i_win, 1}(i_ep, 1)] = robustcov(tmp_p_mw(:));
                end
                
                d2i = nan(nep, 1);
                d2i(bdata.trial, 1) = bdata.dist2ideal;
                tmp_power_mw{i_win, 1}(:, 2) = d2i;
                
                [~, ~, ~, tmp_outliers] = robustcov(bdata.dist2ideal);
                tmp_r_mw = corr(tmp_power_mw{i_win, 1}(~tmp_outliers, :));
                tmp_power_mw{i_win, 2} = tmp_r_mw(1, 2);
            end
            
            sub_power_dist{i_elect, i_freqband} = {tmp_power_mw};
        end
    end
    
    save(fullfile(All_path, [subID, '_corr_traj.mat']), 'sub_power_dist');
    disp([subID, ' finished.']);
    
    All_power_dist{All_sub_i, 1} = subID;
    All_power_dist{All_sub_i, 2} = sub_power_dist;
    All_power_dist{All_sub_i, 3} = rg_all_time(1:n_mv_win);
end

save(fullfile(All_path, ['corr_traj_', num2str(length(All_selected_sub)), 'sub.mat']), 'All_power_dist');

%
% % % plot(rg_all_time(1:n_mv_win), cell2mat(tmp_power_mw(:,end)))