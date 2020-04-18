clear; close all; clc;

%% figure save location
dir_save = '/Users/yen-hsunwu/Dropbox (ASU)/NCML-EEGrasp/for Apr21 2020';

%% Get behavior data: Mcom and peak Roll
% % % beh = load('/Users/yen-hsunwu/Dropbox (ASU)/BIDS_format/NCML-EEGrasp/sub-09/beh/mat/S009_temp_result.mat');
%% Get time freq data
% % % tf = load('/Users/yen-hsunwu/Dropbox (ASU)/NCML-EEGrasp/EEG/eeglab/004 RemoveEye90/linear/RAW/sub-09_timefreq.mat');
%% Get EEG data: voltage for each channel
% load EEG csd data. Raw EEG data are in EEG.dataRaw. CSD EEG data are in EEG.data
% % % EEG = pop_loadset('filename', 'sub-09_eeg_csd.set', 'filepath', '/Users/yen-hsunwu/Dropbox (ASU)/NCML-EEGrasp/EEG/eeglab/004 RemoveEye90/sub-09_onset/');

%{
beh = load('/Users/yen-hsunwu/Dropbox (ASU)/BIDS_format/NCML-EEGrasp/sub-02/beh/mat/S002_temp_result.mat');
tf = load('/Users/yen-hsunwu/Dropbox (ASU)/NCML-EEGrasp/EEG/eeglab/004 RemoveEye90/linear/RAW/sub-02_timefreq.mat');
EEG = pop_loadset('filename', 'sub-02_eeg_csd.set', 'filepath', '/Users/yen-hsunwu/Dropbox (ASU)/NCML-EEGrasp/EEG/eeglab/004 RemoveEye90/sub-02_onset/');
%}

beh = load('/Users/yen-hsunwu/Dropbox (ASU)/BIDS_format/NCML-EEGrasp/sub-09/beh/mat/S009_temp_result.mat');
tf = load('/Users/yen-hsunwu/Dropbox (ASU)/NCML-EEGrasp/EEG/eeglab/004 RemoveEye90/linear/RAW/sub-09_timefreq.mat');
EEG = pop_loadset('filename', 'sub-09_eeg_csd.set', 'filepath', '/Users/yen-hsunwu/Dropbox (ASU)/NCML-EEGrasp/EEG/eeglab/004 RemoveEye90/sub-09_onset/');
%}
%{
beh = load('/Users/yen-hsunwu/Dropbox (ASU)/BIDS_format/NCML-EEGrasp/sub-11/beh/mat/S011_temp_result.mat');
tf = load('/Users/yen-hsunwu/Dropbox (ASU)/NCML-EEGrasp/EEG/eeglab/004 RemoveEye90/linear/RAW/sub-11_timefreq.mat');
EEG = pop_loadset('filename', 'sub-11_eeg_csd.set', 'filepath', '/Users/yen-hsunwu/Dropbox (ASU)/NCML-EEGrasp/EEG/eeglab/004 RemoveEye90/sub-11_onset/');
%}
% correct the channel locations to standard 10-5 cap
EEG = pop_chanedit(EEG, 'lookup', '/Users/yen-hsunwu/Dropbox (Personal)/Programming/Matlab/myLibrary/eeglab2019_1/plugins/dipfit/standard_BESA/standard-10-5-cap385.elp');
subID = EEG.filename(1:6);

%% match trial IDs with trial number
ind_trial_eeg = cellfun(@(x) x{1}, {EEG.epoch.eventtrialID})';
ind_trial_beh = cellfun(@(x) str2double(x(7:9)), {beh.file_list.name})';

%% make EEG data matrix epoch count the same as the trial id
nep_eeg = EEG.trials;
nelec = EEG.nbchan;
ntimes_eeg = length(EEG.times);
data_eeg = nan(nelec, ntimes_eeg, 95);
for i = 1:nep_eeg
    data_eeg(:, :, i) = EEG.dataRaw(:, :, ind_trial_eeg(i));
end

%% compile frequency band complex number into eeg data formate nelectrode x ntimes x ntrials x nfreqBands
freqz = tf.tf_freqs;
rg_freq_band = { {'\theta', '4~8 Hz', 'theta'}, {'\alpha', '8~13 Hz', 'alpha'}, {'\beta_{low}', '13~20 Hz', 'betalow'}, {'\beta_{high}', '20~30 Hz', 'betahigh'}; ...
                 find(freqz >  4 & freqz <=  8), find(freqz >   8 & freqz <= 13), ...
                 find(freqz > 13 & freqz <= 20), find(freqz >  20 & freqz <= 30) };
nfb = length(rg_freq_band);
ntimes_tf = length(tf.tf_times);
data_tf = cell(nfb, 1);
for i_fb = 1:nfb
    for i_elec = 1:nelec
        data_tf{i_fb, 1}(i_elec, :, :) = mean(tf.tf_ersp.Var1{i_elec, 1}(rg_freq_band{2, i_fb}, :, :), 1);
    end
end

%% block trials based on the conditions: IL, TR, PT
ind_trial_condType_eeg = cellfun(@(x) x(1), {EEG.epoch.eventcondType})';
ind_trial_condType = cell(95, 1);
for i = 1:nep_eeg
    ind_trial_condType(ind_trial_eeg(i)) = ind_trial_condType_eeg(i);
end
block_trialID = [strcmpi(ind_trial_condType, 'IL'), strcmpi(ind_trial_condType, 'TR'), strcmpi(ind_trial_condType, 'PT')];
%% match the trial ids between the behavior trials and the eeg trials
ind_trial = nan(95, 2); % the 1st column is the index in beh data, the 2nd is the index in eeg data
for i = 1:length(ind_trial_beh), ind_trial(ind_trial_beh(i), 1) = i; end
for i = 1:length(ind_trial_eeg), ind_trial(ind_trial_eeg(i), 2) = i; end
%% find the time indices of touch onset, lift onset and peak roll
ind_lft_eeg = dsearchn(EEG.times', 0);
ind_lft_tf = dsearchn(tf.tf_times', 0);

tag_time = nan(length(ind_trial), 2);
ind_tag_time_eeg = nan(length(ind_trial), 3);
ind_tag_time_tf = nan(length(ind_trial), 3);
for i = 1:length(ind_trial)
    if ~any(isnan(ind_trial(i, :)))
        tag_time(i, 1) = beh.info_onset_time.tch_time(i) - beh.info_onset_time.lft_time(i); % in second
        tag_time(i, 2) = beh.info_onset_time.pRoll_time(i) - beh.info_onset_time.lft_time(i); % insecond
        
        ind_tag_time_eeg(i, :) = dsearchn(EEG.times', 1000 * [tag_time(i, 1), 0, tag_time(i, 2)]');
        ind_tag_time_tf(i, :) = dsearchn(tf.tf_times', 1000 * [tag_time(i, 1), 0, tag_time(i, 2)]');
    end
end
%% Rearrange EEG raw voltage into blocks of condition types
data_re_eeg = cat( 3, data_eeg(:, :, block_trialID(:, 1)), nan(nelec, ntimes_eeg, 1), ...
                      data_eeg(:, :, block_trialID(:, 2)), nan(nelec, ntimes_eeg, 1), ...
                      data_eeg(:, :, block_trialID(:, 3)) );

ind_tag_time_block_eeg = vertcat( ind_tag_time_eeg(block_trialID(:, 1), :), nan(1, 3), ...
                                  ind_tag_time_eeg(block_trialID(:, 2), :), nan(1, 3), ...
                                  ind_tag_time_eeg(block_trialID(:, 3), :) );
       
data_re_tf = cell(nfb, 1);
for i_fb = 1:nfb
    data_re_tf{i_fb, 1} = cat( 3, data_tf{i_fb, 1}(:, :, block_trialID(:, 1)), nan(nelec, ntimes_tf, 1), ...
                                  data_tf{i_fb, 1}(:, :, block_trialID(:, 2)), nan(nelec, ntimes_tf, 1), ...
                                  data_tf{i_fb, 1}(:, :, block_trialID(:, 3)) );
    
end
ind_tag_time_block_tf = vertcat( ind_tag_time_tf(block_trialID(:, 1), :), nan(1, 3), ...
                                 ind_tag_time_tf(block_trialID(:, 2), :), nan(1, 3), ...
                                 ind_tag_time_tf(block_trialID(:, 3), :) );
                          
%% get time indices for -1000 ms and 2000 ms from lift onset
ind_win_eeg = dsearchn(EEG.times', [-1000, 2000]');
ind_win_tf = dsearchn(tf.tf_times', [-1000, 2000]');

%% baseline signals in each trial with the the avg(touch to lift) in that trial
% voltage
data_bl_sub_eeg = nan(size(data_re_eeg));
for i_elec = 1:nelec
    for i_ep = 1:size(data_re_eeg, 3)
        if ~any(isnan(ind_tag_time_block_eeg(i_ep, 1:2)))
            baseline = mean(data_re_eeg(i_elec, ind_tag_time_block_eeg(i_ep, 1:2), i_ep)); % average from touch to lift onset
            data_bl_sub_eeg(i_elec, :, i_ep) = data_re_eeg(i_elec, :, i_ep) - baseline;
        end
    end
end
% tf
data_bl_sub_tf = cell(size(data_re_tf));
for i_fb = 1:nfb
    for i_elec = 1:nelec
        for i_ep = 1:size(data_re_tf{i_fb, 1}, 3)
            if ~any(isnan(ind_tag_time_block_tf(i_ep, 1:2)))
                baseline = mean(data_re_tf{i_fb, 1}(i_elec, ind_tag_time_block_tf(i_ep, 1:2), i_ep)); % average from touch to lift onset
                data_bl_sub_tf{i_fb, 1}(i_elec, :, i_ep) = data_re_tf{i_fb, 1}(i_elec, :, i_ep) - baseline;
            end
        end
    end
end

%% plotting
baseline_corrected = true;
range_fixed = true;

% plot for voltage
figureSize = 3.3 / nelec;
scale = .6;
if baseline_corrected
    dRange_eeg = nanmean(prctile(prctile(abs(data_bl_sub_eeg), 95), 95));
else
    dRange_eeg = nanmean(prctile(prctile(abs(data_re_eeg), 95), 95));
end
if range_fixed
    dRange_eeg = 30;
end

fig = figure('DefaultAxesFontSize', 18, 'units', 'normalized', 'outerposition', [0, 0, .67, 1]);
for i_elec = 1:nelec
        hh(i_elec) = axes( 'position',...
                           [ sin( deg2rad(EEG.chanlocs(i_elec).theta) ) * EEG.chanlocs(i_elec).radius * scale + .5 - (figureSize / 2), ...
                             cos( deg2rad(EEG.chanlocs(i_elec).theta) ) * EEG.chanlocs(i_elec).radius * scale + .5 - (figureSize / 2), ...
                             figureSize, ...
                             figureSize ]' );
        if baseline_corrected
            imagesc(squeeze(data_bl_sub_eeg(i_elec, ind_win_eeg(1):ind_win_eeg(2), :))', [-dRange_eeg, dRange_eeg]);
        else
            imagesc(squeeze(data_re_eeg(i_elec, ind_win_eeg(1):ind_win_eeg(2), :))', [-dRange_eeg, dRange_eeg]);
        end
        %{
        for i = 1:size(data, 3)
            text(ind_tag_time_block(i, 1) - ind_win(1), i, 't')
            % % %     text(ind_tag_time_block(i, 3) - ind_win(1), i, 'r');
        end
        %}
        vline(ind_lft_eeg - ind_win_eeg(1));
        set(gca, 'YDir', 'normal');
% % %         set(gca, 'YTickLabel', []);
        axis off;
end

if baseline_corrected
    suptitle([subID, ' voltage ', 'w/ baseline corrected, range \pm', num2str(dRange_eeg), ' \muV'])
    filename = [subID, '_voltage_w_baseline_subtraction_', num2str(round(dRange_eeg))];
else
    suptitle([subID, ' voltage ', 'w/o baseline corrected, range \pm', num2str(dRange_eeg), ' \muV'])
    filename = [subID, '_voltage_wo_baseline_subtraction_', num2str(round(dRange_eeg))];
end
saveas(fig, fullfile(dir_save, filename), 'fig')
saveas(fig, fullfile(dir_save, filename), 'png')

% plot for tf
for i_fb = 1:nfb
    if baseline_corrected
        dRange_tf = nanmean(prctile(prctile(abs(data_bl_sub_tf{i_fb, 1}), 95), 95));
    else
        dRange_tf = nanmean(prctile(prctile(abs(data_re_tf{i_fb, 1}), 95), 95));
    end
    if range_fixed
        dRange_tf = 30;
    end
    
    fig = figure('DefaultAxesFontSize', 18, 'units', 'normalized', 'outerposition', [0, 0, .67, 1]);
    for i_elec = 1:nelec
        hh(i_elec) = axes( 'position',...
            [ sin( deg2rad(EEG.chanlocs(i_elec).theta) ) * EEG.chanlocs(i_elec).radius * scale + .5 - (figureSize / 2), ...
            cos( deg2rad(EEG.chanlocs(i_elec).theta) ) * EEG.chanlocs(i_elec).radius * scale + .5 - (figureSize / 2), ...
            figureSize, ...
            figureSize ]' );
        
        if baseline_corrected
            imagesc(squeeze( abs(data_bl_sub_tf{i_fb, 1}(i_elec, ind_win_tf(1):ind_win_tf(2), :)) )', [-dRange_tf, dRange_tf]);
        else
            imagesc(squeeze( abs(data_re_tf{i_fb, 1}(i_elec, ind_win_tf(1):ind_win_tf(2), :)) )', [-dRange_tf, dRange_tf]);
        end
        vline(ind_lft_tf - ind_win_tf(1));
        set(gca, 'YDir', 'normal');
        axis off;
    end
    if baseline_corrected
        suptitle([subID, ' ', rg_freq_band{1, i_fb}{1, 1}, ' ', 'w/ baseline corrected, range \pm', num2str(dRange_tf), ' \muV'])
        filename = [subID, '_', rg_freq_band{1, i_fb}{1, 3}, '_w_baseline_subtraction_', num2str(round(dRange_tf))];
    else
        suptitle([subID, ' ', rg_freq_band{1, i_fb}{1, 1}, ' ', 'w/o baseline corrected, range \pm', num2str(dRange_tf), ' \muV'])
        filename = [subID, '_', rg_freq_band{1, i_fb}{1, 3}, '_wo_baseline_subtraction_', num2str(round(dRange_tf))];
    end
    saveas(fig, fullfile(dir_save, filename), 'fig')
    saveas(fig, fullfile(dir_save, filename), 'png')
end
