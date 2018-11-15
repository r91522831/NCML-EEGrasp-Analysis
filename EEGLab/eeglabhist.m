clear; close all; clc;

% EEGLAB history file generated on the 11-Oct-2018
% ------------------------------------------------
% '/Users/yenhsunw/Dropbox (ASU)/NCML-EEGrasp/EEG/Data/Sxxx_EEG/', '*.vhdr'
[raw_filename, raw_dir, ~] = uigetfile('*.vhdr','Select the raw EEG data file');
[tmp_dir, foldername, ~]= fileparts(fileparts(raw_dir));
sub_id = foldername(1:4);
% '/Users/yenhsunw/Dropbox (ASU)/NCML-EEGrasp/EEG/eeglab/Sxxx/';
[tmp_dir, ~, ~] = fileparts(tmp_dir);
output_dir = fullfile(tmp_dir, 'eeglab', 'archive', sub_id);

%% EEGLab
EEG.etc.eeglabvers = '14.1.2'; % this tracks which version of EEGLAB is being used, you may ignore it
EEG = pop_loadbv(raw_dir, raw_filename);
EEG.setname = [sub_id, '_raw'];
EEG = eeg_checkset( EEG );
EEG = pop_saveset( EEG, 'filename', [EEG.setname, '.set'], 'filepath', output_dir);
EEG = pop_chanedit(EEG, 'load',{'/Users/yenhsunw/Dropbox (ASU)/NCML-EEGrasp/EEG/wg64xyz.xyz' 'filetype' 'xyz'},'settype',{'1:63' 'EEG'},'settype',{'64' 'EOG'});
EEG = eeg_checkset( EEG );
% run insert_behavior_event_in2EEG to put behavior onset into EEG events
behavior_dir = '/Users/yenhsunw/Dropbox (ASU)/NCML-EEGrasp/behavior/matlab data/preliminary results';
behavior_filename = [sub_id, '_info_onset_time.mat'];
EEG = insertEvent2EEG(EEG, behavior_dir, behavior_filename);
EEG = eeg_checkset( EEG );
EEG.setname = [sub_id, '_channel_loc_lift_onset'];
EEG = eeg_checkset( EEG );
EEG = pop_saveset( EEG, 'filename', [EEG.setname, '.set'], 'filepath', output_dir);
% highpass filter at 0.5 Hz
EEG = pop_eegfiltnew(EEG, 'locutoff', 0.5);
EEG = eeg_checkset( EEG );
EEG.setname = [sub_id, '_highpass_halfHz'];
EEG = eeg_checkset( EEG );
EEG = pop_saveset( EEG, 'filename', [EEG.setname, '.set'], 'filepath', output_dir);
% resample at 256 Hz
EEG = pop_resample( EEG, 256);
EEG = eeg_checkset( EEG );
EEG.setname = [sub_id, '_resampled256Hz'];
EEG = eeg_checkset( EEG );
EEG = pop_saveset( EEG, 'filename', [EEG.setname, '.set'], 'filepath', output_dir);
% lowpass filter at 100 Hz
EEG = pop_eegfiltnew(EEG, 'hicutoff', 100);
EEG = eeg_checkset( EEG );
EEG.setname = [sub_id, '_lowpass100Hz'];
EEG = eeg_checkset( EEG );
EEG = pop_saveset( EEG, 'filename', [EEG.setname, '.set'], 'filepath', output_dir);
% rereference at average electrod
EEG = pop_reref( EEG, [],'exclude',64);
EEG = eeg_checkset( EEG );
EEG.setname = [sub_id, '_rereference'];
EEG = eeg_checkset( EEG );
EEG = pop_saveset( EEG, 'filename', [EEG.setname, '.set'], 'filepath', output_dir);
% epoching
EEG = pop_epoch( EEG, {  'onset'  }, [-1.5  2.5], 'newname', [sub_id, '_epochs'], 'epochinfo', 'yes');
EEG = eeg_checkset( EEG );
EEG = pop_saveset( EEG, 'filename', [EEG.setname, '.set'], 'filepath', output_dir);
% run ICA
EEG = pop_runica(EEG, 'chanind', 1:63, 'extended', 1, 'interupt', 'on');
EEG = eeg_checkset( EEG );
EEG.setname = [sub_id, '_ICA'];
EEG = eeg_checkset( EEG );
EEG = pop_saveset( EEG, 'filename', [EEG.setname, '.set'], 'filepath', output_dir);
% remove EOG channel for MARA ICA
EEG = pop_select( EEG,'channel',{'FP1' 'FPZ' 'FP2' 'F7' 'F3' 'FZ' 'F4' 'F8' 'FC5' 'FC1' 'FC2' 'FC6' 'M1' 'T7' 'C3' 'CZ' 'C4' 'T8' 'M2' 'CP5' 'CP1' 'CP2' 'CP6' 'P7' 'P3' 'PZ' 'P4' 'P8' 'POZ' 'O1' 'O2' 'AF7' 'AF3' 'AF4' 'AF8' 'F5' 'F1' 'F2' 'F6' 'FC3' 'FCZ' 'FC4' 'C5' 'C1' 'C2' 'C6' 'CP3' 'CP4' 'P5' 'P1' 'P2' 'P6' 'PO5' 'PO3' 'PO4' 'PO6' 'FT7' 'FT8' 'TP7' 'TP8' 'PO7' 'PO8' 'OZ'});
EEG = eeg_checkset( EEG );
EEG.setname = [sub_id, '_remove_EOG4MARA'];
EEG = eeg_checkset( EEG );
EEG = pop_saveset( EEG, 'filename', [EEG.setname, '.set'], 'filepath', output_dir);

%% Naked eye ICA inspection

% run SASICA

% run ICALabel

% run MARA

%% Manually select ICA component reject 
% load back EOG channel
EEG = pop_loadset('filename', [sub_id, '_ICA.set'], 'filepath', output_dir);
EEG = eeg_checkset( EEG );
% artifact components table
S002_ac = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 25, 26, 27, 28, 31, 33, 34, 35, 36, 37, 38, 41, 42, 43, 44, 47, 48, 51, 55, 57];
S003_ac = [2, 3, 5, 6, 7, 8, 9, 10, 15, 17, 18, 20, 21, 22, 24, 28, 31, 36, 37, 38, 39, 41, 42, 44, 45, 46, 47, 48, 49, 51, 52, 53, 55, 57, 61];
S004_ac = [1, 2, 3, 4, 7, 13, 15, 17, 18, 19, 20, 21, 22, 23, 25, 26, 27  28, 29, 30, 31, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 48, 49, 50, 53, 54, 56, 57];
S005_ac = [1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 14, 17, 18, 22, 25, 27, 28, 29, 30, 31, 32, 33, 35, 36, 37, 38, 39, 40, 43, 45, 47, 48, 49, 50, 51, 52, 53, 55, 57, 59, 60, 62, 63];
S006_ac = [1, 2, 3, 4, 5, 6, 8, 11, 12, 13, 14, 15, 17, 18, 19, 20, 21  22, 23, 24, 25, 27, 28, 29, 30, 32, 33, 34, 36, 37, 38, 39, 40, 41, 42, 43, 44, 46, 47, 48, 49, 53, 54, 55, 60, 62, 63];
S009_ac = [1, 2, 3, 4, 8, 9, 13, 15, 17, 19, 21, 23, 24, 26, 27, 28, 30, 31, 34, 35, 36, 38, 39, 40, 41, 42, 43, 47, 48, 49, 51, 53, 55, 57, 59, 61, 62];
ac_table = cell2table({S002_ac; S003_ac; S004_ac; S005_ac; S006_ac; S009_ac}, 'VariableNames', {'artifact_component'}, 'RowNames', {'S002', 'S003', 'S004', 'S005', 'S006', 'S009'});
% remove artifact components
EEG = pop_subcomp( EEG, ac_table{sub_id, :}{:}, 0);
EEG = eeg_checkset( EEG );
% % % EEG = pop_selectcomps(EEG, 1:63 );
% % % EEG = eeg_checkset( EEG );
EEG.setname = [sub_id, '_pruned_ICA'];
EEG = eeg_checkset( EEG );
EEG = pop_saveset( EEG, 'filename', [EEG.setname, '.set'], 'filepath', output_dir);








%%
EEG = pop_loadset('filename', [sub_id, '_pruned_ICA.set'], 'filepath', output_dir);

% Get experiment conditions
cond_names = {'ALL', 'IL', 'TR', 'PT1', 'PT2', 'PT3'};
% /Users/yenhsunw/Dropbox (ASU)/NCML-EEGrasp/behavior/matlab data/Sxxx
sub_dir = fullfile('/Users/yenhsunw/Dropbox (ASU)/NCML-EEGrasp/behavior/matlab data/', sub_id);
file_list = dir(fullfile(sub_dir, '*.csv'));
tmp = char({file_list.name});
cond = tmp(:, [11:12, 16:18]);

% Get peak roll for each trial
% load peak_roll
load(fullfile('/Users/yenhsunw/Dropbox (ASU)/NCML-EEGrasp/behavior/matlab data/preliminary results/', [sub_id, '_temp_result.mat']));
% sort peak roll according to conditions
pRoll{1, 1} = peak_roll(:, 1);
pRoll{1, 2} = peak_roll(cond(:, 1) == 'I', 1);
pRoll{1, 3} = peak_roll(cond(:, 1) == 'T', 1);
pRoll{1, 4} = peak_roll(and(cond(:, 1) == 'P', cond(:, end) == '1'), 1);
pRoll{1, 5} = peak_roll(and(cond(:, 1) == 'P', cond(:, end) == '2'), 1);
pRoll{1, 6} = peak_roll(and(cond(:, 1) == 'P', cond(:, end) == '3'), 1);
pRoll = cell2table(pRoll, 'VariableNames', cond_names);

% channel id:
% {'Fz', 'FCz', 'C3', 'CP3', 'C1'} = {6, 41, 15, 47, 44}
Fz = 6; FCz = 41; C3 = 15; CP3 = 47; C1 = 44;
electrodes_name = {'Fz', 'FCz', 'C3', 'CP3', 'C1'};
electrodes = {Fz; FCz; C3; CP3; C1};
tf_ersp = cell(size(electrodes));
tf_itc = cell(size(electrodes));
tf_powbase = cell(size(electrodes));
tf_times = cell(size(electrodes));
tf_freqs = cell(size(electrodes));
tf_data = cell(size(electrodes));
for i = 1:length(electrodes)
    figure; [tf_ersp{i}, tf_itc{i}, tf_powbase{i}, tf_times{i}, tf_freqs{i}, ~, ~, tf_data{i}] = pop_newtimef( EEG, 1, electrodes{i}, [-1500, 2500], [3, 0.5] , 'topovec', 15, 'elocs', EEG.chanlocs, 'chaninfo', EEG.chaninfo, 'caption', electrodes_name{i}, 'baseline',[-600, -100], 'basenorm', 'on', 'trialbase', 'full', 'plotitc' , 'off', 'plotphase', 'off', 'padratio', 1, 'trialbase', 'full', 'winsize', 512);
end

% sort the tf_data according to conditions
for i = 1:length(electrodes)
    tf_data{i, 2} = tf_data{i, 1}(:, :, cond(:, 1) == 'I');
    tf_data{i, 3} = tf_data{i, 1}(:, :, cond(:, 1) == 'T');
    tf_data{i, 4} = tf_data{i, 1}(:, :, and(cond(:, 1) == 'P', cond(:, end) == '1'));
    tf_data{i, 5} = tf_data{i, 1}(:, :, and(cond(:, 1) == 'P', cond(:, end) == '2'));
    tf_data{i, 6} = tf_data{i, 1}(:, :, and(cond(:, 1) == 'P', cond(:, end) == '3'));
end

tf_data = cell2table(tf_data, 'VariableNames', cond_names, 'RowNames', electrodes_name);

% Get 150-350 ms for 4-8 Hz (theta); 450-650 ms for 9-12 Hz (alpha); 450-650 ms for 20-30 Hz (beta) for individual trial
t_tf4theta = cell(size(tf_data));
t_tf4alpha = cell(size(tf_data));
t_tf4beta = cell(size(tf_data));
for i = 1:length(electrodes)
    for j = 1:length(cond_names)
        f_range = and(tf_freqs{i} >= 4, tf_freqs{i} <= 8);
        t_range = and(tf_times{i} >= 150, tf_times{i} <= 350);
        for k = 1:length(tf_data{i, j}{:}(1, 1, :))
            tmp_avg = squeeze(mean(mean(tf_data{i, j}{:}(f_range, t_range, k), 1), 2));
            t_tf4theta{i, j}(k, 1) = tmp_avg;
        end
        f_range = and(tf_freqs{i} >= 9, tf_freqs{i} <= 12);
        t_range = and(tf_times{i} >= 450, tf_times{i} <= 650);
        for k = 1:length(tf_data{i, j}{:}(1, 1, :))
            tmp_avg = squeeze(mean(mean(tf_data{i, j}{:}(f_range, t_range, k), 1), 2));
            t_tf4alpha{i, j}(k, 1) = tmp_avg;
        end
        f_range = and(tf_freqs{i} >= 20, tf_freqs{i} <= 30);
        t_range = and(tf_times{i} >= 450, tf_times{i} <= 650);
        for k = 1:length(tf_data{i, j}{:}(1, 1, :))
            tmp_avg = squeeze(mean(mean(tf_data{i, j}{:}(f_range, t_range, k), 1), 2));
            t_tf4beta{i, j}(k, 1) = tmp_avg;
        end
    end
end
t_tf4theta = cell2table(t_tf4theta, 'VariableNames', cond_names, 'RowNames', electrodes_name);
t_tf4alpha = cell2table(t_tf4alpha, 'VariableNames', cond_names, 'RowNames', electrodes_name);
t_tf4beta = cell2table(t_tf4beta, 'VariableNames', cond_names, 'RowNames', electrodes_name);

% compute average for conditions, channels, and brainwave bands for each trial
t_tf4theta_Fz_FCz = cell(1, length(cond_names));
t_tf4alpha_C3_CP3 = cell(1, length(cond_names));
t_tf4beta_C3_C1 = cell(1, length(cond_names));
for j = 1:length(cond_names)
    t_tf4theta_Fz_FCz{1, j} = mean([t_tf4theta{{'Fz'}, j}{:}, t_tf4theta{{'FCz'}, j}{:}], 2);    
    t_tf4alpha_C3_CP3{1, j} = mean([t_tf4alpha{{'C3'}, j}{:}, t_tf4alpha{{'CP3'}, j}{:}], 2);
    t_tf4beta_C3_C1{1, j} = mean([t_tf4beta{{'C3'}, j}{:}, t_tf4beta{{'C1'}, j}{:}], 2);
end
cellfun(@abs, t_tf4theta_Fz_FCz, 'UniformOutput', false);

t_tf_power = [cellfun(@abs, t_tf4theta_Fz_FCz, 'UniformOutput', false); cellfun(@abs, t_tf4alpha_C3_CP3, 'UniformOutput', false); cellfun(@abs, t_tf4beta_C3_C1, 'UniformOutput', false)];
t_tf_power = array2table(t_tf_power, 'VariableNames', cond_names, 'RowNames', {'theta', 'alpha', 'beta'});
t_tf_power_complex = [t_tf4theta_Fz_FCz; t_tf4alpha_C3_CP3; t_tf4beta_C3_C1];
t_tf_power_complex = array2table(t_tf_power_complex, 'VariableNames', cond_names, 'RowNames', {'theta', 'alpha', 'beta'});

% Get 150-350 ms for 4-8 Hz (theta); 450-650 ms for 9-12 Hz (alpha); 450-650 ms for 20-30 Hz (beta)
tf4theta = cell(size(tf_data));
tf4alpha = cell(size(tf_data));
tf4beta = cell(size(tf_data));
for i = 1:length(electrodes)
    for j = 1:length(cond_names)
        f_range = and(tf_freqs{i} >= 4, tf_freqs{i} <= 8);
        t_range = and(tf_times{i} >= 150, tf_times{i} <= 350);
        tmp_avg = squeeze(mean(mean(mean(tf_data{i, j}{:}(f_range, t_range, :), 1), 2), 3));
        tf4theta{i, j} = tmp_avg;
        f_range = and(tf_freqs{i} >= 9, tf_freqs{i} <= 12);
        t_range = and(tf_times{i} >= 450, tf_times{i} <= 650);
        tmp_avg = squeeze(mean(mean(mean(tf_data{i, j}{:}(f_range, t_range, :), 1), 2), 3));
        tf4alpha{i, j} = tmp_avg;
        f_range = and(tf_freqs{i} >= 20, tf_freqs{i} <= 30);
        t_range = and(tf_times{i} >= 450, tf_times{i} <= 650);
        tmp_avg = squeeze(mean(mean(mean(tf_data{i, j}{:}(f_range, t_range, :), 1), 2), 3));
        tf4beta{i, j} = tmp_avg;
    end
end
tf4theta = cell2table(tf4theta, 'VariableNames', cond_names, 'RowNames', electrodes_name);
tf4alpha = cell2table(tf4alpha, 'VariableNames', cond_names, 'RowNames', electrodes_name);
tf4beta = cell2table(tf4beta, 'VariableNames', cond_names, 'RowNames', electrodes_name);

% compute average for conditions, channels, and brainwave bands
tf4theta_Fz_FCz = mean(tf4theta{{'Fz', 'FCz'}, :}, 1);
tf4alpha_C3_CP3 = mean(tf4alpha{{'C3', 'CP3'}, :}, 1);
tf4beta_C3_C1 = mean(tf4beta{{'C3', 'C1'}, :}, 1);
tf_power = [abs(tf4theta_Fz_FCz); abs(tf4alpha_C3_CP3); abs(tf4beta_C3_C1)];
tf_power = array2table(tf_power, 'VariableNames', cond_names, 'RowNames', {'theta', 'alpha', 'beta'});
tf_power_complex = [tf4theta_Fz_FCz; tf4alpha_C3_CP3; tf4beta_C3_C1];
tf_power_complex = array2table(tf_power_complex, 'VariableNames', cond_names, 'RowNames', {'theta', 'alpha', 'beta'});

%'/Users/yenhsunw/Dropbox (ASU)/NCML-EEGrasp/EEG/eeglab/power_chage'
[tmp_dir, ~, ~] = fileparts(output_dir);
tmp_filename = [sub_id, '_tf_power.mat'];
power_change_filename = fullfile(tmp_dir, 'power_change', tmp_filename);
save(power_change_filename, 'tf_power', 'tf_power_complex', 't_tf_power', 't_tf_power_complex', 'pRoll');

