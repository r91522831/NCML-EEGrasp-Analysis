% EEGLAB history file generated on the 11-Oct-2018
% ------------------------------------------------

EEG.etc.eeglabvers = '14.1.2'; % this tracks which version of EEGLAB is being used, you may ignore it
EEG = pop_loadbv('/Users/yenhsunw/Dropbox (ASU)/NCML-EEGrasp/EEG/Data/S009_EEG/', 'S_J_2018-09-21_11-02-02.vhdr', [1 5246488], [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64]);
EEG = eeg_checkset( EEG );
EEG = pop_saveset( EEG, 'filename','S009_raw.set','filepath','/Users/yenhsunw/Dropbox (ASU)/NCML-EEGrasp/EEG/eeglab/S009/');
EEG = eeg_checkset( EEG );
EEG=pop_chanedit(EEG, 'load',{'/Users/yenhsunw/Dropbox (ASU)/NCML-EEGrasp/EEG/wg64xyz.xyz' 'filetype' 'xyz'});
EEG = eeg_checkset( EEG );
EEG = pop_loadset('filename','S009_channel_loc_lift_onset.set','filepath','/Users/yenhsunw/Dropbox (ASU)/NCML-EEGrasp/EEG/eeglab/S009/');
EEG = eeg_checkset( EEG );
EEG = pop_resample( EEG, 256);
EEG = eeg_checkset( EEG );
EEG = pop_eegfiltnew(EEG, 1,55,846,0,[],0);
EEG.setname=' bandpassed1to55Hz';
EEG = eeg_checkset( EEG );
pop_eegplot( EEG, 1, 1, 1);
pop_prop( EEG, 1, 1, NaN, {'freqrange' [2 50] });
pop_eegplot( EEG, 1, 1, 1);
EEG = pop_reref( EEG, []);
EEG.setname=' re_reference';
EEG = pop_loadset('filename','S009_rereference.set','filepath','/Users/yenhsunw/Dropbox (ASU)/NCML-EEGrasp/EEG/eeglab/S009/');
EEG = eeg_checkset( EEG );
EEG = pop_epoch( EEG, {  'onset'  }, [-2  6], 'newname', ' re_reference epochs', 'epochinfo', 'yes');
EEG = eeg_checkset( EEG );
pop_eegplot( EEG, 1, 1, 1);
EEG = eeg_checkset( EEG );
EEG = pop_runica(EEG, 'extended',1,'interupt','on');
EEG = eeg_checkset( EEG );
EEG = pop_loadset('filename','S009_ICA.set','filepath','/Users/yenhsunw/Dropbox (ASU)/NCML-EEGrasp/EEG/eeglab/S009/');
EEG = eeg_checkset( EEG );
pop_selectcomps(EEG, [1:63] );
EEG = eeg_checkset( EEG );
EEG = eeg_checkset( EEG );
pop_selectcomps(EEG, [1:63] );
EEG = eeg_checkset( EEG );
EEG = pop_subcomp( EEG, [1   2   3   5  25  27  31  33  35  43  46  59], 0);
EEG.setname='epochs pruned with ICA';
EEG = pop_loadset('filename','S009_pruned_ICA.set','filepath','/Users/yenhsunw/Dropbox (ASU)/NCML-EEGrasp/EEG/eeglab/S009/');
EEG = eeg_checkset( EEG );

%% Get experiment conditions
sub_dir = uigetdir('', 'select subject folder for behvior matlab raw data.');
file_list = dir(fullfile(sub_dir, '*.csv'));
tmp = char({file_list.name});
cond = tmp(:, [11:12, 16:18]);

%% channel id:
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
    figure; [tf_ersp{i}, tf_itc{i}, tf_powbase{i}, tf_times{i}, tf_freqs{i}, ~, ~, tf_data{i}] = pop_newtimef( EEG, 1, electrodes{i}, [-2000, 5996], [3, 0.5] , 'topovec', 15, 'elocs', EEG.chanlocs, 'chaninfo', EEG.chaninfo, 'caption', electrodes_name{i}, 'baseline',[-600, -100], 'basenorm', 'on', 'trialbase', 'full', 'plotitc' , 'off', 'plotphase', 'off', 'padratio', 1, 'trialbase', 'full', 'winsize', 512);
end

%% sort the tf_data according to conditions
for i = 1:length(electrodes)
    tf_data{i, 2} = tf_data{i, 1}(:, :, cond(:, 1) == 'I');
    tf_data{i, 3} = tf_data{i, 1}(:, :, cond(:, 1) == 'T');
    tf_data{i, 4} = tf_data{i, 1}(:, :, and(cond(:, 1) == 'P', cond(:, end) == '1'));
    tf_data{i, 5} = tf_data{i, 1}(:, :, and(cond(:, 1) == 'P', cond(:, end) == '2'));
    tf_data{i, 6} = tf_data{i, 1}(:, :, and(cond(:, 1) == 'P', cond(:, end) == '3'));
end

%% Get 150-350 ms for 4-8 Hz (theta); 450-650 ms for 9-12 Hz (alpha); 450-650 ms for 20-30 Hz (beta)
f_range = and(tf_freqs{i} >= 4, tf_freqs{i} <= 8);
t_range = and(tf_times{i} >= 150, tf_times{i} <= 350);
tf4theta{1, 1} = tf_data{1, 2}(f_range, t_range);