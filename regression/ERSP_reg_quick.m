close all; clear; clc;
All_path = uigetdir;
All_filelist = dir(fullfile(All_path, '*_timefreq.mat'));

%% 
disp([num2cell((1:length(All_filelist))'), {All_filelist.name}']);
selected_sub = input('Which subject(s) to run regression? ');
if isempty(selected_sub)
    selected_sub = 1:length(All_filelist);
end

nelectrode = 1;
coeff_name = {'\beta_0', '\beta_1', '\beta_2', '\beta_3', '\beta_4', '\beta_5'};
ncoeff = length(coeff_name);


for All_i = selected_sub% 1:length(All_dirlist)
    clearvars -except All_*; close all;
    filepath = fullfile(All_path, All_filelist(All_i).name);
    filelist = dir(fullfile(filepath, '*.set'));
    keystr = 'eeg_csd';
    filename = {filelist(contains({filelist.name}, keystr)).name};
    filename = filename{1}; % get the first file named as key string.
    subID = filename(1:6);
    
    % load data set
% % %     tf = load('/Users/yenhsunw/Dropbox (ASU)/NCML-EEGrasp/EEG/eeglab/002 ProcessAgain/linear/RAW/sub-09_timefreq.mat');
    tf = load(fullfile(All_path, All_filelist(All_i).name));

    timerstamps = tf.tf_times';
    freqz = tf.tf_freqs';
    
    ntime = length(timerstamps);
    nfreq = length(freqz);
    
    
    % delta: 2 ~ 4 Hz; theta: 5 ~ 8 Hz; alpha: 9 ~ 12 Hz; low beta: 13 ~ 20 Hz; high beta: 21 ~ 30 Hz; low gama: 31 ~ 35 Hz
    rg_freq_band = { {'\delta', '2-4 Hz'}, {'\theta', '5-8 Hz'}, {'\alpha', '9-12 Hz'}, {'\beta_{low}', '13-20 Hz'}, {'\beta_{high}', '21-30 Hz'}, {'\gamma_{low}', '31-35 Hz'}; ...
        find(freqz >= 2 & freqz <=  4), find(freqz >  4 & freqz <=  8), ...
        find(freqz >  8   & freqz <= 13), find(freqz > 13 & freqz <= 20), ...
        find(freqz > 20   & freqz <= 30), find(freqz > 30 & freqz <= 35) };
    rg_time_win = { {'-50 to 150 ms'}, {'150 to 350 ms'}, {'350 to 450 ms'}, {'450 to 650 ms'}, {'650 to 850 ms'}; ...
        find(timerstamps >= -50 & timerstamps < 150), find(timerstamps >= 150 & timerstamps < 350), ...
        find(timerstamps >= 350 & timerstamps < 450), find(timerstamps >= 450 & timerstamps < 650), ...
        find(timerstamps >= 650 & timerstamps < 850) };

            
            
end            