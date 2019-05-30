close all; clear; clc;
All_dirpath = uigetdir();
All_dirlist = dir(fullfile(All_dirpath, 'sub*'));

All_figpath = fullfile(All_dirpath, 'erp_figs');
if ~exist(All_figpath, 'dir')
    mkdir(All_figpath);
end
% save result?
switch input('Save result figures?(y/N) ', 's')
    case {'y', 'Y'}
        All_is_saving = true;
    otherwise
        All_is_saving = false;
end
% 
disp([num2cell((1:length(All_dirlist))'), {All_dirlist.name}']);
selected_sub = input('Which subject(s) to plot erpimage? ');
% choose which data set to use b4 or after CSD
switch input('Use data after CSD?(y/N) ', 's')
    case {'y', 'Y'}
        All_isCSDapplied = true;
    otherwise
        All_isCSDapplied = false;
end


for All_i = selected_sub% 1:length(All_dirlist)
    clearvars -except All_*; close all;
    filepath = fullfile(All_dirpath, All_dirlist(All_i).name);
    filelist = dir(fullfile(filepath, '*.set'));
    keystr = 'eeg_csd';
    filename = {filelist(contains({filelist.name}, keystr)).name};
    filename = filename{1}; % get the first file named as key string.
    subID = filename(1:6);
    
    % load data set
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
    EEG = pop_loadset('filename', filename, 'filepath', filepath);
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
    
    % EEG.data is the data applied csd
    % EEG.dataRaw is the data b4 appling csd
    tfEEG = EEG;
    if All_isCSDapplied
        tfEEG.data = EEG.dataRaw;
    end
    
    % channel id:
    % {'Fz', 'FCz', 'C3', 'CP3', 'C1'} = {6, 41, 15, 47, 44}
    Fz = 6; FCz = 41; C3 = 15; CP3 = 47; C1 = 44;
    electrodes_name = {'Fz', 'FCz', 'C3', 'CP3', 'C1'};
    electrodes = {Fz; FCz; C3; CP3; C1};
    nb_epoch = length(EEG.epoch);
    tf_ersp = cell(nb_epoch, length(electrodes));
    tf_itc = tf_ersp; tf_powbase = tf_ersp;
    tf_times = tf_ersp; tf_freqs = tf_ersp;
    tf_data = tf_ersp;
    
    % time frequency analysis for each channel and each epoch
    for i = 1:length(electrodes)
        for j = 1:nb_epoch
            figure;
            [tf_ersp{j, i}, tf_itc{j, i}, tf_powbase{j, i}, tf_times{j, i}, tf_freqs{j, i}, ~, ~, tf_data{j, i}] = ...
                       newtimef(tfEEG.data(electrodes{i}, :, j), size(tfEEG.data, 2), [tfEEG.times(1), tfEEG.times(end)], tfEEG.srate, [3, 0.5], ...
                                    'maxfreq', 35, ...
                                    'topovec', 15, 'elocs', EEG.chanlocs, 'chaninfo', EEG.chaninfo, 'caption', electrodes_name{i}, ...
                                    'baseline', [-600, -100], 'basenorm', 'on', 'trialbase', 'full', 'padratio', 1, 'winsize', 512, ...
                                    'plotitc' , 'off', 'plotphase', 'off'); % 'plotersp', 'off'
        end
    end
    
    % the continuous variables: roll angle (deg)
    err = nan(nb_epoch, size(tf_ersp{1, 1}, 2));
    for i = 1:nb_epoch
        err(i, :) = spline(EEG.behavior.obj_epoch{i, 1}.time, EEG.behavior.obj_epoch{i, 1}.roll, tf_times{i, 1});
    end
    % the categorical variables: condition (IL, TR, PT)
    dummy_cond = dummyvar( grp2idx({EEG.epoch.condType}'));
    
    
    % perform regression using power as depend variable; dummy conditions
    % and err as independ variables
    
    % fit = fitlm(ft, 'ERSP~Err*Cond')
    % E(ERSP) = beta_0 + beta_1 * Err + beta_2 * I[IL] + beta_3 * I[TR] + beta_4 * I[PT] + beta_5 * Err * I[IL] + beta_6 * Err * I[TR] + beta_7 * Err * I[PT]
    % each subject results in 8 beta images in the time frequency plane
    
    
    
    
    
    
    
end