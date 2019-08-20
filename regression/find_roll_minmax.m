close all; clear; clc;
All_dirpath = uigetdir();
All_dirlist = dir(fullfile(All_dirpath, 'sub*'));

All_linearmodel_path = fullfile(All_dirpath, 'linear');
if ~exist(All_linearmodel_path, 'dir')
    mkdir(All_linearmodel_path);
end

% 
disp([num2cell((1:length(All_dirlist))'), {All_dirlist.name}']);
selected_sub = input('Which subject(s) to plot erpimage? ');
if isempty(selected_sub)
    selected_sub = 1:length(All_dirlist);
end

load(fullfile(All_dirpath, 'linear', 'misc.mat'));
All_timestamp = squeeze(tf_times{1, 1}(:, :, 1));

All_roll = nan(length(All_dirlist), length(All_timestamp), 2);
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

    nb_epoch = length(EEG.epoch);
    roll_ang = nan(nb_epoch, length(All_timestamp));
    for i = 1:nb_epoch
        roll_ang(i, :) = spline(EEG.behavior.obj_epoch{i, 1}.time, EEG.behavior.obj_epoch{i, 1}.roll, All_timestamp);
    end
    [All_roll(All_i, :, 1), All_roll(All_i, :, 2)] = bounds(roll_ang);
end


err_range = [max(All_roll(:, :, 1), [], 1); min(All_roll(:, :, 2), [], 1)];

save(fullfile(All_dirpath, 'err_range.mat'), 'err_range');
%%


