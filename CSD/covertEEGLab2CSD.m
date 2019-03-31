clear; close all; clc;
%% Select eeg data set folder
eeg_dataset_folder = uigetdir();
eeg_result_folder = fullfile(fileparts(eeg_dataset_folder), 'eeg_csd');
if ~exist(eeg_result_folder, 'dir')
    mkdir(eeg_result_folder);
end
[project_root, ~, ~] = fileparts(fileparts(fileparts(fileparts(eeg_dataset_folder))));
behavior_folder = fullfile(project_root, 'behavior', 'matlab data', 'preliminary results');
eeg_dataset_list = dir(fullfile(eeg_dataset_folder, '*.set'));

%% turn on EEGLab and select eeg dataset .set file
[ALLEEG, ~, ~, ALLCOM] = eeglab;
% Create G and H for CSD; only need to run one time for the whole project (across subjects and trials)
% load eeg dataset
EEG = pop_loadset('filename', eeg_dataset_list(1).name, 'filepath', eeg_dataset_folder);
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
% get EEG montage
E = {EEG.chanlocs.labels}'; 
% convert ANTneuro EEG montage to CSD standarded nomenclature
E{strcmp(E, 'M1')} = 'TP9';
E{strcmp(E, 'M2')} = 'TP10';
% extract montage spherical coordinates for CSD
M = ExtractMontage('10-5-System_Mastoids_EGI129.csd', E);
% plot electrodes
MapMontage(M);
% Generate Transformation Matrices G and H
[G,H] = GetGH(M);
% clear dataset
STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG = []; CURRENTSET = [];

%% for individual subjects
for sub_i = 1:length(eeg_dataset_list)
    eeg_dataset_file = eeg_dataset_list(sub_i).name;
    %% load behavior info for all trials of one subject
    sub_id = eeg_dataset_file(1:4);
    % load lift onset
    lft_onset = load(fullfile(behavior_folder, [sub_id, '_info_onset_time.mat']));
    tmp = fieldnames(lft_onset);
    lft_onset = lft_onset.(tmp{:}); % lift onset for each epoch in seconds
    % load behavior
    behavior = load(fullfile(behavior_folder, [sub_id, '_for_plot.mat']));
    obj_behavior = cell(length(behavior.obj), 1);
    ind_lft_onset = nan(length(behavior.obj), 1);
    for i = 1:length(behavior.obj)
        tmp_time = array2table(behavior.info_time_trigger{i, 1} - lft_onset(i) * 1000, 'VariableNames', {'time'});
        tmp_time.Properties.VariableUnits = {'ms'};
        tmp_obj = behavior.obj{i};
        tmp_obj.Properties.VariableUnits = {'mm', 'mm', 'mm', 'deg', 'deg', 'deg'};
        obj_behavior{i} = [tmp_time, tmp_obj];
    end
    
    %% load eeg dataset
    EEG = pop_loadset('filename', eeg_dataset_file, 'filepath', eeg_dataset_folder);
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
    obj_epoch = cell(size(obj_behavior));
    for i = 1:length(obj_behavior)
        EEG.behavior.obj_epoch{i, 1} = obj_behavior{i, 1}(obj_behavior{i, 1}.time >= EEG.times(1) & obj_behavior{i, 1}.time <= EEG.times(end), :);
        EEG.behavior.behavior_srate = 1000 / (EEG.behavior.obj_epoch{1, 1}{2, 'time'} - EEG.behavior.obj_epoch{1, 1}{1, 'time'});
        % Prepare the Input potentials for each trial
        D = EEG.data(:, :, i); % Channels x samples of one epoch
        
        % Apply the CSD transform
        EEG.dataCSD(:, :, i) = CSD(D, G, H);
    end
    
    EEG.setname = [sub_id, '_eeg_csd'];
    EEG = eeg_checkset( EEG );
    EEG = pop_saveset( EEG, 'filename', [EEG.setname, '.set'], 'filepath', eeg_result_folder);
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET); 
end

