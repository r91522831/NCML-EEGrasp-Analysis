clear; close all; clc;
%% Select eeg data set folder
eeg_dataset_folder = uigetdir();
eeg_dataset_list = dir(fullfile(eeg_dataset_folder, 'sub-*'));
[project_root, ~, ~] = fileparts(fileparts(fileparts(eeg_dataset_folder)));
behavior_folder = fullfile(project_root, 'behavior', 'preliminary results');
if ~exist(behavior_folder, 'dir')
    disp('Behavior folder is not correct!!')
end

%% turn on EEGLab and select eeg dataset .set file
[ALLEEG, ~, ~, ALLCOM] = eeglab;

disp([num2cell((1:length(eeg_dataset_list))'), {eeg_dataset_list.name}']);
selected_sub = input('Which subject(s) to combine behavior to EEG? ');
if isempty(selected_sub)
    selected_sub = 1:length(eeg_dataset_list);
end

%% for individual subjects
for sub_i = selected_sub%1:length(eeg_dataset_list)
    tmp_dataset_folder = fullfile(eeg_dataset_folder, eeg_dataset_list(sub_i).name);
    tmp_dataset_list = dir(fullfile(tmp_dataset_folder, '*_pruned_ICA.set')); % Should be only one file
    if ~exist(fullfile(tmp_dataset_folder, tmp_dataset_list(1).name), 'file')
        disp('EEG data set name is not correct!')
    end
    
    %% load behavior info for all trials of one subject
    sub_id = tmp_dataset_list(1).name(1:6);
    beh_sub_id = ['S0', tmp_dataset_list(1).name(5:6)];
    % load lift onset
    lft_onset = load(fullfile(behavior_folder, [beh_sub_id, '_info_onset_time.mat']));
    tmp = fieldnames(lft_onset);
    lft_onset = lft_onset.(tmp{:}); % lift onset for each epoch in seconds
    % load behavior
    behavior = load(fullfile(behavior_folder, [beh_sub_id, '_for_plot.mat']));
    obj_behavior = cell(length(behavior.obj), 1);
    ind_lft_onset = nan(length(behavior.obj), 1);
    for i = 1:length(behavior.obj)
        tmp_time = array2table(behavior.info_time_trigger{i, 1} - lft_onset(i) * 1000, 'VariableNames', {'time'});
        tmp_time.Properties.VariableUnits = {'ms'};
        tmp_obj = behavior.obj{i};
        tmp_obj.Properties.VariableUnits = {'mm', 'mm', 'mm', 'rad', 'rad', 'rad'};
        obj_behavior{i} = [tmp_time, tmp_obj];
    end
    
    %% load eeg dataset
    EEG = pop_loadset('filename', tmp_dataset_list(1).name, 'filepath', tmp_dataset_folder);
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
    obj_epoch = cell(size(obj_behavior));
    for i = 1:length(obj_behavior)
        EEG.behavior.obj_epoch{i, 1} = obj_behavior{i, 1}(obj_behavior{i, 1}.time >= EEG.times(1) & obj_behavior{i, 1}.time <= EEG.times(end), :);
        EEG.behavior.behavior_srate = 1000 / (EEG.behavior.obj_epoch{1, 1}{2, 'time'} - EEG.behavior.obj_epoch{1, 1}{1, 'time'});
    end
    
    EEG.setname = [sub_id, '_combine_behavior'];
    EEG = eeg_checkset( EEG );
    EEG = pop_saveset( EEG, 'filename', [EEG.setname, '.set'], 'filepath', tmp_dataset_folder);
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET); 
end

