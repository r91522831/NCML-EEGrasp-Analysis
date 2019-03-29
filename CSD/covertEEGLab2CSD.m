clear; close all; clc; 
%% turn on EEGLab and select eeg dataset .set file
[ALLEEG, ~, ~, ALLCOM] = eeglab;
[eeg_dataset_file, eeg_dataset_folder] = uigetfile('*.set');

%% Create G and H for CSD; only need to run one time for the whole project (across subjects and trials)
% load eeg dataset
EEG = pop_loadset('filename', eeg_dataset_file, 'filepath', eeg_dataset_folder);
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
%% Generate Transformation Matrices G and H
[G,H] = GetGH(M);
% clear dataset
STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG = []; CURRENTSET = [];






%% load behavior info for all trials of one subject
[project_root, ~, ~] = fileparts(fileparts(fileparts(fileparts(fileparts(eeg_dataset_folder)))));
behavior_folder = fullfile(project_root, 'behavior', 'matlab data', 'preliminary results');

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
X = cell(size(obj_behavior));
for i = 1:length(obj_behavior)
    obj_epoch{i, 1} = obj_behavior{i, 1}(obj_behavior{i, 1}.time >= EEG.times(1) & obj_behavior{i, 1}.time <= EEG.times(end), :);
    
    % Prepare the Input potentials for each trial
    D = EEG.data(:, :, i); % Channels x samples of one epoch
    
    % Apply the CSD transform
    X{i, 1} = CSD(D, G, H);
end






%%
%{
figure
for i = 1:length(X)
    subplot 221
    topoplot(X(:, i), EEG.chanlocs, 'style', 'map', 'electrodes', 'off');%, 'maplimits', 'maxmin');
    title('CSD')
    subplot 222
    topoplot(D(:, i), EEG.chanlocs, 'style', 'map', 'electrodes', 'off');%, 'maplimits', 'maxmin');
    title('potential')
    subplot 223
    
    
    pause(.1)
end
%}