clear; close all; clc;
%% Select eeg data set folder
eeg_dataset_folder = uigetdir();
eeg_dataset_list = dir(fullfile(eeg_dataset_folder, 'sub-*'));

%% turn on EEGLab and select eeg dataset .set file
[ALLEEG, ~, ~, ALLCOM] = eeglab;

disp([num2cell((1:length(eeg_dataset_list))'), {eeg_dataset_list.name}']);
selected_sub = input('Which subject(s) to run CSD? ');
if isempty(selected_sub)
    selected_sub = 1:length(eeg_dataset_list);
end

%% for individual subjects
for sub_i = selected_sub%1:length(eeg_dataset_list)
    tmp_dataset_folder = fullfile(eeg_dataset_folder, eeg_dataset_list(sub_i).name);
    tmp_dataset_list = dir(fullfile(tmp_dataset_folder, '*_combine_behavior.set')); % Should be only one file
    if ~exist(fullfile(tmp_dataset_folder, tmp_dataset_list(1).name), 'file')
        disp('EEG data set name is not correct!')
    end
    
    if sub_i == selected_sub(1)
        %% Create G and H for CSD; only need to run one time for the whole project (across subjects and trials)
        % use CSD Toolbox
        % ===========
        % © 2003-2010 by Jürgen Kayser
        % Version 1.1 (July 21, 2010)
        % http://psychophysiology.cpmc.columbia.edu/Software/CSDtoolbox/
        
        % load eeg dataset
        EEG = pop_loadset('filename', tmp_dataset_list(1).name, 'filepath', tmp_dataset_folder);
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
    end
    
    %% load eeg dataset
    EEG = pop_loadset('filename', tmp_dataset_list(1).name, 'filepath', tmp_dataset_folder);
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
    for i = 1:size(EEG.data, 3)
        % Prepare the Input potentials for each trial
        D = EEG.data(:, :, i); % Channels x samples of one epoch
        
        % Apply the CSD transform
        EEG.dataRaw = EEG.data;
        EEG.data(:, :, i) = CSD(D, G, H);
    end
    sub_id = tmp_dataset_list(1).name(1:6);
    EEG.setname = [sub_id, '_eeg_csd'];
    EEG = eeg_checkset( EEG );
    EEG = pop_saveset( EEG, 'filename', [EEG.setname, '.set'], 'filepath', tmp_dataset_folder);
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET); 
end

