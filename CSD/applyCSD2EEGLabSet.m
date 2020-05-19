

%% Create G and H for CSD
% use CSD Toolbox
% ===========
% © 2003-2010 by Jürgen Kayser
% Version 1.1 (July 21, 2010)
% http://psychophysiology.cpmc.columbia.edu/Software/CSDtoolbox/

% load eeg dataset
EEG = pop_loadset('filename', tmp_dataset_list(1).name, 'filepath', tmp_dataset_folder);
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
%% Apply CSD and save result
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