clear; close all; clc; 
%% turn on EEGLab and select eeg dataset .set file
[ALLEEG, ~, ~, ALLCOM] = eeglab;
[eeg_dataset_file, eeg_dataset_folder] = uigetfile('*.set');

%% load eeg dataset
EEG = pop_loadset('filename', eeg_dataset_file, 'filepath', eeg_dataset_folder);
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
% get EEG montage
E = {EEG.chanlocs.labels}'; 
% convert ANTneuro EEG montage to CSD standarded nomenclature
E{strcmp(E, 'M1')} = 'TP9';
E{strcmp(E, 'M2')} = 'TP10';
% extract montage spherical coordinates for CSD
M = ExtractMontage('10-5-System_Mastoids_EGI129.csd', E);
% plot result
MapMontage(M);

%% Generate Transformation Matrices G and H
[G,H] = GetGH(M);

%% Prepare the Input potentials
D = EEG.data(:, :, 1); % Channels x samples of one epoch

%% Apply the CSD transform
X = CSD(D, G, H);
%%
figure
plot(X')
xlim([0, 200])
figure
plot(D')
xlim([0, 200])
%%
figure
for i = 1:length(X)
    subplot 121
    topoplot(X(:, i), EEG.chanlocs, 'style', 'map', 'electrodes', 'off');%, 'maplimits', 'maxmin');
    title('CSD')
    subplot 122
    topoplot(D(:, i), EEG.chanlocs, 'style', 'map', 'electrodes', 'off');%, 'maplimits', 'maxmin');
    title('potential')
    pause(.1)
end