function newEEG = applyCSD2EEGset(EEG)
%% Create G and H for CSD
% use CSD Toolbox
% ===========
% © 2003-2010 by Jürgen Kayser
% Version 1.1 (July 21, 2010)
% http://psychophysiology.cpmc.columbia.edu/Software/CSDtoolbox/

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
EEG.dataRaw = EEG.data; % keep the original data
for i = 1:EEG.trials
    % Prepare the Input potentials for each trial
    D = EEG.data(:, :, i); % Channels x samples of one epoch
    % Apply the CSD transform
    EEG.data(:, :, i) = CSD(D, G, H);
end

EEG = eeg_checkset( EEG );
newEEG = EEG;
end