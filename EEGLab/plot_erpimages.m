close all; clear; clc;
filepath = uigetdir();

filelist = dir(fullfile(filepath, '*.set'));
keystr = 'pruned_ICA';
filename = {filelist(contains({filelist.name}, keystr)).name};
filename = filename{1}; % get the first file named as key string.

%%
% load data set
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
EEG = pop_loadset('filename', filename, 'filepath', filepath);
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );

%%
% plot erpimage w/o filtering


%%
% plot erpimage at theta bands 4 - 8 Hz
EEG = pop_eegfiltnew(EEG, 'locutoff', 4, 'hicutoff', 8);

[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1, 'setname', [EEG.setname, '_theta_4to8Hz'], 'savenew', fullfile(filepath, [EEG.setname, '_theta_4to8Hz.set']), 'gui', 'off'); 


freqName = '\theta';
colID = 3; % column 1, 2, or 3
figure(2);
chanName = 'Cz'; % 'Fz' FCz'
erpimage_subplot(EEG, freqName, chanName, colID);
colID = 2; % column 1, 2, or 3
chanName = 'Fz'; % 'Fz' FCz'
erpimage_subplot(EEG, freqName, chanName, colID);
colID = 1; % column 1, 2, or 3
chanName = 'FCz'; % 'Fz' FCz'
erpimage_subplot(EEG, freqName, chanName, colID);

%{
figure(3);
chanName = 'FCz'; % 'Fz' FCz'
erpimage_subplot(EEG, freqName, chanName, colID);

% plot erpimage at alpha bands 9 - 12 Hz
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2, 'retrieve', 1, 'study', 0); 
EEG = pop_eegfiltnew(EEG, 'locutoff', 9, 'hicutoff', 12);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1, 'setname', [EEG.setname, '_alpha_9to12Hz'], 'savenew', fullfile(filepath, [EEG.setname, '_alpha_9to12Hz.set']), 'gui', 'off'); 


freqName = '\alpha';
colID = 1; % 1, 2, or 3
figure(2);
chanName = 'C3'; % 'C3' 'CP3'
erpimage_subplot(EEG, freqName, chanName, colID);
figure(3);
chanName = 'CP3'; % 'C3' 'CP3'
erpimage_subplot(EEG, freqName, chanName, colID);

% plot erpimage at beta bands 20 - 30 Hz
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 3, 'retrieve', 1, 'study', 0); 
EEG = pop_eegfiltnew(EEG, 'locutoff', 20, 'hicutoff', 30);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1, 'setname', [EEG.setname, '_beta_20to30Hz'], 'savenew', fullfile(filepath, [EEG.setname, '_beta_20to30Hz.set']), 'gui', 'off'); 

freqName = '\beta';
colID = 2; % 1, 2, or 3
figure(2);
chanName = 'C1'; %'C1' 'C3'
erpimage_subplot(EEG, freqName, chanName, colID);
figure(3);
chanName = 'C3'; %'C1' 'C3'
erpimage_subplot(EEG, freqName, chanName, colID);
%}