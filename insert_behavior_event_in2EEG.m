% Use the EEGLab to open the raw eeg data before running the script

% put behavior lift onset into EEG event
[filename, pathname, ~] = uigetfile;
% load('/Users/yenhsunw/Dropbox (ASU)/NCML-EEGrasp/behavior/matlab data/sandbox/S009_info_onset_time.mat')

EEG = insertEvent2EEG(EEG, pathname, filename);
