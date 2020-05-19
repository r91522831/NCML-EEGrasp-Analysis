% prepare data for HERMES
% raw data
EEG = pop_loadset('filename','sub-09_SouceLocalized.set','filepath','/Users/yen-hsunwu/Dropbox (ASU)/BIDS_format/NCML-EEGrasp/sub-09/eeg/set/');
e_type = cellfun(@(x) x(1), {EEG.epoch.eventcondType}');

dataIL = EEG.data(:, :, strcmpi(e_type, 'IL'));
save(fullfile('/Users/yen-hsunwu/Dropbox (ASU)/HERMES/ushape', [EEG.setname, ' IL.mat']), 'dataIL');
dataTR = EEG.data(:, :, strcmpi(e_type, 'TR'));
save(fullfile('/Users/yen-hsunwu/Dropbox (ASU)/HERMES/ushape', [EEG.setname, ' TR.mat']), 'dataTR');
dataPT = EEG.data(:, :, strcmpi(e_type, 'PT'));
save(fullfile('/Users/yen-hsunwu/Dropbox (ASU)/HERMES/ushape', [EEG.setname, ' PT.mat']), 'dataPT');

%% remove eye artifacts
EEG = pop_loadset('filename','sub-09_rmeye.set','filepath','/Users/yen-hsunwu/Dropbox (ASU)/BIDS_format/NCML-EEGrasp/sub-09/eeg/set/');
e_type = cellfun(@(x) x(1), {EEG.epoch.eventcondType}');

dataIL = EEG.data(:, :, strcmpi(e_type, 'IL'));
save(fullfile('/Users/yen-hsunwu/Dropbox (ASU)/HERMES/ushape', [EEG.setname, ' IL.mat']), 'dataIL');
dataTR = EEG.data(:, :, strcmpi(e_type, 'TR'));
save(fullfile('/Users/yen-hsunwu/Dropbox (ASU)/HERMES/ushape', [EEG.setname, ' TR.mat']), 'dataTR');
dataPT = EEG.data(:, :, strcmpi(e_type, 'PT'));
save(fullfile('/Users/yen-hsunwu/Dropbox (ASU)/HERMES/ushape', [EEG.setname, ' PT.mat']), 'dataPT');

%% Only meaningful ICs
EEG = pop_loadset('filename','sub-09_myguess.set','filepath','/Users/yen-hsunwu/Dropbox (ASU)/BIDS_format/NCML-EEGrasp/sub-09/eeg/set/');
e_type = cellfun(@(x) x(1), {EEG.epoch.eventcondType}');

dataIL = EEG.data(:, :, strcmpi(e_type, 'IL'));
save(fullfile('/Users/yen-hsunwu/Dropbox (ASU)/HERMES/ushape', [EEG.setname, ' IL.mat']), 'dataIL');
dataTR = EEG.data(:, :, strcmpi(e_type, 'TR'));
save(fullfile('/Users/yen-hsunwu/Dropbox (ASU)/HERMES/ushape', [EEG.setname, ' TR.mat']), 'dataTR');
dataPT = EEG.data(:, :, strcmpi(e_type, 'PT'));
save(fullfile('/Users/yen-hsunwu/Dropbox (ASU)/HERMES/ushape', [EEG.setname, ' PT.mat']), 'dataPT');



