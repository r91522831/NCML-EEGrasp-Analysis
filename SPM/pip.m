spm('defaults', 'eeg');

S = [];
S.dataset = '/Users/yenhsunw/Dropbox (ASU)/NCML-EEGrasp/EEG/eeglab/S009/S009_pruned_ICA.set';
S.outfile = 'spmeeg_S009_pruned_ICA';
S.channels = 'all';
S.timewin = [];
S.blocksize = 3276800;
S.checkboundary = 1;
S.eventpadding = 0;
S.saveorigheader = 0;
S.conditionlabels = {'Undefined'};
S.inputformat = [];
S.mode = 'epoched';
D = spm_eeg_convert(S);


S = [];
S.task = 'settype';
S.D = '/Users/yenhsunw/Dropbox (ASU)/NCML-EEGrasp/EEG/SPM/S009/ICA pruned/spmeeg_S009_pruned_ICA.mat';
S.ind = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63];
S.type = 'EEG';
S.save = 1;
D = spm_eeg_prep(S);

%%
load(fullfile('/Users/yenhsunw/Dropbox (ASU)/NCML-EEGrasp/EEG/SPM', 'eeg_sensor_info.mat'));
D.sensors = eeg;
load(fullfile('/Users/yenhsunw/Dropbox (ASU)/NCML-EEGrasp/EEG/SPM', 'eeg_channel_info.mat'));
D.channels = channel_info;
sub_id = 'S005';
save(fullfile('/Users/yenhsunw/Dropbox (ASU)/NCML-EEGrasp/EEG/SPM', sub_id, ['spmeeg_', sub_id, '_pruned_ICA.mat']), 'D');

