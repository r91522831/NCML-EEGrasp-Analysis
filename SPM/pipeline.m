spm('defaults', 'eeg');

S = [];
S.dataset = '/Users/yenhsunw/Dropbox (ASU)/NCML-EEGrasp/EEG/Data/S009_EEG/S_J_2018-09-21_11-02-02.eeg';
S.mode = 'continuous';
S.channels = {'all'};
S.outfile = 'S009_raw';
S.eventpadding = 0;
S.blocksize = 3276800;
S.checkboundary = 1;
S.saveorigheader = 0;
S.timewin = [];
S.conditionlabels = {'Undefined'};
S.inputformat = [];
D = spm_eeg_convert(S);


S = [];
S.D = '/Users/yenhsunw/Dropbox (ASU)/NCML-EEGrasp/EEG/SPM/S009/S009_raw.mat';
S.mode = 'write';
S.blocksize = 655360;
S.prefix = 'M';
S.montage = '/Users/yenhsunw/Dropbox (ASU)/NCML-EEGrasp/EEG/SPM/S009/avref_eog.mat';
S.keepothers = 0;
S.keepsensors = 1;
S.updatehistory = 1;
D = spm_eeg_montage(S);


S = [];
S.D = '/Users/yenhsunw/Dropbox (ASU)/NCML-EEGrasp/EEG/SPM/S009/MS009_raw.mat';
S.type = 'butterworth';
S.band = 'high';
S.freq = 0.5;
S.dir = 'twopass';
S.order = 5;
S.prefix = 'f';
D = spm_eeg_filter(S);


S = [];
S.D = '/Users/yenhsunw/Dropbox (ASU)/NCML-EEGrasp/EEG/SPM/S009/fMS009_raw.mat';
S.fsample_new = 256;
S.method = 'resample';
S.prefix = 'd';
D = spm_eeg_downsample(S);


S = [];
S.D = '/Users/yenhsunw/Dropbox (ASU)/NCML-EEGrasp/EEG/SPM/S009/dfMS009_raw.mat';
S.type = 'butterworth';
S.band = 'low';
S.freq = 55;
S.dir = 'twopass';
S.order = 5;
S.prefix = 'f';
D = spm_eeg_filter(S);


S = [];
S.D = '/Users/yenhsunw/Dropbox (ASU)/NCML-EEGrasp/EEG/SPM/S009/fdfMS009_raw.mat';
S.trialdef.conditionlabel = 'onset';
S.trialdef.eventtype = 'Behavior';
S.trialdef.eventvalue = {'onset'};
S.trialdef.trlshift = 0;
S.timewin = [-1500
             2500];
S.bc = 1;
S.prefix = 'e';
S.eventpadding = 0;
D = spm_eeg_epochs(S);


S = [];
S.task = 'settype';
S.D = '/Users/yenhsunw/Dropbox (ASU)/NCML-EEGrasp/EEG/SPM/S009/efdfMS009_raw.mat';
S.ind = 64;
S.type = 'EOG';
S.save = 1;
D = spm_eeg_prep(S);


S = [];
S.D = '/Users/yenhsunw/Dropbox (ASU)/NCML-EEGrasp/EEG/SPM/S009/efdfMS009_raw.mat';
S.mode = 'reject';
S.badchanthresh = 0.2;
S.prefix = 'a';
S.append = true;
S.methods.channels = {'all'};
S.methods.fun = 'threshchan';
S.methods.settings.threshold = 80;
S.methods.settings.excwin = 1000;
D = spm_eeg_artefact(S);


