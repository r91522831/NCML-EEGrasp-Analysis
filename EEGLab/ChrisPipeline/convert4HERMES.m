
% prepare data for HERMES
e_type = cellfun(@(x) x(1), {EEG.epoch.eventcondType}');

dataIL = EEG.data(:, :, strcmpi(e_type, 'IL'));
save(fullfile('/Users/yen-hsunwu/Dropbox (ASU)/HERMES/ushape', 'sub-09 SourceLocalized IL.mat'), 'dataIL');
dataTR = EEG.data(:, :, strcmpi(e_type, 'TR'));
save(fullfile('/Users/yen-hsunwu/Dropbox (ASU)/HERMES/ushape', 'sub-09 SourceLocalized TR.mat'), 'dataTR');
dataPT = EEG.data(:, :, strcmpi(e_type, 'PT'));
save(fullfile('/Users/yen-hsunwu/Dropbox (ASU)/HERMES/ushape', 'sub-09 SourceLocalized PT.mat'), 'dataPT');