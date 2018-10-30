% Use the EEGLab to open the raw eeg data before running the script

% put behavior lift onset into EEG event
[filename, pathname, ~] = uigetfile;
% load('/Users/yenhsunw/Dropbox (ASU)/NCML-EEGrasp/behavior/matlab data/sandbox/S009_info_onset_time.mat')
load(fullfile(pathname, filename));

tmp_latency = num2cell(info_onset_time .* EEG.srate + [EEG.event(2:5:end-1).latency]'); % add 's9' time to the event time
for i = 1:length(tmp_latency)
    tmp_latency{i, 2} = 1; % duration
    tmp_latency{i, 3} = 0; % channel
    tmp_latency{i, 4} = []; % bvtime
    tmp_latency{i, 5} = []; % bvmknum
    tmp_latency{i, 6} = 'onset'; % type
    tmp_latency{i, 7} = 'Behavior'; % code
    tmp_latency{i, 8} = []; % urevent
end
tmp_onset = array2table(tmp_latency, 'VariableNames', {'latency', 'duration', 'channel', 'bvtime', 'bvmknum', 'type', 'code', 'urevent'});
tmp_onset = table2struct(tmp_onset);


tmp_event = [EEG.event, tmp_onset'];
[~, tmp_ind] = sort([tmp_event.latency]);
event_with_onset = tmp_event(tmp_ind);
EEG.event = event_with_onset;

save(fullfile(pathname, [filename(1:4), '_info_onset_event4EEG_EEGLab.mat']), 'event_with_onset');
