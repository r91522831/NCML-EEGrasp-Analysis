function [EEG] = insertEvent2EEG(EEG, pathname, filename)
% insertEvent2EEG Summary of this function goes here
%   Detailed explanation goes here
behavior = load(fullfile(pathname, filename));

tmp_latency = num2cell(behavior.info_onset_time .* EEG.srate + [EEG.event(2:5:end-1).latency]'); % add 's9' time to the event time

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
for i = 1:length(tmp_event)
    event_with_onset(i).bvmknum = i;
    event_with_onset(i).urevent = i;
end

without_bad_behavior_trial = ~isnan([event_with_onset.latency]);
EEG.event = event_with_onset(without_bad_behavior_trial);
EEG.urevent = rmfield(event_with_onset(without_bad_behavior_trial), 'urevent');
save(fullfile(pathname, [filename(1:4), '_info_onset_event4EEG_EEGLab.mat']), 'event_with_onset');
end
