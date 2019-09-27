function [EEG] = insertEvent2EEG(EEG, pathname, filename, behavior_BIDS_dir)
% insertEvent2EEG Summary of this function goes here
%   Detailed explanation goes here
behavior = load(fullfile(pathname, filename));
tmp_ready = [EEG.event(strcmp({EEG.event(:).type}', 's9')).latency]';

if length(tmp_ready) < 95
    disp(['The number of ready cues is ', num2str(length(tmp_ready)),' less than 95'])
    tmp_leftright = [EEG.event(strcmp({EEG.event(:).type}', 's17')).latency]';
    if length(tmp_leftright) < 95
        disp(['The number of left/right cues is ', num2str(length(tmp_ready)),' less than 95'])
    elseif length(tmp_leftright) > 95
        disp(['The number of left/right cues is ', num2str(length(tmp_ready)),' more than 95'])
    else
        % estimate ready from left/right cue: ready + 3000 ms -> left/right
        disp('Estimate ready cue from left/right cue!')
        disp('Warning!')
        % please think if it is necessary to keep original tmp_ready
        tmp_ready = tmp_leftright - (3000 * EEG.srate / 1000);
    end
elseif size(tmp_ready) > 95
    disp(['The number of ready cues is ', num2str(length(tmp_ready)),' more than 95'])
end

trial_filelist = dir(fullfile(behavior_BIDS_dir, '*.csv'));
trial_list = cellfun(@(x) str2double(x(7:9)), {trial_filelist.name}');
nb_trial = length(tmp_ready);

% for onset time
tmp_onset = behavior.info_onset_time .* EEG.srate;
tmp_behavior = nan(nb_trial, 1);
for i = 1:length(trial_list)
    tmp_behavior(trial_list(i), 1) = tmp_onset(i);
end
tmp_latency = num2cell(tmp_behavior + tmp_ready); % add 's9'(ready) time to the event time

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

% for hold time
tmp_hold = num2cell(0.5 * ([EEG.event(strcmp({EEG.event.type}, 's65')).latency] + [EEG.event(strcmp({EEG.event.type}, 's33')).latency])');
for i = 1:length(tmp_hold)
    tmp_hold{i, 2} = 1; % duration
    tmp_hold{i, 3} = 0; % channel
    tmp_hold{i, 4} = []; % bvtime
    tmp_hold{i, 5} = []; % bvmknum
    tmp_hold{i, 6} = 'hold'; % type
    tmp_hold{i, 7} = 'Behavior'; % code
    tmp_hold{i, 8} = []; % urevent
end
tmp_hold = array2table(tmp_hold, 'VariableNames', {'latency', 'duration', 'channel', 'bvtime', 'bvmknum', 'type', 'code', 'urevent'});
tmp_hold = table2struct(tmp_hold);

tmp_event = [EEG.event, tmp_onset', tmp_hold'];
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

