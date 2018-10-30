clear;
%% Insert behavior onset into EEG evenets for SPM
[filename, pathname, ~] = uigetfile;
% load('/Users/yenhsunw/Dropbox (ASU)/NCML-EEGrasp/EEG/SPM/S009_spm_raw.mat')
load(fullfile(pathname, filename));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% put behavior lift onset into EEG event
[filename_behavior, pathname_behavior, ~] = uigetfile;
% load('/Users/yenhsunw/Dropbox (ASU)/NCML-EEGrasp/behavior/matlab data/preliminary results/S009_info_onset_time.mat')
load(fullfile(pathname_behavior, filename_behavior));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmp_time(:, 4) = num2cell(info_onset_time + [D.trials.events(strcmp([strtrim(string(char(D.trials(1).events.value)))], 's9')).time]'); % add 's9' time to the event time
for i = 1:length(tmp_time)
    tmp_time{i, 1} = 'Behavior'; % type
    tmp_time{i, 2} = 'onset'; % value
    tmp_time{i, 3} = 1; % duration
    
    tmp_time{i, 5} = 0; % offset
end
tmp_onset = array2table(tmp_time, 'VariableNames', {'type', 'value', 'duration', 'time', 'offset'});
tmp_onset = table2struct(tmp_onset);

tmp_event = [D.trials.events, tmp_onset'];
[~, tmp_ind] = sort([tmp_event.time]);
event_with_onset = tmp_event(tmp_ind);
D.trials.events = event_with_onset;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save(fullfile(pathname_behavior, [filename_behavior(1:4), '_info_onset_event4EEG_SPM.mat']), 'event_with_onset');
save(fullfile(pathname, filename), 'D');
