function [EEG] = insertEvent2EEG_v1(EEG, pathname, filename, behavior_BIDS_dir)
% insertEvent2EEG Summary of this function goes here
%   Detailed explanation goes here
behavior = load(fullfile(pathname, filename));
n_trial_nominal = 95;
trial_filelist = dir(fullfile(behavior_BIDS_dir, '*.csv'));
trial_list = cellfun(@(x) str2double(x(7:9)), {trial_filelist.name}');

tmp_event = struct2table(EEG.event);
tmp_event{:, 'triggerID'} = nan;
tmp_event{:, 'triggerState'} = {'NA'};
tgr_code = {'s9', 's17', 's33', 's65', 's129'};
tgr_state = {'ready', 'left/right', 'hold', 'relax', 'finish'};
ntgr = length(tgr_code);
for i = 1:ntgr
    tmp_event{strcmp(tmp_event{:, 'type'}, tgr_code{i}), 'triggerID'} = i;
    tmp_event{strcmp(tmp_event{:, 'type'}, tgr_code{i}), 'triggerState'} = tgr_state(i);
end

% EEG.event might have duplicated entries
[~, uni_id, nonu_id] = unique([tmp_event.latency]);
dup_id = setdiff(nonu_id, uni_id);
if ~isempty(dup_id)
    % remove duplicated entries in EEG.event
    tmp = tmp_event;
    tmp(dup_id(:), :) = [];
    nonu_id(dup_id(:), :) = [];
    tmp(:, 'bvmknum') = table(nonu_id);
    tmp(:, 'urevent') = table(nonu_id);
    tmp_event = tmp';
end

n_var_event = width(tmp_event);
name_var_event = tmp_event.Properties.VariableNames;

% go through file list and create a table
%!!!!!!!!! condType: IL, TR, PT; cond: IL, TR, PT1, PT2, PT3; condID: ; condTypeID: ; trialID: ; handleSide: L, R;
n_beh_trial = length(trial_list);
name_var_beh = {'cond', 'condID', 'condCount', 'condType', 'condTypeID', 'condTypeCount', 'trialID', 'handleSide'};
n_var_beh = length(name_var_beh);
tmp_trial_tbl = cell2table( cell(n_trial_nominal, n_var_beh), 'VariableNames', name_var_beh);
% use trial_list, i.e. trial ID, to create the correct trial info
for i = 1:n_beh_trial
    tmp_trial_id = trial_list(i);
    tmp_trial_name = trial_filelist(i).name;
    tmp_trial_tbl{tmp_trial_id, 'trialID'} = {tmp_trial_id};
    tmp_trial_tbl{tmp_trial_id, 'handleSide'} = cellstr(tmp_trial_name(20));
    tmp_cond_type = tmp_trial_name(11:12);
    
    if ~strcmp(tmp_cond_type, 'PT')
        tmp_cond = tmp_trial_name(11:12);
        if strcmp(tmp_cond_type, 'IL')
            tmp_cond_id = 1;
            tmp_cond_count = tmp_trial_name(17:18);
        else
            tmp_cond_id = 2;
            tmp_cond_count = num2str(floor(str2double(tmp_trial_name(13:14)) / 2), '%02d');
        end
        tmp_cond_type_id = tmp_cond_id;
        tmp_cond_type_count = tmp_cond_count;
    else
        tmp_cond = tmp_trial_name([11:12, 18]);
        tmp_cond_id = 2 + str2double(tmp_trial_name(18));
        tmp_cond_count = num2str(floor(str2double(tmp_trial_name(13:14)) / 2), '%02d');
        tmp_cond_type_id = 3;
        tmp_cond_type_count = num2str((floor(str2double(tmp_trial_name(13:14)) / 2) - 1)  * 3 + str2double(tmp_trial_name(16:18)), '%02d');
    end
    tmp_trial_tbl{tmp_trial_id, 'condType'} = cellstr(tmp_cond_type);
    tmp_trial_tbl{tmp_trial_id, 'cond'} = cellstr(tmp_cond);
    tmp_trial_tbl{tmp_trial_id, 'condCount'} = cellstr(tmp_cond_count);
    tmp_trial_tbl{tmp_trial_id, 'condTypeCount'} = cellstr(tmp_cond_type_count);
    tmp_trial_tbl{tmp_trial_id, 'condID'} = num2cell(tmp_cond_id);
    tmp_trial_tbl{tmp_trial_id, 'condTypeID'} =  num2cell(tmp_cond_type_id);
    
end

% fill in tmp_event with corresponding trial info
% get all ready events
n_tgr_state = length(tgr_state);
tmp_raw_event = cell(n_tgr_state, 1);
for i = 1:n_tgr_state
    tmp_raw_event{i, 1} = tmp_event(strcmp(tmp_event.triggerState, tgr_state{i}), :);
    tmp_n_trg = height(tmp_raw_event{i, 1});
    if tmp_n_trg < n_trial_nominal
        tmp_missing = n_trial_nominal - tmp_n_trg;
        tmp_empty = cell2table( cell(tmp_missing, n_var_event), 'VariableNames', name_var_event);
        tmp_empty.latency = nan(tmp_missing, 1);
        tmp_empty.duration = nan(tmp_missing, 1);
        tmp_empty.channel = nan(tmp_missing, 1);
        tmp_empty.bvmknum = nan(tmp_missing, 1);
        tmp_empty.urevent = nan(tmp_missing, 1);
        tmp_empty.triggerID = nan(tmp_missing, 1);

        tmp_raw_event{i, 1} = vertcat( tmp_raw_event{i, 1},  tmp_empty);
    elseif tmp_n_trg > n_trial_nominal
        disp('trigger number is larger than trial numbers.')
    end
end

for i = 1:n_tgr_state
    tmp_raw_event{i, 1} = [tmp_raw_event{i, 1}, tmp_trial_tbl];
end

% the ready cue
tmp_ready = tmp_event{strcmp(tmp_event.triggerState, 'ready'), 'latency'};

n_tgr_trial = length(tmp_ready);
if n_beh_trial ~= n_tgr_trial
    disp('trial number and trigger number mismatched!!')
end

% for onset time
tmp_onset = behavior.info_onset_time .* EEG.srate; % info_onset_time is in second, convert onset into frame number for EEG 2048 Hz sampling rate
tmp_behavior = nan(n_tgr_trial, 1);
for i = 1:n_beh_trial
    % trials might be removed
    tmp_behavior(trial_list(i), 1) = tmp_onset(i);
end

tmp_onset_tbl = cell2table( cell(n_tgr_trial, n_var_event), 'VariableNames', name_var_event);
tmp_onset_tbl.latency = tmp_behavior + tmp_ready; % add 's9'(ready) time to the event time
tmp_onset_tbl.type = repmat({'onset'}, n_tgr_trial, 1);
tmp_onset_tbl.code = repmat({'Behavior'}, n_tgr_trial, 1);
tmp_onset_tbl.duration = nan(n_tgr_trial, 1);
tmp_onset_tbl.channel = nan(n_tgr_trial, 1);
tmp_onset_tbl.bvmknum = nan(n_tgr_trial, 1);
tmp_onset_tbl.urevent = nan(n_tgr_trial, 1);
tmp_onset_tbl.triggerID = nan(n_tgr_trial, 1);
tmp_onset_tbl.triggerState = repmat({'NA'}, n_tgr_trial, 1);
tmp_onset_tbl = [tmp_onset_tbl, tmp_trial_tbl(1:n_tgr_trial, :)];

% for hold time
% middle point between hold and relax cues
tmp_hold = mean([tmp_event{strcmp(tmp_event.triggerState, 'relax'), 'latency'}, tmp_event{strcmp(tmp_event.triggerState, 'hold'), 'latency'}], 2);
tmp_hold_tbl = cell2table( cell(n_tgr_trial, n_var_event), 'VariableNames', name_var_event);
tmp_hold_tbl.latency = tmp_hold;
tmp_hold_tbl.type = repmat({'hold'}, n_tgr_trial, 1);
tmp_hold_tbl.code = repmat({'Behavior'}, n_tgr_trial, 1);
tmp_hold_tbl.duration = nan(n_tgr_trial, 1);
tmp_hold_tbl.channel = nan(n_tgr_trial, 1);
tmp_hold_tbl.bvmknum = nan(n_tgr_trial, 1);
tmp_hold_tbl.urevent = nan(n_tgr_trial, 1);
tmp_hold_tbl.triggerID = nan(n_tgr_trial, 1);
tmp_hold_tbl.triggerState = repmat({'NA'}, n_tgr_trial, 1);
tmp_hold_tbl = [tmp_hold_tbl, tmp_trial_tbl(1:n_tgr_trial, :)];

% combine onset and hold back into the events
tmp_event = vertcat(tmp_raw_event{:}, tmp_onset_tbl, tmp_hold_tbl);
[event_with_onset, ~] = sortrows(tmp_event, 'latency');

% % % for i = 1:length(event_with_onset)
% % %     event_with_onset(i).bvmknum = i;
% % %     event_with_onset(i).urevent = i;
% % % end
without_bad_behavior_trial = ~isnan(event_with_onset.latency);
EEG.event = table2struct(event_with_onset(without_bad_behavior_trial, :));
% % % EEG.urevent = rmfield(event_with_onset(without_bad_behavior_trial), 'urevent');
save(fullfile(pathname, [filename(1:4), '_info_onset_event4EEG_EEGLab.mat']), 'event_with_onset');
end

