clear; close all; clc;

% EEGLAB history file generated on the 11-Oct-2018 by Yen-Hsun Wu @ ASU
% modified after the EEGLab workshop 2018 @ San Diego
% modified 27-Nov-2018 by Yen-Hsun Wu @ ASU
% ------------------------------------------------

%% Section 1: Select raw EEG .vhdr file

% Option if the goal of the process is in IC space or in Electrode space
switch input('Do you want to make data full rank for analyses in IC space?(y/N)', 's')
    case {'y', 'Y'}
        All_reduce2fullrank = true;
    otherwise
        All_reduce2fullrank = false;
end
% Option if the epoch is time-locking at lift onset
switch input('Does the epoch time locked at lift onset (no if at hold)?(y/N)', 's')
    case {'y', 'Y'}
        All_timelocking_type = {  'onset'  };
    otherwise
        All_timelocking_type = {  'hold'  };
end

disp('Select the project folder in the BIDS_format directory:')
All_data = uigetdir;
All_data_list = dir(fullfile(All_data, 'sub-*'));
All_root_dir = fileparts(fileparts(All_data));

disp([num2cell((1:length(All_data_list))'), {All_data_list.name}']);
selected_sub = input('Which subject(s) to plot erpimage? ');
if isempty(selected_sub)
    selected_sub = 1:length(All_data_list);
end

for All_sub_i = selected_sub %1:length(All_data_list) %6:7 % run only sub-09 and sub-10 for meeting on May 7th    
    clearvars -except All_*; close all;
    % get path and file name for raw EEG data in BrainVision format
    disp(['Processing data for ', All_data_list(All_sub_i).name, ' ...'])
    raw_dir = fullfile(All_data_list(All_sub_i).folder, All_data_list(All_sub_i).name, 'eeg', 'bv');
    tmp_filelist = dir(fullfile(raw_dir, '*.vhdr'));
    raw_filename = tmp_filelist.name;
    clear tmp_filelist
    
    %% Step 1:
    % Ex: '/Users/yenhsunw/Dropbox (ASU)/BIDS_format/NCML-EEGrasp/sub-xx/eeg/bv', '*.vhdr'
    % get sub_id
    [~, sub_id, ~]= fileparts(fileparts(fileparts(raw_dir)));
    % set EEG output path
    % output_dir: '/Users/yenhsunw/Dropbox (ASU)/NCML-EEGrasp/EEG/eeglab/001 Process/Sxxx/';
% % %     output_dir = fullfile(All_root_dir, 'NCML-EEGrasp', 'EEG', 'eeglab', '001 Process', [sub_id, '_', All_timelocking_type{:}]);
    output_dir = fullfile(All_root_dir, 'NCML-EEGrasp', 'EEG', 'eeglab', '002 ProcessAgain', [sub_id, '_', All_timelocking_type{:}]);
    
    %% EEGLab
    % Import raw EEG data
    [ALLEEG, EEG, ~, ALLCOM] = eeglab;
    pop_editoptions('option_single', false); % make sure the EEG.data precision is 'double' not 'single'!
    EEG.etc.eeglabvers = 'develop'; % this tracks which version of EEGLAB is being used, you may ignore it
    
    EEG = pop_loadbv(raw_dir, raw_filename);
    
    EEG = eeg_checkset( EEG );
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 0, 'setname', [sub_id, '_raw'], 'gui', 'off');
    if ~isfolder(output_dir)
        mkdir(output_dir);
    end
    EEG = pop_saveset( EEG, 'filename', [EEG.setname, '.set'], 'filepath', output_dir);
    
    %% Section2: Experiment info
    % Load channel locations and insert behavior (lift onset) events
    % Step 1: Load channel locations
    ANTNeuro_montage = { 'Fp1', 'Fpz', 'Fp2', 'F7', 'F3', 'Fz', 'F4', 'F8', ...
                         'FC5', 'FC1', 'FC2', 'FC6', ...
                         'M1', 'T7', 'C3', 'Cz', 'C4', 'T8', 'M2', ...
                         'CP5', 'CP1', 'CP2', 'CP6', 'P7', 'P3', 'Pz', 'P4', 'P8', ...
                         'POz', 'O1', 'O2', 'AF7', 'AF3', 'AF4', 'AF8', ...
                         'F5', 'F1', 'F2', 'F6', 'FC3', 'FCz', 'FC4', ...
                         'C5', 'C1', 'C2', 'C6', 'CP3', 'CP4', ...
                         'P5', 'P1', 'P2', 'P6', 'PO5', 'PO3', 'PO4', 'PO6', ...
                         'FT7', 'FT8', 'TP7', 'TP8', 'PO7', 'PO8', 'Oz', 'HEOG' };
    if ~all(ismember({EEG.chanlocs.labels}, ANTNeuro_montage)) % montage are different, only for S002
        EEG = pop_select(EEG, 'channel', { 'Fp1', 'Fpz', 'Fp2', 'F7', 'F3', 'Fz', 'F4', 'F8', ...
                                           'FC5', 'FC1', 'FC2', 'FC6', ...
                                           'M1', 'T7', 'C3', 'Cz', 'C4', 'T8', 'M2', ...
                                           'CP5', 'CP1', 'CP2', 'CP6', 'P7', 'P3', 'Pz', 'P4', 'P8', ...
                                           'POz', 'O1', 'O2', 'AF7', 'AF3', 'AF4', 'AF8', ...
                                           'F5', 'F1', 'F2', 'F6', 'FC3', 'FCz', 'FC4', ...
                                           'C5', 'C1', 'C2', 'C6', 'CP3', 'CP4', ...
                                           'P5', 'P1', 'P2', 'P6', 'PO5', 'PO3', 'PO4', 'PO6', ...
                                           'FT7', 'FT8', 'TP7', 'TP8', 'PO7', 'PO8', 'Oz', 'BIP1' });
        EEG = pop_chanedit(EEG, 'changefield', {64, 'labels', 'HEOG'});
    end 
    
    EEG = pop_chanedit(EEG, 'load', {fullfile(All_data, 'wg64xyz.xyz'), 'filetype', 'xyz'}, 'settype', {'1:63', 'EEG'}, 'settype',{'64', 'EOG'});
    % Step 2: Insert behavior events
    % get behavior output path
    % behavior_dir: '/Users/yenhsunw/Dropbox (ASU)/NCML-EEGrasp/behavior/';
    behavior_dir = fullfile(All_root_dir, 'NCML-EEGrasp', 'behavior');
    behavior_BIDS_dir = fullfile(All_data_list(All_sub_i).folder, All_data_list(All_sub_i).name, 'beh', 'csv');
    % run insert_behavior_event_in2EEG to put behavior onset into EEG events
% % %     behavior_results_dir = fullfile(behavior_dir, 'preliminary results');
    behavior_results_dir = fullfile(behavior_dir, 'results');
    behavior_filename = ['S0', sub_id(end-1:end), '_info_onset_time.mat'];
    EEG = insertEvent2EEG(EEG, behavior_results_dir, behavior_filename, behavior_BIDS_dir);
    
    EEG.setname = [sub_id, '_channel_loc_lift_onset'];
    EEG = eeg_checkset( EEG );
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = pop_saveset( EEG, 'filename', [EEG.setname, '.set'], 'filepath', output_dir);
    
    %% Section 3: Filtering and downsample
    % Step 1: Lowpass filtering at 512 Hz to remove high freq noise
    EEG = pop_eegfiltnew(EEG, 'hicutoff', 512);
    
    EEG.setname = [sub_id, '_lowpass512Hz'];
    EEG = eeg_checkset( EEG );
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
    % Step 2: Downsample to 256 Hz to reduce computational demand
    EEG = pop_resample( EEG, 256);
    
    EEG.setname = [sub_id, '_resampled256Hz'];
    EEG = eeg_checkset( EEG );
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
    % Step 3: Highpass at 1 Hz to remove slow drift for better ICA
    EEG = pop_eegfiltnew(EEG, 'locutoff', 1);
    
    EEG.setname = [sub_id, '_lp512Hz_resample256Hz_hp1Hz'];
    EEG = eeg_checkset( EEG );
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = pop_saveset( EEG, 'filename', [EEG.setname, '.set'], 'filepath', output_dir);
    
    %% Section 4: Clean data
    % Step 1: Reject bad data use ASR
    % Keep original EEG
    originalEEG = EEG;
    % select only EEG channels
    EEG = pop_select(EEG, 'channel', {EEG.chanlocs(strcmp({EEG.chanlocs.type}, 'EEG')).labels});
%     EEG = clean_rawdata(EEG, 5, [0.25 0.75], 0.8, 4, 5, 0.5);
    EEG = clean_artifacts(EEG, 'WindowCriterion', 0.5);
    % Interpolate channels.
    EEG = pop_interp(EEG, originalEEG.chanlocs, 'spherical');
    % putback EOG channel
    EEG = putback_nonEEG(EEG, originalEEG, EEG.etc.clean_sample_mask);
    
    EEG.setname = [sub_id, '_ASRclean'];
    EEG = eeg_checkset( EEG );
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = pop_saveset( EEG, 'filename', [EEG.setname, '.set'], 'filepath', output_dir);
    % Step 2: Use cleanline plugin to remove 50 or 60 Hz line noise
    
    % Step 3: Apply average reference after adding initial reference channel
    % base on Makoto Miyakoshi's suggestions (https://sccn.ucsd.edu/wiki/Makoto's_preprocessing_pipeline)
    % select only EEG channels
    EEG = pop_select(EEG, 'channel', {EEG.chanlocs(strcmp({EEG.chanlocs.type}, 'EEG')).labels});
    % add 0 to the reference channel
    EEG.nbchan = EEG.nbchan + 1;
    EEG.data(end+1, :) = zeros(1, EEG.pnts);
    EEG.chanlocs(1, EEG.nbchan).labels = 'initialReference';
    EEG = pop_reref(EEG, []);
    % remove the added reference channel
    EEG = pop_select(EEG, 'nochannel', {'initialReference'});
    % putback EOG channel
    EEG = putback_nonEEG(EEG, originalEEG);
    % remove boundary
    EEG.event = EEG.event(~cellfun(@isempty, {EEG.event.urevent}));
    % put back deleted event for 's9' or 's17'; '33' or 's65' or 's129'
    missing_event = find(diff([EEG.event.urevent]) ~= 1);
    for i = 1:length(missing_event)
        switch EEG.event(missing_event(i)).type
            case 's129' % missing 's9'
                % find next 's17'
                tmp_next = missing_event(i) + find(strcmp({EEG.event(missing_event(i):end).type}, 's17'), 1);
                tmp_missing = EEG.event(tmp_next);
                tmp_missing.latency = EEG.event(tmp_next).latency - 3 * EEG.srate; % 3 second earlier than 's17'
                tmp_missing.bvmknum = missing_event(i) + 1;
                tmp_missing.type = 's9';
                tmp_missing.urevent = missing_event(i) + 1;
            case 'hold' % missing 's65'
                tmp_next = missing_event(i) + find(strcmp({EEG.event(missing_event(i):end).type}, 's129'), 1);
                tmp_missing = EEG.event(tmp_next);
                tmp_missing.latency = EEG.event(tmp_next).latency - 3 * EEG.srate; % 3 second earlier than 's129'
                tmp_missing.bvmknum = missing_event(i) + 1;
                tmp_missing.type = 's65';
                tmp_missing.urevent = missing_event(i) + 1;
            case 's9' % missing 's17'
                tmp_next = missing_event(i);
                tmp_missing = EEG.event(tmp_next);
                tmp_missing.latency = EEG.event(tmp_next).latency + 3 * EEG.srate; % 3 second later than 's9'
                tmp_missing.bvmknum = missing_event(i) + 1;
                tmp_missing.type = 's17';
                tmp_missing.urevent = missing_event(i) + 1;
            case 's65' % missing 's129'
                tmp_next = missing_event(i);
                tmp_missing = EEG.event(tmp_next);
                tmp_missing.latency = EEG.event(tmp_next).latency + 3 * EEG.srate; % 3 second later than 's9'
                tmp_missing.bvmknum = missing_event(i) + 1;
                tmp_missing.type = 's129';
                tmp_missing.urevent = missing_event(i) + 1;
            case {'s17', 'onset', 's33'}
                disp('missing event!!!!!!')
            otherwise
        end
        EEG.event = [EEG.event, tmp_missing];
    end
    tmp_event = struct2table(EEG.event);
	sorted_event = sortrows(tmp_event, 'latency');
	EEG.event = table2struct(sorted_event);
    % save dataset
    EEG.setname = [sub_id, '_rereference'];
    EEG = eeg_checkset( EEG );
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = pop_saveset( EEG, 'filename', [EEG.setname, '.set'], 'filepath', output_dir);
    % Step 4: Discard channels to make the data full ranked.
    if All_reduce2fullrank
        % select only EEG channels
        EEG = pop_select(EEG, 'channel', {EEG.chanlocs(strcmp({EEG.chanlocs.type}, 'EEG')).labels});
        dataRank = rank(EEG.data');
        
        % discard channels
        channelSubset = loc_subsets(EEG.chanlocs, dataRank);
        EEG = pop_select(EEG, 'channel', channelSubset{1});
        EEG = pop_chanedit(EEG, 'eval', 'chans = pop_chancenter( chans, [], [] );');
        
        % putback EOG channel
        EEG = putback_nonEEG(EEG, originalEEG);
        
        EEG.setname = [sub_id, '_adjust_rank'];
        EEG = eeg_checkset( EEG );
        [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
        EEG = pop_saveset( EEG, 'filename', [EEG.setname, '.set'], 'filepath', output_dir);
    end
    
    %% Section 5: Epoch around onset
    % Step 1:
    file_list = dir(fullfile(behavior_BIDS_dir, '*.csv'));
%     trial_no = length(file_list);
    trial_no = 95;
    originalEEG_clean = EEG;
    % clean up trials without onset event
    EEG = checkEvent(EEG, 6);
    EEG = eeg_checkset( EEG );
    
    % should re-run epoch before ICA
    % epoch with a window -1500 to 2000 ms around the key event
% %     ind_win = [-1.5, 2.5];
    % epoch with a window -8000 to 4000 ms around the key event to include
    % the left/right cue
    ind_win = [-8, 4];

    tmp_ready = find(strcmp({EEG.event.type}, 's9'));
    tmp_start = tmp_ready';
    use_leftright = false(trial_no, 1);
    if length(tmp_ready) < trial_no
        disp(['The number of ready cues is ', num2str(length(tmp_ready)),' less than ', num2str(trial_no)]);
        tmp_leftright = find(strcmp({EEG.event.type}, 's17'));
        if length(tmp_leftright) < trial_no
            disp(['The number of left/right cues is ', num2str(length(tmp_leftright)),' less than ', num2str(trial_no)]);
        elseif length(tmp_leftright) > trial_no
            disp(['The number of left/right cues is ', num2str(length(tmp_leftright)),' more than ', num2str(trial_no)]);
        else
            tmp_e = [strcmp({EEG.event.type}, 's9')', strcmp({EEG.event.type}, 's17')'];
            tmp_count = 1;
            tmp_start = nan(length(tmp_leftright), 1);
            for i = 1:length(tmp_leftright)
                if tmp_leftright(i) <= 5
                    if tmp_ready(tmp_count) <= 5
                        tmp_start(i) = tmp_ready(tmp_count);
                        tmp_count = tmp_count + 1;
                    end
                    continue;
                end
                tmp_train = tmp_e(tmp_leftright(i) - (1:5), 1);
                if any(tmp_train)
                    tmp_start(i) = tmp_ready(tmp_count);
                    tmp_count = tmp_count + 1;
                else
                    use_leftright(i) = true;
                    tmp_start(i) = tmp_leftright(i);
                end
            end
        end
    elseif length(tmp_ready) > trial_no
        disp(['The number of ready cues is ', num2str(length(tmp_ready)),' more than ', num2str(trial_no)]);
        
    end
    
    tmp_finish = find(strcmp({EEG.event.type}, 's129'));
    tmp_end = tmp_finish';
    use_relax = false(trial_no, 1);
    if length(tmp_finish) < trial_no
        disp(['The number of finish cues is ', num2str(length(tmp_finish)),' less than ', num2str(trial_no)]);
        tmp_relax = find(strcmp({EEG.event.type}, 's65'));
        if length(tmp_relax) < trial_no
            disp(['The number of relax cues is ', num2str(length(tmp_relax)),' less than ', num2str(trial_no)]);
        elseif length(tmp_relax) > trial_no
            disp(['The number of relax cues is ', num2str(length(tmp_relax)),' more than ', num2str(trial_no)]);
        else
            tmp_e = [strcmp({EEG.event.type}, 's65')', strcmp({EEG.event.type}, 's129')'];
            tmp_count = 1;
            tmp_end = nan(length(tmp_relax), 1);
            for i = 1:length(tmp_relax)
                if tmp_relax(i) <= 5
                    if tmp_finish(tmp_count) <= 5
                        tmp_end(i) = tmp_finish(tmp_count);
                        tmp_count = tmp_count + 1;
                    end
                    continue;
                end
                tmp_train = tmp_e(tmp_relax(i) - (1:5), 1);
                if any(tmp_train)
                    tmp_end(i) = tmp_finish(tmp_count);
                    tmp_count = tmp_count + 1;
                else
                    use_relax(i) = true;
                    tmp_end(i) = tmp_relax(i);
                end
            end
        end
    elseif length(tmp_finish) > trial_no
        disp(['The number of finish cues is ', num2str(length(tmp_finish)),' more than ', num2str(trial_no)]);
    end
    
    tmp_ind_event = [tmp_start'; tmp_end'];
    ind_event = nan(length(tmp_ind_event), 4);
    for i = 1:length(tmp_ind_event)
        tmp_event_series = tmp_ind_event(1, i):tmp_ind_event(2, i);
        tmp_epoch_events = {EEG.event(tmp_event_series).type};
        if length(tmp_event_series) >= (4 + 2) % 's9' can be estimated from 's17' {'s17', 's33', 's65', '129'}: EEG triggers + {'onset', 'hold'}: behavior markers
            if use_leftright(i)
                tmp_ready_latency = [EEG.event( tmp_event_series(strcmp(tmp_epoch_events, 's17')) ).latency]' - (3000 * EEG.srate / 1000);
            else
                tmp_ready_latency = [EEG.event( tmp_event_series(strcmp(tmp_epoch_events, 's9')) ).latency]';
            end
            if use_relax(i)
                tmp_finish_latency = [EEG.event( tmp_event_series(strcmp(tmp_epoch_events, 's65')) ).latency]' + (3000 * EEG.srate / 1000);
            else
                tmp_finish_latency = [EEG.event( tmp_event_series(strcmp(tmp_epoch_events, 's129')) ).latency]';
            end
            if ~any(strcmp(tmp_epoch_events, 'onset'))
                continue;
            end
            ind_event(i, :) = round([ tmp_ready_latency, ...  % 's9' or 's17'
                                      [EEG.event( tmp_event_series(strcmp(tmp_epoch_events, 'onset')) ).latency]', ...  % 'onset'
                                      [EEG.event( tmp_event_series(strcmp(tmp_epoch_events, 'hold')) ).latency]', ...  % 'hold'
                                      tmp_finish_latency ]); % 's129'
        end
    end
    
    %{
    % find the tightest window for epoching
    ind_b4afonset = ceil([max(ind_event(:, 1) - ind_event(:, 2)) / EEG.srate, min(ind_event(:, end) - ind_event(:, 2)) / EEG.srate] * 100) / 100;
    ind_b4afhold = ceil([max(ind_event(:, 1) - ind_event(:, 3)) / EEG.srate, min(ind_event(:, end) - ind_event(:, 3)) / EEG.srate] * 100) / 100;
    ind_win = [max([ind_b4afonset(1), ind_b4afhold(1)]), min([ind_b4afonset(2), ind_b4afhold(2)])];
    %}
    EEG = pop_epoch( EEG, All_timelocking_type, ind_win, 'newname', [sub_id, '_epochs'], 'epochinfo', 'yes');
    EEG = eeg_checkset( EEG );
    EEG.etc.epoch_latency = ind_win;
    EEG = eeg_checkset( EEG );
    % Step 2: Get experiment conditions
    
    
    % remove bad trials
    ind_event(any(isnan(ind_event), 2), :) = [];
    
    %!!!!!!!!! condType: IL, TR, PT; cond: IL, TR, PT1, PT2, PT3
    tmp_filename_list = cell(length(ind_event), 1);
    notPT = false(length(ind_event), 1);
    isIL = false(length(ind_event), 1);
    tmp_id = 1;
        
    for i = 1:length(ind_event)
        tmp = char({file_list(tmp_id).name});
        tmp_filename_list{i} = tmp;
        notPT(i) = ~strcmp(tmp(11:12), 'PT');
        isIL(i) = strcmp(tmp(11:12), 'IL');
        tmp_id = tmp_id + 1;
    end
    cond = cell(size(notPT, 1), 1);
    condType = cell(size(notPT, 1), 1);
    condID = cell(size(notPT, 1), 1);
    condTypeID = cell(size(notPT, 1), 1);
    for i = 1:size(notPT, 1)
        tmp_name = tmp_filename_list{i, 1};
        if ~isempty(tmp_name)
            if notPT(i)
                cond(i, 1) = cellstr(tmp_name(11:12));
                if isIL(i)
                    condID(i, 1) = cellstr(tmp_name(17:18));
                else
                    condID(i, 1) = cellstr(num2str(floor(str2double(tmp_name(13:14)) / 2), '%02d'));
                end
                condTypeID(i, 1) = condID(i, 1);
            else
                cond(i, 1) = cellstr(tmp_name([11:12, 18]));
                condID(i, 1) = cellstr(num2str(floor(str2double(tmp_name(13:14)) / 2), '%02d'));
                condTypeID(i, 1) = cellstr((num2str((floor(str2double(tmp_name(13:14)) / 2) - 1)  * 3 + str2double(tmp_name(18)), '%02d')));
            end
            condType(i, 1) = cellstr(tmp_name(11:12));
        end
    end
    for i = 1:length(EEG.epoch)
        EEG.epoch(i).cond = cond{i, 1};
        EEG.epoch(i).condType = condType{i, 1};
        EEG.epoch(i).condID = condID{i, 1};
        EEG.epoch(i).condTypeID = condTypeID{i, 1};
        if ~isempty(tmp_filename_list{i, 1})
            EEG.epoch(i).trialID = tmp_filename_list{i, 1}(7:9);
        end
    end
    for i = 1:length(EEG.event) % for LIMO exp design
        EEG.event(i).cond = EEG.epoch(EEG.event(i).epoch).cond;
        EEG.event(i).condType = EEG.epoch(EEG.event(i).epoch).condType;
        EEG.event(i).condID = EEG.epoch(EEG.event(i).epoch).condID;
        EEG.event(i).condTypeID = EEG.epoch(EEG.event(i).epoch).condTypeID;
        EEG.event(i).trialID = EEG.epoch(EEG.event(i).epoch).trialID;
    end
    
    EEG = eeg_checkset( EEG ); % ISSUE: for some mysterious reason this line is not exacuted!!!!!!!!!!!!!!!!!!!!!!!!!!!
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = pop_saveset( EEG, 'filename', [EEG.setname, '.set'], 'filepath', output_dir);
    disp(['data size', size(EEG.data)])
    
    %% Section 6: ICA
    % Step 1: Use runica() function
    % keep original EEG
    EEG = eeg_checkset( EEG ); % to work around the issue in
    originalEEG_epoched = EEG;
    disp(['data size', size(EEG.data)])
    % select only EEG channels
    EEG = pop_select(EEG, 'channel', {EEG.chanlocs(strcmp({EEG.chanlocs.type}, 'EEG')).labels});
    EEG = pop_runica(EEG, 'chanind', [], 'extended', 1, 'interupt', 'on');
    % putback EOG channel
    EEG = putback_nonEEG(EEG, originalEEG_epoched);
    
    EEG.setname = [sub_id, '_ICA'];
    EEG = eeg_checkset( EEG );
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = pop_saveset( EEG, 'filename', [EEG.setname, '.set'], 'filepath', output_dir);
    %% Step 2: run MARA
    % keep original EEG
    originalEEG_ICA = EEG;
    % select only EEG channels and get MARA ICs rejection suggestion
    EEG = pop_select(EEG, 'channel', {EEG.chanlocs(strcmp({EEG.chanlocs.type}, 'EEG')).labels});
    [ALLEEG, EEG, CURRENTSET] = processMARA( ALLEEG, EEG, CURRENTSET );
    rejIC_MARA = logical(EEG.reject.gcompreject');
    % update original EEG
    EEG = originalEEG_ICA;
    EEG = eeg_checkset( EEG );
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
    % Step 3: run SASICA
    EEG = eeg_SASICA(EEG, 'MARA_enable', 0, 'FASTER_enable', 0, 'FASTER_blinkchanname', 'No channel', ...
                     'ADJUST_enable', 1, 'chancorr_enable', 0, 'chancorr_channames', 'No channel', 'chancorr_corthresh', 'auto 4', ...
                     'EOGcorr_enable', 1, 'EOGcorr_Heogchannames', 'HEOG', 'EOGcorr_corthreshH', 'auto 4', ...
                     'EOGcorr_Veogchannames', 'No channel', 'EOGcorr_corthreshV', 'auto 4', ...
                     'resvar_enable', 0, 'resvar_thresh', 15, ...
                     'SNR_enable', 0, 'SNR_snrcut', 1, 'SNR_snrBL', [-Inf, 0], 'SNR_snrPOI', [0, Inf], ...
                     'trialfoc_enable', 1, 'trialfoc_focaltrialout', 'auto', ...
                     'focalcomp_enable', 1, 'focalcomp_focalICAout', 'auto', ...
                     'autocorr_enable', 1, 'autocorr_autocorrint', 20, 'autocorr_dropautocorr', 'auto', ...
                     'opts_noplot', 1, 'opts_nocompute', 0, 'opts_FontSize',14);
    rejIC_SASICA = EEG.reject.gcompreject';
    % Step 4: run ICLabel
    EEG = pop_iclabel(EEG, 'default');
    [~, class_i] = max(EEG.etc.ic_classification.ICLabel.classifications, [], 2);
    ic_label = cell(size(class_i));
    for i = 1:length(class_i)
        ic_label(i) = EEG.etc.ic_classification.ICLabel.classes(class_i(i));
    end
    % Step 5: combine all suggested ICs rejection (MARA, SASICA, ICLabel)
    % Keep ICs labeled as 'Brain' by ICLabel.
    % When ICs are labeled as 'Other' by ICLabel, check MARA or SASICA classification.
    % Reject the ICs if at least one rejection voted from MARA and SASICA.
    rejIC_ICLabel_nonBrain = ~(strcmp(ic_label, 'Brain') | strcmp(ic_label, 'Other'));
    rejIC_ICLabel_other = strcmp(ic_label, 'Other');
    rejIC_other_vote = (rejIC_MARA | rejIC_SASICA) & rejIC_ICLabel_other;
    rejIC_final = rejIC_other_vote | rejIC_ICLabel_nonBrain;
    % Step 6: Remove the artifact components
    EEG = pop_subcomp( EEG, find(rejIC_final), 0);
    EEG = eeg_checkset( EEG );
    EEG.nbic = size(EEG.icaact, 1);
    % Step 7: Remove HEOG channel
    EEG = pop_select(EEG, 'channel', {EEG.chanlocs(strcmp({EEG.chanlocs.type}, 'EEG')).labels});
    EEG.setname = [sub_id, '_pruned_ICA'];
    EEG = eeg_checkset( EEG );
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = pop_saveset( EEG, 'filename', [EEG.setname, '.set'], 'filepath', output_dir);

    % finish with SXXX_pruned_ICA dataset
    
    %}
    
    
    
    %% Clean up workspace to release memory
    % % % clear originalEEG*
    
    %{
%% Section 7: more behavior info
% Get peak roll for each trial
% load peak_roll
load(fullfile(behavior_dir, 'matlab data/preliminary results/', [sub_id, '_temp_result.mat']));
% sort peak roll according to conditions
pRoll{1, 1} = peak_roll(:, 1);
pRoll{1, 2} = peak_roll(cond(:, 1) == 'I', 1);
pRoll{1, 3} = peak_roll(cond(:, 1) == 'T', 1);
pRoll{1, 4} = peak_roll(and(cond(:, 1) == 'P', cond(:, end) == '1'), 1);
pRoll{1, 5} = peak_roll(and(cond(:, 1) == 'P', cond(:, end) == '2'), 1);
pRoll{1, 6} = peak_roll(and(cond(:, 1) == 'P', cond(:, end) == '3'), 1);
pRoll = cell2table(pRoll, 'VariableNames', cond_names);
%%%%% maybe the roll angle should be plotted on time x trials plot
    %}
    
    %% Section 8: Time-freq analysis on ICs
    
    
    %% Time-freq analysis in tf_transform.m
    %{
%% Section 9: Time-freq analysis on channels
% For each component, plot the three bands (alpha, beta, and theta) power on time x trials (IL: 1~19, TR: 1-19, PT1: 1-19, PT2: 1-19, PT3: 1-19)
cond_names = {'ALL', 'IL', 'TR', 'PT1', 'PT2', 'PT3'};
cond_nb = length(cond_names);
% typeproc - type of processing: 1 process the raw channel data
%                                0 process the ICA component data
typeproc = str2double(input('\nChoose type of processing [1: raw, 0: component]: (default: 0) ', 's'));
if isnan(typeproc)
    typeproc = 0;
    nb_compoent = EEG.nbic;
    topovec_value = EEG.icawinv;
    cap_str = cell(nb_compoent, 1);
    for i = 1:nb_compoent
        cap_str{i} = [' IC ', num2str(i)];
    end
elseif typeproc == 1
    nb_compoent = EEG.nbchan;
    topovec_value = 1:nb_compoent;
    cap_str = {EEG.chanlocs.labels}';
end
tf_ersp = cell(nb_compoent, cond_nb);
tf_itc = cell(nb_compoent, cond_nb);
tf_powbase = cell(nb_compoent, cond_nb);
tf_data = cell(nb_compoent, cond_nb);

fig_dir = fullfile(output_dir, 'figures', sub_id);
if ~isfolder(fig_dir)
    mkdir(fig_dir);
end

EEG_cond = cell(cond_nb, 1);
% EEG_cond = cell(1, 1);
for j = 1:cond_nb
    EEG_cond{j} = EEG;
    if j == 1
        EEG_cond{j}.data = EEG.data;
        EEG_cond{j}.icaact = EEG.icaact;
    else
        EEG_cond{j}.data = EEG.data(:, :, strcmp([EEG.epoch.cond], cond_names{j}));
        EEG_cond{j}.icaact = EEG.icaact(:, :, strcmp([EEG.epoch.cond], cond_names{j}));
    end

    for i = 1:nb_compoent
        h = figure;
        [tf_ersp{i, j}, tf_itc{i, j}, tf_powbase{i, j}, tf_times, tf_freqs, ~, ~, tf_data{i, j}] = ...
            pop_newtimef( EEG_cond{j}, typeproc, i, round(1000 * [EEG.xmin, EEG.xmax]), [3, 0.5], 'topovec', topovec_value, ...
                          'elocs', EEG.chanlocs, 'chaninfo', EEG.chaninfo, 'caption', [cond_names{j}, cap_str{i}], ...
                          'freqs', [0, 35], 'baseline', [-600, -100], 'plotphase', 'off', 'scale', 'abs', 'padratio', 1 ); %'basenorm', 'on', 'trialbase', 'full');

        savefig(h, fullfile(fig_dir, [sub_id, '_fig_', cond_names{j}, '_IC_', num2str(i)]));
        close(h);
        disp([' channel ', num2str(i)]);
    end
end

save(fullfile(output_dir, [sub_id, '_tf_info']), 'tf_data', 'tf_ersp', 'tf_itc', 'tf_powbase', 'tf_times', 'tf_freqs', '-v7.3');
    %}
    
    
    %%
    
    
    
    
    
    
    
    %%
    
    %{
EEG = pop_loadset('filename', [sub_id, '_pruned_ICA.set'], 'filepath', output_dir);



% channel id:
% {'Fz', 'FCz', 'C3', 'CP3', 'C1'} = {6, 41, 15, 47, 44}
Fz = 6; FCz = 41; C3 = 15; CP3 = 47; C1 = 44;
electrodes_name = {'Fz', 'FCz', 'C3', 'CP3', 'C1'};
electrodes = {Fz; FCz; C3; CP3; C1};
tf_ersp = cell(size(electrodes));
tf_itc = cell(size(electrodes));
tf_powbase = cell(size(electrodes));
tf_times = cell(size(electrodes));
tf_freqs = cell(size(electrodes));
tf_data = cell(size(electrodes));
for i = 1:length(electrodes)
    figure; 
    [tf_ersp{i}, tf_itc{i}, tf_powbase{i}, tf_times{i}, tf_freqs{i}, ~, ~, tf_data{i}] = ...
        pop_newtimef( EEG, 1, electrodes{i}, [-1500, 2500], [3, 0.5], ...
                      'topovec', 15, 'elocs', EEG.chanlocs, 'chaninfo', EEG.chaninfo, 'caption', electrodes_name{i}, ...
                      'baseline', [-600, -100], 'basenorm', 'on', 'trialbase', 'full', 'plotitc' , 'off', 'plotphase', 'off', ...
                      'padratio', 1, 'trialbase', 'full', 'winsize', 512);
end

% sort the tf_data according to conditions
for i = 1:length(electrodes)
    tf_data{i, 2} = tf_data{i, 1}(:, :, cond(:, 1) == 'I');
    tf_data{i, 3} = tf_data{i, 1}(:, :, cond(:, 1) == 'T');
    tf_data{i, 4} = tf_data{i, 1}(:, :, and(cond(:, 1) == 'P', cond(:, end) == '1'));
    tf_data{i, 5} = tf_data{i, 1}(:, :, and(cond(:, 1) == 'P', cond(:, end) == '2'));
    tf_data{i, 6} = tf_data{i, 1}(:, :, and(cond(:, 1) == 'P', cond(:, end) == '3'));
end

tf_data = cell2table(tf_data, 'VariableNames', cond_names, 'RowNames', electrodes_name);

% Get 150-350 ms for 4-8 Hz (theta); 450-650 ms for 9-12 Hz (alpha); 450-650 ms for 20-30 Hz (beta) for individual trial
t_tf4theta = cell(size(tf_data));
t_tf4alpha = cell(size(tf_data));
t_tf4beta = cell(size(tf_data));
for i = 1:length(electrodes)
    for j = 1:length(cond_names)
        f_range = and(tf_freqs{i} >= 4, tf_freqs{i} <= 8);
        t_range = and(tf_times{i} >= 150, tf_times{i} <= 350);
        for k = 1:length(tf_data{i, j}{:}(1, 1, :))
            tmp_avg = squeeze(mean(mean(tf_data{i, j}{:}(f_range, t_range, k), 1), 2));
            t_tf4theta{i, j}(k, 1) = tmp_avg;
        end
        f_range = and(tf_freqs{i} >= 9, tf_freqs{i} <= 12);
        t_range = and(tf_times{i} >= 450, tf_times{i} <= 650);
        for k = 1:length(tf_data{i, j}{:}(1, 1, :))
            tmp_avg = squeeze(mean(mean(tf_data{i, j}{:}(f_range, t_range, k), 1), 2));
            t_tf4alpha{i, j}(k, 1) = tmp_avg;
        end
        f_range = and(tf_freqs{i} >= 20, tf_freqs{i} <= 30);
        t_range = and(tf_times{i} >= 450, tf_times{i} <= 650);
        for k = 1:length(tf_data{i, j}{:}(1, 1, :))
            tmp_avg = squeeze(mean(mean(tf_data{i, j}{:}(f_range, t_range, k), 1), 2));
            t_tf4beta{i, j}(k, 1) = tmp_avg;
        end
    end
end
t_tf4theta = cell2table(t_tf4theta, 'VariableNames', cond_names, 'RowNames', electrodes_name);
t_tf4alpha = cell2table(t_tf4alpha, 'VariableNames', cond_names, 'RowNames', electrodes_name);
t_tf4beta = cell2table(t_tf4beta, 'VariableNames', cond_names, 'RowNames', electrodes_name);

% compute average for conditions, channels, and brainwave bands for each trial
t_tf4theta_Fz_FCz = cell(1, length(cond_names));
t_tf4alpha_C3_CP3 = cell(1, length(cond_names));
t_tf4beta_C3_C1 = cell(1, length(cond_names));
for j = 1:length(cond_names)
    t_tf4theta_Fz_FCz{1, j} = mean([t_tf4theta{{'Fz'}, j}{:}, t_tf4theta{{'FCz'}, j}{:}], 2);
    t_tf4alpha_C3_CP3{1, j} = mean([t_tf4alpha{{'C3'}, j}{:}, t_tf4alpha{{'CP3'}, j}{:}], 2);
    t_tf4beta_C3_C1{1, j} = mean([t_tf4beta{{'C3'}, j}{:}, t_tf4beta{{'C1'}, j}{:}], 2);
end
cellfun(@abs, t_tf4theta_Fz_FCz, 'UniformOutput', false);

t_tf_power = [cellfun(@abs, t_tf4theta_Fz_FCz, 'UniformOutput', false); cellfun(@abs, t_tf4alpha_C3_CP3, 'UniformOutput', false); cellfun(@abs, t_tf4beta_C3_C1, 'UniformOutput', false)];
t_tf_power = array2table(t_tf_power, 'VariableNames', cond_names, 'RowNames', {'theta', 'alpha', 'beta'});
t_tf_power_complex = [t_tf4theta_Fz_FCz; t_tf4alpha_C3_CP3; t_tf4beta_C3_C1];
t_tf_power_complex = array2table(t_tf_power_complex, 'VariableNames', cond_names, 'RowNames', {'theta', 'alpha', 'beta'});

% Get 150-350 ms for 4-8 Hz (theta); 450-650 ms for 9-12 Hz (alpha); 450-650 ms for 20-30 Hz (beta)
tf4theta = cell(size(tf_data));
tf4alpha = cell(size(tf_data));
tf4beta = cell(size(tf_data));
for i = 1:length(electrodes)
    for j = 1:length(cond_names)
        f_range = and(tf_freqs{i} >= 4, tf_freqs{i} <= 8);
        t_range = and(tf_times{i} >= 150, tf_times{i} <= 350);
        tmp_avg = squeeze(mean(mean(mean(tf_data{i, j}{:}(f_range, t_range, :), 1), 2), 3));
        tf4theta{i, j} = tmp_avg;
        f_range = and(tf_freqs{i} >= 9, tf_freqs{i} <= 12);
        t_range = and(tf_times{i} >= 450, tf_times{i} <= 650);
        tmp_avg = squeeze(mean(mean(mean(tf_data{i, j}{:}(f_range, t_range, :), 1), 2), 3));
        tf4alpha{i, j} = tmp_avg;
        f_range = and(tf_freqs{i} >= 20, tf_freqs{i} <= 30);
        t_range = and(tf_times{i} >= 450, tf_times{i} <= 650);
        tmp_avg = squeeze(mean(mean(mean(tf_data{i, j}{:}(f_range, t_range, :), 1), 2), 3));
        tf4beta{i, j} = tmp_avg;
    end
end
tf4theta = cell2table(tf4theta, 'VariableNames', cond_names, 'RowNames', electrodes_name);
tf4alpha = cell2table(tf4alpha, 'VariableNames', cond_names, 'RowNames', electrodes_name);
tf4beta = cell2table(tf4beta, 'VariableNames', cond_names, 'RowNames', electrodes_name);

% compute average for conditions, channels, and brainwave bands
tf4theta_Fz_FCz = mean(tf4theta{{'Fz', 'FCz'}, :}, 1);
tf4alpha_C3_CP3 = mean(tf4alpha{{'C3', 'CP3'}, :}, 1);
tf4beta_C3_C1 = mean(tf4beta{{'C3', 'C1'}, :}, 1);
tf_power = [abs(tf4theta_Fz_FCz); abs(tf4alpha_C3_CP3); abs(tf4beta_C3_C1)];
tf_power = array2table(tf_power, 'VariableNames', cond_names, 'RowNames', {'theta', 'alpha', 'beta'});
tf_power_complex = [tf4theta_Fz_FCz; tf4alpha_C3_CP3; tf4beta_C3_C1];
tf_power_complex = array2table(tf_power_complex, 'VariableNames', cond_names, 'RowNames', {'theta', 'alpha', 'beta'});

%'/Users/yenhsunw/Dropbox (ASU)/NCML-EEGrasp/EEG/eeglab/power_chage'
[EEG_dir, ~, ~] = fileparts(output_dir);
tmp_filename = [sub_id, '_tf_power.mat'];
power_change_filename = fullfile(EEG_dir, 'power_change', tmp_filename);
save(power_change_filename, 'tf_power', 'tf_power_complex', 't_tf_power', 't_tf_power_complex', 'pRoll');
    %}

end
