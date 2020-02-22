clear; close all; clc;
% EEGLAB history file generated on the 11-Oct-2018 by Yen-Hsun Wu @ ASU
% modified after the EEGLab workshop 2018 @ San Diego
% modified 27-Nov-2018 by Yen-Hsun Wu @ ASU
% modified 14-Feb-2020 by Yen @ ASU
% ------------------------------------------------

%% Section 1: Select raw EEG .vhdr file
% Option of whether the epoch is time-locking at lift onset or at the middle of holding phase
switch input('Does the epoch time locked at lift onset (no if at holding)?(y/N)', 's')
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
All_selected_sub = input('Which subject(s) to preprocess EEG? ');
if isempty(All_selected_sub)
    All_selected_sub = 1:length(All_data_list);
end

%% Start to prcess individual subjects
for All_sub_i = All_selected_sub
    clearvars -except All_*; close all;
    
    %% Section 1: Load EEG raw data and channel locations with EEGLab
    % Step 1: get path and file name for raw EEG data in BrainVision format
    disp(['Processing data for ', All_data_list(All_sub_i).name, ' ...'])
    raw_dir = fullfile(All_data_list(All_sub_i).folder, All_data_list(All_sub_i).name, 'eeg', 'bv');
    tmp_filelist = dir(fullfile(raw_dir, '*.vhdr'));
    raw_filename = tmp_filelist.name;
    clear tmp_*
    
    % Step 2: get the file directory
    % Ex: '/Users/yenhsunw/Dropbox (ASU)/BIDS_format/NCML-EEGrasp/sub-xx/eeg/bv', '*.vhdr'
    % get sub_id
    [~, sub_id, ~]= fileparts(fileparts(fileparts(raw_dir)));
    % set EEG output path
    % output_dir: '/Users/yenhsunw/Dropbox (ASU)/NCML-EEGrasp/EEG/eeglab/001 Process/Sxxx/';
% % %     output_dir = fullfile(All_root_dir, 'NCML-EEGrasp', 'EEG', 'eeglab', '001 Process', [sub_id, '_', All_timelocking_type{:}]);
% % %     output_dir = fullfile(All_root_dir, 'NCML-EEGrasp', 'EEG', 'eeglab', '002 ProcessAgain', [sub_id, '_', All_timelocking_type{:}]);
    output_dir = fullfile(All_root_dir, 'NCML-EEGrasp', 'EEG', 'eeglab', '003 StartOver', [sub_id, '_', All_timelocking_type{:}]);
    if ~isfolder(output_dir)
        mkdir(output_dir)
    end
    
    % Step 3: call EEGLab and import raw EEG data
    % need to download the bva-io extension in EEGLAB
    [ALLEEG, EEG, ~, ALLCOM] = eeglab;
    pop_editoptions('option_single', false); % make sure the EEG.data precision is 'double' not 'single'!
    EEG.etc.eeglabvers = eeglabUpdater.currentVersionNumber; % this tracks which version of EEGLAB is being used, you may ignore it
    EEG = pop_loadbv(raw_dir, raw_filename);
    EEG = eeg_checkset( EEG );
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 0, 'setname', [sub_id, '_raw'], 'gui', 'off');
    EEG = pop_saveset( EEG, 'filename', [EEG.setname, '.set'], 'filepath', output_dir, 'version', '7.3');
    
    % Step 4: Load channel locations
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
    % channel location provided by ANTNeuro
    EEG = pop_chanedit(EEG, 'load', {fullfile(All_data, 'wg64xyz.xyz'), 'filetype', 'xyz'}, 'settype', {'1:63', 'EEG'}, 'settype',{'64', 'EOG'});
    
    % Step 5: Experiment info
    % insert behavior (lift onset) events
    % Step 1: Insert behavior events and Get experiment conditions
    sub_id = EEG.filename(1:6);
    % get behavior raw file directory, behavior_BIDS_dir: /Users/yenhsunw/Dropbox (ASU)/BIDS_format/NCML-EEGrasp/sub-XX/beh/csv
    behavior_BIDS_dir = fullfile(All_data_list(All_sub_i).folder, All_data_list(All_sub_i).name, 'beh', 'csv');
    % get behavior output path, behavior_dir: '/Users/yenhsunw/Dropbox (ASU)/NCML-EEGrasp/behavior/';
    behavior_dir = fullfile(All_root_dir, 'NCML-EEGrasp', 'behavior');
    % run insertEvent2EEG_v1 to put behavior onset into EEG events
% % %     behavior_results_dir = fullfile(behavior_dir, 'preliminary results');
    behavior_results_dir = fullfile(behavior_dir, 'results');
    % behavior_filename: /Users/yenhsunw/Dropbox (ASU)/NCML-EEGrasp/behavior/results/S0XX_info_onset_time.mat
    behavior_filename = ['S0', sub_id(end-1:end), '_info_onset_time.mat'];
    EEG = insertEvent2EEG_v1(EEG, behavior_results_dir, behavior_filename, behavior_BIDS_dir);
    % combine trial condition to event type
    for i = 1:length(EEG.event), EEG.event(i).type = [EEG.event(i).cond, '_', EEG.event(i).type]; end
    
    EEG.setname = [sub_id, '_channel_loc_lift_onset'];
    EEG = eeg_checkset( EEG );
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = pop_saveset( EEG, 'filename', [EEG.setname, '.set'], 'filepath', EEG.filepath);
    %% Section 2: Filtering and downsample
    % Step 1: Lowpass filtering at 128 Hz to remove high freq noise
    sub_id = EEG.filename(1:6);
    EEG = pop_eegfiltnew(EEG, 'hicutoff', 128);
    
    EEG.setname = [sub_id, '_lowpass128Hz'];
    EEG = eeg_checkset( EEG );
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
    % Step 2: Downsample to 256 Hz to reduce computational demand
    EEG = pop_resample( EEG, 256);
    % Step 3: Highpass at 0.5 Hz to remove slow drift
    EEG = pop_eegfiltnew(EEG, 'locutoff', 0.5);
    
    EEG.setname = [sub_id, '_lp128resampled256_hphalfHz'];
    EEG = eeg_checkset( EEG );
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = pop_saveset( EEG, 'filename', [EEG.setname, '.set'], 'filepath', EEG.filepath);
    
    %% Section 3: remove eye movement artifact
    % Step 1: run ICA on continuos data to identify the eye blinking and eye movement components
    % keep the EEG b4 ICA
    originEEG = EEG;
    % Step 2: Bandpass at 1 to 15 Hz to remove slow drift for better ICA
    EEG = pop_eegfiltnew(EEG, 'locutoff', 1, 'hicutoff', 15);
    % Step 3: select only EEG channels
    EEG = pop_select(EEG, 'channel', {EEG.chanlocs(strcmp({EEG.chanlocs.type}, 'EEG')).labels});
    % Step 4: run ICA
    % cudaica can only run on Nvidia graphics chips
% % %     EEG = pop_runica(EEG, 'extended', 1, 'icatype', 'cudaica', 'chanind', [], 'concatenate', 'off', 'verbose', 'off');
    EEG = pop_runica(EEG, 'extended', 1, 'icatype', 'runica', 'chanind', [], 'concatenate', 'off', 'verbose', 'off');
    % Step 5: copy the ICA result back to the original EEG
    originEEG.icawinv = EEG.icawinv;
    originEEG.icasphere = EEG.icasphere;
    originEEG.icaweights = EEG.icaweights;
    originEEG.icachansind = EEG.icachansind;
    originEEG.icaact = (EEG.icaweights * EEG.icasphere) * originEEG.data(EEG.icachansind, :);
    originEEG.setname = [originEEG.filename(1:6), '_reconstructed'];
    EEG = originEEG;
    
    EEG.setname = [sub_id, '_ICA'];
    EEG = eeg_checkset( EEG );
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = pop_saveset( EEG, 'filename', [EEG.setname, '.set'], 'filepath', EEG.filepath);
    % Step 6: remove the eye blink and eye movement component in each channel
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
% % %     Hcor = corr(EEG.data(strcmp({EEG.chanlocs.labels}, 'HEOG'), :)', EEG.icaact')';
% % %     
% % %     H = find(zscore(abs(Hcor)) > 4 == 1)';
% % %     if isempty(H), [~, H] = max( zscore(abs(Hcor)) ); end
    
% % %     pop_topoplot(EEG, 0, [1 3], 'EOG', 0, 'electrodes', 'on');
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    EEG = iclabel(EEG, 'default');
    tmp_id_eyem_class = strcmp(EEG.etc.ic_classification.ICLabel.classes, 'Eye');
    tmp_eyem = EEG.etc.ic_classification.ICLabel.classifications(:, tmp_id_eyem_class);
    %{
    % choose the two largest eye components to reject
    [~, tmp_id_eyem(1)] = max(tmp_eyem);
    tmp_eyem(tmp_id_eyem(1)) = -Inf;
    [~, tmp_id_eyem(2)] = max(tmp_eyem);
    %}
    % remove all eye element with 90% confidence
    tmp_id_eyem = (tmp_eyem > 0.90);
    
    EEG = pop_subcomp( EEG, tmp_id_eyem, 0);
    EEG = eeg_checkset( EEG );
    EEG.nbic = size(EEG.icaact, 1);
    % Step 7: Remove HEOG channel
    EEG = pop_select(EEG, 'channel', {EEG.chanlocs(strcmp({EEG.chanlocs.type}, 'EEG')).labels});
    
    EEG.setname = [sub_id, '_eyeartifact_removed'];
    EEG = eeg_checkset( EEG );
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = pop_saveset( EEG, 'filename', [EEG.setname, '.set'], 'filepath', EEG.filepath);
    
    %% Section 4: Clean data
    % Step 1: Reject bad data use ASR
    % Keep original EEG
    originalEEG = EEG;
%     EEG = clean_rawdata(EEG, 5, [0.25 0.75], 0.8, 4, 5, 0.5);
    EEG = clean_artifacts(EEG, 'WindowCriterion', 0.5);
    % Interpolate channels.
    EEG = pop_interp(EEG, originalEEG.chanlocs, 'spherical');

    EEG.setname = [sub_id, '_ASRclean'];
    EEG = eeg_checkset( EEG );
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = pop_saveset( EEG, 'filename', [EEG.setname, '.set'], 'filepath', EEG.filepath);
    
    %% Section 5: Epoch around onset
    % Step 1:
    % epoch with a window -1500 to 2000 ms around the key event
% % %     ind_win = [-1.5, 2.5];
    % epoch with a window -8000 to 4000 ms around the key event to include
    % the left/right cue
% % %     ind_win = [-8, 4];
    % epoch with a window -3000 to 6000 ms around the key event to include
    ind_win = [-4, 6];
    tmp_type = {EEG.event.type};
    tmp_onset_typeid = cell2mat(cellfun(@contains, {EEG.event.type}, repmat(All_timelocking_type, size({EEG.event.type})), 'UniformOutput', false));
    tmp_epoch_type = unique(tmp_type(tmp_onset_typeid));
    EEG = pop_epoch( EEG, tmp_epoch_type, ind_win, 'newname', [sub_id, '_epochs'], 'epochinfo', 'yes');
    EEG.etc.epoch_latency = ind_win;
    clear tmp*
    
    EEG.setname = [sub_id, '_epoched'];
    EEG = eeg_checkset( EEG );
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = pop_saveset( EEG, 'filename', [EEG.setname, '.set'], 'filepath', EEG.filepath);
end
