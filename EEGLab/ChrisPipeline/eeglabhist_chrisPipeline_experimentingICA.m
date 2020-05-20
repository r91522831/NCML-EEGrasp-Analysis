clear; close all; clc;
% EEGLAB history file generated on the 11-Oct-2018 by Yen-Hsun Wu @ ASU
% modified after the EEGLab workshop 2018 @ San Diego
% modified 27-Nov-2018 by Yen-Hsun Wu @ ASU
% modified 14-Feb-2020 by Yen @ ASU
% modified 05-May-2020 by Yen @ASU
% ------------------------------------------------

%% Section ini

%{
% Option of whether the epoch is time-locking at lift onset or at the middle of holding phase
switch input('Does the epoch time locked at lift onset (no if at holding)?(y/N)', 's')
    case {'y', 'Y'}
        All_timelocking_type = {  'onset'  };
    otherwise
        All_timelocking_type = {  'hold'  };
end
%}
All_timelocking_type = {  'onset'  }; % epoch at lift onset

disp('Select the project folder in the BIDS_format directory:')
All_data = uigetdir;
All_data_list = dir(fullfile(All_data, 'sub-*'));
All_root_dir = fileparts(fileparts(All_data));
% path for EEGLab dipfit plugins for source localization
All_path_eeglab_dipfit = fullfile(fileparts(All_root_dir), 'Dropbox (Personal)', 'Programming', 'Matlab', 'myLibrary', 'eeglab2019_1', 'plugins', 'dipfit');

disp([num2cell((1:length(All_data_list))'), {All_data_list.name}']);
All_selected_sub = input('Which subject(s) to preprocess EEG? ');
if isempty(All_selected_sub)
    All_selected_sub = 1:length(All_data_list);
end

%% Start to prcess individual subjects
for All_sub_i = All_selected_sub
    clearvars -except All_*; close all;
    
    %% Section 0: prepare EEGLAB
    % Step 1: get path and file name for raw EEG data in BrainVision format
    % Ex: '/Users/yenhsunw/Dropbox (ASU)/BIDS_format/NCML-EEGrasp/sub-xx/eeg/bv', '*.vhdr'
    disp(['Processing data for ', All_data_list(All_sub_i).name, ' ...'])
    raw_dir = fullfile(All_data_list(All_sub_i).folder, All_data_list(All_sub_i).name, 'eeg', 'bv');
    tmp_filelist = dir(fullfile(raw_dir, '*.vhdr'));
    raw_filename = tmp_filelist.name;
    clear tmp_*
    % get sub_id
    [~, sub_id, ~]= fileparts(fileparts(fileparts(raw_dir)));
    % Step 2: set EEG output path
    % output_dir: /Users/yen-hsunwu/Dropbox (ASU)/BIDS_format/NCML-EEGrasp/sub-09/eeg/set/
    output_dir = fullfile(All_data_list(All_sub_i).folder, All_data_list(All_sub_i).name, 'eeg', 'set');
    if ~isfolder(output_dir)
        mkdir(output_dir)
    end
    % Step 3: call EEGLab
    [ALLEEG, EEG, ~, ALLCOM] = eeglab;
    pop_editoptions('option_single', false); % make sure the EEG.data precision is 'double' not 'single'!
    EEG.etc.eeglabvers = eeglabUpdater.currentVersionNumber; % this tracks which version of EEGLAB is being used, you may ignore it
    
    %% Section 1: Load EEG raw data and channel locations with EEGLab
    % Step 1: import raw EEG data
    % need to download the bva-io extension in EEGLAB
    EEG = pop_loadbv(raw_dir, raw_filename);
    EEG = eeg_checkset( EEG );
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 0, 'setname', [sub_id, '_raw'], 'gui', 'off');
    EEG = pop_saveset( EEG, 'filename', [EEG.setname, '.set'], 'filepath', output_dir, 'version', '7.3' );
    %% Step 2: Load channel locations
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
    % standard channel location
    EEG = pop_chanedit(EEG, 'lookup', fullfile(All_path_eeglab_dipfit, 'standard_BESA', 'standard-10-5-cap385.elp'), 'settype', {'1:63', 'EEG'}, 'settype',{'64', 'EOG'});
    EEG.urchanlocs = EEG.chanlocs; % preserve the original locations of the original data channels
    % Step 3: Experiment info
    % Insert behavior events (lift onset, etc.) and Get experiment conditions
    sub_id = EEG.filename(1:6);
    % get behavior raw file directory, behavior_BIDS_dir: /Users/yenhsunw/Dropbox (ASU)/BIDS_format/NCML-EEGrasp/sub-XX/beh/
    behavior_dir = fullfile(All_data_list(All_sub_i).folder, All_data_list(All_sub_i).name, 'beh');
    % behavior_filename: /Users/yenhsunw/Dropbox (ASU)/BIDS_format/NCML-EEGrasp/sub-XX/beh/mat/S0XX_info_onset_time.mat
    behavior_filename = ['S0', sub_id(end-1:end), '_info_onset_time.mat'];
    EEG = insertEvent2EEG_v2(EEG, fullfile(behavior_dir, 'mat'), behavior_filename, fullfile(behavior_dir, 'csv'));
    % combine trial condition to event type
    for i = 1:length(EEG.event), EEG.event(i).type = [EEG.event(i).cond, '_', EEG.event(i).type]; end
    
    EEG.setname = [sub_id, '_channel_loc_lift_onset'];
    EEG = eeg_checkset( EEG );
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = pop_saveset( EEG, 'filename', [EEG.setname, '.set'], 'filepath', EEG.filepath);
    
    %% Section 2: Filtering and downsample
    % Step 0: clean line noise
    EEG = pop_cleanline(EEG, 'bandwidth', 2, 'chanlist', 1:EEG.nbchan, 'computepower', 0, 'linefreqs', [60, 120, 180, 240, 300],...
                             'normSpectrum', 0, 'p', 0.01, 'pad', 2, 'plotfigures', 0, 'scanforlines', 1, 'sigtype', 'Channels', 'tau', 100,...
                             'verb', 1, 'winsize', 4, 'winstep', 4);
    % Step 1: bandpass filtering at 0.5 to 128 Hz to remove slow drift and high freq noise
    EEG = pop_eegfiltnew(EEG, 'locutoff', 0.5, 'hicutoff', 128);
    % Step 2: Downsample to 256 Hz to reduce computational demand
    EEG = pop_resample( EEG, 256);
    
    sub_id = EEG.filename(1:6);
    EEG.setname = [sub_id, '_lp128_hphalfHz_down256Hz'];
    EEG = eeg_checkset( EEG );
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = pop_saveset( EEG, 'filename', [EEG.setname, '.set'], 'filepath', EEG.filepath);
    
    %% Section 3: Clean data
    % Step 1: Reject bad data use ASR
    if any(strcmpi({EEG.chanlocs.type}, 'EOG'))
        % Remove HEOG channel
        EEG = pop_select(EEG, 'channel', {EEG.chanlocs(strcmp({EEG.chanlocs.type}, 'EEG')).labels});
    end
    % Keep original EEG
    originalEEG_b4rmBadChannel = EEG;
    EEG = clean_artifacts(EEG, 'WindowCriterion', 0.5);
    vis_artifacts(EEG, originalEEG_b4rmBadChannel);

    sub_id = EEG.filename(1:6);
    EEG.setname = [sub_id, '_ASRclean'];
    EEG = eeg_checkset( EEG );
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = pop_saveset( EEG, 'filename', [EEG.setname, '.set'], 'filepath', EEG.filepath);
    
    %% clean aggrssively
    EEG = clean_artifacts(EEG, 'WindowCriterion', 0.5, 'BurstCriterion', 1.5);
% % %     most = clean_artifacts(EEG, 'WindowCriterion', 0.5, 'BurstCriterion', 3);
% % %     vis_artifacts(agg, EEG);
    % epoch
    ind_win = [-0.5, 2]; 
    tmp_type = {EEG.event.type};
    tmp_onset_typeid = cell2mat(cellfun(@contains, {EEG.event.type}, repmat(All_timelocking_type, size({EEG.event.type})), 'UniformOutput', false));
    tmp_epoch_type = unique(tmp_type(tmp_onset_typeid));
    EEG = pop_epoch( EEG, tmp_epoch_type, ind_win, 'newname', [EEG.filename(1:6), '_epochs', '_', All_timelocking_type{:}], 'epochinfo', 'yes');
    EEG.etc.epoch_latency = ind_win;
    clear tmp*
    % reduce channel to 39 channels
% % %     EEG = pop_select( EEG, 'channel',{'Fp1' 'Fpz' 'Fp2' 'F7' 'F3' 'Fz' 'F4' 'F8' 'T7' 'C3' 'Cz' 'C4' 'T8' 'P7' 'P3' 'Pz' 'P4' 'P8' 'POz' 'O1' 'O2' 'AF7' 'AF3' 'AF4' 'AF8' 'FC3' 'FCz' 'FC4' 'CP3' 'CP4' 'PO3' 'PO4' 'FT7' 'FT8' 'TP7' 'TP8' 'PO7' 'PO8' 'Oz'});
    % reduce channel to 29 channels
    EEG = pop_select( EEG, 'channel',{'Fp1' 'Fpz' 'Fp2' 'F7' 'F3' 'Fz' 'F4' 'F8' 'T7' 'C3' 'Cz' 'C4' 'T8' 'P7' 'P3' 'Pz' 'P4' 'P8' 'POz' 'O1' 'O2' 'FC3' 'FCz' 'FC4' 'CP3' 'CP4' 'PO3' 'PO4' 'Oz'});
    
    EEG.setname = [sub_id, '_heavilycleaned_shortepoched_reducechannel'];
    EEG = eeg_checkset( EEG );
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = pop_saveset( EEG, 'filename', [EEG.setname, '.set'], 'filepath', EEG.filepath);

    originalEEG = pop_loadset('filename','sub-09_ASRclean.set','filepath', EEG.filepath);
    % epoch
    ind_win = [-0.5, 2]; 
    tmp_type = {originalEEG.event.type};
    tmp_onset_typeid = cell2mat(cellfun(@contains, {originalEEG.event.type}, repmat(All_timelocking_type, size({originalEEG.event.type})), 'UniformOutput', false));
    tmp_epoch_type = unique(tmp_type(tmp_onset_typeid));
    originalEEG = pop_epoch( originalEEG, tmp_epoch_type, ind_win, 'newname', [originalEEG.setname, '_epochs', '_', All_timelocking_type{:}], 'epochinfo', 'yes');
    originalEEG.etc.epoch_latency = ind_win;
    clear tmp*
    % reduce channel to 40 channels
    originalEEG = pop_select( originalEEG, 'channel',{'Fp1' 'Fpz' 'Fp2' 'F7' 'F3' 'Fz' 'F4' 'F8' 'T7' 'C3' 'Cz' 'C4' 'T8' 'P7' 'P3' 'Pz' 'P4' 'P8' 'POz' 'O1' 'O2' 'AF7' 'AF3' 'AF4' 'AF8' 'FC3' 'FCz' 'FC4' 'CP3' 'CP4' 'PO3' 'PO4' 'FT7' 'FT8' 'TP7' 'TP8' 'PO7' 'PO8' 'Oz'});
   
    % Section 5: ICA ICA!!!!
    % remove eye movement artifact
    % run ICA on continuos data to identify the eye blinking and eye movement components
    % Step 2: Highpass at 3 to remove slow drift for better ICA
    EEG = pop_eegfiltnew(EEG, 'locutoff', 3);
    % Step 3: select only EEG channels
    EEG = pop_select(EEG, 'channel', {EEG.chanlocs(strcmp({EEG.chanlocs.type}, 'EEG')).labels});
    % Step 4: run ICA
    % cudaica can only run on Nvidia graphics chips
% % %     EEG = pop_runica(EEG, 'extended', 1, 'icatype', 'cudaica', 'chanind', [], 'concatenate', 'off', 'verbose', 'off');
    EEG = pop_runica(EEG, 'extended', 1, 'icatype', 'runica', 'chanind', [], 'concatenate', 'off', 'verbose', 'off');
    % Step 5: copy the ICA result back to the original EEG
    originalEEG.icawinv = EEG.icawinv;
    originalEEG.icasphere = EEG.icasphere;
    originalEEG.icaweights = EEG.icaweights;
    originalEEG.icachansind = EEG.icachansind;
    originalEEG.icaact = (EEG.icaweights * EEG.icasphere) * originalEEG.data(EEG.icachansind, :);
    EEG = originalEEG;
    
    if isfield(EEG, 'dataRaw'), EEG.setname = [EEG.setname, '_CSD_ICA'];
    else, EEG.setname = [EEG.setname, '_ICA']; end
    EEG = eeg_checkset( EEG );
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = pop_saveset( EEG, 'filename', [EEG.setname, '.set'], 'filepath', EEG.filepath);
    
    % Section 6: check each IC
    % Step 1: label each IC
    EEG = iclabel(EEG, 'default');
    % Step 2: source localization
    COREG = [0.83215, -15.6287, 2.4114, 0.081214, 0.00093739, -1.5732, 1.1742, 1.0601, 1.1485];
    EEG = pop_dipfit_settings( EEG, 'hdmfile', fullfile(All_path_eeglab_dipfit, 'standard_BEM', 'standard_vol.mat'),...
                                    'coordformat', 'MNI', 'mrifile', fullfile(All_path_eeglab_dipfit, 'standard_BEM', 'standard_mri.mat'),...
                                    'chanfile', fullfile(All_path_eeglab_dipfit, 'standard_BEM', 'elec', 'standard_1005.elc'),...
                                    'coord_transform', COREG, 'chansel', 1:EEG.nbchan );
    residual_thres = 40;
    EEG = pop_multifit(EEG, 1:EEG.nbchan, 'threshold', residual_thres, 'plotopt', {'normlen', 'on'});
    pop_dipplot( EEG, 1:EEG.nbchan, 'mri', fullfile(All_path_eeglab_dipfit, 'standard_BEM', 'standard_mri.mat'), 'normlen', 'on');
    
    if isfield(EEG, 'dataRaw'), EEG.setname = [EEG.setname, '_CSD_SouceLocalized'];
    else, EEG.setname = [EEG.setname, '_SouceLocalized']; end
    EEG = eeg_checkset( EEG );
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = pop_saveset( EEG, 'filename', [EEG.setname, '.set'], 'filepath', EEG.filepath);
    
    %% Section 7a: Time-freq decomposition
    tmpEEG = EEG;
    switch input('time-freq analysis on voltage (V) or IC activity (I)?(v/I)', 's')
        case {'v', 'V'}
            icatf = false;
            disp('time-freq on voltage ...')
        otherwise
            icatf = true;
            tmpEEG.data = EEG.icaact;
            disp('time-freq on IC activation ...')
    end
    % No interpolation of channels -> electrodes number equals to ICs number
    electrodes_name = {EEG.chanlocs.labels};
    nelect = length(EEG.chanlocs);
    electrodes = num2cell(1:nelect);
    nepoch = length(EEG.epoch);
    
    % time frequency analysis for each channel and each epoch
    tf_tlim = [tmpEEG.times(1), tmpEEG.times(end)]; % in ms
    tf_ntimesout = 200 * floor(diff(tf_tlim)/4000); % estimating number of time bins to output for time-freq analysis
    tf_ersp = cell(length(electrodes), 1);
    % The progress bar
    ticker = 0;
    h = waitbar(0, 'time frequency analysis.');
    total_iter = nelect * nepoch;
    for i = 1:length(electrodes)        
        % time freq analysis without any baseline normalization or subtraction
% % %         [tf_ersp{i, 1}, tmp_freqs, tmp_times] = timefreq( squeeze(tmpEEG.data(i, :, :)), tmpEEG.srate, ...
% % %                                                                                          'cycles', [3, 0.5], ...
% % %                                                                                          'tlimits', tf_tlim, 'ntimesout', tf_ntimesout, ...
% % %                                                                                          'freqs', [2, 40], 'padratio', 1, ...
% % %                                                                                          'verbose', 'off');
        for j = 1:nepoch
            % tf_ersp unit is power in '10 * log10( \muV^{2}/Hz )'
            [tf_ersp{i, 1}(:, :, j), ~, ~, tmp_times, tmp_freqs, ~, ~, ~] = ...
                newtimef(tmpEEG.data(i, :, j), size(tmpEEG.data, 2), [tmpEEG.times(1), tmpEEG.times(end)], tmpEEG.srate, [3, 0.5], ...
                                               'freqs', [2, 40], 'timesout', tf_ntimesout, ...
                                               'baseline', NaN, 'padratio', 1, ...
                                               'plotitc', 'off', 'plotphase', 'off', 'plotersp', 'off', ...
                                               'scale', 'log', 'verbose', 'off');    
            ticker = ticker + 1;
            progress_precent = 100 * ticker / total_iter;
            waitbar(ticker / total_iter, h, sprintf('time frequency analysis. %2.2f %%', progress_precent));
        end
    end
    close(h);
    
    if icatf, EEG.icatf.tf_times = tmp_times; EEG.icatf.tf_freqs = tmp_freqs; EEG.icatf.tf_ersp = tf_ersp;
    else,     EEG.datatf.tf_times = tmp_times; EEG.datatf.tf_freqs = tmp_freqs; EEG.datatf.tf_ersp = tf_ersp; end
    
    sub_id = EEG.filename(1:6);
    EEG.setname = [sub_id, '_timefreq'];
    EEG = eeg_checkset( EEG );
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = pop_saveset( EEG, 'filename', [EEG.setname, '.set'], 'filepath', EEG.filepath);
    

    
    %% Secion 7b: reject ICs
    % Step 1a: remove eye artifacts 
    eyeIdxBool = EEG.etc.ic_classification.ICLabel.classifications(:, strcmpi(EEG.etc.ic_classification.ICLabel.classes, 'Eye')) >= 0.9; % > 90% eye == eye ICs.
    eyeIdx = find(eyeIdxBool);
    % Perform IC rejection.
    EEG = pop_subcomp(EEG, eyeIdx, 0, 0); % reject eye ICs
    % Post-process to update ICLabel data structure.
    EEG.etc.ic_classification.ICLabel.classifications = EEG.etc.ic_classification.ICLabel.classifications(~eyeIdxBool, :);
    % Post-process to update EEG.icaact.
    EEG.icaact = [];
    EEG = eeg_checkset(EEG, 'ica');
    
    sub_id = EEG.filename(1:6);
    EEG.setname = [sub_id, '_rmeye'];
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = pop_saveset( EEG, 'filename', [EEG.setname, '.set'], 'filepath', EEG.filepath);
    
    
    %%
    % Step 1b: keep only meaningful ICs
    % modified from Makoto's useful EEGLAB code
    brainIdx  = find( EEG.etc.ic_classification.ICLabel.classifications(:, 1) >= 0.7 ); % > 70% brain == good ICs.
 
    % Perform IC rejection using residual variance of the IC scalp maps.
    rvList    = [EEG.dipfit.model.rv];
    goodRvIdx = find(rvList < 0.4)'; % < 40% residual variance == good ICs.
 
    % Perform IC rejection using inside brain criterion.
    load(EEG.dipfit.hdmfile); % This returns 'vol'.
    dipoleXyz = zeros(length(EEG.dipfit.model), 3);
    for icIdx = 1:length(EEG.dipfit.model)
        if isempty(EEG.dipfit.model(icIdx).posxyz), continue; end
        dipoleXyz(icIdx, :) = EEG.dipfit.model(icIdx).posxyz(1, :);
    end
    depth = ft_sourcedepth(dipoleXyz, vol);
    depthThreshold = 1;
    insideBrainIdx = find(depth <= depthThreshold);
    
    % select IC that are meaningful for the UShape object task
    meaningfulIdx = [6, 7, 11, 12, 13, 15, 16, 21, 22, 26, 37];
 
    % Take AND across the three criteria.
    goodIcIdx = intersect(brainIdx, goodRvIdx);
    goodIcIdx = intersect(goodIcIdx, insideBrainIdx);
    goodIcIdx = intersect(goodIcIdx, meaningfulIdx);
 
    % Perform IC rejection.
    EEG = pop_subcomp(EEG, goodIcIdx, 0, 1);
 
    % Post-process to update ICLabel data structure.
    EEG.etc.ic_classification.ICLabel.classifications = EEG.etc.ic_classification.ICLabel.classifications(goodIcIdx,:);
 
    % Post-process to update EEG.icaact.
    EEG.icaact = [];
    EEG = eeg_checkset(EEG, 'ica');
    
    sub_id = EEG.filename(1:6);
    EEG.setname = [sub_id, '_myguess'];
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = pop_saveset( EEG, 'filename', [EEG.setname, '.set'], 'filepath', EEG.filepath);
    
    %% downsample
    
    sub_id = EEG.filename(1:6);
    EEG.setname = [sub_id, '_downsample256Hz'];
    EEG = eeg_checkset( EEG );
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = pop_saveset( EEG, 'filename', [EEG.setname, '.set'], 'filepath', EEG.filepath);
    
    
    % Interpolate channels.
% % %     EEG = pop_interp(EEG, originalEEG_b4rmBadChannel.chanlocs, 'spherical');
    
 
   
 %}
    
    
    
    
    
    
end
