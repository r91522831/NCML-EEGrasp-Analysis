clear; close all; clc;

% EEGLAB history file generated on the 11-Oct-2018
% ------------------------------------------------
% '/Users/yenhsunw/Dropbox (ASU)/NCML-EEGrasp/EEG/Data/Sxxx_EEG/', '*.vhdr'
[raw_filename, raw_dir, ~] = uigetfile('*.vhdr','Select the raw EEG data file');
[tmp_dir, foldername, ~]= fileparts(fileparts(raw_dir));
sub_id = foldername(1:4);
% '/Users/yenhsunw/Dropbox (ASU)/NCML-EEGrasp/EEG/eeglab/Sxxx/';
[EEG_dir, ~, ~] = fileparts(tmp_dir);
output_dir = fullfile(EEG_dir, 'eeglab', sub_id);
[tmp_dir, ~, ~] = fileparts(EEG_dir);
behavior_dir = fullfile(tmp_dir, 'behavior');

% EEGLab
[ALLEEG, EEG, ~, ALLCOM] = eeglab;
pop_editoptions('option_single', false); % make sure the EEG.data precision is 'double' not 'single'!
EEG.etc.eeglabvers = '15.1.1'; % this tracks which version of EEGLAB is being used, you may ignore it
EEG = pop_loadbv(raw_dir, raw_filename);
EEG = eeg_checkset( EEG );
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'setname','S009_raw','gui','off');
if ~isfolder(output_dir)
    mkdir(output_dir);
end
EEG = pop_saveset( EEG, 'filename', [EEG.setname, '.set'], 'filepath', output_dir);
EEG = pop_chanedit(EEG, 'load', {'/Users/yen-hsunwu/Dropbox (ASU)/NCML-EEGrasp/EEG/wg64xyz.xyz', 'filetype', 'xyz'}, 'settype', {'1:63', 'EEG'}, 'settype',{'64', 'EOG'});
% run insert_behavior_event_in2EEG to put behavior onset into EEG events
behavior_results_dir = fullfile(behavior_dir, 'matlab data/preliminary results');
behavior_filename = [sub_id, '_info_onset_time.mat'];
EEG = insertEvent2EEG(EEG, behavior_results_dir, behavior_filename);
EEG.setname = [sub_id, '_channel_loc_lift_onset'];
EEG = pop_saveset( EEG, 'filename', [EEG.setname, '.set'], 'filepath', output_dir);
% lowpass at 512 Hz
EEG = pop_eegfiltnew(EEG, 'hicutoff', 512);
EEG.setname = [sub_id, '_lowpass512Hz'];
EEG = pop_saveset( EEG, 'filename', [EEG.setname, '.set'], 'filepath', output_dir);
% Downsample to 256 Hz
EEG = pop_resample( EEG, 256);
EEG.setname = [sub_id, '_resampled256Hz'];
EEG = pop_saveset( EEG, 'filename', [EEG.setname, '.set'], 'filepath', output_dir);
% highpass at 1Hz
EEG = pop_eegfiltnew(EEG, 'locutoff', 1);
EEG.setname = [sub_id, '_highpass1Hz'];
EEG = pop_saveset( EEG, 'filename', [EEG.setname, '.set'], 'filepath', output_dir);

% reject bad data use ASR
% Keep original EEG
originalEEG = EEG;
% select only EEG channels
EEG = pop_select(EEG, 'channel', {EEG.chanlocs(strcmp({EEG.chanlocs.type}, 'EEG')).labels});
EEG = clean_rawdata(EEG, 5, [0.25 0.75], 0.8, 4, 5, 0.5);
% putback EOG channel
EEG = putback_nonEEG(EEG, originalEEG, EEG.etc.clean_sample_mask);
[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
% eeglab redraw
EEG.setname = [sub_id, '_ASRclean'];
EEG = pop_saveset( EEG, 'filename', [EEG.setname, '.set'], 'filepath', output_dir);

% cleanline

% Apply average reference after adding initial reference
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
[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
% eeglab redraw
EEG.setname = [sub_id, '_rereference'];
EEG = pop_saveset( EEG, 'filename', [EEG.setname, '.set'], 'filepath', output_dir);

% Discard channels to make the data full ranked.
% select only EEG channels
EEG = pop_select(EEG, 'channel', {EEG.chanlocs(strcmp({EEG.chanlocs.type}, 'EEG')).labels});
dataRank = rank(EEG.data');
channelSubset = loc_subsets(EEG.chanlocs, dataRank);
EEG = pop_select(EEG, 'channel', channelSubset{1});
EEG = pop_chanedit(EEG, 'eval', 'chans = pop_chancenter( chans, [], [] );');
% putback EOG channel
EEG = putback_nonEEG(EEG, originalEEG);
[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
% eeglab redraw
EEG.setname = [sub_id, '_adjust_rank'];
EEG = eeg_checkset( EEG );
EEG = pop_saveset( EEG, 'filename', [EEG.setname, '.set'], 'filepath', output_dir);

% epoch
originalEEG_clean = EEG;
% find the tightest window for epoching
ind_event = round([[EEG.event( strcmp({EEG.event.type}, 's9') ).latency]', ...
                   [EEG.event( strcmp({EEG.event.type}, 's17') ).latency]', ...
                   [EEG.event( strcmp({EEG.event.type}, 'onset') ).latency]', ...
                   [EEG.event( strcmp({EEG.event.type}, 's33') ).latency]', ...
                   [EEG.event( strcmp({EEG.event.type}, 's65') ).latency]', ...
                   [EEG.event( strcmp({EEG.event.type}, 's129') ).latency]']);
ind_b4afonset = ceil([max(ind_event(:, 1) - ind_event(:, 3)) / EEG.srate, min(ind_event(:, end) - ind_event(:, 3)) / EEG.srate] * 100) / 100;
EEG = pop_epoch( EEG, {  'onset'  }, ind_b4afonset, 'newname', [sub_id, '_epochs'], 'epochinfo', 'yes');
EEG = eeg_checkset( EEG ); % for some mysterious reason this line is not exacuted!!!!!!!!!!!!!!!!!!!!!!!!!!!
[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
% eeglab redraw
EEG = pop_saveset( EEG, 'filename', [EEG.setname, '.set'], 'filepath', output_dir);
disp(['data size', size(EEG.data)])

% run ICA
% update original EEG
EEG = eeg_checkset( EEG );
originalEEG_epoched = EEG;
disp(['data size', size(EEG.data)])
% select only EEG channels
EEG = pop_select(EEG, 'channel', {EEG.chanlocs(strcmp({EEG.chanlocs.type}, 'EEG')).labels});
EEG = pop_runica(EEG, 'chanind', [], 'extended', 1, 'interupt', 'on');
% putback EOG channel
EEG = putback_nonEEG(EEG, originalEEG_epoched);
[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
eeglab redraw
EEG.setname = [sub_id, '_ICA'];
EEG = eeg_checkset( EEG );
[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = pop_saveset( EEG, 'filename', [EEG.setname, '.set'], 'filepath', output_dir);

% run MARA
% update original EEG
originalEEG_ICA = EEG;
% select only EEG channels and get MARA ICs rejection suggestion
EEG = pop_select(EEG, 'channel', {EEG.chanlocs(strcmp({EEG.chanlocs.type}, 'EEG')).labels});
[ALLEEG, EEG, CURRENTSET] = processMARA( ALLEEG, EEG, CURRENTSET );
rejIC_MARA = logical(EEG.reject.gcompreject');
% update original EEG
EEG = originalEEG_ICA;
EEG = eeg_checkset( EEG );
[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);

% run SASICA
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

% run ICLabel
EEG = pop_iclabel(EEG);
[~, class_i] = max(EEG.etc.ic_classification.ICLabel.classifications, [], 2);
ic_label = cell(size(class_i));
for i = 1:length(class_i)
    ic_label(i) = EEG.etc.ic_classification.ICLabel.classes(class_i(i));
end

% combine all suggested ICs rejection (MARA, SASICA, ICLabel)
% Keep ICLabel brain ICs.
% When ICLabel is other, check MARA or SASICA classification. Reject the IC if at least one rejection from MARA and SASICA.
rejIC_ICLabel_nonBrain = ~(strcmp(ic_label, 'Brain') | strcmp(ic_label, 'Other'));
rejIC_ICLabel_other = strcmp(ic_label, 'Other');
rejIC_other_vote = (rejIC_MARA | rejIC_SASICA) & rejIC_ICLabel_other;
rejIC_final = rejIC_other_vote | rejIC_ICLabel_nonBrain;

%% Naked eye ICA inspection

%% Manually select ICA component reject 
% % load back EOG channel
% EEG = pop_loadset('filename', [sub_id, '_ICA.set'], 'filepath', output_dir);
% EEG = eeg_checkset( EEG );
% % artifact components table
% S002_ac = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 25, 26, 27, 28, 31, 33, 34, 35, 36, 37, 38, 41, 42, 43, 44, 47, 48, 51, 55, 57];
% S003_ac = [2, 3, 5, 6, 7, 8, 9, 10, 15, 17, 18, 20, 21, 22, 24, 28, 31, 36, 37, 38, 39, 41, 42, 44, 45, 46, 47, 48, 49, 51, 52, 53, 55, 57, 61];
% S004_ac = [1, 2, 3, 4, 7, 13, 15, 17, 18, 19, 20, 21, 22, 23, 25, 26, 27  28, 29, 30, 31, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 48, 49, 50, 53, 54, 56, 57];
% S005_ac = [1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 14, 17, 18, 22, 25, 27, 28, 29, 30, 31, 32, 33, 35, 36, 37, 38, 39, 40, 43, 45, 47, 48, 49, 50, 51, 52, 53, 55, 57, 59, 60, 62, 63];
% S006_ac = [1, 2, 3, 4, 5, 6, 8, 11, 12, 13, 14, 15, 17, 18, 19, 20, 21  22, 23, 24, 25, 27, 28, 29, 30, 32, 33, 34, 36, 37, 38, 39, 40, 41, 42, 43, 44, 46, 47, 48, 49, 53, 54, 55, 60, 62, 63];
% S009_ac = [1, 2, 3, 4, 8, 9, 13, 15, 17, 19, 21, 23, 24, 26, 27, 28, 30, 31, 34, 35, 36, 38, 39, 40, 41, 42, 43, 47, 48, 49, 51, 53, 55, 57, 59, 61, 62];
% ac_table = cell2table({S002_ac; S003_ac; S004_ac; S005_ac; S006_ac; S009_ac}, 'VariableNames', {'artifact_component'}, 'RowNames', {'S002', 'S003', 'S004', 'S005', 'S006', 'S009'});

% remove artifact components
% EEG = pop_subcomp( EEG, ac_table{sub_id, :}{:}, 0);
% EEG = eeg_checkset( EEG );
% % % EEG = pop_selectcomps(EEG, 1:63 );
% % % EEG = eeg_checkset( EEG );

%% Remove the artifact components
EEG = pop_subcomp( EEG, find(rejIC_final), 0);
EEG.setname = [sub_id, '_pruned_ICA'];
EEG = eeg_checkset( EEG );
[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = pop_saveset( EEG, 'filename', [EEG.setname, '.set'], 'filepath', output_dir);

% Get experiment conditions
cond_names = {'ALL', 'IL', 'TR', 'PT1', 'PT2', 'PT3'};
% /Users/yenhsunw/Dropbox (ASU)/NCML-EEGrasp/behavior/matlab data/Sxxx
sub_dir = fullfile('/Users/yenhsunw/Dropbox (ASU)/NCML-EEGrasp/behavior/matlab data/', sub_id);
file_list = dir(fullfile(sub_dir, '*.csv'));
tmp = char({file_list.name});
cond = tmp(:, [11:12, 16:18]);

% Get peak roll for each trial
% load peak_roll
load(fullfile('/Users/yenhsunw/Dropbox (ASU)/NCML-EEGrasp/behavior/matlab data/preliminary results/', [sub_id, '_temp_result.mat']));
% sort peak roll according to conditions
pRoll{1, 1} = peak_roll(:, 1);
pRoll{1, 2} = peak_roll(cond(:, 1) == 'I', 1);
pRoll{1, 3} = peak_roll(cond(:, 1) == 'T', 1);
pRoll{1, 4} = peak_roll(and(cond(:, 1) == 'P', cond(:, end) == '1'), 1);
pRoll{1, 5} = peak_roll(and(cond(:, 1) == 'P', cond(:, end) == '2'), 1);
pRoll{1, 6} = peak_roll(and(cond(:, 1) == 'P', cond(:, end) == '3'), 1);
pRoll = cell2table(pRoll, 'VariableNames', cond_names);
%%%%% maybe the roll angle should be plotted on time x trials plot


%% For each component, plot the three bands (alpha, beta, and theta) power on time x trials (IL: 1~19, TR: 1-19, PT1: 1-19, PT2: 1-19, PT3: 1-19)










%{
%%
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
    figure; [tf_ersp{i}, tf_itc{i}, tf_powbase{i}, tf_times{i}, tf_freqs{i}, ~, ~, tf_data{i}] = pop_newtimef( EEG, 1, electrodes{i}, [-1500, 2500], [3, 0.5] , 'topovec', 15, 'elocs', EEG.chanlocs, 'chaninfo', EEG.chaninfo, 'caption', electrodes_name{i}, 'baseline',[-600, -100], 'basenorm', 'on', 'trialbase', 'full', 'plotitc' , 'off', 'plotphase', 'off', 'padratio', 1, 'trialbase', 'full', 'winsize', 512);
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
