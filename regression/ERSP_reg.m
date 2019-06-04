close all; clear; clc;
All_dirpath = uigetdir();
All_dirlist = dir(fullfile(All_dirpath, 'sub*'));

All_linearmodel_path = fullfile(All_dirpath, 'linear');
if ~exist(All_linearmodel_path, 'dir')
    mkdir(All_linearmodel_path);
end

% 
disp([num2cell((1:length(All_dirlist))'), {All_dirlist.name}']);
selected_sub = input('Which subject(s) to plot erpimage? ');
% choose which data set to use b4 or after CSD
switch input('Use data after CSD?(y/N) ', 's')
    case {'y', 'Y'}
        All_isCSDapplied = true;
    otherwise
        All_isCSDapplied = false;
end


for All_i = selected_sub% 1:length(All_dirlist)
    clearvars -except All_*; close all;
    filepath = fullfile(All_dirpath, All_dirlist(All_i).name);
    filelist = dir(fullfile(filepath, '*.set'));
    keystr = 'eeg_csd';
    filename = {filelist(contains({filelist.name}, keystr)).name};
    filename = filename{1}; % get the first file named as key string.
    subID = filename(1:6);
    
    % load data set
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
    EEG = pop_loadset('filename', filename, 'filepath', filepath);
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
    
    % EEG.data is the data applied csd
    % EEG.dataRaw is the data b4 appling csd
    tfEEG = EEG;
    if All_isCSDapplied
        tfEEG.data = EEG.dataRaw;
    end
    
    % channel id:
    % {'Fz', 'FCz', 'C3', 'CP3', 'C1'} = {6, 41, 15, 47, 44}
    Fz = 6; FCz = 41; C3 = 15; CP3 = 47; C1 = 44;
    electrodes_name = {'Fz', 'FCz', 'C3', 'CP3', 'C1'};
    electrodes = {Fz; FCz; C3; CP3; C1};
    nb_epoch = length(EEG.epoch);
    tf_ersp = cell(length(electrodes), 1);
    tf_itc = tf_ersp; tf_powbase = tf_ersp;
    tf_times = tf_ersp; tf_freqs = tf_ersp;
    tf_data = tf_ersp;
    
    % time frequency analysis for each channel and each epoch
    for i = 1:length(electrodes)
        tmp_ersp = cell(nb_epoch, 1);
        tmp_itc = tmp_ersp; tmp_powbase = tmp_ersp;
        tmp_times = tmp_ersp; tmp_freqs = tmp_ersp;
        tmp_data = tmp_ersp;
        for j = 1:nb_epoch
%             figure;
            [tmp_ersp{j, 1}, tmp_itc{j, 1}, tmp_powbase{j, 1}, tmp_times{j, 1}, tmp_freqs{j, 1}, ~, ~, tmp_data{j, 1}] = ...
            newtimef(tfEEG.data(electrodes{i}, :, j), size(tfEEG.data, 2), [tfEEG.times(1), tfEEG.times(end)], tfEEG.srate, [3, 0.5], ...
                                                      'maxfreq', 35, ...
                                                      'topovec', 15, 'elocs', EEG.chanlocs, 'chaninfo', EEG.chaninfo, 'caption', electrodes_name{i}, ...
                                                      'baseline', [-600, -100], 'basenorm', 'on', 'trialbase', 'full', 'padratio', 1, 'winsize', 512, ...
                                                      'plotitc' , 'off', 'plotphase', 'off', 'plotersp', 'off');
        end
        
        % convert
        tf_ersp{i} = cat(3, tmp_ersp{:, 1});
        tf_itc{i} = cat(3, tmp_itc{:, 1});
        tf_powbase{i} = cat(3, tmp_powbase{:, 1});
        tf_times{i} = cat(3, tmp_times{:, 1});
        tf_freqs{i} = cat(3, tmp_freqs{:, 1});
        tf_data{i} = cat(3, tmp_data{:, 1});
    end
    
    % the continuous variables: roll angle (deg)
    roll_ang = nan(nb_epoch, size(tf_ersp{1, 1}, 2));
    for i = 1:nb_epoch
        roll_ang(i, :) = spline(EEG.behavior.obj_epoch{i, 1}.time, EEG.behavior.obj_epoch{i, 1}.roll, squeeze(tf_times{1, 1}(:, :, 1)));
    end
    % the categorical variables: condition (IL, TR, PT)
    dummy_cond = dummyvar( grp2idx({EEG.epoch.condType}'));
    Model_Cond_IL = dummy_cond(:, 1); Model_Cond_TR = dummy_cond(:, 2); Model_Cond_PT = dummy_cond(:, 3);    
    
    % perform regression using power as depend variable; dummy conditions
    % and err as independ variables
    
    % fit = fitlm(ft, 'ERSP ~ Roll * Cond')
    % E(ERSP) = beta_0 * I(IL) + beta_1 * I(TR) + beta_2 * I(PT) + beta_3 * I(IL) * Roll + beta_4 * I(TR) * Roll + beta_5 * I(PT) * Roll
    % each subject results in 8 beta images in the time frequency plane

    Model_Result = cell(length(electrodes), 1);
    Model_coeff_est = Model_Result; Model_coeff_p = Model_Result; Model_Rsquared = Model_Result;
    Model_Cond = categorical( grp2idx({EEG.epoch.condType}') );
    for time_id = 1:size(tf_ersp{1, 1}, 2)
        Model_Roll = roll_ang(:, time_id);
        for freq_id = 1:size(tf_ersp{1, 1}, 1)
            for electrode_id = 1:length(electrodes)
                Model_ERSP = squeeze( tf_ersp{electrode_id, 1}(freq_id, time_id, :) );
                tbl = table(Model_Cond_IL, Model_Cond_TR, Model_Cond_PT, Model_Cond, Model_Roll, Model_ERSP);
                mdl = fitlm(tbl, 'Model_ERSP ~ Model_Cond_IL + Model_Cond_TR + Model_Cond_PT + Model_Cond_IL:Model_Roll + Model_Cond_TR:Model_Roll + Model_Cond_PT:Model_Roll - 1', 'RobustOpts', 'on');
%                 mdl = fitlm(tbl, 'Model_ERSP ~ Model_Roll * Model_Cond', 'RobustOpts', 'on');
                
                Model_Result{electrode_id, 1}{freq_id, time_id} = mdl;
                Model_coeff_est{electrode_id, 1}(freq_id, time_id, :) = mdl.Coefficients{:, 'Estimate'};
                Model_coeff_p{electrode_id, 1}(freq_id, time_id, :) = mdl.Coefficients{:, 'pValue'};
                Model_Rsquared{electrode_id, 1}(freq_id, time_id, :) = mdl.Rsquared.Ordinary;
            end
        end
    end
    
    tmp_filename = fullfile(All_linearmodel_path, [subID, '_LinearModel']);
    save(tmp_filename, 'Model_Result')
    clear Model_Result
    tmp_filename = fullfile(All_linearmodel_path, [subID, '_LinearModel_coeff']);
    save(tmp_filename, 'Model_*')
    tmp_filename = fullfile(All_linearmodel_path, 'misc');
    save(tmp_filename, 'tf_times', 'tf_freqs', 'electrodes_name')
end