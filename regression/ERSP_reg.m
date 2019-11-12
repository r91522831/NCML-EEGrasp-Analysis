close all; clear; clc;
All_dirpath = uigetdir();
All_dirlist = dir(fullfile(All_dirpath, 'sub*'));

%% load behaviore exponential fit
All_beh_expfit = load(fullfile(fileparts(fileparts(fileparts(All_dirpath))), 'behavior', 'preliminary results', 'behavior_pRoll_23subs_pRoll_fit.mat'));

%% 
disp([num2cell((1:length(All_dirlist))'), {All_dirlist.name}']);
selected_sub = input('Which subject(s) to plot erpimage? ');
if isempty(selected_sub)
    selected_sub = 1:length(All_dirlist);
end

% choose which data set to use b4 or after CSD
switch input('Use data after CSD?(y/N) ', 's')
    case {'y', 'Y'}
        All_isCSDapplied = true;
    otherwise
        All_isCSDapplied = false;
end

if All_isCSDapplied
    All_linearmodel_path = fullfile(All_dirpath, 'linear', 'CSD');
else
    All_linearmodel_path = fullfile(All_dirpath, 'linear', 'RAW');
end
if ~exist(All_linearmodel_path, 'dir')
    mkdir(All_linearmodel_path);
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
    %{
    % {'Fz', 'FCz', 'C3', 'CP3', 'C1'} = {6, 41, 15, 47, 44}
    Fz = 6; FCz = 41; C3 = 15; CP3 = 47; C1 = 44;
    electrodes_name = {'Fz', 'FCz', 'C3', 'CP3', 'C1'};
    electrodes = {Fz; FCz; C3; CP3; C1};
    %}
    electrodes_name = {EEG.chanlocs.labels};
    electrodes = num2cell(1:length(EEG.chanlocs));
    
    nb_epoch = length(EEG.epoch);
    tf_ersp = cell(length(electrodes), 1);
    tf_itc = tf_ersp; tf_powbase = tf_ersp;
    tf_times = tf_ersp; tf_freqs = tf_ersp;
    tf_data = tf_ersp;
    
    % find baseline [-600, -100] before left/right cue
    baseline_b4_leftright = cell(nb_epoch, 1);
    tmp_baseline = [-600, -100]; % in milliseconds
    for i = 1:nb_epoch
        if isempty(find(strcmp([EEG.epoch(i).eventtype], 's17'), 1))
            disp([num2str(subID), ' missing left/right cue in trial ', num2str(i)]);
        end
        baseline_b4_leftright{i} = round((EEG.epoch(i).eventlatency{strcmp([EEG.epoch(i).eventtype], 's17')}) * EEG.srate / 1000) + tmp_baseline; % in millisaconds
    end
    
    ticker = 1;
    h = waitbar(0, 'time frequency analysis.');
    total_iter = nb_epoch * length(electrodes);
    % time frequency analysis for each channel and each epoch
    for i = 1:length(electrodes)
        tmp_ersp = cell(nb_epoch, 1);
        tmp_itc = tmp_ersp; tmp_powbase = tmp_ersp;
        tmp_times = tmp_ersp; tmp_freqs = tmp_ersp;
        tmp_data = tmp_ersp;
        for j = 1:nb_epoch
            [tmp_ersp{j, 1}, tmp_itc{j, 1}, tmp_powbase{j, 1}, tmp_times{j, 1}, tmp_freqs{j, 1}, ~, ~, tmp_data{j, 1}] = ...
            newtimef(tfEEG.data(electrodes{i}, :, j), size(tfEEG.data, 2), [tfEEG.times(1), tfEEG.times(end)], tfEEG.srate, [3, 0.5], ...
                                                      'maxfreq', 35, 'timesout', 800, ...
                                                      'topovec', 15, 'elocs', EEG.chanlocs, 'chaninfo', EEG.chaninfo, 'caption', electrodes_name{i}, ...
                                                      'baseline', baseline_b4_leftright{j}, 'basenorm', 'on', 'trialbase', 'full', 'padratio', 1, 'winsize', 512, ...
                                                      'plotitc' , 'off', 'plotphase', 'off', 'plotersp', 'off');
                                                  
                                                  

            %%%%%% resample tf_* to reduce size
            
                                                  
                                                  
                                                  
            
            ticker = ticker + 1;
            progress_precent = 100 * ticker / total_iter;
            waitbar(ticker / total_iter, h, sprintf('time frequency analysis. %2.2f %%', progress_precent));
        end
        
        % convert
        tf_ersp{i} = cat(3, tmp_ersp{:, 1});
        tf_itc{i} = cat(3, tmp_itc{:, 1});
        tf_powbase{i} = cat(3, tmp_powbase{:, 1});
        tf_times{i} = cat(3, tmp_times{:, 1});
        tf_freqs{i} = cat(3, tmp_freqs{:, 1});
        tf_data{i} = cat(3, tmp_data{:, 1});
    end
    close(h);
    
    tmp_filename = fullfile(All_linearmodel_path, [subID, '_timefreq']);
    save(tmp_filename, 'tf_*')
    
    
    % the continuous variables: roll angle (deg)
    % ****************************************************
    % the roll angle is NOT aligned for Left/Right handles
    % ****************************************************
    Model_Roll = nan(nb_epoch, size(tf_ersp{1, 1}, 2));
    for i = 1:nb_epoch
        % something is weird with the last few rows of EEG.behavior.obj_epoch{i, 1}
        tmp_ind = [true; diff(EEG.behavior.obj_epoch{i, 1}.time) == (1000 / EEG.behavior.behavior_srate)];
        
        % the object roll angles are in radian
        Model_Roll(i, :) = spline(EEG.behavior.obj_epoch{i, 1}.time(tmp_ind), EEG.behavior.obj_epoch{i, 1}.roll(tmp_ind), squeeze(tf_times{1, 1}(:, :, 1)));
        
    end
    % the categorical variables: condition (IL, TR, PT)
    % % % dummy_cond = dummyvar( grp2idx({EEG.epoch.condType}')); % unitary dummies
    
    trialID = str2double({EEG.epoch.trialID})';
    available_trials = intersect(All_beh_expfit.pRoll_dummy(:, 1), trialID);
    dummy_cond = All_beh_expfit.pRoll_dummy(available_trials, 2:end); % exponetial dummies
    
    Model_Cond_IL = dummy_cond(:, 1); Model_Cond_TR = dummy_cond(:, 2); Model_Cond_PT = dummy_cond(:, 3);
    
    % perform regression using power as depend variable; dummy conditions
    % and err as independ variables
    
    % fit = fitlm(ft, 'ERSP ~ Roll * Cond')
    % E(ERSP) = beta_0 * I(IL) + beta_1 * I(TR) + beta_2 * I(PT) + beta_3 * I(IL) * Roll + beta_4 * I(TR) * Roll + beta_5 * I(PT) * Roll
    % each subject results in 8 beta images in the time frequency plane

% % %     Model_Result = cell(length(electrodes), 1);
    Model_coeff_est = cell(length(electrodes), 1);
    Model_coeff_p = Model_coeff_est; Model_Rsquared = Model_coeff_est;
    Model_Cond = categorical( grp2idx({EEG.epoch.condType}') );
    
    ticker = 1;
    h = waitbar(0, 'calculating robust regression coefficients.');
    total_iter = size(tf_ersp{1, 1}, 2) * size(tf_ersp{1, 1}, 1) * length(electrodes);
    for time_id = 1:size(tf_ersp{1, 1}, 2)
        for freq_id = 1:size(tf_ersp{1, 1}, 1)
            for electrode_id = 1:length(electrodes)
                Model_ERSP = squeeze( tf_ersp{electrode_id, 1}(freq_id, time_id, :) );
                tbl = table(Model_Cond_IL, Model_Cond_TR, Model_Cond_PT, Model_Cond, Model_Roll(:, time_id), Model_ERSP);
                tbl.Properties.VariableNames{'Var5'} = 'Model_Roll';
                mdl = fitlm(tbl, 'Model_ERSP ~ Model_Cond_IL + Model_Cond_TR + Model_Cond_PT + Model_Cond_IL:Model_Roll + Model_Cond_TR:Model_Roll + Model_Cond_PT:Model_Roll - 1', 'RobustOpts', 'on');
%                 mdl = fitlm(tbl, 'Model_ERSP ~ Model_Roll * Model_Cond', 'RobustOpts', 'on');
                
% % %                 Model_Result{electrode_id, 1}{freq_id, time_id} = mdl;
                Model_coeff_est{electrode_id, 1}(freq_id, time_id, :) = mdl.Coefficients{:, 'Estimate'};
                Model_coeff_p{electrode_id, 1}(freq_id, time_id, :) = mdl.Coefficients{:, 'pValue'};
                Model_Rsquared{electrode_id, 1}(freq_id, time_id, :) = mdl.Rsquared.Ordinary;
                
                ticker = ticker + 1;
                progress_precent = 100 * ticker / total_iter;
                waitbar(ticker / total_iter, h, sprintf('calculating robust regression coefficients. %2.2f %%', progress_precent));
            end
        end
    end
    close(h);
    
% % %     tmp_filename = fullfile(All_linearmodel_path, [subID, '_LinearModel']);
% % %     save(tmp_filename, 'Model_Result')
% % %     clear Model_Result

    tmp_filename = fullfile(All_linearmodel_path, [subID, '_LinearModel_coeff']);
    save(tmp_filename, 'Model_*')
    tmp_filename = fullfile(All_linearmodel_path, 'misc');
    chanlocs = EEG.chanlocs;
    save(tmp_filename, 'tf_times', 'tf_freqs', 'electrodes_name', 'chanlocs')
    
    disp([All_dirlist(All_i).name, ' finished.']);
end