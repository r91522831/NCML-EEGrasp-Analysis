close all; clear; clc;
All_path = uigetdir;
All_filelist = dir(fullfile(All_path, '*_timefreq.mat'));

%% 
disp([num2cell((1:length(All_filelist))'), {All_filelist.name}']);
selected_sub = input('Which subject(s) to run regression? ');
if isempty(selected_sub)
    selected_sub = 1:length(All_filelist);
end

coeff_name = {'\beta_0', '\beta_1', '\beta_2', '\beta_3', '\beta_4', '\beta_5'};
ncoeff = length(coeff_name);

All_power_roll = cell(length(selected_sub), 1);
for All_i = selected_sub% 1:length(All_dirlist)
    clearvars -except All_*; close all;
    subID = All_filelist(All_i).name(1:6);
    filepath_eeg = dir(fullfile(fileparts(fileparts(All_path)), [subID, '*']));
    filename_eeg = dir(fullfile(filepath_eeg.folder, filepath_eeg.name, '*eeg_csd.set'));
    % load eeg data set
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
    EEG = pop_loadset('filename', filename_eeg.name, 'filepath', filename_eeg.folder);
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
    
    % load time frequency data set
% % %     tf = load('/Users/yenhsunw/Dropbox (ASU)/NCML-EEGrasp/EEG/eeglab/002 ProcessAgain/linear/RAW/sub-09_timefreq.mat');
    tf = load(fullfile(All_path, All_filelist(All_i).name));

    timerstamps = tf.tf_times';
    freqz = tf.tf_freqs';
    electrodes = tf.tf_ersp.Properties.RowNames;

    ntime = length(timerstamps);
    nfreq = length(freqz);
    nelectrode = length(electrodes);
    nb_epoch = size(tf.tf_ersp{1, 1}{:}, 3);
    
    rg_time = find(timerstamps >= 400 & timerstamps < 600);
    rg_freq = find(freqz >  4 & freqz <=  8);
    
    for i_electrode = 1:nelectrode
        tmp_power = nan(nb_epoch, 2);
        for i_trial = 1:nb_epoch
            tmp_p = tf.tf_ersp{electrodes(i_electrode), 1}{1}(rg_freq, rg_time, i_trial);
            tmp_power(i_trial, 1) = mean(tmp_p(:));
        end
        tmp_power(:, 2) = table2array(EEG.behavior.obj_roll_peak);
        tmp_r = corr(tmp_power);
        
        All_power_roll{All_i, 1}{i_electrode, 1} = tmp_power;
        All_power_roll{All_i, 1}{i_electrode, 2} = tmp_r(1, 2);
        
        
        tmp_corr = nan(nfreq, ntime);
        for i_f = 1:nfreq
            for i_t = 1:ntime
                tmp_pow_r = corr([squeeze(tf.tf_ersp{electrodes(i_electrode), 1}{1}(i_f, i_t, :)), table2array(EEG.behavior.obj_roll_peak)]);
                tmp_corr(i_f, i_t) = tmp_pow_r(1, 2);
            end
        end
        
        All_power_roll{All_i, 1}{i_electrode, 3} = tmp_corr;
    end
    
    
    
    
    
    %{
    
    % delta: 2 ~ 3 Hz; theta: 4 ~ 7 Hz; alpha: 8 ~ 12 Hz; low beta: 13 ~ 19 Hz; high beta: 20 ~ 29 Hz; low gama: 30 ~ 34 Hz
    rg_freq_band = { {'\delta', '2-4 Hz'}, {'\theta', '5-8 Hz'}, {'\alpha', '9-12 Hz'}, {'\beta_{low}', '13-20 Hz'}, {'\beta_{high}', '21-30 Hz'}, {'\gamma_{low}', '31-35 Hz'}; ...
                     find(freqz >= 2 & freqz <=  4), find(freqz >  4 & freqz <=  8), ...
                     find(freqz >  8   & freqz <= 13), find(freqz > 13 & freqz <= 20), ...
                     find(freqz > 20   & freqz <= 30), find(freqz > 30 & freqz <= 35) };
    rg_time_win = { {'-50 to 150 ms'}, {'150 to 350 ms'}, {'350 to 450 ms'}, {'450 to 650 ms'}, {'650 to 850 ms'}; ...
                    find(timerstamps >= -50 & timerstamps < 150), find(timerstamps >= 150 & timerstamps < 350), ...
                    find(timerstamps >= 350 & timerstamps < 450), find(timerstamps >= 450 & timerstamps < 650), ...
                    find(timerstamps >= 650 & timerstamps < 850) };
                
                
                
    
    % the continuous variables: roll angle (deg)
    % ****************************************************
    % the roll angle is NOT aligned for Left/Right handles
    % ****************************************************
    Model_Roll = nan(nb_epoch, ntime);
    for i = 1:nb_epoch
        % something is weird with the last few rows of EEG.behavior.obj_epoch{i, 1}
        tmp_ind = [true; diff(EEG.behavior.obj_roll_epoch{i, 1}.time) == (1000 / EEG.behavior.behavior_srate)];
        
        % the object roll angles are in radian
        Model_Roll(i, :) = spline(EEG.behavior.obj_roll_epoch{i, 1}.time(tmp_ind), EEG.behavior.obj_roll_epoch{i, 1}.Roll(tmp_ind), timerstamps);
    end
    
    
    % the categorical variables: condition (IL, TR, PT)
    dummy_cond = dummyvar( grp2idx({EEG.epoch.condType}')); % unitary dummies
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
    
    %}
    
    
    
    
    
    
    
    
    
    
    %{
    
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
    
    %}
    
    
    
    disp([subID, ' finished.']);
end

%%
All_r = nan(length(All_power_roll), 1);
for i = 1:length(All_power_roll)
    All_r(i) = All_power_roll{i, 1}{1, 2};
end
disp(mean(All_r))

%%
figure
% plot(All_power_roll{}


