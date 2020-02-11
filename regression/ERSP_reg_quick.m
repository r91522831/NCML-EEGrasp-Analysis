close all; clear; clc;
All_path = uigetdir;
All_filelist = dir(fullfile(All_path, '*_timefreq.mat'));

%% 
disp([num2cell((1:length(All_filelist))'), {All_filelist.name}']);
All_selected_sub = input('Which subject(s) to run regression? ');
if isempty(All_selected_sub)
    All_selected_sub = 1:length(All_filelist);
end

coeff_name = {'\beta_0', '\beta_1', '\beta_2', '\beta_3', '\beta_4', '\beta_5'};
ncoeff = length(coeff_name);

All_power_roll = cell(length(All_selected_sub), 4);
for All_i = All_selected_sub% 1:length(All_dirlist)
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
    
    rg_time = find(timerstamps >= 400 & timerstamps < 600); % 400 to 600 ms after lift onset
    rg_freq = find(freqz >  4 & freqz <=  8); % 4 to 8 Hz: theta band
    
    % define a window size of ~300 ms moving across the range from 0 to 5000 ms
    rg_all_time = timerstamps(timerstamps >= 0 & timerstamps < 5000);
    rg_ind_dt = round(300 / mean(diff(rg_all_time))); % find number of indices corresponding to 300 ms window
    n_mv_win = length(rg_all_time) - rg_ind_dt;
    rg_300_time = cell(n_mv_win, 1);
    for i_win = 1:n_mv_win
        [~, rg_300_time{i_win, 1}] = ismember(rg_all_time(i_win:(i_win + rg_ind_dt)), timerstamps); % find the indices in the timerstamps
    end
    
    % correlation for each electrode
    for i_electrode = 1:nelectrode
        tmp_peak_roll = abs(table2array(EEG.behavior.obj_roll_peak));
        
        tmp_power_mw = cell(n_mv_win, 2);
        % mean with a window size of ~300 ms moving across the range from 0 to 5000 ms for theata band
        for i_win = 1:n_mv_win
            for i_trial = 1:nb_epoch
                tmp_p_mw = tf.tf_ersp{electrodes(i_electrode), 1}{1}(rg_freq, rg_300_time{i_win, 1}, i_trial);
                tmp_power_mw{i_win, 1}(i_trial, 1) = mean(tmp_p_mw(:));
            end
            tmp_power_mw{i_win, 1}(:, 2) = tmp_peak_roll;
            tmp_r_mw = corr(tmp_power_mw{i_win, 1});
            tmp_power_mw{i_win, 2} = tmp_r_mw(1, 2);
        end
        
        tmp_power = nan(nb_epoch, 2);
        for i_trial = 1:nb_epoch
            % mean with 400 to 600 ms window for theta band (4 to 8 Hz)
            tmp_p = tf.tf_ersp{electrodes(i_electrode), 1}{1}(rg_freq, rg_time, i_trial);
            tmp_power(i_trial, 1) = mean(tmp_p(:));
        end
        tmp_power(:, 2) = tmp_peak_roll;
        tmp_r = corr(tmp_power);
        
        tmp_power = array2table(tmp_power, 'VariableNames', {'power', 'peakRoll'});
        tmp_power{:, 'cond'} = {EEG.epoch.cond}';
        
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
        All_power_roll{All_i, 2} = subID;
        All_power_roll{All_i, 3} = tf.tf_ersp{electrodes(i_electrode), 1}{1};
        All_power_roll{All_i, 4} = tmp_power_mw;
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
sub_pow_roll = cell(size(All_power_roll, 1), 2);
sub_pow_roll(:, 1) = cellfun(@(x) x{:, 1}, All_power_roll(:, 1), 'UniformOutput', false);
sub_pow_roll(:, 2) = All_power_roll(:, 2);
sub_r_theta_400to600ms = cellfun(@(x) x{:, 2}, All_power_roll(:, 1), 'UniformOutput', false);
sub_r_time_freq = cellfun(@(x) x{:, 3}, All_power_roll(:, 1), 'UniformOutput', false);
sub_time = tf.tf_times;
sub_freq = tf.tf_freqs;
save(fullfile(All_path, ['corr_', num2str(length(All_selected_sub)), 'sub.mat']), 'sub_*', 'All_power_roll');


%% plot TR1/IL1 peakRoll
for i = 1:length(All_power_roll)
    % find TR1
    i_tr1 = find(strcmp([All_power_roll{1, 1}{1, 1}{:, 'cond'}], 'TR'), 1);
    % find IL1
    i_il1 = find(strcmp([All_power_roll{1, 1}{1, 1}{:, 'cond'}], 'IL'), 1);
    %
    All_power_roll{i, 5} = All_power_roll{i, 1}{1, 1}{i_tr1, 'peakRoll'} / All_power_roll{i, 1}{1, 1}{i_il1, 'peakRoll'};
end
All_power_roll = sortrows(All_power_roll, 5, 'descend');

%%
figure('units', 'normalized', 'outerposition', [0, 0, 1, 1]);
set(0, 'defaultAxesFontSize', 18)
plot(1:size(All_power_roll, 1), [All_power_roll{:, 5}], 'x')
hline([3, 1], {'r', 'b'}, {'3', '1'})
xticks(1:size(All_power_roll, 1))
xticklabels(cellfun(@(x) x(end-1:end), All_power_roll(:, 2), 'UniformOutput', false))
xlabel('Suject ID')
ylabel('ratio (TR_1/IL_1)')
title('pealRoll_{TR_1} / peakRoll_{IL_1} ')




%% power x peak roll plot
figure('units', 'normalized', 'outerposition', [0, 0, 1, 1])
set(0, 'defaultAxesFontSize', 18)
for i = 1:size(sub_pow_roll, 1)
    tmp_IL = strcmp(sub_pow_roll{i, 1}{:, 'cond'}, 'IL');
    tmp_TR = strcmp(sub_pow_roll{i, 1}{:, 'cond'}, 'TR');
    tmp_PT1 = strcmp(sub_pow_roll{i, 1}{:, 'cond'}, 'PT1');
    tmp_PT2 = strcmp(sub_pow_roll{i, 1}{:, 'cond'}, 'PT2');
    tmp_PT3 = strcmp(sub_pow_roll{i, 1}{:, 'cond'}, 'PT3');
    subplot(4, 5, i)
    hold on
    scatter(sub_pow_roll{i, 1}{tmp_IL, 'peakRoll'}, sub_pow_roll{i, 1}{tmp_IL, 'power'}, 'or')
    scatter(sub_pow_roll{i, 1}{tmp_TR, 'peakRoll'}, sub_pow_roll{i, 1}{tmp_TR, 'power'}, 'xb')
    scatter(sub_pow_roll{i, 1}{tmp_PT1, 'peakRoll'}, sub_pow_roll{i, 1}{tmp_PT1, 'power'}, '*k')
    scatter(sub_pow_roll{i, 1}{tmp_PT2, 'peakRoll'}, sub_pow_roll{i, 1}{tmp_PT2, 'power'}, '+y')
    scatter(sub_pow_roll{i, 1}{tmp_PT3, 'peakRoll'}, sub_pow_roll{i, 1}{tmp_PT3, 'power'}, 'xc')
    hold off
    
    xlabel('peak roll ({\circ})')
    ylabel('tf power')
    title(sub_pow_roll{i, 2})
end
legend('IL', 'TR', 'PT1', 'PT2', 'PT3')


savefig(fullfile(All_path, ['roll_power_', num2str(length(All_selected_sub)), 'sub.fig']));

%%
All_r = nan(length(All_power_roll), 1);
for i = 1:size(All_power_roll, 1)
    All_r(i) = All_power_roll{i, 1}{1, 2};
end
disp(mean(All_r))

%% correlation plot
h = figure('units', 'normalized', 'outerposition', [0, 0, 1, 1]);
set(0, 'defaultAxesFontSize', 18)
for i = 1:size(All_power_roll, 1)
    subplot(4, 5, i)
    contourf(sub_time(sub_time >= -1)/1000, sub_freq, All_power_roll{i, 1}{1, 3}(:, sub_time >= -1), 40, 'LineStyle', 'none')
end
mtit('Correlation r')
% % % savefig(h, fullfile(All_path, ['corr_r_', num2str(length(All_selected_sub)), 'sub.fig']));

%% power plot
h = figure('units', 'normalized', 'outerposition', [0, 0, 1, 1]);
set(0, 'defaultAxesFontSize', 18)
for i = 1:size(All_power_roll, 1)
    tmp_p = trimmean(All_power_roll{i, 3}, 10, 3);
    subplot(4, 5, i)
    contourf(sub_time(sub_time >= -1)/1000, sub_freq, tmp_p(:, sub_time >= -1), 40, 'LineStyle', 'none')
end
mtit('Power')
% % % savefig(h, fullfile(All_path, ['raw_p_', num2str(length(All_selected_sub)), 'sub.fig']));




%% plot correlation time trjectory
figure
hold on
for i = 1:length(All_power_roll)
    plot([All_power_roll{i, 4}{:, 2}])
    xlabel('moveing window (width ~300ms)')
    ylabel('corr r')
    title('Mean power of \theta 4~8 Hz moving window ~300 ms vs peark roll corr r time trajectory')
end
