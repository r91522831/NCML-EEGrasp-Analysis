close; clear; clc;

data_dir = uigetdir;
behavior_dir = fullfile(fileparts(fileparts(fileparts(data_dir))), 'behavior', 'matlab data');
file_list = dir(fullfile(data_dir, '*_tf_info.mat'));
sub_id = cellfun(@(x) x(1:4), {file_list.name}', 'UniformOutput', false);
nb_sub = length(sub_id);

cond_names = {'IL', 'TR', 'PT1', 'PT2', 'PT3'};
nb_cond = length(cond_names);

%% Constrained Principal Component Analysis
z_sub = cell(nb_sub, 1);
z_sub_theta = z_sub;
z_sub_alpha = z_sub;
z_sub_beta = z_sub;
g_sub = cell(nb_sub, 1);
dimensions = nan(nb_sub, 4);
for i = 1:nb_sub
    clear tf_* tmp_*
    load(fullfile(data_dir, file_list(i).name));
    
    % tf_data{j, 1}: f x t x epoch; cat(4, tf_data{:, 1}): f x t x epoch x channel;
    % tf_freqs: freq ticks; tf_times: time ticks
    tf_nb_trial = size(tf_data{1, 1}, 3);
    tf_nb_cond = size(tf_data, 2) - 1; % [All, IL, TR, PT1, PT2, PT3]
    tf_all4d = cat(4, tf_data{:, 1});
    
    % get condition for trials
    tmp_cond_list = dir( fullfile(behavior_dir, sub_id{i}, '*.csv') );
    tmp_cond_id = cell(tf_nb_trial, 1);
    for j = 1:length(tmp_cond_list)
        if strcmp(tmp_cond_list(j).name(11:12), 'PT')
            tmp_cond_id{j, 1} = tmp_cond_list(j).name([11:12, 18]);
        else
            tmp_cond_id{j, 1} = tmp_cond_list(j).name(11:12);
        end
    end
    
    % Get indices for 4-8 Hz (theta); 9-12 Hz (alpha); 15-30 Hz (beta)
    ind_theta = tf_freqs >= 4 & tf_freqs <= 8;
    ind_alpha = tf_freqs >= 9 & tf_freqs <= 12;
    ind_beta = tf_freqs >= 13 & tf_freqs <= 30;
    
    % Get indices for [-1000 onset 1000] ms
    ind_t2s = find(tf_times >= -1000 & tf_times <= 1000);
    donwsample = 1;
    
    % average across freq bands
    %{
    tf_bands_theta = abs(mean(tf_all4d(ind_theta, ind_t2s, :, :), 1));
    tf_bands_alpha = abs(mean(tf_all4d(ind_alpha, ind_t2s, :, :), 1));
    tf_bands_beta = abs(mean(tf_all4d(ind_beta, ind_t2s, :, :), 1));
    %}
    % get all without averaging (too large Z matrix)
    tf_bands_theta = tf_all4d(ind_theta, ind_t2s(1:donwsample:end), :, :);
    tf_bands_alpha = tf_all4d(ind_alpha, ind_t2s(1:donwsample:end), :, :);
    tf_bands_beta = tf_all4d(ind_beta, ind_t2s(1:donwsample:end), :, :);
    
    tf_bands_freqs = [tf_bands_theta; tf_bands_alpha; tf_bands_beta]; % truncate 3 Hz and > 30 Hz data
% % %     tf_bands_freqs = tf_all4d(:, ind_t2s(1:donwsample:end), :, :);
    
    dimensions(i, :) = size(tf_bands_freqs);
    ticks_theta = tf_freqs(ind_theta);
    ticks_alpha = tf_freqs(ind_alpha);
    ticks_beta = tf_freqs(ind_beta);
    ticks_time = tf_times(ind_t2s(1:donwsample:end));
    
    %% decimate freqs and times
    %{
    % decimate across 11 freq bins change freq resolution from .5 Hz to 5 Hz
    tf_freq_band = cell(length(tf_data), 1);
    dn_factor_freq = 11;
    for i = 1:length(tf_data)
        [fband, t, ch] = size(tf_data{i, 1});
        dn_f = ceil(fband / dn_factor_freq);
        tmp = nan(dn_f, t, ch);
        for j = 1:t
            for k = 1:ch
                tmp(:, j, k) = decimate(tf_data{i, 1}(:, j, k), dn_factor_freq);
            end
        end
        tf_freq_band{i, 1} = tmp;
    end


    % decimate across 5 time bins ~= 150~200 ms
    tf_ftbands = cell(length(tf_data), 1);
    dn_factor = 5;
    for i = 1:length(tf_freq_band)
        [fband, t, ch] = size(tf_freq_band{i, 1});
        dn_t = ceil(t / dn_factor);
        tmp = nan(fband, dn_t, ch);
        for j = 1:fband
            for k = 1:ch
                tmp(j, :, k) = decimate(tf_freq_band{i, 1}(j, :, k), dn_factor);
            end
        end
        tf_ftbands{i, 1} = tmp;
    end
    %}

    
    %% Z and G matrix
    % z_sub: epoch x (time x freq x channel)
    %{
    %{
    %                    channel 1                        |                    channel 2                        | ...
    %          freq 1       |          freq 2       | ... |          freq 1       |          freq 2       | ... | ...
    % time 1 | time 2 | ... | time 1 | time 2 | ... | ... | time 1 | time 2 | ... | time 1 | time 2 | ... | ... | ...
    %}
    
    % % % z_sub = reshape(permute(cat(4, tf_data{:, 1}), [3, 2, 1, 4]), EEG.trials, []);
    % % % z_sub = reshape(permute(cat(4, tf_freq_band{:, 1}), [3, 2, 1, 4]), EEG.trials, []);
    % % % z_sub = reshape(permute(cat(4, tf_ftbands{:, 1}), [3, 2, 1, 4]), EEG.trials, []);
    
    z_sub{i} = reshape(permute(tf_bands_freqs, [3, 2, 1, 4]), tf_nb_trial, []);
    z_sub_theta{i} = reshape(permute(tf_bands_theta, [3, 2, 1, 4]), tf_nb_trial, []);
    z_sub_alpha{i} = reshape(permute(tf_bands_alpha, [3, 2, 1, 4]), tf_nb_trial, []);
    z_sub_beta{i} = reshape(permute(tf_bands_beta, [3, 2, 1, 4]), tf_nb_trial, []);
    
    % G matrix
    g_sub{i, 1} = nan(tf_nb_trial, nb_cond);
    for j = 1:nb_cond
        g_sub{i, 1}(:, j) = strcmp(tmp_cond_id, cond_names{j})';
    end
    %}
    
    %%
    % z_sub: (epoch x time) x (freq x channel)
    %{
    %                             channel 1                        |                    channel 2                        | ...
    %                   freq 1       |          freq 2       | ... |          freq 1       |          freq 2       | ... | ...
    %          time1
    % epoch1   time2
    %}
    z_sub{i} = reshape(permute(tf_bands_freqs, [2, 3, 1, 4]), tf_nb_trial * length(ticks_time), []);
    z_sub_theta{i} = reshape(permute(tf_bands_theta, [2, 3, 1, 4]), tf_nb_trial * length(ticks_time), []);
    z_sub_alpha{i} = reshape(permute(tf_bands_alpha, [2, 3, 1, 4]), tf_nb_trial * length(ticks_time), []);
    z_sub_beta{i} = reshape(permute(tf_bands_beta, [2, 3, 1, 4]), tf_nb_trial * length(ticks_time), []);
    
    % G matrix
    tmp_sub = false(tf_nb_trial, nb_cond);
    for j = 1:nb_cond
        tmp_sub(:, j) = strcmp(tmp_cond_id, cond_names{j});
    end
    tmp_g_sub = cell(tf_nb_trial, nb_cond);
    for j = 1:nb_cond
        for k = 1:tf_nb_trial
            if tmp_sub(k, j)
                tmp_g_sub{k, j} = ones(length(ticks_time), 1);
            else
                tmp_g_sub{k, j} = zeros(length(ticks_time), 1);
            end
        end
    end
    g_sub{i, 1} = cell2mat(tmp_g_sub);
end

%% combine the huge Z
ALL_z = cat(1, z_sub{:});
ALL_z_theta = cat(1, z_sub_theta{:});
ALL_z_alpha = cat(1, z_sub_alpha{:});
ALL_z_beta = cat(1, z_sub_beta{:});
%% combine the huge G
ALL_g = cell(nb_sub);
for i = 1:nb_sub
    for j = 1:nb_sub
        if i ~= j
            ALL_g{i, j} = zeros(size(g_sub{i, 1}));
        else
            ALL_g{i, j} = g_sub{i, 1};
        end
    end
end
ALL_g = cell2mat(ALL_g);

dimensions = array2table(dimensions, 'VariableNames', {'FreqBands', 'TimeSteps', 'epochs', 'Channels'});
bin_time = dimensions{1, 'TimeSteps'};
bin_freq = dimensions{1, 'FreqBands'};
bin_chan = dimensions{1, 'Channels'};

save(fullfile(data_dir, 'ALL_z'), 'ALL_z', '-v7.3');
save(fullfile(data_dir, 'ALL_z_theta'), 'ALL_z_theta', '-v7.3');
save(fullfile(data_dir, 'ALL_z_alpha'), 'ALL_z_alpha', '-v7.3');
save(fullfile(data_dir, 'ALL_z_beta'), 'ALL_z_beta', '-v7.3');
save(fullfile(data_dir, 'ALL_misc'), 'dimensions', 'ticks_*', 'bin_freq', 'bin_time', 'bin_chan', '-v7.3');
save(fullfile(data_dir, 'ALL_g'), 'ALL_g', '-v7.3');