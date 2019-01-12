clear; clc;

cond_names = {'ALL', 'IL', 'TR', 'PT1', 'PT2', 'PT3'};
cond_nb = length(cond_names);

base_folder = uigetdir('Select base folder');
tf_data_list = dir(fullfile(base_folder, '*_tf_info.mat'));

for sub_i = 1:length(tf_data_list)
    subID = tf_data_list(sub_i).name(1:4);
    
    EEG = pop_loadset('filename', [subID, '_pruned_ICA.set'], 'filepath', fullfile(base_folder, subID));
    
    load(fullfile(base_folder, tf_data_list(sub_i).name));
    
    % tf_data{:, 1}: f x t x epoch; cat(4, tf_data{:, 1}): f x t x epoch x channel;
    % tf_freqs: freq ticks; tf_times: time ticks
    
    % Get 4-8 Hz (theta); 9-12 Hz (alpha); 20-30 Hz (beta)
    ind_theta = tf_freqs >= 4 & tf_freqs <= 8;
    ind_alpha = tf_freqs >= 9 & tf_freqs <= 12;
    ind_beta = tf_freqs >= 20 & tf_freqs <= 30;
    
    %{
    % average across freq bands
    tf_freq_band = cell(length(tf_data), 1);
    dn_f = 3;
    for i = 1:length(tf_data)
        tf_freq_band{i, 1} = cat(1, mean(tf_data{i, 1}(ind_alpha, :, :), 1), mean(tf_data{i, 1}(ind_beta, :, :), 1), mean(tf_data{i, 1}(ind_theta, :, :), 1));
    end
    %}
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
    
    % z_sub: epoch x (time x freq x channel)
    % z_sub = reshape(permute(cat(4, tf_data{:, 1}), [3, 2, 1, 4]), EEG.trials, []);
    % z_sub = reshape(permute(cat(4, tf_freq_band{:, 1}), [3, 2, 1, 4]), EEG.trials, []);
    z_sub = reshape(permute(cat(4, tf_ftbands{:, 1}), [3, 2, 1, 4]), EEG.trials, []);
    
    % G matrix
    g_sub = nan(EEG.trials, length(cond_names) - 1);
    for i = 2:length(cond_names)
        g_sub(:, i - 1) = strcmp([EEG.epoch.cond], cond_names{i})';
    end
    
    %%
    % H matrix
    % % % bin_time = length(tf_times); % without downsample in time
    bin_time = dn_t;
    % % % bin_freq = length(tf_freqs); % without downsample in freq
    bin_freq = dn_f;
    bin_chan = EEG.nbchan;
    
    nb_tf = bin_time * bin_freq;

    %% construct h for channel and freq
    % construct ones for time bins
    h_time = ones(bin_time, 1);
    % construct ones for freq bins
    h_time_freq = cell(bin_freq);
    for i = 1:bin_freq
        for j = 1:bin_freq
            if i == j
                h_time_freq{i, j} = h_time;
            else
                h_time_freq{i, j} = zeros(size(h_time));
            end
        end
    end
    h_time_freq = cell2mat(h_time_freq);
    % construct h for channel and freq
    h_sub = cell(bin_chan); % channel x channel
    for i = 1:bin_chan
        for j = 1:bin_chan
            if i == j
                h_sub{i, j} = h_time_freq;
            else
                h_sub{i, j} = zeros(size(h_time_freq));
            end
        end
    end
    h_sub = cell2mat(h_sub);
    
    %% NG H
    %{
    % construct h for time bins
    h_time = eye(bin_time, bin_time);
    % construct h for freq and time bins
    h_freq = zeros(bin_time * bin_freq, bin_freq + bin_time);
    for i = 1:bin_freq
        tmp_freq = zeros(bin_time, bin_freq);
        tmp_freq(1:bin_time, i) = 1;
        
        h_freq((bin_time * (i - 1) + 1):(bin_time * i), :) = [tmp_freq, h_time];
    end
    % construct h for channel, freq, and time
    q_sub = bin_chan + bin_freq + bin_time;
    h_sub = zeros(size(z_sub, 2), q_sub); % (time * freq * channel) x (channel * freq * time)
    for i = 1:bin_chan
        tmp_sub = zeros(bin_freq * bin_time, bin_chan);
        tmp_sub(1:nb_tf, i) = 1;
        
        h_sub((nb_tf * (i - 1) + 1):(nb_tf * i), :) = [tmp_sub, h_freq];
    end
    %}
    
    %% construct h for channel
    %{
    h_sub = cell(bin_chan); % channel x channel
    for i = 1:bin_chan
        for j = 1:bin_chan
            if i == j
                h_sub{i, j} = ones(nb_tf, 1); % ones(freq x time, 1)
            else
                h_sub{i, j} = zeros(nb_tf, 1);
            end
        end
    end
    h_sub = cell2mat(h_sub);
    %}
    
    %% construct h for freq
    %{
    % construct ones for time bins
    h_time = ones(bin_time, 1);
    % construct ones for freq bins
    h_time_freq = cell(bin_freq);
    for i = 1:bin_freq
        for j = 1:bin_freq
            if i == j
                h_time_freq{i, j} = h_time;
            else
                h_time_freq{i, j} = zeros(size(h_time));
            end
        end
    end
    h_time_freq = cell2mat(h_time_freq);
    
    h_sub = cell(bin_chan, 1); % channel x 1
    for i = 1:bin_chan
        h_sub{i, 1} = h_time_freq; % ones(freq x time, 1)
    end
    h_sub = cell2mat(h_sub);
    %}
    
    %% construct h for time
    %{
    % construct ones for time bins
    h_time = eye(bin_time);
    h_sub = repmat(h_time, bin_freq * bin_chan, 1);
    %}
    
    %% construct h for channel and time
    %{
    % construct ones for time bins
    h_time = eye(bin_time, bin_time);
    % construct ones for freq and time bins
    h_time_freq = repmat(h_time, bin_freq, 1);
    % construct ones for channel and time
    h_sub = cell(bin_chan); % channel x channel
    for i = 1:bin_chan
        for j = 1:bin_chan
            if i == j
                h_sub{i, j} = h_time_freq;
            else
                h_sub{i, j} = zeros(size(h_time_freq));
            end
        end
    end
    h_sub = cell2mat(h_sub);
    %}
    
    %%
    save(fullfile(base_folder, [subID, '_z_sub']), 'z_sub', '-v7.3');
    save(fullfile(base_folder, [subID, '_g_sub']), 'g_sub', '-v7.3');
    save(fullfile(base_folder, [subID, '_h_sub']), 'h_sub', 'bin_freq', 'bin_time', 'bin_chan', '-v7.3');
end