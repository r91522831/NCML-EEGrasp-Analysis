close; clear; clc;

data_dir = uigetdir;
behavior_dir = fullfile(fileparts(fileparts(fileparts(fileparts(data_dir)))), 'behavior', 'matlab data');
% % % output_dir = fullfile(fileparts(data_dir), 'rawVoltage_trial_chanTime');
output_dir = fullfile(fileparts(data_dir), 'rawVoltage_trialTime_chan');
if ~isfolder(output_dir)
    mkdir(output_dir);
end
file_list = dir(fullfile(data_dir, '*_pruned_ICA.set'));
sub_id = cellfun(@(x) x(1:4), {file_list.name}', 'UniformOutput', false);
nb_sub = length(sub_id);

cond_names = {'IL', 'TR', 'PT1', 'PT2', 'PT3'};
nb_cond = length(cond_names);

%% get all raw data across subjects
sub_data = cell(nb_sub, 4); % data, time, srate, trial type
for sub_i = 1:nb_sub
    [ALLEEG, ~, ~, ALLCOM] = eeglab;
    EEG = pop_loadset('filename', file_list(sub_i).name, 'filepath', data_dir);
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );

    sub_data{sub_i, 1} = EEG.data;
    sub_data{sub_i, 2} = EEG.times;
    sub_data{sub_i, 3} = EEG.srate;
    
    % get condition for trials
    tmp_cond_list = dir( fullfile(behavior_dir, sub_id{sub_i}, '*.csv') );
    tmp_cond_id = cell(size(tmp_cond_list, 1), 1);
    for j = 1:length(tmp_cond_list)
        if strcmp(tmp_cond_list(j).name(11:12), 'PT')
            tmp_cond_id{j, 1} = tmp_cond_list(j).name([11:12, 18]);
        else
            tmp_cond_id{j, 1} = tmp_cond_list(j).name(11:12);
        end
    end
    sub_data{sub_i, 4} = tmp_cond_id;
end
%% Trim epoch across subjects
ind_time(:, 1) = cellfun(@(x) find((x == 0)), sub_data(:, 2));
ind_time(:, 2) = cellfun(@(x) length(x) - find((x == 0)), sub_data(:, 2));
for sub_i = 1:nb_sub
    sub_data{sub_i, 1} = sub_data{sub_i, 1}(:, (ind_time(sub_i, 1) - min(ind_time(:, 1)) + 1):(ind_time(sub_i, 1) + min(ind_time(:, 2))), :);
    sub_data{sub_i, 2} = sub_data{sub_i, 2}(1, (ind_time(sub_i, 1) - min(ind_time(:, 1)) + 1):(ind_time(sub_i, 1) + min(ind_time(:, 2))));
end
%% Downsample data
d_factor = 8;
for sub_i = 1:nb_sub
    tmp_data = sub_data{sub_i, 1};
    tmp_dnsampled = nan(size(tmp_data, 1), ceil(size(tmp_data, 2) / d_factor), size(tmp_data, 3));
    for i = 1:size(tmp_data, 1)
        for j = 1:size(tmp_data, 3)
            tmp_dnsampled(i, :, j) = decimate(tmp_data(i, :, j), d_factor);
        end
    end
    sub_data{sub_i, 1} = tmp_dnsampled;
    sub_data{sub_i, 2} = decimate(sub_data{sub_i, 2}, d_factor);
    sub_data{sub_i, 3} = ceil(sub_data{sub_i, 3} / d_factor);
end


%% Constrained Principal Component Analysis
z_sub = cell(nb_sub, 1);
g_sub = cell(nb_sub, 1);
for sub_i = 1:nb_sub
    ind_t0 = find(sub_data{sub_i, 2} == 0);
    dt = 1000 / sub_data{sub_i, 3}; % in ms
    ind_win_half_width = 3000 / dt; % -3000 ~ 3000 ms
    while (ind_t0 - ind_win_half_width) <= 0
        ind_win_half_width = str2double(input('Set time window half width in ms: ', 's')) / dt;
    end
    ind_win = (ind_t0 - ind_win_half_width):(ind_t0 + ind_win_half_width);
    extracted_data = sub_data{sub_i, 1}(:, ind_win, :);    
    
    nb_trial = length(sub_data{sub_i, 4});
    
    % Z and G matrix
    % z_sub: epoch x (time x channel)
    %{
    %{
    %                    channel 1                        |                    channel 2                        | ...
    %          time 1       |          time 2       | ... |          time 1       |          time 2       | ... | ...
    %}
    z_sub{sub_i} = reshape(permute(sub_data{sub_i, 1}, [3, 2, 1]), nb_trial, []);

    % G matrix
    g_sub{sub_i, 1} = nan(nb_trial, nb_cond);
    for j = 1:nb_cond
        g_sub{sub_i, 1}(:, j) = strcmp(sub_data{sub_i, 4}, cond_names{j})';
    end
    %}
    
    ticks_time = sub_data{sub_i, 2};
    % z_sub: (epoch x time) x channel
    %{
    %                             channel 1                        |                    channel 2                        | ...
    %          time1
    % epoch1   time2
    %}
    z_sub{sub_i} = reshape(permute(sub_data{sub_i, 1}, [2, 3, 1]), nb_trial * length(ticks_time), []);
    
    % G matrix
    tmp_sub = false(nb_trial, nb_cond);
    for j = 1:nb_cond
        tmp_sub(:, j) = strcmp(sub_data{sub_i, 4}, cond_names{j})';
    end
    tmp_g_sub = cell(nb_trial, nb_cond);
    for j = 1:nb_cond
        for k = 1:nb_trial
            if tmp_sub(k, j)
                tmp_g_sub{k, j} = ones(length(ticks_time), 1);
            else
                tmp_g_sub{k, j} = zeros(length(ticks_time), 1);
            end
        end
    end
    g_sub{sub_i, 1} = cell2mat(tmp_g_sub);
    %}
end

%% Combine Z and G across subjects
% combine the huge Z
ALL_z = cat(1, z_sub{:});
% combine the huge G
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

bin_time = length(sub_data{1, 2});
bin_chan = size(sub_data{1, 1}, 1);
ticks_times = sub_data{1, 2};

save(fullfile(output_dir, 'ALL_z'), 'ALL_z', '-v7.3');
save(fullfile(output_dir, 'ALL_misc'), 'ticks_times', 'sub_data', 'bin_*', '-v7.3');
save(fullfile(output_dir, 'ALL_g'), 'ALL_g', '-v7.3');