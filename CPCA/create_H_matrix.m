close; clear; clc;

data_dir = uigetdir;
file_list = dir(fullfile(data_dir, 'ALL_*.mat'));

for i = 1:length(file_list)
    if any(strcmp(file_list(i).name, {'ALL_misc.mat'}))
        load(fullfile(data_dir, file_list(i).name));
    end
end

%%
% H matrix
% % % dn = dimensions(1, :);
% % % bin_time = length(tf_times); % without downsample in time
% % % bin_freq = length(tf_freqs); % without downsample in freq
% % % bin_time = dn{1, 'TimeSteps'};
% % % bin_freq = dn{1, 'FreqBands'};
% % % bin_chan = dn{1, 'Channels'};

% % % description = [];
addition_description = [];
% % % bin_freq = 1;
% % % addition_description = '_singleFband';

if exist('bin_freq', 'var')
    nb_tf = bin_time * bin_freq;
else
    nb_tf = bin_time;
end

% construct h for 3 bands x channel (freqs x channel)
%{
if bin_freq == length(ticks_theta) + length(ticks_alpha) + length(ticks_beta)
    iof_ticks = [{ticks_theta}; {ticks_alpha}; {ticks_beta}];
    h_theta_alpha_beta = cell(3, 3);
    for i = 1:3
        for j = 1:3
            if i == j
                h_theta_alpha_beta{i, j} = ones(length(iof_ticks{i}), 1);
            else
                h_theta_alpha_beta{i, j} = zeros(length(iof_ticks{i}), 1);
            end
        end
    end
    h_theta_alpha_beta = cell2mat(h_theta_alpha_beta);
else
    h_theta_alpha_beta = zeros(bin_freq, 3);
    h_theta_alpha_beta(ticks_theta, 1) = 1;
    h_theta_alpha_beta(ticks_alpha, 2) = 1;
    h_theta_alpha_beta(ticks_beta, 3) = 1;
end

ALL_h = cell(bin_chan); % channel x channel
for i = 1:bin_chan
    for j = 1:bin_chan
        if i == j
            ALL_h{i, j} = h_theta_alpha_beta;
        else
            ALL_h{i, j} = zeros(size(h_theta_alpha_beta));
        end
    end
end
ALL_h = sparse(full(cell2mat(ALL_h))); % convert back to full matrix and then convert to sparse matrix can save some workspace run-time memory

description = '3bandChan';
%}

% construct h for channel (freqs x channel)
%{
ALL_h = cell(bin_chan); % channel x channel
for i = 1:bin_chan
    for j = 1:bin_chan
        if i == j
            ALL_h{i, j} = ones(bin_freq, 1);
        else
            ALL_h{i, j} = zeros(bin_freq, 1);
        end
    end
end
ALL_h = sparse(full(cell2mat(ALL_h))); % convert back to full matrix and then convert to sparse matrix can save some workspace run-time memory

description = 'chan';
%}






% construct h for single freq band (theta, alpha, beta) as time x channel
%{
% % % nb_f = length(ticks_theta);
% % % nb_f = length(ticks_alpha);
nb_f = length(ticks_beta);

h_time = speye(bin_time);
h_time_fband = repmat(h_time, nb_f, 1);
ALL_h = cell(bin_chan); % channel x channel
for i = 1:bin_chan
    for j = 1:bin_chan
        if i == j
            ALL_h{i, j} = h_time_fband;
        else
            ALL_h{i, j} = zeros(size(h_time_fband));
        end
    end
end
ALL_h = sparse(full(cell2mat(ALL_h))); % convert back to full matrix and then convert to sparse matrix can save some workspace run-time memory

% % % description = 'TimeChan_theta';
% % % description = 'TimeChan_alpha';
description = 'TimeChan_beta';
%}


% construct h for channel
%{
ALL_h = cell(bin_chan); % channel x channel
for i = 1:bin_chan
    for j = 1:bin_chan
        if i == j
            ALL_h{i, j} = ones(nb_tf, 1); % ones(freq x time, 1)
        else
            ALL_h{i, j} = zeros(nb_tf, 1);
        end
    end
end
ALL_h = cell2mat(ALL_h);
description = 'chan4rawV_chanTime';
%}

% construct h for channel
%{
ALL_h = eye(bin_chan); % channel x channel
description = 'chan';
%}

% construct h for freq bands (theta, alpha, beta)

% construct ones for freq bands with different freq bins: length(theta), length(alpha), length(beta)
h_time_fband = {ones(bin_time * length(ticks_theta), 1), ones(bin_time * length(ticks_alpha), 1), ones(bin_time * length(ticks_beta), 1)};
bin_freq_band = 3;
h_time_freq = cell(3);
for i = 1:bin_freq_band
    for j = 1:bin_freq_band
        if i ~= j
            h_time_freq{i, j} = zeros(size(h_time_fband{i}));
        else
            h_time_freq{i, j} = h_time_fband{i};
        end
    end
end
h_time_freq = cell2mat(h_time_freq);

ALL_h = cell(bin_chan, 1); % channel x 1
for i = 1:bin_chan
    ALL_h{i, 1} = h_time_freq; % ones(freq x time, 1)
end
ALL_h = cell2mat(ALL_h);
description = 'freq_3bands';
%}

% construct h for freq
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

ALL_h = cell(bin_chan, 1); % channel x 1
for i = 1:bin_chan
    ALL_h{i, 1} = h_time_freq; % ones(freq x time, 1)
end
ALL_h = cell2mat(ALL_h);
description = 'freq';
%}

% construct h for time
%{
% construct ones for time bins
h_time = eye(bin_time);
ALL_h = repmat(h_time, bin_freq * bin_chan, 1);
description = 'time';
%}

% construct h for contrast of one channel vs the rest channels
%{
h_theOneChan = (bin_chan - 1) * ones(nb_tf, 1);
ALL_h = cell(bin_chan);
for i = 1:bin_chan
    for j = 1:bin_chan
        if i == j
            ALL_h{i, j} = h_theOneChan;
        else
            ALL_h{i, j} = -1 * ones(size(h_theOneChan));
        end
    end
end
ALL_h = cell2mat(ALL_h);
description = 'channelContrast1vs62';
%}

% construct h for contrast of one time vs the rest comparison
%{
h_theOneTime = -1 * ones(bin_time);
for i = 1:bin_time
    for j = 1:bin_time
        if i == j
            h_theOneTime(i, j) = 39;
        end
    end
end
ALL_h = repmat(h_theOneTime, bin_chan * bin_freq, 1);
description = 'timeContrast1vs39';
%}



% construct h for channel and freq
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
%}

% NG H
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



% construct h for channel and time
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

% construct h for contrast of paired channel comparison
%{
    h_theOneChan = ones(nb_tf, bin_chan - 1);

    h_sub = cell(bin_chan);
    for i = 1:bin_chan
        for j = 1:bin_chan
            h_theOtherChan = zeros(nb_tf, bin_chan - 1);
            if i == j
                h_sub{i, j} = h_theOneChan;
            else
                if i < j
                    h_theOtherChan(:, i) = -1 * ones(nb_tf, 1);
                else
                    h_theOtherChan(:, i - 1) = -1 * ones(nb_tf, 1);
                end
                h_sub{i, j} = h_theOtherChan;
            end
        end
    end
    h_sub = cell2mat(h_sub);
%}

% construct h for contrast of one freq vs the rest comparison
%{
    h_theOneFreq = (bin_freq - 1) * ones(bin_time, 1);
    h_theOtherFreq = -1 * ones(bin_time, 1);
    
    h_freq_time = cell(bin_freq);
    for i = 1:bin_freq
        for j = 1:bin_freq
            if i == j
                h_freq_time{i, j} = h_theOneFreq;
            else
                h_freq_time{i, j} = h_theOtherFreq;
            end
        end
    end
    h_freq_time = cell2mat(h_freq_time);
    h_sub = repmat(h_freq_time, bin_chan, 1);
%}

%
save(fullfile(data_dir, ['ALL_h_', description, addition_description]), 'ALL_h', '-v7.3');