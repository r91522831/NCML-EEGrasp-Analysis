close all; clearvars; clc

%% load aligned data
summarypath = uigetdir;
summaryfile = fullfile(summarypath, 'behavior_resultSummary_w_dist2ideal.mat');
summary = load(summaryfile);

%%
ntrial = 95;
nsub = height(summary.All_info_trial);
mx_x_sub = nan(ntrial, nsub);
dist_x_sub = mx_x_sub;
for sub_i = 1:nsub
    sub_id = summary.All_info_trial.subID{sub_i, 1};
    data = summary.All_info_trial.data{sub_i, 1};
    dt = summary.All_info_trial.dt(sub_i, 1); % in ms
    
    % rotate the odd subs start from L to the other side for averaging
    if mod(str2double(sub_id((end-1):end)), 2)
        mx_x_sub(data.trial, sub_i) = -data.Mcom;
    else
        mx_x_sub(data.trial, sub_i) = data.Mcom;
    end
    dist_x_sub(data.trial, sub_i) = data.dist2ideal;
end

%% compute mean and std
data = summary.All_info_trial.data{end, 1}; % find the complete experiment protocol, IL: 19, 1TR3PT: 19
Mideal = nan(height(data), 1);
for i = 1:height(data)
    switch data.context{i, 1}
        case {'IL', 'PT'}
            Mideal(i, 1) = 395;
        case 'TR'
            Mideal(i, 1) = -395;
    end
end
mcom_x_sub = mx_x_sub - Mideal;
select_sub = 1:nsub;%[1, 5, 6, 8];
avg_mx_x_sub = nanmean(mcom_x_sub(:, select_sub), 2);
std_mx_x_sub = nanstd(mcom_x_sub(:, select_sub), [], 2);
avg_dist_x_sub = nanmean(dist_x_sub(:, select_sub), 2);
std_dist_x_sub = nanstd(dist_x_sub(:, select_sub), [], 2);

%% get block index
data = summary.All_info_trial.data{end, 1}; % find the complete experiment protocol, IL: 19, 1TR3PT: 19
ntrial = height(data);
[~, ~, tmp_block_id] = unique(data.context);
tmp_jump = diff(tmp_block_id);
tmp_b = 1;
b = 1;
for ep = 1:ntrial - 1
    if abs(tmp_jump(ep)) > 0
        tmp_block{b, 1} = tmp_b:ep;
        tmp_b = ep + 1;
        b = b + 1;
    end
end
tmp_block{b, 1} = tmp_b:ntrial;

%% plot
nsub = length(select_sub);
figure('DefaultAxesFontSize', 18, 'units', 'normalized', 'outerposition', [0, 0, 1, 1])
subplot(2, 1, 1)
hold on
nsession = length(tmp_block);
for i = 1:nsession
    switch data.context{tmp_block{i, 1}(1), 1}(1) % the first trial in the block
        case 'I'
            line_spec = '-*r';
        case 'T'
            line_spec = '-ob';
        case 'P'
            line_spec = '-xk';
        otherwise
    end
    errorbar(tmp_block{i, 1}, avg_mx_x_sub(tmp_block{i, 1}, 1), std_mx_x_sub(tmp_block{i, 1}, 1) / sqrt(nsub), line_spec)
end
hold off
legend({'IL', 'TR', 'PT'}, 'Location', 'northwest')
ylabel('Mcom (N-mm)')
% xlabel('trial')
xticklabels([])
xlim([0, 96])
title(['error bars are SE across ', num2str(nsub), ' subjects'])

subplot(2, 1, 2)
hold on
nsession = length(tmp_block);
for i = 1:nsession
    switch data.context{tmp_block{i, 1}(1), 1}(1) % the first trial in the block
        case 'I'
            line_spec = '-*r';
        case 'T'
            line_spec = '-ob';
        case 'P'
            line_spec = '-xk';
        otherwise
    end
    errorbar(tmp_block{i, 1}, avg_dist_x_sub(tmp_block{i, 1}, 1), std_dist_x_sub(tmp_block{i, 1}, 1) / sqrt(nsub), line_spec)
end
hold off
ylabel('dist to ideal surface')
xlabel('trial')
xlim([0, 96])
title([summary.All_info_trial.subID{select_sub}])




%% plot for individual subs
sub_i = 4;
subID = summary.All_info_trial.subID{sub_i, 1}; 
data = summary.All_info_trial.data{sub_i, 1}; % find the complete experiment protocol, IL: 19, 1TR3PT: 19

clear tmp_mcom_err Mideal tmp_block
if mod(str2double(subID((end-1):end)), 2)
    tmp_mcom = -data.Mcom;
else
    tmp_mcom = data.Mcom;
end
for i = 1:height(data)
    switch data.context{i, 1}
        case {'IL', 'PT'}
            Mideal(i, 1) = 395;
        case 'TR'
            Mideal(i, 1) = -395;
    end
end
tmp_mcom_err = tmp_mcom - Mideal;

ntrial = height(data);

[~, ~, tmp_block_id] = unique(data.context);
tmp_jump = diff(tmp_block_id);
tmp_b = 1;
b = 1;
for ep = 1:ntrial - 1
    if abs(tmp_jump(ep)) > 0
        tmp_block{b, 1} = tmp_b:ep;
        tmp_b = ep + 1;
        b = b + 1;
    end
end
tmp_block{b, 1} = tmp_b:ntrial;

figure('DefaultAxesFontSize', 18, 'units', 'normalized', 'outerposition', [0, 0, 1, 1])
subplot(2, 1, 1)
hold on
nsession = length(tmp_block);
for i = 1:nsession
    switch data.context{tmp_block{i, 1}(1), 1}(1) % the first trial in the block
        case 'I'
            line_spec = '-*r';
        case 'T'
            line_spec = '-ob';
        case 'P'
            line_spec = '-xk';
        otherwise
    end
    plot(tmp_block{i, 1}, tmp_mcom_err(tmp_block{i, 1}, 1), line_spec)
end
hold off
legend({'IL', 'TR', 'PT'}, 'Location', 'northwest')
ylabel('Mcom (N-mm)')
% xlabel('trial')
xticklabels([])
xlim([0, 96])
title(subID)

subplot(2, 1, 2)
hold on
nsession = length(tmp_block);
for i = 1:nsession
    switch data.context{tmp_block{i, 1}(1), 1}(1) % the first trial in the block
        case 'I'
            line_spec = '-*r';
        case 'T'
            line_spec = '-ob';
        case 'P'
            line_spec = '-xk';
        otherwise
    end
    plot(tmp_block{i, 1}, data.dist2ideal(tmp_block{i, 1}, 1), line_spec)
end
hold off
ylabel('dist to ideal surface')
xlabel('trial')
xlim([0, 96])


