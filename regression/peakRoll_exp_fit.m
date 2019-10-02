close all; clearvars; clc

%% load aligned data
pathname = uigetdir;
filelist = dir(fullfile(pathname, 'behavior_pRoll_23subs.mat'));

load(fullfile(pathname, filelist(1).name));

% fit exponential for data from single subject
data = -mx_pRoll_cond{4, 2}(:, 2);
x = (1:length(data))';
f1 = fit(x, data, 'a * exp(b * x)', 'StartPoint', [1, 1]);
plot(f1, x, data)

% fit exponential for data from all subjects
cond_id = 2; % TR
tmp = cell(length(mx_pRoll_cond), 1);
for i = 1:length(mx_pRoll_cond)
    tmp{i, 1} = [(1:length(mx_pRoll_cond{i, cond_id}))', mx_pRoll_cond{i, cond_id}(:, 2)];
end
tmp = cell2mat(tmp);
tmp = tmp(~isnan(tmp(:, 2)), :);

f = fit(tmp(:, 1), tmp(:, 2), 'a * exp(b * x)', 'StartPoint', [1, 1]);
plot(f, tmp(:, 1), tmp(:, 2))