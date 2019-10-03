close all; clearvars; clc

%% load aligned data
pathname = uigetdir;
filelist = dir(fullfile(pathname, 'behavior_pRoll_23subs.mat'));

load(fullfile(pathname, filelist(1).name));

%% fill in NaN for less than 19 trials
nsub = length(mx_pRoll_cond);
for sub_i = 1:nsub
    for cond_i = 1:size(mx_pRoll_cond, 2)
        if length(mx_pRoll_cond{sub_i, cond_i}) < 19
            for t_i = (length(mx_pRoll_cond{sub_i, cond_i}) + 1):19
                mx_pRoll_cond{sub_i, cond_i} = [mx_pRoll_cond{sub_i, cond_i}; [NaN, NaN]];
            end
        end
    end
end

%% fit exponential for data from all subjects
% for 1:IL, 2:TR
f = cell(3, 1);
data = cell(3, 1);
for cond_i = 1:2
    tmp = cell(nsub, 1);
    flip_pRoll = 1;
    if cond_i == 2
        flip_pRoll = -1;
    end
    for i = 1:nsub
        tmp{i, 1} = [(1:length(mx_pRoll_cond{i, cond_i}))', flip_pRoll * mx_pRoll_cond{i, cond_i}(:, 2)];
    end
    tmp = cell2mat(tmp);
    data{cond_i} = tmp(~isnan(tmp(:, 2)), :);
    
    f{cond_i} = fit(data{cond_i}(:, 1), data{cond_i}(:, 2), 'a * exp(b * x)', 'StartPoint', [1, 1]);
end

% for PT
pRoll_PT = cell(nsub, 1);
for sub_i = 1:nsub
    tmp = [mx_pRoll_cond{sub_i, 3:5}];
    tmp = tmp(:, 2:2:6)';
    tmp = [repmat([1; 2; 3], length(tmp), 1), reshape(tmp, [], 1)];
    
    tmp = tmp(~isnan(tmp(:, 2)), :);    
    pRoll_PT{sub_i, 1} = tmp;
end
pRoll_PT = cell2mat(pRoll_PT);
f{3} = fit(pRoll_PT(:, 1), pRoll_PT(:, 2), 'a * exp(b * x)', 'StartPoint', [1, 1]);
data{3} = pRoll_PT;

%% plot
for i = 1:3
    figure(i)
    disp(f{i})
    plot(f{i}, data{i}(:, 1), data{i}(:, 2))
end

%% construct dummy variables for IL, TR, PT
x = 1:1:19;
pRoll_fit = nan(length(x), 3);
for i = 1:3
    tmp = feval(f{i}, x);
    pRoll_fit(:, i) = tmp / max(tmp);
end
[~, filename, ~] = fileparts(filelist(1).name);
save(fullfile(pathname, [filename, '_pRoll_fit']), 'pRoll_fit');


%{
% % % %% fit exponential for data from single subject
% % % figure
% % % hold on
% % % cond_id = 2; % 1:IL, 2:TR
% % % f1 = cell(nsub, 1);
% % % flip_pRoll = 1;
% % % if cond_id == 2
% % %     flip_pRoll = -1;
% % % end
% % % for sub_i = 1:nsub
% % %     tmp = flip_pRoll * mx_pRoll_cond{sub_i, cond_id}(:, 2);
% % %     data = [(1:length(tmp))', tmp];
% % %     data = data(~isnan(data(:, 2)), :);
% % % 
% % %     f1{sub_i, :} = fit(data(:, 1), data(:, 2), 'a * exp(b * x)', 'StartPoint', [0, 1]);
% % %     plot(f1{sub_i, :}, data(:, 1), data(:, 2))
% % % end
% % % hold off
% % % legend off
% % % %%
% % % f_PT = cell(nsub, 1);
% % % pRoll_PT = cell(nsub, 1);
% % % figure
% % % hold on
% % % for sub_i = 1:nsub
% % %     tmp = [mx_pRoll_cond{sub_i, 3:5}];
% % %     tmp = tmp(:, 2:2:6)';
% % %     tmp = [repmat([1; 2; 3], length(tmp), 1), reshape(tmp, [], 1)];
% % %     
% % %     tmp = tmp(~isnan(tmp(:, 2)), :);
% % %     f_PT{sub_i, :} = fit(tmp(:, 1), tmp(:, 2), 'a * exp(b * x)', 'StartPoint', [0, 1]);
% % %     plot(f_PT{sub_i, :}, tmp(:, 1), tmp(:, 2))
% % %     
% % %     pRoll_PT{sub_i, 1} = tmp;
% % % end
% % % hold off
% % % fit_coeff = nan(sub_i, 2);
% % % for sub_i = 1:nsub
% % %     fit_coeff(sub_i, :) = coeffvalues(f_PT{sub_i, 1});
% % % end
% % % avg_fit_coeff = mean(fit_coeff, 1);
% % % figure
% % % plot(f_PT, pRoll_PT(:, 1), pRoll_PT(:, 2))
% % % hold on
% % % x = 1:0.2:3;
% % % plot(x, avg_fit_coeff(1, 1) * exp(avg_fit_coeff(1, 2) * x))
% % % 
% % % figure
% % % plot(x, 3.164 * exp(-0.4725 * x), '-b', x, 2.2041 * exp(0.7807 * x), '-r')
%}




