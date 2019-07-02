clear; close all; clc;

All_dirpath = uigetdir();
All_dirlist = dir(fullfile(All_dirpath, 'sub*coeff.mat'));

All_fig_path = fullfile(All_dirpath, 'linear_figs');
if ~exist(All_fig_path, 'dir')
    mkdir(All_fig_path);
end
disp([num2cell((1:length(All_dirlist))'), {All_dirlist.name}']);
All_selected_sub = input('Which subject(s) to plot erpimage? (default: all)');
if isempty(All_selected_sub)
    All_selected_sub = 1:length(All_dirlist);
end

% ERSP = beta0 * IL + beta1 * TR + beta2 * PT + beta3 * IL * Roll + beta4 * TR * Roll + beta5 * PT * Roll
All_stat_plot = {'coeff est', 'coeff pValue', 'Rsquared'};

All_nb_electrode = 5; All_nb_beta = 6;
All_Model_coeff_est = cell(All_nb_electrode, All_nb_beta);
for All_sub_i = All_selected_sub
    clearvars -except All_*; close all;
    sub_id = All_dirlist(All_sub_i).name(1:6);
    load(fullfile(All_dirlist(All_sub_i).folder, All_dirlist(All_sub_i).name));
    load(fullfile(All_dirlist(All_sub_i).folder, 'misc.mat'));
    
    if All_sub_i == All_selected_sub(end) 
        time_stamp = tf_times{1, 1}(:, :, 1)';
        time_tick = cellstr(num2str(time_stamp / 1000, '%.2f'));
        freq_tick = cellstr(num2str(tf_freqs{1, 1}(:, :, 1)', '%.0f'));
    end
    
    for electrode_i = 1:All_nb_electrode
        for beta_i = 1:All_nb_beta
            All_Model_coeff_est{electrode_i, beta_i}(:, :, All_sub_i) = Model_coeff_est{electrode_i, 1}(:, :, beta_i);
        end
    end
end

% plot robust mean
h = gobjects(All_nb_electrode, 1);
for i = 1:All_nb_electrode
    h(i) = figure;
end

trim_percentage = 20; % 5%
for electrode_i = 1:All_nb_electrode
    figure(h(electrode_i))

    for beta_i = 1:All_nb_beta
        ax = subplot(3, 2, beta_i);
        imagesc(trimmean(All_Model_coeff_est{electrode_i, beta_i}, trim_percentage, 3));
        colorbar;
%         colorbar('east');
        set(gca, 'YDir', 'normal');
        
%         if beta_i == nb_beta
        xticklabels(time_tick(xticks));
        xlabel('time(s)');
%         else
%             ax.XTickLabels = [];
%         end
        onset = find(diff(sign(time_stamp))) + 0.5;
        vline(onset, '--r');
        
        yticklabels(freq_tick(yticks));
        ylabel('freq(Hz)')
%         ylabel(['\beta_', num2str(beta_i - 1), ' freq(Hz)'])
        title(['\beta_', num2str(beta_i - 1)])
    end

    suptitle([num2str(trim_percentage), ' % trimmeaned', ' ', electrodes_name{electrode_i}, ' ', 'coeff est', ' ', ' (ERSP=\beta_0\cdotIL+\beta_1\cdotTR+\beta_2\cdotPT+\beta_3\cdotIL\cdotErr+\beta_4\cdotTR\cdotErr+\beta_5\cdotPT\cdotErr+\epsilon)'])
    
    filename_fig = ['All_sub', '_', num2str(trim_percentage), '_percent_trimmed_', electrodes_name{electrode_i}, '_', 'coeff est', '.fig'];
    if ~exist(fullfile(All_fig_path, 'coeff est'), 'dir')
        mkdir(fullfile(All_fig_path, 'coeff est'));
    end
    savefig(h(electrode_i), fullfile(All_fig_path, 'coeff est', filename_fig),'compact');
end