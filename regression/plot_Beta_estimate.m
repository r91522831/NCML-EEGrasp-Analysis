clear; close all; clc;

All_dirpath = uigetdir();
All_dirlist = dir(fullfile(All_dirpath, 'sub*coeff.mat'));

All_fig_path = fullfile(All_dirpath, 'linear_figs');
if ~exist(All_fig_path, 'dir')
    mkdir(All_fig_path);
end
disp([num2cell((1:length(All_dirlist))'), {All_dirlist.name}']);
selected_sub = input('Which subject(s) to plot erpimage? (default: all)');
if isempty(selected_sub)
    selected_sub = 1:length(All_dirlist);
end

% ERSP = beta0 * IL + beta1 * TR + beta2 * PT + beta3 * IL * Roll + beta4 * TR * Roll + beta5 * PT * Roll
All_stat_plot = {'coeff est', 'coeff pValue', 'Rsquared'};

for All_sub_i = selected_sub
    clearvars -except All_*; close all;
    sub_id = All_dirlist(All_sub_i).name(1:6);
    load(fullfile(All_dirlist(All_sub_i).folder, All_dirlist(All_sub_i).name));
    load(fullfile(All_dirlist(All_sub_i).folder, 'misc.mat'));
    
    time_stamp = tf_times{1, 1}(:, :, 1)';
    time_tick = cellstr(num2str(time_stamp / 1000, '%.2f'));
    freq_tick = cellstr(num2str(tf_freqs{1, 1}(:, :, 1)', '%.0f'));
    
    nb_electrode = 5; nb_beta = 6;
    stat_image = {Model_coeff_est, Model_coeff_p, Model_Rsquared};
    for k = 1:3
        h = gobjects(nb_electrode, 1);
        for i = 1:nb_electrode
            h(i) = figure;
        end

        for electrode_i = 1:nb_electrode
            figure(h(electrode_i))
            if k ~= 3
                for beta_i = 1:nb_beta
                    ax = subplot(3, 2, beta_i);
                    imagesc(stat_image{k}{electrode_i, 1}(:, :, beta_i));
                    colorbar;
%                     colorbar('east');
                    set(gca, 'YDir', 'normal');
                    
%                     if beta_i == nb_beta
                        xticklabels(time_tick(xticks));
                        xlabel('time(s)');
%                     else
%                         ax.XTickLabels = [];
%                     end
                    onset = find(diff(sign(time_stamp))) + 0.5;
                    vline(onset, '--r');
                    
                    yticklabels(freq_tick(yticks));
                    ylabel('freq(Hz)')
%                     ylabel(['\beta_', num2str(beta_i - 1), ' freq(Hz)'])
                    title(['\beta_', num2str(beta_i - 1)])
                end
            else
                imagesc(stat_image{k}{electrode_i, 1}(:, :));
                colorbar;
                set(gca, 'YDir', 'normal');
                    
                xticklabels(time_tick(xticks));
                xlabel('time(s)');
                onset = find(diff(sign(time_stamp))) + 0.5;
                vline(onset, '--r');
                
                yticklabels(freq_tick(yticks));
                ylabel('freq(Hz)')
            end
            suptitle([sub_id, ' ', electrodes_name{electrode_i}, ' ', All_stat_plot{k}, ' ', ' (ERSP=\beta_0\cdotIL+\beta_1\cdotTR+\beta_2\cdotPT+\beta_3\cdotIL\cdotErr+\beta_4\cdotTR\cdotErr+\beta_5\cdotPT\cdotErr+\epsilon)'])

            filename_fig = [sub_id, '_', electrodes_name{electrode_i}, '_', All_stat_plot{k}, '.fig'];
            if ~exist(fullfile(All_fig_path, All_stat_plot{k}), 'dir')
                mkdir(fullfile(All_fig_path, All_stat_plot{k}));
            end
            savefig(h(electrode_i), fullfile(All_fig_path, All_stat_plot{k}, filename_fig),'compact');
        end
    end
end