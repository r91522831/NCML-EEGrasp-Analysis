clear; close all; clc;

disp('Select the voltage folder in the EEG directory:')
All_data = uigetdir;
All_data_list = dir(fullfile(All_data, 'sub-*'));
All_root_dir = fileparts(fileparts(All_data));
All_chanlocs = load(fullfile(All_data, 'chanlocs.mat'));
All_n_elect = length(All_chanlocs.chanlocs);

disp([num2cell((1:length(All_data_list))'), {All_data_list.name}']);
All_selected_sub = input('Which subject(s) to plot? ');
if isempty(All_selected_sub)
    All_selected_sub = 1:length(All_data_list);
end
All_n_sub = length(All_selected_sub);

%%
All_f = figure('units', 'normalized', 'outerposition', [0, 0, 1, 1]);
All_f_reg = figure('units', 'normalized', 'outerposition', [0, 0, 1, 1]);
All_plot_i = 1;

All_peak_ratio = cell(All_n_sub, 3);
for All_i = All_selected_sub
    clearvars -except All_*;
    sub_id = All_data_list(All_i).name(1:6);
    
    % load peak roll
    filepath = fullfile(fileparts(All_data_list(All_i).folder), [sub_id, '_onset']);
    filename = [sub_id, '_combine_behavior.set'];
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
    EEG = pop_loadset('filename', filename, 'filepath', filepath);
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
    
    
    % load voltage
    vol = load(fullfile(All_data_list(All_i).folder, All_data_list(All_i).name));
    sig_t2l = nan(All_n_elect, 1);
    sig_l2w = nan(All_n_elect, 1);
    for e_i = 1:All_n_elect
        sig_t2l(e_i, 1) = robustcov(vol.vol_t2l(:, e_i));
        sig_l2w(e_i, 1) = robustcov(vol.vol_l2w(:, e_i));
    end
    
    All_peak_ratio{All_plot_i, 1} = sub_id;
    [All_peak_ratio{All_plot_i, 2}, peak_ratio_e]= max(sig_l2w ./ sig_t2l);
    All_peak_ratio{All_plot_i, 3} = All_chanlocs.chanlocs(peak_ratio_e).labels;
    
    
    clim = max([sig_l2w; sig_t2l]);
    figure(All_f)
    subplot(All_n_sub, 3, 3 * (All_plot_i - 1) + 1, 'FontSize', 20)
    topoplot(sig_t2l, All_chanlocs.chanlocs, 'maplimits', [-clim, clim]);
    if All_plot_i == 1, title('contact to lift onset'); end
    text(-2, 0, sub_id, 'FontSize', 20)
    subplot(All_n_sub, 3, 3 * (All_plot_i - 1) + 2, 'FontSize', 20)
    topoplot(sig_l2w, All_chanlocs.chanlocs, 'maplimits', [-clim, clim]);
    if All_plot_i == 1, title('200 ms window after lift onset'); end
    colorbar('FontSize', 20);
    
    subplot(All_n_sub, 3, 3 * (All_plot_i - 1) + 3, 'FontSize', 20)
    topoplot( (sig_l2w - sig_t2l), All_chanlocs.chanlocs, 'maplimits', [-20, 20]);%0.2 * [-clim, clim]);
    if All_plot_i == 1, title('difference'); end
    colorbar('FontSize', 20);
    
    
    
    figure(All_f_reg)
    hold on
    x = robustcov(abs(EEG.behavior.obj_roll_peak.peakRoll));
    [~, y] = robustcov(sig_l2w - sig_t2l);
    plot(x, y, 'o')
    text(x, y, sub_id, 'FontSize', 18)
    hold off
    xlabel('variance of peak roll')
    ylabel('mean variance difference')
    set(gca, 'FontSize', 18)
    
    
    All_plot_i = All_plot_i + 1;
end