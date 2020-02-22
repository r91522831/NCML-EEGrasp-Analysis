disp('Select the voltage folder in the EEG directory:')
All_data = uigetdir;
All_data_list = dir(fullfile(All_data, 'sub-*'));
All_root_dir = fileparts(fileparts(All_data));
All_chanlocs = load(fullfile(All_data, 'chanlocs.mat'));
All_n_elect = length(All_chanlocs.chanlocs);

disp([num2cell((1:length(All_data_list))'), {All_data_list.name}']);
All_selected_sub = input('Which subject(s) to preprocess EEG? ');
if isempty(All_selected_sub)
    All_selected_sub = 1:length(All_data_list);
end
All_n_sub = length(All_selected_sub);

%%
All_f = figure('units', 'normalized', 'outerposition', [0, 0, 1, 1]);
All_plot_i = 1;
for All_i = All_selected_sub
    clearvars -except All_*;
    sub_id = All_data_list(All_i).name(1:6);
    vol = load(fullfile(All_data_list(All_i).folder, All_data_list(All_i).name));
    sig_t2l = nan(All_n_elect, 1);
    sig_l2w = nan(All_n_elect, 1);
    for e_i = 1:All_n_elect
        sig_t2l(e_i, 1) = robustcov(vol.vol_t2l(:, e_i));
        sig_l2w(e_i, 1) = robustcov(vol.vol_l2w(:, e_i));
    end
    
    clim = max([sig_l2w; sig_t2l]);
    figure(All_f)
    subplot(All_n_sub, 3, 3 * (All_plot_i - 1) + 1, 'FontSize', 20)
    topoplot(sig_t2l, All_chanlocs.chanlocs, 'maplimits', [-clim, clim]);
    if All_plot_i == 1, title('contact to lift onset'); end
    text(-1.5, 0, sub_id, 'FontSize', 20)
    subplot(All_n_sub, 3, 3 * (All_plot_i - 1) + 2, 'FontSize', 20)
    topoplot(sig_l2w, All_chanlocs.chanlocs, 'maplimits', [-clim, clim]);
    if All_plot_i == 1, title('200 ms window after lift onset'); end
    colorbar('FontSize', 20);
    
    subplot(All_n_sub, 3, 3 * (All_plot_i - 1) + 3, 'FontSize', 20)
    topoplot( (sig_l2w - sig_t2l), All_chanlocs.chanlocs, 'maplimits', [-20, 20]);%0.2 * [-clim, clim]);
    if All_plot_i == 1, title('difference'); end
    colorbar('FontSize', 20);
    
    All_plot_i = All_plot_i + 1;
end