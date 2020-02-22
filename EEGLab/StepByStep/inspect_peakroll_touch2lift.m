disp('Select the voltage folder in the EEG directory:')
All_data = uigetdir;
All_data_list = dir(fullfile(All_data, 'sub-*'));
All_root_dir = fileparts(fileparts(All_data));
All_chanlocs = load(fullfile(All_data, 'chanlocs.mat'));

disp([num2cell((1:length(All_data_list))'), {All_data_list.name}']);
All_selected_sub = input('Which subject(s) to preprocess EEG? ');
if isempty(All_selected_sub)
    All_selected_sub = 1:length(All_data_list);
end
All_n_sub = length(All_selected_sub);

All_f = figure;
All_plot_i = 1;
for All_i = All_selected_sub
    clearvars -except All_*;
    vol = load(fullfile(All_data_list(All_i).folder, All_data_list(All_i).name));
    
    clim = max([std(vol.vol_l2w), std(vol.vol_t2l)]);
    figure(All_f)
    subplot(All_n_sub, 3, 3 * (All_plot_i - 1) + 1)
    topoplot(std(vol.vol_t2l), All_chanlocs.chanlocs, 'maplimits', [-clim, clim])
    subplot(All_n_sub, 3, 3 * (All_plot_i - 1) + 2)
    topoplot(std(vol.vol_l2w), All_chanlocs.chanlocs, 'maplimits', [-clim, clim])
    cbar
    subplot(All_n_sub, 3, 3 * (All_plot_i - 1) + 3)
    topoplot(std(vol.vol_l2w)-std(vol.vol_t2l), All_chanlocs.chanlocs, 'maplimits', [-2, 2])
    cbar
    All_plot_i = All_plot_i + 1;
end