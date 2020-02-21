close all; clearvars; clc

%% load aligned data
All_path = uigetdir;
All_filelist = dir(fullfile(All_path, '*_temp_result.mat'));

disp([num2cell((1:length(All_filelist))'), {All_filelist.name}']);
All_selected_sub = input('Which subject(s) to process? ');

if isempty(All_selected_sub)
    All_selected_sub = 1:length(All_filelist);
end

%%
for All_i = All_selected_sub
    clearvars -except All_*; close all;
    sub_id = All_filelist(All_i).name(1:4);
    disp(['Start processing ', sub_id, ' ...']);
    
    load(fullfile(All_path, All_filelist(All_i).name));
    
    f = figure('units', 'normalized', 'outerposition', [0, 0, 1, 1]);
    nb_trailing_zeros = 5;
    for i = 1:height(ind_lft_onset) + 1
        if i > height(ind_lft_onset), subplot(10, 10, i, 'Visible', 'off'); break; end
        subplot(10, 10, i)
        hold on
        yyaxis left
        plot(0.001 * info_time_trigger{i, 1}(1:end-nb_trailing_zeros) - info_onset_time(i, 1), finger_Th{i, 1}.fz(1:end-nb_trailing_zeros))
        plot(0.001 * info_time_trigger{i, 1}(1:end-nb_trailing_zeros) - info_onset_time(i, 1), finger_V{i, 1}.fz(1:end-nb_trailing_zeros))
        ylim([-50, 50])
        yyaxis right
        plot(0.001 * info_time_trigger{i, 1}(1:end-nb_trailing_zeros) - info_onset_time(i, 1), angTilt2R{i, 1}(1:end-nb_trailing_zeros))
        ylim([-20, 20])
        xlim([-3, 3])
        vline(0, '-r');
        vline(info_onset_time(i, 2) - info_onset_time(i, 1), '-k')
        title(['trial ', num2str(i), ' ', file_list(i).name(11:12)])
    end
    
    mtit( f, file_list(1).name(1:4), 'fontsize', 18, 'xoff', 0, 'yoff', 0.025 )
    
    
    figname = [sub_id, '_grip_force_roll'];
    savefig(f, fullfile(All_path, figname))
    saveas(f, fullfile(All_path, [figname, '.png']))
    close(f)
end