corr_summary = load('/Users/yen-hsunwu/Dropbox (ASU)/NCML-EEGrasp/EEG/eeglab/004 RemoveEye90/linear/RAW/corr_traj_7sub.mat');

%%
nsub = length(corr_summary.All_power_dist);
subID = corr_summary.All_power_dist(:, 1);
elec = corr_summary.All_power_dist{1, 2}.Properties.RowNames;
nelec = length(elec);
fband = {'\theta', '\alpha', '\beta_{low}', '\beta_{high}'}; % corr_summary.All_power_dist{1, 2}.Properties.VariableNames;

figure('DefaultAxesFontSize', 18, 'units', 'normalized', 'outerposition', [0, 0, 1, 1])
for i_elec = 41%1:nelec
    for i_fband = 1:4
        subplot(4, 1, 5 - i_fband)
        hold on
        for i_sub = 1:nsub
            time = 0.001 * corr_summary.All_power_dist{i_sub, 3}; % in s
            corr = cell2mat(corr_summary.All_power_dist{i_sub, 2}{elec{i_elec}, i_fband}{:}(:, 2));
            plot(time, corr)
        end
        if i_fband == 4
            legend(subID)
            title(elec{i_elec})
        elseif i_fband == 1
            xlabel('time (s)')
        end
        ylabel({fband{i_fband}, 'r'})
        vline(0, '--r', 'lft')
        hold off
    end
end
