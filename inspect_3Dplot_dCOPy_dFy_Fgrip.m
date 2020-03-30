close all; clearvars; clc

%% load aligned data
summarypath = uigetdir;
summaryfile = fullfile(summarypath, 'behavior_resultSummary_w_dist2ideal.mat');
summary = load(summaryfile);

%% alinged with lift onset
nsub = height(summary.All_info_trial);
for sub_i = 1:nsub
    clear tmp_block
    sub_id = summary.All_info_trial.subID{sub_i, 1};
    data = summary.All_info_trial.data{sub_i, 1};
    dt = summary.All_info_trial.dt(sub_i, 1); % in ms
    data_lft_onset_aligned = summary.All_info_trial.data{sub_i, 1}.onset_aligned;
    [~, lft_onset] = min(abs(data_lft_onset_aligned{1, 1}.time)); % lift onset is at 0 s
    nep = height(data);
    
    % plot 3D plot: Fgrip, dFy, dCOPy
    % get block index
    [~, ~, tmp_block_id] = unique(data.context);
    tmp_jump = diff(tmp_block_id);
    tmp_b = 1;
    b = 1;
    for ep = 1:nep - 1
        if abs(tmp_jump(ep)) > 0
            tmp_block{b, 1} = tmp_b:ep;
            tmp_b = ep + 1;
            b = b + 1;
        end
    end
    tmp_block{b, 1} = tmp_b:nep; % the last block
    
    
    figure('DefaultAxesFontSize', 18, 'units', 'normalized', 'outerposition', [0, 0, 1, 1])
    hold on
    % ideal Mcom surfaces
    [dCOPy, dFy] = meshgrid(-50:50, -10:0.2:10);
    hwidth = -29.6 * 2;
    dCOPz = hwidth; Mcom = 395;
    Fg = (dCOPz * dFy + Mcom) ./ dCOPy;
    surf(dCOPy, dFy, Fg)
    Mcom = -395;
    Fg = (dCOPz * dFy + Mcom) ./ dCOPy;
    surf(dCOPy, dFy, Fg)
    
    for block_ep = 1:length(tmp_block)
        eps = tmp_block{block_ep, 1};
        switch data.context{eps(1), 1}
            case 'IL'
                linespec = 'ro'; tcc = 'r'; linestyl = '-r';
            case 'TR'
                linespec = 'bx'; tcc = 'b'; linestyl = '-b';
            case 'PT'
                linespec = 'k^'; tcc = 'k'; linestyl = '-k';
        end
        
        tmp_data = data_lft_onset_aligned(eps, 1);

        x = cellfun(@(x) x.dCOPy(lft_onset), tmp_data);
        y = cellfun(@(x) x.dFy(lft_onset), tmp_data);
        z = cellfun(@(x) x.GripF(lft_onset), tmp_data);
        xs = data.dCOPy_ideal(eps, 1);
        ys = data.dFy_ideal(eps, 1);
        zs = data.gripF_ideal(eps, 1);
        
        % data points
        plot3(x, y, z, 'LineStyle', 'none', 'Marker', 'none')
        for i = 1:length(eps)
            if (i == 1)
                fw = 'bold';
                fs = 12;
            elseif (i == length(eps))
                fw = 'bold';
                fs = 5;
            else
                fw = 'normal';
                fs = 5; %length(eps) - i;
            end
            x_t = x(i);
            y_t = y(i);
            z_t = z(i);
            text(x_t, y_t, z_t, num2str( data.trial(eps(i), 1) ), 'Color', tcc, 'FontWeight', fw, 'FontSize', fs)
            
            plot3([x_t, xs(i)], [y_t, ys(i)], [z_t, zs(i)], linestyl, 'Marker', 'none')
        end
        % ideal points
        plot3(xs, ys, zs, linespec);
    end
    grid on
    axis equal
    xlim([-30, 50])
    ylim([-20, 20])
    zlim([5, 35])
    set(gca, 'View', [70, 35])
    
    xlabel('{\Delta}COPy_{TH-VF} (mm)');
    ylabel('{\Delta}Fy_{TH-VF} (N)');
    zlabel('Fn_{TH} (N)');
    hold off
    title(sub_id)
    
    savefig(fullfile(summarypath, 'plots', [sub_id, '_3D_dist']))
end