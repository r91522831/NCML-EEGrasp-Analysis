%% plot voltage
sub_id = EEG.filename(1:6);
% electrode_id = find(strcmpi({EEG.chanlocs.labels}, 'CP4'));
cond = {'IL', 'TR', 'PT1'};
smooth = 2;

for electrode_i = 1:length({EEG.chanlocs.labels})
    % get color limit
    data = cell(length(cond), 1);
    for i = 1:length(cond)
        subplot(3, 1, i)
        data{i} = pop_erpimage(EEG, 1, [electrode_i], [[]], [EEG.chanlocs(electrode_i).labels, ' ', cond{i}], smooth, 1, {[cond{i}, '_onset']}, [], 'condID', 'yerplabel', '\muV', 'limits', [-2000 3000 NaN NaN NaN NaN NaN NaN], 'cbar', 'on', 'NoShow', 'on');
        get(gca, 'CLim')
    end
    data = cell2mat(data);
    clim = max(max(max(data), abs(min(min(data)))));
    colorlimit = [-clim, clim];
    % plot
    h = figure('units', 'normalized', 'outerposition', [0, 0, 1, 1]);
    for i = 1:length(cond)
        subplot(3, 1, i)
        pop_erpimage(EEG, 1, [electrode_i], [[]], [EEG.chanlocs(electrode_i).labels, ' ', cond{i}], smooth, 1, {[cond{i}, '_onset']}, [], 'condID', 'yerplabel', '\muV', 'limits', [-2000 3000 NaN NaN NaN NaN NaN NaN], 'cbar', 'on', 'caxis', colorlimit);
    end
    
    electrode_label = EEG.chanlocs(electrode_i).labels;
    figdir = '/Users/yenhsunw/Dropbox (ASU)/NCML-EEGrasp/for Feb18 2020/';
    figname = ['erpimage_', electrode_label, '_', sub_id];
    savefig(h, fullfile(figdir, figname))
    saveas(h, fullfile(figdir, [figname, '.png']))
    close(h)
end

%% plot power amplitude
close all
sub_id = EEG.filename(1:6);
% electrode_id = find(strcmpi({EEG.chanlocs.labels}, 'CP4'));
cond = {'IL', 'TR', 'PT1'};
smooth = 2;

% try to find [-250, -50] befroe touch onset
baselinedb = nan; % [-250, -50]; % baseline window for power amplitude plots, time to lift onset in ms

% theta band: 4 to 8 Hz; low beta: 13 ~ 19 Hz; high beta: 20 ~ 30 Hz
freqband = {[4, 8, 0.01], [13, 20, 0.01], [20, 30, 0.01]};
fb_name = {'theta', 'low_beta', 'high_beta'};
fb_title = {'\theta', 'low\beta', 'high\beta'};
for fb_i = 1:numel(freqband)
    for electrode_i = 1:length({EEG.chanlocs.labels})
        % get color limit
        data = cell(length(cond), 1);
        for i = 1:length(cond)
            [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, data{i}, ~, ~, ~, ~] = ...
                pop_erpimage( EEG, 1, [electrode_i], [[]], [EEG.chanlocs(electrode_i).labels, ' ', cond{i}, ' ', fb_title{fb_i}], ...
                smooth, 1, {[cond{i}, '_onset']}, [], ...
                'condID', 'yerplabel', '\muV', 'limits', [-2000 3000 NaN NaN NaN NaN NaN NaN], 'cbar', 'on', ...
                'plotamps', 'on', 'coher', freqband{fb_i}, 'baselinedb', baselinedb, 'NoShow', 'on' );
            get(gca, 'CLim')
        end
        data = cell2mat(data);
        clim = max(max(max(data), abs(min(min(data)))));
        colorlimit = [-clim, clim];
        % plot
        f_final = figure('units', 'normalized', 'outerposition', [0, 0, 1, 1]);
        for i = 1:length(cond)
            f_tmp = figure;
            pop_erpimage( EEG, 1, [electrode_i], [[]], [EEG.chanlocs(electrode_i).labels, ' ', cond{i}, ' ', fb_title{fb_i}], ...
                smooth, 1, {[cond{i}, '_onset']}, [], ...
                'condID', 'yerplabel', '\muV', 'limits', [-2000 3000 NaN NaN NaN NaN NaN NaN], 'cbar', 'on', 'caxis', colorlimit, ...
                'plotamps', 'on', 'coher', freqband{fb_i}, 'baselinedb', baselinedb );
            h_tmp = findobj(f_tmp, 'type', 'axes');
            
            tmp_h = copyobj(h_tmp(5), f_final); % the 5th ax is the power erpimage in the pop_erpimage plot
            colormap(tmp_h, 'jet')
            xticklabels(tmp_h, 'auto')
            xlabel(tmp_h, 'time (ms)')
            subplot(3, 1, i, tmp_h)
            colorbar(tmp_h)
            close(f_tmp)
        end
        
        electrode_label = EEG.chanlocs(electrode_i).labels;
        figdir = ['/Users/yenhsunw/Dropbox (ASU)/NCML-EEGrasp/for Feb18 2020/', fb_name{fb_i}];
        figname = ['erpimage_', fb_name{fb_i}, '_', electrode_label, '_', sub_id];
        savefig(f_final, fullfile(figdir, figname))
        saveas(f_final, fullfile(figdir, [figname, '.png']))
        close(f_final)
    end
end




