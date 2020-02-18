%%
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
    figdir = '/Users/yenhsunw/Dropbox (ASU)/NCML-EEGrasp/for Feb18 2020';
    figname = ['erpimage_', electrode_label, '_', sub_id];
    savefig(h, fullfile(figdir, figname))
    saveas(h, fullfile(figdir, [figname, '.png']))
    close(h)
end