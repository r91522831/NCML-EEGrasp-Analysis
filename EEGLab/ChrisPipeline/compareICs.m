close all; clear; clc;
%%
projpath = '/Users/yen-hsunwu/Dropbox (ASU)/BIDS_format/NCML-EEGrasp/';
basesubid = 'sub-09';
basefilepath = fullfile(projpath, basesubid, 'eeg', 'set');
basefilename = [basesubid, '_timefreq.set'];
baseEEG = pop_loadset('filename', basefilename, 'filepath', basefilepath);
% 'sub-02', 'sub-11', 'sub-14', 'sub-21'
compsubid = 'sub-02';
compfilepath = fullfile(projpath, compsubid, 'eeg', 'set');
compfilename = [compsubid, '_timefreq.set'];
compEEG = pop_loadset('filename', compfilename, 'filepath', compfilepath);

figfolder = fullfile('/Users/yen-hsunwu/Dropbox (ASU)/NCML-EEGrasp/for Jun01 2020', compEEG.filename(1:6));
if ~isfolder(figfolder), mkdir(figfolder); end

%%
baseidx = [1, 2, 3, 4, 5, 6, 11, 12, 14, 21, 22,26, 39, 47]; % this ICs are for sub-09
dist = nan(compEEG.nbchan, length(baseidx));
selected = nan(3, length(baseidx));
selected_dist = selected;
for b_idx = 1:length(baseidx)
    for idx = 1:compEEG.nbchan
        if ~isempty(compEEG.dipfit.model(idx).posxyz)
            tmp_dist = sqrt(sum((baseEEG.dipfit.model(baseidx(b_idx)).posxyz - compEEG.dipfit.model(idx).posxyz).^2));
            dist(idx, b_idx) = tmp_dist;
        end
    end
    tmp_dist = dist(:, b_idx);
    [selected_dist(:, b_idx), selected(:, b_idx)] = mink(tmp_dist, 3);
end

%%
for b_idx = 1:length(baseidx)
    fig = figure('DefaultAxesFontSize', 18, 'units', 'normalized', 'outerposition', [0, 0, 1, 1]);
    row = 3; col = 4;
    subplot(row, col, [1, 5, 9]);
    basecomp = baseidx(b_idx);
% % %     topoplot(baseEEG.dipfit.model(basecomp).datapot, baseEEG.chanlocs);
    topoplot(baseEEG.icawinv(:, basecomp), baseEEG.chanlocs);
    [~, kkkkk] = max( baseEEG.etc.ic_classification.ICLabel.classifications(basecomp, :) );
    tt = { [ num2str(basecomp), ' ', baseEEG.etc.ic_classification.ICLabel.classes{kkkkk}, ' ', num2str(100 * baseEEG.etc.ic_classification.ICLabel.classifications(basecomp, kkkkk), '%2.1f') ], ...
        baseEEG.dipfit.model(basecomp).areadk, ...
        baseEEG.setname(1:6) };
    title(tt, 'Units', 'normalized', 'Position', [0.5, -0.1, 0]);
    plotICidx = selected(:, b_idx);
    plotICdist = selected_dist(:, b_idx);
    for idx = 1:length(plotICidx)
        comp = plotICidx(idx);
        subplot(row, col, (idx - 1) * col + 2)
% % %         topoplot(compEEG.dipfit.model(comp).datapot, compEEG.chanlocs);
        topoplot(compEEG.icawinv(:, comp), compEEG.chanlocs);
        [~, kkkkk] = max( compEEG.etc.ic_classification.ICLabel.classifications(comp, :) );
        tt = { [ num2str(comp), ' ', compEEG.etc.ic_classification.ICLabel.classes{kkkkk}, ' ', num2str(100 * compEEG.etc.ic_classification.ICLabel.classifications(comp, kkkkk), '%2.1f') ], ...
            compEEG.dipfit.model(comp).areadk, ...
            compEEG.setname(1:6) };
        title(tt, 'Units', 'normalized', 'Position', [0.5, -0.1, 0]);
        subplot(row, col, (idx - 1) * col + 3);
        topoplot(compEEG.dipfit.model(comp).sourcepot, compEEG.chanlocs);
        tt = {'source', ['dist ', num2str(plotICdist(idx), '%2.1f'), ' mm']};
        title(tt, 'Units', 'normalized', 'Position', [0.5, -0.1, 0]);
        subplot(row, col, (idx - 1) * col + 4);
        topoplot(compEEG.dipfit.model(comp).diffmap, compEEG.chanlocs);
        title(['rv ', num2str(compEEG.dipfit.model(comp).rv, '%2.2f')], 'Units', 'normalized', 'Position', [0.5, -0.1, 0]);
    end
    figfilename = [compEEG.filename(1:6), '_ICs2', baseEEG.filename(1:6), '_IC', num2str(basecomp)];
    savefig(fig, fullfile(figfolder, figfilename));
    close(fig);
end