clear; close all; clc;

%% Get behavior data: Mcom and peak Roll
beh = load('/Users/yen-hsunwu/Dropbox (ASU)/BIDS_format/NCML-EEGrasp/sub-09/beh/mat/S009_temp_result.mat');

%% Get EEG data: voltage for each channel
% load EEG csd data. Raw EEG data are in EEG.dataRaw. CSD EEG data are in EEG.data
% % % eegfilename = 'sub-09_timefreq.set';
eegfilename = 'sub-09_epoched_ICA_SouceLocalized.set';
EEG = pop_loadset('filename', eegfilename, 'filepath', '/Users/yen-hsunwu/Dropbox (ASU)/BIDS_format/NCML-EEGrasp/sub-09/eeg/set/');
% % % % Get time freq data
% % % tf = EEG.icatf;

subID = EEG.filename(1:6);
% block ep into IL1, IL2~19, TR1, TR2~19, PT1, or PT2~57
% epBlock = defineBlocks(EEG);
tmp_epb = defineBlocksSeperatePT123(EEG);
epBlock = tmp_epb(1:6);
nepb = length(epBlock);

%% behavior computation
nep = length(beh.resultantF);
dt_beh = diff(beh.info_time_trigger{1, 1}(1:2, 1)); % sampling rate for behavior data
t_win_beh = -3000:dt_beh:3000; % -3000 to 3000 ms
ntwin = length(t_win_beh);
mCom = nan(ntwin, nep);
pRoll = mCom;
for i_ep = 1:nep
    ind_win_beh = beh.info_onset_time.lft_ind(i_ep, 1) + round(t_win_beh ./ dt_beh); 
    mCom(:, i_ep) = beh.resultantF{i_ep, 1}.mx(ind_win_beh, 1);
    pRoll(:, i_ep) = beh.angTilt2R{i_ep, 1}(ind_win_beh, 1);
end

% compute beh average through IL1, IL2~19, TR1, TR2~19, PT1, or PT2~57
mComBlock = nan(ntwin, nepb);
mComBlock_std = mComBlock;
pRollBlock = mComBlock;
pRollBlock_std = mComBlock;
for i_epb = 1:nepb
    mComBlock(:, i_epb) = mean(mCom(:, epBlock{1, i_epb}), 2);
    mComBlock_std(:, i_epb) = std(mCom(:, epBlock{1, i_epb}), 0, 2);
    pRollBlock(:, i_epb) = mean(pRoll(:, epBlock{1, i_epb}), 2);
    pRollBlock_std(:, i_epb) = std(pRoll(:, epBlock{1, i_epb}), 0, 2);
end

%% EEG computation
% get the time and index for the window from -3000 to 3000 ms
ind_win = dsearchn(EEG.times', [-3000, 3000]');
t_win_eeg = EEG.times(ind_win(1):ind_win(2))';
dt_eeg = mean(diff(t_win_eeg));
ntwin_eeg = length(t_win_eeg);

% compute EEG IC activation voltage average through IL1, IL2~19, TR1, TR2~19, PT1, or PT1 2~19
nIC = size(EEG.icaact, 1);
% nIC x ntwin_eeg x nepb
actBlock = nan(nIC, ntwin_eeg, nepb);
actBlock_std = actBlock;
for i_epb = 1:nepb
    actBlock(:, :, i_epb) = mean(EEG.icaact(:, :, epBlock{1, i_epb}), 3);
    actBlock_std(:, :, i_epb) = std(EEG.icaact(:, :, epBlock{1, i_epb}), 0, 3) ./ sqrt(18);
end

actBlock(:, :, nepb + 1) = mean(EEG.icaact(:, :, :), 3);
actBlock_std(:, :, nepb + 1) = std(EEG.icaact(:, :, epBlock{1, i_epb}), 0, 3) ./ sqrt(95);
%{
%% tf computation
% get the time and index for the window from -3000 to 3000 ms
ind_win = dsearchn(EEG.icatf.tf_times', [-3000, 3000]');
t_win_tf = EEG.icatf.tf_times(ind_win(1):ind_win(2))';
dt_tf = mean(diff(t_win_tf));
ntwin_tf = length(t_win_tf);

% compute freq power average through IL1, IL2~19, TR1, TR2~19, PT1, or PT2~57
freqz = EEG.icatf.tf_freqs;
rg_freq_band = { {'\theta', '4~8 Hz'}, {'\alpha', '8~13 Hz'}, {'\beta_{low}', '13~20 Hz'}, {'\beta_{high}', '20~30 Hz'}; ...
                 find(freqz >  4 & freqz <=  8), find(freqz >   8 & freqz <= 13), ...
                 find(freqz > 13 & freqz <= 20), find(freqz >  20 & freqz <= 30) };
nIC = size(EEG.icaact, 1);
nfb = length(rg_freq_band);

pow = cell(nIC, 1);
for i_IC = 1:nIC
    for i_fb = 1:nfb
        pow{i_IC, 1}(i_fb, :, :) = mean(EEG.icatf.tf_ersp{i_IC, 1}(rg_freq_band{2, i_fb}, :, :), 1);
    end
end

powBlock = cell(nIC, 1);
powBlock_std = powBlock;
for i_IC = 1:nIC
    for i_fb = 1:nfb
        for i_epb = 1:nepb
            powBlock{i_IC, 1}(i_fb, :, i_epb) = mean(pow{i_IC, 1}(i_fb, :, epBlock{1, i_epb}), 3);
            powBlock_std{i_IC, 1}(i_fb, :, i_epb) = std(pow{i_IC, 1}(i_fb, :, epBlock{1, i_epb}), 0, 3);
        end
    end
end
%}

%% plot
% % % dir_save = '/Users/yen-hsunwu/Dropbox (ASU)/NCML-EEGrasp/for Apr10 2020';

%{
if mod(str2double(subID((end-1):end)), 2) ~= 0
    side = 1;
else
    side = -1;
end
%}

%{
linspec = {'-r', '--r', '-b', '--b', '-k', ':k'};
fig = figure('DefaultAxesFontSize', 18, 'units', 'normalized', 'outerposition', [0, 0, 1, 1]);

subplot(3, 2, 1) % mCom
hold on
for i_epb = 1:nepb
    h(i_epb) = shadedErrorBar(t_win_beh, mComBlock(:, i_epb), mComBlock_std(:, i_epb), linspec{i_epb}, 1);
end
hold off
vline(0, '--k')
legend([h(:).mainLine], {'IL1', 'IL2~19', 'TR1', 'TR2~19', 'PT1', 'PT2~57'}, 'Location', 'northwest')

subplot(3, 2, 3) % pRoll
hold on
for i_epb = 1:nepb
    h(i_epb) = shadedErrorBar(t_win_beh, pRollBlock(:, i_epb), pRollBlock_std(:, i_epb), linspec{i_epb}, 1);
end
hold off
vline(0, '--k')

subplot(3, 2, 5) % activation
%}

% % % selected = [1, 4, 6, 7, 9, 12, 13, 15, 16, 21, 25, 36, 44, 48];
selected = [1, 2, 3, 4, 5, 6, 11, 12, 14, 21, 22,26, 39, 47];
% % % 'LineStyle', 'none', 'Color' 
linspec = {'-', 'none', '-', 'none', '-', 'none', 'none'};
lincolor = {'r', 'r', 'b', 'b', 'k', 'k', 'k'};
context = {'IL1', 'IL2-19', 'TR1', 'TR2-19', 'PT1', 'PT1 2-19'};
for idx = selected
    fig = figure('DefaultAxesFontSize', 18, 'units', 'normalized', 'outerposition', [0, 0, 1, 1]);
    % hold on
    dRange = nanmax(abs(actBlock(idx, :, :)), [], 'all');
    subplot(3, 2, 1)
    hold on
    for i_epb = 1:2:6
        h(i_epb) = shadedErrorBar(0.001 * t_win_eeg, actBlock(idx, :, i_epb), actBlock_std(idx, :, i_epb), {'LineStyle', linspec{i_epb}, 'Color', lincolor{i_epb}}, 1);
    end
    hold off
    ylim([-dRange, dRange]);
    vline(0, '--k');
    ylabel('\muV');
    legend([h(1:2:6).mainLine], {'IL1', 'TR1', 'PT1'}, 'Location', 'best')

    subplot(3, 2, 3)
    hold on
    for i_epb = 2:2:6
        h(i_epb) = shadedErrorBar(0.001 * t_win_eeg, actBlock(idx, :, i_epb), actBlock_std(idx, :, i_epb), {'LineStyle', linspec{i_epb}, 'Color', lincolor{i_epb}}, 1);
    end
    hold off
% % %     ylim([-dRange, dRange]);
    vline(0, '--k');
    ylabel('\muV');
    xlabel('time (s)')
    tmp_shade = [h(2:2:6).edge];
    legend(tmp_shade(1:2:end), {'IL 2~19', 'TR 2~19', 'PT1 2~19'}, 'Location', 'best')
    
    subplot(3, 2, 5)
    ep_all = nepb + 1;
    h(ep_all) = shadedErrorBar(0.001 * t_win_eeg, actBlock(idx, :, ep_all), actBlock_std(idx, :, ep_all), {'LineStyle', linspec{ep_all}, 'Color', lincolor{ep_all}}, 1);
% % %     ylim([-dRange, dRange]);
    vline(0, '--k');
    ylabel({'all eps', '\muV'});
    xlabel('time (s)')
    
    
    subplot(2, 2, [2, 4])
    topoplot(EEG.icawinv(:, idx), EEG.chanlocs(EEG.icachansind));
    [~, kkkkk] = max( EEG.etc.ic_classification.ICLabel.classifications(idx, :) );
    tt = {[ num2str(idx), ' ', EEG.etc.ic_classification.ICLabel.classes{kkkkk}, ' ', num2str(100 * EEG.etc.ic_classification.ICLabel.classifications(idx, kkkkk), '%2.1f') ], ...
        EEG.dipfit.model(idx).areadk };
    title(tt, 'Units', 'normalized', 'Position', [0.5, -0.1, 0])
end

%% compute correlation between IC pairs
nselected = length(selected);
r = nan(nselected, nselected);
ep_all = nepb + 1;

for idx = 1:nselected
    for idy = 1:nselected
        if idy > idx
            tmp_r = corrcoef(actBlock(selected(idx), :, ep_all), actBlock(selected(idy), :, ep_all));
            r(idx, idy) = tmp_r(1, 2);
        end
    end
end
r = array2table(r, 'VariableNames', cellstr(num2str(selected')), 'RowNames', cellstr(num2str(selected')));

%%
fig = figure('DefaultAxesFontSize', 18, 'units', 'normalized', 'outerposition', [0, 0, 1, 1]);
for idx = 1:nselected
    for idy = 1:nselected
        if idy > idx
            subplot(nselected, nselected, (idx - 1) * nselected + idy)
            plot(0.001 * t_win_eeg, normalize(actBlock(selected(idx), :, ep_all)), 0.001 * t_win_eeg, normalize(actBlock(selected(idy), :, ep_all)))
            axis off
            
            subplot(nselected, nselected, (idy - 1) * nselected + idx)
            plot(0.001 * t_win_eeg, actBlock(selected(idx), :, ep_all), 0.001 * t_win_eeg, actBlock(selected(idy), :, ep_all))
            axis off
        end
    end
end
%{
subplot(2, 1, 2)
for idx = selected
    hold on
    plot(0.001 * t_win_eeg, normalize(actBlock(idx, :, ep_all)))
    hold off
end
%}
