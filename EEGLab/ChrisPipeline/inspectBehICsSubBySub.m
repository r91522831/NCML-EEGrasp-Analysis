clear; close all; clc;

%% Get behavior data: Mcom and peak Roll
% '/Users/yen-hsunwu/Dropbox (ASU)/BIDS_format/NCML-EEGrasp/sub-09/beh/mat/S009_temp_result.mat'
subID = 'sub-09';
behpath = fullfile('/Users/yen-hsunwu/Dropbox (ASU)/BIDS_format/NCML-EEGrasp/', subID, 'beh/mat');
behfilename = ['S0', subID((end-1):end), '_temp_result.mat'];
beh = load(fullfile(behpath, behfilename));

%% Get EEG data: voltage for each channel
% load EEG csd data. Raw EEG data are in EEG.dataRaw. CSD EEG data are in EEG.data
eegpath = fullfile('/Users/yen-hsunwu/Dropbox (ASU)/BIDS_format/NCML-EEGrasp/', subID, 'eeg/set');
% % % eegfilename = [subID, '_epoched_ICA_SouceLocalized.set'];
% % % eegfilename = 'sub-09_timefreq.set';
eegfilename = [subID, '_timefreq.set'];

EEG = pop_loadset('filename', eegfilename, 'filepath', eegpath);
% % % % Get time freq data
% % % tf = EEG.icatf;

subID = EEG.filename(1:6);
% block ep into IL1, IL2~19, TR1, TR2~19, PT1 1, PT1 2~18, PT2 1, PT2 2~18, PT3 1, PT3 2~18
% epBlock = defineBlocks(EEG);
epBlock = defineBlocksSeperatePT123(EEG);
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

% compute EEG IC activation voltage average through IL1, IL2~19, TR1, TR2~19, PT1 1, PT1 2~18, PT2 1, PT2 2~18, PT3 1, PT3 2~18
nIC = size(EEG.icaact, 1);
% nIC x ntwin_eeg x nepb
actBlock = nan(nIC, ntwin_eeg, nepb);
actBlock_std = actBlock;
for i_epb = 1:nepb
    actBlock(:, :, i_epb) = mean(EEG.icaact(:, :, epBlock{1, i_epb}), 3);
    actBlock_std(:, :, i_epb) = std(EEG.icaact(:, :, epBlock{1, i_epb}), 0, 3) ./ sqrt(18);
end

actBlock(:, :, nepb + 1) = mean(EEG.icaact(:, :, :), 3);
actBlock_std(:, :, nepb + 1) = std(EEG.icaact(:, :, :), 0, 3) ./ sqrt(95);

%% tf computation
% get the time and index for the window from -3000 to 3000 ms
ind_win = dsearchn(EEG.icatf.tf_times', [-3000, 3000]');
t_win_tf = EEG.icatf.tf_times(ind_win(1):ind_win(2))';
dt_tf = mean(diff(t_win_tf));
ntwin_tf = length(t_win_tf);

% compute freq power average through IL1, IL2~19, TR1, TR2~19, PT1 1, PT1 2~18, PT2 1, PT2 2~18, PT3 1, PT3 2~18
freqz = EEG.icatf.tf_freqs;
rg_freq_band = { {'\theta', '4~8 Hz'}, {'\alpha', '8~13 Hz'}, {'\beta_{low}', '13~20 Hz'}, {'\beta_{high}', '20~30 Hz'}; ...
                 find(freqz >  4 & freqz <=  8), find(freqz >   8 & freqz <= 13), ...
                 find(freqz > 13 & freqz <= 20), find(freqz >  20 & freqz <= 30) };
nIC = size(EEG.icaact, 1);
nfb = length(rg_freq_band);

pow = cell(nfb, 1);
for i_IC = 1:nIC
    for i_fb = 1:nfb
        pow{i_fb, 1}(i_IC, :, :) = mean(EEG.icatf.tf_ersp{i_IC, 1}(rg_freq_band{2, i_fb}, :, :), 1);
    end
end

powBlock = cell(nfb, 1);
powBlock_std = powBlock;
for i_IC = 1:nIC
    for i_fb = 1:nfb
        for i_epb = 1:nepb
            powBlock{i_fb, 1}(i_IC, :, i_epb) = mean(pow{i_fb, 1}(i_IC, :, epBlock{1, i_epb}), 3);
            powBlock_std{i_fb, 1}(i_IC, :, i_epb) = std(pow{i_fb, 1}(i_IC, :, epBlock{1, i_epb}), 0, 3);
        end
    end
end

for i_fb = 1:nfb
    powBlock{i_fb, 1}(:, :, nepb + 1) = mean(pow{i_fb, 1}(:, :, :), 3);
    powBlock_std{i_fb, 1}(:, :, nepb + 1) = std(pow{i_fb, 1}(:, :, :), 0, 3) ./ sqrt(95);
end

%% plot
% % % dir_save = '/Users/yen-hsunwu/Dropbox (ASU)/NCML-EEGrasp/for Apr10 2020';

%{
if mod(str2double(subID((end-1):end)), 2) ~= 0
    side = 1;
else
    side = -1;
end
%}

%% behavior plots
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



%% 
% % % selected = [1, 4, 6, 7, 9, 12, 13, 15, 16, 21, 25, 36, 44, 48];
selected = [1, 2, 3, 4, 5, 6, 11, 12, 14, 21, 22,26, 39, 47];
% plot the selected ICs
printTrialMapsAxes(EEG, [-3000, 3000], 'ICA', [], [4, 4], selected);

%%
% blocks: 'IL_1', 'IL_{2-19}', 'TR_1', 'TR_{2-19}', 'PT1_1', 'PT1_{2-19}', 'PT2_1', 'PT2_{2-19}', 'PT3_1', 'PT3_{2-19}', 'all'
% nIC x time x blocks (10 + 1)
fqtext = {'\theta', '\alpha', '\beta_{low}', '\beta_{high}'};
fqfilename = {'theta', 'alpha', 'lowbeta', 'highbeta'};
figfolder = fullfile('/Users/yen-hsunwu/Dropbox (ASU)/NCML-EEGrasp/for May25 2020', EEG.filename(1:6));
for i_IC = selected
    close all;
    data_plot = EEG.icaact(i_IC, :, :);
    actBlock_plot = actBlock(i_IC, :, :);
    actBlock_std_plot = actBlock_std(i_IC, :, :);
    fig_act = contextICinspectPlot(EEG, data_plot, epBlock, t_win_eeg, actBlock_plot, actBlock_std_plot, i_IC, 'icaact');
    figfilename = [EEG.filename(1:6), '_IC', num2str(i_IC), '_IC_activity'];
    savefig(fig_act, fullfile(figfolder, figfilename))
    
    for i_fq = 1:4
        data_plot = pow{i_fq, 1}(i_IC, :, :); % unit is power in '10 * log10( \muV^{2}/Hz )'
        blockavg = powBlock{i_fq, 1}(i_IC, :, :);
        blockstd = powBlock_std{i_fq, 1}(i_IC, :, :);
        fig_tf_alpha = contextICinspectPlot(EEG, data_plot, epBlock, t_win_tf, blockavg, blockstd, i_IC, fqtext{i_fq});
        figfilename = [EEG.filename(1:6), '_IC', num2str(i_IC), '_IC_', fqfilename{i_fq}];
        savefig(fig_tf_alpha, fullfile(figfolder, figfilename))
    end
    
end



%% compute correlation between IC pairs
%{
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
%}

%{
subplot(2, 1, 2)
for idx = selected
    hold on
    plot(0.001 * t_win_eeg, normalize(actBlock(idx, :, ep_all)))
    hold off
end
%}
