clear; close all; clc;

%% Get behavior data: Mcom and peak Roll
beh = load('/Users/yen-hsunwu/Dropbox (ASU)/BIDS_format/NCML-EEGrasp/sub-09/beh/mat/S009_temp_result.mat');

%% Get EEG data: voltage for each channel
% load EEG csd data. Raw EEG data are in EEG.dataRaw. CSD EEG data are in EEG.data
eegfilename = 'sub-09_timefreq.set';
EEG = pop_loadset('filename', eegfilename, 'filepath', '/Users/yen-hsunwu/Dropbox (ASU)/BIDS_format/NCML-EEGrasp/sub-09/eeg/set/');
% Get time freq data
tf = EEG.icatf;
subID = EEG.filename(1:6);
% block ep into IL1, IL2~19, TR1, TR2~19, PT1, or PT2~57
epBlock = defineBlocks(EEG);
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
%{
% define the number of topoplot slices
nslice = 15; % 15 slices: each slice is about 250 ms
ind_slice = round(linspace(1, n_t_win, nslice + 1)); % find the intervals
ind_win_slice = cell(1, nslice);
t_win_slice = cell(1, nslice);
for i = 1:nslice
    fst_ind = ind_slice(i);
    if i ~= nslice, end_ind = ind_slice(i + 1) - 1; 
    else, end_ind = ind_slice(i + 1);  
    end
    ind_win_slice{1, i} = ind_win_eeg(fst_ind:end_ind, 1);
    t_win_slice{1, i} = t_win_eeg(fst_ind:end_ind, 1);
end
%}

% compute EEG IC activation voltage average through IL1, IL2~19, TR1, TR2~19, PT1, or PT2~57
nIC = size(EEG.icaact, 1);
% nIC x ntwin_eeg x nepb
actBlock = nan(nIC, ntwin_eeg, nepb);
actBlock_std = actBlock;
for i_epb = 1:nepb
    actBlock(:, :, i_epb) = mean(EEG.icaact(:, :, epBlock{1, i_epb}), 3);
    actBlock_std(:, :, i_epb) = std(EEG.icaact(:, :, epBlock{1, i_epb}), 0, 3);
end

%% tf computation
% get the time and index for the window from -3000 to 3000 ms
ind_win = dsearchn(EEG.icatf.tf_times', [-3000, 3000]');
t_win_tf = EEG.icatf.tf_times(ind_win(1):ind_win(2))';
dt_tf = mean(diff(t_win_tf));
ntwin_tf = length(t_win_tf);

%{
% define the number of topoplot slices
nslice_tf = 15; % 15 slices: each slice is about 250 ms
ind_slice_tf = round(linspace(1, n_t_win_tf, nslice_tf + 1)); % find the intervals
ind_win_slice_tf = cell(1, nslice_tf);
t_win_slice_tf = cell(1, nslice_tf);
for i = 1:nslice_tf
    fst_ind = ind_slice_tf(i);
    if i ~= nslice_tf, end_ind = ind_slice_tf(i + 1) - 1; 
    else, end_ind = ind_slice_tf(i + 1);  
    end
    ind_win_slice_tf{1, i} = ind_win_tf(fst_ind:end_ind, 1);
    t_win_slice_tf{1, i} = t_win_tf(fst_ind:end_ind, 1);
end
%}

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

%% plot
% % % dir_save = '/Users/yen-hsunwu/Dropbox (ASU)/NCML-EEGrasp/for Apr10 2020';

%{
if mod(str2double(subID((end-1):end)), 2) ~= 0
    side = 1;
else
    side = -1;
end
%}

linspec = {'-r', '--r', '-b', '--b', '-k', ':k'};
fig = figure('DefaultAxesFontSize', 18, 'units', 'normalized', 'outerposition', [0, 0, 1, 1]);
%{
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
i_IC = 26;
% hold on
for i_epb = 1:nepb
    subplot(6, 1, i_epb)
    h(i_epb) = shadedErrorBar(t_win_eeg, actBlock(i_IC, :, i_epb), actBlock_std(i_IC, :, i_epb), linspec{i_epb}, 1);
%     axis off
end
% hold off
% vline(0, '--k')


%%







for i_contrast = 1:4
    vol_plot = vol{:, i_contrast};
    pow_plot = pow{:, i_contrast};
    
    fig = figure('DefaultAxesFontSize', 18, 'units', 'normalized', 'outerposition', [0, 0, 1, 1]);
    max_vol = max(abs(vol_plot(:)));
    for i_win = 1:nslice
        subplot(2 + nfb, nslice + 1, i_win)
        topoplot(vol_plot(:, i_win), EEG.chanlocs, 'electrodes', 'off', 'maplimits', [-max_vol, max_vol]);
        
        if i_win == 1, text(-1, 0, vol_text{i_contrast}, 'HorizontalAlignment', 'center', 'Rotation', 90, 'FontWeight', 'bold'); end
        if i_win ~= nslice
            title(round(t_win_slice{1, i_win}(round(0.5 * length(t_win_slice{1, i_win})))))
        else
            title([num2str(round(t_win_slice{1, i_win}(round(0.5 * length(t_win_slice{1, i_win}))))), ' ms'])
        end
    end
    subplot(2 + nfb, nslice + 1, i_win + 1)
    caxis([-max_vol, max_vol])
    set(gca, 'visible','off');
    colorbar;
    
    for i_fb = 1:nfb
        tmp_pow = pow_plot(:, :, i_fb);
        max_fb = max(abs(tmp_pow(:)));
        for i_win = 1:nslice
            subplot(2 + nfb, nslice + 1, i_fb * (nslice + 1) + i_win)
            topoplot(pow_plot(:, i_win, i_fb), EEG.chanlocs, 'electrodes', 'off', 'maplimits', [-max_fb, max_fb]);
            
            if i_win == 1, text(-1, 0, {rg_freq_band{1, i_fb}{1, 1}, pow_text{i_contrast}}, 'HorizontalAlignment', 'center', 'Rotation', 90, 'FontWeight', 'bold'); end
            if i_fb == 1
                if i_win ~= nslice
                    title(round(t_win_slice_tf{1, i_win}(round(0.5 * length(t_win_slice_tf{1, i_win})))))
                else
                    title([num2str(round(t_win_slice_tf{1, i_win}(round(0.5 * length(t_win_slice_tf{1, i_win}))))), ' ms'])
                end
            end
        end
        subplot(2 + nfb, nslice + 1, i_fb * (nslice + 1) + i_win + 1)
        caxis([-max_fb, max_fb])
        set(gca, 'visible','off');
        colorbar;
    end
    
    subplot(2 + nfb, nslice + 1, [(1 + nfb) * (nslice + 1) + 1, (2 + nfb) * (nslice + 1) - 1])
    yyaxis left
    hold on
    h1 = shadedErrorBar(0.001 * t_win_beh, -nanmean(mCom(:, 3:19), 2) * side, nanstd(mCom(:, 3:19), [], 2), '-.b', 1);
    switch i_contrast
        case 1
            h2 = plot(0.001 * t_win_beh, mCom(:, epBlock{i_contrast}) * side, '-.r');
        case 2
            h2 = shadedErrorBar(0.001 * t_win_beh, nanmean(mCom(:, epBlock{i_contrast}), 2) * side, nanstd(mCom(:, epBlock{i_contrast}), [], 2), '-.r', 1);
        case 3
            h2 = plot(0.001 * t_win_beh, -mCom(:, epBlock{i_contrast}) * side, '-.r');
        case 4
            h2 = shadedErrorBar(0.001 * t_win_beh, -nanmean(mCom(:, epBlock{i_contrast}), 2) * side, nanstd(mCom(:, epBlock{i_contrast}), [], 2), '-.r', 1);
    end
    ylim([-300, 500])
    ylabel('Mcom (Nmm)')
    hold off
    yyaxis right
    hold on
    h3 = shadedErrorBar(0.001 * t_win_beh, -nanmean(pRoll(:, 3:19), 2) * side, nanstd(pRoll(:, 3:19), [], 2), '-b', 1);
    switch i_contrast
        case 1
            h4 = plot(0.001 * t_win_beh, -pRoll(:, epBlock{i_contrast}) * side, '-r');
        case 2
            h4 = shadedErrorBar(0.001 * t_win_beh, nanmean(pRoll(:, epBlock{i_contrast}), 2) * side, nanstd(pRoll(:, epBlock{i_contrast}), [], 2), '-r', 1);
        case 3
            h4 = plot(0.001 * t_win_beh, pRoll(:, epBlock{i_contrast}) * side, '-r');
        case 4
            h4 = shadedErrorBar(0.001 * t_win_beh, -nanmean(pRoll(:, epBlock{i_contrast}), 2) * side, nanstd(pRoll(:, epBlock{i_contrast}), [], 2), '-r', 1);
    end
    ylim([-10, 50])
    vline(0, ':k', 'lft')
    ylabel('Roll angle ({\circ})')
    hold off
    xlabel([subID, ' time (s)'])
    if mod(i_contrast, 2) ~= 0
        legend([h1.mainLine, h2, h3.mainLine, h4], {'Mcom_{IL_{3~19}}', beh_text{i_contrast}{1}, 'Roll_{IL_{3~19}}', beh_text{i_contrast}{2}}, 'Location', 'east')
    else
        legend([h1.mainLine, h2.mainLine, h3.mainLine, h4.mainLine], {'Mcom_{IL_{3~19}}', beh_text{i_contrast}{1}, 'Roll_{IL_{3~19}}', beh_text{i_contrast}{2}}, 'Location', 'east')
    end
    
    filename = [subID, fig_text{i_contrast}];
    saveas(fig, fullfile(dir_save, filename), 'fig')
    saveas(fig, fullfile(dir_save, filename), 'png')
end