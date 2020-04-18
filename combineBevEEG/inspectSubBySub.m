clear; close all; clc;

%% Get behavior data: Mcom and peak Roll
% % % beh = load('/Users/yen-hsunwu/Dropbox (ASU)/BIDS_format/NCML-EEGrasp/sub-09/beh/mat/S009_temp_result.mat');
beh = load('/Users/yen-hsunwu/Dropbox (ASU)/BIDS_format/NCML-EEGrasp/sub-17/beh/mat/S017_temp_result.mat');
%% Get time freq data
% % % tf = load('/Users/yen-hsunwu/Dropbox (ASU)/NCML-EEGrasp/EEG/eeglab/004 RemoveEye90/linear/RAW/sub-09_timefreq.mat');
tf = load('/Users/yen-hsunwu/Dropbox (ASU)/NCML-EEGrasp/EEG/eeglab/004 RemoveEye90/linear/RAW/sub-17_timefreq.mat');
%% Get EEG data: voltage for each channel
% load EEG csd data. Raw EEG data are in EEG.dataRaw. CSD EEG data are in EEG.data
% % % EEG = pop_loadset('filename', 'sub-09_eeg_csd.set', 'filepath', '/Users/yen-hsunwu/Dropbox (ASU)/NCML-EEGrasp/EEG/eeglab/004 RemoveEye90/sub-09_onset/');
EEG = pop_loadset('filename', 'sub-17_eeg_csd.set', 'filepath', '/Users/yen-hsunwu/Dropbox (ASU)/NCML-EEGrasp/EEG/eeglab/004 RemoveEye90/sub-17_onset/');
% correct the channel locations to standard 10-5 cap
EEG = pop_chanedit(EEG, 'lookup', '/Users/yen-hsunwu/Dropbox (Personal)/Programming/Matlab/myLibrary/eeglab2019_1/plugins/dipfit/standard_BESA/standard-10-5-cap385.elp');
subID = EEG.filename(1:6);

%% behavior computation
nep = length(beh.resultantF);
dt_beh = diff(beh.info_time_trigger{1, 1}(1:2, 1)); % sampling rate for behavior data
t_win_beh = -1000:dt_beh:3000; % -1000 to 3000 ms

mCom = nan(length(t_win_beh), 20);
pRoll = mCom;
for i_ep = 1:nep
    ind_win_beh = beh.info_onset_time.lft_ind(i_ep, 1) + round(t_win_beh ./ dt_beh); 
    mCom(:, i_ep) = beh.resultantF{i_ep, 1}.mx(ind_win_beh, 1);
    pRoll(:, i_ep) = beh.angTilt2R{i_ep, 1}(ind_win_beh, 1);
end

%% EEG computation
% get the time and index for the window from -1000 to 3000 ms
[~, ind_min_t] = min(abs(EEG.times - (-1000))); % index closet to -1 s
[~, ind_max_t] = min(abs(EEG.times - 3000)); % index closet to 3 s
ind_win_eeg = (ind_min_t:ind_max_t)';
n_t_win = length(ind_win_eeg);
t_win_eeg = EEG.times(ind_win_eeg)';
dt_eeg = mean(diff(t_win_eeg));

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

% define trial blocks into IL, TR and PT
nep = length(EEG.epoch);
context = cell(nep, 1);
for i = 1:nep
    context(i, 1) = EEG.epoch(i).eventcondType(1);
end
[~, ~, tmp_block_id] = unique(context);
tmp_jump = diff(tmp_block_id);
tmp_b = 1;
b = 1;
for i_ep = 1:(nep - 1)
    if abs(tmp_jump(i_ep)) > 0
        block{b, 1} = tmp_b:i_ep;
        tmp_b = i_ep + 1;
        b = b + 1;
    end
end
block{b, 1} = tmp_b:nep;
% block ep into TR1, TR3~TR19, PT1, or PT2~PT57
ep_block = {block{2, 1}, [block{6:2:end, 1}], block{3, 1}(:, 1), [block{3, 1}(:, 3), block{5:2:end, 1}]};

% compute EEG voltage average through IL3~IL19, TR1, TR3~TR19, PT1, or PT2~PT57
nelec = EEG.nbchan;
vol = table(nan(nelec, nslice), nan(nelec, nslice), nan(nelec, nslice), nan(nelec, nslice), ...
            'VariableNames', {'TR1minusIL3to19', 'TR3to19minusIL3to19', 'PT1minusIL3to19', 'PT3to57minusIL3to19'});
for i_win = 1:nslice
    tmp_IL13to19 = mean(mean(EEG.dataRaw(:, ind_win_slice{1, i_win}, block{1, 1}(:, 3:end)), 3), 2);
    tmp_TR1 = mean(mean(EEG.dataRaw(:, ind_win_slice{1, i_win}, ep_block{1}), 3), 2);
    tmp_TR3to19 = mean(mean(EEG.dataRaw(:, ind_win_slice{1, i_win}, ep_block{2}), 3), 2);
    tmp_PT1 = mean(mean(EEG.dataRaw(:, ind_win_slice{1, i_win}, ep_block{3}), 3), 2);
    tmp_PT3to57 = mean(mean(EEG.dataRaw(:, ind_win_slice{1, i_win}, ep_block{4}), 3), 2);
    vol.TR1minusIL3to19(:, i_win) = tmp_TR1 - tmp_IL13to19;
    vol.TR3to19minusIL3to19(:, i_win) = tmp_TR3to19 - tmp_IL13to19;
    vol.PT1minusIL3to19(:, i_win) = tmp_PT1 - tmp_IL13to19;
    vol.PT3to57minusIL3to19(:, i_win) = tmp_PT3to57 - tmp_IL13to19;
end

%% tf computation
% get the time and index for the window from -1000 to 3000 ms
[~, ind_min_t_tf] = min(abs(tf.tf_times - (-1000))); % index closet to -1 s
[~, ind_max_t_tf] = min(abs(tf.tf_times - 3000)); % index closet to 3 s
ind_win_tf = (ind_min_t_tf:ind_max_t_tf)';
n_t_win_tf = length(ind_win_tf);
t_win_tf = tf.tf_times(ind_win_tf)';
dt_tf = mean(diff(t_win_tf));

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

% compute freq power average through IL3~IL19
freqz = tf.tf_freqs;
rg_freq_band = { {'\theta', '4~8 Hz'}, {'\alpha', '8~13 Hz'}, {'\beta_{low}', '13~20 Hz'}, {'\beta_{high}', '20~30 Hz'}; ...
                 find(freqz >  4 & freqz <=  8), find(freqz >   8 & freqz <= 13), ...
                 find(freqz > 13 & freqz <= 20), find(freqz >  20 & freqz <= 30) };
nelec = EEG.nbchan;
nfb = length(rg_freq_band);

pow_TR1minusIL3to19 = nan(nelec, nslice_tf, nfb);
pow = table(nan(nelec, nslice_tf, nfb), nan(nelec, nslice_tf, nfb), nan(nelec, nslice_tf, nfb), nan(nelec, nslice_tf, nfb), ...
            'VariableNames', {'TR1minusIL3to19', 'TR3to19minusIL3to19', 'PT1minusIL3to19', 'PT3to57minusIL3to19'});
for i_win = 1:nslice_tf
    for i_elec = 1:nelec
        for i_fb = 1:nfb
            tmp_IL13to19 = tf.tf_ersp.Var1{i_elec, 1}(rg_freq_band{2, i_fb}, ind_win_slice_tf{1, i_win}, block{1, 1}(:, 3:end));
            tmp_IL13to19 = abs(mean(tmp_IL13to19(:)));
            tmp_TR1 = tf.tf_ersp.Var1{i_elec, 1}(rg_freq_band{2, i_fb}, ind_win_slice_tf{1, i_win}, ep_block{1});
            tmp_TR1 = abs(mean(tmp_TR1(:)));
            tmp_TR3to19 = tf.tf_ersp.Var1{i_elec, 1}(rg_freq_band{2, i_fb}, ind_win_slice_tf{1, i_win}, ep_block{2});
            tmp_TR3to19 = abs(mean(tmp_TR3to19(:)));
            tmp_PT1 = tf.tf_ersp.Var1{i_elec, 1}(rg_freq_band{2, i_fb}, ind_win_slice_tf{1, i_win}, ep_block{3});
            tmp_PT1 = abs(mean(tmp_PT1(:)));
            tmp_PT3to57 = tf.tf_ersp.Var1{i_elec, 1}(rg_freq_band{2, i_fb}, ind_win_slice_tf{1, i_win}, ep_block{4});
            tmp_PT3to57 = abs(mean(tmp_PT3to57(:)));
            
            pow.TR1minusIL3to19(i_elec, i_win, i_fb) = tmp_TR1 - tmp_IL13to19;
            pow.TR3to19minusIL3to19(i_elec, i_win, i_fb) = tmp_TR3to19 - tmp_IL13to19;
            pow.PT1minusIL3to19(i_elec, i_win, i_fb) = tmp_PT1 - tmp_IL13to19;
            pow.PT3to57minusIL3to19(i_elec, i_win, i_fb) = tmp_PT3to57 - tmp_IL13to19;
        end
    end
end

%% plot
vol_text = {'V_{TR_1-IL_{3~19}}', 'V_{TR_{3~19}-IL_{3~19}}', 'V_{PT_1-IL_{3~19}}', 'V_{PT_{3~57}-IL_{3~19}}'};
pow_text = {'mag_{TR_1-IL_{3~19}}', 'mag_{TR_{3~19}-IL_{3~19}}', 'mag_{PT_1-IL_{3~19}}', 'mag_{PT_{3~57}-IL_{3~19}}'};
beh_text = {{'Mcom_{TR_1}', 'Roll_{TR_1}'}, {'Mcom_{TR_{3~19}}', 'Roll_{TR_{3~19}}'}, {'Mcom_{PT_1}', 'Roll_{PT_1}'}, {'Mcom_{PT_{3~57}}', 'Roll_{PT_{3~57}}'}};
fig_text = {'_TR1vsIL3-19', '_TR3-19vsIL3-19', '_PT1vsIL3-19', '_PT3-57vsIL3-19'};
dir_save = '/Users/yen-hsunwu/Dropbox (ASU)/NCML-EEGrasp/for Apr10 2020';

if mod(str2double(subID((end-1):end)), 2) ~= 0
    side = 1;
else
    side = -1;
end

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
            h2 = plot(0.001 * t_win_beh, mCom(:, ep_block{i_contrast}) * side, '-.r');
        case 2
            h2 = shadedErrorBar(0.001 * t_win_beh, nanmean(mCom(:, ep_block{i_contrast}), 2) * side, nanstd(mCom(:, ep_block{i_contrast}), [], 2), '-.r', 1);
        case 3
            h2 = plot(0.001 * t_win_beh, -mCom(:, ep_block{i_contrast}) * side, '-.r');
        case 4
            h2 = shadedErrorBar(0.001 * t_win_beh, -nanmean(mCom(:, ep_block{i_contrast}), 2) * side, nanstd(mCom(:, ep_block{i_contrast}), [], 2), '-.r', 1);
    end
    ylim([-300, 500])
    ylabel('Mcom (Nmm)')
    hold off
    yyaxis right
    hold on
    h3 = shadedErrorBar(0.001 * t_win_beh, -nanmean(pRoll(:, 3:19), 2) * side, nanstd(pRoll(:, 3:19), [], 2), '-b', 1);
    switch i_contrast
        case 1
            h4 = plot(0.001 * t_win_beh, -pRoll(:, ep_block{i_contrast}) * side, '-r');
        case 2
            h4 = shadedErrorBar(0.001 * t_win_beh, nanmean(pRoll(:, ep_block{i_contrast}), 2) * side, nanstd(pRoll(:, ep_block{i_contrast}), [], 2), '-r', 1);
        case 3
            h4 = plot(0.001 * t_win_beh, pRoll(:, ep_block{i_contrast}) * side, '-r');
        case 4
            h4 = shadedErrorBar(0.001 * t_win_beh, -nanmean(pRoll(:, ep_block{i_contrast}), 2) * side, nanstd(pRoll(:, ep_block{i_contrast}), [], 2), '-r', 1);
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