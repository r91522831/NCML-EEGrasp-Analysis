
timelim = [-500, 500]; % -500 to 500 ms
slid_win = 50; % 50 ms
overlap = 25;
% generate sliding window of 50 ms overlap 25 ms from -500 to 500 ms
swins = [timelim(1):overlap:(timelim(2) - slid_win); (timelim(1) + slid_win):overlap:timelim(2)]';

rr_max = cell(60, 60);
rr_all = cell(60, 60);
for c1 = 1:60
    for c2 = 1:60
        if c2 <= c1, continue; end % comparison order doesn't matter.
        srr_max = nan(length(swins), length(swins));
        for t1 = 1:length(swins)
            for t2 = 1:length(swins)
                nnn = dsearchn(EEG.times', swins(t1, :)');
                x1 = squeeze(squeeze( mean(EEG.icaact(c1, nnn(1):nnn(2), :), 2) ));
                nnn = dsearchn(EEG.times', swins(t2, :)');
                x2 = squeeze(squeeze( mean(EEG.icaact(c2, nnn(1):nnn(2), :), 2) ));
                [rr, ~] = corr(x1, x2, 'Type', 'Spearman');
                
                srr_max(t1, t2) = rr;
            end
        end
        [tmp_rrmax, tmp_rrmax_lind]= max( abs(srr_max(:)) );
        [tmp_rrmax_row, tmp_rrmax_col] = ind2sub(size(srr_max), tmp_rrmax_lind);
        
        rr_max{c1, c2} = [srr_max(tmp_rrmax_ind), swins(tmp_rrmax_row, :), swins(tmp_rrmax_col, :)];
        rr_all{c1, c2} = srr_max;
    end
end

%%
big_rr = [];
i = 1;
for c1 = 1:60
    for c2 = 1:60
        if isempty(rr_max{c1, c2}), continue; end
        
        if rr_max{c1, c2}(1) > 0.1
            big_rr{i, 1} = [c1, c2];
            i = i + 1;
        end
    end
end

%%
%{
% recompute icaact using ICs extra from CSDed set to original set
EEG.data = EEG.dataRaw;
EEG.icaact = [];
EEG = eeg_checkset( EEG, 'ica' );
%}
%% Sort icaact into context blocks
epBlock = defineBlocks(EEG);
nepb = length(epBlock);

EEG.icaact = cat( 3, EEG.icaact(:, :, [epBlock{1, 1:2}]), nan(size(EEG.icaact, 1), size(EEG.icaact, 2), 1), ...
                     EEG.icaact(:, :, [epBlock{1, 3:4}]), nan(size(EEG.icaact, 1), size(EEG.icaact, 2), 1), ...
                     EEG.icaact(:, :, [epBlock{1, 5:6}]) );

%%
printTrialMapsAxes(EEG, [-200, 500], 'ICA', [], [1, 1], 21);
%%
printTrialMapsAxes(EEG, [-3000, 3000], 'ICA', [], [7, 6], 1:39);
printTrialMapsAxes(EEG, [-3000, 3000], 'ICA', [], [6, 5], 1:29);
%%
printTrialMapsAxes(EEG, [-3000, 3000], 'ICA', [], [9, 7], 1:61);
disp('finished!')
%%
% sub-09 selected ICs:  1    2   3   4  5    6  11  12  14  21  22      26      39      47
% sub-21 selected ICs: 41	 3	 9	11	5	 2	14	 7	 4	41	24       8      11      37
% sub-11 selected ICs:  7	17	 7	 9	1	25	 6	18	 6	 7	> 50 mm	> 50 mm	 5      24
% sub-14 selected ICs: 26	42	 2	 8	1	12	21	14	43	28	> 50 mm	 8      > 50 mm	12
% sub-02 selected ICs: 28	40	43	25	4	21	23	24	23	41	> 50 mm	23      31      31
% % % selected = [1, 2, 3, 4, 5, 6, 11, 12, 14, 21, 22, 26, 39, 47]; % sub-09
% % % selected = [41,	 3,	 9,	11,	5,	 2,	14,	 7,	 4,	41,	 24,	  8,     11,    37]; % sub-21
% % % selected = [ 7,	17,	 7,	 9,	1,	25,	 6,	18,	 6,	 7, nan,    nan,      5,    24]; % sub-11
% % % selected = [26,	42,	 2,	 8,	1,	12, 21,	14,	43,	28, nan,	  8,    nan,	12]; % sub-14
selected = [28,	40,	43,	25,	4,	21,	23, 24,	23, 41, nan,     23,     31,    31]; % sub-02

printTrialMapsAxes(EEG, [-3000, 3000], 'ICA', [], [4, 4], selected);

%% Sorted the tf_ersp into context blocks
for idx = 1:length(EEG.icatf.tf_ersp)
    tmp = EEG.icatf.tf_ersp{idx, 1};
    sorted = cat( 3, tmp(:, :, [epBlock{1, 1:2}]), nan(size(tmp, 1), size(tmp, 2), 1), ...
                     tmp(:, :, [epBlock{1, 3:4}]), nan(size(tmp, 1), size(tmp, 2), 1), ...
                     tmp(:, :, [epBlock{1, 5:6}]) );
    EEG.icatf.tf_ersp{idx, 1} = sorted;
end

%%
% theta band: 4 to 8 Hz; alpha: 9 ~ 13 Hz; low beta: 14 ~ 20 Hz; high beta: 21 ~ 30 Hz
% % % selected = [1, 4, 6, 7, 9, 12, 13, 15, 16, 21, 25, 36, 44, 48]; % this ICs are for sub-09 version 00
% % % selected = [1, 2, 3, 4, 5,  6, 11, 12, 14, 21, 22, 26, 39, 47]; % this ICs are for sub-09 version 01
selected = [1, 2, 3, 4, 5, 6, 11, 12, 14, 21, 22,26, 39, 47];
% layout = [10, 6];
layout = [4, 4];
printTrialMapsAxes(EEG, [-3000, 3000], 'ICA', [4, 8], layout, selected, [-600, -200]);
suptitle('\theta 4-8 Hz')
printTrialMapsAxes(EEG, [-3000, 3000], 'ICA', [9, 13], layout, selected, [-600, 200]);
suptitle('\alpha 8-13 Hz')
printTrialMapsAxes(EEG, [-3000, 3000], 'ICA', [14, 20], layout, selected, [-600, -200]);
suptitle('low\beta 13-20')
printTrialMapsAxes(EEG, [-3000, 3000], 'ICA', [21, 30], layout, selected, [-600, -200]);
suptitle('high\beta 20-30 Hz')




%%
% plot r > 0.1
disp(['Plotting ', num2str(length(big_rr)), ' pair(s) of ICs'])
for i = 1:20%length(big_rr)
    c1 = big_rr{i}(1);
    c2 = big_rr{i}(2);
    compareICA2(EEG, c1, c2, rr_max{c1, c2}(2:3), rr_max{c1, c2}(4:5))
end

%%
nselected = length(selected);
big_rr = [];
i = 1;
for c1 = 1:nselected
    for c2 = 1:nselected
        if isnan(r{c1, c2}), continue; end
        
        if r{c1, c2} > 0.3
            big_rr{i, 1} = [c1, c2];
            i = i + 1;
        end
    end
end
%%
disp(['Plotting ', num2str(length(big_rr)), ' pair(s) of ICs'])
for i = 1%:length(big_rr)
    c1 = selected(big_rr{i}(1));
    c2 = selected(big_rr{i}(2));
    compareICA2(EEG, c1, c2)
end
