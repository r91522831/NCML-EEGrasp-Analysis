
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
printTrialMapsAxes(EEG, [-3000, 3000], 'ICA', [], [10, 6], 1);
disp('finished!')

%%
% theta band: 4 to 8 Hz; alpha: 9 ~ 13 Hz; low beta: 14 ~ 20 Hz; high beta: 21 ~ 30 Hz
printTrialMapsAxes(EEG, [-3000, 3000], 'ICA', [4, 8], [10, 6], 1, [-600, -200]);
suptitle('\theta 4-8 Hz')
printTrialMapsAxes(EEG, [-3000, 3000], 'ICA', [9, 13], [1, 1], 1, [-400, 0]);
suptitle('\alpha 8-13 Hz')
printTrialMapsAxes(EEG, [-3000, 3000], 'ICA', [14, 20], [10, 6], 1, [-600, -200]);
suptitle('low\beta 13-20 Hz')
printTrialMapsAxes(EEG, [-3000, 3000], 'ICA', [21, 30], [10, 6], 1, [-600, -200]);
suptitle('high\beta 20-30 Hz')


%%
% plot r > 0.1
disp(['Plotting ', num2str(length(big_rr)), ' pair(s) of ICs'])
for i = 1:20%length(big_rr)
    c1 = big_rr{i}(1);
    c2 = big_rr{i}(2);
    compareICA2(EEG, c1, c2, rr_max{c1, c2}(2:3), rr_max{c1, c2}(4:5))
end
