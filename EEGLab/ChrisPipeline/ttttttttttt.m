
% load EEG data

% compute the 


c1 = mean(EEG.icaact(1, :, 1:19), 3);
c2 = mean(EEG.icaact(2, :, 1:19), 3);
[c1c2_r, c1c2_lag] = xcorr(c1, c2, 'coeff');




printTrialMapsAxes(EEG, [-3000, 3000], 'ICA', [], [10, 6], 1)




for e = 1:95
    c_trial(e) = corr(squeeze(EEG.icaact(7, :, e))', squeeze(EEG.icaact(27, :, e))');
end