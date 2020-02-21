

% get data
nepoch = size(EEG.data, 3);
electrode_i = 41;
vol_t2l = nan(nepoch, 1);
vol_l2w = nan(nepoch, 1);
for ep_i = 1:nepoch
    tmp_touch_onset = contains(EEG.epoch(ep_i).eventtype, 'touch');
    [ ~, tmp_touch_id ] = min( abs( EEG.times - EEG.epoch(ep_i).eventlatency{tmp_touch_onset} ) );
    
    tmp_lft_onset = contains(EEG.epoch(ep_i).eventtype, 'onset');
    [ ~, tmp_lft_id ] = min( abs( EEG.times - EEG.epoch(ep_i).eventlatency{tmp_lft_onset} ) );
    
    contact2lft_time_id = tmp_touch_id:tmp_lft_id;
    
    [~, vol_t2l(ep_i, 1)] = robustcov(EEG.data(electrode_i, contact2lft_time_id, ep_i));
    
    [ ~, tmp_win_id ] = min( abs( EEG.times - 200 ) );
    lft2win_time_id = tmp_lft_id:tmp_win_id;
    [~, vol_l2w(ep_i, 1)] = robustcov(EEG.data(electrode_i, lft2win_time_id, ep_i));
end

%%
histogram(EEG.behavior.obj_roll_peak.peakRoll)
xlabel('peak roll ({\circ})')
title(['histogram for ', EEG.setname(1:6)])

%%
figure
subplot 121
plot(EEG.behavior.obj_roll_peak.peakRoll, vol_t2l, '*')
ylim([min([vol_t2l; vol_l2w]), max([vol_t2l; vol_l2w])])
xlabel('peak roll ({\circ})')
ylabel('voltage averaged from contact to lift ({\mu}V)')
subplot 122
plot(EEG.behavior.obj_roll_peak.peakRoll, vol_l2w, '*')
ylim([min([vol_t2l; vol_l2w]), max([vol_t2l; vol_l2w])])
ylabel('voltage averaged from lift to 200 ms ({\mu}V)')
xlabel('peak roll ({\circ})')
suptitle(['voltage against peak roll for ', EEG.setname(1:6)])