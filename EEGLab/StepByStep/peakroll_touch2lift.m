

% get data
electrode_i = 41;

for ep_i = 1:size(EEG.datat, 3)
    tmp_touch_id = contains(EEG.epoch(ep_i).eventtype, 'touch');
    contact2lft_time_id = EEG.times(EEG.times == EEG.epoch(ep_i).eventlatency{tmp_touch_id});
    
    vol{ep_i, 1} = EEG.data(electrode_i, contact2lft_time_id:, ep_i);
end