% get data
nepoch = size(EEG.data, 3);
sub_id = EEG.setname(1:6);
%%
for electrode_i = 1:63
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
    
    electrode_label = EEG.chanlocs(electrode_i).labels;
    
    %%
    f = figure;
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
    suptitle([electrode_label, ' voltage against peak roll for ', sub_id])
    
    figdir = ['/Users/yenhsunw/Dropbox (ASU)/NCML-EEGrasp/for Feb21 2020/', sub_id];
    figname = ['peakRoll_voltage_', '_', electrode_label, '_', sub_id];
    savefig(f, fullfile(figdir, figname))
    saveas(f, fullfile(figdir, [figname, '.png']))
    close(f)
end
%%
f = figure;
histogram(EEG.behavior.obj_roll_peak.peakRoll);
xlabel('peak roll ({\circ})')
title(['histogram for ', EEG.setname(1:6)])

figdir = ['/Users/yenhsunw/Dropbox (ASU)/NCML-EEGrasp/for Feb21 2020/', sub_id];
figname = ['peakRoll_hist', '_', sub_id];
savefig(f, fullfile(figdir, figname))
saveas(f, fullfile(figdir, [figname, '.png']))
close(f)