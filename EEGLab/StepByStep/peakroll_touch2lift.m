disp('Select the project folder in the EEG directory:')
All_data = uigetdir;
All_data_list = dir(fullfile(All_data, 'sub-*'));
All_root_dir = fileparts(fileparts(All_data));

disp([num2cell((1:length(All_data_list))'), {All_data_list.name}']);
All_selected_sub = input('Which subject(s) to preprocess EEG? ');
if isempty(All_selected_sub)
    All_selected_sub = 1:length(All_data_list);
end

for All_i = All_selected_sub
    clearvars -except All_*;
    filepath = fullfile(All_data_list(All_i).folder, All_data_list(All_i).name);
    filelist = dir(fullfile(filepath, '*.set'));
    keystr = 'combine_behavior';
    filename = {filelist(contains({filelist.name}, keystr)).name};
    filename = filename{1}; % get the first file named as key string.
    subID = filename(1:6);
    
    % load 
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
    EEG = pop_loadset('filename', filename, 'filepath', filepath);
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
    
    % get data
    nepoch = size(EEG.data, 3);
    sub_id = EEG.setname(1:6);
    %%
    vol_t2l = nan(nepoch, 63);
    vol_l2w = nan(nepoch, 63);
    for electrode_i = 1:63
        for ep_i = 1:nepoch
            tmp_touch_onset = contains(EEG.epoch(ep_i).eventtype, 'touch');
            [ ~, tmp_touch_id ] = min( abs( EEG.times - EEG.epoch(ep_i).eventlatency{tmp_touch_onset} ) );
            
            tmp_lft_onset = contains(EEG.epoch(ep_i).eventtype, 'onset');
            [ ~, tmp_lft_id ] = min( abs( EEG.times - EEG.epoch(ep_i).eventlatency{tmp_lft_onset} ) );
            
            contact2lft_time_id = tmp_touch_id:tmp_lft_id;
            [~, vol_t2l(ep_i, electrode_i)] = robustcov(EEG.data(electrode_i, contact2lft_time_id, ep_i));
            
            [ ~, tmp_win_id ] = min( abs( EEG.times - 200 ) );
            lft2win_time_id = tmp_lft_id:tmp_win_id;
            [~, vol_l2w(ep_i, electrode_i)] = robustcov(EEG.data(electrode_i, lft2win_time_id, ep_i));
        end
        
        electrode_label = EEG.chanlocs(electrode_i).labels;
        
        %{
    %%
    f = figure;
    subplot 121
    plot(abs(EEG.behavior.obj_roll_peak.peakRoll), vol_t2l, '*')
    ylim([min([vol_t2l; vol_l2w]), max([vol_t2l; vol_l2w])])
    xlabel('peak roll ({\circ})')
    ylabel('voltage averaged from contact to lift ({\mu}V)')
    subplot 122
    plot(abs(EEG.behavior.obj_roll_peak.peakRoll), vol_l2w, '*')
    ylim([min([vol_t2l; vol_l2w]), max([vol_t2l; vol_l2w])])
    ylabel('voltage averaged from lift to 200 ms ({\mu}V)')
    xlabel('peak roll ({\circ})')
    suptitle([electrode_label, ' voltage against peak roll for ', sub_id])
    
    figdir = ['/Users/yen-hsunwu/Dropbox (ASU)/NCML-EEGrasp/for Feb21 2020/', sub_id];
    figname = ['peakRoll_voltage_', '_', electrode_label, '_', sub_id];
    savefig(f, fullfile(figdir, figname))
    saveas(f, fullfile(figdir, [figname, '.png']))
    close(f)
        %}
    end
    
    tmp_path = fullfile(All_data, 'vol');
    if ~exist(tmp_path, 'dir')
        mkdir(tmp_path);
    end
    tmp_filename = fullfile(tmp_path, [subID, '_vol']);
    save(tmp_filename, 'vol_t2l', 'vol_l2w');
end
chanlocs = EEG.chanlocs;
save(fullfile(tmp_path, 'chanlocs'), 'chanlocs')

%%
%{
f = figure;
histogram(EEG.behavior.obj_roll_peak.peakRoll);
xlabel('peak roll ({\circ})')
title(['histogram for ', EEG.setname(1:6)])

figdir = ['/Users/yen-hsunwu/Dropbox (ASU)/NCML-EEGrasp/for Feb21 2020/', sub_id];
figname = ['peakRoll_hist', '_', sub_id];
savefig(f, fullfile(figdir, figname))
saveas(f, fullfile(figdir, [figname, '.png']))
close(f)
%}