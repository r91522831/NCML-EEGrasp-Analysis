close all; clear; clc;
All_dirpath = uigetdir();
All_dirlist = dir(fullfile(All_dirpath, 'sub*'));

All_figpath = fullfile(All_dirpath, 'erp_figs');
if ~exist(All_figpath, 'dir')
    mkdir(All_figpath);
end

for All_i = 1:length(All_dirlist)
    clearvars -except All_*; close all;
    filepath = fullfile(All_dirpath, All_dirlist(All_i).name);
    filelist = dir(fullfile(filepath, '*.set'));
    keystr = 'pruned_ICA';
    filename = {filelist(contains({filelist.name}, keystr)).name};
    filename = filename{1}; % get the first file named as key string.
    subID = filename(1:6);
    
    % load data set
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
    EEG = pop_loadset('filename', filename, 'filepath', filepath);
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
    
    %%
    % plot erpimage w/o filtering
    
    
    %%
    % plot erpimage at theta bands 4 - 8 Hz
    EEG = pop_eegfiltnew(EEG, 'locutoff', 4, 'hicutoff', 8);
    
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1, 'setname', [EEG.setname, '_theta_4to8Hz'], 'savenew', fullfile(filepath, [EEG.setname, '_theta_4to8Hz.set']), 'gui', 'off');
    
    % setup figures
    h_fig1 = figure;
    h_fig2 = figure;
    h_fig3 = figure;
    
    freqName3 = '\theta';
    colID = 3; % column 1, 2, or 3
    figure(h_fig1);
    chanName3 = 'Cz'; % 'Fz' FCz'
    erpimage_subplot(EEG, freqName3, chanName3, colID);
    colID = 2; % column 1, 2, or 3
    chanName2 = 'Fz'; % 'Fz' FCz'
    erpimage_subplot(EEG, freqName3, chanName2, colID);
    colID = 1; % column 1, 2, or 3
    chanName1 = 'FCz'; % 'Fz' FCz'
    erpimage_subplot(EEG, freqName3, chanName1, colID);
    suptitle([subID, ' ', freqName3, ' ', chanName1, ' ', chanName2, ' ', chanName3]);
    filename_fig = [subID, '_', chanName1, '_', chanName2, '_', chanName3, '.fig'];
    savefig(h_fig1, fullfile(All_figpath, filename_fig),'compact');
    
    figure(h_fig2);
    colID = 3; % column 1, 2, or 3
    chanName32 = 'Fz'; % 'Fz' FCz'
    erpimage_subplot(EEG, freqName3, chanName32, colID);
    figure(h_fig3);
    chanName33 = 'FCz'; % 'Fz' FCz'
    erpimage_subplot(EEG, freqName3, chanName33, colID);
    
    % plot erpimage at alpha bands 9 - 12 Hz
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2, 'retrieve', 1, 'study', 0);
    EEG = pop_eegfiltnew(EEG, 'locutoff', 9, 'hicutoff', 12);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1, 'setname', [EEG.setname, '_alpha_9to12Hz'], 'savenew', fullfile(filepath, [EEG.setname, '_alpha_9to12Hz.set']), 'gui', 'off');
    
    freqName1 = '\alpha';
    colID = 1; % 1, 2, or 3
    figure(h_fig2);
    chanName12 = 'C3'; % 'C3' 'CP3'
    erpimage_subplot(EEG, freqName1, chanName12, colID);
    figure(h_fig3);
    chanName13 = 'CP3'; % 'C3' 'CP3'
    erpimage_subplot(EEG, freqName1, chanName13, colID);
    
    % plot erpimage at beta bands 20 - 30 Hz
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 3, 'retrieve', 1, 'study', 0);
    EEG = pop_eegfiltnew(EEG, 'locutoff', 20, 'hicutoff', 30);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1, 'setname', [EEG.setname, '_beta_20to30Hz'], 'savenew', fullfile(filepath, [EEG.setname, '_beta_20to30Hz.set']), 'gui', 'off');
    
    freqName2 = '\beta';
    colID = 2; % 1, 2, or 3
    figure(h_fig2);
    chanName22 = 'C1'; %'C1' 'C3'
    erpimage_subplot(EEG, freqName2, chanName22, colID);
    figure(h_fig3);
    chanName23 = 'C3'; %'C1' 'C3'
    erpimage_subplot(EEG, freqName2, chanName23, colID);
    
    figure(h_fig2);
    suptitle([subID, ' ', freqName1, ' ', chanName12, ' ', freqName2, ' ', chanName22, ' ', freqName3, ' ', chanName32]);
    figure(h_fig3);
    suptitle([subID, ' ', freqName1, ' ', chanName13, ' ', freqName2, ' ', chanName23, ' ', freqName3, ' ', chanName33]);
    
    filename_fig = [subID, '_', chanName12, '_', chanName22, '_', chanName32, '.fig'];
    savefig(h_fig2, fullfile(All_figpath, filename_fig),'compact');
    filename_fig = [subID, '_', chanName13, '_', chanName23, '_', chanName33, '.fig'];
    savefig(h_fig3, fullfile(All_figpath, filename_fig),'compact');
end

