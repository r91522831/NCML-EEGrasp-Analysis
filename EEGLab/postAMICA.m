

[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;
EEG = pop_loadset('filename','S009_AMICA.set','filepath','/Users/yen-hsunwu/GitHub/testfolder/');
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
EEG = pop_loadmodout(EEG,'/Users/yen-hsunwu/GitHub/testfolder/amicaout/');
[ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = eeg_checkset( EEG );

%%
model_no = 10;
for i = 1:model_no
    pop_topohistplot(EEG, 0, 1:62, ['S009_Highpass1Hz, Model ', num2str(i)], [8, 8], 0, 'electrodes', 'on', 'showhist', 1, 'use_block', 1, 'model', i);
end

%%
EEG = pop_modprobplot(EEG, 1:10, 1, 2, 0, 'onset', 's129', 's17', 's33', 's65', 's9');

%%
EEG = pop_epoch( EEG, {'onset'}, [-7.02, 6.05], 'newname', 'S009_epochs', 'epochinfo', 'yes');

%%
figure; pop_moderp(EEG, 1, 2, [], 'FP1', 0, 1, {}, [], '', 'yerplabel', '\muV', 'cbar', 'on', 'caxis', [0.5 10.5] );

%%
EEG = pop_iclabel(EEG);

%%
ind_event = [round([EEG.event( strcmp({EEG.event.type}, 's9') ).latency]'), ...
             round([EEG.event( strcmp({EEG.event.type}, 's17') ).latency]'), ...
             round([EEG.event( strcmp({EEG.event.type}, 'onset') ).latency]'), ...
             round([EEG.event( strcmp({EEG.event.type}, 's33') ).latency]'), ...
             round([EEG.event( strcmp({EEG.event.type}, 's65') ).latency]'), ...
             round([EEG.event( strcmp({EEG.event.type}, 's129') ).latency]')];
ind_b4afonset = ceil([max(ind_event(:, 3) - ind_event(:, 1)) / EEG.srate, max(ind_event(:, end) - ind_event(:, 3)) / EEG.srate] * 100) / 100;

% win_width = 