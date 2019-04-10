close all; clear; clc;
base_folder = uigetdir();
file_list = dir(fullfile(base_folder, '*.set'));

%%
[ALLEEG, ~, ~, ~] = eeglab;
for i = 1:length(file_list)
    sub_id = file_list(i).name(1:4);
    
    EEG = pop_loadset('filename', file_list(i).name, 'filepath', base_folder);
    [ALLEEG, EEG, ~] = eeg_store( ALLEEG, EEG, 0 );
    EEG = eeg_checkset( EEG );
    EEG = pop_epoch( EEG, {  'onset'  }, [-1.5, 2.5], 'newname', [sub_id, '_eeg_csd epochs'], 'epochinfo', 'yes');
    [~, EEG, ~] = pop_newset(ALLEEG, EEG, 1, 'savenew', fullfile(base_folder, [sub_id, '_eeg_csd_epochs.set']), 'gui', 'off');
    EEG = eeg_checkset( EEG );
end