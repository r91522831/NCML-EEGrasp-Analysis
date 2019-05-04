% Set study design info for LIMO
orgSTUDY = STUDY;
STUDY = std_maketrialinfo(STUDY, ALLEEG);
EEG = eeg_checkset( EEG );
[STUDY, EEG] = pop_savestudy( STUDY, EEG, 'savemode','resave');

