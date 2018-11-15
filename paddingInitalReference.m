% Apply average reference after adding initial reference
EEG.nbchan = EEG.nbchan + 1;
EEG.data(end + 1, :) = zeros(1, EEG.pnts); % add an all zero channel as the reference channel
EEG.chanlocs(1, EEG.nbchan).labels = 'initialReference'; % label initial reference channel
EEG = pop_reref(EEG, []); % do re-referncing
EEG = pop_select(EEG, 'nochannel', {'initialReference'}); % remove the initial reference channel