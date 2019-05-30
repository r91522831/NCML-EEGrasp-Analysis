smooth_factor = 5;

[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
EEG = pop_loadset('filename','sub-02_eeg_csd.set','filepath','/Users/yenhsunw/Dropbox (ASU)/NCML-EEGrasp/EEG/eeglab/001 Process/sub-02_onset/');
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
EEG = pop_eegfiltnew(EEG, 'locutoff',4,'hicutoff',8,'chantype',{'EEG'});
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname','sub-02_eeg_csd_theta4to8Hz','savenew','/Users/yenhsunw/Dropbox (ASU)/NCML-EEGrasp/EEG/eeglab/001 Process/sub-02_onset/sub-02_eeg_csd_theta4to8Hz.set','gui','off'); 
EEG = eeg_checkset( EEG );
figure; pop_erpimage(EEG,1, [6],[[]],'FZ', smooth_factor, 1, {}, [], '','yerplabel','\muV','erp','on','cbar','on','topo', { [6] EEG.chanlocs EEG.chaninfo } );
EEG = eeg_checkset( EEG );
figure; pop_erpimage(EEG,1, [41],[[]],'FCz', smooth_factor, 1, {}, [], '','yerplabel','\muV','erp','on','cbar','on','topo', { [41] EEG.chanlocs EEG.chaninfo } );
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'retrieve',1,'study',0); 
EEG = pop_eegfiltnew(EEG, 'locutoff',9,'hicutoff',12,'chantype',{'EEG'});
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname','sub-02_eeg_csd_alpha9to12Hz','savenew','/Users/yenhsunw/Dropbox (ASU)/NCML-EEGrasp/EEG/eeglab/001 Process/sub-02_onset/sub-02_eeg_csd_alpha9to12Hz.set','gui','off'); 
EEG = eeg_checkset( EEG );
figure; pop_erpimage(EEG,1, [15],[[]],'C3',smooth_factor,1,{},[],'' ,'yerplabel','\muV','erp','on','cbar','on','topo', { [15] EEG.chanlocs EEG.chaninfo } );
EEG = eeg_checkset( EEG );
figure; pop_erpimage(EEG,1, [47],[[]],'CP3',smooth_factor,1,{},[],'' ,'yerplabel','\muV','erp','on','cbar','on','topo', { [47] EEG.chanlocs EEG.chaninfo } );
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 3,'retrieve',1,'study',0); 
EEG = pop_eegfiltnew(EEG, 'locutoff',20,'hicutoff',30);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname','sub-02_eeg_csd_beta20to30Hz','savenew','/Users/yenhsunw/Dropbox (ASU)/NCML-EEGrasp/EEG/eeglab/001 Process/sub-02_onset/sub-02_eeg_csd_beta20to30Hz.set','gui','off'); 
EEG = eeg_checkset( EEG );
figure; pop_erpimage(EEG,1, [15],[[]],'C3',smooth_factor,1,{},[],'' ,'yerplabel','\muV','erp','on','cbar','on','topo', { [15] EEG.chanlocs EEG.chaninfo } );
EEG = eeg_checkset( EEG );
figure; pop_erpimage(EEG,1, [44],[[]],'C1',smooth_factor,1,{},[],'' ,'yerplabel','\muV','erp','on','cbar','on','topo', { [44] EEG.chanlocs EEG.chaninfo } );
EEG = eeg_checkset( EEG );



hline(find(strcmp({EEG.epoch.cond}, 'TR')), '-b')