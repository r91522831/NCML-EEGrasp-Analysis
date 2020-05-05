
% Load EEG data
% update channel info
% remove 60Hz line noise

for i=1:2
    if i==1
        EEG = pop_loadbv('/R/BLAIS/NCML-SensoryGating/RAW/sub-01/eeg/SEP60Trial/', 'S_M_2019-11-20_14-08-40.vhdr', [],[]);
    else
        EEG = pop_loadbv('/R/BLAIS/NCML-SensoryGating/RAW/sub-01/eeg/SenoryGating50Trial/', 'S_M_2019-11-20_15-36-08.vhdr', [],[]);
    end

    
    %? why add AFZ and not use it??
    EEG=pop_chanedit(EEG, 'append',1,'changefield',{2 'labels' 'AFZ'});
    EEG=pop_chanedit(EEG, 'lookup','/R/MATLABPATH/eeglab_latest/plugins/dipfit3.0/standard_BESA/standard-10-5-cap385.elp');
    %EEG = pop_reref( EEG, [13 19] ,'refloc',struct('labels',{'AFZ'},'type',{''},'theta',{180},'radius',{0.12662},'X',{-32.9279},'Y',{-4.0325e-15},'Z',{78.363},'sph_theta',{-180},'sph_phi',{67.208},'sph_radius',{85},'urchan',{65},'ref',{''},'datachan',{0}),'exclude',64,'keepref','on');

    EEG = pop_cleanline(EEG, 'bandwidth', 2,'chanlist', [1:EEG.nbchan], 'computepower', 0, 'linefreqs', [60 120 180 240 300],...
    'normSpectrum', 0, 'p', 0.01, 'pad', 2, 'plotfigures', 0, 'scanforlines', 1, 'sigtype', 'Channels', 'tau', 100,...
    'verb', 1, 'winsize', 4, 'winstep', 4);

    pop_saveset(EEG,[pwd '/merged/ica/' EEG.comments([16:18 30:38]) '.set']);

end

% Load the two files
% merge
% filter 3-100Hz
% find bad channels
% epoch -0.5 to 1.0 s
    EEG = pop_loadset('filename',{'S_M_14-08-40.set' 'S_M_15-36-08.set'},'filepath','/R/BLAIS/NCML-SensoryGating/merged/ica/');
    EEG = pop_mergeset( ALLEEG, [1  2], 0);
    EEG  = pop_basicfilter( EEG,  1:EEG.nbchan , 'Cutoff', [ 3 100], 'Design', 'butter', 'Filter', 'bandpass', 'Order',  2, 'RemoveDC', 'on' ); % GUI: 10-Dec-2019 17:04:13
    EEG.BadChannels = getBadChannelIndices(EEG);
    
    
    
    EEG = pop_epoch( EEG, {  's4'  }, [-0.5           1], 'newname', '3-200Hz epochs', 'epochinfo', 'yes');
 
 % You'll need a NVIDIA graphics card and the following software
 % Linux/OSX  github.com/fraimondo/cudaica
 % Windows    gihhub.com/yhz-1995/cudaica_win
 % To run this extremely fast implementation of infoMAX
 % Alternatively, you could use binica (25X longer) or runica (~300+x
 % longer)
     EEG = pop_runica(EEG, 'extended', 1, 'icatype', 'cudaica', 'chanind',...
         setdiff(1:EEG.nbchan-1,EEG.BadChannels),'concatenate' ,'off', 'verbose', 'off');
     
 % label ICs
 % source localization
     EEG=pop_iclabel(EEG,'default');
 
 
     COREG=[0.83215 -15.6287 2.4114 0.081214 0.00093739 -1.5732 1.1742 1.0601 1.1485];
    EEG = pop_dipfit_settings( EEG, 'hdmfile','/R/MATLABPATH/eeglab2019_1/plugins/dipfit/standard_BEM/standard_vol.mat',...
                                'coordformat','MNI','mrifile','/R/MATLABPATH/eeglab2019_1/plugins/dipfit/standard_BEM/standard_mri.mat',...
                                'chanfile','/R/MATLABPATH/eeglab2019_1/plugins/dipfit/standard_BEM/elec/standard_1005.elc',...
                                'coord_transform',COREG ,'chansel',[1:EEG.nbchan-1] );
    EEG = pop_multifit(EEG, [1:EEG.nbchan-1] ,'threshold',40,'plotopt',{'normlen' 'on'});
    pop_dipplot( EEG, [1:EEG.nbchan-1] ,'mri','/R/MATLABPATH/eeglab2019_1/plugins/dipfit/standard_BEM/standard_mri.mat','normlen','on');
 
 
 %EEG = pop_rmbase( EEG, [-500    0]);
 pop_saveset(EEG,[pwd '/merged/ica/MarcoMergedICA_noHEO.set']);



end
