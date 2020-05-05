function badChannels=getBadChannelIndices(EEG)
% Helper function:
% computes SD for each electrodes
% considers electrode "bad" if the |z-score| of that SD is >3.5
badChannels=0;
if EEG.nbchan < 40
    values=1:30;
else
    values=1:EEG.nbchan-1;
end

elect=setdiff(values,badChannels);
go=1;
TESTDATA=EEG.data(elect,:)';
while go == 1
    x=std(TESTDATA);
    x(isnan(x))=0;
    index=find(zscore(x)>3.5);
    if numel(index) > 0
        badChannels=[badChannels index];
        TESTDATA(:,index)=1;
    else
        go=0;
    end
end
badChannels=unique(badChannels(2:end));
end