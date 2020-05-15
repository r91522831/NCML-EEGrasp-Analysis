% This example code compares PSD in dB (left) vs. uV^2/Hz (right) rendered as scalp topography (setfile must be loaded.)
lowerFreq  = 4; % Hz
higherFreq = 8; % Hz
meanPowerDb     = zeros(60, 1);%EEG.nbchan,1);
meanPowerMicroV = zeros(60, 1);%EEG.nbchan,1);
for channelIdx = 1:60%EEG.nbchan
    [psdOutDb(channelIdx,:), freq] = spectopo(reshape(EEG.icaact(channelIdx, :, :), 60, []), 0, EEG.srate, 'plot', 'off');
    lowerFreqIdx    = find(freq==lowerFreq);
    higherFreqIdx   = find(freq==higherFreq);
    meanPowerDb(channelIdx) = mean(psdOutDb(channelIdx, lowerFreqIdx:higherFreqIdx));
    meanPowerMicroV(channelIdx) = mean(10.^((psdOutDb(channelIdx, lowerFreqIdx:higherFreqIdx))/10), 2);
end
 
figure
subplot(1,2,1)
topoplot(meanPowerDb, EEG.chanlocs)
title('Theta band (4-8Hz) power distribution')
cbarHandle = colorbar;
set(get(cbarHandle, 'title'), 'string', '(dB)')
 
subplot(1,2,2)
topoplot(meanPowerMicroV, EEG.chanlocs)
title('Theta band (4-8Hz) power distribution')
cbarHandle = colorbar;
set(get(cbarHandle, 'title'), 'string', '(uV^2/Hz)')