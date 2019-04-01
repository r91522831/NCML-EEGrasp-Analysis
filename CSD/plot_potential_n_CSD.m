close all; clear; clc;

[eeg_dataset_file, eeg_dataset_folder, ~]= uigetfile('*.set');

EEG = pop_loadset('filename', eeg_dataset_file, 'filepath', eeg_dataset_folder);
%{
for i = 1%:nb_epoch
    figure
    for t = 1:100:length(EEG.dataCSD(:, :, i))
        subplot 231
        topoplot(EEG.dataCSD(:, t, i), EEG.chanlocs, 'style', 'map', 'electrodes', 'off');%, 'maplimits', 'maxmin');
        title('CSD')
        subplot 232
        topoplot(EEG.data(:, t, i), EEG.chanlocs, 'style', 'map', 'electrodes', 'off');%, 'maplimits', 'maxmin');
        title('potential')
        subplot 233
        plot(EEG.behavior.obj_epoch{i, 1}.time, EEG.behavior.obj_epoch{i, 1}.roll)
        vline(EEG.behavior.obj_epoch{i, 1}.time(round(t * EEG.behavior.behavior_srate / EEG.srate)), '-b')
        vline(0, ':r')
        xlabel('time (ms)')
        ylabel('object roll ({\circ})')
        %             pause(.1)
        
        pause
    end
end
%}
%%
id_epoch = [1, 2, 3, 17, 18, 19, 20, 21, 22, 24];
nb_epoch = length(id_epoch);
figure
for t = 1:100:length(EEG.dataCSD(:, :, 1))
    for i = 1:length(id_epoch)
        subplot(nb_epoch, 3, 3 * (i - 1) + 1)
        topoplot(EEG.dataCSD(:, t, id_epoch(i)), EEG.chanlocs, 'style', 'map', 'electrodes', 'off');%, 'maplimits', 'maxmin');
        if i == 1
            title('CSD')
        end
        subplot(nb_epoch, 3, 3 * (i - 1) + 2)
        topoplot(EEG.data(:, t, id_epoch(i)), EEG.chanlocs, 'style', 'map', 'electrodes', 'off');%, 'maplimits', 'maxmin');
        if i == 1
            title('potential')
        end
        subplot(nb_epoch, 3, 3 * (i - 1) + 3);
        plot(EEG.behavior.obj_epoch{id_epoch(i), 1}.time, EEG.behavior.obj_epoch{id_epoch(i), 1}.roll);
        ylim([-0.5, 0.5])
        vline(EEG.behavior.obj_epoch{id_epoch(i), 1}.time(round(t * EEG.behavior.behavior_srate / EEG.srate)), '-b');
        vline(0, ':r');
        if i == nb_epoch
            xlabel('time (ms)');
        else
            xticklabels([]);
        end
        ylabel('object roll ({\circ})');
        
    end
    pause
    
end