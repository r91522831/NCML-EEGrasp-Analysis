%% Plot peak mx and peak roll around lift onset in all trials
close all; clearvars; clc

%% load aligned data
pathname_plot = uigetdir;
filelist = dir(fullfile(pathname_plot, '*.mat'));

for sub = 1:length(filelist)
    clearvars -except sub filelist pathname_plot
    close all
    load(fullfile(pathname_plot, filelist(sub).name));
    
    tmp = [str2double(file_list(1).name(7:9)), mx_onset(1), peak_roll{1, 'peakRoll'}];%, peak_mx{1, 'peakMx'}, peak_roll{1, 'peakRoll'}];
    session = cell(33, 2);
    j = 1;
    for i = 2:length(file_list)
        time_of_plot = ind_lft_onset(i, 1);
        if (file_list(i).name(:, 11) ~= file_list(i - 1).name(:, 11))
            session(j, :) = {file_list(i - 1).name(:, 11), tmp};
            j = j + 1;
            tmp = [];
        end

        tmp = [tmp; str2double(file_list(i).name(7:9)), mx_onset(i), peak_roll{i, 'peakRoll'}];%, peak_mx{i, 'peakMx'}, peak_roll{i, 'peakRoll'}];
    end
    session(j, :) = {file_list(end).name(:, 11), tmp};

    figure(1)
    subplot 411
    hold on
    for i = 1:length(session)
        switch session{i, 1}
            case 'I'
                line_spec = '-*r';
            case 'T'
                line_spec = '-ob';
            case 'P'
                line_spec = '-xk';
            otherwise
        end
        plot(session{i, 2}(:, 1), session{i, 2}(:, 2), line_spec)
    end
    hold off

    subplot 412
    hold on
    for i = 1:length(session)
        switch session{i, 1}
            case 'I'
                line_spec = '-*r';
            case 'T'
                line_spec = '-ob';
            case 'P'
                line_spec = '-xk';
            otherwise
        end
        plot(session{i, 2}(:, 1), session{i, 2}(:, 3), line_spec)
    end
    hold off

    subplot 413
    hold on
    for i = 1:length(session)
        switch session{i, 1}
            case 'I'
                line_spec = '-*r';
            case 'T'
                line_spec = '-ob';
            case 'P'
                line_spec = '-xk';
            otherwise
        end
        plot(session{i, 2}(:, 1), abs(session{i, 2}(:, 2)), line_spec)
    end
    hold off

    subplot 414
    hold on
    for i = 1:length(session)
        switch session{i, 1}
            case 'I'
                line_spec = '-*r';
            case 'T'
                line_spec = '-ob';
            case 'P'
                line_spec = '-xk';
            otherwise
        end
        plot(session{i, 2}(:, 1), abs(session{i, 2}(:, 3)), line_spec)
    end
    hold off
    % ylim([0, 15])
    legend({'IL', 'TR', 'PT'})
    mtit(filelist(sub).name(1:4))
    savefig(fullfile(pathname_plot, [filelist(sub).name(1:4), '_onset']))
end