close; clear; clc;

data_dir = uigetdir;
file_list = dir(fullfile(data_dir, '*_pruned_ICA.set'));
sub_id = cellfun(@(x) x(1:4), {file_list.name}', 'UniformOutput', false);
nb_sub = length(sub_id);

cond_names = {'ALL', 'IL', 'TR', 'PT1', 'PT2', 'PT3'};
cond_nb = length(cond_names);

% typeproc - type of processing: 1 process the raw channel data
%                                0 process the ICA component data
typeproc = str2double(input('\nChoose type of processing [1: raw, 0: component]: (default: 0) ', 's'));

for sub_i = 1:nb_sub
    [ALLEEG, ~, ~, ALLCOM] = eeglab;
    EEG = pop_loadset('filename', file_list(sub_i).name, 'filepath', data_dir);
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
    
    %% Time-freq analysis on channels    
    % typeproc - type of processing: 1 process the raw channel data
    %                                0 process the ICA component data
    if isnan(typeproc)
        typeproc = 0;
        nb_compoent = EEG.nbic;
        topovec_value = EEG.icawinv;
        cap_str = cell(nb_compoent, 1);
        for i = 1:nb_compoent
            cap_str{i} = [' IC ', num2str(i)];
        end
    elseif typeproc == 1
        nb_compoent = EEG.nbchan;
        topovec_value = 1:nb_compoent;
        cap_str = {EEG.chanlocs.labels}';
    end
    tf_ersp = cell(nb_compoent, 1);
    tf_itc = cell(nb_compoent, 1);
    tf_powbase = cell(nb_compoent, 1);
    tf_data = cell(nb_compoent, 1);
    
% % %     fig_dir = fullfile(output_dir, 'figures', sub_id);
% % %     if ~isfolder(fig_dir)
% % %         mkdir(fig_dir);
% % %     end
    
    EEG_cond = EEG;
    time_window = [-3000, 3000];
    for i = 1:nb_compoent
        h = figure;
        [tf_ersp{i, 1}, tf_itc{i, 1}, tf_powbase{i, 1}, tf_times, tf_freqs, ~, ~, tf_data{i, 1}] = ...
            pop_newtimef( EEG_cond, typeproc, i, time_window, [3, 0.5], 'topovec', topovec_value, ...
                                     'elocs', EEG.chanlocs, 'chaninfo', EEG.chaninfo, 'caption', [cond_names{1}, cap_str{i}], ...
                                     'freqs', [0, 35], 'baseline', [-600, -100], 'plotphase', 'off', 'scale', 'abs', 'padratio', 1 ); %'basenorm', 'on', 'trialbase', 'full');
        
% % %         savefig(h, fullfile(fig_dir, [sub_id, '_fig_', cond_names{1}, '_IC_', num2str(i)]));
        close(h);
        disp([' channel ', num2str(i)]);
    end
    
    save(fullfile(fileparts(data_dir), [sub_id{sub_i}, '_tf_info']), 'tf_data', 'tf_ersp', 'tf_itc', 'tf_powbase', 'tf_times', 'tf_freqs', '-v7.3');
end