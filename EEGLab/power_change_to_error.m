close; clear; clc;

data_dir = uigetdir;
file_list = dir(fullfile(data_dir, '*.mat'));
sub_label = cellfun(@(x) str2double(x(2:4)), {file_list.name}');

sub_no = length(file_list);

t_tf_with_error = cell(sub_no, 1);
t_tf_complex_with_error = cell(sub_no, 1);
for i = 1:sub_no
    load(fullfile(data_dir, file_list(i).name));
    
    tmp = cell(size(t_tf_power));
    tmp_complex = cell(size(t_tf_power));
    for c = 1:width(pRoll)
        for f = 1:height(t_tf_power)
            tmp{f, c} = array2table([table2array(pRoll{1, c}{:}), t_tf_power{f, c}{:}], 'VariableNames', {'pRoll', 'powerChange'});
            tmp_complex{f, c} = array2table([table2array(pRoll{1, c}{:}), t_tf_power_complex{f, c}{:}], 'VariableNames', {'pRoll', 'powerChange'});
        end
    end
    
    t_tf_with_error{i, 1} = cell2table(tmp, 'VariableNames', t_tf_power.Properties.VariableNames, 'RowNames', t_tf_power.Properties.RowNames);
    t_tf_complex_with_error{i, 1} = cell2table(tmp_complex, 'VariableNames', t_tf_power_complex.Properties.VariableNames, 'RowNames', t_tf_power_complex.Properties.RowNames);
end
% t_tf_with_error = cell2table(tmp_tf_with_error, 'VariableNames', t_tf_power.Properties.VariableNames, 'RowNames', t_tf_power.Properties.RowNames);
% t_tf_complex_with_error = cell2table(tmp_tf_complex_with_error, 'VariableNames', t_tf_power_complex.Properties.VariableNames, 'RowNames', t_tf_power_complex.Properties.RowNames);



%% plot
for sub_id = 1%:length(sub_label)
    figure(sub_id)
    set(0,'defaultAxesFontSize', 18)
    
    color_maps = {winter, summer, copper, copper, copper};
    cond = {'IL', 'TR', 'PT1', 'PT2', 'PT3'};
    marker_style = {'x', 'o', '^', 's', 'd'};
    freq_band = {'alpha', 'theta', 'beta'};
    title_array = {'\alpha', '\theta', '\beta'};
    
    invert_obj_direct = 1;
    if mod(sub_label(sub_id), 2) == 0
        invert_obj_direct = -1;
    end
    
    for f = 1:length(freq_band)
        subplot(2, 2, f)
        hold on;
        
        for c = 1:length(cond)
            x = t_tf_complex_with_error{sub_id, 1}{freq_band{f}, cond{c}}{:}{:, 'pRoll'};
            y = abs(t_tf_complex_with_error{sub_id, 1}{freq_band{f}, cond{c}}{:}{:, 'powerChange'});
            % create color array
            color_map_source = flip(color_maps{c});
            color_step = floor(length(color_map_source)/length(x));
            color_map = color_map_source(1:color_step:length(color_map_source), :);
            
            for i = 1:length(x)
                plot(invert_obj_direct * x(i), y(i), marker_style{c}, 'MarkerEdgeColor', color_map(i, :));
%                 text(invert_obj_direct * x(i) + 0.3, y(i), num2str(i));
            end
        end
% % % % % % %         lsline
        hold off;
        ylabel('Power change (%)');
        xlabel('Peak roll ({\circ})');
        xlim([-25, 15]);
        ylim([-2, 8]);
        title(title_array{f}, 'FontWeight', 'bold');
    end
    
    mtit(['sub00', num2str(sub_label(sub_id))], 'yoff', .5)
    
    savefig(fullfile(data_dir, ['S00', num2str(sub_id)]));
end