close all; clearvars; clc

%% load aligned data
[filename, pathname, ~] = uigetfile;
load(fullfile(pathname, filename));

coord_table_x = [1, 0, 0];
coord_table_y = [0, 1, 0];
coord_table_z = [0, 0, 1];

%% lowpass filter all data
dt = diff(data{1}{1:2, 1}) * 0.001; % in second
cutoff = 30; % in Hz
data_filtered = data;
for i = 1:length(file_list)
    data_filtered{i}{:, 3:end} = filtmat_class( dt, cutoff, data{i}{:, 3:end} );
end

%%
resultant_mx = cell(size(file_list));
angTilt = cell(size(file_list));
coord_obj_y = cell(size(file_list));
coord_obj_z = cell(size(file_list));
ind_lft_onset = zeros(size(file_list));
obj_height = cell(size(file_list));

% assign filtered or raw data to analyze
input = data_filtered; % or data
for i = 1:length(file_list)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% get time
    time = input{i}{:, {'time_ms'}};
    % get trigger indices
    audio_trigger = input{i}{:, {'trigger'}};
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Get object coordinate before reaching initiation
    % get PS Marker6: the center of the object
    % object coordinate before audio go cue and will keep until object
    % lift*
    ind_b4go = (audio_trigger == 1);
    obj_marker0_b4go = mean(input{i}{ind_b4go, var_PS{1}}, 1);
    obj_marker1_b4go = mean(input{i}{ind_b4go, var_PS{2}}, 1);
    obj_marker2_b4go = mean(input{i}{ind_b4go, var_PS{3}}, 1);
    obj_marker4_b4go = mean(input{i}{ind_b4go, var_PS{5}}, 1);
    obj_marker5_b4go = mean(input{i}{ind_b4go, var_PS{6}}, 1);
    obj_marker6_b4go = mean(input{i}{ind_b4go, var_PS{7}}, 1);
    obj_marker7_b4go = mean(input{i}{ind_b4go, var_PS{8}}, 1);
    
    % coordinate fixed on the object for each frame!!!
    coord_obj_origin = obj_marker6_b4go;
    
    obj_vector1_b4go = obj_marker1_b4go - obj_marker4_b4go;
    obj_vector2_b4go = obj_marker0_b4go - obj_marker4_b4go;
    coord_obj_x_b4go = cross(obj_vector1_b4go, obj_vector2_b4go);
    
    obj_vector3_b4go = obj_marker5_b4go - obj_marker2_b4go;
    obj_vector4_b4go = obj_marker7_b4go - obj_marker2_b4go;
    % coord_obj_y_b4go = cross(obj_vector3_b4go, obj_vector4_b4go);
    
    % coord_obj_z_b4go = cross(coord_obj_x_b4go, coord_obj_y_b4go);
    
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    coord_obj_y_b4go = cross(obj_vector1_b4go, coord_obj_x_b4go);
    coord_obj_z_b4go = obj_vector1_b4go; % for now!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    % Angle b/w line (object z) and plane (table xz)
    angTilt_b4go = asind( abs(dot(coord_table_y, coord_obj_z_b4go)) / (sqrt(sum(coord_table_y.^2)) * sqrt(sum(coord_obj_z_b4go.^2))) );
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% compute frame by frame
    for time_id = 1:height(input{i})
        %% compute Torque compensation
        obj_side = file_list(i).name(end-4);
        switch obj_side
            case 'R'
                finger_Th = input{i}{time_id, var_ATI{1}};
                finger_V  = input{i}{time_id, var_ATI{3}};
            case 'L'
                finger_Th = input{i}{time_id, var_ATI{2}};
                finger_V  = input{i}{time_id, var_ATI{4}};
            otherwise
                disp('something wrong with the file name!!')
        end

        % find th_cop_y and v_cop_y
        % since the ATI reference origin is on tool surface: dy = (mx - fy * dz) / fz
        th_cop_z2surface = -5; % in mm
        v_cop_z2surface = 5; % in mm
        th_cop_y = (finger_Th(:, 4) - finger_Th(:, 2) .* th_cop_z2surface )./ finger_Th(:, 3);
        v_cop_y = (finger_V(:, 4) - finger_V(:, 2) .* v_cop_z2surface )./ finger_V(:, 3);

        % width b/w two nano 25 tool surface in mm
        width_obj = 2 * (21.6 + 3);

        diff_finger = finger_Th - finger_V;
        resultant_mx{i}(time_id, 1) = diff_finger(:, 2) .* width_obj - (diff_finger(:, 3) .* (th_cop_y - v_cop_y));
    
        %% compute object tilt
        obj_marker0 = input{i}{time_id, var_PS{1}};
        obj_marker1 = input{i}{time_id, var_PS{2}};
        obj_marker2 = input{i}{time_id, var_PS{3}};
        obj_marker4 = input{i}{time_id, var_PS{5}};
        obj_marker5 = input{i}{time_id, var_PS{6}};
        obj_marker7 = input{i}{time_id, var_PS{8}};
        
        obj_vector1 = obj_marker1 - obj_marker4;
        obj_vector2 = obj_marker0 - obj_marker4;
        coord_obj_x = cross(obj_vector1, obj_vector2);
        
        obj_vector3 = obj_marker5 - obj_marker2;
        obj_vector4 = obj_marker7 - obj_marker2;
        % coord_obj_y = cross(obj_vector3, obj_vector4);
        % coord_obj_z = cross(coord_obj_x, coord_obj_y);
        
        %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        coord_obj_y{i}(time_id, :) = cross(obj_vector1, coord_obj_x);
        coord_obj_z{i}(time_id, :) = obj_vector1; % for now!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        % Angle b/w line (object z) and plane (table xz)
        angTilt{i}(time_id, 1) = abs(angTilt_b4go - asind( abs(dot(coord_table_y, coord_obj_z{i}(time_id, :))) / (sqrt(sum(coord_table_y.^2)) * sqrt(sum(coord_obj_z{i}(time_id, :).^2))) ));
        
    end
    
    disp(i);

    %% compute finger tip coordinate without missing frames
    
end

save(fullfile(pathname, [filename(1:4), '_temp_result.mat']), 'resultant_mx', 'angTilt', 'ind_lft_onset', 'file_list', 'pathname', 'obj_height')

%% plot
mx = resultant_mx;
dt = diff(data{1}{1:2, 1}) * 0.001;
cutoff = 5; % in Hz
mx_filtered = mx;
obj_height_filtered = obj_height;
for i = 1:length(mx)
    mx_filtered{i} = filtmat_class( dt, cutoff, mx{i} );
    obj_height_filtered{i} = filtmat_class( dt, cutoff, obj_height{i} );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% find lift onset
for i = 1:length(file_list)
    %{
    lift_marker0 = zeros(height(input{i}), 1);
    for j = 1:height(input{i})
        obj_marker0 = input{i}{j, var_PS{1}};
        lift_marker0(j, 1) = sqrt(sum((obj_marker0 - obj_marker0_b4go).^2));
    end
    avg_lft = mean(lift_marker0(ind_b4go));
    std_lft = std(lift_marker0(ind_b4go));
    %     tmp = find(abs(lift_marker0 - avg_lft) > 5 * std_lft);
    obj_height{i} = abs(lift_marker0 - avg_lft);
    tmp = find(abs(lift_marker0 - avg_lft) > 10); % larger than 10 mm
    ind_lft_onset(i, 1) = tmp(1) - 1;
    %}
    
    
    ind_hold = find(data{i, 1}{:, 2} == 3);
    for j = ind_hold(1, 1):-1:1
        if obj_height_filtered{i, :}(j, 1) < 30 % 10 mm
            ind_lft_onset(i, 1) = j;
            break;
        end
    end
end

for i = 1:length(file_list)
    figure(2)
    subplot 211
    plot(obj_height_filtered{i, :})
    hold on
    vline(ind_lft_onset(i, 1));
    hold off
    subplot 212
    plot(mx_filtered{i, :})
    hold on
    vline(ind_lft_onset(i, 1));
    hold off
    
%     disp(i)
%     pause
end





%%

time_of_plot = ind_lft_onset(1, 1);
tmp = [str2double(file_list(1).name(7:9)), mx{1}(time_of_plot, 1), angTilt{1}(time_of_plot,1)];
session = cell(33, 2);
j = 1;
for i = 2:length(file_list)
    time_of_plot = ind_lft_onset(i, 1);
    if (file_list(i).name(:, 11) ~= file_list(i - 1).name(:, 11))
        session(j, :) = {file_list(i - 1).name(:, 11), tmp};
        j = j + 1;
        tmp = [];
    end
    
    tmp = [tmp; str2double(file_list(i).name(7:9)), mx{i}(time_of_plot, 1), angTilt{i}(time_of_plot, 1)];
end
session(j, :) = {file_list(end).name(:, 11), tmp};

figure(1)
subplot 211
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
subplot 212
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
ylim([0, 15])
legend({'IL', 'TR', 'PT'})

