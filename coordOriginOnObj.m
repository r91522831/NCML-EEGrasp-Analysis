function [output_coord] = coordOriginOnObj(markers, ids, basis, holdSide)
%coordOnObj Summary of this function goes here
% coordinate fixed on the object for each frame!!!

% rigid body object sensor locations with respect to the right handle center
offset_R = [ 17  , -49.4256,    0  ; ...
             17  ,  49.4256,    0  ; ...
             11.5,  62     ,   14  ; ...
             17  , -49.4256, -115  ; ...
             17  ,  49.4256, -115  ; ...
             12  ,  62     , -128.5; ...
              0  , -58     , - 57.5; ...
            -12  ,  62     , -101.5 ];
% rigid body object sensor locations with respect to the left handle center
offset_L = [ 17  , -49.4256,  115  ; ...
             17  ,  49.4256,  115  ; ...
             11.5,  62     ,  129  ; ...
             17  , -49.4256,    0  ; ...
             17  ,  49.4256,    0  ; ...
             12  ,  62     , - 13.5; ...
              0  , -58     ,   57.5; ...
            -12  ,  62     ,   13.5 ];
% rotate the offsets back to lab frame
lab_frame = transpose( table2array(basis) );
offset_R_lab = offset_R * lab_frame;
offset_L_lab = offset_L * lab_frame;       

% compute origin from all good sensors for later average
origin_r = zeros(length(ids), 3);
origin_l = zeros(length(ids), 3);
for i = 1:length(ids)
    raw = markers{ids(i), 1};
    origin_r(i, :) = raw - offset_R_lab(ids(i), :);
    origin_l(i, :) = raw - offset_L_lab(ids(i), :);
end
origin_r_avg = mean(origin_r, 1);
origin_l_avg = mean(origin_l, 1);

switch holdSide
    case 'R'
        origin = origin_r_avg;
    case 'L'
        origin = origin_l_avg;
    otherwise
        disp("Something wrong with holding side!")
end

output_coord = array2table(origin', 'RowNames', {'x', 'y', 'z'}, 'VariableNames', {'origin'});
end

