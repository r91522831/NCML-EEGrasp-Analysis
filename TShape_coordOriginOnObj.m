function [output_coord] = TShape_coordOriginOnObj(markers, ids, basis)
%coordOnObj Summary of this function goes here
% coordinate fixed on the object for each frame!!!

% rigid body object sensor locations with respect to the right handle center
offset_R = [ 17  , -49.4256,    0  ; ...
             17  ,  49.4256,    0  ; ...
            -12  ,  62     ,   13.5; ...
             12  ,  62     , - 13.5 ];

% rotate the offsets back to lab frame
lab_frame = transpose( table2array(basis) );
offset_R_lab = offset_R * lab_frame;   

% compute origin from all good sensors for later average
origin_r = zeros(length(ids), 3);
for i = 1:length(ids)
    raw = markers{ids(i), 1};
    origin_r(i, :) = raw - offset_R_lab(i, :);
end
origin = mean(origin_r, 1);

output_coord = array2table([origin', origin', origin'], 'RowNames', {'x', 'y', 'z'}, 'VariableNames', {'origin', 'RCenter', 'LCenter'});
end

