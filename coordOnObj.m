function [output_coord] = coordOnObj(markers, holdSide)
%coordOnObj Summary of this function goes here
% coordinate fixed on the object for each frame!!!

obj_vector1 = markers{4, 1} - markers{1, 1}; % PS_3 - PS_0
obj_vector2 = markers{5, 1} - markers{1, 1}; % PS_4 - PS_0
tmp = cross(obj_vector1, obj_vector2);
coord_obj_x = tmp / norm(tmp);

obj_vector3 = markers{6, 1} - markers{3, 1}; % PS_5 - PS_2
obj_vector4 = markers{8, 1} - markers{3, 1}; % PS_7 - PS_2
tmp = cross(obj_vector3, obj_vector4);
coord_obj_y = tmp / norm(tmp);

tmp = cross(coord_obj_x, coord_obj_y);
coord_obj_z = tmp / norm(tmp);

% x coordinate of origin
tmp = cell2mat(markers(6:8, 1)); % PS_5, PS_6, PS_7
origin(1, 1) = mean(tmp(:, 1));
% y coordinate of origin
tmp = cell2mat(markers([1, 2, 4, 5], 1)); % PS_0, PS_1, PS_3, PS_4
origin(1, 2) = mean(tmp(:, 2));

% z coordinate of origin
% If the obj_side is R, the right handle is grasped.
% Set the handle coordinate at the center of the right handle
switch holdSide
    case 'R'
        tmp = cell2mat(markers([1, 2], 1)); % PS_0, PS_1
    case 'L'
        tmp = cell2mat(markers([4, 5, 6, 8], 1)); % PS_3, PS_4, PS_5, PS_7
    otherwise
        disp("Something wrong with holding side!")
end
origin(1, 3) = mean(tmp(:, 3));
tmp = [coord_obj_x; coord_obj_y; coord_obj_z; origin];
output_coord = array2table(tmp, 'RowNames', {'x_axis', 'y_axis', 'z_axis', 'origin'}, 'VariableNames', {'x', 'y', 'z'});
end

