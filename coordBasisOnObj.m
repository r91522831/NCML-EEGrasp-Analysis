function [output_coord] = coordBasisOnObj(markers, yz_id, xz_id)
%coordBasisOnObj Summary of this function goes here
% coordinate fixed on the object for each frame!!!

% use makers on yz plan to get x axis unit vector
obj_vector1 = markers{yz_id(2), 1} - markers{yz_id(1), 1}; % PS_3 - PS_0 or PS_3 - PS_1 or PS_1 - PS_0
obj_vector2 = markers{yz_id(3), 1} - markers{yz_id(1), 1}; % PS_4 - PS_0 or PS_4 - PS_1 or PS_3 - PS_0
tmp = cross(obj_vector1, obj_vector2);
coord_obj_x = tmp / norm(tmp);

% use makers on xz plan to get y axis unit vector
obj_vector3 = markers{xz_id(2), 1} - markers{xz_id(1), 1}; % PS_5 - PS_2
obj_vector4 = markers{xz_id(3), 1} - markers{xz_id(1), 1}; % PS_7 - PS_2
tmp = cross(obj_vector3, obj_vector4);
coord_obj_y = tmp / norm(tmp);

% use x and y axes to get z axis unit vector
tmp = cross(coord_obj_x, coord_obj_y);
coord_obj_z = tmp / norm(tmp);

tmp = [coord_obj_x; coord_obj_y; coord_obj_z];
output_coord = array2table(tmp, 'RowNames', {'x_axis', 'y_axis', 'z_axis'}, 'VariableNames', {'x', 'y', 'z'});
end

