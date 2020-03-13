function [cop] = centerOfPressure(fgr, z_offset, force_resolution)
%centerOfPressure Summary of this function goes here
%   Detailed explanation goes here

cop_z = z_offset * ones(height(fgr), 1);
% COPy = (fy * COPz + mx) / fz
cop_y = ((fgr.fy .* cop_z) + fgr.mx) ./ fgr.fz;
% COPx = (my - fx * COPz) / fz
cop_x = (fgr.my - (fgr.fx .* cop_z)) ./ fgr.fz;

% remove cop estimation before contact
cop_y(abs(fgr.fz) < force_resolution) = 0;
cop_x(abs(fgr.fz) < force_resolution) = 0;
cop_z(abs(fgr.fz) < force_resolution) = 0;

cop = array2table([cop_x, cop_y, cop_z], 'VariableNames', {'COPx', 'COPy', 'COPz'});
end

