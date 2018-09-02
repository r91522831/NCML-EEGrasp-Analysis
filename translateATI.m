function [alignedATI] = translateATI(rotatedATI, newZ)
%translateATI Summary of this function goes here
% aligned to the new z location with respective to the original location
% translate all rotated ATI to their respective handle center
% [1, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0; 0, 0, 1, 0, 0, 0; 0, newZ, 0, 1, 0, 0; -newZ, 0, 0, 0, 1, 0; 0, 0, 0, 0, 0, 0]
% rotatedATI's columns are fx, fy, fz, mx, my, mz
Tr = eye(6, 6);
Tr(4, 2) = newZ; % in the handle coordinate
Tr(5, 1) = -newZ;
alignedATI = (Tr * rotatedATI')';
end

