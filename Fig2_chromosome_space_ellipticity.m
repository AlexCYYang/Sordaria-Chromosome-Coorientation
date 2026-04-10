% Calculate ellipticity of chromosome space
% Ellipticity is defined as the long axis length divided by the average 
% of the two short axis lengths of the fitted ellipsoid.
%
% External dependency:
% This script requires 'ellipsoid_fit_new.m' (a third-party function).
% Please ensure it is downloaded and added to your MATLAB path.
% Yury (2026). Ellipsoid fit 
% (https://www.mathworks.com/matlabcentral/fileexchange/24693-ellipsoid-fit), 
% MATLAB Central File Exchange. Retrieved April 5, 2026.
clear;
close all;

% Load image dimensions for converting pixel indices to physical coordinates
load('obl_info.mat');

% Specify the input TIFF stack
% Example: 'TCV.tif', 'RCV.tif', 'LCV.tif', or 'DCV.tif'
tif_filename = 'LCV.tif';

% Read the 3D binary stack
A = [];

for j = 1:tz
    A(:, :, j) = imread(tif_filename, j);
end

% Extract voxel coordinates from the binary mask
[row_idx, col_idx, z_idx] = ind2sub(size(A), find(A > 1));

x = col_idx * px;
y = row_idx * py;
z = (z_idx - 1) * pz;

% Fit an ellipsoid to the 3D voxel coordinates
[center, radii, evecs, v, chi2] = ellipsoid_fit_new([x, y, z]);

% Sort semi-axis lengths from longest to shortest
r = sort(radii(:), 'descend');

% Ellipticity: long axis / average of the two short axes
elp = r(1) / ((r(2) + r(3)) / 2);

fprintf('Ellipticity = %.4f\n', elp);