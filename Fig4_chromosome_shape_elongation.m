% Calculate chromosome elongation (EL)
% EL is defined as the mean long-axis length of the minimum bounding boxes
% fitted to the traced chromosome paths in a nucleus.
%
% In this analysis, each traced chromosome path is represented by a set of
% 3D coordinates. A minimum bounding box is fitted to each chromosome path,
% and the longest axis of that box is taken as the chromosome long-axis length.
% EL for the nucleus is then calculated as the mean of these long-axis lengths
% across all traced chromosomes.
%
% External dependency:
% This script requires 'minboundbox.m' (a third-party function).
% Please ensure 'minboundbox.m' is downloaded and added to your MATLAB path.
% Johannes Korsawe (2026). Minimal Bounding Box 
% (https://www.mathworks.com/matlabcentral/fileexchange/18264-minimal-bounding-box), 
% MATLAB Central File Exchange. Retrieved April 5, 2026.

clearvars -except C

% C is a cell array.
% Each C{i} contains an N-by-3 matrix of XYZ coordinates
% for one traced chromosome path in a nucleus.

% Expansion factor used to rescale measurements from expanded samples
% back to original-scale units.
expansion_factor = 7.33;

% Number of traced chromosomes in the nucleus
nChr = length(C);

% Preallocate an array to store the long-axis length of the minimum
% bounding box for each chromosome.
long_l = nan(nChr,1);

for i = 1:nChr
    % Extract the 3D coordinates of chromosome i
    c1 = C{i};

    % Fit a minimum bounding box to the 3D chromosome coordinates.
    % 'cornerpoints' contains the 8 vertices of the fitted box.
    [~, cornerpoints, ~, ~, ~] = minboundbox(c1(:,1), c1(:,2), c1(:,3));

    % Define the three box-edge vectors originating from one corner.
    % These correspond to the three axes of the minimum bounding box.
    ax1 = cornerpoints(4,:) - cornerpoints(1,:);
    ax2 = cornerpoints(5,:) - cornerpoints(1,:);
    ax3 = cornerpoints(2,:) - cornerpoints(1,:);

    % Compute the lengths of the three box axes.
    edge_lengths = [norm(ax1), norm(ax2), norm(ax3)];

    % Take the longest box axis as the chromosome long-axis length.
    % Divide by the expansion factor to convert to original-scale units.
    long_l(i) = max(edge_lengths) / expansion_factor;
end

% Calculate chromosome elongation (EL) for the nucleus as the mean
% long-axis length across all traced chromosomes.
EL = mean(long_l)