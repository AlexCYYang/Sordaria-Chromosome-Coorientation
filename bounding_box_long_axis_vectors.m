% Calculate chromosome long-axis vectors (ve) for coorientation analysis
%
% This script computes one normalized long-axis direction vector for each
% traced chromosome path in a nucleus. For each chromosome, a minimum
% bounding box is fitted to its 3D coordinates, and the longest box axis
% is taken as the chromosome orientation vector. The output variable 've'
% is then used for downstream coorientation analysis.
%
% Input:
%   C : cell array
%       Each C{i} is an N-by-3 matrix containing the XYZ coordinates
%       of one traced chromosome path in a nucleus.
%
% Output:
%   ve : nChr-by-3 matrix
%       Each row of ve is the normalized long-axis vector of the minimum
%       bounding box fitted to one chromosome.
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


nChr = length(C);% Number of traced chromosomes in this nucleus
ve = nan(nChr, 3);

for i = 1:nChr
    c1 = C{i};    % Extract the 3D coordinates of chromosome i

% Fit a minimum bounding box to the 3D chromosome coordinates
    % cornerpoints contains the 8 vertices of the fitted box
    [~, cornerpoints, ~, ~, ~] = minboundbox(c1(:,1), c1(:,2), c1(:,3));

    % Define the three box-edge vectors originating from one corner
    % These represent the three principal box axes
    ax1 = cornerpoints(4,:) - cornerpoints(1,:);
    ax2 = cornerpoints(5,:) - cornerpoints(1,:);
    ax3 = cornerpoints(2,:) - cornerpoints(1,:);
    
    
    edge_lengths = [norm(ax1), norm(ax2), norm(ax3)];% Compute the lengths of the three box axes
    axes_vec = [ax1; ax2; ax3];% Store the three candidate axis vectors together


    [~, idx_long] = max(edge_lengths); % Identify the longest axis of the bounding box

    % Normalize the longest axis vector to unit length
    % This removes length information and keeps only direction
    ve(i,:) = axes_vec(idx_long,:) / norm(axes_vec(idx_long,:));
end

% 've' can then be used to calculate pairwise inter-chromosome angles
% for coorientation analysis.