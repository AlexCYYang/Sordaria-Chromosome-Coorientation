clearvars -except out;

% Load two 7x3 matrices:
% each row is the long-axis direction vector of one chromosome,
% obtained from the minimum bounding box.
% The two matrices correspond to the two partner nuclei.
load('R_cube_long_only.mat');%load; name of nucleus 1
chromosome_axis_vectors_R = ve;

load('L_cube_long_only.mat');%load; name of nucleus 2
chromosome_axis_vectors_L = ve;

% Normalize each chromosome axis vector to unit length
axis_lengths_R = sqrt(sum(chromosome_axis_vectors_R.^2, 2));
axis_lengths_L = sqrt(sum(chromosome_axis_vectors_L.^2, 2));

chromosome_axis_unit_vectors_R = chromosome_axis_vectors_R ./ axis_lengths_R;
chromosome_axis_unit_vectors_L = chromosome_axis_vectors_L ./ axis_lengths_L;

% Compute pairwise acute angles between chromosomes across the two nuclei
% Opposite directions are treated as equivalent by using abs(...)
cos_alpha_across = abs(chromosome_axis_unit_vectors_R * chromosome_axis_unit_vectors_L');

% Prevent numerical errors from giving values slightly outside [-1, 1]
cos_alpha_across = min(max(cos_alpha_across, -1), 1);

alpha_matrix_across = acosd(cos_alpha_across);

% Extract all 49 cross-nuclear pairwise angles
alpha_values_across = alpha_matrix_across(:);

% Mean pairwise angle across nuclei
mean_alpha_across = mean(alpha_values_across)

