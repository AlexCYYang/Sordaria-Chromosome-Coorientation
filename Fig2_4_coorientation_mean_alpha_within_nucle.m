clearvars -except out;

load('B_cube_long_only.mat'); % Load a 7x3 matrix (ve):
% each row is the long-axis direction vector of one chromosome,
% obtained from the minimum bounding box.

chromosome_axis_vectors = ve;

% Normalize each chromosome axis vector to unit length
axis_lengths = sqrt(sum(chromosome_axis_vectors.^2, 2));
chromosome_axis_unit_vectors = chromosome_axis_vectors ./ axis_lengths;

% Compute pairwise acute angles (alpha) between chromosome axes
cos_alpha = abs(chromosome_axis_unit_vectors * chromosome_axis_unit_vectors');
cos_alpha = min(max(cos_alpha, -1), 1);
alpha_matrix = acosd(cos_alpha);

% Extract the 21 unique pairwise angles
alpha_values = [];

for j = 1:6
    alpha_values = [alpha_values; alpha_matrix(j+1:7, j)];
end

% Mean pairwise angle
mean_alpha = mean(alpha_values);

fprintf('Mean alpha = %.2f degrees\n', mean_alpha);

