% =========================================================================
% Random 3D Vector Coorientation Simulation
% Generates sets of random 3D vectors and calculates the mean pairwise 
% acute angles to establish a random expectation for coorientation.
% =========================================================================
clear; close all; clc;

% =========================================================================
% CONFIGURATION
% =========================================================================
CFG.n_simulations = 10000;   % Monte Carlo iterations
CFG.n_vectors     = 7;       % Number of vectors per set

% Font and plot settings
CFG.font_name  = 'Arial';
CFG.font_size  = 20;
CFG.title_size = 20;
CFG.line_width = 2;

% Plot mode: 
%   'mean_per_set'    - Histogram of the mean of pairwise angles per set
%   'pairwise_angles' - Histogram of all pairwise angles across all sets
CFG.plot_mode = 'mean_per_set';

% Histogram bin control ('auto', 'binWidth', or 'numBins')
CFG.bin_mode  = 'auto';
CFG.bin_width = 0.25;      % Used if bin_mode = 'binWidth'
CFG.num_bins  = 50;        % Used if bin_mode = 'numBins'

% Axis limits (leave [] for auto)
CFG.x_lim = [];            
CFG.y_lim = [];            

% =========================================================================
% SIMULATION
% =========================================================================
n_simulations = CFG.n_simulations;
n_vectors     = CFG.n_vectors;

% Pre-allocate storage
if strcmp(CFG.plot_mode, 'mean_per_set')
    data = zeros(n_simulations, 1);
elseif strcmp(CFG.plot_mode, 'pairwise_angles')
    n_pairs = n_vectors * (n_vectors - 1) / 2;
    data = zeros(n_simulations * n_pairs, 1);
    idx = 1;
else
    error('CFG.plot_mode must be ''mean_per_set'' or ''pairwise_angles''.');
end

% Monte Carlo loop
for i = 1:n_simulations
    % 1) Generate random isotropic 3D vectors and normalize
    vecs = randn(n_vectors, 3);
    vecs = vecs ./ vecnorm(vecs, 2, 2);

    % 2) Calculate pairwise dot products
    dot_matrix = vecs * vecs';

    % 3) Extract unique pairs (upper triangular, excluding diagonal)
    mask = triu(true(n_vectors), 1);
    pairwise_dots = dot_matrix(mask);

    % Numeric safety clamp to avoid complex numbers in acosd
    pairwise_dots = max(-1, min(1, pairwise_dots));

    % 4) Convert to acute angles (0-90 degrees)
    pairwise_angles = acosd(abs(pairwise_dots));

    % Store results based on selected mode
    if strcmp(CFG.plot_mode, 'mean_per_set')
        data(i) = mean(pairwise_angles);
    else
        data(idx:idx+numel(pairwise_angles)-1) = pairwise_angles;
        idx = idx + numel(pairwise_angles);
    end
end

% =========================================================================
% STATISTICS
% =========================================================================
mu = mean(data);
sd = std(data);

fprintf('=== Simulation Results ===\n');
fprintf('Mode: %s\n', CFG.plot_mode);
fprintf('Mean: %.4f degrees\n', mu);
fprintf('SD  : %.4f degrees\n\n', sd);

% =========================================================================
% VISUALIZATION
% =========================================================================
figure;

% Generate histogram
switch CFG.bin_mode
    case 'auto'
        h = histogram(data, 'Normalization', 'pdf');
    case 'binWidth'
        h = histogram(data, 'BinWidth', CFG.bin_width, 'Normalization', 'pdf');
    case 'numBins'
        h = histogram(data, CFG.num_bins, 'Normalization', 'pdf');
    otherwise
        error('CFG.bin_mode must be ''auto'', ''binWidth'', or ''numBins''.');
end

h.FaceColor = [0.7 0.7 0.7];
hold on;

% Add Mean ± 1 SD indicators
xline(mu,      'r',  'LineWidth', CFG.line_width);
xline(mu - sd, 'r--', 'LineWidth', 1.5);
xline(mu + sd, 'r--', 'LineWidth', 1.5);

% Formatting labels
if strcmp(CFG.plot_mode, 'mean_per_set')
    xlabel_txt = 'Mean pairwise acute angle per set (degrees)';
else
    xlabel_txt = 'Pairwise acute angle (degrees)';
end

xlabel(xlabel_txt, 'FontName', CFG.font_name, 'FontSize', CFG.font_size);
ylabel('Probability density', 'FontName', CFG.font_name, 'FontSize', CFG.font_size);

title(sprintf('%s (N=%d): %.2f \\pm %.2f^\\circ', ...
    strrep(CFG.plot_mode, '_', ' '), n_simulations, mu, sd), ...
    'FontName', CFG.font_name, 'FontSize', CFG.title_size);

set(gca, 'FontName', CFG.font_name, 'FontSize', CFG.font_size, 'LineWidth', 1);

if ~isempty(CFG.x_lim); xlim(CFG.x_lim); end
if ~isempty(CFG.y_lim); ylim(CFG.y_lim); end

% =========================================================================
% OUTPUT DISTRIBUTION DATA
% =========================================================================
% Extract bin centers and probability density for potential export
hist_x = h.BinEdges(1:end-1) + diff(h.BinEdges)/2;   
hist_y = h.Values;                                   

hist_xy = table(hist_x(:), hist_y(:), 'VariableNames', {'BinCenter_deg', 'ProbDensity'});
disp(head(hist_xy)); % Display only the top few rows in command window to avoid clutter