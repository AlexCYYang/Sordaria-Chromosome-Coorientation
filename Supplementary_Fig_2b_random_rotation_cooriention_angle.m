% =========================================================================
% Within-Nucleus Randomization (Random Rotation of Traced Chromosomes)
% Keeps chromosome shapes fixed but independently rotates each chromosome
% in 3D around its arc-length-weighted centroid. Computes the mean of 
% all pairwise acute angles to establish a randomized baseline.
%
% Note: Ensure the cell array 'C' (containing the Nx3 polylines for the 
% 7 chromosomes) is loaded in the workspace before running this script.
% =========================================================================

% =========================================================================
% CONFIGURATION
% =========================================================================
CFG_fontName  = 'Arial';
CFG_fontSize  = 20;
CFG_titleSize = 20;

% Histogram bin control ('auto' | 'binWidth' | 'numBins')
CFG_binMode   = 'binWidth';     
CFG_binWidth  = 1;              % degrees, used if binMode='binWidth'
CFG_numBins   = 30;             % used if binMode='numBins'

% Axis limits (leave [] for auto)
CFG_xLim      = [22.9 69.1];    
CFG_yLim      = [0 0.1];        

% =========================================================================
% INPUT & INITIALIZATION
% =========================================================================
% CC should be a 1x7 (or 7x1) cell array, each cell is an Nx3 polyline
CC = C;                 
nChr = numel(CC);
nRep = 10000;            
axisMethod = 'pca';     % 'pca' or 'end2end'

meanAngles = zeros(nRep, 1);

% =========================================================================
% SIMULATION LOOP
% =========================================================================
for rep = 1:nRep
    U = zeros(nChr, 3);  % One unit direction vector per chromosome
    
    for i = 1:nChr
        P = CC{i};                  % Nx3 trace points
        CM = arclen_centroid(P);    % 1x3 arc-length-weighted centroid
        
        % Random 3D rotation (isotropic)
        R = randrot3();             % 3x3 rotation matrix, det=+1
        
        % Rigid rotation about centroid
        Prot = (R * (P - CM).').' + CM;
        
        % Extract direction vector for this chromosome
        u = chr_axis(Prot, axisMethod);   
        U(i, :) = u / norm(u);
    end
    
    % Pairwise angles (fold to 0-90 degrees)
    D = U * U.';                           % Dot products
    dots = D(triu(true(nChr), 1));         % Extract 21 unique pairs
    dots = max(-1, min(1, dots));          % Numeric safety clamp
    ang  = acosd(abs(dots));               % Fold to [0, 90]
    
    meanAngles(rep) = mean(ang);
end

% =========================================================================
% STATISTICS
% =========================================================================
mu = mean(meanAngles);
sd = std(meanAngles, 0);

fprintf('=== Simulation Results (Independent random rotations, nChr=%d) ===\n', nChr);
fprintf('Mean of mean-angles = %.4f deg\n', mu);
fprintf('SD across repeats   = %.4f deg\n\n', sd);

% =========================================================================
% VISUALIZATION 1: Errorbar Plot
% =========================================================================
figure('Name', 'Mean and SD'); hold on;
errorbar(1, mu, sd, 'o', 'LineWidth', 1.5);
xlim([0.5 1.5]);
xticks(1);
xticklabels({'Independent rotation'});
ylabel('Mean pairwise angle (deg, folded 0-90)', 'FontName', CFG_fontName, 'FontSize', CFG_fontSize);
set(gca, 'FontName', CFG_fontName, 'FontSize', CFG_fontSize, 'LineWidth', 1);
grid on;

% =========================================================================
% VISUALIZATION 2: Histogram Distribution
% =========================================================================
ctrl_mean = mean(meanAngles);
ctrl_sd   = std(meanAngles, 0);

figure('Name', 'Randomized Distribution'); hold on;

% Plot histogram (PDF-normalized)
switch CFG_binMode
    case 'auto'
        histogram(meanAngles, 'Normalization', 'pdf', 'FaceColor', [0.7 0.7 0.7]);
    case 'binWidth'
        histogram(meanAngles, 'BinWidth', CFG_binWidth, 'Normalization', 'pdf', 'FaceColor', [0.7 0.7 0.7]);
    case 'numBins'
        histogram(meanAngles, CFG_numBins, 'Normalization', 'pdf', 'FaceColor', [0.7 0.7 0.7]);
    otherwise
        error('CFG_binMode must be ''auto'', ''binWidth'', or ''numBins''.');
end

% Reference lines: mean ± 1 SD
xline(ctrl_mean, 'r',  'LineWidth', 2);
xline(ctrl_mean - ctrl_sd, 'r--', 'LineWidth', 1.5);
xline(ctrl_mean + ctrl_sd, 'r--', 'LineWidth', 1.5);

% Labels & Formatting
xlabel('Mean pairwise angle (deg, folded 0-90)', 'FontName', CFG_fontName, 'FontSize', CFG_fontSize);
ylabel('Probability density', 'FontName', CFG_fontName, 'FontSize', CFG_fontSize);

title(sprintf('Randomized (N=%d): %.2f \\pm %.2f^\\circ', numel(meanAngles), ctrl_mean, ctrl_sd), ...
    'FontName', CFG_fontName, 'FontSize', CFG_titleSize);

set(gca, 'FontName', CFG_fontName, 'FontSize', CFG_fontSize, 'LineWidth', 1);

if ~isempty(CFG_xLim); xlim(CFG_xLim); end
if ~isempty(CFG_yLim); ylim(CFG_yLim); end

grid off;
box on;

% =========================================================================
% HELPER FUNCTIONS
% =========================================================================

function CM = arclen_centroid(P)
% Arc-length-weighted centroid of a polyline P (Nx3)
% CM = sum_k (d_k * midpoint_k) / sum_k d_k
    dP = diff(P, 1, 1);                 % (N-1)x3
    d  = sqrt(sum(dP.^2, 2));           % (N-1)x1
    M  = (P(1:end-1, :) + P(2:end, :)) / 2; % (N-1)x3
    CM = sum(M .* d, 1) / sum(d);       % 1x3
end

function u = chr_axis(P, method)
% Return a 1x3 direction vector for a chromosome polyline
    switch lower(method)
        case 'end2end'
            u = P(end, :) - P(1, :);
        case 'pca'
            % PCA major axis (direction only)
            X = P - mean(P, 1);
            [~, ~, V] = svd(X, 'econ');  % columns of V are principal directions
            u = V(:, 1).';               % 1x3
        otherwise
            error('Unknown axis method: %s', method);
    end
    if norm(u) < 1e-12
        error('Axis vector is near-zero. Check input trace.');
    end
end

function R = randrot3()
% Isotropic random 3D rotation matrix (det=+1), Shoemake quaternion method
    u1 = rand();
    u2 = rand();
    u3 = rand();
    
    q1 = sqrt(1 - u1) * sin(2 * pi * u2);
    q2 = sqrt(1 - u1) * cos(2 * pi * u2);
    q3 = sqrt(u1)     * sin(2 * pi * u3);
    q4 = sqrt(u1)     * cos(2 * pi * u3);   % (q1, q2, q3, q4)
    
    % Convert quaternion to rotation matrix
    x = q1; y = q2; z = q3; w = q4;
    R = [1 - 2*(y^2 + z^2),     2*(x*y - z*w),     2*(x*z + y*w);
             2*(x*y + z*w), 1 - 2*(x^2 + z^2),     2*(y*z - x*w);
             2*(x*z - y*w),     2*(y*z + x*w), 1 - 2*(x^2 + y^2)];
end