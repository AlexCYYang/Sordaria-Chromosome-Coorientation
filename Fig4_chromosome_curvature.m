% Calculate average chromosome curvature (kap) for a single nucleus
%
% Input: 
% 'C' - A 1x7 cell array.  
% Each cell contains an N-by-3 matrix of XYZ coordinates 
% for one traced chromosome path in a nucleus.
%
% External function:
% This script requires 'curvature.m' and 'circumcenter.m' files (third-party functions).
% We use 'curvature.m' to calculate curvature, 
% and 'curvature.m' needs 'circumcenter.m' to help calculate curvature.
% Please ensure both are added to your MATLAB path.

clearvars -except C

% --- 1. Define Parameters ---
bead = 0.4;              % Inter-spot distance in microns (400 nm)
exp = 7.33;              % Expansion factor (spatial calibration)
nChr = length(C);        % Number of chromosomes (7 for pre-karyogamy)

% Initialize a master pool to store every single local curvature value
all_kap_nucleus = []; 

% --- 2. Loop through each chromosome ---
for i = 1:nChr
    
    % Extract the 3D coordinates
    c1 = C{i};
    x = c1(:, 1);
    y = c1(:, 2);
    z = c1(:, 3);
    
    % Calculate cumulative distance (B) along the traced path
    d = sqrt(diff(x).^2 + diff(y).^2 + diff(z).^2);
    B = cumsum(d);
    f_length = max(B); 
    
    % Calculate the number of 400 nm segments
    l = bead * exp; 
    seg = round(f_length / l); 
    
    % Skip if the chromosome is too short to have at least 2 segments
    if seg < 2
        continue;
    end
    
    % Initialize the resampled path with the starting point
    pos = [x(1), y(1), z(1)];
    spot = zeros(seg-1, 3);
    
    % --- Step 2.1: assign equal-distant spots along each chromosome ---
    for m = 1:(seg-1)
        [M_val, I_idx] = min(abs(B - m*l));
        
        if M_val == 0
            % If the target distance matches an existing point exactly
            spot(m, :) = [x(I_idx+1), y(I_idx+1), z(I_idx+1)];
        else
            % Calculate the direction vector toward the next/previous point
            r_val = -abs(B(I_idx) - m*l) / (B(I_idx) - m*l);
            
            % Direction vector v
            v_vec = [x(I_idx+1+r_val) - x(I_idx+1), ...
                     y(I_idx+1+r_val) - y(I_idx+1), ...
                     z(I_idx+1+r_val) - z(I_idx+1)];
            
            % Project the new spot at exactly distance m*l
            spot(m, :) = v_vec .* (abs(M_val) / sqrt(sum(v_vec.*v_vec))) + ...
                         [x(I_idx+1), y(I_idx+1), z(I_idx+1)];
        end
        pos = [pos; spot(m, :)];
    end
    % Append the final terminal point
    pos = [pos; x(end), y(end), z(end)];
    
    % --- Step 2.2: Local Curvature Calculation by the 3rd part curvature function---
    % Compute curvature on the equidistant resampled nodes
    [~, R, ~] = curvature(pos);
    RR = rmmissing(R);
    kap = 1 ./ RR;
    
    if isempty(kap)
        continue;
    end
    
    % Store all local kappa values into the nucleus master pool
    scaled_kap = kap * exp; %scale back to pre-expansion
    all_kap_nucleus = [all_kap_nucleus; scaled_kap(:)]; 
end

% --- 3. Final Calculation ---
% The average curvature is the mean of all local curvature values pooled together
nucleus_avg_curvature = mean(all_kap_nucleus, 'omitnan');

fprintf('Average curvature for this nucleus: %.4f\n', nucleus_avg_curvature);