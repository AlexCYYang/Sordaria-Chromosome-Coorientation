% Plot Inter-homolog Distance Profile for a homolog Pair
% This script demonstrates the calculation of inter-homolog distances
%
% Input: 
% 'C' - A 1x14 cell array containing 3D coordinates for 7 homologous pairs.
%       (e.g., C{1} and C{2} are pair 1; C{3} and C{4} are pair 2, etc.)

clearvars -except C
close all;

% === 1. Parameter Settings ===
targetChr = 2;          % Select which chromosome pair to plot (1 to 7)
expand = 7.33;          % Spatial expansion factor

% Determine cell array indices for the selected pair
idx1 = (targetChr * 2) - 1;
idx2 = (targetChr * 2);

% Extract the coordinates for the chosen pair
profile = {C{idx1}, C{idx2}};

if isempty(profile{1}) || isempty(profile{2})
    error('Selected chromosome pair is empty.');
end

% === 2. Assign equaldistant points along homologs ===
% This section resamples the raw traced paths into a set of 
% evenly spaced points (200 nm; pre-expansion scale). 
% This is required to compare distance of corresponding position between two homologs.
% Since the number of points need to be the same while the 
% length of homolog pairs might be slightly different,
% the point-to-point distance of homolog 2 needs to be slightly adjust 
% for the purpose of the same amount of total spots between the homologs

% ===== homolog 1 =====
v   = profile{1};
pos = []; spot= [];
homolog1 = v;
x = homolog1(:,1); y = homolog1(:,2); z = homolog1(:,3);

% Calculate the Euclidean distance between consecutive points in the raw trace
d = sqrt(diff(x).^2 + diff(y).^2 + diff(z).^2);

% Calculate the cumulative length (B) along the traced path
B = cumsum(d);
f  = max(B);   % Total physical length of homolog 1 
f1 = f;

% Determine fixed distance step sizes
l    = 0.2 * expand;   % Target segment length
seg  = round(f / l);   % Total number of equal segments for homolog 1
perl = l / expand;     % Physical length of each segment in micrometers
       
pos = [x(1), y(1), z(1)]; % Initialize the new path with the first raw point

% Interpolate new spots at exact intervals of 'l'
for ii = 1:seg-1
    % Find the closest existing raw point 'I' to our target cumulative distance (ii * l)
    % M is the absolute difference between the closest raw point and the target distance
    [M,I] = min(abs(B - ii * l));  
    
    if M == 0
        % If the target distance perfectly matches a raw point, use it directly
        spot(ii,:) = [x(ii+1), y(ii+1), z(ii+1)];
    else
        % Otherwise, interpolate a new point between the raw points.
        % 'r' calculates the fractional step direction to the adjacent point
        r = -abs(B(I) - ii*l) / (B(I) - ii*l); 
        
        % 'vv' is the 3D direction vector between the bounding raw points
        vv = [x(I+1+r) - x(I+1), y(I+1+r) - y(I+1), z(I+1+r) - z(I+1)];
        
        % Project the new 'spot' by scaling the direction vector 'vv' by the 
        % required distance 'M', and adding it to the base coordinates
        spot(ii,:) = vv .* (abs(M) / sqrt(sum(vv.*vv))) + [x(I+1), y(I+1), z(I+1)];
    end
    pos = [pos; spot(ii,:)]; % Append the newly interpolated spot
end
pos = [pos; x(end), y(end), z(end)]; % Append the final terminal point
h1  = pos; % Resampled coordinates for homolog 1

% ===== homolog 2 =====
v   = profile{2};
pos = []; spot= [];
homolog2 = v;
x = homolog2(:,1); y = homolog2(:,2); z = homolog2(:,3);

% Calculate cumulative distances for the second homolog trace
d = sqrt(diff(x).^2 + diff(y).^2 + diff(z).^2);
B = cumsum(d);
f  = max(B);   % Total physical length of homolog 2
f2 = f;
avg_length = (f1 + f2) / 2;    

% CRITICAL STEP: Re-calculate the step size 'l' based on the length of homolog 2
% divided by the EXACT SAME number of segments ('seg') derived from homolog 1.
% This forces homolog 1 and homolog 2 to have an identical number of points,
% enabling a 1-to-1 array comparison for inter-homolog distances.
l = f / seg;
pos  = [x(1), y(1), z(1)];

% The interpolation geometric logic here is identical to homolog 1
for ii = 1:seg-1
    [M,I] = min(abs(B - ii * l));
    if M == 0
        spot(ii,:) = [x(ii+1), y(ii+1), z(ii+1)];
    else
        r  = -abs(B(I) - ii*l) / (B(I) - ii*l);
        vv = [x(I+1+r) - x(I+1), y(I+1+r) - y(I+1), z(I+1+r) - z(I+1)];
        spot(ii,:) = vv .* (abs(M) / sqrt(sum(vv.*vv))) + [x(I+1), y(I+1), z(I+1)];
    end
    pos = [pos; spot(ii,:)];
end
pos = [pos; x(end), y(end), z(end)];
h2  = pos; % Resampled coordinates for homolog 2


% === 3. Inter-distance & Orientation Calculation ===
r1 = h1;
r2 = h2;
r2f = flip(r2); % flipped homolog 2 to test for reversed tracing direction

% Calculate 3D point-by-point distances for Original Orientation (dist1)
dr = r1 - r2;
vv_dist = dr .* dr;
lvv = sum(vv_dist, 2);
lv = sqrt(lvv);
dist1 = 1000 * lv / expand; % rescale back to pre-expansion and convert to nm

% Calculate 3D point-by-point distances for Reversed Orientation (dist2)
drf = r1 - r2f;
vvf = drf .* drf;
lvvf = sum(vvf, 2);
lvf = sqrt(lvvf);
dist2 = 1000 * lvf / expand; %rescale back to pre-expansion and convert to nm

% --- Orientation Selection Logic ---
% Tracing start/end points might be arbitrary. This corrects the alignment.
if targetChr == 2
    % For chromosome 2, the orientation is always the same: NU ends = the
    % first point in the imported coordinates.
    dist = dist1;
else
    % For other chromosomes, pick the orientation that yields the smaller 
    % overall mean distance, ensuring homologous ends are properly paired
    if mean(dist1) < mean(dist2)
        dist = dist1;
    else
        dist = dist2;
    end
end


% === 4. Plotting the Profile ===
% Generate the X-axis array
xp = 0 : perl : perl*(length(dist)-1);

% Create Plot
figure('Name', sprintf('Chromosome %d Profile', targetChr), 'Position', [100, 100, 800, 350]);
p = plot(xp, dist);
p.LineWidth = 4;
p.Color     = [0 0 0];

% Formatting
box off;
set(gca, 'TickDir', 'in', 'TickLength', [0.03 0.03], ...
         'LineWidth', 2.5, 'FontSize', 20, 'FontName', 'Arial');
     
xlabel('Position (\mum)', 'FontSize', 22, 'FontName', 'Arial');
ylabel('Inter-homolog distance (nm)', 'FontSize', 22, 'FontName', 'Arial');
title(sprintf('Chromosome %d Distance Profile', targetChr), 'FontSize', 24, 'FontName', 'Arial');

% Apply typical Y limits
ylim([-200 3500]); 

fprintf('Mean distance: %.2f nm\n', mean(dist));