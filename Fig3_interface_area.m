% Final interface-area pipeline used for the paper
% This script:
%   1) reads two traced NE volumes in the current folder
%   2) converts them to 3D point clouds
%   3) defines interface points using a single threshold from the combined
%      bidirectional nearest-neighbor distance distribution
%   4) estimates interface area on each side using clustered planar
%      triangulation (DBSCAN + PCA + Delaunay + alpha-shape mask)
%   5) reports A-side area, B-side area, and their mean

clear;
close all;

load('obl_info.mat');

%% 1. Read NE volumes and convert them to point clouds
files = dir(fullfile(pwd, '*NE.tif'));
names = string({files.name});

% Assign the two surfaces
A_filename = names(2);
B_filename = names(1);

A = [];
B = [];
for j = 1:tz
    A(:,:,j) = imread(A_filename, j);
    B(:,:,j) = imread(B_filename, j);
end

% Downsampling factor (kept as 1 for final analysis)
f = 1;

Asm = imresize3(A, f);
Bsm = imresize3(B, f);

[rA, cA, vA] = ind2sub(size(Asm), find(Asm > 1));
[rB, cB, vB] = ind2sub(size(Bsm), find(Bsm > 1));

px_use = px / f;
py_use = px_use;
zf = size(A,3) / size(Asm,3);
pz_use = pz * zf;

xA = cA * px_use;
yA = rA * py_use;
zA = (vA - 1) * pz_use;

xB = cB * px_use;
yB = rB * py_use;
zB = (vB - 1) * pz_use;

ptsA = [xA(:), yA(:), zA(:)];
ptsB = [xB(:), yB(:), zB(:)];

%% 2. Compute bidirectional nearest-neighbor distances
pct_limit = 15.0;
num_bins = 100;
skip_first_peak_if_alternatives = true;

fprintf('Computing bidirectional nearest-neighbor distances...\n');
[~, D_AtoB] = knnsearch(ptsB, ptsA);
[~, D_BtoA] = knnsearch(ptsA, ptsB);

%% 3. Define a single interface threshold from the combined distance distribution
D_all = [D_AtoB(:); D_BtoA(:)];

[counts_all, edges_all] = histcounts(D_all, num_bins);
bin_centers_all = edges_all(1:end-1) + diff(edges_all)/2;

T_limit = prctile(D_all, pct_limit);

counts_all = counts_all(:);

internal_peak_idx = find(counts_all(2:end-1) > counts_all(1:end-2) & ...
                         counts_all(2:end-1) > counts_all(3:end)) + 1;

start_peak_idx = [];
if numel(counts_all) > 1 && counts_all(1) > counts_all(2)
    start_peak_idx = 1;
end

end_peak_idx = [];
if numel(counts_all) > 1 && counts_all(end) > counts_all(end-1)
    end_peak_idx = numel(counts_all);
end

peak_idx = unique([start_peak_idx; internal_peak_idx; end_peak_idx]);
peak_distances = bin_centers_all(peak_idx);

valid_peak_idx = peak_idx(peak_distances <= T_limit);

if skip_first_peak_if_alternatives && any(valid_peak_idx > 1)
    valid_peak_idx = valid_peak_idx(valid_peak_idx > 1);
end

if isempty(valid_peak_idx)
    error('No valid peak was found within the lower %.1f%% of the combined distance distribution.', pct_limit);
end

final_peak_idx = valid_peak_idx(1);
T_interface = edges_all(final_peak_idx + 1);

fprintf('Interface threshold = %.4f\n', T_interface);

%% 4. Identify interface points on both sides
ptsA_interface = ptsA(D_AtoB <= T_interface, :);
ptsB_interface = ptsB(D_BtoA <= T_interface, :);

fprintf('A-side interface points: %d\n', size(ptsA_interface, 1));
fprintf('B-side interface points: %d\n', size(ptsB_interface, 1));

%% 5. Estimate interface area by clustered planar triangulation
exm_factor = 7.33;
area_scale = exm_factor^2;

areaA = clustered_interface_area(ptsA_interface, T_interface, area_scale);
areaB = clustered_interface_area(ptsB_interface, T_interface, area_scale);

%% 6. Report results
fprintf('A-side interface area: %.6f\n', areaA);
fprintf('B-side interface area: %.6f\n', areaB);

if areaA > 0 && areaB > 0
    interface_area_pair = mean([areaA, areaB]);
    fprintf('Pair-averaged interface area: %.6f\n', interface_area_pair);
else
    interface_area_pair = NaN;
    fprintf('Pair-averaged interface area: NaN\n');
end


function total_area = clustered_interface_area(interface_pts, dist_threshold, area_scale)
% Estimate interface area from clustered interface points.
% Each spatial cluster is projected onto its own best-fit plane.
% Area is measured in the projected plane using Delaunay triangulation,
% restricted to the cluster footprint by an alpha-shape mask.

    total_area = 0;

    if isempty(interface_pts) || size(interface_pts, 1) < 3
        return;
    end

    eps_dbscan = max(1e-6, 1.5 * dist_threshold);
    min_pts = max(15, round(0.005 * size(interface_pts, 1)));

    labels = dbscan(interface_pts, eps_dbscan, min_pts);
    cluster_ids = unique(labels(labels > 0));

    if isempty(cluster_ids)
        return;
    end

    tri_area_sum_2d = @(XY, tri_idx) ...
        sum(abs( ...
            (XY(tri_idx(:,2),1) - XY(tri_idx(:,1),1)) .* (XY(tri_idx(:,3),2) - XY(tri_idx(:,1),2)) - ...
            (XY(tri_idx(:,3),1) - XY(tri_idx(:,1),1)) .* (XY(tri_idx(:,2),2) - XY(tri_idx(:,1),2)) ...
        ) / 2);

    for k = 1:numel(cluster_ids)
        Pk = interface_pts(labels == cluster_ids(k), :);

        if size(Pk, 1) < 3
            continue;
        end

        [~, score_k] = pca(Pk);
        XY = score_k(:, 1:2);

        try
            shp = alphaShape(XY(:,1), XY(:,2));
            shp.Alpha = max(eps, 0.9 * criticalAlpha(shp));
        catch
            shp = alphaShape(XY(:,1), XY(:,2));
        end

        tri = delaunay(XY(:,1), XY(:,2));
        if isempty(tri)
            continue;
        end

        tri_centers = (XY(tri(:,1),:) + XY(tri(:,2),:) + XY(tri(:,3),:)) / 3;
        keep_tri = inShape(shp, tri_centers(:,1), tri_centers(:,2));
        tri_in = tri(keep_tri, :);

        if isempty(tri_in)
            continue;
        end

        cluster_area = tri_area_sum_2d(XY, tri_in) / area_scale;
        total_area = total_area + cluster_area;
    end
end