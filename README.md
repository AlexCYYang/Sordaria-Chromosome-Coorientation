Supplementary Code for Data Analysis

This file describes the MATLAB scripts used for quantitative analyses in this study.

Software: MATLAB (Release 2024b, MathWorks).

File descriptions:
---------------------------------------------------------------------------

1. bounding_box_long_axis_vectors.m
– Fit each chromosome with minimum bounding box and extract the long-axis vectors of all chromosomes for coorientation calculation 
– The output of this script is the input for Fig2_4_coorientation_mean_alpha_within_nucle.m
- Note: This script requires the external function 'minboundbox.m' described below.


2. Fig2_4_coorientation_mean_alpha_within_nucle.m
- Computes the intra-nuclear chromosome coorientation (mean pairwise angle alpha) for Fig. 2 and Fig. 4.

3. Fig2_chromosome_space_ellipticity.m
- Calculates the ellipticity of the chromosome space for Fig. 2.
- Note: This script requires the external function 'ellipsoid_fit_new.m' described below.

4. Fig3_coorientation_mean_alpha_across_nuclei.m
- Computes the inter-nuclear chromosome coorientation across pair nuclei (mean pairwise angle alpha) for Fig. 3.

5. Fig3_interface_area.m
- Quantifies the nuclear interface area in Fig. 3. 

6. Fig4_chromosome_shape_elongation.m
– Quantifies the chromosome shape elongation (EL) in Fig. 4. 
- Note: This script requires the external function 'minboundbox.m' described below.

7. Fig4_chromosome_curvature.m
– Quantifies the chromosome curvature (average per nucleus) in Fig. 4. 
- Note: This script requires the external function 'curvature.m' and 'circumcenter.m' (both from the same source)described below.

8. Fig5_interhomolog_distance_profile
– Measure and plot the distance between homologous chromosomes along their lengths.


External functions:
---------------------------------------------------------------------------
The custom scripts above utilize the following third-party functions. These files are NOT included in this repository/archive. To execute the analyses, please download the required external functions listed below from the MATLAB Central File Exchange. Ensure that all custom scripts and downloaded external functions are placed within the same active MATLAB directory or added to your MATLAB path before running. For ellipsoid fitting, please rename the downloaded file 'ellipsoid_fit.m' to 'ellipsoid_fit_new.m' to match the function call used in this workflow.

1. ellipsoid_fit.m
- Function: ellipsoid fitting.
- Source: 
Yury (2026). Ellipsoid fit (https://www.mathworks.com/matlabcentral/fileexchange/24693-ellipsoid-fit), MATLAB Central File Exchange. Retrieved April 5, 2026.
- Note: In our workflow, the original file 'ellipsoid_fit.m' was renamed to 'ellipsoid_fit_new.m'. Please download the original function from the source above and rename it as 'ellipsoid_fit_new.m' before running the script.

2. minboundbox.m
- Function: Minimum bounding box calculation.
- Source: Johannes Korsawe (2026). Minimal Bounding Box (https://www.mathworks.com/matlabcentral/fileexchange/18264-minimal-bounding-box), MATLAB Central File Exchange. Retrieved April 6, 2026.
- Note: 'minboundbox.m' must be downloaded from the source package and added to the MATLAB path.


3. curvature.m and circumcenter.m
- Function: 3D local curvature calculation.
- Source: Are Mjaavatten (2026). Curvature of a 1D curve in a 2D or 3D space (https://www.mathworks.com/matlabcentral/fileexchange/69452-curvature-of-a-1d-curve-in-a-2d-or-3d-space), MATLAB Central File Exchange. Retrieved April 10, 2026.
- Note: Both 'curvature.m' and its helper function 'circumcenter.m' must be downloaded from the source package and added to the MATLAB path.

Other scripts:
---------------------------------------------------------------------------
This repository may also contain additional custom scripts in the 'other' folder (e.g., for 3D chromosome visualization, rendering, or movie generation). These files were not used for the quantitative analyses reported in this study.
