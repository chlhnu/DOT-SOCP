function show_presets_dot2d(varargin)
%SHOW_PRESETS_DOT2D Display and demonstrate available problem and algorithm presets
%
%   show_presets_dot2d()              - Show all available presets
%   show_presets_dot2d('problems')    - Show problem types only
%   show_presets_dot2d('algorithms')  - Show algorithms only
%   show_presets_dot2d('usage')       - Show usage examples
%
% SEE ALSO: demo_dot2d, get_example, solver_dotsocp2d

if nargin == 0 || strcmp(varargin{1}, 'problems')
    show_problems();
end

if nargin == 0 || strcmp(varargin{1}, 'algorithms')
    show_algorithms();
end

if nargin == 0 || strcmp(varargin{1}, 'usage')
    show_usage();
end

end

function show_problems()
fprintf('\n');
fprintf('════════════════════════════════════════════════════════════════\n');
fprintf('Problem presets (set in demo_dot2d.m: Problem = "...")\n');
fprintf('════════════════════════════════════════════════════════════════\n\n');

fprintf('Basic:\n');
fprintf('  "example1"   - Single Gaussian to Single Gaussian\n');
fprintf('  "example2"   - Single Gaussian splits into four Gaussians\n');
fprintf('  "example3"   - Central Laplacian to four Gaussians\n');
fprintf('  "example4"   - Quartic polynomial to four Gaussians\n');
fprintf('  "example5"   - Centaur to Man\n');
fprintf('  "example6"   - Stitch 4 images from DOTmark dataset\n');
fprintf('  "example7"   - Gaussian to Dirac\n');
fprintf('  "circle"     - Circular densities\n\n');

fprintf('  Usage:\n');
fprintf('      [rho0, rho1] = get_example(Problem, nx, ny, lowerBound);\n');
fprintf('          lowerBound - Optional. Add lower bound to densities (default: 0)\n\n');

fprintf('────────────────────────────────────────────────────────────────\n\n');

fprintf('Additional notes on DOTmark dataset (example6):\n');

fprintf('  Description:\n');
fprintf('      Constructs rho0 and rho1 by stitching 4 images in 2×2 grid.\n');
fprintf('      Images are selected from DOTmark benchmark datasets.\n\n');

fprintf('  Extra parameters:\n');
fprintf('      DOTmark_type      - "ClassicImages" or "Shapes"\n');
fprintf('      stitch1_indices   - [i,j,k,l] where i,j,k,l ∈ [1,8] for rho0\n');
fprintf('      stitch2_indices   - [i,j,k,l] where i,j,k,l ∈ [1,8] for rho1\n\n');

fprintf('  Usage:\n');
fprintf('      [rho0, rho1] = get_example("example6", nx, ny, 0, ...\n');
fprintf('                     "DOTmark_type", "ClassicImages", ...\n');
fprintf('                     "stitch1_indices", [1,2,3,4], ...\n');
fprintf('                     "stitch2_indices", [5,6,7,8]);\n\n');
end

function show_algorithms()
fprintf('════════════════════════════════════════════════════════════════\n');
fprintf('Available algorithms\n');
fprintf('════════════════════════════════════════════════════════════════\n\n');

fprintf('  "inpALM"        - inexact proximal ALM (recommended)\n');
fprintf('  "pALM"          - proximal ALM\n');
fprintf('  "ALG2"          - Augmented Lagrangian Algorithm (ALG2) splitting scheme\n');
fprintf('  "sGS-inpALM"    - sGS-based inexact proximal ALM\n');
fprintf('  "acc-pADMM"      - Accelerated pADMM\n');
fprintf('  "acc-sGS-pADMM"  - Accelerated sGS-based pADMM\n\n');
end

function show_usage()
fprintf('════════════════════════════════════════════════════════════════\n');
fprintf('Quick usage guide\n');
fprintf('════════════════════════════════════════════════════════════════\n\n');

fprintf('1. Edit demo_dot2d.m, change these lines:\n');
fprintf('   Problem = "example1";   %% Try different problems\n');
fprintf('   algo = "inpALM";        %% Try different algorithms\n\n');

fprintf('2. Run:   >> demo_dot2d\n\n');

fprintf('3. Common parameters to adjust:\n');
fprintf('       nt = 33;       %% Time resolution (2^k+1 for multilevel coarsening)\n');
fprintf('       nx = 129;      %% Spatial resolution (2^k+1 for multilevel coarsening)\n');
fprintf('       levelN = 3;    %% Multilevel depth\n');
fprintf('   Note: Grid sizes of form 2^k+1 (k >= levelN) are only convenient for multilevel coarsening,\n');
fprintf('         which can be overcome with more flexible grid handling.\n\n');
end
